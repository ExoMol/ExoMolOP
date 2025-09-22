
import numpy as np
try:
    import pickle
except:
    import cPickle as pickle
import os
import sys
import imp
import glob
import multiprocessing
from multiprocessing import Pool
import argparse
import h5py

def get_specgrid( R=10000, lambda_min=0.3, lambda_max=15.0):
    #generating wavelength grid with uniform binning in log(lambda)
    #lambda min and max are in microns, R is the spectral resolution
    #R = lambda/delta_lamdba
    specgrid = []
    delta_lambda =[]
    specgrid.append(lambda_min)
    run = True
    i=1
    while run:
        dlam= specgrid[i-1]/R
        specgrid.append(specgrid[i-1]+dlam)
        delta_lambda.append(dlam)

        if specgrid[i] >= lambda_max:
            run=False
        i+=1
    return np.asarray(specgrid),np.asarray(delta_lambda)

def get_specbingrid(wavegrid, specgrid, binwidths=None):
    #function calculating the bin boundaries for the data
    #this is used to bin the internal spectrum to the data in fitting module

    if not isinstance(binwidths, (np.ndarray, np.generic)):
        bingrid =[]
        bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
        for i in range(len(wavegrid)-1):
            bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
        bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
        bingrid_idx = np.digitize(specgrid,bingrid) #getting the specgrid indexes for bins
    else:

        # this bingrid is actually useless, as it doesn't allow for gaps in the data
        bingrid = []
        for i in range(len(wavegrid)):
            bingrid.append(wavegrid[i]-binwidths[i]/2.)
            bingrid.append(wavegrid[i]+binwidths[i]/2.)

        # build bin grid index array (an index for each model datapoint. If the point is outside the input
        # spectrum bins, the idx is -99999 (invalid integer...)
        bingrid_idx = np.empty(len(specgrid))
        bingrid_idx[:] = np.NaN

        for i in range(len(specgrid)):
            for j in range(len(wavegrid)):
                if specgrid[i] >= (wavegrid[j]-binwidths[j]/2.) and specgrid[i] < (wavegrid[j]+binwidths[j]/2.):
                    bingrid_idx[i] = j+1
                    break

    return bingrid, bingrid_idx



#loading parameter file parser
parser = argparse.ArgumentParser()
parser.add_argument('--dict', '--dictionary',
                    dest='dictionary',
                    default=False
                   )
parser.add_argument('--resolution',
                    dest='resolution',
                    type=int,
                    default=False
                    )
parser.add_argument('--wl',
                    dest='wlrange',
                    type=str,
                    default=False
                   )
parser.add_argument('--key_iso_ll',
                    dest='key_iso_ll',
                    type=str,                                       
                    default=False                                   
                    )
parser.add_argument('--mname',
                    dest='mname',
                    type=str,
                    default=False
                    )
parser.add_argument('--mol_mass',
                    dest='mol_mass',
                    type=int,
                    default=1
                    ) 
options = parser.parse_args()


## Import input dictionary
try:
    mol = imp.load_source('molecule', options.dictionary)
    define = mol.define
except:
    print(' Cannot import dictionary file')
    exit()

mol_name = define['molecule_name']
ncores = 24 # equivalent to number of temps

# determine lambda range and create interpolation grid

lambdamin_str, lambdamax_str = options.wlrange.split(',')
lambdamin = float(lambdamin_str)
lambdamax = float(lambdamax_str)
wlgrid, dlamb_grid = get_specgrid(R=options.resolution,lambda_min=lambdamin,lambda_max=lambdamax)
wngrid = np.sort(10000./wlgrid)

# input and output folders
folder_xsec = os.path.join(define['working_folder'], 'xsec_combined_taurex')
output_folder = os.path.join(define['working_folder'], 'xsec_hdf5_sampled_R%i_%s-%s_final_v1' % (options.resolution, lambdamin_str, lambdamax_str))

# create dir
os.system('mkdir -p %s' % output_folder)

def worker(temp_idx):
    # worker for a given temperature idx
    temp_val = define['temp_list'][temp_idx]
    print('Doing T', temp_val)
    xsec_file =  os.path.join(folder_xsec, '%s_T%i.TauREx.pickle' % (mol_name, temp_val))

    sigma_in = pickle.load(open(xsec_file, 'rb'), encoding='latin1')

    sigma_array = np.zeros((len(define['press_list']), 1, len(wngrid)))

    for press_idx, press_val in enumerate(np.sort(define['press_list'])):
        sigma_array[press_idx, 0] =  np.interp(wngrid, sigma_in['wno'],
                                                 sigma_in['xsecarr'][press_idx, 0],
                                                 left=0, right=0)
    sigma_out = {
        'name': sigma_in['name'],
        'p': sigma_in['p'],
        't': sigma_in['t'],
        'wno': wngrid,
        'xsecarr': sigma_array,
    }
    pickle.dump(sigma_out, open(os.path.join(output_folder, '%s_T%i.TauREx.pickle' %
                                             (mol_name, temp_val)), 'wb'), protocol=2)
    print('Done T', temp_val)

pool = Pool(processes=int(ncores)) #setting number of cores on which to run
pool_result = pool.map(worker,[temp_idx for temp_idx, temp_val in enumerate(define['temp_list'])])
pool.close()                            #closing the pool
pool.join()                             #closing the pool. Needs to be done otherwise createob.reset has no effect.

# aggregate xsec into a single file

print('Aggregate xsec into a single file')

sigma_array = np.zeros((len(define['press_list']), len(define['temp_list']), len(wngrid)))

for temp_idx, temp_val in enumerate(np.sort(define['temp_list'])):
    xsec_file =  os.path.join(output_folder, '%s_T%i.TauREx.pickle' % (mol_name, temp_val))
    sigma_in = pickle.load(open(xsec_file, 'rb'), encoding='latin1')['xsecarr']
    sigma_array[:, temp_idx, :] = sigma_in[:, 0, :]
    os.system('rm %s' % xsec_file)

sigma_out = {
    'name': mol_name,
    'p': np.sort(define['press_list']).astype(float),
    't': np.sort(define['temp_list']).astype(float),
    'wno': wngrid,
    'xsecarr': sigma_array,
}


hdf5_out = h5py.File(os.path.join(output_folder, '%s.R%i_%s-%smu.xsec.TauREx.h5' %
                                         (options.key_iso_ll, options.resolution,
                                           lambdamin_str, lambdamax_str)),'w')

str_type = h5py.new_vlen(str)
dset = hdf5_out.create_dataset("mol_name",(1,), dtype=str_type)
dset2 = hdf5_out.create_dataset("key_iso_ll",(1,), dtype=str_type)
dset3 = hdf5_out.create_dataset("DOI",(1,), dtype=str_type)
dset4 = hdf5_out.create_dataset("Date_ID",(1,), dtype=str_type)

dset[0] = options.mname
dset2[0] = options.key_iso_ll
hdf5_out['p'] = np.sort(define['press_list']).astype(float)
hdf5_out['p'].attrs['units'] = 'bar'
hdf5_out['t'] = np.sort(define['temp_list']).astype(float)
hdf5_out['t'].attrs['units'] = 'kelvin'
hdf5_out['bin_edges'] = wngrid
hdf5_out['bin_edges'].attrs['units'] = 'wavenumbers'
hdf5_out['xsecarr'] = sigma_array
hdf5_out['xsecarr'].attrs['units'] = 'cm^2/molecule'
dset3[0] = 'rep4'
dset4[0] = 'rep5'
hdf5_out['mol_mass'] = options.mol_mass


hdf5_out.close()
#pickle.dump(sigma_out, open(os.path.join(output_folder, '%s.R%i.TauREx.pickle' %
#                                         (mol_name, options.resolution)), 'wb'), protocol=2)
print('End')
