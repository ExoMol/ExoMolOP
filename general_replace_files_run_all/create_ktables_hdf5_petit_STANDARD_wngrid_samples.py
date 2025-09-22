'''

Create ktables from high resolution cross sections

'''

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

def get_specgrid( R=5000, lambda_min=0.1, lambda_max=20.0):
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
parser.add_argument('--dict', '--dictionary', # definition dictionary, usually MOL_input.py
                    dest='dictionary',
                    default=False
                   )
parser.add_argument('--input_grid',
                    dest='input_grid',
                    type=str,
                    default=False
                    )
parser.add_argument('--title',
                    dest='title',
                    type=str,
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
parser.add_argument('--ncores',  # number of cores to use
                    dest='ncores',
                    type=int,
                    default=1
                   )
parser.add_argument('--key_iso_ll',
                    dest='key_iso_ll',
                    type=str,
                    default=False
                    )
parser.add_argument('--mol_mass',
                    dest='mol_mass',
                    type=int,
                    default=1
                    )  
parser.add_argument('--mname',
                    dest='mname',
                    type=str,
                    default=False
                    )
parser.add_argument('--DOI',
                    dest='DOI',
                    type=str,
                    default=False
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
ncores = options.ncores # equivalent to number of temps

if options.resolution:

    # determine lambda range and create binning grid
    lambdamin_str, lambdamax_str = options.wlrange.split(',')
    lambdamin = float(lambdamin_str)
    lambdamax = float(lambdamax_str)
    wl_grid = get_specgrid(int(options.resolution), float(lambdamin), float(lambdamax))[0] # build bin grid in wavelength
    wn_grid = np.sort(10000./wl_grid) # convert to wavenumber and sort
    bincentres = np.sort(10000./(0.5*(wl_grid[1:] + wl_grid[:-1]))) # get the bin centres in wl space

elif options.input_grid:

    wn_grid = np.sort(np.loadtxt(options.input_grid)) # this are the bin edges
    bincentres = np.sort((0.5*(wn_grid[1:] + wn_grid[:-1]))) # bin centres in wn
    wl_grid = np.sort(10000./wn_grid) # bin edges in wn
    lambdamin = np.min(wl_grid)
    lambdamax = np.max(wl_grid)


    #wl_grid = np.sort(np.loadtxt(options.input_grid)) # this are the bin edges
    #bincentres = np.sort(10000./(0.5*(wl_grid[1:] + wl_grid[:-1]))) # bin centres in wn
    #wn_grid = np.sort(10000./wl_grid) # bin edges in wn
    #lambdamin = np.min(wl_grid)
    #lambdamax = np.max(wl_grid)

# get gaussian quadrature points points
#ngauss = int(options.ngauss)
#gauss = np.polynomial.legendre.leggauss(ngauss) # get the legendre-gauss sample points and weights
#samples = gauss[0]
#weights = gauss[1]



samples = np.array([-0.964260871,
-0.81699983,
-0.572979169,
-0.265091178,
0.065091178,
0.372979169,
0.61699983,
0.764260871,
0.803971014,
0.820333352,
0.847446759,
0.881656536,
0.918343464,
0.952553241,
0.979666648,
0.996028986])

ngauss=16

weights = np.array([0.091105683,
0.200142931,
0.282335981,
0.326415405,
0.326415405,
0.282335981,
0.200142931,
0.091105683,
0.010122854,
0.022238103,
0.031370665,
0.036268378,
0.036268378,
0.031370665,
0.022238103,
0.010122854])

gauss = (samples,weights)

# input and output folders
folder_xsec = os.path.join(define['working_folder'], 'xsec_combined')

if options.title:
    output_folder = os.path.join(define['working_folder'], 'xsec_ktable_%s_%ig' % (options.title, ngauss))
else:
    output_folder = os.path.join(define['working_folder'], 'xsec_ktable_R%i_%s-%smu_%ig' % (options.resolution, lambdamin,
                                                                                            lambdamax,ngauss))


# create dir
os.system('mkdir -p %s' % output_folder)

def worker(press_idx):
    press_val = define['press_list'][press_idx]

    for temp_idx, temp_val in enumerate(define['temp_list']):

        #if options.title:
        ktable_file = '%s_T%i_P%.4e_%s.ktable.h5' % (options.key_iso_ll, temp_val, press_val,
                                                             options.title)
       # else:
       #     ktable_file = '%s_T%i_P%.4e_R%i_%s-%s.ktable.h5' % (options.key_iso_ll, temp_val, press_val,
       #                                                             options.resolution, lambdamin,
       #                                                             lambdamax)

        if os.path.isfile(os.path.join(output_folder, ktable_file)):
            print('Skip pressure ', press_val, ' Temperature ', temp_val)
            continue

        print('Doing pressure ', press_val, ' Temperature ', temp_val)


        xsec_file = '%s_T%i_P%.4e.xsec.pickle' % (mol_name, temp_val, press_val)

        xsec_in = pickle.load(open(os.path.join(folder_xsec, xsec_file), 'rb'))
        bingrid_idx = np.digitize(xsec_in[:,0], wn_grid) # get the indexes for bins

        ktable = np.zeros((len(wn_grid)-1, ngauss))
        for i in range(1, len(wn_grid)):

            x = xsec_in[:,0][bingrid_idx == i]
            y = xsec_in[:,1][bingrid_idx == i]

            if len(x) > 0:

                # reinterpolate to finer grid
                x_new = np.linspace(np.min(x), np.max(x), len(x)*10)
                y_new = np.interp(x_new, x, y)

                sort_bin = np.sort(y_new)
                norm_x = ((x_new-np.min(x_new))/(np.max(x_new) - np.min(x_new)))*2 - 1
                kcoeff = np.interp(gauss[0], norm_x, sort_bin)
                ktable[i-1,:]  = kcoeff

        # dumping ktable to file
        #pickle.dump(ktable, open(os.path.join(output_folder, ktable_file), 'wb'), protocol=2)
        hdf5_out_T = h5py.File(os.path.join(output_folder, ktable_file))
        hdf5_out_T['kcoeffs'] = ktable

pool = Pool(processes=int(ncores))
pool_result = pool.map(worker,[press_idx for press_idx, press_val in enumerate(define['press_list'])])
pool.close()
pool.join()

# aggreagate ktables to TauREx format

print('Aggregate ktables into a single file in TauREx format')

kcoeff = np.zeros((len(define['press_list']), len(define['temp_list']), len(bincentres), ngauss))

for temp_idx, temp_val in enumerate(np.sort(define['temp_list'])):
    for press_idx, press_val in enumerate(np.sort(define['press_list'])):

        #if options.title:
        #    ktable_file = '%s_T%i_P%.4e_%s.ktable.pickle' % (mol_name, temp_val, press_val,
        #                                                            options.title)
        #else:
        #    ktable_file = '%s_T%i_P%.4e_R%i_%s-%s.ktable.pickle' % (mol_name, temp_val, press_val,
        #                                                            options.resolution, lambdamin_str,
        #                                                            lambdamax_str)

        #if options.title:
        ktable_file = '%s_T%i_P%.4e_%s.ktable.h5' % (options.key_iso_ll, temp_val, press_val,
                                                             options.title)
        #else:
        #    ktable_file = '%s_T%i_P%.4e_R%i_%s-%s.ktable.h5' % (options.key_iso_ll, temp_val, press_val,
        #                                                            options.resolution, lambdamin,
        #                                                            lambdamax)
        ktab = h5py.File(os.path.join(output_folder, ktable_file), 'r')['kcoeffs']

        kcoeff[press_idx, temp_idx, :, :] = ktab[:,:]
        os.system('rm %s' % open(os.path.join(output_folder, ktable_file)))

kdist_out = {
    'bin_centers' : bincentres,
    'bin_edges' : wn_grid,
    'wlrange' : (lambdamin, lambdamax),
    'wnrange' : (10000./lambdamax, 10000./lambdamin),
    'weights' :  weights/2.,
    'samples' : (samples+1.)/2.,
    'ngauss' : ngauss,
    'method' : 'petitRADTRANS samples',
    'kcoeff' : kcoeff,
    'name' : mol_name,
    't' : np.sort(define['temp_list']).astype(float),
    'p' : np.sort(define['press_list']).astype(float),
}

if options.resolution:
    kdist_out['resolution'] = float(options.resolution)


hdf5_out = h5py.File(os.path.join(output_folder, '%s.R1000_0.3-50mu.ktable.petitRADTRANS.h5' %
                                         (options.key_iso_ll)),'w')

mass_mol  = options.mol_mass

str_type = h5py.new_vlen(str)
dset = hdf5_out.create_dataset("method",(1,), dtype=str_type)
dset2 = hdf5_out.create_dataset("mol_mass",(1,), dtype=int)
dset3 = hdf5_out.create_dataset("DOI",(1,), dtype=str_type)
dset4 = hdf5_out.create_dataset("Date_ID",(1,), dtype=str_type)
dset5 = hdf5_out.create_dataset("mol_name",(1,), dtype=str_type)

hdf5_out['bin_centers'] = bincentres
hdf5_out['bin_edges'] = wn_grid
hdf5_out['wlrange'] = (lambdamin, lambdamax)
hdf5_out['wnrange'] = (10000./lambdamax, 10000./lambdamin)
hdf5_out['weights'] = weights/2.
hdf5_out['samples'] = (samples+1.)/2.
#ngauss and method don't output
hdf5_out['ngauss'] = 16
#hdf5_out['method'] = 'polynomial.legendre.leggauss'
dset[0] = 'petit_samples'
hdf5_out['kcoeff'] = kcoeff
hdf5_out['kcoeff'].attrs['units'] = 'cm^2/molecule'
hdf5_out['t'] = np.sort(define['temp_list']).astype(float)
hdf5_out['p'] = np.sort(define['press_list']).astype(float)
hdf5_out['p'].attrs['units'] = 'bar'
dset2[0] = options.mol_mass
dset3[0] = 'rep4'
dset4[0] = 'rep5'
dset5[0] = options.mname

#str_type = h5py.new_vlen(str)
#dset = hdf5_out.create_dataset("name",(1,), dtype=str_type)
#dset2 = hdf5_out.create_dataset("key_iso_ll",(1,), dtype=str_type)
#dset3 = hdf5_out.create_dataset("resolution",(1,))
#dset4 = hdf5_out.create_dataset("lambdamin",(1,))
#dset5 = hdf5_out.create_dataset("lambdamax",(1,))
#dset[0] = mol_name
#dset2[0] = options.key_iso_ll
#dset3[0] = options.resolution
#dset4[0] = lambdamin
#dset5[0] = lambdamax





#res = options.resolution
#hdf5_out['resolution'] = res
hdf5_out.close()


    #pickle.dump(kdist_out, open(os.path.join(output_folder, '%s_R%i_%s-%smu.ktable.TauREx.pickle' %
    #                                         (mol_name, options.resolution,
    #                                          lambdamin_str, lambdamax_str)), 'wb'), protocol=2)

print('End')
