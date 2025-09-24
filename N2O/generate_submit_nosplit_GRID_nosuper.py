import numpy as np
try:
    import cPickle as pickle
except:
    import pickle
import os
import sys
import imp

# Define some variables
AMU = 1.660538921e-27
KBOLTZ = 1.380648813e-16
AU = 1.49e11
C = 29979245800.00
NA = 6.022140857e23
ATM = 0.986923

## Define functions
def get_gamma_V(T, P, nu, n=0.5, gamma=0.07):
    Tref = 296.25
    n = 0.5
    m =  18.01528 / NA
    P = P*ATM # in atmosphere
    gamma_G = np.sqrt(2*KBOLTZ*T/m) * (nu/C)
    gamma_L = (Tref/T)**n * P * gamma
    gamma_V =  0.05346 * (2*gamma_L) + np.sqrt(0.2166 * (2*gamma_L)**2 + gamma_G**2)
    return gamma_V

## Import input dictionary
try:
    mol = imp.load_source('molecule', sys.argv[1])
    define = mol.define
except:
    print(' Cannot import dictionary file')
    exit()

## Set some variable names
mol_name = define['molecule_name']
temp_list = define['temp_list']
press_list = define['press_list']
ngroup = define['ngroup']
nvsample = define['nvsample']
min_gridspacing = define['min_gridspacing']
max_cutoff = define['max_cutoff']
nthreads = define['nthreads']
folder_tools = define['tools_folder']

## Make some folders: input, output, submit and xsec
os.system('mkdir -p ' + define['working_folder'])
folder_inp = os.path.join(define['working_folder'], 'input')
os.system('mkdir -p ' + folder_inp)
folder_out = os.path.join(define['working_folder'], 'output')
os.system('mkdir -p ' + folder_out)
#folder_submit = os.path.join(define['working_folder'], 'submit')
#os.system('mkdir -p ' + folder_submit)
folder_xsec = os.path.join(define['working_folder'], 'xsec')
os.system('mkdir -p ' + folder_xsec)
folder_comb = os.path.join(define['working_folder'], 'xsec_combined')
os.system('mkdir -p ' + folder_comb)

# Submission script preamble for clusters
#if  define['platform'] == 'cobweb':
#    submit_pre = '#! /bin/bash \n'
#    submit_pre += '#PBS -S /bin/bash \n'
#    if 'queue' in define:
#        submit_pre += '#PBS -q %s \n' % define['queue']
#    else:
#        submit_pre += '#PBS -q compute \n'
#
#    submit_pre += '#PBS -V \n'
#    submit_pre += '#PBS -l nodes=1:ppn=%i \n' % nthreads
#    submit_pre += '#PBS -l walltime=100:00:00 \n\n'
#    submit_pre += 'module unload compiler/intel/l_ics_2013.1.046+mpi \n'
#    submit_pre += 'module load compiler/intel/15.1/133 \n\n'
#elif define['platform'] == 'legion':
#    submit_pre = '#!/bin/bash -l \n'
#    submit_pre += '#$ -S /bin/bash \n'
#    submit_pre += '#$ -l h_rt=12:0:0 \n'
#    submit_pre += '#$ -l mem=%iG \n' % np.int(define['memory'])
#    submit_pre += '#$ -pe smp %i \n' % nthreads
 #   submit_pre += '#$ -wd %s \n' % os.path.join(define['working_folder'], 'submit')

# interpolate partition function
if 'partition_function' in define:
    p = np.loadtxt(os.path.join(define['linelist_folder'] , define['partition_function']))
    from scipy.interpolate import InterpolatedUnivariateSpline
    part_func = InterpolatedUnivariateSpline(p[:,0], p[:,1], k=2) # quadratic spline extrapolation for T outside range

# Generate input and submit files.

submitsh = ''

batch_submit_all = open(os.path.join(define['working_folder'], 'batch_submit_all.sh'), "w")

#for rng_idx, rng_val in enumerate(define['ranges']):

 #   print('Range ', rng_val)

    #batch_submit = open(os.path.join(define['working_folder'], 'batch_submit_%s-%s.sh' % (str(rng_val[0]).zfill(5), str(rng_val[1]).zfill(5))), "w")

    # get the line lists for this range if we are on legion


for press_idx, press_val in enumerate(press_list):

    for temp_part_idx in range(int(len(temp_list)/ngroup)):

        out_file = '%s_T%i_P%.4e' % (mol_name,
                                                                          temp_list[temp_part_idx*ngroup],
                                                                              press_val)

        if not os.path.exists(os.path.join(folder_inp, '%s.done' % out_file)):


            if define['software'] == 'exocross':

                # temperatures are NOT grouped with exocross. One computation per temperature
                for i in range(ngroup):
                    xsec_file = '%s_T%i_P%.4e' % (mol_name,
                                              temp_list[temp_part_idx*ngroup+i],
                                              press_val)

                        # partition function?
                    outinp = ''
                    if 'dbtype' in define:
                    	outinp = '%s \n' % define['dbtype']
                    if 'iso' in define:
                    	outinp += 'iso %i \n' % define['iso']
                    outinp += 'absorption \n'
                    outinp += 'voigt \n'
                    outinp += 'verbose 3 \n'
                    outinp += 'mass %f \n' % define['mean-mass']
                    outinp += 'temperature %f \n' %  temp_list[temp_part_idx*ngroup+i]
                    outinp += 'pressure %f \n' %  press_val
                    outinp += 'pffile %s \n' %  os.path.join(define['linelist_folder'], define['ptfile'])
                    outinp += 'output %s \n' % os.path.join(folder_xsec, xsec_file)
                    outinp += 'states %s\n' % os.path.join(define['linelist_folder'], define['state-file'])


                    for rng_idx2, rng_val2 in enumerate(define['ranges']):

                       # determing average gamma_V for group of temperatures in this wavenumber range
                       gamma_V = np.average([ get_gamma_V(temp_list[temp_part_idx*ngroup+i],
                                               press_val,
                                               np.linspace(rng_val2[0], rng_val2[1], 10),
                                               n=define['gamma-n'], gamma=define['gamma']) for i in range(ngroup)])

                       grid_spacing = gamma_V/nvsample # determine grid spacing
                       if grid_spacing > min_gridspacing:
                         grid_spacing = min_gridspacing

                       nsamples =  int(np.round((rng_val2[1]-rng_val2[0])/grid_spacing)) # total number of samples

                       # wavenumber range
                       range_start = rng_val2[0]
                       if rng_idx2 < (len(define['ranges'])-1):
                           range_end = rng_val2[1]-grid_spacing # remove last point
                       else:
                           range_end = rng_val2[1]

                       if 'fixed_cutoff' in define:
                           cutoff = define['fixed_cutoff']
                       else:
                           # cutoff
                           cutoff = (gamma_V*500.)
                           if cutoff > max_cutoff:
                               cutoff = max_cutoff

                    outinp += 'Range  %s %s \n' % (rng_val2[0], rng_val2[1])
                    outinp += 'R 1000000 \n'
                    outinp += 'offset %.1f \n' % cutoff
                    outinp += 'transitions \n'
                    outinp += '  "%s" \n' %  os.path.join(define['linelist_folder'], define['trans-files'][0][0])
                    outinp += 'end \n'

                    outinp += 'species \n'
                    outinp += '  0 gamma %f n %f t0 296 ratio 0.86  \n' % \
                                  ( define['gamma'], define['gamma-n'])
                    outinp += '  1 gamma %f n %f t0 296 ratio 0.14  \n' % \
                                  ( define['gamma-He'], define['gamma-n-He'])
                    outinp += 'end \n'

            else:

                outinp = 'EXOMOL \n'
                outinp += 'MEAN-MASS %f \n' % define['mean-mass']

            text_file = open(os.path.join(folder_inp, '%s.inp' % out_file), "w")
            text_file.write(outinp)
            text_file.close()

            # run only if .done file is not found
            submitsh =  'if [ ! -f %s ]; then \n' % os.path.join(folder_inp, '%s.done ' % out_file)

            # run file
            submitsh += '  %s < %s > %s \n' % (define['executable'], os.path.join(folder_inp, '%s.inp' % out_file),
                                                  os.path.join(folder_out, '%s.out' % out_file))

            if define['software'] != 'exocross':
                for i in range(ngroup):
                    # cross sections for different temperatures are all in the .out file. cat them and store
                    # them in seprate .xsec files. Only for CEXSY
                    xsec_file = '%s_T%i_P%.4e' % (mol_name,
                                              temp_list[temp_part_idx*ngroup+i],
                                              press_val)
                    submitsh +=  ' cat %s | grep "\[%i\]" | sed \'s/\[%i\]//g\' > %s \n' % \
                                 (os.path.join(folder_out, '%s.out' % out_file),
                                  i, i, os.path.join(folder_xsec, '%s.xsec' % xsec_file))

            for i in range(ngroup):
                xsec_file = '%s_T%i_P%.4e' % (mol_name,
                                          temp_list[temp_part_idx*ngroup+i],
                                          press_val)

                # convert .xsec to pickle file (use python script located in [tools_folder]/pickle.xsec.py)
                submitsh +=  ' %s %s %s \n' % (define['python'],
                                               os.path.join(folder_tools, 'pickle_xsec.py'),
                                               os.path.join(folder_xsec, '%s.xsec' % xsec_file))

                # transfer files to cobweb
                if define['platform'] == 'legion':
                    submitsh +=  ' scp -r %s %s@cobweb-login.star.ucl.ac.uk:%s \n' % (os.path.join(folder_xsec, '%s.*' % xsec_file),
                                                                                   define['cobweb_username'],
                                                                                   os.path.join(define['cobweb_working_folder'], 'xsec_in'))

            # create .done file
            submitsh +=  ' touch %s \n' % os.path.join(folder_inp, '%s.done \n' % out_file)

            #remove output
            submitsh += '  rm %s \n' % os.path.join(folder_out, '%s.out' % out_file)

            # transfer files to cobweb
            if define['platform'] == 'legion':
                submitsh +=  ' scp -r %s %s@cobweb-login.star.ucl.ac.uk:%s \n' % (os.path.join(folder_inp, '%s.done' % out_file),
                                                                               define['cobweb_username'],
                                                                               os.path.join(define['cobweb_working_folder'], 'input_in'))

            submitsh += 'else \n'
            submitsh += '  echo "Skip file, already done. " \n'
            submitsh += 'fi \n'

                #submit_file = open(os.path.join(folder_submit, '%s.sh' % out_file), "w")
                #submit_file.write(submit_pre)
                #submit_file.write(submitsh)
                #submit_file.close()

              #  if define['platform'] == 'emerald':
             #      batch_submit.write('bsub < %s \n' % os.path.join(folder_submit, '%s.sh' % out_file))
            #       batch_submit_all.write('bsub < %s \n' % os.path.join(folder_submit, '%s.sh' % out_file))
           #    if define['platform'] == 'cobweb' or define['platform'] == 'legion':
                #batch_submit.write('qsub %s \n' % os.path.join(folder_submit, '%s.sh' % out_file))
                #batch_submit_all.write('qsub %s \n' % os.path.join(folder_submit, '%s.sh' % out_file))
                #batch_submit.write('./b_trove_sup2.sh %s 1\n' % out_file)
        batch_submit_all.write('./b_trove_split_grid.sh %s 4\n' % out_file)

    #batch_submit.close()
batch_submit_all.close()

