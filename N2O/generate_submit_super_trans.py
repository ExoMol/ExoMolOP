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
#folder_tools = define['tools_folder']

## Make some folders: input, output, submit and xsec
os.system('mkdir -p ' + define['working_folder'])
folder_inp = os.path.join(define['working_folder'], 'input')
os.system('mkdir -p ' + folder_inp)
folder_out = os.path.join(define['working_folder'], 'output')
os.system('mkdir -p ' + folder_out)
folder_submit = os.path.join(define['working_folder'], 'submit')
os.system('mkdir -p ' + folder_submit)
folder_xsec = os.path.join(define['working_folder'], 'xsec')
os.system('mkdir -p ' + folder_xsec)
folder_done = os.path.join(define['working_folder'], 'done')
os.system('mkdir -p ' + folder_done)
folder_comb = os.path.join(define['working_folder'], 'xsec_combined')
os.system('mkdir -p ' + folder_comb)

# Submission script preamble for clusters
#if  define['platform'] == 'cobweb':
submit_pre = '#! /bin/bash \n'
submit_pre += '#PBS -S /bin/bash \n'
if 'queue' in define:
    submit_pre += '#PBS -q %s \n' % define['queue']
else:
    submit_pre += '#PBS -q compute \n'

submit_pre += '#PBS -V \n'
submit_pre += '#PBS -l nodes=1:ppn=%i \n' % nthreads
submit_pre += '#PBS -l walltime=100:00:00 \n\n'
submit_pre += 'module unload compiler/intel/l_ics_2013.1.046+mpi \n'
submit_pre += 'module load compiler/intel/15.1/133 \n\n'
#elif define['platform'] == 'legion':
#    submit_pre = '#!/bin/bash -l \n'
#    submit_pre += '#$ -S /bin/bash \n'
#    submit_pre += '#$ -l h_rt=12:0:0 \n'
#    submit_pre += '#$ -l mem=%iG \n' % np.int(define['memory'])
#    submit_pre += '#$ -pe smp %i \n' % nthreads
#    submit_pre += '#$ -wd %s \n' % os.path.join(define['working_folder'], 'submit')

# interpolate partition function
if 'partition_function' in define:
    p = np.loadtxt(os.path.join(define['linelist_folder'] , define['partition_function']))
    from scipy.interpolate import InterpolatedUnivariateSpline
    part_func = InterpolatedUnivariateSpline(p[:,0], p[:,1], k=2) # quadratic spline extrapolation for T outside range

# Generate input and submit files.

#submitsh1 = ''
#submitsh2 = ''

batch_submit_super1 = open(os.path.join(define['working_folder'], 'batch_submit_super1.sh'), "w")
batch_submit_super2 = open(os.path.join(define['working_folder'], 'batch_submit_super2.sh'), "w")

for rng_idx, rng_val in enumerate(define['ranges']):

    print('Range ', rng_val)

    #batch_submit = open(os.path.join(define['working_folder'], 'batch_submit_%s-%s.sh' % (str(rng_val[0]).zfill(5), str(rng_val[1]).zfill(5))), "w")


#superlines batch 1 - only for each temperature, not pressure
    for temp_part_idx in range(int(len(temp_list)/ngroup)):
        out_file1 = '%s_T%i_sup1' % (mol_name, temp_list[temp_part_idx*ngroup])

        if not os.path.exists(os.path.join(folder_done, '%s.done' % out_file1)):

   # wavenumber range
            #range_start = rng_val[0]
            #range_end = rng_val[1]
               # wavenumber range 
            sup1_rng_start = define['fullrange'][0][0]
            sup1_rng_end = define['fullrange'][0][1]

                    # temperatures are NOT grouped with exocross. One computation per temperature
            for i in range(ngroup):
                xsec1_file = '%s_T%i_sup1' % (mol_name, temp_list[temp_part_idx*ngroup+i])

                outinp1 = ''
                if 'dbtype' in define:
                    outinp1 += '%s \n' % define['dbtype']
                if 'iso' in define:
                    outinp1 += 'iso %i \n' % define['iso']
                outinp1 += 'absorption \n'
                outinp1 += 'bin \n'
                outinp1 += 'verbose 3 \n'
                if 'molname' in define:
                    outinp1 += 'molname %s \n' % define['molname']
                #outinp1 += 'offset %.1f \n' % cutoff
                outinp1 += 'mass %f \n' % define['mean-mass']
                outinp1 += 'temperature %f \n' %  temp_list[temp_part_idx*ngroup+i]
                #outinp1 += 'pressure %f \n' %  press_val
                outinp1 += 'range %f %f \n' % (sup1_rng_start, sup1_rng_end)
#always use R 1000000 for superlines1
                outinp1 += 'R 1000000 \n'
                if 'partition_function' in define:
                    outinp1 += 'pf %f \n' % part_func(temp_list[temp_part_idx * ngroup + i])
                if 'ptfile' in define:
                    outinp1 += 'pffile %s \n' %  os.path.join(define['linelist_folder'], define['ptfile'])

                outinp1 += 'output %s \n' % os.path.join(folder_xsec, xsec1_file)
                if 'state-file' in define:
                    outinp1 += 'states %s\n' % os.path.join(define['linelist_folder'], define['state-file'])
                outinp1 += 'transitions \n'
                for i in range(len(define['full-trans'][0])):
                    outinp1 += '  "%s" \n' %  os.path.join(define['linelist_folder'], define['full-trans'][0][i])
                outinp1 += 'end \n'



            text_file1 = open(os.path.join(folder_inp, '%s.inp' % out_file1), "w")
            text_file1.write(outinp1)
            text_file1.close()
            #submitsh1 = ''
            #submitsh1 += ' set OMP_NUM_THREADS=%i \n' % define['nthreads']
            #submitsh1 +=  'if [ ! -f %s ]; then \n' % os.path.join(folder_done, '%s.done ' % out_file1)
            #submitsh1 += '  %s < %s > %s \n' % (define['executable'], os.path.join(folder_inp, '%s.inp' % out_file1),
            #                                      os.path.join(folder_out, '%s.out' % out_file1))

            #submitsh1 +=  ' touch %s \n' % os.path.join(folder_done, '%s.done \n' % out_file1)
            #submitsh1 += 'else \n'
            #submitsh1 += '  echo "Skip file, already done. " \n'
            #submitsh1 += 'fi \n'


            #submit_file1 = open(os.path.join(folder_submit, '%s.sh' % out_file1), "w")
           # submit_file1.write(submit_pre)
          #  submit_file1.write(submitsh1)
         #   submit_file1.close()


            batch_submit_super1.write('./b_trove_sup1.sh %s.inp 24\n' % out_file1)

    batch_submit_super1.close()

#superlines batch 2 - for both pressure and temperature lists 
    for press_idx, press_val in enumerate(press_list):

        for temp_part_idx in range(int(len(temp_list)/ngroup)):

            out_file2 = '%s_T%i_P%.4e_super2' % (mol_name, temp_list[temp_part_idx*ngroup], press_val)

            if not os.path.exists(os.path.join(folder_done, '%s.done' % out_file2)):

                # determing average gamma_V for group of temperatures in this wavenumber range
                gamma_V = np.average([ get_gamma_V(temp_list[temp_part_idx*ngroup+i],
                                                   press_val,
                                                   np.linspace(rng_val[0], rng_val[1], 100),
                                                   n=define['gamma-n'], gamma=define['gamma']) for i in range(ngroup)])

                grid_spacing = gamma_V/nvsample # determine grid spacing
                if grid_spacing > min_gridspacing:
                    grid_spacing = min_gridspacing

                nsamples =  int(np.round((rng_val[1]-rng_val[0])/grid_spacing)) # total number of samples

                # wavenumber range
                range_start = rng_val[0]
                if rng_idx < (len(define['ranges'])-1):
                    range_end = rng_val[1]-grid_spacing # remove last point
                else:
                    range_end = rng_val[1]

                if 'fixed_cutoff' in define:
                    cutoff = define['fixed_cutoff']
                else:
                    # cutoff
                    cutoff = (gamma_V*500.)
                    if cutoff > max_cutoff:
                        cutoff = max_cutoff


                if define['software'] == 'exocross':

                    # temperatures are NOT grouped with exocross. One computation per temperature
                    for i in range(ngroup):
                        xsec2_file = '%s_T%i_P%.4e' % (mol_name, temp_list[temp_part_idx*ngroup+i],press_val)

                        outinp2 = ''
                        if 'dbtype' in define:
                            outinp2 += '%s \n' % define['dbtype']
                        if 'iso' in define:
                            outinp2 += 'iso %i \n' % define['iso']
                        outinp2 += 'super-lines \n'
                        outinp2 += 'absorption \n'
                        outinp2 += 'voigt \n'
                        outinp2 += 'verbose 3 \n'
                        if 'molname' in define:
                            outinp2 += 'molname %s \n' % define['molname']
                        outinp2 += 'offset %.1f \n' % cutoff
                        outinp2 += 'mass %f \n' % define['mean-mass']
                        outinp2 += 'temperature %f \n' %  temp_list[temp_part_idx*ngroup+i]
                        outinp2 += 'pressure %f \n' %  press_val
                        outinp2 += 'range %f %f \n' % (range_start, range_end)
                        outinp2 += 'npoints %i \n' % nsamples
                        if 'partition_function' in define:
                            outinp2 += 'pf %f \n' % part_func(temp_list[temp_part_idx * ngroup + i])
                        if 'ptfile' in define:
                            outinp2 += 'pffile %s \n' %  os.path.join(define['linelist_folder'], define['ptfile'])

                        outinp2 += 'output %s \n' % os.path.join(folder_xsec, xsec2_file)
                        if 'state-file' in define:
                            outinp2 += 'states %s\n' % os.path.join(define['linelist_folder'], define['state-file'])

                        xsec1_super = '%s_T%i_sup1.xsec' % (mol_name,temp_list[temp_part_idx*ngroup+i])
                        outinp2 += 'transitions %s\n' % os.path.join(folder_xsec, xsec1_super)



                        outinp2 += 'species \n'
                        outinp2 += '  0 gamma %f n %f t0 296 ratio 0.86  \n' % \
                                      ( define['gamma'], define['gamma-n'])
                        outinp2 += '  1 gamma %f n %f t0 296 ratio 0.14  \n' % \
                                      ( define['gamma-He'], define['gamma-n-He'])
                        outinp2 += 'end \n'

                       # if 'broadeners' in define:
                      #      outinp2 += 'species \n'
#                            for broad_idx in define['broadeners']:
#
#                                if isinstance(define['broadeners'][broad_idx][1], list):
#                                    outinp2 += '  %i gamma %f n %f t0 296 ratio %f  \n' %  \
#                                              (broad_idx,
#                                               define['broadeners'][broad_idx][1][0], # gamma
#                                               define['broadeners'][broad_idx][1][1], # n
#                                               define['broadeners'][broad_idx][0]) # ratio
#                                else:
#                                    outinp2 += '  %i gamma 0.05 n 0.07 t0 296 file %s model J ratio %f  \n' %  \
#                                              (broad_idx,
#                                               os.path.join(define['linelist_folder'], define['broadeners'][broad_idx][1]),
#                                               define['broadeners'][broad_idx][0])
##                            outinp2 += 'end \n'
#                        else:
#                            outinp2 += 'species \n'
#                            outinp2 += '  0 gamma %f n %f t0 296 ratio 1  \n' % \
#                                      ( define['gamma'], define['gamma-n'])
#
                text_file2 = open(os.path.join(folder_inp, '%s.inp' % out_file2), "w")
                text_file2.write(outinp2)
                text_file2.close()

                #submitsh2 = ''

                #if define['software'] == 'exocross':
                #    submitsh2 += ' set OMP_NUM_THREADS=%i \n' % define['nthreads']

                # run only if .done file is not found
                #submitsh2 +=  'if [ ! -f %s ]; then \n' % os.path.join(folder_done, '%s.done ' % out_file2)

                # run file
                #submitsh2 += '  %s < %s > %s \n' % (define['executable'], os.path.join(folder_inp, '%s.inp' % out_file2),
                #                                  os.path.join(folder_out, '%s.out' % out_file2))

                for i in range(ngroup):
                    xsec_file2 = '%s_T%i_P%.4e' % (mol_name,
                                              temp_list[temp_part_idx*ngroup+i],
                                              press_val)

                    # convert .xsec to pickle file (use python script located in [tools_folder]/pickle.xsec.py)

                #    submitsh2 +=  ' %s %s %s \n' % (define['python'],
                #                                   os.path.join(folder_tools, 'pickle_xsec.py'),
                #                                   os.path.join(folder_xsec, '%s.xsec' % xsec_file2))

                # create .done file
                #submitsh2 +=  ' touch %s \n' % os.path.join(folder_done, '%s.done \n' % out_file2)

                #remove output - maybe add back this line if very large
#                submitsh += '  rm %s \n' % os.path.join(folder_out, '%s.out' % out_file)

                #submitsh2 += 'else \n'
                #submitsh2 += '  echo "Skip file, already done. " \n'
                #submitsh2 += 'fi \n'


                #submit_file2 = open(os.path.join(folder_submit, '%s.sh' % out_file2), "w")
                #submit_file2.write(submit_pre)
                #submit_file2.write(submitsh2)
                #submit_file2.close()

                batch_submit_super2.write('./b_trove_sup2.sh %s.inp 4\n' % out_file2)

    batch_submit_super2.close()
