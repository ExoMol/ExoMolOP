######## Definition dictionary for CH4 -- Optimised for cobweb

define = {}

# using exocross
define['software'] = 'exocross'
define['molecule_name'] = 'N2O'
# define['dbtype'] = 'exomol' # hitran format
# define['iso'] = 261   # molecule n 26, isotope 1

# for HITEMP xsec use custom exocross, and specify 'molname'. E.g. for CO2
#define['molname'] = 'N2O'

define['platform'] = 'cobweb'
define['queue'] = 'gpu'

# needed if you work on legion:
# define['legion_username'] = 'zcapfa9'
#define['cobweb_username'] = 'kchubb'
#define['cobweb_working_folder'] = '/scratch/dp060/dc-chub1/ExoCross/superlines_LiOH/'
#define['cobweb_linelists_folder'] = '/scratch/dp060/dc-chub1/ExoCross/linelists/LiOH/'

# if define['platform'] == 'legion':
#     define['python'] = '/shared/ucl/apps/python/bundles/python2/venv/bin/python'
#     define['executable'] = '/home/%s/Scratch/exocross/j-xsec_1206_C.x' % define['legion_username']
#     define['working_folder'] = '/home/%s/Scratch/CO2/' % define['legion_username']
#     define['linelist_folder'] = '/home/%s/Scratch/linelists/CO2/' % define['legion_username']
#     define['tools_folder'] = '/home/%s/Scratch/CEXSY/tools/' % define['legion_username']

if define['platform'] == 'cobweb':
    define['python'] = '/cm/shared/apps/python/intelpython3/bin/python'
    define['executable'] = '/scratch/dp060/dc-chub1/ExoCross/j-xcross_1606.x'
    define['working_folder'] = 'N2O'
    define['linelist_folder'] = 'linelists/'
    define['tools_folder'] = '/scratch/dp060/dc-chub1/ExoCross/tools/'

else:
    print('Only cobweb available')
    exit()

# Download line list on the fly from cobweb (only valid for legion!)
define['copy_linelists'] = False


# You can leave these as they are
define['ngroup'] = 1
define['nvsample'] = 4
define['min_gridspacing'] = 0.01
define['max_cutoff'] = 25
define['fixed_cutoff'] = 25

if define['platform'] == 'legion':
    define['nthreads'] = 4
    define['memory'] = 2

# on cobweb 1 core per xsec, memory 5 GB
elif define['platform'] == 'cobweb':
    define['nthreads'] = 1
    define['memory'] = 5

# For the standard TauREx T,P list use these:

# temperatures in kelvin
define['temp_list'] = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600,
                        1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400]

# pressures in bar
define['press_list'] = [  1.00000000e-05,   2.15443469e-05,   4.64158883e-05,   1.00000000e-04,
                           2.15443469e-04,   4.64158883e-04,   1.00000000e-03,   2.15443469e-03,
                           4.64158883e-03,   1.00000000e-02,   2.15443469e-02,   4.64158883e-02,
                           1.00000000e-01,   2.15443469e-01,   4.64158883e-01,   1.00000000e+00,
                           2.15443469e+00,   4.64158883e+00,   1.00000000e+01,   2.15443469e+01,
                           4.64158883e+01,   1.00000000e+02]

define['mean-mass'] = 44.0


# Recent version of Exocross supports partition file input (ptfile input). Use ptfile for HITRAN / HITEMP
define['ptfile'] = '14N2-16O__TYM.pf'

# Old version of exocross didn't support partition function files, you had to specify the exact value of pf
# for each run. Use this option if you are using old v. of exocross, but you want to use external pf file.
# The Script will interpolate and provide exocross with the right pf value
#
# *** You might want to use old v. of exocross for HITEMP linelists data which have weird format and are not supported
# by new v. of exocross
#
# define['partition_function'] = 'partition_func.dat'

define['gamma'] = 0.06
define['gamma-n'] = 0.3
define['gamma-He'] = 0.0301
define['gamma-n-He'] = 0.5


define['ranges'] = [(100, 20000)] # 0

define['fullrange'] = [(100,20000)] # 0


define['full-trans'] = [(
'14N2-16O__TYM__00000-01000.trans',
'14N2-16O__TYM__01000-02000.trans',
'14N2-16O__TYM__02000-03000.trans',
'14N2-16O__TYM__03000-04000.trans',
'14N2-16O__TYM__04000-05000.trans',
'14N2-16O__TYM__05000-06000.trans',
'14N2-16O__TYM__06000-07000.trans',
'14N2-16O__TYM__07000-08000.trans',
'14N2-16O__TYM__08000-09000.trans',
'14N2-16O__TYM__09000-10000.trans',
'14N2-16O__TYM__10000-11000.trans',
'14N2-16O__TYM__11000-12000.trans',
'14N2-16O__TYM__12000-13000.trans',
'14N2-16O__TYM__13000-14000.trans',
'14N2-16O__TYM__14000-15000.trans',
'14N2-16O__TYM__15000-16000.trans',
'14N2-16O__TYM__16000-17000.trans',
'14N2-16O__TYM__17000-18000.trans',
'14N2-16O__TYM__18000-19000.trans',
'14N2-16O__TYM__19000-20000.trans',
'14N2-16O__TYM__20000-21000.trans'
)]


define['trans-files'] = {}
define['trans-files'][0] = [
'14N2-16O__TYM__00000-01000.trans',
'14N2-16O__TYM__01000-02000.trans',
'14N2-16O__TYM__02000-03000.trans',
'14N2-16O__TYM__03000-04000.trans',
'14N2-16O__TYM__04000-05000.trans',
'14N2-16O__TYM__05000-06000.trans',
'14N2-16O__TYM__06000-07000.trans',
'14N2-16O__TYM__07000-08000.trans',
'14N2-16O__TYM__08000-09000.trans',
'14N2-16O__TYM__09000-10000.trans',
'14N2-16O__TYM__10000-11000.trans',
'14N2-16O__TYM__11000-12000.trans',
'14N2-16O__TYM__12000-13000.trans',
'14N2-16O__TYM__13000-14000.trans',
'14N2-16O__TYM__14000-15000.trans',
'14N2-16O__TYM__15000-16000.trans',
'14N2-16O__TYM__16000-17000.trans',
'14N2-16O__TYM__17000-18000.trans',
'14N2-16O__TYM__18000-19000.trans',
'14N2-16O__TYM__19000-20000.trans',
'14N2-16O__TYM__20000-21000.trans',
]



define['state-file'] = '14N2-16O__TYM.states'

#define['broadeners'] = {}
#define['broadeners'][0] = (0.86, 'LiOH_H2.broad')
#define['broadeners'][1] = (0.14, 'LiOH_He.broad')

