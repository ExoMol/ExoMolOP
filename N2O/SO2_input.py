######## Definition dictionary for SO2 -- Optimised for cobweb

define = {}

# using exocross
define['software'] = 'exocross'
define['molecule_name'] = 'SO2'
# define['dbtype'] = 'exomol' # hitran format
# define['iso'] = 261   # molecule n 26, isotope 1

# for HITEMP xsec use custom exocross, and specify 'molname'. E.g. for CO2
#define['molname'] = 'CO2'

define['platform'] = 'cobweb'
define['queue'] = 'gpu'

# needed if you work on legion:
# define['legion_username'] = 'zcapfa9'
#define['cobweb_username'] = 'kchubb'
#define['cobweb_working_folder'] = '/scratch/dp060/dc-chub1/ExoCross/superlines_SO2/'
#define['cobweb_linelists_folder'] = '/scratch/dp060/dc-chub1/ExoCross/linelists/SO2/'

# if define['platform'] == 'legion':
#     define['python'] = '/shared/ucl/apps/python/bundles/python2/venv/bin/python'
#     define['executable'] = '/home/%s/Scratch/exocross/j-xsec_1206_C.x' % define['legion_username']
#     define['working_folder'] = '/home/%s/Scratch/CO2/' % define['legion_username']
#     define['linelist_folder'] = '/home/%s/Scratch/linelists/CO2/' % define['legion_username']
#     define['tools_folder'] = '/home/%s/Scratch/CEXSY/tools/' % define['legion_username']

if define['platform'] == 'cobweb':
    define['python'] = '/cm/shared/apps/python/intelpython3/bin/python'
    define['executable'] = '/scratch/dp060/dc-chub1/ExoCross/j-xsec_2211_i17.x'
    define['working_folder'] = '/scratch/dp060/dc-chub1/ExoCross/files_ExoMolOP/'
    define['linelist_folder'] = '/scratch/dp060/dc-chub1/ExoCross/linelists/SO2/'
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
                           4.64158883e+01,   1.00000000e+02,]

define['mean-mass'] = 63.8664


# Recent version of Exocross supports partition file input (ptfile input). Use ptfile for HITRAN / HITEMP
define['ptfile'] = 'SO2.pf'

# Old version of exocross didn't support partition function files, you had to specify the exact value of pf
# for each run. Use this option if you are using old v. of exocross, but you want to use external pf file.
# The Script will interpolate and provide exocross with the right pf value
#
# *** You might want to use old v. of exocross for HITEMP linelists data which have weird format and are not supported
# by new v. of exocross
#
# define['partition_function'] = 'partition_func.dat'

define['gamma'] = 0.1408
define['gamma-n'] = 0.75
define['gamma-He'] = 0.07
define['gamma-n-He'] = 0.64

define['ranges'] = [(100, 8000)] # 0

define['fullrange'] = [(100,8000)] # 0
define['full-trans'] = [(
'SO2_00000-00100.trans',
'SO2_00100-00200.trans',
'SO2_00200-00300.trans',
'SO2_00300-00400.trans',
'SO2_00400-00500.trans',
'SO2_00500-00600.trans',
'SO2_00600-00700.trans',
'SO2_00700-00800.trans',
'SO2_00800-00900.trans',
'SO2_00900-01000.trans',
'SO2_01000-01100.trans',
'SO2_01100-01200.trans',
'SO2_01200-01300.trans',
'SO2_01300-01400.trans',
'SO2_01400-01500.trans',
'SO2_01500-01600.trans',
'SO2_01600-01700.trans',
'SO2_01700-01800.trans',
'SO2_01800-01900.trans',
'SO2_01900-02000.trans',
'SO2_02000-02100.trans',
'SO2_02100-02200.trans',
'SO2_02200-02300.trans',
'SO2_02300-02400.trans',
'SO2_02400-02500.trans',
'SO2_02500-02600.trans',
'SO2_02600-02700.trans',
'SO2_02700-02800.trans',
'SO2_02800-02900.trans',
'SO2_02900-03000.trans',
'SO2_03000-03100.trans',
'SO2_03100-03200.trans',
'SO2_03200-03300.trans',
'SO2_03300-03400.trans',
'SO2_03400-03500.trans',
'SO2_03500-03600.trans',
'SO2_03600-03700.trans',
'SO2_03700-03800.trans',
'SO2_03800-03900.trans',
'SO2_03900-04000.trans',
'SO2_04000-04100.trans',
'SO2_04100-04200.trans',
'SO2_04200-04300.trans',
'SO2_04300-04400.trans',
'SO2_04400-04500.trans',
'SO2_04500-04600.trans',
'SO2_04600-04700.trans',
'SO2_04700-04800.trans',
'SO2_04800-04900.trans',
'SO2_04900-05000.trans',
'SO2_05000-05100.trans',
'SO2_05100-05200.trans',
'SO2_05200-05300.trans',
'SO2_05300-05400.trans',
'SO2_05400-05500.trans',
'SO2_05500-05600.trans',
'SO2_05600-05700.trans',
'SO2_05700-05800.trans',
'SO2_05800-05900.trans',
'SO2_05900-06000.trans',
'SO2_06000-06100.trans',
'SO2_06100-06200.trans',
'SO2_06200-06300.trans',
'SO2_06300-06400.trans',
'SO2_06400-06500.trans',
'SO2_06500-06600.trans',
'SO2_06600-06700.trans',
'SO2_06700-06800.trans',
'SO2_06800-06900.trans',
'SO2_06900-07000.trans',
'SO2_07000-07100.trans',
'SO2_07100-07200.trans',
'SO2_07200-07300.trans',
'SO2_07300-07400.trans',
'SO2_07400-07500.trans',
'SO2_07500-07600.trans',
'SO2_07600-07700.trans',
'SO2_07700-07800.trans',
'SO2_07800-07900.trans',
'SO2_07900-08000.trans',
)]

define['trans-files'] = {}
define['trans-files'][0] = [
'SO2_00000-00100.trans',
'SO2_00100-00200.trans',
'SO2_00200-00300.trans',
'SO2_00300-00400.trans',
'SO2_00400-00500.trans',
'SO2_00500-00600.trans',
'SO2_00600-00700.trans',
'SO2_00700-00800.trans',
'SO2_00800-00900.trans',
'SO2_00900-01000.trans',
'SO2_01000-01100.trans',
'SO2_01100-01200.trans',
'SO2_01200-01300.trans',
'SO2_01300-01400.trans',
'SO2_01400-01500.trans',
'SO2_01500-01600.trans',
'SO2_01600-01700.trans',
'SO2_01700-01800.trans',
'SO2_01800-01900.trans',
'SO2_01900-02000.trans',
'SO2_02000-02100.trans',
'SO2_02100-02200.trans',
'SO2_02200-02300.trans',
'SO2_02300-02400.trans',
'SO2_02400-02500.trans',
'SO2_02500-02600.trans',
'SO2_02600-02700.trans',
'SO2_02700-02800.trans',
'SO2_02800-02900.trans',
'SO2_02900-03000.trans',
'SO2_03000-03100.trans',
'SO2_03100-03200.trans',
'SO2_03200-03300.trans',
'SO2_03300-03400.trans',
'SO2_03400-03500.trans',
'SO2_03500-03600.trans',
'SO2_03600-03700.trans',
'SO2_03700-03800.trans',
'SO2_03800-03900.trans',
'SO2_03900-04000.trans',
'SO2_04000-04100.trans',
'SO2_04100-04200.trans',
'SO2_04200-04300.trans',
'SO2_04300-04400.trans',
'SO2_04400-04500.trans',
'SO2_04500-04600.trans',
'SO2_04600-04700.trans',
'SO2_04700-04800.trans',
'SO2_04800-04900.trans',
'SO2_04900-05000.trans',
'SO2_05000-05100.trans',
'SO2_05100-05200.trans',
'SO2_05200-05300.trans',
'SO2_05300-05400.trans',
'SO2_05400-05500.trans',
'SO2_05500-05600.trans',
'SO2_05600-05700.trans',
'SO2_05700-05800.trans',
'SO2_05800-05900.trans',
'SO2_05900-06000.trans',
'SO2_06000-06100.trans',
'SO2_06100-06200.trans',
'SO2_06200-06300.trans',
'SO2_06300-06400.trans',
'SO2_06400-06500.trans',
'SO2_06500-06600.trans',
'SO2_06600-06700.trans',
'SO2_06700-06800.trans',
'SO2_06800-06900.trans',
'SO2_06900-07000.trans',
'SO2_07000-07100.trans',
'SO2_07100-07200.trans',
'SO2_07200-07300.trans',
'SO2_07300-07400.trans',
'SO2_07400-07500.trans',
'SO2_07500-07600.trans',
'SO2_07600-07700.trans',
'SO2_07700-07800.trans',
'SO2_07800-07900.trans',
'SO2_07900-08000.trans'
]

define['state-file'] = 'SO2.states'

#define['broadeners'] = {}
#define['broadeners'][0] = (0.86, 'SO2_H2.broad')
#define['broadeners'][1] = (0.14, 'SO2_He.broad')

