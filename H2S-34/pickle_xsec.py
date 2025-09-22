try:
    import cPickle as pickle
except:
    import pickle

import sys
import numpy as np

xs = np.loadtxt(sys.argv[1])
pickle.dump(xs, open('%s.pickle' % sys.argv[1], 'wb'), protocol=2)
