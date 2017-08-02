'''
Created on Jul 30, 2017
@author: lixi729
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

Test if the user has the required Python environment for Xanthos v1.0
'''


import sys
print "# Python version: " + sys.version

# List all the locally installed Python packages
import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
for package in installed_packages_list:
    print package
print "\nRequired Packages for Xanthos v1.0:"

#numpy
try:
    import numpy as np
    print "# numpy version", np.__version__
except:
    print "!!! numpy is not installed."
    
#scipy
try:
    import scipy
    from scipy import io as spio
    print "# scipy version", scipy.__version__
except:
    print "!!! scipy is not installed."
    
# Matplotlib
try:
    import matplotlib
    print "# matplotlib version", matplotlib.__version__
except:
    print "!!! matplotlib is not installed."

# Pandas
try:
    import pandas
    print "# pandas version", pandas.__version__
except:
    print "!!! pandas is not installed."
    
# configobj
try:
    import configobj
    print "# configobj version", configobj.__version__
except:
    print "!!! configobj is not installed."
    
