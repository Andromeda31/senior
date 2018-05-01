##lots of imports. Some are unnecessary but I left a lot just to be safe...
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
from pathlib import Path
import math
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.patches import Ellipse
import numpy.random as rnd
from matplotlib import patches
import sys
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable


drpall = t.Table.read('/home/celeste/Documents/astro_research/drpall-v2_3_1.fits')
drpall = fits.open('/home/celeste/Documents/astro_research/drpall-v2_3_1.fits')

print(drpall.info())

print(drpall[1].header)
print(drpall[2].data['nsa_sersic_flux'])


filename = '/home/celeste/Documents/astro_research/thesis_git/Good_Galaxies_SPX_3_N2S2.txt'
file_names = np.genfromtxt(filename, usecols=(0), skip_header=1, dtype=str, delimiter=',')

plate_num = []
fiber_num = []
split = []

header = drpall[1].header
data = drpall[2].data

print(header)

print(data)

mass_total = math.log10(drpall['nsa_elpetro_mass'])-np.log10(.49)
colorg = obj(['GFWHM'])
colorr = obj(['RFWHM'])
colori = obj(['IFWHM'])

#file_open = open("error_files.txt", "w")

mass_new = []
colorg_new = []
colorr_new = []
colori_new = []

##Goes through all files in the folder
for ii in range(0, len(file_names)):
    ##Removes all non alphanumeric characters and only leaves numbers and periods
    file_names[ii] = re.sub("[^0-9-]", "", file_names[ii])
    #print(file_names[ii])
    #print(file_names[ii][4:])
    #print(file_names[ii][:4])
    ##splits the two numbers into a plate number and fiber number
    one, two = (str(file_names[ii]).split('-'))
    ##splits the two numbers into a plate number and fiber number
    plate_num.insert(ii, one)
    fiber_num.insert(ii, two)
    
    
    
for i in range(0, len(plate_num)):
    plateifu = (str(plate_num[i]) + '-' + str(fiber_num[i]))
    obj = drpall[drpall['plateifu']==plateifu][0]
    mass_new.insert(i, math.log10(obj['nsa_elpetro_mass'])-np.log10(.49))
    colorg_new.insert(i, obj(['GFWHM']))
    colorr_new.insert(i, obj(['RFWHM']))
    colori_new.insert(i, obj(['IFWHM']))
    
    
plt.plot(mass_total, colorg)
plt.plot(mass_total, colorr)
plt.plot(mass_total, colori)
plt.plot(mass_new, colorg_new)
plt.plot(mass_new, colorr_new)
plt.plot(mass_new, colori_new)
plt.show()

