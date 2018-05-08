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
import re


drpall=fits.open('/home/celeste/Documents/astro_research/drpall-v2_3_1.fits')

mags=drpall[1].data['NSA_ELPETRO_ABSMAG']

u=mags[:,2]
g=mags[:,3]
r=mags[:,4]
i=mags[:,5]

mass_total = np.log10(drpall[1].data['nsa_elpetro_mass'])-np.log10(.49)

plt.scatter(mass_total, g-i, c = 'g', alpha = 0.5)
#plt.scatter(mass_total, r-i, c = 'r', alpha = 0.5)
#plt.scatter(mass_total, u-r, c = 'b', alpha = 0.5)
plt.xlabel("$Log_{10}$ Stellar Mass")
xmin, xmax = 8.5, 12
ymin, ymax = -1, 2
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)


filename = '/home/celeste/Documents/astro_research/thesis_git/Good_Galaxies_SPX_3_N2S2.txt'
file_names = np.genfromtxt(filename, usecols=(0), skip_header=1, dtype=str, delimiter=',')

plate_num = []
fiber_num = []
split = []


#file_open = open("error_files.txt", "w")

mass_new = []
colorg_new = []
colorr_new = []
colori_new = []
coloru_new = []

tbdata = drpall[1].data

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
    
    
  
for x in range(0, len(plate_num)):
    plateifu = (str(plate_num[x]) + '-' + str(fiber_num[x]))

    ind = np.where(tbdata['plateifu'] == plateifu)
    mass_new.insert(x, np.log10(tbdata['nsa_elpetro_mass'][ind][0])-np.log10(.49))
    coloru_new.insert(x, tbdata['NSA_ELPETRO_ABSMAG'][ind][:,2])
    colorg_new.insert(x, tbdata['NSA_ELPETRO_ABSMAG'][ind][:,3])
    colorr_new.insert(x, tbdata['NSA_ELPETRO_ABSMAG'][ind][:,4])
    colori_new.insert(x, tbdata['NSA_ELPETRO_ABSMAG'][ind][:,5])
    
print(mass_new)
print(np.asarray(coloru_new)-np.asarray(colorr_new))
plt.scatter(mass_new, np.asarray(colorg_new)-np.asarray(colori_new), c='black')
plt.ylabel("g-i color")
plt.savefig("g-i_color.png")
plt.close()


