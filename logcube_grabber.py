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

plt.rcParams['axes.facecolor'] = 'white'

#fun galaxy: 8332-12701


#from marvin.tools.cube import Cube
'''
import urllib.request
import urllib.parse
import urllib.request
'''
import re
import csv

import os

import requests

import numpy as np
from scipy.stats import chi2
#import pylab as mp

filename = '/home/celeste/Documents/astro_research/thesis_git/Good_Galaxies_SPX_3_N2S2.txt'
file_names = np.genfromtxt(filename, usecols=(0), skip_header=1, dtype=str, delimiter=',')
"""
with open('/home/celeste/Documents/astro_research/thesis_git/Good_Galaxies_SPX_2_N2S2.txt') as f:
    file_names=[]
    for line in f:
        file_names.append(line)
"""
with open('/home/celeste/Documents/astro_research/thesis_git/mass_data.txt') as f:
    mass_data=[]
    for line in f:
        mass_data.append(line)
#print(file_names)
##creates the empty arrays to append the names of the files in the folder
plate_num = []
fiber_num = []
split = []

#file_open = open("error_files.txt", "w")

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
    
#print(plate_num[0] + "-" + fiber_num[0])
#print(file_names[0])


    ##Main loop over all the plates
    
"""
Bad Plots?

8445-3701
8332-1902
8309-3703

"""
    
#plate_num=['9183']
#fiber_num = ['9102']

import shutil
from urllib.request import urlretrieve

plate_num = ['8718']
fiber_num = ['3703']

for i in range(0, len(plate_num)): ##len(plate_num)
##for j in range(0, len(fiber_num)):
        print(plate_num[i] + '-' + fiber_num[i])
        print("Index: " + str(i))
        '''We gon get all the logcubes
        
        '''
        r = requests.get('https://data.sdss.org/sas/mangawork/manga/spectro/redux/v2_4_3/' + plate_num[i] + '/stack/manga-' + plate_num[i] + '-' + fiber_num[i] + '-LOGCUBE.fits.gz', auth=('sdss', '2.5-meters'), stream = True)

        ##Saves the file
        #handle = open('/home/celeste/Documents/astro_research/logcube_files/manga-'+ str(plate_num[i])+ '-' + str(fiber_num[i]) + '-LOGCUBE.fits.gz', 'wb')
        handle = open('/media/celeste/Hypatia/MPL7/LOGCUBES/manga-'+ str(plate_num[i])+ '-' + str(fiber_num[i]) + '-LOGCUBE.fits.gz', 'wb')
        for chunk in r.iter_content(chunk_size=1024*100):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)
            #shutil.copyfileobj(r.raw, f)

        print("Downloaded " + str(plate_num[i]) + '-' + str(fiber_num[i]))
        print("----------------------------------")
        

        if i == 20:
            blargh
