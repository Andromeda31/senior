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
import re
import csv

import os

import requests

filename = '/home/celeste/Documents/astro_research/thesis_git/Good_Galaxies_SPX_2_N2S2.txt'

file_names = np.genfromtxt(filename, usecols=(0), skip_header=1, dtype=str, delimiter=',')

plate_num = []
fiber_num = []
split = []


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
    
print(len(plate_num))
print(len(fiber_num))

    
for c in range(0, len(plate_num)):
    #for cc in range(0, len(plate_num))
        r = requests.get('https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-6/SPX-GAU-MILESHC/' + str(plate_num[c]) + '/' + str(fiber_num[c]) + '/manga-' + str(plate_num[c]) + '-' + str(fiber_num[c]) + '-MAPS-SPX-GAU-MILESHC.fits.gz', auth=('sdss', '2.5-meters'))

        ##Saves the file
        with open('/home/celeste/Documents/astro_research/downloaded_data/MPL-6/manga-' + str(plate_num[c]) + '-' + str(fiber_num[c]) + '-MAPS-SPX-GAU-MILESHC.fits.gz', 'wb') as fd:
	        for chunk in r.iter_content(chunk_size=128):
		        fd.write(chunk)
        print("Downloaded " + str(plate_num[c]) + '-' + str(fiber_num[c]))
        
        r = requests.get('https://data.sdss.org/sas/mangawork/manga/spectro/redux/v2_3_1/' + str(plate_num[c]) + '/stack/images/' + str(fiber_num[c]) + '.png', auth=('sdss', '2.5-meters'))

        ##Saves the image
        with open('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + str(plate_num[c]) + '-' + str(fiber_num[c]) + '.png', 'wb') as fd:
	        for chunk in r.iter_content(chunk_size=128):
		        fd.write(chunk)
        fd.close()
		        
        print("Saved picture " + str(plate_num[c]) + "-" + str(fiber_num[c]))
        print("----------------------------------")
