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
#from marvin.tools.cube import Cube

drpall = t.Table.read('/home/celeste/Documents/astro_research/drpall-v2_3_1.fits')
obj = drpall['plateifu']

plate_num = []
fiber_num = []
split = []

'''

for ii in range(0, len(obj)):
    ##Removes all non alphanumeric characters and only leaves numbers and periods
    obj[ii] = re.sub("[^0-9]", "", obj[ii])
    #print(file_names[ii])
    #print(file_names[ii][4:])
    #print(file_names[ii][:4])
    ##splits the two numbers into a plate number and fiber number
    split.insert(1, obj[ii][4:])
    split.insert(0, obj[ii][:4])
    plate_num.insert(ii, split[0])
    fiber_num.insert(ii, split[1])
    ##gets rid of the old split so the new one can be inserted at the 0 and 1 position
    del split[1]
    del split[0]
    
'''

with open('/home/celeste/Documents/astro_research/thesis_git/adam_galaxies.txt') as f:
    file_names=[]
    for line in f:
        file_names.append(line)
    
for ii in range(0, len(file_names)):
    ##Removes all non alphanumeric characters and only leaves numbers and periods
    file_names[ii] = re.sub("[^0-9]", "", file_names[ii])
    #print(file_names[ii])
    #print(file_names[ii][4:])
    #print(file_names[ii][:4])
    ##splits the two numbers into a plate number and fiber number
    split.insert(1, file_names[ii][4:])
    split.insert(0, file_names[ii][:4])
    plate_num.insert(ii, split[0])
    fiber_num.insert(ii, split[1])
    ##gets rid of the old split so the new one can be inserted at the 0 and 1 position
    del split[1]
    del split[0]
    
print(len(plate_num))
print(len(fiber_num))

import os
import re
import requests

for x in range(0, len(plate_num)):
        

        r = requests.get('https://data.sdss.org/sas/mangawork/manga/spectro/redux/v2_3_1/' + str(plate_num[x]) + '/stack/images/' + str(fiber_num[x]) + '.png', auth=('sdss', '2.5-meters'))

        ##Saves the image
        with open('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + str(plate_num[x]) + '-' + str(fiber_num[x]) + '.png', 'wb') as fd:
	        for chunk in r.iter_content(chunk_size=128):
		        fd.write(chunk)
        fd.close()
		        
        print("Saved " + str(plate_num[x]) + "-" + str(fiber_num[x]))
