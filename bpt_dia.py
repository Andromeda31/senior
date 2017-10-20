##lots of imports. Some are unnecessary but I left a lot just to be safe...
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
import math
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

    ##stop():

    ##ON:
    ##7977-12705

    #Check later: 
    #7977-9102

    #Good examples of boring, normal galaxies:
    #7990-6101

    #Place-IFU (plate # - number of fibers, number bundle that is out
    #of the number of bundles with that number of fibers)
    ##Code from Zach, do not touch!!
def get_line_ew(maps,key,sn=3.):
    sew_hdu=maps['EMLINE_GFLUX']
    sew_ivar_hdu=maps['EMLINE_GFLUX_IVAR']

    # get a mapping from eline key to channel key
    v2k={v:k for(k,v)in sew_hdu.header.items()}
    # get a mapping from channel key to channel
    cstring2ix=lambda s:int(s[1:])-1

    ix=cstring2ix(v2k[key])

    ew=sew_hdu.data[ix,...]
    ew_ivar=sew_ivar_hdu.data[ix,...]
    snr= ew*np.sqrt(ew_ivar)
    ew_mask=(snr<sn)

    return np.ma.array(ew,mask=ew_mask)


    #plate_num = ['8454', '9041', '7990', '8259', '8318', '9026']
    #fiber_num = ['12703', '12701', '6104', '9101', '9101', '3703']

    ##plate_num = ['8455']
    ##fiber_num = ['3701']
    ##fiber_num = ['1901', '1902', '3701', '3702', '3703', '3704', '6101', '6102', '6103', '6104', '9101', '9102', '12701', '12702', '12704', '12705']

    ##After hoefully downloading all the required fits files, this will read all the names
file_names = os.listdir("/home/celeste/Documents/astro_research/keepers")

    ##creates the empty arrays to append the names of the files in the folder
plate_num = []
fiber_num = []
split = []

    ##Goes through all files in the folder
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

    ##Main loop over all the plates
for i in range(0, len(plate_num)): ##len(plate_num)
    ##for j in range(0, len(fiber_num)):
        print(plate_num[i] + '-' + fiber_num[i])
        print("Index: " + str(i))
        ##some black magic
        hdulist = fits.open('/home/celeste/Documents/astro_research/keepers/manga-' + plate_num[i] + '-' + fiber_num[i] + '-MAPS-SPX-GAU-MILESHC.fits.gz')

        ##gets the hydrogen alpha and hydrogen beta data
        hd2 = hdulist['SPX_ELLCOO'].data
        Ha = hdulist['EMLINE_GFLUX'].data[7,...] ##21 element array
        Hb = hdulist['EMLINE_GFLUX'].data[1,...]

        ##assigns the plate id based on what is in the data cube
        plate_id = hdulist['PRIMARY'].header['PLATEIFU']

        ##gets official plate number
        plate_number = hdulist['PRIMARY'].header['PLATEID']
        fiber_number = hdulist['PRIMARY'].header['IFUDSGN']

        #print(hdulist['PRIMARY'].header)

        #quit()

        ##Imports the thingy we need to get the images of the galaxy without having to download directly all the pictures. This also bypasses the password needed
		
	

        ##Also need OIII 5007 and NII 6584

        ##Gets the correct lines
        H_alpha = get_line_ew(hdulist, 'Ha-6564')
        H_beta = get_line_ew(hdulist, 'Hb-4862')


        OIII = get_line_ew(hdulist, 'OIII-5008')
        NII = get_line_ew(hdulist, 'NII-6585')

        ##line ratios
        #x: n2/Halpha (plot the log)
        O_B = OIII/H_beta

        N_A = NII/H_alpha
        
        zeros = False
        
        for element in range(0, len(O_B)):
            for item in range(0, len(O_B[element])):
                if O_B[element][item] >= 0:
                    zeros = True
                
        if zeros == True:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.set_xscale("log")
            #ax.set_yscale("log")
            minx=1000000
            maxx=-1000000
            miny=1000000
            maxy=-1000000
            ax.set_aspect(1)
            ax.set_title("BPT Diagram for " + str(plate_number) + "-" + str(fiber_number))
            for element in range(0, len(O_B)):
                for item in range(0, len(O_B[element])):
                    if (math.log10(O_B[element][item])<0.61/(math.log10(N_A[element][item])-0.47)+1.19) and (math.log10(O_B[element][item])<0.61/(math.log10(N_A[element][item])-0.05)+1.3) and (math.log10(N_A[element][item])<0.0):
                        #print("red dot")
                        ax.plot(math.log10(N_A[element][item]), math.log10(O_B[element][item]), color = "red", marker = ".", ls = "None")
                        """
                        if minx>(math.log10(N_A[element][item])):
                            minx=math.log10(N_A[element][item])
                        if maxx<(math.log10(N_A[element][item])):
                            maxx=math.log10(N_A[element][item])
                        if miny>(math.log10(O_B[element][item])):
                            miny=math.log10(O_B[element][item])
                        if maxy<(math.log10(O_B[element][item])):
                            maxy=math.log10(O_B[element][item])
                        """
                    else:
                        ax.plot(math.log10(N_A[element][item]), math.log10(O_B[element][item]), color = "gray", marker = ".", ls = "None")
                        """
                        if minx>(math.log10(N_A[element][item])):
                            minx=math.log10(N_A[element][item])
                        if maxx<(math.log10(N_A[element][item])):
                            maxx=math.log10(N_A[element][item])
                        if miny>(math.log10(O_B[element][item])):
                            miny=math.log10(O_B[element][item])
                        if maxy<(math.log10(O_B[element][item])):
                            maxy=math.log10(O_B[element][item])
                        """
            #Kewley
            X = np.linspace(-1.5, 0.3)
            Y = ((0.61/(X-0.47))+1.19)
            
            #Kauffmann
            Xk = np.linspace(-1.5,0.)
            Yk= (0.61/(Xk-0.05)+1.3)
            
            
            ax.plot(X, Y, '--', color = "red", lw = 1, label = "Kewley+01")
            ax.plot(Xk, Yk, '-', color = "blue", lw = 1, label = "Kauffmann+03")
            ax.set_xlim(-1.2, 1.2)
            ax.set_ylim(-1.5, 1)
            plt.draw()
            #plt.show()
            plt.savefig('/home/celeste/Documents/astro_research/astro_images/bpt_diagrams/bpt_' + str(plate_number) + "-" + str(fiber_number))
            plt.clf()
