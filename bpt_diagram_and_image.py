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


#from https://casper.berkeley.edu/astrobaki/index.php/Plotting_Ellipses_in_Python

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
    

def pp04_o3n2_w_errs(n2,ha,o3,hb,n2_err,ha_err,o3_err,hb_err):
    o3n2=(o3/hb) / (n2/ha)
    lo3n2=np.log10(o3n2)
    o3hb_err=(o3/hb)*np.sqrt((o3_err/o3)**2 + (hb_err/hb)**2)
    n2ha_err=(n2/ha)*np.sqrt((n2_err/n2)**2 + (ha_err/ha)**2)
    o3n2_err=o3n2*np.sqrt((o3hb_err/(o3/hb))**2 + ((n2ha_err)/(n2/ha))**2)
    lo3n2_err=o3n2_err/(o3n2*np.log(10))
    met=8.73-0.32*lo3n2
    met_err=0.32*lo3n2_err
    return met,met_err
    
def is_sf_array(n2,ha,o3,hb):
    '''
    Checks whether arrays of line fluxes come from star formation based on BPT diagnostics
    returns 1 if spaxel is star-forming, 0 if non-star forming and nan if not determinable.
    '''
    issf=np.zeros(n2.shape)
    x=np.log10(n2/ha)
    y=np.log10(o3/hb)
    issf[np.where(x>0)]=0
    goodpix=np.where((y<(0.61/(x-0.47))+1.19) & (y<(0.61/(x-0.05))+1.3) & (x<0.0))
    badpix=np.where((np.isnan(x)==True) | (np.isnan(y)==True) | (np.isinf(x)==True) | (np.isinf(y)==True))
    issf[badpix]=np.nan
    issf[goodpix]=1
    return issf
    
def n2s2_dopita16_w_errs(ha,n22,s21,s22,ha_err,n22_err,s21_err,s22_err):
    '''
    N2S2 metallicity diagnostic from Dopita et al. (2016)
    includes a calculation of the errors
    '''
    y=np.log10(n22/(s21+s22))+0.264*np.log10(n22/ha)
    s2=s21+s22
    s2_err=np.sqrt(s21_err**2 + s22_err**2)
    met_err=(1.0/np.log(10)) * np.sqrt( (1+0.264**2)*(n22_err/n22)**2 + (s2_err/s2)**2 + (ha_err/ha)**2 )
    met=8.77+y
    return met, met_err



##Takes an input array and returns the values if a condition is met. Basically a glorified call to numpy.where

def array_if(array,condition=None):
    array1=np.array(array)
    if condition is None:
        return array1
    ret=array1[np.where(condition)]
    return ret
    
#It will take an image and radius array (distarr; same size as the image) as arguments. Because of the way it was coded up, it works best when the distarr array is in units of pixels, i.e. don't use the R_re array.


def radial_profile(image,centre=None,distarr=None,mask=None,binwidth=2,radtype='unweighted'):
    '''
    image=2D array to calculate RP of.
    centre=centre of image in pixel coordinates. Not needed if distarr is given.
    distarr=2D array giving distance of each pixel from the centre.
    mask = 2D array, 1 if you want to include given pixels, 0 if not.
    binwidth= width of radial bins in pixels.
    radtype='weighted' or 'unweighted'. Weighted will give you the average radius of pixels in the bin. Unweighted will give you the middle of each radial bin.
    '''
    distarr=distarr/binwidth
    if centre is None:
        centre=np.array(image.shape,dtype=float)/2
    if distarr is None:
        y,x=np.indices(image.shape)
        distarr=np.sqrt((x-centre[0])**2 + (y-centre[1])**2)
    if mask is None:
        mask=np.ones(image.shape)
    rmax=int(np.max(distarr))
    r_edge=np.linspace(0,rmax,rmax+1)
    rp=np.zeros(len(r_edge)-1)*np.nan
    nums=np.zeros(len(r_edge)-1)*np.nan
    sig=np.zeros(len(r_edge)-1)*np.nan
    r=np.zeros(len(r_edge)-1)*np.nan
    for i in range(0,len(r)):
        rp[i]=np.nanmean(image[np.where((distarr>=r_edge[i]) & (distarr<r_edge[i+1]) & (mask==1.0) & (np.isinf(image)==False))])
        nums[i]=len(np.where((distarr>=r_edge[i]) & (distarr<r_edge[i+1]) & (mask==1.0) & (np.isinf(image)==False) & (np.isnan(image)==False))[0])
        sig[i]=np.nanstd((image[np.where((distarr>=r_edge[i]) & (distarr<r_edge[i+1]) & (mask==1.0) & (np.isinf(image)==False))]))
        if radtype=='unweighted':
            r[i]=(r_edge[i]+r_edge[i+1])/2.0
        elif radtype=='weighted':
            r[i]=np.nanmean(distarr[np.where((distarr>=r_edge[i]) & (distarr<r_edge[i+1]) & (mask==1.0) & (np.isinf(image)==False))])
    r=r*binwidth
    return r,rp,nums,sig


"""
scatter_if:
Creates a scatter plot from two arrays, but only plots points if a condition is met. All other plt.scatter functionality should be retained.
"""

def scatter_if(x_in,y_in,condition=None,**kwargs):
    x_ret=array_if(x_in,condition)
    y_ret=array_if(y_in,condition)
    if 'c' in kwargs and type(kwargs['c'])==np.ndarray:
        kwargs['c']=array_if(kwargs['c'],condition)
    ax=plt.scatter(x_ret,y_ret,**kwargs)
    return ax


#plate_num = ['8454', '9041', '7990', '8259', '8318', '9026']
#fiber_num = ['12703', '12701', '6104', '9101', '9101', '3703']

##plate_num = ['8455']
##fiber_num = ['3701']
##fiber_num = ['1901', '1902', '3701', '3702', '3703', '3704', '6101', '6102', '6103', '6104', '9101', '9102', '12701', '12702', '12704', '12705']

##After hoefully downloading all the required fits files, this will read all the names
#file_names = os.listdir("/home/celeste/Documents/astro_research/keepers")
#file_names = os.listdir("/home/celeste/Documents/astro_research/fits_files")
filename = '/home/celeste/Documents/astro_research/thesis_git/Good_Galaxies_SPX_3_N2S2.txt'
file_names = np.genfromtxt(filename, usecols=(0), skip_header=1, dtype=str, delimiter=',')

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
    
#plate_num=['9194']
#fiber_num = ['12701']
#plate_num = ['8252', '8338', '8568', '9865']
#fiber_num = ['12705', '12701', '12702', '12705']

count_continue1 = 0
count_continue2 = 0
count_continue3 = 0
failed_maps = "failed maps\n"
failed_logcube = "failed_logcube\n"
failed_other = "failed_TYPERROR\n"

for i in range(0, len(plate_num)): ##len(plate_num)
##for j in range(0, len(fiber_num)):
        print(plate_num[i] + '-' + fiber_num[i])
        print("Index: " + str(i))

        
        ##some black magic
        #hdulist = fits.open('/home/celeste/Documents/astro_research/keepers/manga-' + plate_num[i] + '-' + fiber_num[i] + '-MAPS-SPX-GAU-MILESHC.fits.gz')
        #hdulist = fits.open('/home/celeste/Documents/astro_research/downloaded_data/MPL-7/manga-' + plate_num[i] + '-' + fiber_num[i] + '-HYB-SPX-GAU-MILESHC.fits.gz')

        try:
            hdulist = fits.open('/media/celeste/Hypatia/MPL7/HYB/allmaps/manga-' + plate_num[i] + '-' + fiber_num[i] + '-MAPS-HYB10-GAU-MILESHC.fits.gz')
        except FileNotFoundError:
            failed_maps = failed_maps + str(plate_num[i]) + "-" + str(fiber_num[i]) + "\n"
            print("failed on the MAPS file.")
            #print(failed_maps)
            print("------------------------------------------")
            continue
        
        
        #logcube = fits.open('/home/celeste/Documents/astro_research/logcube_files/manga-'+ str(plate_num[i])+ '-' + str(fiber_num[i]) + '-LOGCUBE.fits.gz')
        try:
            logcube = fits.open('/media/celeste/Hypatia/MPL7/LOGCUBES/manga-'+ str(plate_num[i])+ '-' + str(fiber_num[i]) + '-LOGCUBE.fits.gz')
        except FileNotFoundError:
            failed_logcube = failed_logcube + str(plate_num[i]) + "-" + str(fiber_num[i]) + "\n"
            print("failed on the LOGCUBE file.")
            print(failed_logcube)
            print("------------------------------------------")
            continue

        ##assigns the plate id based on what is in the data cube
        plate_id = hdulist['PRIMARY'].header['PLATEIFU']

        ##gets official plate number
        plate_number = hdulist['PRIMARY'].header['PLATEID']
        fiber_number = hdulist['PRIMARY'].header['IFUDSGN']
        
        ##gets the hydrogen alpha and hydrogen beta data
        Ha = hdulist['EMLINE_GFLUX'].data[18,...]
        Hb = hdulist['EMLINE_GFLUX'].data[1,...]
        snmap = hdulist['SPX_SNR'].data
        fluxes = hdulist['EMLINE_GFLUX'].data
        #errs = hdulist['EMLINE_GFLUX_ERR'].data
        errs=(hdulist['EMLINE_GFLUX_IVAR'].data)**-0.5
        H_alpha = fluxes[18,:,:]
        Ha = H_alpha
        Ha_err = errs[18,:,:]
        OIII = fluxes[13,:,:]
        o3_err = errs[13,:,:]
        H_beta = fluxes[11,:,:]
        Hb_err = errs[11,:,:]
        n2_err = errs[19,:,:]
        NII = fluxes[19,:,:]
        s21 = fluxes[20,:,:]
        s22 = fluxes[21,:,:]
        s21_err = errs[20,:,:]
        s22_err = errs[21,:,:]
        
        #I band for contours
        contours_i = logcube['IIMG'].data
        contours_i_same = contours_i
        
        


        velocity = hdulist['EMLINE_GVEL'].data[18,...]
        velocity_err = (hdulist['EMLINE_GVEL_IVAR'].data[18,...])**-0.5
        
        ew_cut = hdulist['EMLINE_GEW'].data[18,...]
        
        #print(hdulist['PRIMARY'].header)

        #quit()

        ##Imports the thingy we need to get the images of the galaxy without having to download directly all the pictures. This also bypasses the password needed
        
        import requests
        """
        
        r = requests.get('https://data.sdss.org/sas/mangawork/manga/spectro/redux/v2_3_1/' + str(plate_num[i]) + '/stack/images/' + str(fiber_num[i]) + '.png', auth=('sdss', '2.5-meters'))
        

        ##Saves the image
        with open('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + str(plate_id) + '.png', 'wb') as fd:
	        for chunk in r.iter_content(chunk_size=128):
		        fd.write(chunk)
		        
		 """
        
		
	

        ##Also need OIII 5007 and NII 6584

        ##Gets the correct lines
        """"
        H_alpha = get_line_ew(hdulist, 'Ha-6564')
        H_beta = get_line_ew(hdulist, 'Hb-4862')


        OIII = get_line_ew(hdulist, 'OIII-5008')
        NII = get_line_ew(hdulist, 'NII-6585')
        """

        ##line ratios
        #x: n2/Halpha (plot the log)
        O_B = OIII/H_beta

        N_A = NII/H_alpha

        R = O_B/N_A	
        logR = np.log10(R)
        
        


        c0 = 0.281
        c1 = (-4.765)
        c2 = (-2.268)

        cs = np.array([c0, c1, c2])


        ##A lambda function, do not touch!
        x2logR = lambda x, cs: np.sum((c*x**p for p,c in zip(np.linspace(0, len(cs)-1, len(cs)), cs)))

        x2logR_zero = lambda x, cs, logR: x2logR(x, cs)-logR-0.001

        ##takes the log of the OH12 array
        #logOH12 = np.ma.array(np.zeros(logR.shape),mask=logR.mask)
        
        logOH12_old = 8.73-0.32*np.log10(R)
        logOH12, logOH12error = n2s2_dopita16_w_errs(H_alpha, NII, s21, s22, Ha_err, n2_err, s21_err, s22_err)
        is_starforming = is_sf_array(NII,H_alpha,OIII, H_beta)
        
        ##Finds the standard deviation and mean for future use	
        std_dev = np.std(Ha)
        mean = np.mean(Ha)	
		
        ##if it deviates too much from the mean it is removed		
        for j in range(len(Ha)):
	        for k in range(len(Ha[0])):
		        if (Ha[j][k] > std_dev*20+mean):
			        Ha[j][k] = np.nan
		

        ##Creates a shape that is the same size as the h-alpha array	
        shape = (Ha.shape[1])
        shapemap = [-.25*shape, .25*shape, -.25*shape, .25*shape]

        ##Changes the font size
        matplotlib.rcParams.update({'font.size': 20})
        #Second image we want?
        fig = plt.figure(figsize=(30,18), facecolor='white')
        #plt.plot(hd2[0], Ha[7])
        #sky coordinates relative to center
        #exent = .5

        ##places text on the plot
        plateifu = plate_id
        
        image = img.imread('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + str(plateifu) + '.png')
        
        drpall = t.Table.read('/home/celeste/Documents/astro_research/drpall-v2_3_1.fits')
        r = hdulist['SPX_ELLCOO'].data[0, ...]
        obj = drpall[drpall['plateifu']==plateifu][0]
        #nsa_sersic_ba for axis ratio
        #axis=drpall['nsa_sersic_ba'].data[0, ...]
        Re = obj['nsa_elpetro_th50_r']
        pa = obj['nsa_elpetro_phi']
        ba = obj['nsa_elpetro_ba']
        #radius of each spec in 
        r_Re = r/Re	
        r_Re = hdulist['SPX_ELLCOO'].data[1]

        print(plateifu)
        mass = math.log10(obj['nsa_elpetro_mass'])-np.log10(.49)
        #petrosian
        print("mass from data", mass)
        axis=obj['nsa_sersic_ba']
        #closest to 1, above .8
        #print("Axis ratio: ", axis) 
        #print("Mass is " + str(mass))
        
        
        if mass > 10.75:
            m = -0.16878698011761817
            b = 8.92174257450408
        if mass > 10.50 and mass <= 10.75:
            m = -0.19145937059393828
            b = 8.898917413495317
        if mass > 10.25 and mass <= 10.50:
            m = -0.16938127151421675
            b = 8.825998835583249
        if mass > 10.00 and mass <= 10.25:
            m = -0.1762907767970223
            b = 8.713865209075324
        if mass > 9.75 and mass <= 10.00:
            m = -0.14756252418062643
            b = 8.59167993089605
        if mass > 9.50 and mass <= 9.75:
            m = -0.07514461331863775
            b = 8.36144939226056
        if mass > 9.25 and mass <= 9.50:
            m = -0.05300368644036175
            b = 8.26602769508888
        if mass <= 9.25:
            m = -0.05059620593888811
            b = 8.147647436306206
        
        #print("Slope: ", m)
        #print("intercept: ", b)

            
        zeros= False
        for element in range(0, len(O_B)):
            for item in range(0, len(O_B[element])):
                if O_B[element][item] >= 0:
                    zeros = True
                    
        ##Adds another subplot with the plateifu
        a = fig.add_subplot(1,2,1)
        print("plate ifu for plotting image" + plateifu)
        #print(plateifu)
        try:
            image = img.imread('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + plateifu + '.png')
        except ValueError:
            print("No image.")
            print("========================================================================================")
        lum_img = image[:,:,0]
        #plt.subplot(121)
        imgplot = plt.imshow(image)
        plt.title("Galaxy "  + str(plate_number) + "-" + str(fiber_number))
        
#################################################################################
#
#   BPT Diagram Creator
#  
#################################################################################

#Sum the fluxes over a 3" center of the galaxy, put into is starforming
        
        ax_bpt = fig.add_subplot(1, 2, 2)
        if zeros == True:
            
            total=0
            sfr=0
            nsfr=0
            
            ax_bpt.set_aspect(1)
            ax_bpt.set_title("BPT Diagram")

        
        #Kewley
        X = np.linspace(-1.5, 0.3)
        Y = ((0.61/(X-0.47))+1.19)
        
        #Kauffmann
        Xk = np.linspace(-1.5,0.)
        Yk= (0.61/(Xk-0.05)+1.3)
       
        
        ax_bpt.plot(X, Y, '--', color = "red", lw = 1, label = "Kewley+01")
        ax_bpt.plot(Xk, Yk, '-', color = "blue", lw = 1, label = "Kauffmann+03")
        
        x=np.linspace(-0.133638005,0.75,100)
        y=2.1445*x+0.465
        ax_bpt.plot(x, y, '--', color = "green", lw = 1, label = "Seyfert/LINER")

        
        bpt_n2ha = np.log10(NII/H_alpha)
        bpt_o3hb = np.log10(OIII/H_beta)
        
        badpix = ((Ha/Ha_err) < 5) | ((H_beta/Hb_err) < 5) | ((OIII/o3_err) < 3) | ((NII/n2_err) < 3) | np.isinf(bpt_n2ha) | np.isinf(bpt_o3hb)
        bpt_n2ha[badpix] = np.nan
        bpt_o3hb[badpix] = np.nan
        
        bpt_o3hb95 = np.nanpercentile(bpt_o3hb, 98)
        bpt_o3hb5 = np.nanpercentile(bpt_o3hb, 2)
        
        xmin = np.nanmin(bpt_n2ha) - 0.1
        xmax = np.nanmax(bpt_n2ha) + 0.1
        #ymin = bpt_o3hb5 - 0.1
        #ymax = bpt_o3hb95 + 0.1
        ymin = np.nanmin(bpt_o3hb) - 0.1
        ymax = np.nanmax(bpt_o3hb) + 0.1
        plt.legend()
        
        
        
        #bad = is_starforming != 1
        #r_Rebpt = r_Re[bad]
        
        scatter_if(bpt_n2ha, bpt_o3hb, (is_starforming == 1) | (is_starforming == 0), c=r_Re, marker = ".", s = 65, alpha = 0.5, cmap = 'jet_r')
        #scatter_if(bpt_n2ha, bpt_o3hb, is_starforming == 0, c=r_Re, marker = ".", s = 65, alpha = 0.5, cmap = 'jet')
        
        ax_bpt.set_xlim(xmin, xmax)
        ax_bpt.set_ylim(ymin, ymax)
        ax_bpt.set_aspect((xmax-xmin)/(ymax-ymin))
        
        #cb_max = math.ceil(np.amax(r_Re))
        
        
        cb_bpt = plt.colorbar(shrink = .7)
        cb_bpt.set_label('r/$R_e$', rotation = 270, labelpad = 25)
        #plt.axes().set_aspect('equal')
        #plt.clim(0,20)
        #first_time = 1
        
        plt.tight_layout()
        
        try:
            plt.tight_layout()
        except ValueError:
            print("all NaN")
            print("==========================================================================================")
            #file_open = open("error_files.txt", "a")
            #first_time = 0
            #file_open.write(plateifu + "\n")
            #file.close()
            print("value error, all NaN")
            count_continue1=count_continue1+1
            continue

        
        Ha[is_starforming==0]=np.nan
        logOH12[is_starforming==0]=np.nan
                
        #print("total", total)
        #print("nsfr", nsfr)
        #print("sfr", sfr)

       

        #plt.show()
        plt.savefig('/home/celeste/Documents/astro_research/paper_plots/bpt_dia/bpt_dia_image_' + str(plate_id) + '.png')
