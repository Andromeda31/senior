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
    
#plate_num=['9183']
#fiber_num = ['9102']

count_continue1 = 0
count_continue2 = 0
count_continue3 = 0
failed_maps = "failed maps\n"
failed_logcube = "failed_logcube\n"

for i in range(230, len(plate_num)): ##len(plate_num)
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
            print(failed_maps)
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
        
        logOH12 = 8.73-0.32*np.log10(R)
        logOH12, logOH12error = n2s2_dopita16_w_errs(H_alpha, NII, s21, s22, Ha_err, n2_err, s21_err, s22_err)
        is_starforming = is_sf_array(NII,H_alpha,OIII, H_beta)
        
        #plt.errorbar(rad.ravel(), met.ravel(), yerr=logOH12error.ravel(), fmt='o', errorevery=10)
        

        ##Not sure what this is either
        ##This is the logOH12 issue here
        #print(*logR.shape)
        '''
        for ii, jj in np.ndindex(*logR.shape):
	        if logR.mask[ii,jj]:continue
	        try:
		        logOH12[ii, jj] = newton(x2logR_zero,x0=.1,args=(cs, logR[ii,jj]))+8.69
	        except (SystemExit, KeyboardInterrupt):
		        raise
	        except:
		        pass
		'''
		
        #print(len(logOH12))
		        
	
        #print("average ", np.ma.median(np.log10(R)))

        #xplus = (-c1 + np.sqrt(c1**2 - 4*(c0-np.log10(R))*(c2)))/(2*c2)
        #xminus =(-c1 - np.sqrt(c1**2 - 4*(c0-np.log10(R))*(c2)))/(2*c2)
        #print((~(xplus.mask)).sum())
        #print((xminus<0).sum())
        #pows = np.linspace(0, 2, 3)
        #x = np.sum([c*(R**p) for p, c in zip(pows,cs)])
        #logOH12 = np.log10(np.absolute(xminus))+8.69
        #print((~(logOH12.mask)).sum())
		
        ##Gets rid of the Ha value if it is ridiculously high		
        
        """
        for l in range(len(Ha)):
	        for m in range(len(Ha[0])):
		        if (Ha[l][m] > 150):
			        Ha[l][m] = np.nan
		        if (Ha[l][m] < 0):
		            Ha[l][m] = np.nan
		"""
	    
		
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
        
        #Add bpt diagram

        #check for pre existence of bpt plot:

        #bpt_file=Path("/home/celeste/Documents/astro_research/astro_images/bpt_diagrams/bpt_" + str(plateifu) + ".png")
        
        image = img.imread('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + str(plateifu) + '.png')
        
        
        """
        if bpt_file.is_file():
            a=fig.add_subplot(2, 3, 4)
            image=img.imread("/home/celeste/Documents/astro_research/astro_images/bpt_diagrams/bpt_" + str(plateifu) + ".png")
            lum_img = image[:,:,0]
            a.axis('off')
            imgplot=plt.imshow(image[90:500, 30:450])
            need_new_plot=False
        else:
            print("Image not found. Plotting BPT manually...")
            print(plateifu)
        """
        
        
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
        
#################################################################################
#
#   Mass of galaxy to get slope of average profile from Belfiore
#  
#################################################################################
        
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

#################################################################################
#
#   BPT Diagram Creator
#  
#################################################################################

#Sum the fluxes over a 3" center of the galaxy, put into is starforming
        
        ax_bpt = fig.add_subplot(2, 3, 4)
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

       




#################################################################################
#
#   Halpha Flux 2-D
#  
#################################################################################

        a = fig.add_subplot(2,3,2)
        
        badpix = ((Ha/Ha_err) < 3)
        Ha_2d = Ha
        Ha_2d[badpix]=np.nan
        contours_i[badpix]=np.nan
        imgplot = plt.imshow(Ha_2d, cmap = "viridis", extent = shapemap, zorder = 1)
        plt.title("H-alpha Flux")
        #cs=plt.gca().contour(Ha_2d, 8, colors='k', alpha = 0.3, linewidths= [1], extent=shapemap, origin='upper', zorder = 3)
        csss=plt.gca().contour(contours_i, 8, colors = 'k', alpha = 0.6, extent = shapemap, origin = 'upper', zorder = 4)
        #plt.gca().clabel(cs, inline=1, fontsize=5)
        css = plt.gca().contour(r_Re*2,[2],extent=shapemap, colors='r', origin = 'upper', zorder = 2, z = 2)
        
        """
        
        plt.close()
        plt.imshow(r_Re)
        plt.colorbar()
        plt.show()
        print(r_Re)
        
        """

        #plt.gca().clabel(css, inline=1)
        plt.gca().invert_yaxis()
        
        #ec = patches.Ellipse(xy=(0,0), width=Re*ba*2, height=Re*2, angle=pa, linewidth = 4, edgecolor='k',fill=False, zorder = 3)
        
        #a.add_patch(ec)
        
        plt.xlabel('Arcseconds')
        plt.ylabel('Arcseconds')
        cb = plt.colorbar(shrink = .7)
        #cb.set_label('F(H$\alpha$)', rotation = 270, labelpad = 25)
        cb.set_label('H-alpha flux [$10^{17} {erg/s/cm^2/pix}$]', rotation = 270, labelpad = 25)
        #plt.xlim, plt.ylim
        
        mask = hdulist['EMLINE_GVEL_MASK'].data[18,:,:]
        
        #make white be zero
        #find min and max of velocity, whichever one's absolute value is larger, use as min and max

        velocity[np.isinf(velocity)==True]=np.nan
        #velocity[(np.isnan(snmap)==True)|(snmap==0)]=np.nan
        velocity[mask != 0]=np.nan
        velocity_err[np.isinf(velocity_err)==True]=np.nan

        size=50
        ##Makes a mask
        #mask = (Ha==0).flatten()
        ##Adds another subplot

#################################################################################
#
#   Velocity Map
#  
#################################################################################
        
        a = fig.add_subplot(2,3,3)
        ##Makes  a scatter plot of the data
        badpix_vel = ((velocity_err) > 25)
        velocity[badpix_vel]=np.nan
        
        vel_min = np.nanpercentile(velocity, 5)
        vel_max = np.nanpercentile(velocity, 95)
        
        if abs(vel_min) > abs(vel_max):
            vel_final = abs(vel_min)
            want = vel_max
        else:
            vel_final = abs(vel_max)
            want = vel_min

            
        print(vel_min)
        print(vel_max)
        
        

        imgplot = plt.imshow(velocity, origin = "lower", cmap = "RdYlBu_r", extent = shapemap, vmin = -vel_final, vmax = vel_final)
        css = plt.gca().contour(r_Re*2,[2],extent=shapemap, colors='green', origin = 'lower', zorder = 2, z = 1)
        csss=plt.gca().contour(contours_i, 8, colors = 'k', alpha = 0.6, extent = shapemap, zorder = 4)
        
        plt.title("Gas Velocity")
        

        cb = plt.colorbar(shrink = .7)
        if ((vel_min <=0) and (vel_max <=0)):
            plt.clim(-vel_final, want)
        else:
            plt.clim(-vel_final,vel_final)
        cb.set_label('km/s', rotation = 270, labelpad = 25)
        a.set_facecolor('white')
        
#################################################################################
#
#   Plots Galaxy Image
#  
#################################################################################


        ##Adds another subplot with the plateifu
        a = fig.add_subplot(2,3,1)
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
#   Metallicity Gradient with fitted Lines
#  
#################################################################################
       
			        

        a = fig.add_subplot(2, 3, 6)
        plt.xlabel("Effective Radii r/$R_e$")
        plt.ylabel('12+log(O/H)')
        #plt.scatter(r_Re.flatten(), logOH12.flatten(), s=size)
        idx = np.isfinite(r_Re.flatten()) & np.isfinite(logOH12.flatten())
        indarr=np.argsort(r_Re.flatten()[idx])
        
        yfit = [b + m * xi for xi in r_Re.flatten()[idx][indarr]]
        try:
            plt.tight_layout()
        except ValueError:
            print("all NaN")
            print("==========================================================================================")
            count_continue2=count_continue2+1
        plt.plot(r_Re.flatten()[idx][indarr], yfit, color = "red", zorder = 3, label = 'Expected profile for stellar mass')
        
        
     
        def func(x, m, b):
            return m*x+b
        
        abundance = logOH12
        radius = r_Re
            
        trialX = np.linspace(np.amin(radius.ravel()), np.amax(radius.ravel()), 1000)
        
        cond_err = logOH12error.ravel()<np.nanpercentile(logOH12error.ravel(), 95)
        max_err = np.nanpercentile(logOH12error.ravel(), 95)
        condition = (logOH12error.flatten() < max_err) & ((Ha/Ha_err).flatten() > 3) & ((s22/s22_err).flatten() >3) & ((s21/s21_err).flatten() > 3) & ((NII/n2_err).flatten() >3)
        scatter_if(r_Re.flatten(), logOH12.flatten(), condition, s= 10, edgecolors = "black", color = "gray", zorder = 1)
        plt.errorbar(r_Re.ravel()[condition], logOH12.ravel()[condition], yerr=logOH12error.ravel()[condition], fmt=None, errorevery = 45, capsize = 15, color = "black", zorder = 2)
        
        condition2 = (logOH12error < max_err) & ((Ha/Ha_err) > 3) & ((s22/s22_err) >3) & ((s21/s21_err) > 3) & ((NII/n2_err) >3)
        logOH12_2=logOH12.copy()
        logOH12_2[condition2==False]=np.nan
        
        rad_pix=hdulist['SPX_ELLCOO'].data[0,:,:]*2.0 #since there are 2 pixels/arcsec
        rad, rp, n, sig =radial_profile(image=logOH12_2,distarr=rad_pix)
        rad=rad/(2*Re) #This is now in units of Re.
        
        
        valid = ~(np.isnan(rad) | np.isnan(rp) | np.isinf(rad) | np.isinf(rp) | ((rad < .5) | (rad > 2) ) | (n < 10))
        
        plt.plot(rad, rp, '.m', label = 'binned median', markersize =7, marker = 'D')
        
        try:
            popt, pcov = curve_fit(func, rad[valid], rp[valid], check_finite = True)
        except TypeError:
            print("Improper input: N=2 must not exceed M=0")
            print("==========================================================================================")
            plt.close('all')
            count_continue=count_continue3+1
            continue
        plt.plot(rad[valid], func(rad[valid], *popt), 'cyan', label = '0.5-2 $R_e$ fit', linewidth = 5)

        
        plt.legend()
        plt.title("Metallicity Gradient")
        plt.xlim(xmin = 0)
        
#################################################################################
#
#   Metallicity Radial Plot 3-D
#  
#################################################################################
     
        

        shape = (logOH12.shape[1])
        shapemap = [-.25*shape, .25*shape, -.25*shape, .25*shape]

        #logOH12_flat = logOH12.flatten()
        

        #print("--------")
			
        matplotlib.rcParams.update({'font.size': 20})
        #Second image we want?
        #plt.plot(hd2[0], Ha[7])
        #sky coordinates relative to center
        #exent = .5
        logOH12[np.isinf(logOH12)==True]=np.nan
        
        minimum = np.nanpercentile(logOH12, 5)
        maximum = np.nanpercentile(logOH12, 95)
        median = np.nanpercentile(logOH12, 50)
        
        badpix = (logOH12error > max_err) | ((Ha/Ha_err) < 3) | ((s22/s22_err) < 3) | ((s21/s21_err) < 3) | ((NII/n2_err) < 3) |  (ew_cut < 3)
        logOH12[badpix]=np.nan
        

        Ha_contour = Ha
        Ha_contour[badpix]=np.nan
        contours_i_same[badpix]=np.nan
        
        a = fig.add_subplot(2,3,5)
        try:
            plt.tight_layout()
        except ValueError:
            print("all NaN")
            print("==========================================================================================")
            
        if ((maximum - minimum) < .2):
            maximum = median + 0.1
            minimum = median - 0.1
        
        
        #badpix=np.where((np.isnan(radius)==True) | (np.isnan(abundance)==True)
        #issf[badpix]=np.nan

            
        imgplot = plt.imshow(logOH12, cmap = "viridis", extent = shapemap, vmin = minimum, vmax = maximum, zorder = 1)

        plt.title("Metallicity Map")
        
            
        #write if statement 
        """
        good=np.where(logOH12!=None)
        total=len(good[1])
        if total >=100:
            cs=plt.gca().contour(logOH12, 8, colors='k', extent = shapemap, origin="upper")
        """
        try:
            #cs=plt.gca().contour(Ha_contour, 8, colors='k', alpha = 0.3, linewidths= [1], extent=shapemap, origin='upper', zorder = 2)
            csss=plt.gca().contour(contours_i_same, 8, colors = 'k', alpha = 0.6, extent = shapemap, origin = 'upper', zorder = 3)
            #plt.contour(logOH12, 20, colors='k')
        except ValueError:
            print("Value error! Skipping the log0H12 contour plotting....")
            print("==========================================================================================")
        #plt.gca().clabel(cs, inline=1, fontsize=5)
        

        css = plt.gca().contour(r_Re*2,[2],extent=shapemap, colors='red', origin = 'upper', alpha = .6, zorder = 2, z = 1, edgecolors = "black")


        #plt.gca().clabel(css)

        plt.gca().invert_yaxis()
        #plt.gca().invert_yaxis()

        plt.xlabel('Arcseconds')
        plt.ylabel('Arcseconds')
        cb = plt.colorbar(shrink = .7)
        cb.set_label('12+log(O/H)', rotation = 270, labelpad = 25)
        #plt.xlim, plt.ylim
        
        """
        
        plt.close()
        
        fig = plt.figure(figsize=(30,18))
        
        a = fig.add_subplot(1, 1, 1)
        imgplot = plt.imshow(r_Re, cmap = "viridis", extent = shapemap, zorder = 1)
        cs=plt.gca().contour(r_Re, 8, colors='k', extent=shapemap, origin='upper', zorder = 2)
        plt.gca().clabel(cs, inline=1, fontsize=5)

        
        ec = patches.Ellipse(xy=(0,0), width=Re*ba*2, height=Re*2, angle=pa+90, linewidth = 4, edgecolor='k',fill=False, zorder = 3)
        
        
        a.add_patch(ec)
        
        plt.gca().invert_yaxis()
        
        plt.show()
        """
     
       

            
        ##Saves the final image
        print("Saving...")
        print(fig.dpi)
        plt.savefig('/home/celeste/Documents/astro_research/manga_images/final_images/TERABYTE/MPL7/logcube_' + plateifu +"_v5.1.png", bbox_inches = 'tight', dpi = 90)
        #plt.savefig('/home/celeste/Documents/astro_research/manga_images/final_images/star_faceon_' + plateifu +"_v4.1.png", bbox_inches = 'tight')
        #plt.savefig('/home/celeste/Documents/astro_research/thesis_git/show_adam/gaalxy_faceon_average_line_' + plateifu +".png", bbox_inches = 'tight')
        #plt.show()
        #plt.close('all')
        plt.close('all')
        print("Done with this one.")
        print("--------------------------------------------------")
print(failed_logcube)
print(failed_maps)
