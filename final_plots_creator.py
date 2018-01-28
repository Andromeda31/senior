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
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams['axes.facecolor'] = 'white'



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

##Takes an input array and returns the values if a condition is met. Basically a glorified call to numpy.where

def array_if(array,condition=None):
    array1=np.array(array)
    if condition is None:
        return array1
    ret=array1[np.where(condition)]
    return ret


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
with open('/home/celeste/Documents/astro_research/thesis_git/adam_galaxies.txt') as f:
    file_names=[]
    for line in f:
        file_names.append(line)
with open('/home/celeste/Documents/astro_research/thesis_git/mass_data.txt') as f:
    mass_data=[]
    for line in f:
        mass_data.append(line)
#print(file_names)
##creates the empty arrays to append the names of the files in the folder
plate_num = []
fiber_num = []
split = []

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
    
#plate_num=['8137']
#fiber_num = ['12704']

for i in range(0, len(plate_num)): ##len(plate_num)
##for j in range(0, len(fiber_num)):
        print(plate_num[i] + '-' + fiber_num[i])
        print("Index: " + str(i))
        ##some black magic
        #hdulist = fits.open('/home/celeste/Documents/astro_research/keepers/manga-' + plate_num[i] + '-' + fiber_num[i] + '-MAPS-SPX-GAU-MILESHC.fits.gz')
        hdulist = fits.open('/home/celeste/Documents/astro_research/downloaded_data/MPL-6/manga-' + plate_num[i] + '-' + fiber_num[i] + '-MAPS-SPX-GAU-MILESHC.fits.gz')
        
        #logcube = fits.open('/home/celeste/Documents/astro_research/logcube_files/manga-10001-12702-LOGCUBE.fits')
        

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
        logOH12, logOH12error = pp04_o3n2_w_errs(NII, H_alpha, OIII, H_beta, n2_err, Ha_err, o3_err, Hb_err)
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

        print(plateifu)
        mass_adam = float(mass_data[i])
        mass = math.log10(obj['nsa_sersic_mass'])-np.log10(.49)
        print("mass from Adam", float(mass_adam))
        print("mass from data", mass)
        axis=obj['nsa_sersic_ba']
        #closest to 1, above .8
        #print("Axis ratio: ", axis) 
        #print("Mass is " + str(mass))
        
        if mass > 10.75:
            m = -0.06576751761269369
            b = 8.826541125371161
        if mass > 10.50 and mass <= 10.75:
            m = -0.08089943006458505
            b = 8.824202580328505
        if mass > 10.25 and mass <= 10.50:
            m = -0.10526083712519829
            b = 8.81995506359879
        if mass > 10.00 and mass <= 10.25:
            m = -0.12427905762033492
            b = 8.791749940700669
        if mass > 9.75 and mass <= 10.00:
            m = -0.08320549783022192
            b = 8.709540598630488
        if mass > 9.50 and mass <= 9.75:
            m = -0.052153896639010065
            b = 8.617252739908452
        if mass > 9.25 and mass <= 9.50:
            m = -0.047762559161765784
            b = 8.502512520775594
        if mass <= 9.25:
            m = 0.011644012238230033
            b = 8.359909467340257
        
        #print("Slope: ", m)
        #print("intercept: ", b)

            
        zeros= False
        for element in range(0, len(O_B)):
            for item in range(0, len(O_B[element])):
                if O_B[element][item] >= 0:
                    zeros = True
        
        ax = fig.add_subplot(2, 3, 4)
        if zeros == True:
            #fig = plt.figure()
            #ax.set_xscale("log")
            #ax.set_yscale("log")
            
            total=0
            sfr=0
            nsfr=0
            
            ax.set_aspect(1)
            ax.set_title("BPT Diagram")

        
        #Kewley
        X = np.linspace(-1.5, 0.3)
        Y = ((0.61/(X-0.47))+1.19)
        
        #Kauffmann
        Xk = np.linspace(-1.5,0.)
        Yk= (0.61/(Xk-0.05)+1.3)
       
        
        
        ax.plot(X, Y, '--', color = "red", lw = 1, label = "Kewley+01")
        ax.plot(Xk, Yk, '-', color = "blue", lw = 1, label = "Kauffmann+03")

        
        bpt_n2ha = np.log10(NII/H_alpha)
        bpt_o3hb = np.log10(OIII/H_beta)
        
        badpix = ((Ha/Ha_err) < 3) | ((OIII/o3_err) < 3) | ((H_beta/Hb_err) < 3) | ((NII/n2_err) < 3)
        bpt_n2ha[badpix] = np.nan
        bpt_o3hb[badpix] = np.nan
        
        xmin = np.nanmin(bpt_n2ha) - 0.1
        xmax = np.nanmax(bpt_n2ha) + 0.1
        ymin = np.nanmin(bpt_o3hb) - 0.1
        ymax = np.nanmax(bpt_o3hb) + 0.1
        
        
        
        scatter_if(bpt_n2ha, bpt_o3hb, is_starforming == 1, c=r_Re, marker = ".", s = 65, alpha = 0.5, cmap = 'jet')
        scatter_if(bpt_n2ha, bpt_o3hb,is_starforming == 0, c=r_Re, marker = ".", s = 65, alpha = 0.5, cmap = 'jet')
        
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        
        cb_max = math.ceil(np.amax(r_Re))
        
        cb = plt.colorbar(shrink = .7)
        cb.set_label('R/R_e', rotation = 270, labelpad = 25)
        #plt.axes().set_aspect('equal')
        plt.clim(0,cb_max)
        plt.tight_layout()
        
        Ha[is_starforming==0]=np.nan
        logOH12[is_starforming==0]=np.nan
                
        #print("total", total)
        #print("nsfr", nsfr)
        #print("sfr", sfr)

       



        ##Adds to the plot
        ##Is 2by3 and this is the second image
        a = fig.add_subplot(2,3,2)
        
        badpix = ((Ha/Ha_err) < 3)
        Ha_2d = Ha
        Ha_2d[badpix]=np.nan
        imgplot = plt.imshow(Ha_2d, cmap = "viridis", extent = shapemap, zorder = 1)
        plt.title("H-alpha Flux")
        cs=plt.gca().contour(Ha_2d, 8, colors='k', extent=shapemap, origin='upper', zorder = 2)
        plt.gca().clabel(cs, inline=1, fontsize=5)
        css = plt.gca().contour(r_Re,[1],extent=shapemap, colors='r', origin = 'upper', zorder = 2, z = 1)


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

        velocity[np.isinf(velocity)==True]=np.nan
        #velocity[(np.isnan(snmap)==True)|(snmap==0)]=np.nan
        velocity[mask != 0]=np.nan
        velocity_err[np.isinf(velocity_err)==True]=np.nan

        size=50
        ##Makes a mask
        #mask = (Ha==0).flatten()
        ##Adds another subplot
        a = fig.add_subplot(2,3,3)
        ##Makes  a scatter plot of the data
        badpix_vel = ((velocity_err) > 25)
        velocity[badpix_vel]=np.nan
        
        vel_min = np.nanpercentile(velocity, 5)
        vel_max = np.nanpercentile(velocity, 95)
        

        imgplot = plt.imshow(velocity, origin = "lower", cmap = "RdYlBu_r", extent = shapemap, vmin = vel_min, vmax = vel_max)
        
        plt.title("Gas Velocity")
        

        cb = plt.colorbar(shrink = .7)
        cb.set_label('km/s', rotation = 270, labelpad = 25)
        a.set_facecolor('white')
        


        ##Adds another subplot with the plateifu
        a = fig.add_subplot(2,3,1)
        print("plate ifu for plotting image" + plateifu)
        #print(plateifu)
        try:
            image = img.imread('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + plateifu + '.png')
        except ValueError:
            print("No image.")
        lum_img = image[:,:,0]
        #plt.subplot(121)
        imgplot = plt.imshow(image)
        plt.title("Galaxy "  + str(plate_number) + "-" + str(fiber_number))
        

        #print((~(logOH12.mask)).sum())

        for a in range(len(logOH12)):
	        for c in range(len(logOH12[0])):
		        if (logOH12[a][c] < 7):
			        logOH12[a][c] = None
        #print("Max R_E")
        #print(np.max(r_Re.flatten()))
       
			        

        a = fig.add_subplot(2, 3, 6)
        plt.xlabel("Effective Radii r/$R_e$")
        plt.ylabel('12+log(O/H)')
        #plt.scatter(r_Re.flatten(), logOH12.flatten(), s=size)
        idx = np.isfinite(r_Re.flatten()) & np.isfinite(logOH12.flatten())
        indarr=np.argsort(r_Re.flatten()[idx])
        
        yfit = [b + m * xi for xi in r_Re.flatten()[idx][indarr]]
        plt.tight_layout()
        plt.plot(r_Re.flatten()[idx][indarr], yfit, color = "red", zorder = 3)
        
        
        #plt.scatter(r_Re.flatten(), logOH12.flatten(), s=10, edgecolors = "black", color = "gray", zorder = 1)
        
        """
        H_alpha = fluxes[18,:,:]
        Ha = H_alpha
        Ha_err = errs[18,:,:]
        OIII = fluxes[13,:,:]
        o3_err = errs[13,:,:]
        H_beta = fluxes[11,:,:]
        Hb_err = errs[11,:,:]
        n2_err = errs[19,:,:]
        NII = fluxes[19,:,:]
        """
                
        cond_err = logOH12error.ravel()<np.nanpercentile(logOH12error.ravel(), 95)
        max_err = np.nanpercentile(logOH12error.ravel(), 95)
        condition = (logOH12error.flatten() < max_err) & ((Ha/Ha_err).flatten() > 3) & ((OIII/o3_err).flatten() >3) & ((H_beta/Hb_err).flatten() > 3) & ((NII/n2_err).flatten() >3)
        scatter_if(r_Re.flatten(), logOH12.flatten(), condition, s= 10, edgecolors = "black", color = "gray", zorder = 1)
        plt.title("Metallicity Gradient")

        plt.errorbar(r_Re.ravel()[condition], logOH12.ravel()[condition], yerr=logOH12error.ravel()[condition], fmt=None, errorevery = 45, capsize = 15, color = "green", zorder = 2)
        #imgplot = plt.scatter(r_Re.flatten(), logOH12.flatten(), s=41, color = "red")
        
        #m=np.polyfit(r_Re.flatten()[idx], logOH12.flatten()[idx], 1)
        #f = np.poly1d(m)
        #m= np.polyfit(r_Re.flatten(), logOH12.flatten(), 2)
        #indarr=np.argsort(r_Re.flatten()[idx])
        #plt.plot(r_Re.flatten()[idx][indarr], f(r_Re.flatten()[idx][indarr]))
        #plt.ylim(8,9)
        #plt.xlim(0, 2.25)
        #plt.show()
        

        shape = (logOH12.shape[1])
        shapemap = [-.5*shape, .5*shape, -.5*shape, .5*shape]

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
        
        badpix = (logOH12error > max_err) | ((Ha/Ha_err) < 3) | ((OIII/o3_err) < 3) | ((H_beta/Hb_err) < 3) | ((NII/n2_err) < 3) |  (ew_cut < 3)
        logOH12[badpix]=np.nan
        

        Ha_contour = Ha
        Ha_contour[badpix]=np.nan
        
        a = fig.add_subplot(2,3,5)
        plt.tight_layout()
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
            cs=plt.gca().contour(Ha_contour, 8, colors='k', extent=shapemap, origin='upper', zorder = 2)
            #plt.contour(logOH12, 20, colors='k')
        except ValueError:
            print("Value error! Skipping the log0H12 contour plotting....")
        plt.gca().clabel(cs, inline=1, fontsize=5)
        

        css = plt.gca().contour(r_Re,[1],extent=shapemap, colors='red', origin = 'upper', zorder = 2, z = 1, edgecolors = "black")


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
        
        plt.savefig('/home/celeste/Documents/astro_research/manga_images/final_images/star_faceon_' + plateifu +".png", bbox_inches = 'tight')
        #plt.savefig('/home/celeste/Documents/astro_research/thesis_git/show_adam/gaalxy_faceon_average_line_' + plateifu +".png", bbox_inches = 'tight')
        #plt.show()
        plt.close()
        print("Done with this one.")
        print("--------------------------------------------------")
        #leedle
        """
        leedle
        if i >5:
            leedle
        """

        #quit()

	    #plt.savefig('/home/celeste/Documents/astro_research/astro_images/six_galaxies/ha_flux_Z_' + plate_num[i] + '_' + fiber_num[i], bbox_inches = 'tight')


