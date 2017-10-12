##lots of imports. Some are unnecessary but I left a lot just to be safe...
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
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
	    import requests
	    r = requests.get('https://data.sdss.org/sas/mangawork/manga/spectro/redux/v2_0_1/' + str(plate_number) + '/stack/images/' + str(fiber_number) + '.png', auth=('sdss', '2.5-meters'))

	    ##Saves the image
	    with open('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + str(plate_id) + '.png', 'wb') as fd:
		    for chunk in r.iter_content(chunk_size=128):
			    fd.write(chunk)
			
		

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

	    R = O_B/N_A	
	    logR = np.log(R)


	    c0 = 0.281
	    c1 = (-4.765)
	    c2 = (-2.268)

	    cs = np.array([c0, c1, c2])
	
	
	    ##A lambda function, do not touch!
	    x2logR = lambda x, cs: np.sum((c*x**p for p,c in zip(np.linspace(0, len(cs)-1, len(cs)), cs)))

	    x2logR_zero = lambda x, cs, logR: x2logR(x, cs)-logR-0.001

	    ##takes the log of the OH12 array
	    logOH12 = np.ma.array(np.zeros(logR.shape),mask=logR.mask)
	    

	    ##Not sure what this is either
	    ##This is the logOH12 issue here
	    print(*logR.shape)
	    for ii, jj in np.ndindex(*logR.shape):
		    if logR.mask[ii,jj]:continue
		    try:
			    logOH12[ii, jj] = newton(x2logR_zero,x0=.1,args=(cs, logR[ii,jj]))+8.69
		    except (SystemExit, KeyboardInterrupt):
			    raise
		    except:
			    pass
	    print(len(logOH12))
			    
		
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
	    for l in range(len(Ha)):
		    for m in range(len(Ha[0])):
			    if (Ha[l][m] > 150):
				    Ha[l][m] = 0
			
	    ##Finds the standard deviation and mean for future use	
	    std_dev = np.std(Ha)
	    mean = np.mean(Ha)	
			
	    ##if it deviates too much from the mean it is removed		
	    for j in range(len(Ha)):
		    for k in range(len(Ha[0])):
			    if (Ha[j][k] > std_dev*20+mean):
				    Ha[j][k] = 0
			
	
	    ##Creates a shape that is the same size as the h-alpha array	
	    shape = (Ha.shape[1])
	    shapemap = [-.5*shape, .5*shape, -.5*shape, .5*shape]

	    ##Changes the font size
	    matplotlib.rcParams.update({'font.size': 20})
	    #Second image we want?
	    fig = plt.figure(figsize=(30,18))
	    #plt.plot(hd2[0], Ha[7])
	    #sky coordinates relative to center
	    #exent = .5

	    ##places text on the plot
	    plateifu = plate_id
	
	    ##Adds to the plot
	    ##Is 2by3 and this is the second image
	    a = fig.add_subplot(2,3,2)
	    imgplot = plt.imshow(Ha, cmap = "viridis", extent = shapemap)
	    plt.gca().invert_yaxis()
	    plt.xlabel('Arcseconds')
	    plt.ylabel('Arcseconds')
	    cb = plt.colorbar(shrink = .7)
	    cb.set_label('H-alpha flux [$10^{17} {erg/s/cm^2/pix}$]',
	    rotation = 270, labelpad = 25)
	    #plt.xlim, plt.ylim


	    drpall = t.Table.read('drpall-v2_0_1.fits')
	    r = hdulist['SPX_ELLCOO'].data[0, ...]
	    obj = drpall[drpall['plateifu']==plateifu][0]
	    Re = obj['nsa_elpetro_th50_r']
	    #radius of each spec in 
	    r_Re = r/Re	
	
	    mass = obj['nsa_sersic_mass']

	
	    fig.text(0.12, 0.45, 'ID: ' + plateifu, fontsize=30)
	    fig.text(0.12, 0.40, 'Mass:', fontsize=30)
	    fig.text(0.12, 0.35, str(mass) + "$\log_{10}{\frac{M}{M_{\odot}}}$", fontsize = 25)
	    fig.text(0.12, 0.30, 'MARVIN Link: ', fontsize=30)
	    fig.text(0.12, 0.28, 'https://sas.sdss.org/marvin2/galaxy/' + plateifu, fontsize=15)

	
	    ##Makes a mask
	    mask = (Ha==0).flatten()
	    ##Adds another subplot
	    a = fig.add_subplot(2,3,3)
	    ##Makes  a scatter plot of the data
	    imgplot = plt.scatter(r_Re.flatten()[~mask], Ha.flatten()[~mask])
	    plt.xlabel("Effective Radii r/$R_e$")
	    plt.ylabel('H-alpha flux [$10^{17} {erg/s/cm^2/pix}$]')
	    ##Finds the maximum and puts it in the plot
	    r_max = np.amax(r_Re.flatten()[~mask])
	    Ha_max = np.amax(Ha.flatten()[~mask])
	    plt.xlim(0,r_max)
	    plt.ylim(0,Ha_max)

	    ##Adds another subplot with the plateifu
	    a = fig.add_subplot(2,3,1)
	    image = img.imread('/home/celeste/Documents/astro_research/astro_images/marvin_images/' + plateifu + '.png')
	    lum_img = image[:,:,0]
	    #plt.subplot(121)
	    imgplot = plt.imshow(image)


	    #print((~(logOH12.mask)).sum())
	
	    for a in range(len(logOH12)):
		    for b in range(len(logOH12[0])):
			    if (logOH12[a][b] < 7):
				    logOH12[a][b] = None
				    

	    a = fig.add_subplot(2, 3, 6)
	    imgplot = plt.scatter(r_Re.flatten(), logOH12.flatten()) #.flatten
	    plt.xlabel("Effective Radii r/$R_e$")
	    plt.ylabel('12+log(O/H)')
	    plt.scatter(r_Re.flatten(), logOH12.flatten())
	    #plt.ylim(8,9)
	    #plt.xlim(0, 2.25)
	    #plt.show()

	    shape = (logOH12.shape[1])
	    shapemap = [-.5*shape, .5*shape, -.5*shape, .5*shape]

	    #logOH12_flat = logOH12.flatten()
	    
	    print(len(logOH12))
	    print(logOH12[30])
	
	    print("--------")
				
	    matplotlib.rcParams.update({'font.size': 20})
	    #Second image we want?
	    #plt.plot(hd2[0], Ha[7])
	    #sky coordinates relative to center
	    #exent = .5
	    a = fig.add_subplot(2,3,5)
	    imgplot = plt.imshow(logOH12, cmap = "viridis", extent = shapemap)
	    plt.gca().invert_yaxis()
	    plt.xlabel('Arcseconds')
	    plt.ylabel('Arcseconds')
	    cb = plt.colorbar(shrink = .7)
	    cb.set_label('log(OH) + 12', rotation = 270, labelpad = 25)
	    #plt.xlim, plt.ylim
	    

	    ##Saves the final image
	    plt.savefig('/home/celeste/Documents/astro_research/astro_images/keeper_images/6_images_' + plateifu, bbox_inches = 'tight')
	    #plt.show()
	    plt.close() 
	    #quit()

	    #plt.savefig('/home/celeste/Documents/astro_research/astro_images/six_galaxies/ha_flux_Z_' + plate_num[i] + '_' + fiber_num[i], bbox_inches = 'tight')

'''
hdulist = fits.open('/home/celeste/Documents/astro_research/manga-10001-12702-LOGCUBE.fits')
lam = hdulist['WAVE'].data
flux = hdulist['FLUX'].data


#print(len(flux[:, 22, 22]))
#print((flux[:, 22, 22]))

#print(len(lam))
#print(lam)

plt.plot(lam, flux[:, 22, 22])
plt.xlabel('wavelength [Angstrom]')
plt.ylabel('flux')
#plt.show()
plt.close('all')

#Finding 5500A

for x in range(0, (len(lam)-1)):
    if lam[x] < 5501 and lam[x] > 5499:
	    print(x)

print(flux[1815])

im_axes = [0, 71, 0, 72]


rband = hdulist['RIMG'].data
im1 = plt.imshow(rband, cmap='viridis')
#plt.xlabel('Radius')
#plt.ylabel('Flux')
plt.colorbar(im1)
plt.savefig('rband')
plt.show()
plt.close()

im2 = plt.imshow(flux[1815], cmap='viridis')
plt.ylabel('Flux')
plt.colorbar(im2)
plt.savefig('flux_5500')
plt.show()
plt.close('all')
'''



