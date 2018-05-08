
from astropy.io import fits

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

#############################################################################################

'''COMMENT  *** Column names ***                                                   COMMENT                                                                         TTYPE1  = 'PLATEIFU'           /                                                TTYPE2  = 'MANGAID '           /                                                TTYPE3  = 'RA      '           /                                                TTYPE4  = 'DEC     '           /                                                TTYPE5  = 'Z       '           /                                                TTYPE6  = 'FNUGRIZ_ABSMAG'     /                                                TTYPE7  = 'LOG_MASS'           /                                                TTYPE8  = 'SUBJECT_ID'         /                                                TTYPE9  = 'NCLASS  '           /                                                TTYPE10 = 'BAD_RE  '           /                                                TTYPE11 = 'BAD_RE_ERR'         /                                                TTYPE12 = 'PA_SHIFT'           /                                                TTYPE13 = 'PA_SHIFT_ERR'       /                                                TTYPE14 = 'KINE_TWIST'         /                                                TTYPE15 = 'KINE_TWIST_ERR'     /                                                TTYPE16 = 'DISTURBED_KINE'     /                                                TTYPE17 = 'DISTURBED_KINE_ERR' /                                                TTYPE18 = 'MERGING '           /                                                TTYPE19 = 'MERGING_ERR'        /                                                TTYPE20 = 'SYMMETRIC_OH'       /                                                TTYPE21 = 'SYMMETRIC_OH_ERR'   /                                                TTYPE22 = 'DISTORTED_OH'       /                                                TTYPE23 = 'DISTORTED_OH_ERR'   /                                                TTYPE24 = 'CHAOTIC_OH'         /                                                TTYPE25 = 'CHAOTIC_OH_ERR'     /                                                TTYPE26 = 'BAD_OH  '           /                                                TTYPE27 = 'BAD_OH_ERR'         /                                                TTYPE28 = 'LOW_KNOTS'          /                                                TTYPE29 = 'LOW_KNOTS_ERR'      /                                                TTYPE30 = 'HIGH_KNOTS'         /                                                TTYPE31 = 'HIGH_KNOTS_ERR'     /                                                TTYPE32 = 'LINEAR_OHGRAD'      /                                                TTYPE33 = 'LINEAR_OHGRAD_ERR'  /                                                TTYPE34 = 'SLOPE_CHANGE'       /                                                TTYPE35 = 'SLOPE_CHANE_ERR'    /                                                TTYPE36 = 'IRREGULAR_OHGRAD'   /                                                TTYPE37 = 'IRREGULAR_OHGRAD_ERR' /                                              TTYPE38 = 'BAD_OHGRAD'         /                                                TTYPE39 = 'BAD_OHGRAD_ERR'     /
'''

hdul = fits.open('/home/celeste/Documents/astro_research/thesis_git/zooinverse_summary_v0.fits')

print(hdul[1].header)
hdr = hdul[1].data


print(hdr[0][0])

interacting = 0
bad_rad = 0
mis_align = 0
twist = 0
disturbed_kin = 0
total = 0

symm_oh = 0
distor_oh = 0
chaotic_oh = 0
bad_oh = 0

symm_arr = []
distor_arr = []
chaotic_arr = []
bad_arr = []
mass_arr = []

lin_grad = 0
slope_ch = 0
irr_grad = 0
bad_grad = 0

m = 0
b = 0

slope = []

for x in range(0, len(hdr)):
    plateid = hdr[x][0]
    
    bad_rad = bad_rad + hdr[x][9]
    mis_align = mis_align + hdr[x][11]
    twist = twist + hdr[x][13]
    disturbed_kin = disturbed_kin + hdr[x][15]
    interacting = interacting + hdr[x][17]
    
    symm_oh = symm_oh + int(round(hdr[x][19]+.01))
    distor_oh = distor_oh + int(round(hdr[x][21]+.01))
    chaotic_oh = chaotic_oh + int(round(hdr[x][23]+.01))
    bad_oh = bad_oh + int(round(hdr[x][25]+.01))
    
    log_mass = hdr[x][6]
    
    mass_arr.insert(x, log_mass)
    symm_arr.insert(x, int(round(hdr[x][19]+.01)))
    distor_arr.insert(x, int(round(hdr[x][21]+.01)))
    chaotic_arr.insert(x, int(round(hdr[x][23]+.01)))
    bad_arr.insert(x, int(round(hdr[x][25]+.01)))
    
    bad_grad = bad_grad + int(round(hdr[x][37]+.01))
    irr_grad = irr_grad + int(round(hdr[x][35]+.01))
    slope_ch = slope_ch + int(round(hdr[x][33]+.01))
    lin_grad = lin_grad + int(round(hdr[x][31]+.01))
    
    mass = log_mass
    
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
            
    slope.insert(x, m)
    
    total = total + 1
    
#####################################################################

### BAR GRAPH FROM THE FIRST ZOONIVERSE QUESTION

#####################################################################


y = [bad_rad, mis_align, twist, disturbed_kin, interacting]
y = [(i/total)*100 for i in y]
N = len(y)
x = range(N)
width = .5
fig, ax = plt.subplots()
rects1 = ax.bar(x, y, width, color = "blue")
ax.set_xticklabels(('0','Bad Radius\nor Inclination', 'Misalignment of\nPosition Angle', 'Kinematic\nTwist', 'Highly\nDisturbed\nKinematics', 'Interacting or\nMerging'))
ax.set_title("Data from First Zooniverse Question")
ax.set_ylabel('Percent of Galaxies with Said Feature')

def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height,
                '%1.1d%%' % height,
                ha='center', va='bottom')
                
autolabel(rects1)
plt.savefig("/home/celeste/Documents/astro_research/paper_plots/analysis_plots/bar_graph_q1.png")
plt.close('all')

#####################################################################

### PIE GRAPH FROM DATA FROM QUESTION 2

#####################################################################

'''
print(total)
print(symm_oh + distor_oh + chaotic_oh + bad_oh)
print(symm_oh)
print(distor_oh)
print(chaotic_oh)
print(bad_oh)
'''



# Data to plot
labels = ['Symmetric O/H Map', 'Distorted O/H Contours', 'Disturbed, Chaotic O/H', 'Bad O/H Map']
sizes = [symm_oh, distor_oh, chaotic_oh, bad_oh]
colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
explode = (0.1, 0, 0, 0)  # explode 1st slice
 
# Plot
plt.pie(sizes, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90)
plt.legend(labels, loc = 'best')
plt.tight_layout()
plt.title('Best Descriptors of O/H Map')
plt.axis('equal')
#plt.show()
plt.savefig("/home/celeste/Documents/astro_research/paper_plots/analysis_plots/pie_graph_q2a.png")
plt.close('all')

#####################################################################

### LINE PLOT FROM DATA FROM QUESTION 2

#UNFINISHED

#####################################################################

'''
print(mass_arr[0])
print(slope[0])

sorted_list_slope_by_list_mass_arr = sorted(mass_arr, key=lambda x: slope.index(x))

print(mass_arr[0])
print(slope[0])

asdfsd
'''

max_mass = 11
min_mass = 8.5
binsize = 0.5

numbin = int((max_mass-min_mass)/(binsize))+1

bins = [[0 for x in range(4)] for x in range(numbin)]



for i in range(0, len(mass_arr)):
    ind = math.floor(2*(mass_arr[i]-min_mass+(binsize/2)))
    bins[ind][0] += symm_arr[i]
    bins[ind][1] += distor_arr[i]
    bins[ind][2] += chaotic_arr[i]
    bins[ind][3] += bad_arr[i]
    
sums = np.sum(bins, axis = 1)/100

x = np.linspace(min_mass, max_mass, numbin)

bins = np.asarray(bins)

#print(bins[:,0])

print(x)

plt.xlim(min_mass-.15, max_mass+.15)
labels = ['Symmetric O/H Map', 'Distorted O/H Contours', 'Disturbed, Chaotic O/H', 'Bad O/H Map']
sym = plt.plot(x, bins[:,0]/sums, c = "lightpink", linewidth = 4, label = 'Symmetric O/H Map')
dist = plt.plot(x, bins[:,1]/sums, c = "coral", linewidth = 4, label = 'Distorted O/H Contours')
chaos = plt.plot(x, bins[:,2]/sums, c = "rebeccapurple", linewidth = 4, label = 'Disturbed, Chaotic O/H')
bad = plt.plot(x, bins[:,3]/sums, c = "deepskyblue", linewidth = 4, label ='Bad O/H Map', alpha = 0.5)

capsize = 7

plt.errorbar(x, bins[:,0]/sums, np.sqrt(bins[:,0])/sums, ecolor = "lightpink", label = 'Symmetric O/H Map', visible = False, capsize = capsize)
plt.errorbar(x, bins[:,1]/sums, np.sqrt(bins[:,1])/sums, ecolor = "coral", label = 'Distorted O/H Contours', visible = False, capsize = capsize)
plt.errorbar(x, bins[:,2]/sums, np.sqrt(bins[:,2])/sums, ecolor = "rebeccapurple", label = 'Disturbed, Chaotic O/H', visible = False, capsize = capsize)
plt.errorbar(x, bins[:,3]/sums, np.sqrt(bins[:,3])/sums, ecolor = "deepskyblue", label ='Bad O/H Map', visible = False, capsize = capsize)

plt.legend(labels, loc = 'best')
plt.ylabel('Percent')
plt.xlabel('$log_{10}$ Stellar Mass')

#plt.show()
plt.savefig("/home/celeste/Documents/astro_research/paper_plots/analysis_plots/line_plot_q2.png")
plt.close('all')

#####################################################################

### PIE CHART FROM QUESTION 3

#####################################################################

labels = ['Linear O/H Gradient', 'Change in Slope of O/H Gradient', 'Irregular Fit to O/H Gradient', 'Bad O/H Gradient']
sizes = [lin_grad, slope_ch, irr_grad, bad_grad]
colors = ['mediumspringgreen', 'darkviolet', 'fuchsia', 'crimson']
explode = (0.1, 0, 0, 0)  # explode 1st slice
 
# Plot
plt.pie(sizes, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90)
plt.legend(labels, loc = 'best')
plt.tight_layout()
plt.axis('equal')
plt.title('Best Descriptors of O/H Radial Profile')
#plt.show()
plt.savefig("/home/celeste/Documents/astro_research/paper_plots/analysis_plots/pie_graph_q3a.png")
plt.close('all')

#####################################################################

### LINE SLOPE VS STELLAR MASS

# UNFINISHED

#####################################################################


