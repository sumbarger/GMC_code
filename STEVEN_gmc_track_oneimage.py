import sys
import matplotlib.pyplot as plt
import numpy as np
import fof_analysis
import pdb
import gizmo_analysis as gizmo
import utilities as ut
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import pickle
import seaborn as sns
import pandas as pd
from astropy.io import ascii
from astropy.table import Table

MsunToGm = 1.99e33 #conversion from mass of sun to the mass times the constant G
KpcToCm = 3.086e21 #conversion from kiloparsecs to centimeters
mp = 1.67e-24 #???
#bin_edge = 10.
bin_edge = 20. #the bin edge value
##gas properties to carry around

xg = [] #x position of gas
yg = [] #y position of gas
zg = [] #z position of gas
vxg = [] #x velocity of gas
vyg = [] #y velocity of gas
vzg = [] #z velocity of gas
mg = [] #mass of gas
rhog = [] #density of gas
tg = [] #temperature of gas
idg = [] #id number of gas

xs = [] #x position of star
ys = [] #y position of star
zs = [] #z position of star
ms = [] #mass of the star
ids = [] #id number of star
    
#star properties to carry around

sxp = [] #star x position?
syp = [] #star y position?
szp = [] #star z position?
sxp_new = [] #new star x position?
syp_new = [] #new star y position?
szp_new = [] #new star z position?
sage = [] #age of the star

gxp = [] #gas x position?
gyp = [] #gas y position?
gzp = [] #gas z position?

vdispx = [] #velocity at the x position?
vdispy = [] #velocity at the y position?
vdispz = [] #veloctiy at the z position?

snaps = np.arange(590,593,1) #set snaps to an arrange of values
cl_id = pickle.load(open('m12m_591_clusters.pkl','rb')) 
#loads a pickle file of m12m cluster and sets it to a variable

sxcm = np.zeros(shape=(len(snaps),len(cl_id))) #star x position in cm
sycm = np.zeros(shape=(len(snaps),len(cl_id))) #star y position in cm
szcm = np.zeros(shape=(len(snaps),len(cl_id))) #star z position in cm
ap_num = np.zeros(shape=(len(snaps),len(cl_id))) #???

r_ap = 0.05 #???

point_gas = [] #gas point

#snap time
tsnap = []
############################################################################
#read in sim files and find relevant particles
############################################################################

for i in range(len(snaps)): #if i is in the range of the length of snaps
    #/home/smbeninc/scr/m12m/m12m_res7100/ #where it is located
    #part = gizmo.io.Read.read_snapshots('all', 'snapshot_index', snaps[i],assign_host_principal_axes=True)
    #read the snapshot

    part = gizmo.io.Read.read_snapshots('all', 'snapshot_index', snaps[i],'../gmc/m12m/',assign_host_principal_axes=True) #read snapshot of m12m and assingn principal axes

    tsnap.append(np.max(part['star'].prop('form.time'))) #add the max star form time to the snap times

    ig = np.where((part['gas'].prop('host.distance.principal.cylindrical')[:,0] <= bin_edge) & (np.fabs(part['gas'].prop('host.distance.principal.cartesian')[:,2]) <= 1.5) & (part['gas']['temperature'] <= 1e4))
    #sets a variable where the gas cylindrical distance is less than or equal to the bin edge and 
    #the cartesian distance is less than or equal to 1.5 and the gas temperature is less than 10000
    
    ist = np.where((part['star'].prop('host.distance.principal.cylindrical')[:,0] <= bin_edge) & (np.fabs(part['star'].prop('host.distance.principal.cartesian')[:,2]) <= 1.5))
    #sets a variable where the star cylindrical distance is less than or equal to the bin edge and
    #the cartesian distance is less than or equal to 1.5
    
    xg.append(part['gas'].prop('host.distance.principal.cartesian')[ig[0],0]) 
    #add the cartesian x coordinate of ig to the x values of the gas
    yg.append(part['gas'].prop('host.distance.principal.cartesian')[ig[0],1])
    #add the cartesian y coordinate of ig to the y values of the gas
    zg.append(part['gas'].prop('host.distance.principal.cartesian')[ig[0],2])
    #add the cartesian z coordinate of ig to the z values of the gas
    
    vxg.append(part['gas'].prop('host.velocity.principal.cartesian')[ig[0],0])
    #add the x velocity of ig to the x velocity values of the gas
    vyg.append(part['gas'].prop('host.velocity.principal.cartesian')[ig[0],1])
    #add the y velocity of ig to the y velocity values of the gas
    vzg.append(part['gas'].prop('host.velocity.principal.cartesian')[ig[0],2])
    #add the z velocity of ig to the z velocity values of the gas
    
    mg.append(part['gas']['mass'][ig])
    #add the mass of ig to the mass values of the gas
    rhog.append(part['gas'].prop('number.density')[ig])
    #add the density of ig to the density values of the gas
    tg.append(part['gas']['temperature'][ig])
    #add the temperature of ig to the temperature values of the gas
    idg.append(part['gas']['id'][ig])
    #add the id of ig to the id values of the gas
    
    xs.append(part['star'].prop('host.distance.principal.cartesian')[ist[0],0])
    #add the cartesian x coordinate of ist to the x values of the stars
    ys.append(part['star'].prop('host.distance.principal.cartesian')[ist[0],1])
    #add the cartesian y coordinate of ist to the y values of the stars
    zs.append(part['star'].prop('host.distance.principal.cartesian')[ist[0],2])
    #add the cartesian z coordinate of ist to the z values of the stars
    ms.append(part['star']['mass'][ist])
    #add the mass of ist to the mass values of the stars
    ids.append(part['star']['id'][ist])
    #add the id of ist to the id values of the stars
    sage.append(part['star'].prop('age')[ist]*1000.)
    #add the age of ist to the age values of the stars
        
    msh = np.array(ms[i])
    #assign variable msh to be the array of the mass of stars of i
    xsh = np.array(xs[i])
    #assign variable xsh to be the array of the x coordinates of the stars of i
    ysh = np.array(ys[i])
    #assign variable ysh to be the array of the y coordinates of the stars of i
    zsh = np.array(zs[i])
    #assign variable zsh to be the array of the z coordinates of the stars of i

    mgh = np.array(mg[i])
    #assign variable mgh to be the array of the mass of the gas of i
    xgh = np.array(xg[i])
    #assign variable xgh to be the array of the x coordinates of the gas of i
    ygh = np.array(yg[i])
    #assign variable ygh to be the array of the y coordinates of the gas of i
    zgh = np.array(zg[i])
    #assign variable zgh to be the array of the z coordinates of the gas of i

    vxgh = np.array(vxg[i])
    #assign variable vxgh to be the array of the x velocity of the gas of i
    vygh = np.array(vyg[i])
    #assign variable vygh to be the array of the y velocity of the gas of i
    vzgh = np.array(vzg[i])
    #assign variable vzgh to be the array of the z velocity of the gas of i
    
    for j in range(len(cl_id)):
        #when it is in the range of the length of the pickle file we opened
        cl = cl_id[j] #assign variable cl to be the points in the range of the pickle file m12m
        ind = [] #???
        check = [] #???
        print(j, ind, check) #print the j, ind, and check
        for k in range(len(cl)):
            #when it is in the range of the length of cl
            ind.append(np.where(cl[k] == ids[i])[0]) 
            #add the value of cl of k is equal to the first value of the id for stars to ind
            check.append(len(ind[k]))
            #add the length of ind of k to the check
        ind = np.array(ind)
        #set variable ind to be an array of ind
        ind2 = np.array(ind[np.where(np.array(check) == 1)[0]])
        #set variable ind2 to be an array where the array check is equal to 1 or if the check is true
        ind2 = ind2.astype(dtype='int') 
        #convert array ind2 to be a int type of array
        sxcm[i,j] = np.sum(msh[ind2]*xsh[ind2])/np.sum(msh[ind2])
        #sets up a variable as a math equation to get the star x position in centimeters
        sycm[i,j] = np.sum(msh[ind2]*ysh[ind2])/np.sum(msh[ind2])
        #sets up a variable as a math equation to get the star y position in centimeters
        szcm[i,j] = np.sum(msh[ind2]*zsh[ind2])/np.sum(msh[ind2])
        #sets up a variable as a math equation to get the star z position in centimeters
        sxp.append(xsh[ind2])
        #add the ind2 of the xsh to the star x position
        syp.append(ysh[ind2])
        #add the ind2 of the ysh to the star y position
        szp.append(zsh[ind2])
        #add the ind2 of the zsh to the star z position

        rh = np.sqrt((xgh-sxcm[i,j])**2 + (ygh-sycm[i,j])**2 + (zgh-szcm[i,j])**2)
        #calculate the density of the gas using the x, y, and z of the gas
        rsh = np.sqrt((xsh-sxcm[i,j])**2 + (ysh-sycm[i,j])**2 + (zsh-szcm[i,j])**2)
        #calculate the density of the gas using the x, y, and z of the stars
        ind = np.where(rh <= r_ap)
        #assigns variable ind to be the values where the density is less than r_ap
        n = np.float(len(ind[0]))
        #assigns variable n to be the floating points of the length of the first value of the array in ind
        indst = np.where((rsh <= r_ap) & (sage[i] <= np.float(i)*2.2))
        #assigns a variable where rsh is less than r_ap and the age is less than the floating points
        gxp.append(xgh[ind])
        #add the ind of the xgh to the gxp values
        gyp.append(ygh[ind])
        #add the ind of the ygh of the gyp values
        sxp_new.append(xsh[indst])
        #add the indst of the xsh to sxp_new
        syp_new.append(ysh[indst])
        #add the indst of the ysh to syp_new
        if (len(ind[0]) > 0): #if statement if it is greater than 0
            vmean = np.average(vxgh[ind],weights=mgh[ind])
            #calculate the average and the weights
            vdispx.append(np.sqrt(np.sum(mgh[ind]*(vxgh[ind]-vmean)**2)/ ((n-1)*np.sum(mgh[ind])/n)))
            #add it to vdispx
            vmean = np.average(vygh[ind],weights=mgh[ind])
            #calculate the average and the weights
            vdispy.append(np.sqrt(np.sum(mgh[ind]*(vygh[ind]-vmean)**2)/ ((n-1)*np.sum(mgh[ind])/n)))
            #add it to vdispy
            vmean = np.average(vzgh[ind],weights=mgh[ind])
            #calculate the average and the weights
            vdispz.append(np.sqrt(np.sum(mgh[ind]*(vzgh[ind]-vmean)**2)/ ((n-1)*np.sum(mgh[ind])/n)))
            #add it to vdispz
        else: #if it is not greater than zero
            vdispx.append(0.) #add it to vdispx
            vdispy.append(0.) #add it to vdispy
            vdispz.append(0.) #add it to vdispz

            

######compare to clouds & track clouds
#which cluster at 591 to focus on
#focus_list = [4, 9, 13, 16, 17, 27, 53, 63, 65, 67, 70, 74, 82, 83]
#focus = 83

#focuses on different things
#focus = 13
#focus_list = [13,13]
focus_list = np.arange(1,2,1) #makes an arrange of numbers of 1 to 2 going up by 1
#focus_list = np.arange(len(cl_id))
#focus_list = [1,38]

alpha_track = np.zeros(shape=(3,len(cl_id))) #creates an array of zeros with the length of cl_id
alpha_s_track = np.zeros(shape=(3,len(cl_id))) #creates an array of zeros with the length of cl_id
vdisp_track = np.zeros(shape=(3,len(cl_id))) #creates an array of zeros with the length of cl_id
mass_track = np.zeros(shape=(3,len(cl_id))) #creates an array of zeros with the length of cl_id
xcm_track = np.zeros(shape=(3,len(cl_id))) #creates an array of zeros with the length of cl_id
ycm_track = np.zeros(shape=(3,len(cl_id))) #creates an array of zeros with the length of cl_id
zcm_track = np.zeros(shape=(3,len(cl_id))) #creates an array of zeros with the length of cl_id

data0 = ascii.read('cloud_props_m12m_590.txt') #read the text file
idata0 = ascii.read('cloud_indices_m12m_590.txt') #read the text file

data = ascii.read('cloud_props_m12m_591.txt') #read the text file
idata = ascii.read('cloud_indices_m12m_591.txt') #read the text file

data2 = ascii.read('cloud_props_m12m_592.txt') #read the text file
idata2 = ascii.read('cloud_indices_m12m_592.txt') #read the text file

for focus in focus_list:
    #focus on the ones in the focus_list
    #load in the cloud catalogs
    ##########find the original cloud
    xh = sxcm[1,focus] #set xh to be the focus of the first array of sxcm
    yh = sycm[1,focus] #set yh to be the focus of the first array of sycm
    zh = szcm[1,focus] #set zh to be the focus of the first array of szcm
    
    distx = np.array(data['xcm']) - xh 
    #set distx to be the array where xcm minus xh
    disty = np.array(data['ycm']) - yh
    #set disty to be the array where ycm minus yh
    distz = np.array(data['zcm']) - zh
    #set distz to be the array where zcm minus zh
    dist = np.sqrt(distx**2 + disty**2 + distz**2)
    #calculate the distance
    focus_cl = np.where(dist == min(dist))[0]
    #set focus_cl to be where distance is the minimum
    focus_cl_m = data['mass'][np.where(data['cloud number'] == focus_cl)]
    #set focus_cl_m to be the mass where the the cloud number is equal to the focus cloud
    
    cluster_gas = idata['id'][np.where(idata['cloud number'] == focus_cl)[0]]
    #set cluster_gas to the id where the cloud number is equal to the first number of the array focus_cl
    
    #use Andrew's particle tracking...eventually
    #find the original GMC
    backward = [] #use tracking
    for i in range(len(cluster_gas)):
        if (len(np.where(idata0['id'] == cluster_gas[i])[0]) == 2):
            #if it is equal to the cluster gas id number
            dists = np.sqrt((idata0['x'][np.where(idata0['id'] == cluster_gas[i])[0]]-xh)**2 + (idata0['y'][np.where(idata0['id'] == cluster_gas[i])[0]]-yh)**2 + (idata0['z'][np.where(idata0['id'] == cluster_gas[i])[0]]-zh)**2) #calculate the distance
            #print(idata0['id'][np.where(idata0['id'] == cluster_gas[i])[0][np.where(dists == np.min(dists))]])
            #print the id of the clouds
            backward.append(idata0['cloud number'][np.where(idata0['id'] == cluster_gas[i])[0][np.where(dists == np.min(dists))]])
            #add the cloud number to backward tracking
            #print(idata0['x'][np.where(idata0['id'] == cluster_gas[i])[0][np.where(dists == np.min(dists))]])
            #print the data of the cluster gas
        if (len(np.where(idata0['id'] == cluster_gas[i])[0]) == 1):
            #if the length of id is equal to the cluster gas
            #print(idata0['x'][np.where(idata0['id'] == cluster_gas[i])[0]])
            #print the data
            backward.append(idata0['cloud number'][np.where(idata0['id'] == cluster_gas[i])[0]])
            #add it to the backward tracking
            
    original_cl = np.int(np.median(backward))
    #set the variable to the median of the backward tracking data
    
    #find out what happened to it
    xh = sxcm[2,focus] #look at the x value of the position
    yh = sycm[2,focus] #look at the y value of the position
    zh = szcm[2,focus] #look at the z value of the position

    evolve = [] #look at the evolution
    evolve_mass = [] #look at the mass of the evolution
    
    for i in range(len(cluster_gas)):
        #i in the range of the length of the cluster_gas
        if (len(np.where(idata2['id'] == cluster_gas[i])[0]) == 1):
            #if the length of of the id is equal to the cluster_gas
            evolve.append(idata2['cloud number'][np.where(idata2['id'] == cluster_gas[i])[0]])
            evolve_mass.append(np.float(idata2['mass'][np.where(idata2['id'] == cluster_gas[i])[0]]))
        if (len(np.where(idata2['id'] == cluster_gas[i])[0]) == 2):
            dists = np.sqrt((idata2['x'][np.where(idata2['id'] == cluster_gas[i])[0]]-xh)**2 + (idata2['y'][np.where(idata2['id'] == cluster_gas[i])[0]]-yh)**2 + (idata2['z'][np.where(idata2['id'] == cluster_gas[i])[0]]-zh)**2)
            evolve.append(idata2['cloud number'][np.where(idata2['id'] == cluster_gas[i])[0][np.where(dists == np.min(dists))]])
            evolve_mass.append(np.float(idata2['mass'][np.where(idata2['id'] == cluster_gas[i])[0][np.where(dists == np.min(dists))]]))

    #look at the whole cloud
        
    evolve = np.array(evolve)
    evolve_mass = np.array(evolve_mass)
    unq = np.unique(evolve)
    unq_mass = np.zeros(len(unq))
    for i in range(len(unq)):
        k = np.where(evolve == unq[i])[0]
        unq_mass[i] = np.sum(evolve_mass[k])

    if (len(unq_mass) > 0):
        match = np.where(unq_mass == np.max(unq_mass))
        evolve_cl = unq[match]
        evolve_cl_m = unq_mass[match][0]
        evolve_cl_mtot = data2['mass'][np.where(data2['cloud number'] == evolve_cl)]
        ratio = evolve_cl_m/focus_cl_m
    else:
        ratio = 0.

    mass_track[0,focus] = data0['mass'][np.where(data0['cloud number'] == original_cl)]
    mass_track[1,focus] = data['mass'][np.where(data['cloud number'] == focus_cl)]
    alpha_track[0,focus] = data0['alpha_vir'][np.where(data0['cloud number'] == original_cl)]
    alpha_track[1,focus] = data['alpha_vir'][np.where(data['cloud number'] == focus_cl)]
    #alpha_s_track[0,focus] = data0['alpha_vir_s'][np.where(data0['cloud number'] == original_cl)]
    #alpha_s_track[1,focus] = data['alpha_vir_s'][np.where(data['cloud number'] == focus_cl)]
    vdisp_track[0,focus] = np.sqrt(data0['vdisp_x'][np.where(data0['cloud number'] == original_cl)]**2 + data0['vdisp_y'][np.where(data0['cloud number'] == original_cl)]**2 + data0['vdisp_z'][np.where(data0['cloud number'] == original_cl)]**2)/np.sqrt(3.)
    vdisp_track[1,focus] = np.sqrt(data['vdisp_x'][np.where(data['cloud number'] == focus_cl)]**2 + data['vdisp_y'][np.where(data['cloud number'] == focus_cl)]**2 + data['vdisp_z'][np.where(data['cloud number'] == focus_cl)]**2)/np.sqrt(3.)
    xcm_track[1,focus] = data['xcm'][np.where(data['cloud number'] == focus_cl)]
    ycm_track[1,focus] = data['ycm'][np.where(data['cloud number'] == focus_cl)]
    zcm_track[1,focus] = data['zcm'][np.where(data['cloud number'] == focus_cl)]
    xcm_track[0,focus] = data0['xcm'][np.where(data0['cloud number'] == original_cl)]
    ycm_track[0,focus] = data0['ycm'][np.where(data0['cloud number'] == original_cl)]
    zcm_track[0,focus] = data0['zcm'][np.where(data0['cloud number'] == original_cl)]

        
    if ((ratio < 0.5) | (evolve_cl_m/evolve_cl_mtot < 0.5)):
        print(focus,evolve_cl,np.float(ratio),"DEAD",np.float(evolve_cl_m/evolve_cl_mtot))
        alpha_track[2,focus] = 0.
        vdisp_track[2,focus] = 0.
        mass_track[2,focus] = 0.
        xcm_track[2,focus] = sxcm[2,focus]
        ycm_track[2,focus] = sycm[2,focus]
        zcm_track[2,focus] = szcm[2,focus]
    elif ((ratio >= 0.5) & (evolve_cl_m/evolve_cl_mtot >= 0.5)):
        print(focus, evolve_cl,np.float(ratio),'ALIVE')
        mass_track[2,focus] = data2['mass'][np.where(data2['cloud number'] == evolve_cl)]
        alpha_track[2,focus] = data2['alpha_vir'][np.where(data2['cloud number'] == evolve_cl)]
        #alpha_s_track[2,focus] = data2['alpha_vir_s'][np.where(data2['cloud number'] == evolve_cl)]
        vdisp_track[2,focus] = np.sqrt(data2['vdisp_x'][np.where(data2['cloud number'] == evolve_cl)]**2 + data2['vdisp_y'][np.where(data2['cloud number'] == evolve_cl)]**2 + data2['vdisp_z'][np.where(data2['cloud number'] == evolve_cl)]**2)/np.sqrt(3.)
        xcm_track[2,focus] = data2['xcm'][np.where(data2['cloud number'] == evolve_cl)]
        ycm_track[2,focus] = data2['ycm'][np.where(data2['cloud number'] == evolve_cl)]
        zcm_track[2,focus] = data2['zcm'][np.where(data2['cloud number'] == evolve_cl)]

print(sxcm[:,38])
print(sycm[:,38])
print(szcm[:,38])
################################################################################
# try to make some nicer images
################################################################################

#for 38
#[-3.70781625 -3.91423683 -3.80127101]
#[-2.62598792 -0.78919105 -1.14497895]
#[-0.37837485 -0.02532101 -0.03734116]


#image some stuff
bin_edge = 0.5
sizebox = bin_edge*2.
stepn = 30

count = 0
htot = 0.
w2d_kernel = []
nkernel = 200
w2d_table = np.zeros(nkernel+3)

vel_array0 = np.zeros(shape=(2,stepn, stepn))
vel_array1 = np.zeros(shape=(2,stepn, stepn))
vel_array2 = np.zeros(shape=(2,stepn, stepn))

#read the kernel file
kernelfile = 'kernel2d'
temp = open(kernelfile,'r').read().split('\n')

for i in range(len(temp)-1):
    if (i%2 == 0):
        w2d_kernel.append(np.float(temp[i]))

w2d_kernel.append(0.)
w2d_kernel = np.array(w2d_kernel)

stepnh= (stepn-1)/2
res= sizebox/stepn
print('stepnh = ', stepnh)
print('res = ', res)

#focus = 38
focus = 1

for t in range(len(snaps)):
    w2darray = np.zeros(shape=(stepn, stepn))
    vel_array = np.zeros(shape=(stepn, stepn))
    xcen = xcm_track[t,focus]
    ycen = ycm_track[t,focus]
    zcen = zcm_track[t,focus]
    
    #read in sim snap
    part = gizmo.io.Read.read_snapshots('all', 'snapshot_index', snaps[t],'/home/smbeninc/scr/m12m/m12m_res7100/',assign_host_principal_axes=True)
    #part = gizmo.io.Read.read_snapshots('all', 'snapshot_index', snaps[t],assign_host_principal_axes=True)
    
    ig = np.where((part['gas'].prop('host.distance.principal.cartesian')[:,0] <= (xcen + bin_edge)) & (part['gas'].prop('host.distance.principal.cartesian')[:,0] > (xcen - bin_edge)) &  (part['gas'].prop('host.distance.principal.cartesian')[:,1] > (ycen - bin_edge)) & (part['gas'].prop('host.distance.principal.cartesian')[:,1] <= (ycen + bin_edge)) & (np.fabs(part['gas'].prop('host.distance.principal.cartesian')[:,2]) <= 1.5))
    ngas = len(ig[0])
    
    xg = part['gas'].prop('host.distance.principal.cartesian')[ig[0],0]
    yg = part['gas'].prop('host.distance.principal.cartesian')[ig[0],1]
    zg = part['gas'].prop('host.distance.principal.cartesian')[ig[0],2]
    mg = part['gas']['mass'][ig]
    hg = part['gas']['smooth.length'][ig]
    
    #set the kernel weights, make sure to normalize properly
    for i in range(ngas):
        htot = htot + hg[i]/ngas
        
        
    ##big loop
    for i in range(ngas):
        xi = (xg[i]-xcen)/res+stepnh
        yi = (yg[i]-ycen)/res+stepnh
        hi = 2*hg[i]/res
        
        ixmin= np.int(np.floor(xi-hi))+1
        if (ixmin < 0):
            ixmin=0
        if (ixmin > stepn-1):
            ixmin=stepn-1
        ixmax= np.int(np.floor(xi+hi))
        if (ixmax > stepn-1):
            ixmax=stepn-1
        if (ixmax < 0):
            ixmax=0
        iymin= np.int(np.floor(yi-hi))+1
        if (iymin < 0):
            iymin=0
        if (iymin > stepn-1):
            iymin=stepn-1
        iymax= np.int(np.floor(yi+hi))
        if (iymax > stepn-1):
            iymax=stepn-1
        if (iymax < 0):
            iymax=0
            
        if (i%1000 == 0):
            print("Doing particle: ",i,ngas)
            
        count = count + 1
        w2d_sum = 0.
        d2max = 4.*hg[i]*hg[i]
        kernfrac = nkernel*.5/hg[i]
        loop = np.arange(ixmin,ixmax+1,1)
        for ii in range(len(loop)):
            r0xx = (ii - stepnh)*res + xcen
            for jj in range(len(loop)):
                r0yy = (jj-stepnh)*res + ycen
                xx = xg[i] - r0xx
                yy = yg[i] - r0yy
                d2 = xx*xx + yy*yy
                if (d2 <= d2max):
                    d = np.sqrt(d2)
                    xii = d*kernfrac
                    xxi = np.int(np.floor(xii))
                    w2d= ((w2d_table[xxi+1] - w2d_table[xxi])*(xii - xxi)+ w2d_table[xxi])
                    w2darray[ii,jj] = w2d
                    w2d_sum = w2d_sum + w2d
                else:
                    w2darray[ii,jj] = 0.
        if (w2d_sum > 0):
            weight = mg[i]/w2dsum
            for ii in range(len(loop)):
                for jj in range(len(loop)):
                    if (w2darray[ii,jj] > 0):
                        vel_array[ii,jj] = vel_array[ii,jj] + w2darray[ii,jj]*weight
        else:
            ii = np.int(np.floor(xi + 0.5))
            jj = np.int(np.floor(yi + 0.5))
            if ((ii >=0) & (ii <= stepn-1) & (jj >= 0) & (jj <= stepn-1)):
                vel_array[ii,jj] = vel_array[ii,jj] + mg[i]

    if (t == 0):
        vel_array0[0,:,:] = vel_array
    if (t == 1):
        vel_array1[0,:,:] = vel_array
    if (t == 2):
        vel_array2[0,:,:] = vel_array
   
        

fig = plt.figure(figsize=(8,2.5))

gs0 = gridspec.GridSpec(1, 4, width_ratios=[1,1,1,0.1])
#gs00 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[0],wspace=0.05)
#gs01 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[1],wspace=0.05)

ax0 = plt.subplot(gs0[0])
ax1 = plt.subplot(gs0[1])
ax2 = plt.subplot(gs0[2])
axcb = plt.subplot(gs0[3])

ax0.imshow(vel_array0[0,:,:]/((res*1000.)**2), interpolation='gaussian',norm=colors.LogNorm(),vmin=1,vmax=1000,extent=(xcm_track[0,1]-0.5,xcm_track[0,1]+0.5,ycm_track[0,1]-0.5,ycm_track[0,1]+0.5))
ax1.imshow(vel_array1[0,:,:]/((res*1000.)**2), interpolation='gaussian',norm=colors.LogNorm(),vmin=1,vmax=1000,extent=(xcm_track[1,1]-0.5,xcm_track[1,1]+0.5,ycm_track[1,1]-0.5,ycm_track[1,1]+0.5))
ax2.imshow(vel_array2[0,:,:]/((res*1000.)**2), interpolation='gaussian',norm=colors.LogNorm(),vmin=1,vmax=1000,extent=(xcm_track[2,1]-0.5,xcm_track[2,1]+0.5,ycm_track[2,1]-0.5,ycm_track[2,1]+0.5))

ax0.plot([5.65,5.65+0.25],[-7.05,-7.05],color='k',linewidth=1.5)
ax0.text(5.625,-7,'250 pc')

#plot the stars
ax1.scatter(sxp[1+len(cl_id)], syp[1+len(cl_id)], s=8,marker='*', edgecolors='k',facecolors='none',linewidth=0.5)
ax1.scatter(sxp[37+len(cl_id)], syp[37+len(cl_id)], s=8,marker='*',edgecolors='k',facecolors='none',linewidth=0.5)
ax1.scatter(sxp[54+len(cl_id)], syp[54+len(cl_id)], s=8,marker='*', edgecolors='k',facecolors='none',linewidth=0.5)
ax1.scatter(sxp[73+len(cl_id)], syp[73+len(cl_id)], s=8,marker='*', edgecolors='k',facecolors='none',linewidth=0.5)

ax2.scatter(sxp[1+2*len(cl_id)], syp[1+2*len(cl_id)], s=8,marker='*', edgecolors='k',facecolors='none',linewidth=0.5)
ax2.scatter(sxp[37+2*len(cl_id)], syp[37+2*len(cl_id)], s=8,marker='*', edgecolors='k',facecolors='none',linewidth=0.5)
ax2.scatter(sxp[54+2*len(cl_id)], syp[54+2*len(cl_id)], s=8,marker='*', edgecolors='k',facecolors='none',linewidth=0.5)
ax2.scatter(sxp[73+2*len(cl_id)], syp[73+2*len(cl_id)], s=8,marker='*', edgecolors='k',facecolors='none',linewidth=0.5)




#ax11.scatter(sxcm[1,1],sycm[1,1],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)
#ax11.scatter(sxcm[1,37],sycm[1,37],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)
#ax11.scatter(sxcm[1,54],sycm[1,54],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)
#ax11.scatter(sxcm[1,73],sycm[1,73],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)

#ax12.scatter(sxcm[2,1],sycm[2,1],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)
#ax12.scatter(sxcm[2,37],sycm[2,37],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)
#ax12.scatter(sxcm[2,54],sycm[2,54],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)
#ax12.scatter(sxcm[2,73],sycm[2,73],marker='o',s=80,edgecolors='k',facecolors='none',linewidth=0.75)

ax0.get_xaxis().set_visible(False)
ax0.get_yaxis().set_visible(False)
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

norm1 = matplotlib.colors.LogNorm(vmin=1,vmax=1000)
cb = matplotlib.colorbar.ColorbarBase(axcb,norm=norm1)
cb.set_label('$\Sigma$ (M$_{\odot}$/pc$^2$)')

ax0.set_title(r'$t_0$')
ax1.set_title(r'$t_1$')
ax2.set_title(r'$t_2$')

plt.savefig('andrew_gmc_cluster.pdf',bbox_inches='tight') 
plt.close(fig)
