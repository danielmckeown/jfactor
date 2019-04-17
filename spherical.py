import numpy as np
import h5py
import scipy.stats as st
from astropy import constants as const
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from pylab import *
import math
from operator import truediv
######################                                                                                                                                                 

def get_data(filename,num_of_file,key1,key2):
    if num_of_file == 1:
        f = h5py.File(filename+'.hdf5', 'r')
        if key1 == 'Header':
            return f[key1].attrs[key2]
        else:
            return f[key1][key2][:]
    else:
        for i in range(0,num_of_file):
            f = h5py.File(filename+'.'+str(i)+'.hdf5', 'r')
            if key1 == 'Header':
                return f[key1].attrs[key2]
            else:
                if ( len(f[key1][key2][:].shape)==1 ):
                    if i==0:
                        result = f[key1][key2][:]
                    else:
                        result = np.hstack( (result,f[key1][key2][:]) )
                else:
                    if i==0:
                        result = f[key1][key2][:]
                    else:
                        result = np.vstack( (result,f[key1][key2][:]) )

        return result


######## ROTATIONS                                                                                                                                                     

def rotateX(lx, ly, lz):
    """ Rotates the point around the X axis by the given angle in degrees. """
    rad = 0.2
    cosa = math.cos(rad)
    sina = math.sin(rad)
    y = ly * cosa - lz * sina
    z = ly * sina + lz * cosa
    return lx, y, z

def rotateY(mx, my, mz):
    """ Rotates the point around the Y axis by the given angle in degrees. """
    rad = 0.2
    cosa = math.cos(rad)
    sina = math.sin(rad)
    z = mz * cosa - mx * sina
    x = mz * sina + mx * cosa
    return x, my, z



def rotateZ(nx, ny, nz):
    """ Rotates the point around the Z axis by the given angle in degrees. """
    rad = 0.02
    cosa = math.cos(rad)
    sina = math.sin(rad)
    x = nx * cosa - ny * sina
    y = nx * sina + ny * cosa
    return x, y, nz



def rotate_X_back(lx, ly, lz):
    """ Rotates the point around the X axis by the given angle in degrees. """
    rad = - 0.2 * index2
    cosa = math.cos(rad)
    sina = math.sin(rad)
    y = ly * cosa - lz * sina
    z = ly * sina + lz * cosa
    return lx, y, z


def rotate_Y_back(lx, ly, lz):
    """ Rotates the point around the X axis by the given angle in degrees. """
    rad = - 0.2 * index1
    cosa = math.cos(rad)
    sina = math.sin(rad)
    y = ly * cosa - lz * sina
    z = ly * sina + lz * cosa
    return lx, y, z




def rotate_Y_negative(lx, ly, lz):
    """ Rotates the point around the X axis by the given angle in degrees. """
    rad = - 0.2 * index1
    cosa = math.cos(rad)
    sina = math.sin(rad)
    y = ly * cosa - lz * sina
    z = ly * sina + lz * cosa
    return lx, y, z



def rotate_X_back_positive(lx, ly, lz):
    """ Rotates the point around the X axis by the given angle in degrees. """
    rad =  0.2 * index2
    cosa = math.cos(rad)
    sina = math.sin(rad)
    y = ly * cosa - lz * sina
    z = ly * sina + lz * cosa
    return lx, y, z



############################################### SWITCHES ##########################################################                                                   \
                                                                                                                                                                       
                                                                                                                                                                      \

#... Elvis (0) or LATTE (1)                                                                                                                                           \
                                                                                                                                                                       
                                                                                                                                                                      \

sim = 1

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):    

#... Elvis (0) or LATTE (1)                                                                                                                                           \
                                                                                                                                                                       
                                                                                                                                                                      \

sim = 1

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):                                                                                                                \
                                                                                                                                                                       
                                                                                                                                                                      \

itype = 1

#...Snapshot                                                                                                                                                          \
                                                                                                                                                                       
                                                                                                                                                                      \

isnap = 600

#... What run?                                                                                                                                                        \
                                                                                                                                                                       
                                                                                                                                                                      \

irun = 'dark'
#...Number of bins in histogram                                                                                                                                       \
                                                                                                                                                                       
                                                                                                                                                                      \

my_bins = 150
#...Out to radius (kpc) (KEEP AS FLOAT)


###############################################  BOOKKEEPING  #####################################################                                                   \
                                                                                                                                                                       

part_list = ['PartType0','PartType1','PartType2','PartType4']
p_type = part_list[itype]

type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
type_tag = type_list[itype]

#...Path to data                                                                                                                                                      \
                                                                                                                                                                       


if sim == 0:
    base = '/data11/home/dmckeown/output/m12i/ELVIS/'
    base1 = '/data11/home/dmckeown/output/m12i/ELVIS/'
if sim == 1:
    base = '/data18/brunello/dmckeown/output'
    base1 = '/data18/brunello/dmckeown/plots/'
fname = base + '/' + irun + '/snapshot_'+ str(isnap)
fname1 = base + '/' + irun + '/halo_' + str(isnap)

f = h5py.File(fname1 + '.hdf5', 'r')


print f.keys()
print "these are the keys"                                                                                                                                            \




h = get_data(fname, 1,'Header','HubbleParam')
print 'h: '+ str(h)
dm_xyz = get_data(fname, 1, p_type, 'Coordinates')/h
velocities  = get_data(fname, 1, p_type, 'Velocities')
print "printing len of velocities"
print len(velocities)


radius = np.sqrt(dm_xyz[:,0]*dm_xyz[:,0] + dm_xyz[:,1]*dm_xyz[:,1] + dm_xyz[:,2]*dm_xyz[:,2])

dm_mass = get_data(fname, 1, p_type, 'Masses')*(10**10)/h
print 'DM data loaded...'
print "now printing dm mass"
print dm_mass
f = h5py.File(fname1 + '.hdf5', 'r')
halo_pos = np.array(f['position'])#/h#/(10**3)                                                                                                                         


virial = np.array(f['mass.vir'])
host_velocity = np.array(f['host.velocity'])

import operator
index, value = max(enumerate(virial), key=operator.itemgetter(1))
print "printing index and value"
print index,value
virial_radius = np.array(f['radius'])
print "printing virial radius"
vir =  virial_radius[index]
print vir
host_v = host_velocity[index]
print "printing all velocities"
print host_velocity
host_cent = halo_pos[index]

my_rad = vir
print "printing host vel"
print host_v

dm_xyz = dm_xyz - host_cent
thetaa1 = 0.2
phii1 = 0.2
###### Defining the coordinates                                                                                                                                        

index1 = 0
counter = 0
F_J = []
Th = []
Ph = []
print "printing first double check"
print dm_xyz[:,1][0]
while index1 < 4:



    index2 = 0
    while index2 < 4 :
        r = np.sqrt(( dm_xyz[:,0] )*(  dm_xyz[:,0] ) + (dm_xyz[:,1] )* (dm_xyz[:,1] ) + (dm_xyz[:,2] + 8.0 )*(dm_xyz[:,2] + 8.0 ) )

        Thetta = np.arctan2(( dm_xyz[:,1] ) ,(dm_xyz[:,0]  ))

        u = (dm_xyz[:,2] + 8.0)  / ( r )

        Phi = np.arccos(u)


        index_rad = ( r < 200.0 ) & (Phi < 0.2 )
        dm_xyz1 = dm_xyz[index_rad]

        print "printing len r after coords"
        print len(dm_xyz1)

        radius = np.sqrt(dm_xyz1[:,0]*dm_xyz1[:,0] + dm_xyz1[:,1]*dm_xyz1[:,1] + dm_xyz1[:,2]*dm_xyz1[:,2])
        n_bins = 10 #Number of bins you want                                                                                                                           

        part_per_bin = np.zeros(n_bins)

        n_bins_frac = np.linspace(1.0/n_bins, 1.0, n_bins)
#... Leave this as is. It only sets the percentage of                                                                                                                  
#... particles within each bin given your number of bins                                                                                                               

#... Binning Function                                                                                                                                                  
#... Where r is the distance from the cone vertex                                                                                                                      
#... of all the particles within the cone                                                                                                                              


        theta = [0.2]

        print "printing theta"
        print theta

        dl = []

        edges = st.mstats.mquantiles(radius,n_bins_frac)
        for j in range(len(edges)):
            if j == 0:
                mask = ( radius > 0) & ( radius < edges[j] )
            else:
                mask = (radius > edges[j - 1]) & ( radius < edges[j])
#...This array now tells you the number of particles in each bin                                                                                                       
            part_per_bin[j] = len(radius[mask])
            print "printing parts per bin"
            print part_per_bin[j]
            print edges[j]
            dl.insert(j,edges[j])


        dl.insert(0,-8.0)
	print "printing dl ########################## "
        dl.insert(11,200)
        print dl


        height = [j-i for i, j in zip(dl[:-1], dl[1:])]
        print "printing height"
        print height
        lengths = [sum(height[:i+1]) for i in range(len(height))]

        base = []
        for j in range(len(edges)):
            bases =  lengths[j] * tan(theta)
            base.insert(j,bases)

        print "printing base"
        print base

        j_factor = []
        for j in range(len(edges) ):
            if j == 0:
                volume = (1.0/3.0) *np.pi* height[j] * (base[j] * base[j] )

                mass = dm_mass[j] * part_per_bin[j]
                

                density = mass / volume
                density_squared =  np.power(density, 2)
		j_fact = density_squared[j] * height[j]

                j_factor.insert(j,j_fact)
            else:

                volume =(1.0/3.0)*np.pi* height[j]  * (  base[j] * base[j]  +  base[ j - 1] * base[j] + base[j - 1] * base[j - 1] )
		mass = dm_mass[j] * part_per_bin[j]
                density= mass / volume
                density_squared =  np.power(density, 2)
                j_fact = density_squared * height[j]

		j_factor.insert(j,j_fact)
        print "printing final j factor"

	thetaa =  thetaa1* index2
        phii =  phii1 * index1
        final_J =  np.sum(j_factor)
        print final_J
        F_J.insert(index,final_J)
	print "printing theta and phi"
        Th.insert(counter,thetaa)
	Ph.insert(counter,phii)
        print thetaa
        print phii
        print "printing y coord"
        print dm_xyz[:,1][0]
        # Rotate coordinates about x ( same as moving cone downward in the x direction)                                                                                
        dm_xyz[:,0],dm_xyz[:,1],dm_xyz[:,2] = rotateX(dm_xyz[:,0],dm_xyz[:,1],dm_xyz[:,2])
        index2 = index2 + 1
        counter = counter + 1

        print dm_xyz[:,1][0]

    dm_xyz[:,0],dm_xyz[:,1],dm_xyz[:,2] =  rotate_X_back(dm_xyz[:,0],dm_xyz[:,1],dm_xyz[:,2])

    index1 = index1 + 1
    print "double check coordinates"
    print dm_xyz[:,1][0]
    print dm_xyz[:,0][0]
    dm_xyz[:,0],dm_xyz[:,1],dm_xyz[:,2] = rotateY(dm_xyz[:,0],dm_xyz[:,1],dm_xyz[:,2])
    print "hi"
    print "printing x coord"
    print dm_xyz[:,0][0]



                                                                                                                                                                       
                                                                                                                                                                       
                                                                                                                                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
