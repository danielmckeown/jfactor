import numpy as np
import h5py
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
############################################### SWITCHES ##########################################################                                         \


#... Elvis (0) or LATTE (1)                                                                                                                                 \
                                                                                                                                                             
                                                                                                                                                            \

sim = 1

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):                                                                                                      \
                                                                                                                                                             
                                                                                                                                                            \

itype = 1

#...Snapshot                                                                                                                                                \
                                                                                                                                                             
                                                                                                                                                            \

isnap = 600

irun = 'dark'
#...Number of bins in histogram                                                                                                                             \
                                                                                                                                                             
                                                                                                                                                            \

my_bins = 150
#...Out to radius (kpc) (KEEP AS FLOAT)                                                                                                                     \
                                                                                                                                                             
                                                                                                                                                            \

###############################################  BOOKKEEPING  #####################################################                                         \
                                                                                                                                                             

part_list = ['PartType0','PartType1','PartType2','PartType4']
p_type = part_list[itype]


type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
type_tag = type_list[itype]

#...Path to data                                                                                                                                            \
                                                                                                                                                             


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
print "these are the keys"                                                                                                                                  \




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
print host_cent

my_rad = vir
print "printing host vel"
print host_v

dm_xyz = dm_xyz - host_cent
print "printing len xzy1"
radius = np.sqrt(dm_xyz[:,0]*dm_xyz[:,0] + dm_xyz[:,1]*dm_xyz[:,1] + dm_xyz[:,2]*dm_xyz[:,2])
print len(dm_xyz)

xyz = dm_xyz
#index = (radius > radii) & (radius < radii + increment)                                                                                                     
index = (  radius < 210 ) & ( - 2080 * xyz[:,0] <  0 ) &  (    2080 * xyz[:,1] < 0 )

xyz = xyz[index]
print len(xyz)

print "printing j factor length"
print len(xyz)



fig = plt.figure()
plt.scatter(xyz[:,2],xyz[:,1] )



save_fig_file = base1+str(irun)+'/dm/'+str(isnap)+'.png'
#...Report saving:                                                                                                                                           
print "Saving : "
print str(save_fig_file)
#...Save Figure:                                                                                                                                             
fig.savefig(save_fig_file)

A = [  0.0,0.0,-8.0]
B = [0.0,0.0 , 200.0 ]
C = [0.0,- 1.0, 200.0]
D = [-1.0, -1.0, 200.0]
E = [-1.0, 0.0,200.0]

A = np.array(A)
B = np.array(B)
C = np.array(C)
D = np.array(D)
E = np.array(E)
##### Now calculate the planes                                                                                                                              \
                                                                                                                                                             

AB = B - A
print "printing AB"
AC = C - A
print AB
print "printing AD"
print AC
AE = E - A
print "printing AE"
print AE
AD = D - A
BC = C - B
BE = E - B

def plane(AC ,AB,AE,AD,BC,BE):
    import numpy as np
    n1 = [ np.cross(AB, AC) ]
    n2 = [ np.cross(AE, AB) ]
    n3 = [ np.cross(AD, AE) ]
    n4 = [ np.cross(AC, AD) ]
    n5 = [ np.cross(BC, BE) ]
    # Useful for holding the place of a given vector component                                                                                               
    unit = [1,1,1]
    xyz = dm_xyz
    vector = [xyz[:,0],xyz[:,1],xyz[:,2]]
    vector = np.array(vector)
    

    eq1 = [a*b for a,b in zip(n1,unit)]
    eq11 = [a*b for a,b in zip(n1,-B)]
    n1 = np.array(n1)
    n1 = n1.T
    plane11 = []
    iplane1 = n1[0] * vector[0]
    jplane1 = n1[1] * vector[1]
    kplane1 = n1[2] * vector[2]
    if iplane1[0] == 0.0:
        plane11.insert(0, )
    else:
        plane11.insert(0,n1[0] * vector[0])
    if jplane1[0] == 0.0:
        plane11.insert(1,0.0 )
    else:
        plane11.insert(1,n1[1]* vector[1])


    eq11 = np.array(eq11)
    scalar1 = np.sum(eq11)

    eq2 = [a*b for a,b in zip(n2,unit)]
    eq22 = [a*b for a,b in zip(n2,- E)]
    eq222 = [a + b for a,b in zip(eq2,eq22)]


    eq3 = [a*b for a,b in zip(n3,unit)]
    eq33 = [a*b for a,b in zip(n3,- D)]
    eq333 = [a + b for a,b in zip(eq3,eq33)]


    eq4 = [a*b for a,b in zip(n4,unit)]
    eq44 = [a*b for a,b in zip(n4,- C)]
    eq444 = [a + b for a,b in zip(eq4,eq44)]

    eq5 = [a*b for a,b in zip(n5,unit)]
    eq55 = [a*b for a,b in zip(n5,- B)]
    eq555 = [a + b for a,b in zip(eq5,eq55)]





    index = (plane11[0] + plane11[1] + plane11[2] + scalar1 < 0 )



    return n1,n1[0],eq11,jplane1[1],plane11[0],len(plane11[0]),plane11[1],scalar1


print "printing plane"
print plane(AC,AB,AE,AD,BC,BE)



                      
