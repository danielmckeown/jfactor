
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
############################################### SWITCHES ##########################################################                                                             

#... Elvis (0) or LATTE (1)                                                                                                                                                     

sim = 1

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):                                                                                                                          

itype = 1


#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):                                                                                                                          

itype = 1

#...Snapshot                                                                                                                                                                    

isnap = 600

#... What run?                                                                                                                                                                  

irun = 'dark'
#...Number of bins in histogram                                                                                                                                                 

my_bins = 150
#...Out to radius (kpc) (KEEP AS FLOAT)                                                                                                                                         

###############################################  BOOKKEEPING  #####################################################                                                             

part_list = ['PartType0','PartType1','PartType2','PartType4']
p_type = part_list[itype]

type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
type_tag = type_list[itype]

#...Path to data                                    




if sim == 0:
    base = '/data11/home/dmckeown/output/m12i/ELVIS/'
    base1 = '/data11/home/dmckeown/output/m12i/ELVIS/'
if sim == 1:
    base = '/data18/brunello/dmckeown/output'
    base1 = '/data18/brunello/dmckeown/plots/'
fname = base + '/' + irun + '/snapshot_'+ str(isnap)
fname1 = base + '/' + irun + '/halo_' + str(isnap)

###############################################  LOAD DATA  #####################################################                                                               



def Diff(li1, li2):
    return (list(set(li1) - set(li2)))
    
    
f = h5py.File(fname1 + '.hdf5', 'r')


#print f.keys()                                                                                                                                                                 
#print "these are the keys"                                                                                                                                                     


THETTA = 2.754074388
PHII = 2.754074388
J = []
n = 0
counter = 0
while n < 10:
    m = 0

    while m < 10:
        h = get_data(fname, 1,'Header','HubbleParam')
       # print 'h: '+ str(h)                                                                                                                                                    
        dm_xyz = get_data(fname, 1, p_type, 'Coordinates')/h
        velocities  = get_data(fname, 1, p_type, 'Velocities')
       # print "printing len of velocities"                                                                                                                                     
       # print len(velocities)                                                                                                                                                  


        radius = np.sqrt(dm_xyz[:,0]*dm_xyz[:,0] + dm_xyz[:,1]*dm_xyz[:,1] + dm_xyz[:,2]*dm_xyz[:,2])

        dm_mass = get_data(fname, 1, p_type, 'Masses')*(10**10)/h
       # print 'DM data loaded...'                                                                                                                                              
       # print "now printing dm mass"                                                                                                                                           
       # print dm_mass                                                                                                                                                          
        f = h5py.File(fname1 + '.hdf5', 'r')
        halo_pos = np.array(f['position'])#/h#/(10**3)                                                                                                                          



        virial = np.array(f['mass.vir'])
        host_velocity = np.array(f['host.velocity'])

        import operator
        index, value = max(enumerate(virial), key=operator.itemgetter(1))
      #  print "printing index and value"                                                                                                                                       
      #  print index,value                                                                                                                                                      
        virial_radius = np.array(f['radius'])
      #  print "printing virial radius"                                                                                                                                         
        vir =  virial_radius[index]
      #  print vir                                                                                                                                                              
        host_v = host_velocity[index]
      #  print "printing all velocities"                                                                                                                                        
      #  print host_velocity                                                                                                                                                    
        host_cent = halo_pos[index]
       # print host_cent                                                                                                                                                        

        my_rad = vir
       # print "printing host vel"                                                                                                                                              
       # print host_v                                                                                                                                                           
        dm_xyz = dm_xyz - host_cent
       # print "printing len xzy1"                                                                                                                                              
        radius = np.sqrt(dm_xyz[:,0]*dm_xyz[:,0] + dm_xyz[:,1]*dm_xyz[:,1] + dm_xyz[:,2]*dm_xyz[:,2])
        #print len(dm_xyz)                                                                                                                                                      

        xyz = dm_xyz

        R = [0,0,-8]
        P = [0.5 + (n*1.0),0.5 + (m*1.0),200]
        Q = [0.5 + (n*1.0), -0.5 + (m*1.0),200]
        S = [-0.5 + (n*1.0), 0.5 + (m*1.0),200]
        T = [ -0.5 + (n*1.0),-0.5 + (m*1.0),200]
        Center = [(n*1.0),(m*1.0),200]
       # print "printing coordinates"                                                                                                                                           
       # print P,Q,S,T                                                                                                                                                          
        RP = [a_i - b_i for a_i, b_i in zip(P, R)]
        RQ  = [a_i - b_i for a_i, b_i in zip(Q, R)]
        RS = [a_i - b_i for a_i, b_i in zip(S, R)]
        RT = [a_i - b_i for a_i, b_i in zip(T, R)]
        QP = [a_i - b_i for a_i, b_i in zip(P, Q)]
        SP = [a_i - b_i for a_i, b_i in zip(P, S)]

      #  print "printing differences"                                                                                                                                           
      #  print RP,RQ,RS,RT                                                                                                                                                      
        TQ_plane = np.cross(RT, RQ)
      #  print "printing TQ plane"                                                                                                                                              
      #  print RT, RQ                                                                                                                                                           
     #   print TQ_plane                                                                                                                                                         
        QP_plane = np.cross(RQ,RP)
        SP_plane = np.cross(RP,RS)
        ST_plane = np.cross(RS,RT)
        PQS_plane = np.cross(QP,SP)

        equation1 = TQ_plane[0] * - R[0] +  TQ_plane[1] * (- R[1]) +  TQ_plane[2] * (- R[2])
  #      print "printing equation1"                                                                                                                                             
  #      print equation1                                                                                                                                                        
        plane_equation1 =  TQ_plane[0] * dm_xyz[:,0] +  TQ_plane[1] * dm_xyz[:,1] +  TQ_plane[2] * dm_xyz[:,2]
  #      print plane_equation1                                                                                                                                                  

        equation2 = QP_plane[0] * - R[0] +  QP_plane[1] * (- R[1]) +  QP_plane[2] * (- R[2])

        plane_equation2 =  QP_plane[0] * dm_xyz[:,0] +  QP_plane[1] * dm_xyz[:,1] +  QP_plane[2] * dm_xyz[:,2]

        equation3 = SP_plane[0] * - R[0] +  SP_plane[1] * (- R[1]) +  SP_plane[2] * (- R[2])

        plane_equation3 =  SP_plane[0] * dm_xyz[:,0] +  SP_plane[1] * dm_xyz[:,1] +  SP_plane[2] * dm_xyz[:,2]


        equation4 = ST_plane[0] * - R[0] +  ST_plane[1] * (- R[1]) +  ST_plane[2] * (- R[2])

        plane_equation4 =  ST_plane[0] * dm_xyz[:,0] +  ST_plane[1] * dm_xyz[:,1] +  ST_plane[2] * dm_xyz[:,2]



        index2 = (radius < 200) & (plane_equation1 > - equation1) & (plane_equation2 > - equation2) &  (plane_equation3 > - equation3) &  (plane_equation4 > - equation4)

        xyz = xyz[index2]
        
        
                n_bins = 10 #Number of bins you want                                                                                                                                    

        part_per_bin = np.zeros(n_bins)

        n_bins_frac = np.linspace(1.0/n_bins, 1.0, n_bins)
#... Leave this as is. It only sets the percentage of                                                                                                                           
#... particles within each bin given your number of bins                                                                                                                        

#... Binning Function                                                                                                                                                           
        edges = st.mstats.mquantiles(xyz[:,2], n_bins_frac)
#... Where z is the distance from the cone vertex                                                                                                                               
#... of all the particles within the cone                                                                                                                                       

#... For loop that separates all particles in each bin                                                                                                                          

        theta = np.arctan([5.0/208.0])

        print " ######################### printing theta"

      #  print theta                                                                                                                                                            
        dl = []
#xyz[:,2] = 20.8 * xyz[:,0] - 8.0                                                                                                                                               

        for j in range(len(edges)):
            if j == 0:
                mask = (xyz[:,2] > - 8.0) & (xyz[:,2] < edges[j])
            else:
                mask = (xyz[:,2] > edges[j - 1]) & (xyz[:,2] < edges[j])
#...This array now tells you the number of particles in each bin                                                                                                                
            part_per_bin[j] = len(xyz[:,2][mask])
            
            
                   #     print part_per_bin[j]                                                                                                                                              
        #    print edges[j]                                                                                                                                                     
            dl.insert(j,edges[j])
        dl.insert(0,-8.0)
       # print "printing dl "                                                                                                                                                   
       # print dl                                                                                                                                                               
        height1 = [j-i for i, j in zip(dl[:-1], dl[1:])]
       # print "printing height"                                                                                                                                                
       # print height1                                                                                                                                                          
        height = [x * 3.086 * 10**21 for x in height1]
       # print height                                                                                                                                                           
        lengths = [sum(height[:i+1]) for i in range(len(height))]
       # print "printing lengths"                                                                                                                                               
      #  print lengths                                                                                                                                                          
        base = []
        for j in range(len(edges)):
            bases = 2.0 * lengths[j] * tan(theta)
            base.insert(j,bases)
       # print "printing base1"                                                                                                                                                 
       # print base                                                                                                                                                             
       # print "printing base"                                                                                                                                                  



        j_factor = []
        for j in range(len(edges) ):
            if j == 0:
                volume = (1.0/3.0) * height[j] * (base[j] * base[j] )
        #        print "printing triangle volume"                                                                                                                               
         #       print volume[j]                                                                                                                                                
                mass = dm_mass[j] * part_per_bin[j] * (5.609588 * 10**26) * ( 1.998 * 10**30)
          #      print "mass of dm particle"                                                                                                                                    
           #     print dm_mass[j]                                                                                                                                               
                density = mass / volume
                density_squared =  np.power(density, 2)
                j_fact = density_squared[j] * height[j]
            #    print "printing jfactor"                                                                                                                                       
             #   print j_fact                                                                                                                                                   
                j_factor.insert(j,j_fact)
            else:
              #  print "prints and height"                                                                                                                                      
               # print base[j],base[j - 1],height[j]                                                                                                                            
                volume =(1.0/3.0)* height[j]  * (  base[j] * base[j]  +  base[ j - 1] * base[j] + base[j - 1] * base[j - 1] )

                #print "printing volume of chopped pyramind"                                                                                                                    
                #print volume                                                                                                                                                   
                #print "printing dm mass"                                                                                                                                       
                #print dm_mass[j]                                                                                                                                               
                mass = dm_mass[j] * part_per_bin[j] * (5.609588 * 10**26) * ( 1.998 * 10**30)
                density= mass / volume
               # print "printing density"                                                                                                                                       
               # print density                        
               
               
                density_squared =  np.power(density, 2)
                j_fact = density_squared * height[j]
                #print "printing j factor"                                                                                                                                      
                #print j_fact                                                                                                                                                   
                j_factor.insert(j,j_fact)
       # print "printing final j factor"                                                                                                                                        

        #for j in range(len(edges) ):                                                                                                                                           
        #    print "printing j factor"                                                                                                                                          
         #   print j_factor[j]                                                                                                                                                  
       # print "printing the total jfactor"                                                                                                                                     
        print "###################################"
        J_fact = sum(j_factor)
        J.insert(counter, J_fact)
        m = m + 1
        counter = counter + 1
    n = n + 1
print J                                                                                                                                                                        


        
        
 
