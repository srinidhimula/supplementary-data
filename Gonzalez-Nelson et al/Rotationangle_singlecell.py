import linecache
import numpy as np
import re
import os
import glob


"""Function to get the path and filename of .xyz CP2K trajectory file from MD simulations.
This file will be used for calculating the ring rotation angles. Returns path and filename."""
def load_data_fromxyz(base):
    for file in glob.glob(base+ "*.xyz"):
        filename=file
        my_path=os.path.dirname(filename)
        print (filename)
    return filename,my_path


"""start, end are starting and ending time step numbers in the MD simulation."""
start=0
end=101
cell_param=np.array([6.6490,16.382,13.320]) #Single cell lattice parameters of NO2-MIL53 [1] 
num_atoms=84 #Number of atoms in unit cell

"""Function returns the time steps of MD simulation when CP2K trajectory file is given as input file."""
def get_num_steps(filename):
    line_number=2
    num_line=[]
    for i in range(start,end):
        line_number=line_number+(i*(num_atoms+2))
        line=linecache.getline(filename,line_number)
        spl_line=re.split('[,\s]\s*',line)
        num_line.append(spl_line[3])
        line_number=2
    time_steps=np.array(num_line,dtype=int)
    return time_steps

"""Function to read ring1 cartesian coordinates at each timestep and calculates the angle between benzene ring plane 
and reference plane.Takes .xyz trajectory files from CP2K as input. Returns ring rotation angle in degree"""
def ring1_rotation(filename):
    ring1_carbon=np.array([30,72,39,28,43,61]) #ring1 carbon labels 
    ring1coord=[]
    ref_plane=[0,cell_param[2],cell_param[1]] # Normal vector of (011) plane
    rotation_ring1=np.zeros(end-start,dtype=float)
    for j in range(start,end):
        k = j - start
        for i in range(len(ring1_carbon)):
            line_number=ring1_carbon[i]+(j*(num_atoms+2))+2 
            line=linecache.getline(filename,line_number)
            rline=line.split()
            ring1coord.append(rline[1:]) 
        ring1coord=np.array(ring1coord,dtype=float) # Cartesian coordinates of the benzene carbons
        req1_coord=np.array([ring1coord[1,:],ring1coord[0,:],ring1coord[1,:]])
        rotation_ring1[k]=angle_calculate(ring1coord,ref_plane,req1_coord)
        ring1coord=[]
    np.savetxt('ring1.dat',np.column_stack([time_steps,rotation_ring1]),fmt='%d  %1.5f')
    return rotation_ring1

"""Function to read ring2 cartesian coordinates, correct for periodicity (i.e., obtain a complete ring) at each 
timestep and calculates the angle between benzene ring plane and reference plane. Takes .xyz trajectory files from CP2K as input.
Returns ring rotation angle in degree"""
def ring2_rotation(filename):
    ring2_carbon=np.array([29,71,40,27,44,62]) #ring2 carbon labels
    ring2coord=[]
    rotation_ring2=np.zeros(end-start,dtype=float)
    ref_plane=[0,cell_param[2],-cell_param[1]] #Normal vector of (0-11) plane
    for j in range(start,end):
        k = j - start
        for i in range(len(ring2_carbon)):
            line_number=ring2_carbon[i]+(j*(num_atoms+2))+2
            line=linecache.getline(filename,line_number)
            rline=line.split()
            ring2coord.append(rline[1:]) 
        ring2coord=np.array(ring2coord,dtype=float) # Cartesian coordinates of the benzene carbons
        ring2coord[3:6,1]=ring2coord[3:6,1]+cell_param[1] # Unit cell parameters are added to complete the ring to account for periodic boundary conditions
        req2_coord=np.array([ring2coord[0,:],ring2coord[1,:],ring2coord[4,:]])
        rotation_ring2[k]=angle_calculate(ring2coord,ref_plane,req2_coord)
        ring2coord=[]
    np.savetxt('ring2.dat',np.column_stack([time_steps,rotation_ring2]),fmt='%d  %1.5f')
    return rotation_ring2


"""Function to read ring3 cartesian coordinates, correct for periodicity (i.e., obtain a complete ring) at each 
timestep and calculates the angle between benzene ring plane and reference plane. Takes .xyz trajectory files from CP2K as input.
Returns ring rotation angle in degree"""
def ring3_rotation(filename):
    ring3_carbon=np.array([32,70,37,26,41,63]) # ring3 carbon labels
    ring3coord=[]
    rotation_ring3=np.zeros(end-start,dtype=float)
    ref_plane=[0,cell_param[2],-cell_param[1]] # Normal vector of (0-11) plane
    for j in range(start,end):
        k = j - start
        for i in range(len(ring3_carbon)):
            line_number=ring3_carbon[i]+(j*(num_atoms+2))+2
            line=linecache.getline(filename,line_number)
            rline=line.split()
            ring3coord.append(rline[1:]) 
        ring3coord=np.array(ring3coord,dtype=float) # Cartesian coordinates of the benzene carbons
        ring3coord[2,0]=ring3coord[2,0]+cell_param[0] #Unit cell parameters are added to complete the ring to account for periodic boundary conditions
        ring3coord[3,0]=ring3coord[3,0]+cell_param[0]
        ring3coord[3:6,2]=ring3coord[3:6,2]+cell_param[2]
        req3_coord=np.array([ring3coord[0,:],ring3coord[1,:],ring3coord[4,:]])
        rotation_ring3[k]=angle_calculate(ring3coord,ref_plane,req3_coord)
        ring3coord=[]
    np.savetxt('ring3.dat',np.column_stack([time_steps,rotation_ring3]),fmt='%d  %1.5f')
    return rotation_ring3

"""Function to read ring4 cartesian coordinates, correct for periodicity (i.e., obtain a complete ring) at each 
timestep and calculates the angle between benzene ring plane and reference plane. Takes .xyz trajectory files from CP2K as input.
Returns ring rotation angle in degree"""
def ring4_rotation(filename):
    ring4_carbon=np.array([31,69,38,25,42,64]) #ring4 carbon labels
    ring4coord=[]
    ref_plane=[0,cell_param[2],cell_param[1]] #Normal vector of (011) plane
    rotation_ring4=np.zeros(end-start,dtype=float)
    for j in range(start,end):
        k = j - start
        for i in range(len(ring4_carbon)):
            line_number=ring4_carbon[i]+(j*(num_atoms+2))+2
            line=linecache.getline(filename,line_number)
            rline=line.split()
            ring4coord.append(rline[1:]) 
        ring4coord=np.array(ring4coord,dtype=float) # Cartesian coordinates of the benzene carbons
        ring4coord[2:4,0]=ring4coord[2:4,0]-cell_param[0] #Unit cell parameters are added to complete the ring to account for periodic boundary conditions
        ring4coord[3:6,1]=ring4coord[3:6,1]+cell_param[1]
        ring4coord[3:6,2]=ring4coord[3:6,2]-cell_param[2]
        req4_coord=np.array([ring4coord[1,:],ring4coord[0,:],ring4coord[1,:]])
        rotation_ring4[k]=angle_calculate(ring4coord,ref_plane,req4_coord)
        ring4coord=[]
    np.savetxt('ring4.dat',np.column_stack([time_steps,rotation_ring4]),fmt='%d  %1.5f')
    return rotation_ring4

"""Function to calculate the rotation angle of benzene rings as per notation in Figure2 of the article.
Angle between normal of benzene ring plane and normal of (011)/(0-11) plane is calculated."""
def angle_calculate(ring_coord,ref_plane,req_coord):
    average=np.average(ring_coord,axis=0)
    ring_coord=ring_coord-average
    req_coord=req_coord-average
    ring_coord=ring_coord.T
    svd=np.linalg.svd(ring_coord) 
    svd [1] [2] # The normal vector of the best-fitting plane is the left singular vector corresponding to the least singular value.
    normal = np.transpose(svd [0]) [2] # normal vector of the plane of benzene carbons
    ring_coord=ring_coord.T
    updir=np.cross(req_coord[0,:],req_coord[1,:])
    test_direction=np.dot(updir,normal) # Consistent definition of ring normal direction
    dir_correction=test_direction/np.absolute(test_direction) 
    normal_corrected=normal*dir_correction # This normal is  used to calculate the rotation angle
    ref_plane = ref_plane / np.linalg.norm(ref_plane)
    dot=np.dot(normal_corrected,ref_plane) # Dot product gives the angle between the reference planes and ring plane
    
    rotcross=np.cross(normal_corrected,ref_plane)
    pos_neg=np.sign(np.dot(rotcross,req_coord[2,:])) # Sign for the angle based on the notation in Figure 2
    angle=pos_neg*np.arccos(np.clip(dot,-1,1))
    angle_degree=np.degrees(angle)
    return angle_degree
    
#[1]Biswas,S., Ahnfeldt, T. & Stock, N. Inorg. Chem. 50, 9518–9526 (2011).


"""for example, if .xyz file is in input_demo directory, following can be used to read the files and calculate the
rotation angles. In this case start, end variables have to be changed to 0 and 101 since .xyz files has only 100 MD timesteps.
.dat files with time step and rotation angle are saved."""
filename,my_path=load_data_fromxyz('input_demo/') 
time_steps= get_num_steps(filename)
rotation_ring1=ring1_rotation(filename)
rotation_ring2=ring2_rotation(filename)
rotation_ring3=ring3_rotation(filename)
rotation_ring4=ring4_rotation(filename)
