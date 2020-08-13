import numpy as np
import math
import matplotlib.pyplot as plt
import os.path as path

def Check_Number_breaks(line):
    b=np.zeros(4)
    c=0
    a=0
    for i in line:
        #print(i)
        if i==' ':
            b[c]=int(a)
            c+=1
        a+=1
    b[3]=len(line)
    #print(b)
    return(b.astype(int))

def Frame_Array(frame):
    b=np.array([2+27*(frame),27*(frame+1)])
    return(b.astype(int))

def Frame_Iterate(myarray,frame):

    data_lines=np.zeros(2)
    splitted=np.empty((25,3))
    
    data_lines=Frame_Array(frame)
    data_lines=data_lines.astype(int)    
    splitted= np.split(
        myarray,[data_lines[0],data_lines[1]],axis=0)
    return(splitted[1])

def Open_File(filename):
    f = open(filename,"r")  #open the xyz file
    contents=f.readlines()
    myarray = np.array(contents)
    return(contents)

def Data_Array_xyz(myarray,current_frame):

    data=Frame_Iterate(myarray,current_frame)
    #break up the array into just the data which in this case is lines 3-203 
    length=len(data)

    a=0   #initialize everything
    c=0  
    x_coordinate = np.array([])
    y_coordinate = np.array([])
    z_coordinate = np.array([])
    b=np.zeros(4)

        
    while a<length:  #iterate over the length of the saved data and save each column as float numbers
        for line in data:
            index=Check_Number_breaks(line)
            x_coordinate = np.append(x_coordinate,float(data[a][index[0]+1:index[1]]))
            y_coordinate = np.append(y_coordinate,float(data[a][index[1]+1:index[2]]))
            #print(data[a][index[2]+1:index[3]])
            z_coordinate = np.append(z_coordinate,float(data[a][index[2]+1:index[3]]))
            a+=1
      #transform numbers into integers
    data=np.vstack((x_coordinate.T,y_coordinate.T,z_coordinate.T))
    return(data.T)

def XYZ_To_NParray(filename):
    myarray=Open_File(filename)
    Frame_Total=int(len(myarray)/27)
    data=np.zeros((Frame_Total,25,3))
    current_frame=0
    while current_frame<Frame_Total:
        data[current_frame]=Data_Array_xyz(myarray,current_frame)
        current_frame+=1
    return(data)