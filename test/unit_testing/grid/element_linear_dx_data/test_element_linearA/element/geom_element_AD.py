from __future__ import division

import torch
import torch.autograd as autograd
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
import sys
import os
import time


#
# TORCH INSTALLATION: refer to https://pytorch.org/get-started/locally/
#



def update_progress(job_title, progress):
    length = 20 # modify this to change the length
    block = int(round(length*progress))
    msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
    if progress >= 1: msg += " DONE\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()

def cls():
    os.system('cls' if os.name=='nt' else 'clear')


cls()


################################################################################################################

# Initialize torch tensor for coordiantes 
coords_data  = [[0.0,0.0,0.0],
                [5.0,0.0,0.0],
                [0.0,1.0,0.0],
                [5.0,1.0,0.0],
                [0.0,0.0,1.0],
                [5.0,0.0,1.0],
                [0.0,1.0,1.0],
                [5.0,1.0,1.0],
                ] 
coords = torch.tensor(coords_data,requires_grad=True,dtype=torch.float64)


nnodes_r   = coords.size(0)
nnodes_ie  = 8 
nnodes_if  = 4 
nterms_s   = 8
ndirs      = 3
coord_sys  = 'CARTESIAN'



# Define matrix of polynomial basis terms at support nodes
val_r_data = [[ 1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0,-1.0],
              [ 1.0,-1.0,-1.0, 1.0,-1.0,-1.0, 1.0, 1.0],
              [ 1.0, 1.0,-1.0,-1.0,-1.0, 1.0,-1.0, 1.0],
              [ 1.0, 1.0,-1.0, 1.0, 1.0,-1.0,-1.0,-1.0],
              [ 1.0,-1.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0],
              [ 1.0,-1.0, 1.0, 1.0,-1.0, 1.0,-1.0,-1.0],
              [ 1.0, 1.0, 1.0,-1.0,-1.0,-1.0, 1.0,-1.0],
              [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
              ] 
val_r = torch.tensor(val_r_data,requires_grad=False,dtype=torch.float64)


# Define matrices at interpolation nodes (quadrature, level = 1)
val_i_data = [[ 1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0, 1.0/3.0,-1.0/3.0*np.sqrt(1.0/3.0)],
              [ 1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0, 1.0/3.0, 1.0/3.0*np.sqrt(1.0/3.0)],
              [ 1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0,-1.0/3.0, 1.0/3.0*np.sqrt(1.0/3.0)],
              [ 1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0,-1.0/3.0,-1.0/3.0*np.sqrt(1.0/3.0)],
              [ 1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0,-1.0/3.0, 1.0/3.0*np.sqrt(1.0/3.0)],
              [ 1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0,-1.0/3.0,-1.0/3.0*np.sqrt(1.0/3.0)],
              [ 1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0, 1.0/3.0,-1.0/3.0*np.sqrt(1.0/3.0)],
              [ 1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0*np.sqrt(1.0/3.0)],
              ]
val_i = torch.tensor(val_i_data,requires_grad=False,dtype=torch.float64)


ddxi_i_data  = [[ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                [ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                [ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                [ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                ]
ddxi_i = torch.tensor(ddxi_i_data,requires_grad=False,dtype=torch.float64)

ddeta_i_data = [[ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0), 1.0/3.0],
                [ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0), 1.0/3.0],
                ]
ddeta_i = torch.tensor(ddeta_i_data,requires_grad=False,dtype=torch.float64)

ddzeta_i_data= [[ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                [ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                ]
ddzeta_i = torch.tensor(ddzeta_i_data,requires_grad=False,dtype=torch.float64)

# Define element interpolation nodes weights for linear element
weights_e_data = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0] 
weights_e = torch.tensor(weights_e_data,requires_grad=False,dtype=torch.float64)


# Define val_f for each face
# Face 1, XI_MIN
val_1_data   = [[ 1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
               [ 1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
               [ 1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
               [ 1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
               ]
val_1 = torch.tensor(val_1_data,requires_grad=False,dtype=torch.float64)

# Face 2, XI_MAX
val_2_data   = [[ 1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                [ 1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                [ 1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                [ 1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                ] 
val_2 = torch.tensor(val_2_data,requires_grad=False,dtype=torch.float64)

# Face 3, ETA_MIN
val_3_data   = [[ 1.0,-1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0),-1.0/3.0],
                [ 1.0,-1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0), 1.0/3.0],
                [ 1.0,-1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 1.0,-1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                ] 
val_3 = torch.tensor(val_3_data,requires_grad=False,dtype=torch.float64)

# Face 4, ETA_MAX
val_4_data   = [[ 1.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 1.0,1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                [ 1.0,1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0),-1.0/3.0],
                [ 1.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0), 1.0/3.0],
                ]
val_4 = torch.tensor(val_4_data,requires_grad=False,dtype=torch.float64)

# Face 5, ZETA_MIN
val_5_data   = [[ 1.0,-np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                [ 1.0,-np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                [ 1.0, np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 1.0, np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                ]
val_5 = torch.tensor(val_5_data,requires_grad=False,dtype=torch.float64)

# Face 6, ZETA_MAX
val_6_data   = [[ 1.0,-np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                [ 1.0,-np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                [ 1.0, np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                [ 1.0, np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                ]
val_6 = torch.tensor(val_6_data,requires_grad=False,dtype=torch.float64)





#--------------------------------------------------------------------

# Matrix modes_to_nodes
val_r_inv = torch.inverse(val_r)


# Computes coordiantes modes
coords_modes = torch.mm(val_r_inv,coords)


# Initialized coordiantes
interp_coords = torch.mm(val_i,coords_modes)


# Initialized jacobian
jacobian = torch.empty(3,3,nnodes_ie, dtype=torch.float64) 
for inode in range(0,nnodes_ie):
    jacobian[0,0,inode] = torch.dot(ddxi_i[inode,:]   , coords_modes[:,0])
    jacobian[0,1,inode] = torch.dot(ddeta_i[inode,:]  , coords_modes[:,0])
    jacobian[0,2,inode] = torch.dot(ddzeta_i[inode,:] , coords_modes[:,0])
    jacobian[1,0,inode] = torch.dot(ddxi_i[inode,:]   , coords_modes[:,1])
    jacobian[1,1,inode] = torch.dot(ddeta_i[inode,:]  , coords_modes[:,1])
    jacobian[1,2,inode] = torch.dot(ddzeta_i[inode,:] , coords_modes[:,1])
    jacobian[2,0,inode] = torch.dot(ddxi_i[inode,:]   , coords_modes[:,2])
    jacobian[2,1,inode] = torch.dot(ddeta_i[inode,:]  , coords_modes[:,2])
    jacobian[2,2,inode] = torch.dot(ddzeta_i[inode,:] , coords_modes[:,2])
    update_progress("Computing Jacobian                 ", inode/(nnodes_ie-1))
if coord_sys == 'CYLINDRICAL':
    scaling_factor = torch.mm(val_i,coords_modes[:,0])
    for inode in range(0,nnodes_ie):
        jacobian[1,0,inode] = jacobian[1,0,inode] * scaling_factor[inode]
        jacobian[1,1,inode] = jacobian[1,1,inode] * scaling_factor[inode]
        jacobian[1,2,inode] = jacobian[1,2,inode] * scaling_factor[inode]



# Matrics and Determinant
metrics = torch.empty(3,3,nnodes_ie, dtype=torch.float64) 
jinv    = torch.empty(nnodes_ie, dtype=torch.float64) 
for inode in range(0,nnodes_ie):
    ijacobian = torch.empty(3,3, dtype=torch.float64)
    imetric   = torch.empty(3,3, dtype=torch.float64)
    for irow in range(0,3):
        for icol in range(0,3):
            ijacobian[irow,icol] = jacobian[irow,icol,inode]
    # Compute jacobian for the ith node
    update_progress("Computing Jinv and Metric          ", inode/(nnodes_ie-1))
    jinv[inode] = torch.det(ijacobian) 
    imetric     = torch.inverse(ijacobian)
    for irow in range(0,3):
        for icol in range(0,3):
            metrics[irow,icol,inode] = imetric[irow,icol]

 
# Compute inverse Mass matrix
invmass = torch.empty(nterms_s,nterms_s,nnodes_ie, dtype=torch.float64)
mass    = torch.empty(nterms_s,nterms_s,nnodes_ie, dtype=torch.float64)
val_tmp = torch.empty(nterms_s,nnodes_ie, dtype=torch.float64)
i = 1
for iterm in range(0,nterms_s):
    for inode in range(0,nnodes_ie):
        val_tmp[inode,iterm] = val_i[inode,iterm] * weights_e[inode] * jinv[inode]
        update_progress("Computing invmass                  ", i/(nterms_s*nnodes_ie))
        i += 1
mass    = torch.mm(torch.t(val_tmp),val_i)
invmass = torch.inverse(mass)


# Compute BR2_VOL for each face
br2_vol_face1 = torch.mm(val_i,torch.mm(invmass,torch.t(val_1)))
br2_vol_face2 = torch.mm(val_i,torch.mm(invmass,torch.t(val_2)))
br2_vol_face3 = torch.mm(val_i,torch.mm(invmass,torch.t(val_3)))
br2_vol_face4 = torch.mm(val_i,torch.mm(invmass,torch.t(val_4)))
br2_vol_face5 = torch.mm(val_i,torch.mm(invmass,torch.t(val_5)))
br2_vol_face6 = torch.mm(val_i,torch.mm(invmass,torch.t(val_6)))
update_progress("Computing br2_vol                  ", 1)


# Compute BR2_FACE for each face
br2_face_face1 = torch.mm(val_1,torch.mm(invmass,torch.t(val_1)))
br2_face_face2 = torch.mm(val_2,torch.mm(invmass,torch.t(val_2)))
br2_face_face3 = torch.mm(val_3,torch.mm(invmass,torch.t(val_3)))
br2_face_face4 = torch.mm(val_4,torch.mm(invmass,torch.t(val_4)))
br2_face_face5 = torch.mm(val_5,torch.mm(invmass,torch.t(val_5)))
br2_face_face6 = torch.mm(val_6,torch.mm(invmass,torch.t(val_6)))
update_progress("Computing br2_face                 ", 1)


# Grad1, Grad2, and Grad3
grad1 = torch.empty(nnodes_ie,nterms_s, dtype=torch.float64)
grad2 = torch.empty(nnodes_ie,nterms_s, dtype=torch.float64)
grad3 = torch.empty(nnodes_ie,nterms_s, dtype=torch.float64)
i = 1
for iterm in range(0,nterms_s):
    for inode in range(0,nnodes_ie):
        grad1[inode,iterm] = metrics[0,0,inode] * ddxi_i[inode,iterm] + metrics[1,0,inode] * ddeta_i[inode,iterm] + metrics[2,0,inode] * ddzeta_i[inode,iterm] 
        grad2[inode,iterm] = metrics[0,1,inode] * ddxi_i[inode,iterm] + metrics[1,1,inode] * ddeta_i[inode,iterm] + metrics[2,1,inode] * ddzeta_i[inode,iterm] 
        grad3[inode,iterm] = metrics[0,2,inode] * ddxi_i[inode,iterm] + metrics[1,2,inode] * ddeta_i[inode,iterm] + metrics[2,2,inode] * ddzeta_i[inode,iterm] 
        update_progress("Computing grad1, grad2, grad3      ", i/(nnodes_ie*nterms_s))
        i += 1







#WRITE_____________________

#
# Metrics
#
f = open("metrics.txt","w")
i = 1
for inode in range (0,nnodes_ie):
    f.write("Metric interpolation node %d \n" % (inode+1))
    array = np.zeros([3, 3])
    for irow in range(0,3):
        for icol in range(0,3):
            array[irow,icol] = metrics[irow,icol,inode].item() 
            update_progress("Writing metrics to file            ", i/(nnodes_ie*9))
            i += 1
    np.savetxt(f,array)
f.close()

#
# jinv
#
f = open("jinv.txt","w")
array = np.zeros([1])
i = 1
for inode in range (0,nnodes_ie):
    f.write("Jinv interpolation node %d \n" % (inode+1))
    array[0] = jinv[inode].item()
    np.savetxt(f,array)
    update_progress("Writing jinv to file               ", i/(nnodes_ie))
    i += 1
f.close()

#
# Grad1
#
f = open("grad1.txt","w")
f.write("Grad1 \n")
array = np.zeros([nnodes_ie,nterms_s])
i = 1
for inode in range (0,nnodes_ie):
    for iterm in range(0,nterms_s):
        array[inode,iterm] = grad1[inode,iterm].item()
        update_progress("Writing grad1 to file              ", i/(nnodes_ie*nterms_s))
        i += 1
np.savetxt(f,array)
f.close()

#
# Grad2
#
f = open("grad2.txt","w")
f.write("Grad2 \n")
array = np.zeros([nnodes_ie,nterms_s])
i = 1
for inode in range (0,nnodes_ie):
    for iterm in range(0,nterms_s):
        array[inode,iterm] = grad2[inode,iterm].item()
        update_progress("Writing grad2 to file              ", i/(nnodes_ie*nterms_s))
        i += 1
np.savetxt(f,array)
f.close()

#
# Grad3
#
f = open("grad3.txt","w")
f.write("Grad3 \n")
array = np.zeros([nnodes_ie,nterms_s])
i = 1
for inode in range (0,nnodes_ie):
    for iterm in range(0,nterms_s):
        array[inode,iterm] = grad3[inode,iterm].item()
        update_progress("Writing grad3 to file              ", i/(nnodes_ie*nterms_s))
        i += 1
np.savetxt(f,array)
f.close()

#
# dmetric_dx
#
f = open("dmetric_dx.txt","w")
i = 1
for inode in range (0,nnodes_ie):
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            array = np.zeros([3,3])
            f.write("dmetric_dx interpolation node %s, diff_node %s, diff_dir %s \n" % (inode+1,inode_diff+1,idir+1))
            for irow in range(0,3):
                for icol in range(0,3):
                    data = metrics[irow,icol,inode]
                    data.backward(retain_graph=True)
                    ddata = coords.grad
                    ddata_np = ddata.numpy()
                    array[irow,icol] = ddata_np[inode_diff,idir]
                    update_progress("Writing dmetric_dx to file         ", i/(nnodes_ie*nnodes_r*ndirs*3*3))
                    # This avoid to accumulate derivatives
                    dummy = coords.grad.data.zero_()
                    i += 1
            np.savetxt(f,array)
f.close()

#
# interp_coords_dx
#
f = open("dinterp_xcoords_dx.txt","w")
i = 1
f.write("xcoord interpolation, coord 1,  row=node, col=nnodes_r*dir \n")
array = np.zeros([nnodes_ie,nnodes_r*ndirs])
for inode in range (0,nnodes_ie):
    data = interp_coords[inode,0]
    data.backward(retain_graph=True)
    ddata = coords.grad
    ddata_np = ddata.numpy()
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            if idir == 0:
                index = inode_diff
            elif idir == 1:
                index = nnodes_r + inode_diff
            elif idir == 2:
                index = 2*nnodes_r + inode_diff
            array[inode,index] = ddata_np[inode_diff,idir]
            update_progress("Writing interp_xcoords_dx to file  ", i/(nnodes_ie*nnodes_r*3))
            i += 1
    # This avoid to accumulate derivatives
    dummy = coords.grad.data.zero_()
np.savetxt(f,array)
f.close()
f = open("dinterp_ycoords_dx.txt","w")
i = 1
f.write("ycoord interpolation, coord 2,  row=node, col=nnodes_r*dir \n")
array = np.zeros([nnodes_ie,nnodes_r*ndirs])
for inode in range (0,nnodes_ie):
    data = interp_coords[inode,1]
    data.backward(retain_graph=True)
    ddata = coords.grad
    ddata_np = ddata.numpy()
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            if idir == 0:
                index = inode_diff
            elif idir == 1:
                index = nnodes_r + inode_diff
            elif idir == 2:
                index = 2*nnodes_r + inode_diff
            array[inode,index] = ddata_np[inode_diff,idir]
            update_progress("Writing interp_ycoords_dx to file  ", i/(nnodes_ie*nnodes_r*3))
            i += 1
    # This avoid to accumulate derivatives
    dummy = coords.grad.data.zero_()
np.savetxt(f,array)
f.close()
f = open("dinterp_zcoords_dx.txt","w")
i = 1
f.write("zcoord interpolation, coord 3,  row=node, col=nnodes_r*dir \n")
array = np.zeros([nnodes_ie,nnodes_r*ndirs])
for inode in range (0,nnodes_ie):
    data = interp_coords[inode,2]
    data.backward(retain_graph=True)
    ddata = coords.grad
    ddata_np = ddata.numpy()
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            if idir == 0:
                index = inode_diff
            elif idir == 1:
                index = nnodes_r + inode_diff
            elif idir == 2:
                index = 2*nnodes_r + inode_diff
            array[inode,index] = ddata_np[inode_diff,idir]
            update_progress("Writing interp_zcoords_dx to file  ", i/(nnodes_ie*nnodes_r*3))
            i += 1
    # This avoid to accumulate derivatives
    dummy = coords.grad.data.zero_()
np.savetxt(f,array)
f.close()

#
# djinv_dx
#
f = open("djinv_dx.txt","w")
i = 1
for inode in range (0,nnodes_ie):
    array = np.zeros([nnodes_r,ndirs])
    f.write("djinv_dx interpolation node %s, row=inode_diff, col=dir \n" % (inode+1))
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
                data = jinv[inode]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[inode_diff,idir] = ddata_np[inode_diff,idir]
                update_progress("Writing djinv_dx to file           ", i/(nnodes_ie*nnodes_r*ndirs))
                dummy = coords.grad.data.zero_()
                i += 1
    np.savetxt(f,array)
f.close()

#
# dmass_dx
# 
f = open("dmass_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dmass_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nterms_s,nterms_s])
        for irow in range(0,nterms_s):
            for icol in range(0,nterms_s):
                data = mass[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dmass_dx to file           ", i/(nterms_s*nnodes_r*ndirs*nterms_s))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()

#
# dinvmass_dx
# 
f = open("dinvmass_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dinvmass_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nterms_s,nterms_s])
        for irow in range(0,nterms_s):
            for icol in range(0,nterms_s):
                data = invmass[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dinvmass_dx to file        ", i/(nterms_s*nnodes_r*ndirs*nterms_s))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()

#
# dbr2_vol_dx
#
# 
f = open("dbr2_vol_face1_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face1_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data = br2_vol_face1[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_vol_face1_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_vol_face2_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face2_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data = br2_vol_face2[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_vol_face2_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_vol_face3_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face3_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data = br2_vol_face3[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_vol_face3_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_vol_face4_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face4_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data = br2_vol_face4[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_vol_face4_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_vol_face5_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face5_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data = br2_vol_face5[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_vol_face5_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_vol_face6_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face6_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data = br2_vol_face6[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_vol_face6_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()

#
# dbr2_face_dx
#
# 
f = open("dbr2_face_face1_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face1_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data = br2_face_face1[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_face_face1_dx to file ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_face_face2_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face2_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data = br2_face_face2[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_face_face2_dx to file ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_face_face3_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face3_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data = br2_face_face3[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_face_face3_dx to file ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_face_face4_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face4_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data = br2_face_face4[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_face_face4_dx to file ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_face_face5_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face5_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data = br2_face_face5[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_face_face5_dx to file ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()
f = open("dbr2_face_face6_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face6_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data = br2_face_face6[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dbr2_face_face6_dx to file ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()

#
# dgrad1_dx
#
f = open("dgrad1_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dgrad1_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nterms_s])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nterms_s):
                data = grad1[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dgrad1_dx to file          ", i/(nnodes_ie*nnodes_r*ndirs*nterms_s))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()

#
# dgrad2_dx
#
f = open("dgrad2_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dgrad2_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nterms_s])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nterms_s):
                data = grad2[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dgrad2_dx to file          ", i/(nnodes_ie*nnodes_r*ndirs*nterms_s))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()

#
# dgrad3_dx
#
f = open("dgrad3_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dgrad3_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_ie,nterms_s])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nterms_s):
                data = grad3[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dgrad3_dx to file          ", i/(nnodes_ie*nnodes_r*ndirs*nterms_s))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()








