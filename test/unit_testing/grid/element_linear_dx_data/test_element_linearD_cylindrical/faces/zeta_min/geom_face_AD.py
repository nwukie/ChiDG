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
coords_data  = [[ 1.0,   0.0,       0.0],
                [ 2.0,   0.0,       0.0],
                [ 2.0,   np.pi/2.0, 0.0],
                [ 1.0,   np.pi/2.0, 0.0],
                [ 1.0,   0.0,       1.0],
                [ 2.0,   0.0,       1.0],
                [ 2.0,   np.pi/2.0, 1.0],
                [ 1.0,   np.pi/2.0, 1.0],
                ]
coords = torch.tensor(coords_data,requires_grad=True,dtype=torch.float64)


nnodes_r   = coords.size(0)
nnodes_ie  = 8 
nnodes_if  = 4 
nterms_s   = 8
ndirs      = 3
coord_sys  = 'CYLINDRICAL'
iface      = 5



# Define matrix of polynomial basis terms at support nodes
val_r_data = [[ 1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0,-1.0],
              [ 1.0,-1.0,-1.0, 1.0,-1.0,-1.0, 1.0, 1.0],
              [ 1.0, 1.0,-1.0, 1.0, 1.0,-1.0,-1.0,-1.0],
              [ 1.0, 1.0,-1.0,-1.0,-1.0, 1.0,-1.0, 1.0],
              [ 1.0,-1.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0],
              [ 1.0,-1.0, 1.0, 1.0,-1.0, 1.0,-1.0,-1.0],
              [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
              [ 1.0, 1.0, 1.0,-1.0,-1.0,-1.0, 1.0,-1.0],
              ] 
val_r = torch.tensor(val_r_data,requires_grad=False,dtype=torch.float64)



# Define matrices at interpolation nodes (quadrature, level = 1)
if iface == 1:
    val_i_data    = [[ 1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
                     [ 1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
                     [ 1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
                     [ 1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
                     ] 
    val_i = torch.tensor(val_i_data,requires_grad=False,dtype=torch.float64)

    ddxi_i_data   = [[ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                     [ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                     ]
    ddxi_i = torch.tensor(ddxi_i_data,requires_grad=False,dtype=torch.float64)

    ddeta_i_data  = [[ 0.0,1.0,0.0,0.0,-1.0,0.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,-1.0,0.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,-1.0,0.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,-1.0,0.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     ] 
    ddeta_i = torch.tensor(ddeta_i_data,requires_grad=False,dtype=torch.float64)

    ddzeta_i_data = [[ 0.0,0.0,1.0,0.0,0.0,-1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,-1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,-1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,-1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     ] 
    ddzeta_i = torch.tensor(ddzeta_i_data,requires_grad=False,dtype=torch.float64)


if iface == 2:
    val_i_data    = [[ 1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                     [ 1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                     [ 1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                     [ 1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                     ] 
    val_i = torch.tensor(val_i_data,requires_grad=False,dtype=torch.float64)

    ddxi_i_data   = [[ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                     [ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0,-1.0/3.0],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),0.0, 1.0/3.0],
                     ] 
    ddxi_i = torch.tensor(ddxi_i_data,requires_grad=False,dtype=torch.float64)

    ddeta_i_data  = [[ 0.0,1.0,0.0,0.0,1.0,0.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,1.0,0.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,1.0,0.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,1.0,0.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     ] 
    ddeta_i = torch.tensor(ddeta_i_data,requires_grad=False,dtype=torch.float64)


    ddzeta_i_data = [[ 0.0,0.0,1.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0)],
                     ] 
    ddzeta_i = torch.tensor(ddzeta_i_data,requires_grad=False,dtype=torch.float64)



if iface == 3:
    val_i_data    = [[ 1.0,-1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 1.0,-1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0), 1.0/3.0],
                     [ 1.0,-1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                     [ 1.0,-1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                     ]
    val_i = torch.tensor(val_i_data,requires_grad=False,dtype=torch.float64)

    ddxi_i_data   = [[ 0.0,0.0,0.0,1.0,-1.0,-np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,-1.0,-np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,-1.0, np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,-1.0, np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0)],
                     ] 
    ddxi_i = torch.tensor(ddxi_i_data,requires_grad=False,dtype=torch.float64)

    ddeta_i_data  = [[ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0), 1.0/3.0],
                     ]
    ddeta_i = torch.tensor(ddeta_i_data,requires_grad=False,dtype=torch.float64)


    ddzeta_i_data = [[ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0)],
                     ]
    ddzeta_i = torch.tensor(ddzeta_i_data,requires_grad=False,dtype=torch.float64)



if iface == 4:
    val_i_data    = [[ 1.0,1.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                     [ 1.0,1.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 1.0,1.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 1.0,1.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0), 1.0/3.0],
                     ] 
    val_i = torch.tensor(val_i_data,requires_grad=False,dtype=torch.float64)

    ddxi_i_data   = [[ 0.0,0.0,0.0,1.0,1.0,-np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,1.0,-np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,1.0, np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,1.0, np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0)],
                     ]
    ddxi_i = torch.tensor(ddxi_i_data,requires_grad=False,dtype=torch.float64)

    ddeta_i_data  = [[ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0), 1.0/3.0],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,-np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0, np.sqrt(1.0/3.0), 1.0/3.0],
                     ]
    ddeta_i = torch.tensor(ddeta_i_data,requires_grad=False,dtype=torch.float64)


    ddzeta_i_data = [[ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0)],
                     ]
    ddzeta_i = torch.tensor(ddzeta_i_data,requires_grad=False,dtype=torch.float64)



if iface == 5:
    val_i_data   = [[ 1.0,-np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                    [ 1.0,-np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                    [ 1.0, np.sqrt(1.0/3.0),-1.0,-np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                    [ 1.0, np.sqrt(1.0/3.0),-1.0, np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                    ]
    val_i = torch.tensor(val_i_data,requires_grad=False,dtype=torch.float64)

    ddxi_i_data   = [[ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-1.0,0.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),-1.0,0.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),-1.0,0.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),-1.0,0.0,-np.sqrt(1.0/3.0)],
                     ]
    ddxi_i = torch.tensor(ddxi_i_data,requires_grad=False,dtype=torch.float64)

    ddeta_i_data  = [[ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,-1.0, np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,-1.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,-1.0, np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,-1.0,-np.sqrt(1.0/3.0)],
                     ]
    ddeta_i = torch.tensor(ddeta_i_data,requires_grad=False,dtype=torch.float64)


    ddzeta_i_data = [[ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                     ]
    ddzeta_i = torch.tensor(ddzeta_i_data,requires_grad=False,dtype=torch.float64)



if iface == 6:
    val_i_data   = [[ 1.0,-np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0), 1.0/3.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                    [ 1.0,-np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0),-1.0/3.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                    [ 1.0, np.sqrt(1.0/3.0),1.0,-np.sqrt(1.0/3.0),-1.0/3.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                    [ 1.0, np.sqrt(1.0/3.0),1.0, np.sqrt(1.0/3.0), 1.0/3.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                    ]
    val_i = torch.tensor(val_i_data,requires_grad=False,dtype=torch.float64)

    ddxi_i_data   = [[ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),1.0,0.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0,-np.sqrt(1.0/3.0),1.0,0.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),1.0,0.0, np.sqrt(1.0/3.0)],
                     [ 0.0,0.0,0.0,1.0, np.sqrt(1.0/3.0),1.0,0.0, np.sqrt(1.0/3.0)],
                     ]
    ddxi_i = torch.tensor(ddxi_i_data,requires_grad=False,dtype=torch.float64)

    ddeta_i_data  = [[ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,1.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,1.0, np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),0.0,1.0,-np.sqrt(1.0/3.0)],
                     [ 0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),0.0,1.0, np.sqrt(1.0/3.0)],
                     ]
    ddeta_i = torch.tensor(ddeta_i_data,requires_grad=False,dtype=torch.float64)


    ddzeta_i_data = [[ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0), 1.0/3.0],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0),-np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,0.0,1.0,0.0,0.0,-np.sqrt(1.0/3.0), np.sqrt(1.0/3.0),-1.0/3.0],
                     [ 0.0,0.0,1.0,0.0,0.0, np.sqrt(1.0/3.0), np.sqrt(1.0/3.0), 1.0/3.0],
                     ]
    ddzeta_i = torch.tensor(ddzeta_i_data,requires_grad=False,dtype=torch.float64)






#--------------------------------------------------------------------

# Matrix modes_to_nodes
val_r_inv = torch.inverse(val_r)


# Computes coordiantes modes
coords_modes = torch.mm(val_r_inv,coords)


# Initialized coordiantes
interp_coords = torch.mm(val_i,coords_modes)


# Initialized jacobian
jacobian     = torch.empty(3,3,nnodes_if, dtype=torch.float64) 
jacobian_tmp = torch.empty(3,3,nnodes_if, dtype=torch.float64) 
for inode in range(0,nnodes_if):
    jacobian_tmp[0,0,inode] = torch.dot(ddxi_i[inode,:]   , coords_modes[:,0])
    jacobian_tmp[0,1,inode] = torch.dot(ddeta_i[inode,:]  , coords_modes[:,0])
    jacobian_tmp[0,2,inode] = torch.dot(ddzeta_i[inode,:] , coords_modes[:,0])
    jacobian_tmp[1,0,inode] = torch.dot(ddxi_i[inode,:]   , coords_modes[:,1])
    jacobian_tmp[1,1,inode] = torch.dot(ddeta_i[inode,:]  , coords_modes[:,1])
    jacobian_tmp[1,2,inode] = torch.dot(ddzeta_i[inode,:] , coords_modes[:,1])
    jacobian_tmp[2,0,inode] = torch.dot(ddxi_i[inode,:]   , coords_modes[:,2])
    jacobian_tmp[2,1,inode] = torch.dot(ddeta_i[inode,:]  , coords_modes[:,2])
    jacobian_tmp[2,2,inode] = torch.dot(ddzeta_i[inode,:] , coords_modes[:,2])
    update_progress("Computing Jacobian                 ", inode/(nnodes_if-1))
if coord_sys == 'CYLINDRICAL':
    jacobian = jacobian_tmp.clone()
    scaling_factor = torch.mv(val_i,coords_modes[:,0])
    for inode in range(0,nnodes_if):
        jacobian[1,0,inode] = jacobian_tmp[1,0,inode] * scaling_factor[inode]
        jacobian[1,1,inode] = jacobian_tmp[1,1,inode] * scaling_factor[inode]
        jacobian[1,2,inode] = jacobian_tmp[1,2,inode] * scaling_factor[inode]
else:
    jacobian = jacobian_tmp.clone()



# Matrics and Determinant
metrics = torch.empty(3,3,nnodes_if, dtype=torch.float64) 
jinv    = torch.empty(nnodes_if, dtype=torch.float64) 
for inode in range(0,nnodes_if):
    ijacobian = torch.empty(3,3, dtype=torch.float64)
    imetric   = torch.empty(3,3, dtype=torch.float64)
    for irow in range(0,3):
        for icol in range(0,3):
            ijacobian[irow,icol] = jacobian[irow,icol,inode]
    # Compute jacobian for the ith node
    update_progress("Computing Jinv and Metric          ", inode/(nnodes_if-1))
    jinv[inode] = torch.det(ijacobian) 
    imetric     = torch.inverse(ijacobian)
    for irow in range(0,3):
        for icol in range(0,3):
            metrics[irow,icol,inode] = imetric[irow,icol]

 

# Normals
normals = torch.empty(nnodes_if,ndirs, dtype=torch.float64) 
if iface == 1 or iface == 2:
    for inode in range(0,nnodes_if):
        normals[inode,0] = jinv[inode] * metrics[0,0,inode]
        normals[inode,1] = jinv[inode] * metrics[0,1,inode]
        normals[inode,2] = jinv[inode] * metrics[0,2,inode]
        update_progress("Computing Normals                  ", inode/(nnodes_if-1))
if iface == 3 or iface == 4:
    for inode in range(0,nnodes_if):
        normals[inode,0] = jinv[inode] * metrics[1,0,inode]
        normals[inode,1] = jinv[inode] * metrics[1,1,inode]
        normals[inode,2] = jinv[inode] * metrics[1,2,inode]
        update_progress("Computing Normals                  ", inode/(nnodes_if-1))
if iface == 5 or iface == 6:
    for inode in range(0,nnodes_if):
        normals[inode,0] = jinv[inode] * metrics[2,0,inode]
        normals[inode,1] = jinv[inode] * metrics[2,1,inode]
        normals[inode,2] = jinv[inode] * metrics[2,2,inode]
        update_progress("Computing Normals                  ", inode/(nnodes_if-1))
if iface == 1 or iface == 3 or iface == 5:
    for inode in range(0,nnodes_if):
        normals[inode,0] = -normals[inode,0]
        normals[inode,1] = -normals[inode,1]
        normals[inode,2] = -normals[inode,2]




# Grad1, Grad2, and Grad3
grad1 = torch.empty(nnodes_if,nterms_s, dtype=torch.float64)
grad2 = torch.empty(nnodes_if,nterms_s, dtype=torch.float64)
grad3 = torch.empty(nnodes_if,nterms_s, dtype=torch.float64)
i = 1
for iterm in range(0,nterms_s):
    for inode in range(0,nnodes_if):
        grad1[inode,iterm] = metrics[0,0,inode] * ddxi_i[inode,iterm] + metrics[1,0,inode] * ddeta_i[inode,iterm] + metrics[2,0,inode] * ddzeta_i[inode,iterm] 
        grad2[inode,iterm] = metrics[0,1,inode] * ddxi_i[inode,iterm] + metrics[1,1,inode] * ddeta_i[inode,iterm] + metrics[2,1,inode] * ddzeta_i[inode,iterm] 
        grad3[inode,iterm] = metrics[0,2,inode] * ddxi_i[inode,iterm] + metrics[1,2,inode] * ddeta_i[inode,iterm] + metrics[2,2,inode] * ddzeta_i[inode,iterm] 
        update_progress("Computing grad1, grad2, grad3      ", i/(nnodes_if*nterms_s))
        i += 1







#WRITE_____________________

#
# Metrics
#
f = open("metrics.txt","w")
i = 1
for inode in range (0,nnodes_if):
    f.write("Metric interpolation node %d \n" % (inode+1))
    array = np.zeros([3, 3])
    for irow in range(0,3):
        for icol in range(0,3):
            array[irow,icol] = metrics[irow,icol,inode].item() 
            update_progress("Writing metrics to file            ", i/(nnodes_if*9))
            i += 1
    np.savetxt(f,array)
f.close()

#
# jinv
#
f = open("jinv.txt","w")
array = np.zeros([1])
i = 1
for inode in range (0,nnodes_if):
    f.write("Jinv interpolation node %d \n" % (inode+1))
    array[0] = jinv[inode].item()
    np.savetxt(f,array)
    update_progress("Writing jinv to file               ", i/(nnodes_if))
    i += 1
f.close()

#
# Grad1
#
f = open("grad1.txt","w")
f.write("Grad1 \n")
array = np.zeros([nnodes_if,nterms_s])
i = 1
for inode in range (0,nnodes_if):
    for iterm in range(0,nterms_s):
        array[inode,iterm] = grad1[inode,iterm].item()
        update_progress("Writing grad1 to file              ", i/(nnodes_if*nterms_s))
        i += 1
np.savetxt(f,array)
f.close()

#
# Grad2
#
f = open("grad2.txt","w")
f.write("Grad2 \n")
array = np.zeros([nnodes_if,nterms_s])
i = 1
for inode in range (0,nnodes_if):
    for iterm in range(0,nterms_s):
        array[inode,iterm] = grad2[inode,iterm].item()
        update_progress("Writing grad2 to file              ", i/(nnodes_if*nterms_s))
        i += 1
np.savetxt(f,array)
f.close()

#
# Grad3
#
f = open("grad3.txt","w")
f.write("Grad3 \n")
array = np.zeros([nnodes_if,nterms_s])
i = 1
for inode in range (0,nnodes_if):
    for iterm in range(0,nterms_s):
        array[inode,iterm] = grad3[inode,iterm].item()
        update_progress("Writing grad3 to file              ", i/(nnodes_if*nterms_s))
        i += 1
np.savetxt(f,array)
f.close()


#
# Normals
#
f = open("normals.txt","w")
array = np.zeros([3])
i = 1
for inode in range (0,nnodes_if):
    f.write("Norm interpolation node %d \n" % (inode+1))
    for idir in range(0,ndirs):
        array[idir] = normals[inode,idir].item()
        update_progress("Writing norm to file               ", i/(nnodes_if*ndirs))
        i += 1
    np.savetxt(f,array)
f.close()

#
# dmetric_dx
#
f = open("dmetric_dx.txt","w")
i = 1
for inode in range (0,nnodes_if):
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
                    update_progress("Writing dmetric_dx to file         ", i/(nnodes_if*nnodes_r*ndirs*3*3))
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
array = np.zeros([nnodes_if,nnodes_r*ndirs])
for inode in range (0,nnodes_if):
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
            update_progress("Writing interp_xcoords_dx to file  ", i/(nnodes_if*nnodes_r*3))
            i += 1
    # This avoid to accumulate derivatives
    dummy = coords.grad.data.zero_()
np.savetxt(f,array)
f.close()
f = open("dinterp_ycoords_dx.txt","w")
i = 1
f.write("ycoord interpolation, coord 2,  row=node, col=nnodes_r*dir \n")
array = np.zeros([nnodes_if,nnodes_r*ndirs])
for inode in range (0,nnodes_if):
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
            update_progress("Writing interp_ycoords_dx to file  ", i/(nnodes_if*nnodes_r*3))
            i += 1
    # This avoid to accumulate derivatives
    dummy = coords.grad.data.zero_()
np.savetxt(f,array)
f.close()
f = open("dinterp_zcoords_dx.txt","w")
i = 1
f.write("zcoord interpolation, coord 3,  row=node, col=nnodes_r*dir \n")
array = np.zeros([nnodes_if,nnodes_r*ndirs])
for inode in range (0,nnodes_if):
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
            update_progress("Writing interp_zcoords_dx to file  ", i/(nnodes_if*nnodes_r*3))
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
for inode in range (0,nnodes_if):
    array = np.zeros([nnodes_r,ndirs])
    f.write("djinv_dx interpolation node %s, row=inode_diff, col=dir \n" % (inode+1))
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
                data = jinv[inode]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[inode_diff,idir] = ddata_np[inode_diff,idir]
                update_progress("Writing djinv_dx to file           ", i/(nnodes_if*nnodes_r*ndirs))
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
        array = np.zeros([nnodes_if,nterms_s])
        for irow in range(0,nnodes_if):
            for icol in range(0,nterms_s):
                data = grad1[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dgrad1_dx to file          ", i/(nnodes_if*nnodes_r*ndirs*nterms_s))
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
        array = np.zeros([nnodes_if,nterms_s])
        for irow in range(0,nnodes_if):
            for icol in range(0,nterms_s):
                data = grad2[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dgrad2_dx to file          ", i/(nnodes_if*nnodes_r*ndirs*nterms_s))
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
        array = np.zeros([nnodes_if,nterms_s])
        for irow in range(0,nnodes_if):
            for icol in range(0,nterms_s):
                data = grad3[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dgrad3_dx to file          ", i/(nnodes_if*nnodes_r*ndirs*nterms_s))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()






#
# dnorm_dx
#
f = open("dnorm_dx.txt","w")
array = np.zeros([3])
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dnorm_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = np.zeros([nnodes_if,ndirs])
        for irow in range(0,nnodes_if):
            for icol in range(0,ndirs):
                data = normals[irow,icol]
                data.backward(retain_graph=True)
                ddata = coords.grad
                ddata_np = ddata.numpy()
                array[irow,icol] = ddata_np[inode_diff,idir]
                update_progress("Writing dnorm_dx to file           ", i/(nnodes_r*ndirs*nnodes_if*ndirs))
                dummy = coords.grad.data.zero_()
                i += 1
        np.savetxt(f,array)
f.close()




