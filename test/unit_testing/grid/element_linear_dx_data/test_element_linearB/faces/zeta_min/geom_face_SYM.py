from __future__ import division

import sys
import os
import time
import numpy
import pickle
from sympy import *


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
print "WARNING: This script is very slow, it might run for hours. It is strongly recommended to watch Netflix in the meanwhile."
################################################################################################################

# Define symbols for each coordinate for support node
x1,y1,z1 = symbols('x1 y1 z1') 
x2,y2,z2 = symbols('x2 y2 z2') 
x3,y3,z3 = symbols('x3 y3 z3') 
x4,y4,z4 = symbols('x4 y4 z4') 
x5,y5,z5 = symbols('x5 y5 z5') 
x6,y6,z6 = symbols('x6 y6 z6') 
x7,y7,z7 = symbols('x7 y7 z7') 
x8,y8,z8 = symbols('x8 y8 z8') 

coords_ = Matrix( [[x1,y1,z1],
                   [x2,y2,z2],
                   [x3,y3,z3],
                   [x4,y4,z4],
                   [x5,y5,z5],
                   [x6,y6,z6],
                   [x7,y7,z7],
                   [x8,y8,z8],
                   ] )


nnodes_r  = coords_.shape[0]
nnodes_if = 4 
nterms_s  = 8
ndirs     = 3
coord_sys = 'CARTESIAN'

# Define face to test 
#XI_MIN    = 1
#XI_MAX    = 2
#ETA_MIN   = 3
#ETA_MAX   = 4
#ZETA_MIN  = 5
#ZETA_MAX  = 6
iface = 5 

# Define coordinate values at support nodes
coords  = Matrix( [[0.0,0.0,0.0],
                   [5.0,0.0,0.0],
                   [0.0,1.0,0.0],
                   [5.0,1.0,0.0],
                   [0.0,0.0,1.0],
                   [5.0,0.0,1.0],
                   [0.0,1.0,1.0],
                   [5.0,1.0,1.0],
                   ] )

# Define matrix of polynomial basis terms at support nodes
val_r   = Matrix( [[ 1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0,-1.0],
                   [ 1.0,-1.0,-1.0, 1.0,-1.0,-1.0, 1.0, 1.0],
                   [ 1.0, 1.0,-1.0,-1.0,-1.0, 1.0,-1.0, 1.0],
                   [ 1.0, 1.0,-1.0, 1.0, 1.0,-1.0,-1.0,-1.0],
                   [ 1.0,-1.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0],
                   [ 1.0,-1.0, 1.0, 1.0,-1.0, 1.0,-1.0,-1.0],
                   [ 1.0, 1.0, 1.0,-1.0,-1.0,-1.0, 1.0,-1.0],
                   [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                   ] )


# Define matrices at interpolation nodes (quadrature, level = 1)
if iface == 1:
    val_i   = Matrix( [[ 1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
                       [ 1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
                       ] )

    ddxi_i  = Matrix( [[ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0, 1.0/3.0],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0,-1.0/3.0],
                       [ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),0.0,-1.0/3.0],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),0.0, 1.0/3.0],
                       ] )

    ddeta_i = Matrix( [[ 0.0,1.0,0.0,0.0,-1.0,0.0,-sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,-1.0,0.0,-sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,-1.0,0.0, sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,-1.0,0.0, sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       ] )

    ddzeta_i= Matrix( [[ 0.0,0.0,1.0,0.0,0.0,-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       ] )

if iface == 2:
    val_i   = Matrix( [[ 1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                       [ 1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                       ] )

    ddxi_i  = Matrix( [[ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0, 1.0/3.0],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0,-1.0/3.0],
                       [ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),0.0,-1.0/3.0],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),0.0, 1.0/3.0],
                       ] )

    ddeta_i = Matrix( [[ 0.0,1.0,0.0,0.0,1.0,0.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,1.0,0.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,1.0,0.0, sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,1.0,0.0, sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       ] )

    ddzeta_i= Matrix( [[ 0.0,0.0,1.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0)],
                       ] )

if iface == 3:
    val_i   = Matrix( [[ 1.0,-1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                       [ 1.0,-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                       [ 1.0,-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                       [ 1.0,-1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                       ] )

    ddxi_i  = Matrix( [[ 0.0,0.0,0.0,1.0,-1.0,-sqrt(1.0/3.0),0.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,-1.0,-sqrt(1.0/3.0),0.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,-1.0, sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,-1.0, sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0)],
                       ] )

    ddeta_i = Matrix( [[ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0), 1.0/3.0],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0, sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0, sqrt(1.0/3.0), 1.0/3.0],
                       ] )

    ddzeta_i= Matrix( [[ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0)],
                       ] )

if iface == 4:
    val_i   = Matrix( [[ 1.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                       [ 1.0,1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                       [ 1.0,1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                       [ 1.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                       ] )

    ddxi_i  = Matrix( [[ 0.0,0.0,0.0,1.0,1.0,-sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,1.0,-sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,1.0, sqrt(1.0/3.0),0.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,1.0, sqrt(1.0/3.0),0.0, sqrt(1.0/3.0)],
                       ] )

    ddeta_i = Matrix( [[ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0), 1.0/3.0],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0, sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0, sqrt(1.0/3.0), 1.0/3.0],
                       ] )

    ddzeta_i= Matrix( [[ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),1.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),1.0, sqrt(1.0/3.0)],
                       ] )

if iface ==5:
    val_i   = Matrix( [[ 1.0,-sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                       [ 1.0,-sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                       ] )

    ddxi_i  = Matrix( [[ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-1.0,0.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-1.0,0.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),-1.0,0.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),-1.0,0.0,-sqrt(1.0/3.0)],
                       ] )

    ddeta_i = Matrix( [[ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,-1.0, sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,-1.0,-sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,-1.0, sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,-1.0,-sqrt(1.0/3.0)],
                       ] )

    ddzeta_i= Matrix( [[ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0],
                       ] )

if iface ==6:
    val_i   = Matrix( [[ 1.0,-sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                       [ 1.0,-sqrt(1.0/3.0),1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                       [ 1.0, sqrt(1.0/3.0),1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                       ] )

    ddxi_i  = Matrix( [[ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),1.0,0.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),1.0,0.0,-sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),1.0,0.0, sqrt(1.0/3.0)],
                       [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),1.0,0.0, sqrt(1.0/3.0)],
                       ] )

    ddeta_i = Matrix( [[ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,1.0,-sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,1.0, sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,1.0,-sqrt(1.0/3.0)],
                       [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,1.0, sqrt(1.0/3.0)],
                       ] )

    ddzeta_i= Matrix( [[ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0],
                       [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0],
                       ] )


#--------------------------------------------------------------------

# Matrix modes_to_nodes
val_r_inv = val_r**(-1)


# Computes coordiantes modes
coords_modes_ = val_r_inv * coords_
coords_modes  = lambdify(coords_,coords_modes_,"numpy")


# Initialized jacobian
jacobian_ = MutableSparseNDimArray.zeros(3, 3, nnodes_if)
for inode in range(0,nnodes_if):
    jacobian_[0,0,inode] = ddxi_i[inode,:]   * coords_modes_[:,0]
    jacobian_[0,1,inode] = ddeta_i[inode,:]  * coords_modes_[:,0]
    jacobian_[0,2,inode] = ddzeta_i[inode,:] * coords_modes_[:,0]
    jacobian_[1,0,inode] = ddxi_i[inode,:]   * coords_modes_[:,1]
    jacobian_[1,1,inode] = ddeta_i[inode,:]  * coords_modes_[:,1]
    jacobian_[1,2,inode] = ddzeta_i[inode,:] * coords_modes_[:,1]
    jacobian_[2,0,inode] = ddxi_i[inode,:]   * coords_modes_[:,2]
    jacobian_[2,1,inode] = ddeta_i[inode,:]  * coords_modes_[:,2]
    jacobian_[2,2,inode] = ddzeta_i[inode,:] * coords_modes_[:,2]
    update_progress("Computing Jacobian                  ", inode/(nnodes_if-1))

if coord_sys == 'CYLINDRICAL':
    scaling_factor = val_i * coords_modes_[:,0]
    for inode in range(0,nnodes_if):
        jacobian_[1,0,inode] = jacobian_[1,0,inode] * scaling_factor[inode]
        jacobian_[1,1,inode] = jacobian_[1,1,inode] * scaling_factor[inode]
        jacobian_[1,2,inode] = jacobian_[1,2,inode] * scaling_factor[inode]




# Matrics and Determinant
metrics_ = MutableSparseNDimArray.zeros(3, 3, nnodes_if)
jinv_    = zeros(nnodes_if)
for inode in range(0,nnodes_if):
    ijacobian = zeros(3,3)
    for irow in range(0,3):
        for icol in range(0,3):
            ijacobian[irow,icol] = jacobian_[irow,icol,inode]
    # Compute jacobian for the ith node
    update_progress("Computing Jinv and Metric           ", inode/(nnodes_if-1))
    jinv_[inode] = ijacobian.det() 
    imetric      = ijacobian**(-1)
    for irow in range(0,3):
        for icol in range(0,3):
            metrics_[irow,icol,inode] = imetric[irow,icol]
    


# Normals
normals_ = MutableSparseNDimArray.zeros(nnodes_if,ndirs)
if iface == 1 or iface == 2:
    for inode in range(0,nnodes_if):
        normals_[inode,0] = jinv_[inode] * metrics_[0,0,inode]
        normals_[inode,1] = jinv_[inode] * metrics_[0,1,inode]
        normals_[inode,2] = jinv_[inode] * metrics_[0,2,inode]
        update_progress("Computing Normals                   ", inode/(nnodes_if-1))
if iface == 3 or iface == 4:
    for inode in range(0,nnodes_if):
        normals_[inode,0] = jinv_[inode] * metrics_[1,0,inode]
        normals_[inode,1] = jinv_[inode] * metrics_[1,1,inode]
        normals_[inode,2] = jinv_[inode] * metrics_[1,2,inode]
        update_progress("Computing Normals                   ", inode/(nnodes_if-1))
if iface == 5 or iface == 6:
    for inode in range(0,nnodes_if):
        normals_[inode,0] = jinv_[inode] * metrics_[2,0,inode]
        normals_[inode,1] = jinv_[inode] * metrics_[2,1,inode]
        normals_[inode,2] = jinv_[inode] * metrics_[2,2,inode]
        update_progress("Computing Normals                   ", inode/(nnodes_if-1))
if iface == 1 or iface == 3 or iface == 5:
    for inode in range(0,nnodes_if):
        normals_[inode,0] = -normals_[inode,0]
        normals_[inode,1] = -normals_[inode,1]
        normals_[inode,2] = -normals_[inode,2]




# Grad1, Grad2, and Grad3
grad1_ = zeros(nnodes_if,nterms_s)
grad2_ = zeros(nnodes_if,nterms_s)
grad3_ = zeros(nnodes_if,nterms_s)
i = 1
for iterm in range(0,nterms_s):
    for inode in range(0,nnodes_if):
        grad1_[inode,iterm] = metrics_[0,0,inode] * ddxi_i[inode,iterm] + metrics_[1,0,inode] * ddeta_i[inode,iterm] + metrics_[2,0,inode] * ddzeta_i[inode,iterm] 
        grad2_[inode,iterm] = metrics_[0,1,inode] * ddxi_i[inode,iterm] + metrics_[1,1,inode] * ddeta_i[inode,iterm] + metrics_[2,1,inode] * ddzeta_i[inode,iterm] 
        grad3_[inode,iterm] = metrics_[0,2,inode] * ddxi_i[inode,iterm] + metrics_[1,2,inode] * ddeta_i[inode,iterm] + metrics_[2,2,inode] * ddzeta_i[inode,iterm] 
        update_progress("Computing grad1, grad2, grad3       ", i/(nnodes_if*nterms_s))
        i += 1




# Differentiate determinant
djinv_dx_ = MutableSparseNDimArray.zeros(nnodes_if, nnodes_r, ndirs) 
i = 1
for inode in range(0,nnodes_if):
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            djinv_dx_[inode,inode_diff,idir] = jinv_[inode].diff(coords_[inode_diff,idir])
            update_progress("Computing djinv_dx                  ", i/(nnodes_if*nnodes_r*ndirs))
            i += 1


# Differentiate metrics
dmetric_dx_ = MutableSparseNDimArray.zeros(3,3,nnodes_if,nnodes_r,ndirs) 
i = 1
for inode in range(0,nnodes_if):
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            for irow in range(0,3):
                for icol in range(0,3):
                    dmetric_dx_[irow,icol,inode,inode_diff,idir] = metrics_[irow,icol,inode].diff(coords_[inode_diff,idir])
                    update_progress("Computing dmetric_dx                ", i/(nnodes_if*nnodes_r*ndirs*9))
                    i += 1


# Differentaite Gradients
dgrad1_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nterms_s,nnodes_r,ndirs)
dgrad2_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nterms_s,nnodes_r,ndirs)
dgrad3_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nterms_s,nnodes_r,ndirs)
i = 1
for inode in range(0,nnodes_if):
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            for inode in range(0,nnodes_if):
                for iterm in range(0,nterms_s):
                    dgrad1_dx_[inode,iterm,inode_diff,idir] = grad1_[inode,iterm].diff(coords_[inode_diff,idir])
                    dgrad2_dx_[inode,iterm,inode_diff,idir] = grad2_[inode,iterm].diff(coords_[inode_diff,idir])
                    dgrad3_dx_[inode,iterm,inode_diff,idir] = grad3_[inode,iterm].diff(coords_[inode_diff,idir])
                    update_progress("Computing dgrad1_dx, dgrad2_dx, ..  ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if*nterms_s))
                    i += 1


# Differentiate Normals
dnorm_dx_ =  MutableSparseNDimArray.zeros(nnodes_if,ndirs,nnodes_r,ndirs)
i = 1
for inode in range(0,nnodes_if):
    for idir in range(0,ndirs):
        for inode_diff in range(0,nnodes_r):
            for idir_diff in range(0,ndirs):
                dnorm_dx_[inode,idir,inode_diff,idir_diff] = normals_[inode,idir].diff(coords_[inode_diff,idir_diff])
                update_progress("Computing dnorm_dx                  ", i/(nnodes_if*ndirs*nnodes_r*ndirs))
                i += 1



#WRITE_____________________

#
# Metrics
#
f = open("metrics.txt","w")
i = 1
for inode in range (0,nnodes_if):
    f.write("Metric interpolation node %d \n" % (inode+1))
    array = numpy.zeros([3, 3])
    for irow in range(0,3):
        for icol in range(0,3):
            data_sym = lambdify(coords_,metrics_[irow,icol,inode],"numpy")
            data_val = data_sym(*flatten(coords))
            array[irow,icol] = data_val
            update_progress("Writing metrics to file             ", i/(nnodes_if*9))
            i += 1
    numpy.savetxt(f,array)
f.close()

#
# jinv
#
f = open("jinv.txt","w")
array = numpy.zeros([1])
i = 1
for inode in range (0,nnodes_if):
    f.write("Jinv interpolation node %d \n" % (inode+1))
    data_sym = lambdify(coords_,jinv_[inode],"numpy")
    data_val = data_sym(*flatten(coords))
    array[0] = data_val
    numpy.savetxt(f,array)
    update_progress("Writing jinv to file                ", i/(nnodes_if))
    i += 1
f.close()

#
# Grad1
#
f = open("grad1.txt","w")
f.write("Grad1 \n")
array = numpy.zeros([nnodes_if,nterms_s])
i = 1
for inode in range (0,nnodes_if):
    for iterm in range(0,nterms_s):
        data_sym  = lambdify(coords_,grad1_[inode,iterm],"numpy")
        data_val = (data_sym(*flatten(coords)))
        array[inode,iterm] = data_val
        update_progress("Writing grad1 to file               ", i/(nnodes_if*nterms_s))
        i += 1
numpy.savetxt(f,array)
f.close()

#
# Grad2
#
f = open("grad2.txt","w")
f.write("Grad2\n")
array = numpy.zeros([nnodes_if,nterms_s])
i = 1
for inode in range (0,nnodes_if):
    for iterm in range(0,nterms_s):
        data_sym  = lambdify(coords_,grad2_[inode,iterm],"numpy")
        data_val = (data_sym(*flatten(coords)))
        array[inode,iterm] = data_val
        update_progress("Writing grad2 to file               ", i/(nnodes_if*nterms_s))
        i += 1
numpy.savetxt(f,array)
f.close()

#
# Grad3
#
f = open("grad3.txt","w")
f.write("Grad3\n")
array = numpy.zeros([nnodes_if,nterms_s])
i = 1
for inode in range (0,nnodes_if):
    for iterm in range(0,nterms_s):
        data_sym  = lambdify(coords_,grad3_[inode,iterm],"numpy")
        data_val = (data_sym(*flatten(coords)))
        array[inode,iterm] = data_val
        update_progress("Writing grad3 to file               ", i/(nnodes_if*nterms_s))
        i += 1
numpy.savetxt(f,array)
f.close()


#
# Normals
#
f = open("normals.txt","w")
array = numpy.zeros([3])
i = 1
for inode in range (0,nnodes_if):
    f.write("Norm interpolation node %d \n" % (inode+1))
    for idir in range(0,ndirs):
        data_sym = lambdify(coords_,normals_[inode,idir],"numpy")
        data_val = data_sym(*flatten(coords))
        array[idir] = data_val
        update_progress("Writing norm to file                ", i/(nnodes_if*ndirs))
        i += 1
    numpy.savetxt(f,array)
f.close()




#
# dmetric_dx
#
f = open("dmetric_dx.txt","w")
i = 1
for inode in range (0,nnodes_if):
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            array = numpy.zeros([3,3])
            f.write("dmetric_dx interpolation node %s, diff_node %s, diff_dir %s \n" % (inode+1,inode_diff+1,idir+1))
            for irow in range(0,3):
                for icol in range(0,3):
                    data_sym = lambdify(coords_,dmetric_dx_[irow,icol,inode,inode_diff,idir],"numpy")
                    data_val = data_sym(*flatten(coords))
                    array[irow,icol] = data_val
                    update_progress("Writing dmetric_dx to file          ", i/(nnodes_if*nnodes_r*ndirs*3*3))
                    i += 1
            numpy.savetxt(f,array)
f.close()

#
# djinv_dx
#
f = open("djinv_dx.txt","w")
i = 1
for inode in range (0,nnodes_if):
    array = numpy.zeros([nnodes_r,ndirs])
    f.write("djinv_dx interpolation node %s, row=inode_diff, col=dir \n" % (inode+1))
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            data_sym = lambdify(coords_,djinv_dx_[inode,inode_diff,idir],"numpy")
            data_val = data_sym(*flatten(coords))
            array[inode_diff,idir] = data_val
            update_progress("Writing djinv_dx to file            ", i/(nnodes_if*nnodes_r*ndirs))
            i += 1
    numpy.savetxt(f,array)
f.close()

#
# dgrad1_dx
#
f = open("dgrad1_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dgrad1_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nterms_s])
        for irow in range(0,nnodes_if):
            for icol in range(0,nterms_s):
                data_sym  = lambdify(coords_,dgrad1_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dgrad1_dx to file           ", i/(nnodes_if*nnodes_r*ndirs*nterms_s))
                i += 1
        numpy.savetxt(f,array)
f.close()

#
# dgrad2_dx
#
f = open("dgrad2_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dgrad2_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nterms_s])
        for irow in range(0,nnodes_if):
            for icol in range(0,nterms_s):
                data_sym  = lambdify(coords_,dgrad2_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dgrad2_dx to file           ", i/(nnodes_if*nnodes_r*ndirs*nterms_s))
                i += 1
        numpy.savetxt(f,array)
f.close()

#
# dgrad3_dx
#
f = open("dgrad3_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dgrad3_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nterms_s])
        for irow in range(0,nnodes_if):
            for icol in range(0,nterms_s):
                data_sym  = lambdify(coords_,dgrad3_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dgrad3_dx to file           ", i/(nnodes_if*nnodes_r*ndirs*nterms_s))
                i += 1
        numpy.savetxt(f,array)
f.close()






#
# Normals
#
f = open("dnorm_dx.txt","w")
array = numpy.zeros([3])
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dnorm_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,ndirs])
        for irow in range(0,nnodes_if):
            for icol in range(0,ndirs):
                data_sym  = lambdify(coords_,dnorm_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dnorm_dx to file            ", i/(nnodes_r*ndirs*nnodes_if*ndirs))
                i += 1
        numpy.savetxt(f,array)
f.close()






