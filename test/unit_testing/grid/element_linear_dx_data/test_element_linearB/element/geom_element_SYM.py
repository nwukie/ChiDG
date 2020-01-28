from __future__ import division

import sys
import os
import time
import numpy
import pickle
from sympy import *
from sympy.tensor.array import MutableSparseNDimArray


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
nnodes_ie = 8 
nnodes_if = 4 
nterms_s  = 8
ndirs     = 3



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
val_i   = Matrix( [[ 1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0, 1.0/3.0,-1.0/3.0*sqrt(1.0/3.0)],
                   [ 1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0, 1.0/3.0, 1.0/3.0*sqrt(1.0/3.0)],
                   [ 1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0,-1.0/3.0, 1.0/3.0*sqrt(1.0/3.0)],
                   [ 1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0,-1.0/3.0,-1.0/3.0*sqrt(1.0/3.0)],
                   [ 1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0,-1.0/3.0, 1.0/3.0*sqrt(1.0/3.0)],
                   [ 1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0,-1.0/3.0,-1.0/3.0*sqrt(1.0/3.0)],
                   [ 1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0, 1.0/3.0,-1.0/3.0*sqrt(1.0/3.0)],
                   [ 1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0*sqrt(1.0/3.0)],
                   ] )

ddxi_i  = Matrix( [[ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0, 1.0/3.0],
                   [ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0, 1.0/3.0],
                   [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0,-1.0/3.0],
                   [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),0.0,-1.0/3.0],
                   [ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),0.0,-1.0/3.0],
                   [ 0.0,0.0,0.0,1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),0.0,-1.0/3.0],
                   [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),0.0, 1.0/3.0],
                   [ 0.0,0.0,0.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),0.0, 1.0/3.0],
                   ] )

ddeta_i = Matrix( [[ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0), 1.0/3.0],
                   [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0), 1.0/3.0],
                   [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0,-sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0, sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0, sqrt(1.0/3.0), 1.0/3.0],
                   [ 0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),0.0, sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,1.0,0.0,0.0, sqrt(1.0/3.0),0.0, sqrt(1.0/3.0), 1.0/3.0],
                   ] )

ddzeta_i= Matrix( [[ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0],
                   [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0],
                   [ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0],
                   [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,0.0,1.0,0.0,0.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0],
                   [ 0.0,0.0,1.0,0.0,0.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0],
                   ] )

# Define element interpolation nodes weights for linear element
weights_e = Matrix( [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0] )


# Define val_f for each face
# Face 1, XI_MIN
val_1   = Matrix( [[ 1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
                   [ 1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, 1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-1.0/3.0],
                   ] )

# Face 2, XI_MAX
val_2   = Matrix( [[ 1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                   [ 1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, 1.0/3.0],
                   ] ) 

# Face 3, ETA_MIN
val_3   = Matrix( [[ 1.0,-1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                   [ 1.0,-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                   [ 1.0,-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                   [ 1.0,-1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                   ] ) 

# Face 4, ETA_MAX
val_4   = Matrix( [[ 1.0,1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                   [ 1.0,1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                   [ 1.0,1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                   [ 1.0,1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                   ] )

# Face 5, ZETA_MIN
val_5   = Matrix( [[ 1.0,-sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                   [ 1.0,-sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0),-1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0),-1.0, sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                   ] )

# Face 6, ZETA_MAX
val_6   = Matrix( [[ 1.0,-sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0), sqrt(1.0/3.0),-1.0/3.0,-sqrt(1.0/3.0), 1.0/3.0],
                   [ 1.0,-sqrt(1.0/3.0),1.0, sqrt(1.0/3.0),-sqrt(1.0/3.0), 1.0/3.0,-sqrt(1.0/3.0),-1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0),1.0,-sqrt(1.0/3.0),-sqrt(1.0/3.0),-1.0/3.0, sqrt(1.0/3.0),-1.0/3.0],
                   [ 1.0, sqrt(1.0/3.0),1.0, sqrt(1.0/3.0), sqrt(1.0/3.0), 1.0/3.0, sqrt(1.0/3.0), 1.0/3.0],
                   ] )





#--------------------------------------------------------------------

# Matrix modes_to_nodes
val_r_inv = val_r**(-1)


# Computes coordiantes modes
coords_modes_ = val_r_inv * coords_
coords_modes  = lambdify(coords_,coords_modes_,"numpy")


# Initialized coordiantes
interp_coords_ = MutableSparseNDimArray.zeros(nnodes_ie,3)
for inode in range(0,nnodes_ie):
    for idir in range(0,3):
        interp_coords_[inode,idir] = val_i[inode,:] * coords_modes_[:,idir]


# Initialized jacobian
jacobian_ = MutableSparseNDimArray.zeros(3, 3, nnodes_ie)
for inode in range(0,nnodes_ie):
    jacobian_[0,0,inode] = ddxi_i[inode,:]   * coords_modes_[:,0]
    jacobian_[0,1,inode] = ddeta_i[inode,:]  * coords_modes_[:,0]
    jacobian_[0,2,inode] = ddzeta_i[inode,:] * coords_modes_[:,0]
    jacobian_[1,0,inode] = ddxi_i[inode,:]   * coords_modes_[:,1]
    jacobian_[1,1,inode] = ddeta_i[inode,:]  * coords_modes_[:,1]
    jacobian_[1,2,inode] = ddzeta_i[inode,:] * coords_modes_[:,1]
    jacobian_[2,0,inode] = ddxi_i[inode,:]   * coords_modes_[:,2]
    jacobian_[2,1,inode] = ddeta_i[inode,:]  * coords_modes_[:,2]
    jacobian_[2,2,inode] = ddzeta_i[inode,:] * coords_modes_[:,2]
    update_progress("Computing Jacobian             ", inode/(nnodes_ie-1))


# Matrics and Determinant
metrics_ = MutableSparseNDimArray.zeros(3, 3, nnodes_ie)
jinv_    = zeros(nnodes_ie)
for inode in range(0,nnodes_ie):
    ijacobian = zeros(3,3)
    for irow in range(0,3):
        for icol in range(0,3):
            ijacobian[irow,icol] = jacobian_[irow,icol,inode]
    # Compute jacobian for the ith node
    update_progress("Computing Jinv and Metric      ", inode/(nnodes_ie-1))
    jinv_[inode] = ijacobian.det() 
    imetric      = ijacobian**(-1)
    for irow in range(0,3):
        for icol in range(0,3):
            metrics_[irow,icol,inode] = imetric[irow,icol]

 
# Compute inverse Mass matrix
invmass_ = zeros(nterms_s,nterms_s)
mass_    = zeros(nterms_s,nterms_s)
i = 1
val_tmp = val_i
for iterm in range(0,nterms_s):
    for inode in range(0,nnodes_ie):
        val_tmp[inode,iterm] = val_tmp[inode,iterm] * weights_e[inode] * jinv_[inode]
        update_progress("Computing invmass              ", i/(nterms_s*nnodes_ie))
        i += 1
mass_    = transpose(val_tmp)*val_i
invmass_ = (mass_)**(-1)


# Compute BR2_VOL for each face
br2_vol_face1_ = zeros(nnodes_ie,nnodes_if)
br2_vol_face2_ = zeros(nnodes_ie,nnodes_if)
br2_vol_face3_ = zeros(nnodes_ie,nnodes_if)
br2_vol_face4_ = zeros(nnodes_ie,nnodes_if)
br2_vol_face5_ = zeros(nnodes_ie,nnodes_if)
br2_vol_face6_ = zeros(nnodes_ie,nnodes_if)
br2_vol_face1_ = val_i*(invmass_*transpose(val_1))
br2_vol_face2_ = val_i*(invmass_*transpose(val_2))
br2_vol_face3_ = val_i*(invmass_*transpose(val_3))
br2_vol_face4_ = val_i*(invmass_*transpose(val_4))
br2_vol_face5_ = val_i*(invmass_*transpose(val_5))
br2_vol_face6_ = val_i*(invmass_*transpose(val_6))
update_progress("Computing br2_vol              ", 1)


# Compute BR2_FACE for each face
br2_face_face1_ = zeros(nnodes_if,nnodes_if)
br2_face_face2_ = zeros(nnodes_if,nnodes_if)
br2_face_face3_ = zeros(nnodes_if,nnodes_if)
br2_face_face4_ = zeros(nnodes_if,nnodes_if)
br2_face_face5_ = zeros(nnodes_if,nnodes_if)
br2_face_face6_ = zeros(nnodes_if,nnodes_if)
br2_face_face1_ = val_1*(invmass_*transpose(val_1))
br2_face_face2_ = val_2*(invmass_*transpose(val_2))
br2_face_face3_ = val_3*(invmass_*transpose(val_3))
br2_face_face4_ = val_4*(invmass_*transpose(val_4))
br2_face_face5_ = val_5*(invmass_*transpose(val_5))
br2_face_face6_ = val_6*(invmass_*transpose(val_6))
update_progress("Computing br2_face             ", 1)


## Grad1, Grad2, and Grad3
#grad1_ = zeros(nnodes_ie,nterms_s)
#grad2_ = zeros(nnodes_ie,nterms_s)
#grad3_ = zeros(nnodes_ie,nterms_s)
#i = 1
#for iterm in range(0,nterms_s):
#    for inode in range(0,nnodes_ie):
#        grad1_[inode,iterm] = metrics_[0,0,inode] * ddxi_i[inode,iterm] + metrics_[1,0,inode] * ddeta_i[inode,iterm] + metrics_[2,0,inode] * ddzeta_i[inode,iterm] 
#        grad2_[inode,iterm] = metrics_[0,1,inode] * ddxi_i[inode,iterm] + metrics_[1,1,inode] * ddeta_i[inode,iterm] + metrics_[2,1,inode] * ddzeta_i[inode,iterm] 
#        grad3_[inode,iterm] = metrics_[0,2,inode] * ddxi_i[inode,iterm] + metrics_[1,2,inode] * ddeta_i[inode,iterm] + metrics_[2,2,inode] * ddzeta_i[inode,iterm] 
#        update_progress("Computing grad1, grad2, grad3       ", i/(nnodes_ie*nterms_s))
#        i += 1



# Differentiate coordinates at interpolation points
interp_coords_dx_ = MutableSparseNDimArray.zeros(nnodes_ie, 3, nnodes_r, ndirs)
i = 1
for inode in range(0,nnodes_ie):
    for direct in range(0,3):
        for inode_diff in range(0,nnodes_r):
            for idir in range(0,ndirs):
                interp_coords_dx_[inode,direct,inode_diff,idir] = interp_coords_[inode,direct].diff(coords_[inode_diff,idir])
                update_progress("Computing interp_coords_dx          ", i/(nnodes_ie*nnodes_r*ndirs*3))
                i += 1



# Differentiate determinant
djinv_dx_ = MutableSparseNDimArray.zeros(nnodes_ie, nnodes_r, ndirs) 
i = 1
for inode in range(0,nnodes_ie):
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            djinv_dx_[inode,inode_diff,idir] = jinv_[inode].diff(coords_[inode_diff,idir])
            update_progress("Computing djinv_dx                  ", i/(nnodes_ie*nnodes_r*ndirs))
            i += 1


# Differentiate metrics
dmetric_dx_ = MutableSparseNDimArray.zeros(3,3,nnodes_ie,nnodes_r,ndirs) 
i = 1
for inode in range(0,nnodes_ie):
    for inode_diff in range(0,nnodes_r):
        for idir in range(0,ndirs):
            for irow in range(0,3):
                for icol in range(0,3):
                    dmetric_dx_[irow,icol,inode,inode_diff,idir] = metrics_[irow,icol,inode].diff(coords_[inode_diff,idir])
                    update_progress("Computing dmetric_dx                ", i/(nnodes_ie*nnodes_r*ndirs*9))
                    i += 1


# Differentiate invmass
dinvmass_dx_ = MutableSparseNDimArray.zeros(nterms_s,nterms_s,nnodes_r,ndirs) 
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        for irow in range(0,nterms_s):
            for icol in range(0,nterms_s):
                dinvmass_dx_[irow,icol,inode_diff,idir] = invmass_[irow,icol].diff(coords_[inode_diff,idir])
                update_progress("Computing dinvmass_dx               ", i/(nnodes_r*ndirs*nterms_s*nterms_s))
                i += 1


# Differentiate BR2_vol
dbr2_vol_face1_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nnodes_if,nnodes_r,ndirs) 
dbr2_vol_face2_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nnodes_if,nnodes_r,ndirs) 
dbr2_vol_face3_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nnodes_if,nnodes_r,ndirs) 
dbr2_vol_face4_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nnodes_if,nnodes_r,ndirs) 
dbr2_vol_face5_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nnodes_if,nnodes_r,ndirs) 
dbr2_vol_face6_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nnodes_if,nnodes_r,ndirs) 
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                dbr2_vol_face1_dx_[irow,icol,inode_diff,idir] = br2_vol_face1_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_vol_face2_dx_[irow,icol,inode_diff,idir] = br2_vol_face2_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_vol_face3_dx_[irow,icol,inode_diff,idir] = br2_vol_face3_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_vol_face4_dx_[irow,icol,inode_diff,idir] = br2_vol_face4_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_vol_face5_dx_[irow,icol,inode_diff,idir] = br2_vol_face5_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_vol_face6_dx_[irow,icol,inode_diff,idir] = br2_vol_face6_[irow,icol].diff(coords_[inode_diff,idir])
                update_progress("Computing dbr2_vol_faces_dx         ", i/(nnodes_r*ndirs*nnodes_ie*nnodes_if))
                i += 1

# Differentiate BR2_face
dbr2_face_face1_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nnodes_if,nnodes_r,ndirs) 
dbr2_face_face2_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nnodes_if,nnodes_r,ndirs) 
dbr2_face_face3_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nnodes_if,nnodes_r,ndirs) 
dbr2_face_face4_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nnodes_if,nnodes_r,ndirs) 
dbr2_face_face5_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nnodes_if,nnodes_r,ndirs) 
dbr2_face_face6_dx_ = MutableSparseNDimArray.zeros(nnodes_if,nnodes_if,nnodes_r,ndirs) 
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                dbr2_face_face1_dx_[irow,icol,inode_diff,idir] = br2_face_face1_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_face_face2_dx_[irow,icol,inode_diff,idir] = br2_face_face2_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_face_face3_dx_[irow,icol,inode_diff,idir] = br2_face_face3_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_face_face4_dx_[irow,icol,inode_diff,idir] = br2_face_face4_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_face_face5_dx_[irow,icol,inode_diff,idir] = br2_face_face5_[irow,icol].diff(coords_[inode_diff,idir])
                dbr2_face_face6_dx_[irow,icol,inode_diff,idir] = br2_face_face6_[irow,icol].diff(coords_[inode_diff,idir])
                update_progress("Computing dbr2_face_faces_dx        ", i/(nnodes_r*ndirs*nnodes_if*nnodes_if))
                i += 1



## Differentaite Gradients
#dgrad1_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nterms_s,nnodes_r,ndirs)
#dgrad2_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nterms_s,nnodes_r,ndirs)
#dgrad3_dx_ = MutableSparseNDimArray.zeros(nnodes_ie,nterms_s,nnodes_r,ndirs)
#i = 1
#for inode in range(0,nnodes_ie):
#    for inode_diff in range(0,nnodes_r):
#        for idir in range(0,ndirs):
#            for inode in range(0,nnodes_ie):
#                for iterm in range(0,nterms_s):
#                    dgrad1_dx_[inode,iterm,inode_diff,idir] = grad1_[inode,iterm].diff(coords_[inode_diff,idir])
#                    dgrad2_dx_[inode,iterm,inode_diff,idir] = grad2_[inode,iterm].diff(coords_[inode_diff,idir])
#                    dgrad3_dx_[inode,iterm,inode_diff,idir] = grad3_[inode,iterm].diff(coords_[inode_diff,idir])
#                    update_progress("Computing dgrad1_dx, dgrad2_dx, ..  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_ie*nterms_s))
#                    i += 1






#WRITE_____________________

##
## Metrics
##
#f = open("metrics.txt","w")
#i = 1
#for inode in range (0,nnodes_ie):
#    f.write("Metric interpolation node %d \n" % (inode+1))
#    array = numpy.zeros([3, 3])
#    for irow in range(0,3):
#        for icol in range(0,3):
#            data_sym = lambdify(coords_,metrics_[irow,icol,inode],"numpy")
#            data_val = data_sym(*flatten(coords))
#            array[irow,icol] = data_val
#            update_progress("Writing metrics to file             ", i/(nnodes_ie*9))
#            i += 1
#    numpy.savetxt(f,array)
#f.close()
#
##
## jinv
##
#f = open("jinv.txt","w")
#array = numpy.zeros([1])
#i = 1
#for inode in range (0,nnodes_ie):
#    f.write("Jinv interpolation node %d \n" % (inode+1))
#    data_sym = lambdify(coords_,jinv_[inode],"numpy")
#    data_val = data_sym(*flatten(coords))
#    array[0] = data_val
#    numpy.savetxt(f,array)
#    update_progress("Writing jinv to file                ", i/(nnodes_ie))
#    i += 1
#f.close()

##
## Grad1
##
#f = open("grad1.txt","w")
#f.write("Grad1 \n")
#array = numpy.zeros([nnodes_ie,nterms_s])
#i = 1
#for inode in range (0,nnodes_ie):
#    for iterm in range(0,nterms_s):
#        data_sym  = lambdify(coords_,grad1_[inode,iterm],"numpy")
#        data_val = (data_sym(*flatten(coords)))
#        array[inode,iterm] = data_val
#        update_progress("Writing grad1 to file              ", i/(nnodes_ie*nterms_s))
#        i += 1
#numpy.savetxt(f,array)
#f.close()
#
##
## Grad2
##
#f = open("grad2.txt","w")
#f.write("Grad2\n")
#array = numpy.zeros([nnodes_ie,nterms_s])
#i = 1
#for inode in range (0,nnodes_ie):
#    for iterm in range(0,nterms_s):
#        data_sym  = lambdify(coords_,grad2_[inode,iterm],"numpy")
#        data_val = (data_sym(*flatten(coords)))
#        array[inode,iterm] = data_val
#        update_progress("Writing grad2 to file              ", i/(nnodes_ie*nterms_s))
#        i += 1
#numpy.savetxt(f,array)
#f.close()
#
##
## Grad3
##
#f = open("grad3.txt","w")
#f.write("Grad3\n")
#array = numpy.zeros([nnodes_ie,nterms_s])
#i = 1
#for inode in range (0,nnodes_ie):
#    for iterm in range(0,nterms_s):
#        data_sym  = lambdify(coords_,grad3_[inode,iterm],"numpy")
#        data_val = (data_sym(*flatten(coords)))
#        array[inode,iterm] = data_val
#        update_progress("Writing grad3 to file              ", i/(nnodes_ie*nterms_s))
#        i += 1
#numpy.savetxt(f,array)
#f.close()

##
## dmetric_dx
##
#f = open("dmetric_dx.txt","w")
#i = 1
#for inode in range (0,nnodes_ie):
#    for inode_diff in range(0,nnodes_r):
#        for idir in range(0,ndirs):
#            array = numpy.zeros([3,3])
#            f.write("dmetric_dx interpolation node %s, diff_node %s, diff_dir %s \n" % (inode+1,inode_diff+1,idir+1))
#            for irow in range(0,3):
#                for icol in range(0,3):
#                    data_sym = lambdify(coords_,dmetric_dx_[irow,icol,inode,inode_diff,idir],"numpy")
#                    data_val = data_sym(*flatten(coords))
#                    array[irow,icol] = data_val
#                    update_progress("Writing dmetric_dx to file         ", i/(nnodes_ie*nnodes_r*ndirs*3*3))
#                    i += 1
#            numpy.savetxt(f,array)
#f.close()

#
# interp_coords_dx
#
f = open("interp_coords_dx.txt","w")
i = 1
for inode in range (0,nnodes_ie):
    for direct in range (0,3):
        array = numpy.zeros([nnodes_r,ndirs])
        f.write("coord interpolation node %s, coord %s,  row=inode_diff, col=dir \n" % (inode+1,direct+1))
        for inode_diff in range(0,nnodes_r):
            for idir in range(0,ndirs):
                data_sym = lambdify(coords_,interp_coords_dx_[inode,direct,inode_diff,idir],"numpy")
                data_val = data_sym(*flatten(coords))
                array[inode_diff,idir] = data_val
                update_progress("Writing interp_coords_dx to file   ", i/(nnodes_ie*nnodes_r*ndirs*3))
                i += 1
        numpy.savetxt(f,array)
f.close()

##
## djinv_dx
##
#f = open("djinv_dx.txt","w")
#i = 1
#for inode in range (0,nnodes_ie):
#    array = numpy.zeros([nnodes_r,ndirs])
#    f.write("djinv_dx interpolation node %s, row=inode_diff, col=dir \n" % (inode+1))
#    for inode_diff in range(0,nnodes_r):
#        for idir in range(0,ndirs):
#            data_sym = lambdify(coords_,djinv_dx_[inode,inode_diff,idir],"numpy")
#            data_val = data_sym(*flatten(coords))
#            array[inode_diff,idir] = data_val
#            update_progress("Writing djinv_dx to file           ", i/(nnodes_ie*nnodes_r*ndirs))
#            i += 1
#    numpy.savetxt(f,array)
#f.close()

#
# dinvmass_dx
#
f = open("dinvmass_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dinvmass_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nterms_s,nterms_s])
        for irow in range(0,nterms_s):
            for icol in range(0,nterms_s):
                data_sym  = lambdify(coords_,dinvmass_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dinvmass_dx to file        ", i/(nterms_s*nnodes_r*ndirs*nterms_s))
                i += 1
        numpy.savetxt(f,array)
f.close()

#
# dbr2_vol_dx
#
f = open("dbr2_vol_face1_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face1_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_vol_face1_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_vol_face1_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_vol_face2_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face2_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_vol_face2_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_vol_face2_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_vol_face3_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face3_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_vol_face3_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_vol_face3_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_vol_face4_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face4_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_vol_face4_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_vol_face4_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_vol_face5_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face5_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_vol_face5_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_vol_face5_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_vol_face6_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_vol_face6_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_ie,nnodes_if])
        for irow in range(0,nnodes_ie):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_vol_face6_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_vol_face6_dx to file  ", i/(nnodes_ie*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()

#
# dbr2_face_dx
#
f = open("dbr2_face_face1_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face1_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_face_face1_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_face_face1_dx to file  ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_face_face2_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face2_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_face_face2_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_face_face2_dx to file  ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_face_face3_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face3_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_face_face3_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_face_face3_dx to file  ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_face_face4_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face4_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_face_face4_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_face_face4_dx to file  ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_face_face5_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face5_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_face_face5_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_face_face5_dx to file  ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()
f = open("dbr2_face_face6_dx.txt","w")
i = 1
for inode_diff in range(0,nnodes_r):
    for idir in range(0,ndirs):
        f.write("dbr2_face_face6_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
        array = numpy.zeros([nnodes_if,nnodes_if])
        for irow in range(0,nnodes_if):
            for icol in range(0,nnodes_if):
                data_sym  = lambdify(coords_,dbr2_face_face6_dx_[irow,icol,inode_diff,idir],"numpy")
                data_val = (data_sym(*flatten(coords)))
                array[irow,icol] = data_val
                update_progress("Writing dbr2_face_face6_dx to file  ", i/(nnodes_if*nnodes_r*ndirs*nnodes_if))
                i += 1
        numpy.savetxt(f,array)
f.close()

##
## dgrad1_dx
##
#f = open("dgrad1_dx.txt","w")
#i = 1
#for inode_diff in range(0,nnodes_r):
#    for idir in range(0,ndirs):
#        f.write("dgrad1_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
#        array = numpy.zeros([nnodes_ie,nterms_s])
#        for irow in range(0,nnodes_ie):
#            for icol in range(0,nterms_s):
#                data_sym  = lambdify(coords_,dgrad1_dx_[irow,icol,inode_diff,idir],"numpy")
#                data_val = (data_sym(*flatten(coords)))
#                array[irow,icol] = data_val
#                update_progress("Writing dgrad1_dx to file          ", i/(nnodes_ie*nnodes_r*ndirs*nterms_s))
#                i += 1
#        numpy.savetxt(f,array)
#f.close()
#
##
## dgrad2_dx
##
#f = open("dgrad2_dx.txt","w")
#i = 1
#for inode_diff in range(0,nnodes_r):
#    for idir in range(0,ndirs):
#        f.write("dgrad2_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
#        array = numpy.zeros([nnodes_ie,nterms_s])
#        for irow in range(0,nnodes_ie):
#            for icol in range(0,nterms_s):
#                data_sym  = lambdify(coords_,dgrad2_dx_[irow,icol,inode_diff,idir],"numpy")
#                data_val = (data_sym(*flatten(coords)))
#                array[irow,icol] = data_val
#                update_progress("Writing dgrad2_dx to file          ", i/(nnodes_ie*nnodes_r*ndirs*nterms_s))
#                i += 1
#        numpy.savetxt(f,array)
#f.close()
#
##
## dgrad3_dx
##
#f = open("dgrad3_dx.txt","w")
#i = 1
#for inode_diff in range(0,nnodes_r):
#    for idir in range(0,ndirs):
#        f.write("dgrad3_dx => diff_node %s, diff_dir %s \n" % (inode_diff+1,idir+1))
#        array = numpy.zeros([nnodes_ie,nterms_s])
#        for irow in range(0,nnodes_ie):
#            for icol in range(0,nterms_s):
#                data_sym  = lambdify(coords_,dgrad3_dx_[irow,icol,inode_diff,idir],"numpy")
#                data_val = (data_sym(*flatten(coords)))
#                array[irow,icol] = data_val
#                update_progress("Writing dgrad3_dx to file          ", i/(nnodes_ie*nnodes_r*ndirs*nterms_s))
#                i += 1
#        numpy.savetxt(f,array)
#f.close()








