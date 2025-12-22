import numpy as np
import math as m

# Euler angles given in degrees: Z-Y-X
psi, theta, phi = np.deg2rad([45, 0, 5])  # yaw, pitch, roll


def derive_rotation_matrix(psi, theta, phi):
    Rx = np.array([
        [1, 0, 0],
        [0, m.cos(phi), -m.sin(phi)],
        [0, m.sin(phi),  m.cos(phi)]
    ])

    Ry = np.array([
        [m.cos(theta), 0, m.sin(theta)],
        [0, 1, 0],
        [-m.sin(theta), 0, m.cos(theta)]
    ])

    Rz = np.array([
        [m.cos(psi), -m.sin(psi), 0],
        [m.sin(psi),  m.cos(psi), 0],
        [0, 0, 1]
    ])

    R = Rz @ Ry @ Rx
    print(R)
    return R

def calculate_global_inertia(R):
    # when the local inertia is given, and the rotation matrix is calculated the formula is I = RI´R_T
    R_t = R.T

    #local inertia:
    I_local = np.array([[4, 0, 0],
                 [0, 4, 0],
                 [0, 0, 2]])
    
    I = R @ I_local @ R_t

    print(I)
    return(I)

R = derive_rotation_matrix(psi, theta, phi)

calculate_global_inertia(R)
