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
    #print(R)
    return R

def calculate_global_inertia(R):
    # when the local inertia is given, and the rotation matrix is calculated the formula is I = RI´R_T
    R_t = R.T

    #local inertia:
    I_local = np.array([[4, 0, 0],
                 [0, 4, 0],
                 [0, 0, 2]])
    
    I = R @ I_local @ R_t

    #print(I)
    return(I)

R = derive_rotation_matrix(psi, theta, phi)

#calculate_global_inertia(R)

drone_mass = 6  # kg

#motor forces:
f1 = np.array([[0],[0],[10]]) #N
f2 = np.array([[0],[0],[25]]) #N
f3 = np.array([[0],[0],[5]]) #N
f4 = np.array([[0],[0],[30]]) #N

#payload force:
payload = np.array([[6], [-3], [-30]]) #N

#gravitational force:

g_force = np.array([[0], [0], [-9.81 * drone_mass]])

force_local = f1 + f2 + f3 + f4 


force_world = R @ force_local

force_world = force_world + g_force + payload 


acceleration = force_world / drone_mass

#print(acceleration)

def newton_euler_total_force(sum_force, mass, acc_c):
    pass


def omega_dot(omega, forces, positions, inertia):
    """
    Beregner vinkelacceleration omega_dot vha. Newton–Euler.

    Parameters
    ----------
    omega : np.array (3,)
        Angular velocity [rad/s] i body-frame
    forces : list of np.array (3,)
        Rotor-kræfter (body-frame)
    positions : list of np.array (3,)
        Rotor-positioner ift. COM (body-frame)
    inertia : np.array (3,3)
        Lokal inertimatrix (body-frame)

    Returns
    -------
    omega_dot : np.array (3,)
        Angular acceleration [rad/s^2]
    """

    # 1) Beregn samlet moment τ = Σ (r × F)
    tau = np.zeros(3)
    for r, f in zip(positions, forces):
        tau += np.cross(r, f)

    # 2) Gyroskopisk led: ω × (Iω)
    gyro = np.cross(omega, inertia @ omega)

    # 3) Newton–Euler rotation
    omega_dot = np.linalg.inv(inertia) @ (tau - gyro)

    return omega_dot

# Inerti (body-frame)
I = np.array([
    [4, 0, 0],
    [0, 4, 0],
    [0, 0, 2]
])

# Vinkelhastighed
omega = np.array([0.1, 0.1, 0.7])

# Arm-længde
L = 0.1

# Rotor-positioner (body-frame)
positions = [
    np.array([ L,  0, 0]),
    np.array([ 0,  L, 0]),
    np.array([-L,  0, 0]),
    np.array([ 0, -L, 0])
]

# Rotor-kræfter (body-frame)
forces = [
    np.array([0, 0, 10]),
    np.array([0, 0, 25]),
    np.array([0, 0, 5]),
    np.array([0, 0, 30])
]

# Beregn omega_dot
omega_dot_result = omega_dot(omega, forces, positions, I)
print(omega_dot_result)


