# ============================================================================
# 3D DRONE DYNAMICS SIMULATION
# ============================================================================
# Dette program beregner kræfter, momenter og accelerationer for en drone
# i 3D rummet ved hjælp af Newton-Euler ligninger.
# ============================================================================

import numpy as np
import math as m

# ============================================================================
# SYSTEM PARAMETRE
# ============================================================================

# Euler vinkler i radianer (Z-Y-X konvention)
# roll (x): rotation omkring x-aksen
# pitch (y): rotation omkring y-aksen  
# yaw (z): rotation omkring z-aksen
x, y, z = np.deg2rad([30, 5, 3])  # roll, pitch, yaw

# Rotorkræfter i lokal koordinatsystem [Fx, Fy, Fz] (N)
# Antaget at rotorerne kun genererer kraft i z-retning
f_1 = np.array([[0],[0],[10]])  # Rotor 1
f_2 = np.array([[0],[0],[25]])  # Rotor 2
f_3 = np.array([[0],[0],[10]])   # Rotor 3
f_4 = np.array([[0],[0],[20]])  # Rotor 4

# Dronens masse (kg)
drone_mass = 4

# Afstand fra massecentrum til hver rotor (m)
L = 0.4

# Payload kraft i globalt koordinatsystem [Fx, Fy, Fz] (N)
# Dette kunne være en ekstern kraft fra f.eks. vind eller last
p = np.array([[6], [-3], [-30]])

# Vinkelhastighed i globalt koordinatsystem [wx, wy, wz] (rad/s)
omega = np.array([0.1, 0.5, 0.1])

# ============================================================================
# FUNKTIONER
# ============================================================================

def derive_rotation_matrix(roll, pitch, yaw, XYZ):
    """
    Beregner rotationsmatricen fra lokal til global koordinatsystem.
    
    Rotationsmatricen transformerer vektorer fra dronens lokal koordinatsystem
    til det globale (inertielle) koordinatsystem.
    
    Parameters:
    -----------
    roll : float
        Rotation omkring x-aksen (radianer)
    pitch : float
        Rotation omkring y-aksen (radianer)
    yaw : float
        Rotation omkring z-aksen (radianer)
    XYZ : bool
        True: bruger X-Y-Z konvention (R = Rx @ Ry @ Rz)
        False: bruger Z-Y-X konvention (R = Rz @ Ry @ Rx)
    
    Returns:
    --------
    R : ndarray (3x3)
        Rotationsmatrice
    """
    # Rotation omkring x-aksen (roll)
    Rx = np.array([
        [1, 0, 0],
        [0, m.cos(roll), -m.sin(roll)],
        [0, m.sin(roll),  m.cos(roll)]
    ])

    # Rotation omkring y-aksen (pitch)
    Ry = np.array([
        [m.cos(pitch), 0, m.sin(pitch)],
        [0, 1, 0],
        [-m.sin(pitch), 0, m.cos(pitch)]
    ])

    # Rotation omkring z-aksen (yaw)
    Rz = np.array([
        [m.cos(yaw), -m.sin(yaw), 0],
        [m.sin(yaw),  m.cos(yaw), 0],
        [0, 0, 1]
    ])
    
    # Sammensæt rotationsmatrice baseret på konvention
    if XYZ == True:
        R = Rx @ Ry @ Rz  # X-Y-Z konvention
    if XYZ == False:
        R = Rz @ Ry @ Rx  # Z-Y-X konvention (ofte brugt i luftfart)

    print("================")
    print("Rotation matrix:")
    print(R)
    print("================")
    return R

def calculate_global_inertia(R):
    """
    Transformerer inertitensoren fra lokal til global koordinatsystem.
    
    Når dronen roterer, ændres inertimomentets repræsentation i det
    globale koordinatsystem. Formlen er: I_global = R @ I_local @ R^T
    
    Parameters:
    -----------
    R : ndarray (3x3)
        Rotationsmatrice
    
    Returns:
    --------
    I : ndarray (3x3)
        Inertitensor i globalt koordinatsystem
    """
    # Transponeret rotationsmatrice
    R_t = R.T

    # Lokal inertitensor (symmetrisk, diagonal for symmetrisk drone)
    # Ixx, Iyy: inertimoment omkring x og y akser
    # Izz: inertimoment omkring z-aksen (lodret)
    I_local = np.array([[3, 0, 0],
                        [0, 3, 0],
                        [0, 0, 1]])
    
    # Transformer til globalt koordinatsystem
    I = R @ I_local @ R_t
    print("================")
    print("Global inertia tensor")
    print(I)
    print("================")
    return(I)

def total_force(R, f1, f2, f3, f4, drone_mass, payload_f):
    """
    Beregner den totale kraft på dronen i globalt koordinatsystem.
    
    Samler alle kræfter: rotorkræfter (transformeret til globalt),
    tyngdekraft og eventuel payload kraft.
    
    Parameters:
    -----------
    R : ndarray (3x3)
        Rotationsmatrice
    f1, f2, f3, f4 : ndarray (3x1)
        Rotorkræfter i lokal koordinatsystem
    drone_mass : float
        Dronens masse (kg)
    payload_f : ndarray (3x1) eller None
        Payload kraft i globalt koordinatsystem
    
    Returns:
    --------
    f_total : ndarray (3x1)
        Total kraft i globalt koordinatsystem
    """
    # Tyngdekraft i globalt koordinatsystem (peger nedad)
    mg = drone_mass * -9.82  # g = 9.82 m/s²
    gravity = np.array([[0],[0],[mg]])

    # Summer alle rotorkræfter i lokalt koordinatsystem
    f_local = f1 + f2 + f3 + f4
    
    # Transformer rotorkræfter til globalt koordinatsystem
    f_global = R @ f_local

    # Tilføj payload hvis den findes (payload er allerede i globalt koordinatsystem)
    if payload_f is not None:
        f_total = f_global + payload_f + gravity
        print("=====================")
        print("total force")
        print(f_total)
        print("=====================")
        return f_total
    else:
        # Kun rotorkræfter og tyngdekraft
        f_total = f_global + gravity
        print("=====================")
        print("total force")
        print(f_total)
        print("=====================")
        return f_total
    

    

def total_torque(R, f1, f2, f3, f4, L):
    """
    Beregner det totale moment (torque) på dronen.
    
    Momentet beregnes ved krydsproduktet af positionsvektorer og kræfter
    for hver rotor: τ = r × F. Beregnes først i lokalt koordinatsystem,
    derefter transformeret til globalt.
    
    Parameters:
    -----------
    R : ndarray (3x3)
        Rotationsmatrice
    f1, f2, f3, f4 : ndarray (3x1)
        Rotorkræfter i lokal koordinatsystem
    L : float
        Afstand fra massecentrum til rotor (m)
    
    Returns:
    --------
    tau_global : ndarray (3,)
        Totalt moment i globalt koordinatsystem
    """
    # Rotorpositioner i lokalt koordinatsystem (jævnt fordelt i + konfiguration)
    # Rotor 1: foran (+x)
    # Rotor 2: højre (+y)
    # Rotor 3: bag (-x)
    # Rotor 4: venstre (-y)
    r1 = np.array([ L, 0, 0])
    r2 = np.array([ 0, L, 0])
    r3 = np.array([-L, 0, 0])
    r4 = np.array([ 0,-L, 0])

    # Beregn moment for hver rotor: τ = r × F (krydsproduktet)
    # Flatten() konverterer fra (3,1) til (3,) for np.cross
    tau1 = np.cross(r1, f1.flatten())
    tau2 = np.cross(r2, f2.flatten())
    tau3 = np.cross(r3, f3.flatten())
    tau4 = np.cross(r4, f4.flatten())

    # Totalt moment i lokalt koordinatsystem
    tau_local = tau1 + tau2 + tau3 + tau4

    print("=====================")
    print("tau local:", tau_local)
    
    # Transformer til globalt koordinatsystem
    tau_global = R @ tau_local
    
    print("tau global", tau_global)
    print("=====================")
    return tau_global


# ============================================================================
# HOVEDPROGRAM - BEREGNING AF DYNAMIK
# ============================================================================

# 1. Beregn rotationsmatrice (Z-Y-X konvention)
R = derive_rotation_matrix(x, y, z, True)

# 2. Beregn global inertitensor
I_global = calculate_global_inertia(R)

# 3. Beregn total kraft (inkl. payload)
F_total = total_force(R, f_1, f_2, f_3, f_4, drone_mass, None)

# 4. Beregn totalt moment
tau_global = total_torque(R, f_1, f_2, f_3, f_4, L)

# ============================================================================
# LINEÆR ACCELERATION (Newton's 2. lov: F = ma → a = F/m)
# ============================================================================
v_dot = F_total / drone_mass
print("=====================")
print("Linear acceleration:")
print(v_dot)
print("=====================")

# ============================================================================
# VINKELACCELERATION (Euler's rotationsligning)
# ============================================================================
# Euler's ligning: τ = I·ω̇ + ω × (I·ω)
# Omskrevet: ω̇ = I⁻¹ · (τ - ω × (I·ω))

# Beregn I·ω (angular momentum)
Iomega = I_global @ omega

# Beregn gyroskopisk led (ω × (I·ω))
# Dette led skyldes at koordinatsystemet roterer
gyroscopic = np.cross(omega, Iomega)

# Beregn vinkelacceleration
omega_dot = np.linalg.inv(I_global) @ (tau_global - gyroscopic)
print("=====================")
print("Angular acceleration:")
print(omega_dot)
print("=====================")
