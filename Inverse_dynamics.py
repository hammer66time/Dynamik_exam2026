import sympy as sp
import math as m

t = sp.symbols("t")

theta1 = 3 * sp.sin(3.14 * t + 1.1)
theta2 = 2 * sp.sin(3.14 * t + 0.4)

# ω1 = dθ1/dt
omega1 = sp.diff(theta1, t)

# ω2 = d(θ1 + θ2)/dt = ω1 + ω2_individual
omega2 = sp.diff(theta1 + theta2, t)

print("vel link 1:", omega1)
print("vel link 2:", omega2)

acc_link1 = sp.diff(omega1, t)

acc_link2 = sp.diff(omega2, t)

print("acc link 1:", acc_link1)
print("acc link 2:", acc_link2)

#------------------------------------------------------
#insert number for t:

time = 0.1

om1 = 9.42*m.cos(3.14*time + 1.1)
om2 = 6.28*m.cos(3.14*time + 0.4) + 9.42*m.cos(3.14*time + 1.1)

acc1 = -29.5788*m.sin(3.14*time + 1.1)
acc2 = -19.7192*m.sin(3.14*time + 0.4) - 29.5788*m.sin(3.14*time + 1.1)

print("------------------------------------------------")
print("with t = 0.1")
print("velocity link 1:", om1)
print("velocity link 2:", om2)
print("acceleration link 1:", acc2)
print("acceleration link 2:", acc2)
print("------------------------------------------------")

#------------------------------------------------------
#Velocity at link midpoints:
def compute_velotity_for_c():

    l1 = 0.5
    l2 = 0.3

    theta1_t = 3 * sp.sin(3.14 * time + 1.1)
    theta2_t = 2 * sp.sin(3.14 * time + 0.4)

    c1_vector = sp.Matrix([[-(l1/2) * m.sin(theta1_t)],
                            [(l1/2) * m.cos(theta1_t)]])

    velocity_c2 = sp.Matrix([[-om1 * l1 * m.sin(theta1_t) - om2 * (l2/2) * m.sin(theta1_t + theta2_t)],
                            [om1 * l1 * m.cos(theta1_t) + om2 * (l2/2) * m.cos(theta1_t + theta2_t)]])


    velocity_c1 = om1 * c1_vector

    print("------------------------------------------------")
    print("Velocity at c1:", velocity_c1)
    print("Velocity at c2:", velocity_c2)
    print("------------------------------------------------")


def compute_accelerations(L1=0.5, L2=0.3, t=0.1):
    pi = m.pi

    # angles
    arg1 = pi * t + 1.1
    arg2 = pi * t + 0.4

    theta1 = 3.0 * m.sin(arg1)
    theta2 = 2.0 * m.sin(arg2)
    theta12 = theta1 + theta2

    # angular velocities and accelerations (exact with pi)
    omega1 = 3.0 * pi * m.cos(arg1)        # d/dt theta1
    omega2_ind = 2.0 * pi * m.cos(arg2)    # d/dt theta2 (relative)
    omega_total = omega1 + omega2_ind        # total for link2

    alpha1 = -3.0 * (pi**2) * m.sin(arg1)   # d2/dt2 theta1
    alpha2_ind = -2.0 * (pi**2) * m.sin(arg2) # d2/dt2 theta2
    alpha_total = alpha1 + alpha2_ind

    # r_C1 (from joint1 to C1)
    rC1_x = (L1/2.0) * m.cos(theta1)
    rC1_y = (L1/2.0) * m.sin(theta1)

    # Tangential and centripetal for C1
    tang_C1_x = -alpha1 * rC1_y
    tang_C1_y =  alpha1 * rC1_x
    cent_C1_x  = - (omega1**2) * rC1_x
    cent_C1_y  = - (omega1**2) * rC1_y

    aC1_x = tang_C1_x + cent_C1_x
    aC1_y = tang_C1_y + cent_C1_y

    # r_C2: position of C2 relative to base
    rC2_x = L1 * m.cos(theta1) + (L2/2.0) * m.cos(theta12)
    rC2_y = L1 * m.sin(theta1) + (L2/2.0) * m.sin(theta12)

    # For C2 we sum contributions:
    #  - contribution due to rotation of link1 acting on the vector to the joint at link2:
    r_joint2_x = L1 * m.cos(theta1)
    r_joint2_y = L1 * m.sin(theta1)
    # tangential and centripetal for that first piece (rotation omega1, alpha1)
    tang_j1_x = -alpha1 * r_joint2_y
    tang_j1_y =  alpha1 * r_joint2_x
    cent_j1_x  = - (omega1**2) * r_joint2_x
    cent_j1_y  = - (omega1**2) * r_joint2_y

    #  - contribution due to rotation of link2 about joint2 (use omega_total, alpha_total on r from joint2 to C2)
    r_link2_half_x = (L2/2.0) * m.cos(theta12)
    r_link2_half_y = (L2/2.0) * m.sin(theta12)

    tang_l2_x = -alpha_total * r_link2_half_y
    tang_l2_y =  alpha_total * r_link2_half_x
    cent_l2_x  = - (omega_total**2) * r_link2_half_x
    cent_l2_y  = - (omega_total**2) * r_link2_half_y

    # total acceleration at C2 = (contrib from link1 endpoint) + (contrib from link2 half)
    aC2_x = (tang_j1_x + cent_j1_x) + (tang_l2_x + cent_l2_x)
    aC2_y = (tang_j1_y + cent_j1_y) + (tang_l2_y + cent_l2_y)

    # Return a dictionary with intermediate values if you want to inspect
    return {
        "theta1": theta1, "theta2": theta2, "theta12": theta12,
        "omega1": omega1, "omega2_ind": omega2_ind, "omega_total": omega_total,
        "alpha1": alpha1, "alpha2_ind": alpha2_ind, "alpha_total": alpha_total,
        "rC1": (rC1_x, rC1_y),
        "aC1": (aC1_x, aC1_y),
        "rC2": (rC2_x, rC2_y),
        "aC2": (aC2_x, aC2_y),
    }

if __name__ == "__main__":
    out = compute_accelerations(L1=0.5, L2=0.3, t=0.1)
    print("theta1, theta2:", out["theta1"], out["theta2"])
    print("omega1, omega2_ind, omega_total:", out["omega1"], out["omega2_ind"], out["omega_total"])
    print("alpha1, alpha2_ind, alpha_total:", out["alpha1"], out["alpha2_ind"], out["alpha_total"])
    print("rC1:", out["rC1"])
    print("aC1 (x,y):", out["aC1"])
    print("rC2:", out["rC2"])
    print("aC2 (x,y):", out["aC2"])
