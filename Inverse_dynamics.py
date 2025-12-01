import sympy as sp
import math as m

def compute_for_joints():
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

    time = sp.symbols('time', real=True)
    l1, l2 = 0.5, 0.3

    theta1 = 3 * sp.sin(3.141*time + 1.1)
    theta2 = 2 * sp.sin(3.141*time + 0.4)

    # positioner
    r_c1 = sp.Matrix([-(l1/2)*sp.sin(theta1),
                    ( l1/2)*sp.cos(theta1)])
    r_c2 = sp.Matrix([-l1*sp.sin(theta1) - (l2/2)*sp.sin(theta1+theta2),
                    l1*sp.cos(theta1) + (l2/2)*sp.cos(theta1+theta2)])

    # hastigheder = tidsafledte
    v_c1 = sp.diff(r_c1, time)
    v_c2 = sp.diff(r_c2, time)

    # acceleration:
    a_c1 = sp.diff(v_c1, time)
    a_c2 = sp.diff(v_c2, time)

    # hvis du vil have numeriske værdier ved t=0.1:
    v_c1_num = sp.N(v_c1.subs(time, 0.1))
    v_c2_num = sp.N(v_c2.subs(time, 0.1))

    a_c1_num = sp.N(a_c1.subs(time, 0.1))
    a_c2_num = sp.N(a_c2.subs(time, 0.1))
    
    print("------------------------------------------------")
    print("v_c1 =", v_c1_num)
    print("v_c2 =", v_c2_num)
    print("a_c1 =", a_c1_num)
    print("a_c2 =", a_c2_num)
    print("------------------------------------------------")

compute_velotity_for_c()

