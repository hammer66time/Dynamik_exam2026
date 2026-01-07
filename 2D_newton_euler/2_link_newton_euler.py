import sympy as sp

#-------------------------- motion at links-----------------------------
# ---------- symbols ----------
t = sp.symbols("t")
l1 = 0.5
l2 = 0.3
m1 = 3
m2 = 2
g = 9.82
t_eval = 0.1

# ---------- motion ----------
#----------- link1 -----------
theta1 = 3*sp.sin(3.141*t+1.1)
omega1 = sp.diff(theta1, t)
alpha1 = sp.diff(omega1, t)

theta0_1 = float(theta1.subs(t, t_eval))
omega0_1 = float(omega1.subs(t, t_eval))
alpha0_1 = float(alpha1.subs(t, t_eval))

print("\n===== Kinematics =====")
print("theta1")
print(f"θ(0.1) = {theta0_1:.3f} rad")
print(f"ω(0.1) = {omega0_1:.2f} rad/s")
print(f"alpha(0.1) = {alpha0_1:.2f} rad/s²")

print("\n===== Symbolske ligninger Link 1 =====")
print("ω₁(t) = dθ₁/dt:")
sp.pprint(omega1)
print("\nα₁(t) = dω₁/dt:")
sp.pprint(alpha1)

#----------- link2 -----------
theta2_rel = 2*sp.sin(3.141*t+0.4)  # relative angle
theta2 = theta1 + theta2_rel  # absolute angle (link2 depends on link1)
omega2 = sp.diff(theta2, t)
alpha2 = sp.diff(omega2, t)

theta0_2 = float(theta2.subs(t, t_eval))
omega0_2 = float(omega2.subs(t, t_eval))
alpha0_2 = float(alpha2.subs(t, t_eval))

print("\n=================================")
print("theta2")
print(f"θ(0.1) = {theta0_2:.3f} rad")
print(f"ω(0.1) = {omega0_2:.2f} rad/s")
print(f"alpha(0.1) = {alpha0_2:.2f} rad/s²")

print("\n===== Symbolske ligninger Link 2 =====")
print("ω₂(t) = dθ₂/dt:")
sp.pprint(omega2)
print("\nα₂(t) = dω₂/dt:")
sp.pprint(alpha2)
print("\n=================================")

#-------------------------- Center of mass positions and accelerations ----------
# Link 1 - center of mass C1
rc1 = sp.Matrix([
    (l1/2)*sp.cos(theta1),
    (l1/2)*sp.sin(theta1)
])
vc1 = sp.diff(rc1, t)
ac1 = sp.diff(vc1, t)

vc1_0 = vc1.subs(t, t_eval)
ac1_0 = ac1.subs(t, t_eval)

print("\n===== Velocities =====")
print("Link 1 - Center of mass velocity:")
print(f"v_C1 = [{float(vc1_0[0]):.3f}, {float(vc1_0[1]):.3f}] m/s")

# Link 2 - center of mass C2
# Position: first go to end of link1, then halfway through link2
rc2 = sp.Matrix([
    l1*sp.cos(theta1) + (l2/2)*sp.cos(theta2),
    l1*sp.sin(theta1) + (l2/2)*sp.sin(theta2)
])
vc2 = sp.diff(rc2, t)
ac2 = sp.diff(vc2, t)

vc2_0 = vc2.subs(t, t_eval)
ac2_0 = ac2.subs(t, t_eval)

print("Link 2 - Center of mass velocity:")
print(f"v_C2 = [{float(vc2_0[0]):.3f}, {float(vc2_0[1]):.3f}] m/s")

print("\n===== Accelerations =====")
print("Link 1 - Center of mass acceleration:")
print(f"a_C1 = [{float(ac1_0[0]):.2f}, {float(ac1_0[1]):.2f}] m/s²")

print("Link 2 - Center of mass acceleration:")
print(f"a_C2 = [{float(ac2_0[0]):.2f}, {float(ac2_0[1]):.2f}] m/s²")

#-------------------------- Inertias ----------
# Using given inertias (about center of mass)
I1 = 0.5      # given directly in assignment
I2 = 0.3      # given directly in assignment

print("\n===== Inertias =====")
print(f"I_1 = {I1:.4f} kg·m²")
print(f"I_2 = {I2:.4f} kg·m²")

#-------------------------- Torque calculations (Moment method) ----------
# 2D cross product z-component helper
def cross2d(a, b):
    """Calculate z-component of cross product for 2D vectors"""
    return a[0]*b[1] - a[1]*b[0]

# gravity vector
g_vec = sp.Matrix([0, -g])

# Position vectors at evaluation time
rc1_0 = rc1.subs(t, t_eval)
rc2_0 = rc2.subs(t, t_eval)

# Evaluate angular accelerations
alpha1_0 = float(alpha1.subs(t, t_eval))
alpha2_0 = float(alpha2.subs(t, t_eval))

#-------------------------- Torque τ2 (moment about joint A) ----------
# Joint A is at the end of link 1
# Vector from A to C2 (center of mass of link 2)
rAC2 = sp.Matrix([
    (l2/2)*sp.cos(theta2),
    (l2/2)*sp.sin(theta2)
])
rAC2_0 = rAC2.subs(t, t_eval)

# Moment equation about A:
# τ2 = I2*α2 + r_AC2 × (m2*aC2 + m2*g)
tau2 = I2*alpha2_0 \
       + cross2d(rAC2_0, m2*ac2_0) \
       + cross2d(rAC2_0, m2*g_vec)

#-------------------------- Torque τ1 (moment about base O) ----------
# Moment equation about O:
# τ1 = I1*α1 + r_C1 × (m1*aC1 + m1*g) + r_C2 × (m2*aC2 + m2*g)
tau1 = I1*alpha1_0 \
       + cross2d(rc1_0, m1*ac1_0) \
       + cross2d(rc2_0, m2*ac2_0) \
       + cross2d(rc1_0, m1*g_vec) \
       + cross2d(rc2_0, m2*g_vec)

print("\n===== Joint Torques =====")
print(f"τ_1 (base joint) = {float(tau1):.3f} N·m")
print(f"τ_2 (elbow joint) = {float(tau2):.3f} N·m")

#-------------------------- Reaction Forces ----------
# Unknown reaction forces at joints
Ox, Oy, Ax, Ay = sp.symbols('Ox Oy Ax Ay')

# Force balance equations (Newton's 2nd law for translation)
# Link 2: sum of forces = m2*aC2
# Forces on link 2: -Ax (from joint A in x), -Ay (from joint A in y), m2*g (gravity)
eq2x = -Ax - 0 + m2*g_vec[0] - m2*ac2_0[0]
eq2y = -Ay + m2*g_vec[1] - m2*ac2_0[1]

# Link 1: sum of forces = m1*aC1
# Forces on link 1: Ox (from base O), Oy (from base O), Ax (reaction from link 2), Ay (reaction from link 2), m1*g (gravity)
eq1x = Ox + Ax + m1*g_vec[0] - m1*ac1_0[0]
eq1y = Oy + Ay + m1*g_vec[1] - m1*ac1_0[1]

# Solve system of equations for reaction forces
solutions = sp.solve((eq1x, eq1y, eq2x, eq2y), (Ox, Oy, Ax, Ay))

Ox_val = float(solutions[Ox])
Oy_val = float(solutions[Oy])
Ax_val = float(solutions[Ax])
Ay_val = float(solutions[Ay])

print("\n===== Reaction Forces =====")
print(f"O_x = {Ox_val:.3f} N")
print(f"O_y = {Oy_val:.3f} N")
print(f"A_x = {Ax_val:.3f} N")
print(f"A_y = {Ay_val:.3f} N")
