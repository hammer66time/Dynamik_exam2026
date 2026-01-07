import sympy as sp

# -------------------------- symbols & constants --------------------------
t = sp.symbols("t")

l1 = 0.5
l2 = 0.3

m1 = 3
m2 = 2

I1 = 0.5     # about C1
I2 = 0.3     # about C2

g = 9.82
t_eval = 0.1

# -------------------------- joint motion --------------------------
theta1 = 3*sp.sin(sp.pi*t + 1.1)
theta2_rel = 2*sp.sin(sp.pi*t + 0.4)

# absolute angle for link 2
theta2 = theta1 + theta2_rel

omega1 = sp.diff(theta1, t)
omega2 = sp.diff(theta2, t)

alpha1 = sp.diff(omega1, t)
alpha2 = sp.diff(omega2, t)

# evaluate at time
th1_0 = float(theta1.subs(t, t_eval))
th2_0 = float(theta2.subs(t, t_eval))

om1_0 = float(omega1.subs(t, t_eval))
om2_0 = float(omega2.subs(t, t_eval))

al1_0 = float(alpha1.subs(t, t_eval))
al2_0 = float(alpha2.subs(t, t_eval))

print("\n===== Joint kinematics at t = 0.1 s =====")
print(f"theta1 = {th1_0:.4f} rad,  omega1 = {om1_0:.4f} rad/s,  alpha1 = {al1_0:.4f} rad/s²")
print(f"theta2 = {th2_0:.4f} rad,  omega2 = {om2_0:.4f} rad/s,  alpha2 = {al2_0:.4f} rad/s²")

# -------------------------- positions of mass centers --------------------------
# IMPORTANT: angles measured from y-axis
# x = r*sin(theta), y = r*cos(theta)

# C1 position
rc1 = sp.Matrix([
    (l1/2)*sp.sin(theta1),
    (l1/2)*sp.cos(theta1)
])

# C2 position
rc2 = sp.Matrix([
    l1*sp.sin(theta1) + (l2/2)*sp.sin(theta2),
    l1*sp.cos(theta1) + (l2/2)*sp.cos(theta2)
])

# velocities and accelerations
vc1 = sp.diff(rc1, t)
ac1 = sp.diff(vc1, t)

vc2 = sp.diff(rc2, t)
ac2 = sp.diff(vc2, t)

rc1_0 = rc1.subs(t, t_eval)
rc2_0 = rc2.subs(t, t_eval)

ac1_0 = ac1.subs(t, t_eval)
ac2_0 = ac2.subs(t, t_eval)

print("\n===== Accelerations of mass centers =====")
print(f"a_C1 = [{float(ac1_0[0]):.4f}, {float(ac1_0[1]):.4f}] m/s²")
print(f"a_C2 = [{float(ac2_0[0]):.4f}, {float(ac2_0[1]):.4f}] m/s²")

# -------------------------- helpers --------------------------
def cross2d(a, b):
    return a[0]*b[1] - a[1]*b[0]

g_vec = sp.Matrix([0, -g])

# -------------------------- Torques --------------------------

# --- PARALLEL AXIS TERMS ---

# link 1 about O
I1_O = I1 + m1*(l1/2)**2

# vectors needed
rOC2 = rc2
rOC2_0 = rOC2.subs(t, t_eval)

# vector from A to C2
rAC2 = sp.Matrix([
    (l2/2)*sp.sin(theta2),
    (l2/2)*sp.cos(theta2)
])
rAC2_0 = rAC2.subs(t, t_eval)

# inertias moved to correct rotation points
I2_A = I2 + m2*(l2/2)**2
I2_O = I2 + m2*(rOC2_0.norm())**2

# ---- torque τ2 about joint A ----
tau2 = I2_A*al2_0 \
       + cross2d(rAC2_0, m2*ac2_0 + m2*g_vec)

# ---- torque τ1 about base O ----
tau1 = I1_O*al1_0 \
       + I2_O*al2_0 \
       + cross2d(rc1_0, m1*ac1_0 + m1*g_vec) \
       + cross2d(rc2_0, m2*ac2_0 + m2*g_vec)

print("\n===== Joint torques =====")
print(f"tau1 = {float(tau1):.4f} N·m")
print(f"tau2 = {float(tau2):.4f} N·m")

# -------------------------- Reaction forces --------------------------
Ox, Oy, Ax, Ay = sp.symbols("Ox Oy Ax Ay")

# link 2 translation
eq2x = Ax - m2*ac2_0[0]
eq2y = Ay - m2*ac2_0[1] - m2*g

# link 1 translation
eq1x = Ox - Ax - m1*ac1_0[0]
eq1y = Oy - Ay - m1*ac1_0[1] - m1*g

solutions = sp.solve((eq1x, eq1y, eq2x, eq2y), (Ox, Oy, Ax, Ay))

print("\n===== Reaction forces =====")
print(f"O_x = {float(solutions[Ox]):.4f} N")
print(f"O_y = {float(solutions[Oy]):.4f} N")
print(f"A_x = {float(solutions[Ax]):.4f} N")
print(f"A_y = {float(solutions[Ay]):.4f} N")
