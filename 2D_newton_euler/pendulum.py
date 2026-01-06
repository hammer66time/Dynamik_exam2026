import sympy as sp

# ---------- symbols ----------
t = sp.symbols("t")
l = 0.8
m = 2
g = 9.82
t_eval = 0.25

# ---------- motion ----------
theta = 3*sp.sin(sp.pi*t)
omega = sp.diff(theta, t)
alpha = sp.diff(omega, t)

theta0 = float(theta.subs(t, t_eval))
omega0 = float(omega.subs(t, t_eval))
alpha0 = float(alpha.subs(t, t_eval))

print("\n===== PART (1): Kinematics =====")
print(f"θ(0.25) = {theta0:.3f} rad")
print(f"ω(0.25) = {omega0:.2f} rad/s")
print(f"α(0.25) = {alpha0:.2f} rad/s²")

# ---------- position of C ----------
rc = sp.Matrix([
    (l/2)*sp.cos(theta),
    (l/2)*sp.sin(theta)
])

vc = sp.diff(rc, t)
ac = sp.diff(vc, t)

vc0 = vc.subs(t, t_eval)
ac0 = ac.subs(t, t_eval)

print("\nVelocity of C at t=0.25 s:")
print(f"v_C = [{float(vc0[0]):.3f}, {float(vc0[1]):.3f}] m/s")

print("\nAcceleration of C at t=0.25 s:")
print(f"a_C = [{float(ac0[0]):.2f}, {float(ac0[1]):.2f}] m/s²")

# ---------- inertia ----------
IC = (1/12)*m*l**2
print("\n===== PART (2): Inertia =====")
print(f"I_C = {float(IC):.3f} kg·m²")

# ---------- torque ----------
def cross2d(ax, ay, bx, by):
    return ax*by - ay*bx

# numbers at time t
theta_num = theta0

# positions of C and B
rC = sp.Matrix([
    (l/2)*sp.cos(theta_num),
    (l/2)*sp.sin(theta_num)
])

rB = sp.Matrix([
    l*sp.cos(theta_num),
    l*sp.sin(theta_num)
])

# forces
Fx = 5
Fy = -2

# torques from forces about A
tau_g = cross2d(rC[0], rC[1], 0, -m*g)
tau_F = cross2d(rB[0], rB[1], Fx, Fy)

# inertia about A
IA = IC + m*(l/2)**2

tau_drive = IA*alpha0 - tau_g - tau_F

print("\n===== PART (5): Driving torque at A =====")
print(f"τ = {float(tau_drive):.2f} N·m")
