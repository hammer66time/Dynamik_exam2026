import sympy as sp

# Define time variable:
t = sp.symbols("t")

# Define function:
theta = 3 * sp.sin(3.14 * t)

# Time derivative, which is the velocity
omega = sp.diff(theta, t)

alpha = sp.diff(omega, t)

#print(f"Velocity: {omega}")
#print(f"acceleration: {alpha}")

# Velocity result
def velocity(t):
    result = 9.42 * sp.cos(3.14 * t)
    return result

def acceleration(t):
    result = -29.58 * sp.sin(3.14 * t)
    return result

acc = acceleration(0.25)

vel = velocity(0.25)
#print(vel)
#print(acc)

#cal = -0.4 * sp.sin(2.12) * (6.66 ** 2) + 0.4 * sp.cos(2.12) * -20.

#cal = (2*(0.8**2))/12

cal = 0.1 + 2 * (0.4**2)
#print(cal)

def vector_cross_2d(Ax, Ay, Bx, By):
    cross = Ax * By - Ay * Bx
    return cross

theta = 2.12
m = 2
g = 9.82

tau_g = vector_cross_2d(
    0.4*sp.cos(theta),
    0.4*sp.sin(theta),
    0,
    -m*g
)

print("tau_g =", float(tau_g))

tau_F = vector_cross_2d(
    0.8*sp.cos(theta),
    0.8*sp.sin(theta),
    5,
    -2
)

print("tau_F =", float(tau_F))

I_A = 0.42
alpha = -20.9

tau = I_A*alpha - tau_g - tau_F
print("tau =", float(tau))
