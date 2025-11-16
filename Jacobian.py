import sympy as sp
from sympy import latex
import math as m

theta1, theta2 , s = sp.symbols("theta1 theta2 s")

y1 = 5 * sp.cos(theta1) + 3 * sp.cos(theta1 + theta2)
y2 = 5 * sp.sin(theta1) + 3 * sp.sin(theta1 + theta2)
y3 = 0.5 * s + 3

#--------------------------------------------------------------------
#determine the jacobian

dy1_dtheta1  =   sp.diff(y1, theta1)
dy2_dtheta1 =   sp.diff(y2, theta1)
dy3_dtheta1 =   sp.diff(y3, theta1)
dy1_dtheta2 =   sp.diff(y1, theta2)
dy2_dtheta2 =   sp.diff(y2, theta2)
dy3_dtheta2 =   sp.diff(y3, theta2)
dy1_ds      =   sp.diff(y1, s)
dy2_ds      =   sp.diff(y2, s)
dy3_ds      =   sp.diff(y3, s)

jacobian = sp.Matrix([  [dy1_dtheta1, dy1_dtheta2, dy1_ds],
                        [dy2_dtheta1, dy2_dtheta2, dy2_ds],
                        [dy3_dtheta1, dy3_dtheta2, dy3_ds]])

print("------------------------------------------------------------------")
print("Jacobian")
sp.pprint(jacobian)
print("------------------------------------------------------------------")
print("in latex form:")
print(latex(jacobian))
print("------------------------------------------------------------------")

#-----------------------------------------------------------------------
# calculate singularity

detJ = jacobian.det()
print("------------------------------------------------------------------")
print("det(J):")
sp.pprint(detJ)

# simplify determinant
detJ_simplified = sp.simplify(detJ)
print("------------------------------------------------------------------")
print("Simplified det(J):")
sp.pprint(detJ_simplified)

# solve det(J)=0
print("------------------------------------------------------------------")
print("Solve det(J)=0:")
solutions = sp.solve(sp.Eq(detJ_simplified, 0), theta2)
sp.pprint(solutions)

#----------------------------------------------------------------------
# Illustrate graphicaly singular config:


