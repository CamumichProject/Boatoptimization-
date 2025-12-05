from ship import Ship as sp
from physics import to_meters
import numpy as np
import cvxpy as cp

#Our inputs are (self, length: float, draught: float, beam: float, speed: float,
#                 CB: float, midship_coefficient: float, 
#                 waterplane_coefficient: float, lcb: float, prop_diameter: float)

#We define our optimizer as optimizing over basic inputs BLT 

L = to_meters(406.7)
D = cp.Variable(pos=True)
B = cp.Variable(pos=True)
CB = cp.Variable(pos=True)
k = 20 #Using a constant speed for convenient
Prismatic = 0.75 #We use a constant prismatic coefficient
Displacement = 10000 #We assume we must have a constant ship size
WPCF= 0.8 #Again assume this is constant (will grab from series 60)
MSCF = 0.8 #Constant again
LCB = 0.5 * L 
PropD = 6 

x = sp(L, D, B, k, CB, MSCF, WPCF, LCB, PropD)

costWf = 10 * x.calc_wake_fraction()
Costthrr = 10 * x.calc_thrust_reduction()
Costfrc= 10 * fictional_resistance_coeff(x.speed,x.length)
Costr = 10000 * x.resistance

objective = cp.Minimize(cost)
constraints = 200-L
constraints = constraints.append(50-B)
constraints = constraints.append(50-D)
prob = cp.problem(objective, constraints)
prob.solve()
