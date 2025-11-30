from ship import Ship
from physics import froude_number, reynolds_number
import numpy as np

# Testing out the ship class

ship1 = Ship(205.00, 10.00, 32.00, 25.0, 37500, 0.98, 0.750, -0.75, 8)

print(f"S: {ship1.calc_S()}")
print(f"B/TA: {ship1.beam / ship1.draught}")
print(f"c9: {ship1.calc_c9()}")
print(f"CV: {ship1.calc_CV()}")
print(f"c11: {ship1.calc_c11()}")
Re = reynolds_number(ship1.length, ship1.speed) 
print(f"1 + k: {ship1.onePlusk1}")
print(f"CF: {0.075 / ((np.log10(Re) - 2)) ** 2}")
print(f"CA: {0.006 * (ship1.length + 100) ** (-0.16) - 0.00205}")
print(f"w: {ship1.wake_fraction()}")

print(f"t: {ship1.calc_thrust_deduction()}")
