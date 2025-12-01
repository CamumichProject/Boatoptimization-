from ship import Ship
from physics import froude_number, reynolds_number, frictional_resistance_coef, to_meters
import numpy as np

# Testing out the ship class

ship1 = Ship(205.00, 10.00, 32.00, 25.0, 0.5716, 0.98, 0.750, -0.75, 8)

print(f"S: {ship1.calc_S()}")
print(f"B/TA: {ship1.beam / ship1.draught}")
print(f"c9: {ship1.calc_c9()}")
print(f"CV: {ship1.calc_CV()}")
print(f"c11: {ship1.calc_c11()}")
Re = reynolds_number(ship1.length, ship1.speed) 
print(f"1 + k: {ship1.onePlusk1}")
print(f"CF: {0.075 / ((np.log10(Re) - 2)) ** 2}")
print(f"CA: {0.006 * (ship1.length + 100) ** (-0.16) - 0.00205}")
print(f"w: {ship1.calc_wake_fraction()}")
print(f"t: {ship1.calc_thrust_deduction()}")
print(f"eta_h: {ship1.eta_h}")
print("All hull efficiency terms match with paper!")
print(f"CF calculated: {frictional_resistance_coef(ship1.speed, ship1.length)}")
print(f"Resistance {ship1.resistance / 1000}kN")

print("\n\n\nNow Working on Wave Resistance:")


print("\nNow working on Resistance\n")


# Comparing with PPP
L = to_meters(406.7)
B = to_meters(53.33)
T = to_meters(21.33)
disp = 7807 * 2240 / 64 * 0.3048**3

CB = 0.6
CM = 0.977
CWP = 0.624

LCB = -1.5

V = 16.79
D = to_meters(26.62)

ship2 = Ship(L, T, B, V, CB, CM, CWP, LCB, D)

print(f"calculated tonnage: {ship2.displacement / 0.3048**3 * 64 / 2240 :.1f}")
print(f"calculated resistance: {ship2.resistance:.1f}")
print(f"Circle C: {ship2.circle_C():.3f}")
print(f"Circle K: {ship2.circle_K():.1f}")
print(f"wake fraction: {ship2.calc_wake_fraction():.4f}")
print(f"thrust deduction: {ship2.calc_thrust_deduction():.4f}")
print(f"eta_h: {ship2.eta_h:.3f}")


L = to_meters(406.7)
B = to_meters(55.17)
T = to_meters(22.09)

CB = 0.650
CM = 0.982
CWP = 0.746

LCB = -0.5

V = 17.25
D = to_meters(26.40)
print("\n\n\n Starting Ship 3 \n \n\n")
ship3 = Ship(L, T, B, V, CB, CM, CWP, LCB, D)
print("\n\n\n")
print("Graph Info: ")
print(f"Circle K: {ship3.circle_K():.2f}")
print(f"B/H Ratio: {ship3.beam / ship3.draught :.1f}")
print(f"CB: {ship3.CB:.3f}")
print(f"L/B: {ship3.length / ship3.beam :.2f}")

print("\nResults:")

print(f"calculated tonnage: {ship3.displacement / 0.3048**3 * 64 / 2240 :.1f}")
print(f"calculated resistance: {ship3.resistance:.1f}")
print(f"calculated resistance coeff: {ship3.resistance_coeff()*1000:.3f}x10^3")
print(f"Circle C: {ship3.circle_C():.3f}")
print(f"wake fraction: {ship3.calc_wake_fraction():.4f}")
print(f"thrust deduction: {ship3.calc_thrust_deduction():.4f}")
print(f"eta_h: {ship3.eta_h:.3f}")