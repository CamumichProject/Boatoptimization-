import numpy as np
from ship import Ship
import matplotlib.pyplot as plt
from physics import to_meters


def get_series60_properties(target_Cb: float) -> tuple[float, float, float, float]:
    # 1. Define Arrays
    parents_Cb  = np.array([0.60,  0.65,  0.70,  0.75,  0.80])
    parents_Cm  = np.array([0.977, 0.982, 0.986, 0.992, 0.994])
    parents_Cwp = np.array([0.718, 0.747, 0.778, 0.817, 0.864])
    parents_lcb = np.array([-1.5, -0.5, 0.5,  1.5,  2.5]) 

    # 2. Interpolate and Convert to Native Python Float using .item()
    # np.interp returns numpy.float64 by default
    cm = np.interp(target_Cb, parents_Cb, parents_Cm).item()
    cwp = np.interp(target_Cb, parents_Cb, parents_Cwp).item()
    lcb = np.interp(target_Cb, parents_Cb, parents_lcb).item()
    
    # 3. Calculate Cp
    # Ensure this is also a native float
    cp = float(target_Cb / cm)
    
    return cm, cp, cwp, lcb

def get_V(circle_K: float, CB: float, length:float, beam: float, draught: float) -> float:
    displacement = CB * length / 1.01675 * beam * draught
    disp = displacement / 0.3048**3 * 64 / 2240
    
    return circle_K * disp**(1/6) / 0.5834




L = to_meters(406.7)
BT = 2.5
circle_K = 2.2
D = to_meters(25.3)

cbs = np.arange(.6, .7, .001)
lbs = np.arange(8.5, 5.5, -.01)

res_values = np.empty((len(lbs), len(cbs)))
eta_values = np.empty((len(lbs), len(cbs)))
res_coeffs = np.empty((len(lbs), len(cbs)))
RES_values = np.empty((len(lbs), len(cbs)))
cost_value = np.empty((len(lbs), len(cbs)))

for x, cb in enumerate(cbs):
    for y, lb in enumerate(lbs):
        CB = float(cb)
        LB = float(lb)

        B = L / LB
        T = B / BT
        
        CM, CP, CWP, LCB = get_series60_properties(CB)
        V = get_V(circle_K, CB, L, B, T)

        shp = Ship(L, T, B, V, CB, CM, CWP, LCB, D)

        res_values[y, x] = shp.circle_C()
        res_coeffs[y, x] = shp.resistance_coeff()
        eta_values[y, x] = shp.eta_h
        RES_values[y, x] = shp.resistance
        cost_value[y, x] = shp.circle_C() / shp.eta_h




CB_GRID, LB_GRID = np.meshgrid(cbs, lbs)
fig, axes = plt.subplots(nrows=1, ncols=3)


im1 = axes[0].contour(CB_GRID, LB_GRID, res_values, 12)
fig.colorbar(im1, ax=axes[0])

im2 = axes[1].contour(CB_GRID, LB_GRID, eta_values, 12)
fig.colorbar(im2, ax=axes[1])

axes[0].set_title("Curves of Froude Resistance Coefficient: 10 m/s")
axes[0].set_ylabel("L/B Ratio")
axes[0].set_xlabel("Block Coefficient")


axes[1].set_title("Curves of Hull Efficiency: 10 m/s")
axes[1].set_ylabel("L/B Ratio")
axes[1].set_xlabel("Block Coefficient")


im3 = axes[2].contour(CB_GRID, LB_GRID, cost_value, 12)
axes[2].set_ylabel("L/B Ratio")
axes[2].set_xlabel("Block Coefficient")
axes[2].set_title("Cost: resistance_coeff/eta_h")
fig.colorbar(im3, ax=axes[2])

plt.show()