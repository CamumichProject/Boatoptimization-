from physics import residual_resistance_coef, froude_number, reynolds_number, frictional_resistance_coef
import numpy as np


class Ship:
    """
    Class of ship object, can be initialize with zero argument.
    """

    def __init__(self, length: float, draught: float, beam: float, speed: float,
                 slenderness_coefficient: float, prismatic_coefficient: float,
                 CM: float, lcb: float, prop_diameter: float) -> None:
        """
        Assign values for the main dimension of a ship.

        :param length: metres length of the vehicle
        :param draught: metres draught of the vehicle
        :param beam: metres beam of the vehicle
        :param speed: m/s speed of the vehicle
        :param slenderness_coefficient: Slenderness coefficient dimensionless :math:`L/(∇^{1/3})` where L is length of ship,
            ∇ is displacement
        :param prismatic_coefficient: Prismatic coefficient dimensionless :math:`∇/(L\cdot A_m)` where L is length of ship,
            ∇ is displacement Am is midsection area of the ship
        :param midship coefficient, :math:`C_M = A_M / BT`
        :param lcb, LCB as a percentage forward of midship
        """
        self.length = length
        self.draught = draught
        self.beam = beam
        self.speed = speed
        self.slenderness_coefficient = slenderness_coefficient
        self.CP = prismatic_coefficient
        self.displacement = (self.length / self.slenderness_coefficient) ** 3
        self.CB = self.displacement / self.length / self.beam / self.draught
        self.CM = CM
        self.surface_area = 1.7 * self.length * self.draught + self.displacement / self.draught
        self.lcb = lcb
        self.LR = (1 - self.CP + 0.06 * self.CP * self.lcb / (4 * self.CP - 1)) * self.length
        self.c12 = self.calc_c12()
        self.onePlusk1 = 0.93 + self.c12 * (self.beam / self.LR )**0.92497 * (0.95 - self.CP) ** -0.521448 * (1 - self.CP + 0.0225 * self.lcb) ** 0.6906
        self.prop_diameter = prop_diameter




    @property
    def resistance(self) -> float:
        """
        Return resistance of the vehicle.

        :return: newton the resistance of the ship
        """
        total_resistance_coef = frictional_resistance_coef(self.length, self.speed) + \
                                residual_resistance_coef(self.slenderness_coefficient,
                                                         self.CP,
                                                         froude_number(self.speed, self.length))
        return 1 / 2 * total_resistance_coef * 1025 * self.surface_area * self.speed ** 2

    def resistance_coeff(self):
        return frictional_resistance_coef(self.length, self.speed) + \
                                residual_resistance_coef(self.slenderness_coefficient,
                                                         self.CP,
                                                         froude_number(self.speed, self.length))

    def maximum_deck_area(self, water_plane_coef: float = 0.88) -> float:
        """
        Return the maximum deck area of the ship

        :param water_plane_coef: optional water plane coefficient
        :return: Area of the deck
        """
        return self.beam * self.length * water_plane_coef

    @property
    def reynold_number(self) -> float:
        """
        Return Reynold number of the ship

        :return: Reynold number of the ship
        """
        return reynolds_number(self.length, self.speed)

    def propulsion_power(self, propulsion_eff: float = 0.7, sea_margin: float = 0.2) -> float:
        """
        Total propulsion power of the ship.

        :param propulsion_eff: Shaft efficiency of the ship
        :param sea_margin: Sea margin take account of interaction between ship and the sea, e.g. wave
        :return: Watts shaft propulsion power of the ship
        """
        return (1 + sea_margin) * self.resistance * self.speed / propulsion_eff

    def wake_fraction(self) -> float:
        c9 = self.calc_c9()
        CV = self.calc_CV()
        c11 = self.calc_c11()
        CP1 = 1.45 * self.CP - 0.315 - 0.0225 * self.lcb

        term1 = c9 * CV * self.length / self.draught * (0.661875 + 1.21756 * c11 * CV / (1 - CP1))
        term2 = 0.24588 * np.sqrt(self.beam / (self.length * (1 - CP1)))
        term3 = -0.09726 / (0.95 - self.CP)
        term4 = 0.11434 / (0.95 - self.CB)
        
        return term1 + term2 + term3 + term4

    def calc_c12(self) -> float:
        if self.draught / self.length > 0.05:
            return (self.draught / self.length) ** 0.2228446
        
        elif self.draught / self.length > 0.02:
            return 48.2 * (self.draught / self.length - 0.02) ** 2.078 + 0.479948
        
        else:
            return 0.479948
        
    def calc_CV(self) -> float:
        Re = reynolds_number(self.length, self.speed) 
        # From ITTC 1957 Friction Line
        CF = 0.075 / ((np.log10(Re) - 2)**2)
        CA = 0.006 * (self.length + 100) ** (-0.16) - 0.00205
        return (self.onePlusk1 * CF) + CA

    def calc_c9(self) -> float:
        if self.beam / self.draught < 5:
            c8 = self.beam * self.surface_area / (self.length * self.prop_diameter * self.draught)
        
        else:
            c8 = self.surface_area * (7 * self.beam / self.draught - 25) / (self.length * self.prop_diameter * (self.beam / self.draught - 3))

        if c8 < 28:
            return c8
        else:
            return 32 - 16 / (c8 - 24)
        
    def calc_c11(self) -> float:
        if self.draught / self.prop_diameter < 2:
            return self.draught / self.prop_diameter
        else:
            return 0.0833333 * (self.draught / self.prop_diameter) + 1.333333


