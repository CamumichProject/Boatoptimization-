from physics import residual_resistance_coef, froude_number, reynolds_number, frictional_resistance_coef
import numpy as np


class Ship:
    """
    Class of ship object, can be initialize with zero argument.
    """

    def __init__(self, length: float, draught: float, beam: float, speed: float,
                 CB: float, midship_coefficient: float, 
                 waterplane_coefficient: float, lcb: float, prop_diameter: float) -> None:
        """
        Assign values for the main dimension of a ship.

        :param length: metres length of the vehicle
        :param draught: metres draught of the vehicle
        :param beam: metres beam of the vehicle
        :param speed: kts speed of the vehicle
        :param diaplacement: m^3 moulded displacement 
        :param prismatic_coefficient: Prismatic coefficient 
        """

        # Assumption: Ship is sitting EVEN KEEL (T_A = T_F)
        self.length = length
        self.draught = draught
        self.beam = beam
        self.CB = CB
        self.displacement = CB * length / 1.01675 * beam * draught
        print(f"displacement: {self.displacement}")
        self.speed = speed * 0.5144444 # now in m/s
        self.CS = self. length / self.displacement ** 0.33333 
        self.CM = midship_coefficient
        self.CP = self.CB / self.CM
        self.CWP = waterplane_coefficient
        self.ABT = 0 # 20 # REMOVE FOR LATER
        self.CA = 0.006 * (self.length + 100) ** (-0.16) - 0.00205
        self.surface_area = self.calc_S()
        self.lcb = lcb
        self.LR = (1 - self.CP + 0.06 * self.CP * self.lcb / (4 * self.CP - 1)) * self.length
        self.c12 = self.calc_c12()
        self.c13 = 1 # + 0.003 * 10 # REMOVE FOR LATER
        self.onePlusk1 = self.c13 * (0.93 + self.c12 * (self.beam / self.LR )**0.92497 * (0.95 - self.CP) ** -0.521448 * (1 - self.CP + 0.0225 * self.lcb) ** 0.6906)
        self.prop_diameter = prop_diameter
        self.eta_h = (1-self.calc_thrust_deduction()) / (1 - self.calc_wake_fraction())
        self.hb = 0. # 4.0 # remove for later
        self.resistance = self.calc_resistance()
    
        
# For reference: residual_resistance_coef(self.CS, self.CP, froude_number(self.speed, self.length))


    def calc_resistance(self) -> float:
        """
        Return resistance of the vehicle.

        :return: newton the resistance of the ship
        """
        total_resistance_coef = frictional_resistance_coef(self.length, self.speed) + \
                                self.wave_resistance_coef() + self.CA
        
        # print(f"Frictional Resistance Coeff: {frictional_resistance_coef(self.length, self.speed):.6f}")
        # print(f"Residual Resistance Coeff: {self.wave_resistance_coef():.6f}")
        # print(f"Correlation Allowance: {self.CA:.6f}")
        # print(f"Total Resistance Coeff: {total_resistance_coef:.6f}")

        return 1 / 2 * total_resistance_coef * 1025 * self.surface_area * self.speed ** 2

    def resistance_coeff(self):
        return frictional_resistance_coef(self.length, self.speed) + \
                                self.wave_resistance_coef() + self.CA

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

    
    def calc_wake_fraction(self) -> float:
        c9 = self.calc_c9()
        CV = self.calc_CV()
        c11 = self.calc_c11()
        c19 = self.calc_c19()
        c20 = 1
        CP1 = 1.45 * self.CP - 0.315 - 0.0225 * self.lcb

        term1 = c9 * c20 * CV * self.length / self.draught * (0.050776 + 0.93405 * c11 * CV / (1 - CP1))
        term2 = 0.27915 * c20 * np.sqrt(self.beam / (self.length * (1 - CP1)))
        term3 = c19 * c20
        
        return term1 + term2 + term3 

    def calc_c12(self) -> float:
        return (self.draught / self.length) ** 0.2228446
        
        
    def calc_CV(self) -> float:
        Re = reynolds_number(self.length, self.speed) 
        # From ITTC 1957 Friction Line
        CF = frictional_resistance_coef(self.length, self.speed)
        CA = self.CA
        return (self.onePlusk1 * CF) + CA

    def calc_c9(self) -> float:
        c8 = self.beam * self.surface_area / (self.length * self.prop_diameter * self.draught)
        return c8
        
    def calc_c11(self) -> float:

        return self.draught / self.prop_diameter



    def calc_c19(self) -> float:
        if self.CP < 0.7:
            return 0.12997/(0.95 - self.CB) - 0.11056 / (0.95 - self.CP)
        else:
            return 0.18567/(1.3571 - self.CM) - 0.71276 + 0.38648 * self.CP

    def calc_S(self) -> float:
        mult_factor = self.length * (2 * self.draught + self.beam) * np.sqrt(self.CM)
        group1 = 0.453 + 0.4425 * self.CB
        group2 = -0.2862 * self.CM - 0.003467 * self.beam / self.draught + 0.3696 * self.CWP 

        return mult_factor * (group1 + group2) #+ 2.38 * self.ABT / self.CB # REMOVE FOR LATER
    
    def calc_thrust_deduction(self):
        term1 = 0.25014 * (self.beam / self.length) ** 0.28956
        term2 = (np.sqrt(self.beam * self.draught) / self.prop_diameter) ** 0.2642
        term3 = (1 - self.CP + 0.0225 * self.lcb) ** 0.01762 
        term4 = 0.0015 * 10 # remove for later

        return term1 * term2 / term3 # + term4
    
    def circle_C(self):
        R = self.resistance / 4.448 # Resistance in pounds
        V = self.speed / 0.5144444   # Speed in knots
        EHP = R * V / 326            # Effective horsepower

        disp = self.displacement / 0.3048**3 * 64 / 2240 #displacement in tons
        return EHP / disp**(2/3) / V**3 * 427.1
    
    def circle_K(self):
        disp = self.displacement / 0.3048**3 * 64 / 2240 #displacement in tons
        return 0.5834 * self.speed / .51444 / disp**(1/6) 

    def wave_resistance_coef(self):

            c7 = self.beam / self.length

            eterm1 = -(self.length / self.beam)**0.80856 * (1-self.CWP)**0.30484
            eterm2 = (1-self.CP - 0.0225*self.lcb)**0.6367 * (self.LR / self.beam)**0.34574
            eterm3 = (100 * self.displacement / self.length**3) **0.16302
            ie = 1 + 89 * np.exp(eterm1 * eterm2 * eterm3)

            c1 = 2223105 * c7 ** 3.78613 * (self.draught / self.beam) ** 1.07961 * (90 - ie) ** (-1.37565)
            c3 = 0.56 * self.ABT ** 1.5 / (self.beam * self.draught * (0.31 * np.sqrt(self.ABT) + self.draught - self.hb)) 
            c2 = np.exp(-1.89 * np.sqrt(c3))

            AT = 0 #16 # remove for later
            c5 = 1 - 0.8 * AT / (self.beam * self.draught * self.CM)

            L = self.length
            B = self.beam
            T = self.draught
            

            lam = 1.446 * self.CP - 0.03 * L/B

            
            if self.CP < 0.8:
                c16 = 8.07981 * self.CP - 13.8673 * self.CP**2 + 6.984388 * self.CP**3
            else:
                c16 = 1.73014 - 0.7067*self.CP

            m1 = 0.0140407 * L/T - 1.75254 * self.displacement**(1/3) / L - 4.79323 * B/L - c16
            


            c15 = -1.69385 + (L/self.displacement**(1/3) - 8.0)/2.36

            Fn = froude_number(self.speed, L)
            m2 = c15 * self.CP**2 * np.exp(-0.1 * Fn**(-2))
            d = -0.9
            rho = 1025
            g = 9.81
            Rw = c1 * c2 * c5 * self.displacement * rho * g * np.exp(m1 * Fn**d + m2 * np.cos(lam * Fn**-2))
            return Rw / (0.5 * rho * self.speed**2 * self.surface_area)


