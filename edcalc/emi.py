"""
Functions calculating common-mode and differential-mode EMI filters.
"""

from typing import Optional, NamedTuple
import math

from .format import format_F, format_C, format_L, format_R, format_V


# noinspection PyPep8Naming
def dm_filter_frequency(F_sw, D, I_sw, ESR120, tcross):
    """
    Calculated DM filter's pole frequency.

    :param F_sw: Switching frequency of the power supply
    :param D: Duty cycle of the power supply
    :param I_sw: Center of the ramp of the current flowing through
        the main switches
    :param ESR120: ESR of the input bulk capacitor measured at 120Hz
    :param tcross: time it takes for the main switches to turn ON or OFF
    """
    nbreak1 = 1 / (math.pi * D)
    # nbreak2 = 1 / (math.pi * tcross * F_sw)
    fbreak1 = F_sw / (math.pi * D)
    # fbreak2 = 1 / (math.pi * tcross)
    c_estimated1 = 2 * I_sw / ((nbreak1 if fbreak1 > F_sw else 1) * math.pi)
    c_estimated2 = 2 * I_sw / (2 * math.pi)
    ESR = ESR120 / 2.25
    V_dm1 = 20 * math.log10((I_sw * c_estimated1 * ESR / 2) / 1e-6)
    V_dm2 = 20 * math.log10((I_sw * c_estimated2 * ESR / 2) / 1e-6)
    # print(V_dm1, V_dm2)
    attenuation = min(V_dm1, V_dm2) - 64
    fpole = 10 ** (math.log10(2 * F_sw) - (attenuation / 40))
    return fpole


# noinspection PyPep8Naming
def dm_dampening_RC_network(L: float, C: float, Q: float = 1.0):
    """
    Calculating the RC chain at the output of the DM filter, reducing its Q
    (and eliminating potential instability problems).
    watch https://youtu.be/k3gCiL6SFSE for more info.
    """
    C_damp = 5 * C
    R_damp = math.sqrt(L / C) / Q
    return R_damp, C_damp


# noinspection PyPep8Naming
class DMFilterParams(NamedTuple):
    # corner_frequency
    F_sw: float
    F: float
    V_in: float
    LC: float
    Z_out: float
    C_min: float

    L_max: float

    def print(self):
        print(f'DM filter at {format_V(self.V_in)} should:')
        print(f' - have corner frequency about {format_F(self.F)}')
        print(f' - have impedance of less than {format_R(self.Z_out)}')
        print(f' - use X capacitor of more than {format_C(self.C_min)}')
        print(f' - use total inductance (Lup + Ldown)'
              f'of less than {format_L(self.L_max)}')

    # noinspection PyPep8Naming
    def check_LC(self, L: float, C: float, L_source: float = 100e-6):
        print(f'Checking {format_L(L)}-{format_C(C)} '
              f'DM filter at {format_V(self.V_in)}')
        F = 1 / (2 * math.pi * math.sqrt(L * C))
        print(f'  corner frequency is {format_F(F)} vs calculated '
              f'{format_F(self.F)}')

        Z_out = math.sqrt(L / C)
        print('  The resulting impedance is', format_R(Z_out))
        if self.Z_out < Z_out:
            print('    which is more than maximum allowed',
                  format_R(self.Z_out), 'which can cause instability')
        else:
            print('    which is good, since it is less than maximum allowed',
                  format_R(self.Z_out))
        print('  optimum capacitance with the given L is',
              format_C(self.LC / L), 'vs', format_C(C), 'selected')
        print('  optimum inductance with the given C is',
              format_L(self.LC / C), 'vs', format_L(L), 'selected')
        R_damp, C_damp = dm_dampening_RC_network(L, C)
        print(f'  recommended dampening RC chain: {format_R(R_damp)} '
              f'resistor and {format_C(C_damp)} capacitor')
        # Input X-cap calculation, based on https://youtu.be/8M8B8GytW78
        L_total = 1 / (1 / L + 1 / L_source)
        C_pi_ub1 = C / 5  # ub stands for upper boundary
        C_pi_ub2 = 1 / (L_total * (4 * math.pi * self.F_sw)**2)
        C_pi_lb1 = 1 / (L_total * (math.pi * self.F_sw)**2)
        C_pi_lb2 = 1 / (10 * math.pi * self.F_sw)
        possible = 'NOT POSSIBLE' if C_pi_ub1 < C_pi_lb2 else 'check'
        print(f'  ({possible}) select Cpi < {format_C(C_pi_ub1)} '
              f'AND (Cpi << {format_C(C_pi_ub2)} '
              f'OR Cpi >> {format_C(C_pi_lb1)}) '
              f'AND Cpi > {format_C(C_pi_lb2)}')

    @staticmethod
    def create(F_sw, efficiency, P_out, ESR120,
               V_in, I_sw, D, T_rise) -> 'DMFilterParams':
        """
        :param F_sw: power supply's switching frequency
        :param efficiency: expected efficiency of the power supply.
        :param P_out: total output power of the power supply
        :param ESR120: ESR of the PSU's input capacitor measured at 120Hz
        :param V_in: Input voltage to verify (usually V_in_min or V_in_max)
        :param I_sw: Center of the ramp of the current flowing through
            the main switches.
        :param D: the lowest duty cycle of the power supply
        :param T_rise: time it takes for the main switches to turn ON or OFF
        """
        F_corner = dm_filter_frequency(F_sw, D, I_sw, ESR120, T_rise)
        LC = 1 / (2 * math.pi * F_corner) ** 2
        Z_in = efficiency * V_in**2 / P_out  # PSU input impedance
        Z_out = 0.1 * Z_in  # Filter output impedance
        C_min = 1 / (2 * math.pi * F_corner * Z_out)
        L_max = Z_out / (2 * math.pi * F_corner)
        return DMFilterParams(
            F_sw=F_sw,
            F=F_corner, V_in=V_in, LC=LC, Z_out=Z_out,
            C_min=C_min, L_max=L_max)


# noinspection PyPep8Naming
def cm_filter_corner_frequency(F_sw, D, V_amp, attenuation=None):
    Cp = 1e-10  # parasitic capacitance to earth
    Lqp = 65  # QP limit (65 dBµV) of IEC/EN61000-6-3
    # Breakpoints
    nbreak1 = 1 / (math.pi * D)
    # print('nbreak1:', nbreak1)
    fbreak1 = F_sw / (math.pi * D)
    # print('fbreak1:', fbreak1)
    c1 = 2 * V_amp / (nbreak1 * math.pi)
    c2 = 2 * V_amp / (2 * math.pi)
    # print('C1 =', format_V(c1))
    # print('C2 =', format_V(c2))
    # EMI spectrum with no filter
    I_cm1_mag = (2 * math.pi * fbreak1 * Cp * c1
                 / math.sqrt((50 * math.pi * Cp)**2 + 1))
    I_cm2_mag = (4 * math.pi * F_sw * Cp * c2
                 / math.sqrt((100 * math.pi * Cp)**2 + 1))
    # print('I_cm1_mag:', format_I(I_cm1_mag))
    # print('I_cm2_mag:', format_I(I_cm2_mag))
    assert abs(I_cm1_mag - I_cm2_mag) < 1e-5
    Vcm = 20 * math.log10(25 * I_cm1_mag/1e-6)
    # print(f'Vcm is {round(Vcm, 2)}dBuV')
    if attenuation is None:
        attenuation = Vcm - Lqp - 3
    f_c = F_sw * 10**(-attenuation / 40)
    return f_c


# noinspection PyPep8Naming
def cm_filter_check_LC(L: float, C: float, D: float, V_amp: float,
                       attenuation: Optional[float] = None):
    f_c = cm_filter_corner_frequency(D, V_amp, attenuation)
    LC_desired = 1 / (2 * math.pi * f_c)**2
    F_LC = 1 / (2 * math.pi * math.sqrt(L * C))
    print(f'Checking {format_L(L)}-{format_C(C)} CM filter:')
    print('  - the desired corner frequency is', format_F(f_c))
    print('  - the LC filter has frequency', format_F(F_LC),
          '(must be close to the desired)')
    print(f'  - the optimal C for the L={format_L(L)} '
          f'is {format_C(LC_desired / L)}')
    print(f'  - the optimal L for the C={format_C(C)} '
          f'is {format_L(LC_desired / C)}')

