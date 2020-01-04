"""
Data structures describing various components such as magnetic cores,
MOSFETs, voltage references, opto-couplers, etc.
"""
import math
from dataclasses import dataclass, replace
from typing import Optional, NamedTuple

from .format import block_of_values_plain, format_value_with_warning, format_W, \
    PrintableValues


class Core:
    """
    Describes a magnetic core (for a transformer or an inductor),
    given as many datasheet parameters as can be provided and trying
    to fill in (calculate) the rest.
    """

    # noinspection PyPep8Naming
    def __init__(
            self,
            # permeability
            mu: Optional[float] = None,
            # Effective magnetic path length, mm
            l_e: Optional[float] = None,
            # Effective core cross-sectional area, mm^2
            A_e: Optional[float] = None,
            # Effective core volume, mm^3
            V_e: Optional[float] = None,
            # nH / turns^2
            A_L: Optional[float] = None,
            # Max flux density before saturation, Tesla
            B_sat: Optional[float] = 0.28,
            # winding area of the bobbin, mm^2
            W_a: Optional[float] = None,
            note: str = ''):
        """
        :param mu: permeability of the material
        :param float l_e: Effective magnetic path, mm
        :param A_e: Effective cross-sectional area, mm^2
        :param V_e: Effective core volume, mm^3
        :param A_L: Inductance per square turn, nH / turns^2
        :param B_sat: Max flux density before the core saturates, Tesla
        :param W_a: Winding area of the bobbin, mm^2
        :param note: An optional comment
        """
        self.mu = mu
        self.l_e = l_e
        self.A_e = A_e
        self.V_e = V_e
        self.A_L = A_L
        self.B_sat = B_sat
        self.W_a = W_a
        self.note = note
        self.fill_in_blanks()

    def fill_in_blanks(self):
        for _ in range(2):
            if (self.V_e is None
                    and self.A_e is not None
                    and self.l_e is not None):
                self.V_e = self.A_e * self.l_e
            if (self.A_e is None
                    and self.V_e is not None
                    and self.l_e is not None):
                self.A_e = self.V_e / self.l_e
            if (self.l_e is None
                    and self.V_e is not None
                    and self.A_e is not None):
                self.l_e = self.V_e / self.A_e
            if (self.mu is None
                    and self.A_L is not None
                    and self.A_e is not None
                    and self.l_e is not None):
                self.mu = (self.A_L_mks * self.l_e_mks
                           / (4 * math.pi * 1e-7 * self.A_e_mks))
            if (self.A_L is None
                    and self.mu is not None
                    and self.A_e is not None
                    and self.l_e is not None):
                self.A_L = (
                        1e9 * (4 * math.pi * 1e-7 * self.mu * self.A_e_mks)
                        / self.l_e_mks)

    # noinspection PyPep8Naming
    @property
    def A_L_mks(self):
        """
        Inductance per square turn in MKS units (in henries / n^2)
        """
        return self.A_L / 1e9

    # noinspection PyPep8Naming
    @property
    def A_e_mks(self):
        """
        Effective cross-sectional area in MKS units (square meters)
        """
        return self.A_e / 1e6

    @property
    def l_e_mks(self):
        """
        Effective magnetic path in MKS units (meters)
        """
        return self.l_e / 1e3

    # noinspection PyPep8Naming
    @property
    def V_e_mks(self):
        """
        Effective core volume in MKS units (cubic meters)
        """
        return self.V_e / 1e9

    def __repr__(self):
        return (
            f'Core(mu={self.mu!r}, l_e={self.l_e!r}, A_e={self.A_e!r}, '
            f'V_e={self.V_e!r}, A_L={self.A_L!r}, B_sat={self.B_sat!r}, '
            f'W_a={self.W_a!r}, note={self.note!r})')

    def _repr_markdown_(self):
        def value_or_unknown(value, unit):
            return (
                f'{round(value, 2)} {unit}'
                if value is not None else 'unknown')

        return block_of_values_plain(
            ('Permeability', value_or_unknown(self.mu, '')),
            ('Effective magnetic path', value_or_unknown(self.l_e, '$mm$')),
            ('Effective cross-sectional area',
             value_or_unknown(self.A_e, '$mm^2$')),
            ('Effective core volume', value_or_unknown(self.V_e, '$mm^3$')),
            ('Inductance per square turn',
             value_or_unknown(self.A_L, '$nH / turns^2$')),
            ('Saturation flux density', value_or_unknown(self.B_sat, "T")),
            ('Winding area',
             value_or_unknown(self.W_a, '$mm^2$')),
            ('Note', self.note)
        )


class OptoCoupler(NamedTuple):
    # LED's forward voltage
    V_f: float
    # Collector-emitter saturation voltage
    V_sat: float
    # Minimum Current Transfer Ratio (CTR). For example,
    # at the end of the optocoupler's life
    CTR_min: float
    # Output parasitic capacitance, F
    # You can read on how to measure it from this presentation:
    # [The TL431 in the Control of Switching Power Supplies]
    # (https://www.onsemi.com/pub/Collateral/TND381-D.PDF)
    C_out: float


class ShuntReference(NamedTuple):
    """
    Shunt voltage reference, such as an adjustable TL431 or even a zener diode
    """
    # Recommended bias current necessary to guarantee accurate regulation, A
    I_bias: float
    # Minimum regulated cathode voltage, V
    V_min: float
    # Voltage at the reference pin (only for adjustable references like TL431)
    V_ref: Optional[float] = None


# Some tests


def test_core_math():
    from pytest import approx
    core = Core(mu=2300, A_e=119, l_e=46.3, V_e=5490)
    assert core.mu == 2300
    assert core.A_L == approx(7428.54)

    core = Core(A_L=5800, A_e=119, l_e=46.3, V_e=5490)
    assert core.mu == approx(1795.776)


# noinspection PyPep8Naming
@dataclass(frozen=True)
class HeatingComponent:
    # Junction-to-Case Thermal Resistance, °C/W
    Rt_JC: float
    # Case-to-Ambient Thermal Resistance (JA - JC), °C/W
    Rt_CA: float
    # Maximum junction temperature, °C
    T_j: float

    def thermal_analysis(self, I_rms: float, Rt_CS: float = None,
                         Rt_SA: float = None,
                         T_ambient: float = 25) -> PrintableValues:
        result = []
        dissipated_power = self.dissipated_power(I_rms)
        result.append(('Power being dissipated, per component',
                       format_W(dissipated_power)))
        T_no_sink = T_ambient + dissipated_power * (self.Rt_JC + self.Rt_CA)
        result.append(
            ('Junction temperature without a heatsink',
             format_value_with_warning(T_no_sink, '°C', self.T_j,
                                       'Too hot (>= {max_value})!')))
        if Rt_CS is not None and Rt_SA is not None:
            T_heatsink = (
                T_ambient
                + dissipated_power * (self.Rt_JC + Rt_CS + Rt_SA))
            result.append(
                ('Junction temperature with the heatsink',
                 format_value_with_warning(T_heatsink, '°C', self.T_j,
                                           'Too hot (>= {max_value})!')))
        return result

    def dissipated_power(self, I_rms):
        raise NotImplementedError()


# noinspection PyPep8Naming
@dataclass(frozen=True)
class MOSFET(HeatingComponent):
    # Reverse Transfer Capacitance, F
    C_rss: float
    # Input Capacitance, F
    C_iss: float
    # Static Drain-to-Source On-Resistance, Ohms
    R_ds: float

    def gate_drive_charge(self, V_d: float, V_gs: float) -> float:
        """
        Calculates the charge we need to bring the gate to V_gs,
        assuming the drain has potential V_d and given C_rss/C_iss
        parameters from the datasheet.
        Knowing the charge we can later calculate the current we need
        to deliver this charge in time given.
        The theory is explained here:
        https://youtu.be/of_v2N5f788
        http://www.ti.com/lit/pdf/slua618

        :param V_d: initial potential of the drain
        :param V_gs: target voltage between the gate and the source
            we need to reach.
        """
        C_gd = self.C_rss
        C_gs = self.C_iss - self.C_rss
        C_equiv = C_gs + C_gd * (1 + V_d / V_gs)
        Q_gs = V_gs * C_equiv
        return Q_gs

    def dissipated_power(self, I_rms: float):
        return I_rms**2 * self.R_ds


# noinspection PyPep8Naming
@dataclass(frozen=True)
class Diode(HeatingComponent):
    # Maximum DC Blocking Voltage, V
    V_max: float
    # Forward Voltage drop on each of the output rectifying diodes, V
    V_drop: float
    # Maximum average forward current, A
    I_avg: float
    # Peak forward surge current, A
    I_peak: float
    # Number of diodes connected in parallel
    num_in_parallel: int = 1

    def dissipated_power(self, I_rms: float):
        return (I_rms / self.num_in_parallel) * self.V_drop

    def parallel(self, number: int):
        return replace(self, num_in_parallel=number)

    def check_conditions(self, V_max: float,
                         I_avg: float,
                         I_pk: float) -> PrintableValues:
        return [
            (f'Max reverse voltage on each of the diodes',
             format_value_with_warning(V_max, 'V', self.V_max)),
            (f'Peak current through each diode',
             format_value_with_warning(
                 I_pk, 'A', self.I_peak / self.num_in_parallel)),
            ('Average current per diode',
             format_value_with_warning(
                 I_avg, 'A',
                 self.I_avg / self.num_in_parallel)),
        ]
