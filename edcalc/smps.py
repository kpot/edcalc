"""
This is a bunch of functions especially useful for designing various power 
supplies using Jupyter.
"""

import math
from typing import Optional, List, Union, Tuple, NamedTuple

from lcapy import Circuit

from .component_types import Core, OptoCoupler, ShuntReference
from .eseries import nearest_e12, nearest_e24
from .format import (
    format_I, format_markdown_table, format_F, format_V, PrintableValues,
    format_C, format_R, format_L, format_W, warning_message,
    format_flux_density, block_of_values)


class OutputParams(NamedTuple):
    # Output voltage, not taking the diode drop into account
    V: float
    # Max output current
    I: float
    # Type of rectification: center-tapped, full-bridge or None
    rectifier: Optional[str] = 'center-tapped'
    # Voltage drop on each of the output rectifying diodes, V
    V_D = 1
    # Combined (junction-to-case-to-heatsink-to-ambient) temperature
    # resistances of the output diodes
    Rt_JA_D: Optional[float] = None
    # Magnetic core that is to be used for the output filter
    filter_core: Optional[Core] = None
    # Output filter capacitors. Use None if the value is to be calculated
    # automatically
    C: Optional[float] = None
    # Output filter inductance. Use None if the value is to be calculated
    # automatically
    L: Optional[float] = None
    # ESRs of the output capacitors. Use None if the value is to be
    # determined automatically
    ESR: Optional[float] = None
    # Desired voltage ripple on output capacitors (portion of the input
    # voltage)
    C_ripple: float = 0.001

    def full_diode_drop(self):
        if self.rectifier == 'center-tapped':
            return self.V_D
        elif self.rectifier == 'full-bridge':
            return 2 * self.V_D
        elif self.rectifier is None:
            return 0
        else:
            raise ValueError(
                f'Unknown type of rectification: {self.rectifier!r}')

    @property
    def is_center_tapped(self):
        return self.rectifier == 'center-tapped'


class CoilChoice(NamedTuple):
    wire_dia: float
    is_too_thick: bool
    copper_dia: float
    strands: int
    total_copper_area: float
    winding_area_fill: int


class CoilSelectionInfo(NamedTuple):
    """
    Basic info regarding the winding selection
    """
    # Required switching frequency
    frequency: float
    # Number of turns necessary
    turns: int
    # Expected RMS current, Amperes
    rms_current: float
    # Maximum current density, Amperes / mm^2
    current_density: float
    # How deep from the surface of the wire the current is still flowing, mm
    skin_depth: float
    # Copper area required for the current density, mm^2
    required_copper_area: float
    # Maximum strand dia given the chosen frequency and the skin effect, mm
    optimal_strand_dia: float
    # Maximum strand's copper area given the chosen frequency, mm^2
    optimal_strand_area: float
    # Possible choices for different wires
    choices: List[CoilChoice]

    def _repr_markdown_(self):
        basic_info = (
            f'With the RMS current of {format_I(self.rms_current)} '
            f'and max current density of {self.current_density} $A/mm^2$, '
            f'we need at least {round(self.required_copper_area, 2)} $mm^2$ '
            f'of cross-sectional copper area available for it. \n\n'
            f'Thus the optimum wire dia is: '
            f'{round(2 * self.skin_depth, 2)} $ mm $ '
            f'({round(math.pi * self.skin_depth**2, 2)} $mm^2$ each strand)'
        )
        table_header = (
            'wire Ø', 'copper Ø', 'fully utilized', 'strands', 'copper total',
            f'{self.turns} turns occupy'
        )
        table_rows = [
            [f'{c.wire_dia} $mm$',
             f'{c.copper_dia} $mm$',
             f'{"TOO THICK" if c.is_too_thick else "OK"}',
             f'{c.strands}',
             f'{round(c.total_copper_area, 4)} $mm^2$',
             (f'{c.winding_area_fill}%' if c.winding_area_fill is not None
              else '')]
            for c in self.choices
        ]
        return (basic_info
                + '\n\n'
                + format_markdown_table(table_header, table_rows))


def coil_selection(
            frequency: float,
            rms_current: float,
        wires_dia: List[Union[float, Tuple[float, float]]],
        current_density: float,
        turns: int,
        winding_area: Optional[float] = None,
        window_utilization: float = 0.3) -> CoilSelectionInfo:
    """
    Analyze different options of how to make a winding and summarizes them
    into a table.

    :param frequency: switching frequency, Hz
    :param rms_current: average current flowing through the coil, A
    :param wires_dia: a list of all diameters of wires to choose from, mm
       Also, for heavily insulated wires a tuple of
       (<wire diameter>, <copper core diameter>) can be used.
    :param current_density: recommended current density, A
    :param winding_area: Winding area of the bobbin
    :param window_utilization: window utilization factor, usually around 0.3.
         Reflects how much of the winding area is actually occupied by copper.
    :param turns: the number of turns needed
    """
    skin_depth_mm = 66.1 / math.sqrt(frequency)
    required_copper_area = rms_current / current_density
    choices = []
    for wd in wires_dia:
        if isinstance(wd, tuple):
            # dividing into wire dia and copper diameters
            wire_dia, copper_dia = wd
        else:
            wire_dia = copper_dia = wd
        # noinspection PyPep8Naming
        A_cu = (math.pi * copper_dia**2 / 4)
        # noinspection PyPep8Naming
        A_w = (math.pi * wire_dia**2 / 4)
        # How much of the copper is unavailable due to the skin effect
        unavailable_area = (
            (math.pi * (copper_dia/2 - skin_depth_mm)**2)
            if copper_dia > 2 * skin_depth_mm
            else 0)
        copper_area_available = A_cu - unavailable_area
        wire_strands = max(
            1,
            math.ceil(required_copper_area / copper_area_available))
        winding_area_occupied = A_w * wire_strands * turns / window_utilization
        winding_choice = CoilChoice(
            wire_dia=wire_dia,
            is_too_thick=copper_dia > 2 * skin_depth_mm,
            copper_dia=copper_dia,
            strands=wire_strands,
            total_copper_area=wire_strands * copper_area_available,
            winding_area_fill=(
                math.ceil(100 * winding_area_occupied / winding_area)
                if winding_area is not None
                else None)
        )
        choices.append(winding_choice)

    result = CoilSelectionInfo(
        frequency=frequency,
        turns=turns,
        rms_current=rms_current,
        current_density=current_density,
        skin_depth=skin_depth_mm,
        required_copper_area=required_copper_area,
        optimal_strand_dia=2 * skin_depth_mm,
        optimal_strand_area=math.pi * skin_depth_mm**2,
        choices=choices)
    return result


# noinspection PyPep8Naming
class CurrentTransformerNetwork(NamedTuple):
    I_in: float
    I_in_rms: float
    V_out: float
    primary_turns: int
    secondary_turns: int
    diode_drop: float
    R_burden_split: float
    P_burden_split: float
    B_pk: float
    core: Core

    # noinspection PyPep8Naming
    @staticmethod
    def create(I_in: float, V_out: float, I_in_rms: float, core: Core,
               diode_drop: float = 0.5, primary_turns: int = 1,
               secondary_turns: int = 100) -> 'CurrentTransformerNetwork':
        """
        Calculates all components of a transformer-based current sensing
        circuit.

        :param I_in: peak input current
        :param V_out: peak output voltage the input current translates into
        :param I_in_rms: input RMS current, so we could estimate power
            dissipated by the burden resistors.
        :param core: magnetic core being used (to check if the core saturates)
        :param diode_drop: voltage drop on each of the rectifying diodes
            (this function assumes that a full-bridge rectifier is used).
        :param primary_turns: the number of turns in the primary winding
            of the current transformer
        :param secondary_turns: the number of turns in the secondary winding
            of the current transformer
        """
        transformer_ratio = secondary_turns / primary_turns
        secondary_current = (I_in / transformer_ratio)
        R_burden = (V_out + 2 * diode_drop) / secondary_current
        B_pk = I_in * core.A_L_mks * primary_turns / core.A_e_mks
        R_burden_split = nearest_e24(R_burden * 2)
        V_out_real = R_burden * secondary_current - 2 * diode_drop
        return CurrentTransformerNetwork(
            I_in=I_in,
            I_in_rms=I_in_rms,
            V_out=V_out_real,
            primary_turns=primary_turns,
            secondary_turns=secondary_turns,
            diode_drop=diode_drop,
            R_burden_split=R_burden_split,
            P_burden_split=((I_in_rms / transformer_ratio)**2 * R_burden) / 2,
            B_pk=B_pk,
            core=core)

    def markdown(self):
        transformer_ratio = self.secondary_turns / self.primary_turns
        V_out_rms = (
            (self.R_burden_split / 2)
            * (self.I_in_rms / transformer_ratio)
            - 2 * self.diode_drop)
        diode_rating = (
            (self.R_burden_split * self.I_in / transformer_ratio) / 2)
        return block_of_values(
            ('Primary winding turns', str(self.primary_turns)),
            ('Secondary winding turns', str(self.secondary_turns)),
            ('Primary side inductance',
             format_L(self.primary_turns**2 * self.core.A_L_mks)),
            ('Secondary side inductance',
             format_L(self.secondary_turns**2 * self.core.A_L_mks)),
            ("First and second burden resistors",
             format_R(self.R_burden_split)),
            ("Power dissipated by each resistor",
             format_W(self.P_burden_split)),
            ("Peak magnetic flux",
             format_flux_density(self.B_pk, self.core.B_sat)),
            ('Peak input current', format_I(self.I_in)),
            ('Peak output voltage', format_V(self.V_out)),
            ('Output voltage matching RMS current', format_V(V_out_rms)),
            ('Max voltage the diodes should be able to withstand',
             format_V(diode_rating)),
        )

    def draw(self):
        """
        Draws in Jupyter schematics of a TL431-based feedback loop circuit.
        """
        cct = Circuit(
            f"Prectified Vin1 Vin2; down, l=input\n"
            f"TF1 t1 t2 t3 t4 core; right\n"
            f"W Vin1 t3; right, i=$I_{{in}}$\n"
            f"W Vin2 t4; right\n"
            f"W t1 t1e; up=0.8\n"
            f"W t2 t2e; down=0.8\n"
            f"W t1e rb1_1; right\n"
            f"W t2e rb1_2; right\n"
            f"Rb1 rb1_1 rb1_2 {self.R_burden_split}; down\n"
            f"W rb1_1 dbin1; right=2\n"
            f"W rb1_2 dbin2; right=2\n"
            f"D1 dbin1 dboutp; rotate=225, scale=0.5\n"
            f"D2 dboutn dbin1; rotate=135, scale=0.5\n"
            f"D3 dbin2 dboutp; rotate=135, scale=0.5\n"
            f"D4 dboutn dbin2; rotate=225, scale=0.5\n"
            f"Rb2 dboutp dboutn {self.R_burden_split}; right, v=$V_{{out}}$\n"
        )
        cct.draw(scale=0.5, style='american', draw_nodes='connections',
                 label_nodes=False)


def tan2esr(tan_delta: float, capacitance: float, frequency: float = 120):
    """
    For capacitors: converts tan(sigma) to ESR
    """
    return tan_delta / (2 * math.pi * capacitance * frequency)


class ReferenceFeedbackBase(NamedTuple):
    """
    Describes the foundation of a TL431-based feedback network without
    the integrator-related components.
    """
    # The rail voltage the output of the optocoupler is being pulled to
    V_dd: float
    # A resistance pulling optocoupler's output to the rail (V_dd or ground).
    R_opto_pull: float
    # resistor providing current to the fixed reference
    R_fixed_ref: float
    # Limits the current through both the LED and the main reference
    R_LED: float
    # Biasing resistor around the LED
    R_bias: float
    # max power dissipated by the fixed reference
    P_z: float
    # The optocoupler used
    opto: OptoCoupler
    # The main voltage reference used
    reference: ShuntReference
    # A regulator providing a fixed DC point to power the main reference
    fixed_reference: ShuntReference

    # noinspection PyPep8Naming
    @classmethod
    def create(cls, supply_voltage: float,
               opto: OptoCoupler,
               main_reference: ShuntReference,
               fixed_reference: ShuntReference,
               V_dd: float = 5.0,
               R_opto_pull: float = 10e3) -> 'ReferenceFeedbackBase':
        """
        Calculates base elements of the TL431-like-based feedback circuit
        setting all necessary currents, all without the fast-lane.
        For more details, read ["The TL431 in Switch-Mode Power Supplies loops"
        by C. Basso](https://cbasso.pagesperso-orange.fr/Downloads/Papers/
        The%20TL431%20in%20loop%20control.pdf)

        :param supply_voltage: Max voltage the whole circuit is feeding from.
            Usually it's the output voltage of the power supply.
            If this voltage can vary, use its lowest value but check that
            the components can survive the stress of the highest.
        :param opto: Optocoupler being used to link the output with
            the PWM controller.
        :param main_reference: The main reference controlling the output
            voltage of the converter.
        :param fixed_reference: A regulator providing a fixed DC point
            to power both the main reference and the optocoupler, thus removing
            the fast-lane problem described by C. Basso.
        :param V_dd: The rail voltage the output of the optocoupler is being
            pulled to.
        :param R_opto_pull: A resistance pulling optocoupler's output
            to the rail (V_dd or ground).
        """
        # 'R_bias' is the biasing resistor around the LED
        R_bias = nearest_e12(opto.V_f / main_reference.I_bias, 'lower')
        I_bias_actual = opto.V_f / R_bias
        I_LED_max = (R_opto_pull - opto.V_sat) / (R_opto_pull * opto.CTR_min)
        # x10 to lower the impedance of the zener-based voltage source
        I_R_z = fixed_reference.I_bias + 10 * (I_bias_actual + I_LED_max)
        # 'R_z' is the resistor providing current to the fixed reference
        R_z = nearest_e12((supply_voltage - fixed_reference.V_min) / I_R_z)
        R_LED_reduction = 0.5
        # 'R_LED' limits the current through both the LED
        # and the main reference
        R_LED = nearest_e12(
            R_LED_reduction
            * ((fixed_reference.V_min - opto.V_f - main_reference.V_min)
               / (V_dd
                  - opto.V_sat
                  + main_reference.I_bias * opto.CTR_min * R_opto_pull))
            * R_opto_pull * opto.CTR_min)
        return ReferenceFeedbackBase(
            V_dd=V_dd,
            R_opto_pull=R_opto_pull,
            R_fixed_ref=R_z,
            R_LED=R_LED,
            R_bias=R_bias,
            # max power dissipated by the fixed reference
            P_z=I_R_z * fixed_reference.V_min,
            opto=opto,
            reference=main_reference,
            fixed_reference=fixed_reference)

    # noinspection PyPep8Naming
    def buck_type3_compensator(
            self, V_in, V_o, L_o, C_o,
            ESR, F_sw, V_ramp, I_ref_in=250e-6,
            split_poles_and_zeros=False) -> Tuple['ReferenceFeedback',
                                                  PrintableValues]:
        """
        Calculates all components for a Type III compensator of
        a buck converter, usually needed for voltage-mode control.
        The procedure is based on "The TL431 in Switch-Mode Power Supplies
        loops: part IV" by C. Basso and on "AN-1162 Compensator Design
        Procedure for Buck Converter with Voltage-Mode Error-Amplifier"
        by Amir M. Rahimi et al.

        :param V_in: Input voltage of the converter as it's seen by its
            main stage. For example, for a forward converter
            we can ignore the transformer and pretend that the input voltage
            is what comes out of it.
        :param V_o: Target output voltage of the converter we need to achieve.
        :param L_o: Output filter inductor.
        :param C_o: Output filter capacitor.
        :param ESR: ESR of the output filter capacitor
        :param F_sw: switching frequency
        :param V_ramp: amplitude of the control signal of the PWM controller
        :param I_ref_in: How much current the main voltage divider should
            provide to the reference pin of the main voltage reference (TL431)
            to guarantee accurate regulation.
        :param split_poles_and_zeros: Whether we should co-locate the two poles
            and two zeroes (False) or try to assign separate values
            to each of them (True).

        :returns: a tuple of two values - complete description of the network
            and a block of additional information for printing.
        """
        if self.reference.V_ref is None:
            raise ValueError('Not a programmable reference')
        # Based on "The TL431 in Switch-Mode Power Supplies loops: part IV"
        # by C. Basso
        F_LC = 1 / (2 * math.pi * math.sqrt(L_o * C_o))
        F_ESR = 1 / (2 * math.pi * C_o * ESR)
        fc = F_cross = F_sw / 6
        Rlower_fb = nearest_e24(self.reference.V_ref / I_ref_in)
        R1_fb = nearest_e24((V_o - self.reference.V_ref) / I_ref_in)
        # Static gain from the basic TL431-optocoupler network
        G0 = self.R_opto_pull * self.opto.CTR_min / self.R_LED
        # DC gain of the whole converter
        G = V_ramp / V_in
        # DC gain of the compensator itself
        G1 = G / G0
        if split_poles_and_zeros:
            if F_ESR < F_sw / 2:
                compensator_type = 'Type III A'
                fz1 = 0.75 * F_LC
                fz2 = F_LC
                fp1 = F_ESR
                fp2 = F_sw / 2
            else:
                compensator_type = 'Type III B'
                phase_lead = 70 * math.pi / 180
                fz2 = F_cross * math.sqrt(
                    (1 - math.sin(phase_lead)) / (1 + math.sin(phase_lead)))
                fz1 = 0.5 * fz2
                fp1 = F_cross * math.sqrt(
                    (1 + math.sin(phase_lead)) / (1 - math.sin(phase_lead)))
                fp2 = F_sw / 2
            C3b = nearest_e12(1 / (2 * math.pi * R1_fb * fz1))
            R3b = nearest_e12(1 / (2 * math.pi * fp2 * C3b))
            # Formula (3-50) from the book "SMPS SPICE Simulations
            # and Practical Designs" by C. Basso.
            R2b = nearest_e12(
                math.sqrt((fp1 ** 2 + fc ** 2) * (fp2 ** 2 + fc ** 2) / (
                            (fz1 ** 2 + fc ** 2) * (fz2 ** 2 + fc ** 2)))
                * (G1 * fc * R3b / fp1))
            C1b = nearest_e12(1 / (2 * math.pi * fz2 * R2b))
            C2b = nearest_e12(
                max(0,
                    (1 / (2 * math.pi * fp1 * self.R_opto_pull))
                    - self.opto.C_out))

        else:
            compensator_type = 'Type III'
            fz = fz1 = fz2 = F_LC
            fp = fp1 = fp2 = F_ESR
            C3b = nearest_e12(1 / (2 * math.pi * R1_fb * fz))
            R3b = nearest_e12(1 / (2 * math.pi * fp * C3b))
            R2b = nearest_e12(
                ((fp ** 2 + fc ** 2) / (fz ** 2 + fc ** 2)) * (
                            G1 * fc * R3b / fp))
            C1b = nearest_e12(1 / (2 * math.pi * fz * R2b))
            C2b = nearest_e12(
                max(0,
                    (1 / (2 * math.pi * fp * self.R_opto_pull))
                    - self.opto.C_out))
        V_out1_max = self.reference.V_ref * (Rlower_fb + R1_fb) / Rlower_fb

        result = ReferenceFeedback(
            base=self,
            C_out=C_o,
            L_out=L_o,
            R_divider_upper=R1_fb,
            R_divider_lower=Rlower_fb,
            # Feedback loop RC network
            C1=C1b, R2=R2b,
            # RC network parallel to the divider's upper resistor
            # Note: the seeming confusion between C2 and C3 is not a mistake,
            # just Basso and Maniktala identify components differently
            C2=C3b, R3=R3b,
            # Parallel to the optocoupler's output
            C3=C2b)
        extra_info = [
            ('Crossover frequency', format_F(F_cross)),
            ('LC filter frequency', format_F(F_LC)),
            ('ESR frequency', format_F(F_ESR)),
            ('DC plant gain', str(round(G, 4))),
            ('Max regulated output voltage', format_V(V_out1_max)),
            ('Compensator type', compensator_type),
            ('$F_{z1}$', format_F(fz1)),
            ('$F_{z2}$', format_F(fz2)),
            ('$F_{p1}$', format_F(fp1)),
            ('$F_{p2}$', format_F(fp2)),
        ]
        return result, extra_info


class ReferenceFeedback(NamedTuple):
    """
    A complete shunt reference-based feedback network, relying on TL431
    or one of its analogs.
    """
    base: ReferenceFeedbackBase
    C_out: float
    L_out: float
    R_divider_upper: float
    R_divider_lower: float
    # RC network between the cathode and the ref pins
    C1: float
    # RC network between the cathode and the ref pins
    R2: float
    # RC network parallel to the divider's upper resistor
    C2: float
    # RC network parallel to the divider's upper resistor
    R3: float
    # Parallel to the optocoupler's output
    C3: float

    def draw(self):
        """
        Draws in Jupyter schematics of a TL431-based feedback loop circuit.
        """
        cct = Circuit(
            f"Prectified Vin1 Vin2; down, l=V_{{IN}}\n"
            f"W Vin1 1_30; right\n"
            f"W Vin2 0_7; down, sground\n"
            f"W 0_cout 0_cout2; down=0.25,sground\n"
            f"Lout 1_30 1_31 {self.L_out}; right\n"
            f"W 1_31 1_3; right=8\n"
            f"Cout 1_31 0_cout {self.C_out}; down\n"
            f"U1 chip1313; pindefs={{ref=out,cathode=t1,anode=b1}}, "
            f"pinlabels={{ref=REF,cathode=Cathode,anode=Anode}},l=TL431\n"
            f"DTL431 U1.anode U1.cathode; up,size=2,scale=0.5,kind=zener,l=,"
            f"  color=gray,fixed\n"
            f"W U1.anode 0_6; down=0.1, sground\n"
            f"W U1.cathode u1cat_j; up\n"
            f"U2 chip2121; pindefs={{anode=r1,cathode=r2,c=l1,e=l2}}, "
            f"  pinlabels={{anode=Anode,cathode=Cathode,"
            f"    c=Collector,e=Emitter}},"
            f"  l=Optocoupler\n"
            f"DOpto u2d_1 u2d_2; down, scale=0.5, l=,color=gray,kind=led,"
            f"fixed\n"
            f"W u2d_2 U2.cathode; right=0.25, dashed, fixed, color=gray\n"
            f"W U2.anode u2d_1; left=0.25, dashed, fixed, color=gray\n"
            f"W U2.anode u2a_j1; right=0.5\n"
            f"W u2a_j1 u2a_j2; up\n"
            f"W u2a_j2 u2a_j3; right=0.5\n"
            f"RLED u2a_j3 rled2 {self.base.R_LED}; up\n"
            f"W rled2 zen1; left=1.5\n"
            f"W zen1 zcap1; left=1.5\n"
            f"DZ zen2 zen1; up,kind=zener,"
            f"  l={self.base.fixed_reference.V_min}V\n"
            f"CZ zcap1 zcap2 100e-9; down\n"
            f"W zen2 zcap2; left\n"
            f"W zen2 zen2g; down=0.1,sground\n"
            f"RZ rled2 rup2 {self.base.R_fixed_ref}; right\n"
            f"W u1cat_j U2.cathode; left\n"
            f"Rbias u1cat_j u2a_j3 {self.base.R_bias}; up\n"
            f"CF1 u1cat_j cf1_2 {self.C1}; right\n"
            f"RF2 cf1_2 r2_j {self.R2}; right\n"
            f"W r2_j rlj; down\n"
            f"W U1.ref rlj; right\n"
            f"Rlower rlj rl2 {self.R_divider_lower}; down\n"
            f"W rl2 0_5; down=0.1, sground\n"
            f"RF1 r2_j rup2 {self.R_divider_upper}; up=2\n"
            f"W rup2 1_3; up\n"
            f"W r2_j fastlane1; right\n"
            f"CF2 fastlane1 fastlane2 {self.C2}; up\n"
            f"RF3 fastlane2 fastlane3 {self.R3}; up\n"
            f"W fastlane3 rup2; left\n"
            f"W U2.e u2e_j; left\n"
            f"Rpulldown u2e_j r_optopull2 {self.base.R_opto_pull}; down\n"
            f"U3 chip4141; "
            f"  pindefs={{ref=r1,inplus=r2,inminus=r3,fb=r4,gnd=vss}},"
            f"  pinlabels={{ref=REF,inplus=IN+,inminus=IN-,fb=FB,gnd=GND}},"
            f"  l=TL494\n"
            f"W U2.c cf3_1; left\n"
            f"W cf3_1 U3.ref; left\n"
            f"W u2e_j u2e_j2; left\n"
            f"W u2e_j2 u2e_j3; up=0.5\n"
            f"W u2e_j3 U3.inplus; left\n"
            f"W r_optopull2 pwmc_gnd1; down\n"
            f"W U3.gnd pwmc_gnd2; down\n"
            f"W pwmc_gnd1 pwmc_gnd2; left\n"
            f"W U3.inminus inmin1; right=0.5\n"
            f"W inmin1 inmin2; down=0.5\n"
            f"W inmin2 U3.fb; left=0.5\n"
            f"CF3 cf3_1 u2e_j {self.C3}; down\n"
        )
        cct.draw(scale=0.5, style='american', draw_nodes='connections',
                 label_nodes=False)


class BuckOutputParams(NamedTuple):
    """
    Recommended output parameters for a buck converter
    """
    L: float
    N: float
    C: float
    ESR: float
    core_fits: bool


# noinspection PyPep8Naming
def buck_output_params(
        F_sw: float, V_in_min: float, V_in_max: float, I_pft: float,
        D_min: float, D_max: float, V_o: float,
        I_o: float, ripple: float, core: Core,
        C_ripple: float) -> Tuple[BuckOutputParams, PrintableValues]:
    """
    For a buck converter, calculates parameters of the inductor
    (the inductance, the number of turns and so on), and the output capacitor.
    The same function can be used for any other buck-derived topology such
    as forward, push-pull or half-bridge converters.

    :param F_sw: Switching frequency
    :param V_in_min: Minimum input voltage as it is seen by the inductor.
        For isolated converter the voltage must be scaled by the turn ratio
        of the transformer.
    :param V_in_max: Maximum input voltage as it is seen by the inductor.
        Similar to V_in_min.
    :param I_pft: peak current flowing through the inductor.
    :param D_min: minimum duty cycle
    :param D_max: maximum duty cycle
    :param V_o: (regulated) output voltage
    :param I_o: output current
    :param ripple: current ripple ratio (dI / I_L)
    :param core: the magnetic core that is supposed to be used for the inductor
    :param C_ripple: maximum output voltage ripple ratio (dV / V)
    :return: returns a structure describing the component values,
        as well as some extra information to be displayed.
    """
    # these are desired values assuming the core can handle
    # the output current without saturation and thus works in CCM
    dI = ripple * I_o
    dV_out = V_o * C_ripple
    L_out = V_in_max * (1 - D_min) / (I_o * ripple * F_sw)
    N = round(math.sqrt(L_out / core.A_L_mks))
    B_pk_pft = L_out * I_pft / (N * core.A_e_mks)
    # Everything considering peak flux density and potential problems with it
    extra_info = []
    if B_pk_pft > core.B_sat:
        core_fits = False
        # max possible number of turns without saturation
        B_pk_max = 0.8 * core.B_sat
        N = math.floor(core.A_e_mks * B_pk_max / (core.A_L_mks * I_pft))
        # max inductance without saturation
        L_out = core.A_L_mks * N ** 2
        # realistic current swing
        dI = (V_in_min - V_o) * D_max / (F_sw * L_out)
        ripple_real = dI / I_o
        if ripple_real < 2:
            # Continuous mode
            mode = f'continuous (r = {round(ripple_real, 2)})'
            C_out = dI / (8 * F_sw * dV_out)
            I_C_out_RMS = dI / math.sqrt(12)
            I_pk = I_o + dI / 2
        else:
            # Discontinuous mode
            # p.23 on http://www.ti.com/lit/an/slva057/slva057.pdf and
            # "Analog Circuit Design Volume 2: Immersion in the Black Art of
            # Analog Design" p.108
            mode = f'discontinuous (r = {round(ripple_real, 2)})'
            I_pk = dI
            C_out = I_o * (1 - I_o / dI) ** 2 / (F_sw * dV_out)
            I_C_out_RMS = math.sqrt(
                2 * I_o * V_o * (V_in_max - V_o) / (L_out * F_sw * V_in_min))

        dt_off_max = L_out * ripple * I_o / V_o  # T_off necessary to be
        # fine with this inductance
        F_sw_target = (1 - D_min) / dt_off_max
        D_min_target = 1 - dt_off_max * F_sw
        extra_info += [
            ('',
             warning_message(
                 f"We can't use this core to achieve the current ripple r = "
                 f"{ripple} we need. "
                 f"Please, use a larger core and/or a core with lower $A_L$. "
                 f"Worst case, we can still use this core for filtering "
                 f"by reducing the number of turns and compensating increased "
                 f"pulsations with a larger output capacitor "
                 f"(below you can see the numbers for that scenario). "
                 f"Alternatively, you can change the switching frequency "
                 f"to {format_F(F_sw_target)} or D_min to "
                 f"{round(D_min_target, 2)}, or do a bit of both").data)]
    else:
        core_fits = False
        mode = 'continuous'
        I_pk = I_o * (1 + ripple / 2)
        # Derivation of this is on figure 13.1 of "Switching power supplies
        # A-Z" by S. Maniktala
        C_out = dI / (8 * F_sw * dV_out)
        L_out = core.A_L_mks * N ** 2
        I_C_out_RMS = dI / math.sqrt(12)
    ESR_max = dV_out / dI
    P_C_out = I_C_out_RMS ** 2 * ESR_max
    B_pk = L_out * I_pft / (N * core.A_e_mks)
    extra_info += [
        ('Output voltage', format_V(V_o)),
        ('Output current', format_I(I_o)),
        (f'Output voltage ripple ({round(100 * C_ripple, 1)}%)',
         format_V(dV_out)),
        ('Output inductor', f'{format_L(L_out)}, {N} turns'),
        ('Inductor mode', mode),
        ('Inductor\'s peak flux density',
         format_flux_density(B_pk, core.B_sat)),
        ('Recommended minimum capacitor:', format_C(C_out)),
        ('Max ESR of this capacitor should be', format_R(ESR_max)),
        ('Output capacitor RMS current', format_I(I_C_out_RMS)),
        ('Such capacitor will dissipate', format_W(P_C_out)),
        (f'Desirable inductor', f'{format_L(L_out)}, {N} turns'),
        ('Inductor\'s peak current:', format_I(I_pk)),
    ]
    return (
        BuckOutputParams(L=L_out, N=N, C=C_out, ESR=ESR_max,
                         core_fits=core_fits),
        extra_info)


class NonDissipativeIsolatedClamp(NamedTuple):
    pass


# noinspection PyPep8Naming
def non_dissipative_isolated_clamp(
        V_in_min: float, V_in_max: float, I_pri_pk: float,
        L_leakage: float, T_clamp_on: float, T_clamp_off: float,
        V_D_clamp: float = 1.0,
        V_ripple: float = 0.1) -> Tuple[NonDissipativeIsolatedClamp,
                                        PrintableValues]:
    """
    Calculates a non-dissipating clamp suitable for forward, flyback
    and push-pull topologies. It utilizes two diodes, a capacitor
    and an inductor to absorb energy spike produced by the leakage inductance
    and dump it back to the input rail.

    :param V_in_min: Minimum input voltage
    :param V_in_max: Maximum input voltage
    :param I_pri_pk: Peak current flowing through the primary winding
        when the switch goes off
    :param L_leakage: leakage inductance
    :param T_clamp_on: for how long the clamp's transistor stays on
    :param T_clamp_off: for how long the clamp's transistor stays off
    :param V_D_clamp: voltage drop on the clamp's diodes
    :param V_ripple: voltage ripple on the clamp's capacitor
    """
    V_clamp_base = V_in_min
    V_clamp_max = (1 + V_ripple) * V_clamp_base
    # how much energy is stored in the leakage inductance when
    # the switch turns off
    E_leakage = L_leakage * I_pri_pk ** 2 / 2
    # the energy stored in the leakage inductance must be spent
    # by raising the energy
    # of the capacitor, thus raising its voltage from V_base to V_max.
    # In other words E_leak = C*V_max^2/2 - C*V_base^2/2,
    # from where we can find C
    C_clamp = nearest_e24(
        2 * E_leakage / (V_clamp_max ** 2 - V_clamp_base ** 2), 'higher')
    # We have a capacitor supposedly charged to V_clamp_max
    # Now we need to discharge it down to V_clamp_base into the inductor
    # in a time T_on, thus effectively dumping all absorbed energy
    # into the inductor.
    # However, because both V_in and T_on are floating, it's better to choose
    # an inductance large enough so that Cs will not be discharged
    # significantly. That will lead to a steady growth V_clamp_base,
    # until it reaches some equilibrium, but this will also help to keep
    # V_clamp_base above V_in, which is important since we don't want,
    # for example, the clamp to absorb the energy reflected from the other
    # half of the primary winding in push-pull topology. So we empirically
    # choose an inductance that normally will drop voltage across Cs
    # to V_clamp_base in T_clamp_on + T_clamp_off.
    # This will also guarantee that the inductor works in continuous mode and
    # T_on_clamp + T_off_clamp, which will simplify further calculation.
    # This LC circuit with an initially charged capacitor follows equation
    # V(t) = V_clamp_max * cos(t / sqrt(Ls * Cs)), so it V(t) = V_clamp_base,
    # we can derive Ls:
    LC_clamp = ((T_clamp_on + T_clamp_off)
                / math.acos(V_clamp_base / V_clamp_max)) ** 2
    L_clamp = LC_clamp / C_clamp
    # How much power we loose and (ideally) recover each cycle
    P_R_clamp = E_leakage / (T_clamp_on + T_clamp_off)

    # eventually min/max currents in Ls will settle down around specific values
    # allowing the capacitor to fully dump a fixed amount
    # (E_leakage, assuming 100% efficiency) of joules
    # into the inductor during transistor's T_clamp_on = D / F_sw and
    # then allowing the inductor to fully pump the same amount of energy
    # into the input rail  during T_clamp_off
    # (for push-pull this interval is equal to (2 - D) / F_sw).
    # For Ls (works in continuous mode) when it's pumping its energy back
    # to the input: V_in - 2 * V_diode = Ls * dI_Ls / T_clamp_off
    # From where we can find the current swing
    dI_L_clamp = (V_clamp_base - 2 * V_D_clamp) * T_clamp_off / L_clamp
    # would be the amplitute of current in Ls-Cs circuit
    # if allowed to oscillate freely with E_leakage energy in it
    I_L_clamp_amp = math.sqrt(E_leakage * 2 / L_clamp)
    # Roughly the center of the ramp of the current in Ls
    I_L_clamp = I_L_clamp_amp / 2

    V_sw_peak = (2 + V_ripple) * V_in_max
    # Duration of the leakage inductance surge
    Q_D_upper = 2 * E_leakage / V_clamp_base
    I_D_upper_pk = I_pri_pk + dI_L_clamp
    I_D_upper_rms = Q_D_upper / (T_clamp_on + T_clamp_off)

    components = NonDissipativeIsolatedClamp()
    return components, [
        ('Peak current through the upper diode $D_{SB}$',
         format_I(I_D_upper_pk)),
        ('RMS current through the upper diode $D_{SB}$',
         format_I(I_D_upper_rms)),
        ('Leakage inductance', format_L(L_leakage)),
        ('Max drain-source voltage:', format_V(V_sw_peak)),
        # ('Max voltage on $C_S$', format_V(V_clamp_max_new)),
        ('Min pre-charged voltage on $C_S$', format_V(V_clamp_base)),
        ('Minimum $L_S$', format_L(L_clamp)),
        ('Minimum $C_S$', format_C(C_clamp)),
        # ('RMS current through $L_S$ (and the lower diode $D_{SA})$',
        # format_I(I_L_clamp_rms)),
        ('Peak current through $C_S$', format_I(I_pri_pk)),
        ('Power recovered by one clamps is', format_W(P_R_clamp)),
        ('$L_s$ peak-to-peak ripple $\\Delta I$', format_I(dI_L_clamp)),
        ('$L_s$ RMS (center-ramp) current $I_{Ls}$', format_I(I_L_clamp)),
        ('$L_s$ peak current', format_I(I_L_clamp + dI_L_clamp / 2)),
        ('Max drain-source voltage/reverse voltage across each of the diodes:',
         format_V(V_sw_peak)),
    ]


class NonDissipativeCoupledClamp(NamedTuple):
    C: float
    I_L_rec_rms: float


# noinspection PyPep8Naming
def non_dissipative_coupled_clamp(
        V_in_min: float, V_in_max: float, I_pri_pk: float,
        L_leakage: float, T_clamp_on: float, T_clamp_off: float,
        V_D_clamp: float = 1.0,
        V_ripple: float = 0.1) -> Tuple[NonDissipativeCoupledClamp,
                                        PrintableValues]:
    """
    Calculates a non-dissipating inductor-less clamp for a push-pull converter.
    It utilizes (per side) one diode, a capacitor and an extra winding
    on the main transformer to absorb energy spike produced by the leakage
    inductance and dump it back to the input rail.

    :param V_in_min: Minimum input voltage
    :param V_in_max: Maximum input voltage
    :param I_pri_pk: Peak current flowing through the primary winding
        when the switch goes off
    :param L_leakage: leakage inductance
    :param T_clamp_on: for how long the clamp's transistor stays on
    :param T_clamp_off: for how long the clamp's transistor stays off
    :param V_D_clamp: voltage drop on the clamp's diodes
    :param V_ripple: voltage ripple on the clamp's capacitor
    """
    V_clamp_base = V_in_min
    V_clamp_max = (1 + V_ripple) * V_clamp_base
    # how much energy is stored in the leakage inductance when
    # the switch turns off
    E_leakage = L_leakage * I_pri_pk ** 2 / 2
    # How much charge we dump into the capacitor
    Q_dump = E_leakage / (V_clamp_base - V_D_clamp)
    I_D_rms = Q_dump / (T_clamp_on + T_clamp_off)
    # the energy stored in the leakage inductance must be spent
    # by raising the energy
    # of the capacitor, thus raising its voltage from V_base to V_max.
    # In other words E_leak = C*V_max^2/2 - C*V_base^2/2,
    # from where we can find C
    # Note: The article "Nondissipative Clamping Benefits DC-DC Converters"
    # uses the same formula, just written a bit differently.
    C_clamp = nearest_e24(
        2 * E_leakage / (V_clamp_max ** 2 - V_clamp_base ** 2), 'higher')
    D = T_clamp_on / (T_clamp_on + T_clamp_off)
    I_L_rec_rms = I_pri_pk * math.sqrt(D * (1 - D)) / 2
    F_res = 0.5 * math.pi * math.sqrt(L_leakage * C_clamp)
    components = NonDissipativeCoupledClamp(C=C_clamp, I_L_rec_rms=I_L_rec_rms)
    return components, [
        ('Minimum $C_S$', format_C(C_clamp)),
        ('RMS current of the recovery winding', format_I(I_L_rec_rms)),
        ('Resonant frequency (must be $<< F_{sw}$)', format_F(F_res)),
        ('Max drain-source voltage:', format_V(2 * V_in_max)),
        ('$D_S$ reverse voltage', format_V(2 * V_in_max)),
        ('$D_S$ peak current', format_V(I_pri_pk)),
        ('$D_S$ RMS current', format_I(I_D_rms)),
    ]
