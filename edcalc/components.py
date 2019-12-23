"""
A small library of components that simplifies entering the design parameters
"""
import pathlib
import os.path

from .component_types import Core, OptoCoupler, ShuntReference
from .format import format_markdown_table


# A bunch of magnetic cores (for transformer or inductors) to avoid copying
# their parameters between the designs.
# Each core can be described using only few basic datasheet parameters,
# the rest will be calculated. This allows having just a caliper and an
# LCR meter quickly add even second-hand unknown cores to the list.
cores = {
    'PQ2620-PC40': Core(
        A_L=5800, A_e=117.67, l_e=46.32, V_e=5450.4, B_sat=0.3,
        note='PQ26/20 core, material: PC40'),
    'R12.5x7.5x5-N87': Core(
        mu=2200, A_L=1120, l_e=30.09, A_e=12.23, V_e=368, B_sat=0.3,
        note='Ring 12.5x7.5x5, material: N87 (Epcos)'),
    'R33x17.8x11.1-MM52': Core(
        A_e=80.5, V_e=6410, mu=75, A_L=95, B_sat=1.6,
        note='Green-blue ring 33x17.8x11.1, material: Micrometals -52'),
    'R22x14x6-N87': Core(
        A_L=1340, mu=2200, l_e=54.15, A_e=26.17, V_e=1417, B_sat=0.3,
        note='Ring R22x14x6, material: N87 (Epcos)'),
    'EE40-PC40': Core(
        A_L=4000, l_e=78.943, A_e=145.95, V_e=11521.7, W_a=94, B_sat=0.3,
        note='EE40 core, material: PC40 or equivalent, unknown brand'),
    'PC21/13-USSR': Core(
        A_e=57.7, V_e=1862, l_e=32.27, A_L=100,
        note='soviet pot core 21/13 from some old equipment')
}


optocouplers = {
    'PC817': OptoCoupler(V_f=1.0, V_sat=0.2, CTR_min=0.5, C_out=4.5e-9),
}


references = {
    'TL431': ShuntReference(I_bias=1e-3, V_min=2.5, V_ref=2.5),
    'BZX55C6V2': ShuntReference(I_bias=1e-3, V_min=6.2),
}


# Relative path to this very file
LIBRARY_PATH = os.path.join(*pathlib.Path(__file__).parts[-2:])


def table_of_cores():
    """
    Returns a Markdown table containing all cores from this library.
    """
    return (
        '### Magnetic core library '
        f'(from [{LIBRARY_PATH}]({LIBRARY_PATH}))\n' +
        'Variable `cores[\'<some key\']` can be used for adding the cores ' +
        'below to the design parameters.\n\n' +
        format_markdown_table(
            ['key', 'permeability $ \\mu $', '$ A_e $', '$ l_e $', '$ V_e $',
             '$ B_{sat} $', 'note'],
            [
                [key, str(round(c.mu, 2)), str(round(c.A_e, 2)),
                 str(round(c.l_e, 2)), str(round(c.V_e, 2)),
                 str(c.B_sat), str(c.note)]
                for key, c in cores.items()
            ]
        ))
