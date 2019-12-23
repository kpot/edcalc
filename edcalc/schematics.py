import os.path

from lcapy import Circuit


def static_file_path(*path_parts: str) -> str:
    dir_name = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(dir_name, *path_parts)


def static_image_file(image_file_name: str) -> str:
    return static_file_path('static/images', image_file_name)


def half_bridge_topology_image(cache: bool = True):
    cct = Circuit()
    # the input (switches-capacitors-transformer)
    cct.add('''
    P1 p1a p1b; down,v=V_{in}
    C1 c1a c1b; down=3
    C2 c1b c2b; down=3
    R1 r1a r1b; down=3,l=Rb
    R2 r1b r2b; down=3,l=Rb
    W r1b c1b; right
    W p1a r1a; right
    W p1b r2b; right
    W r1a c1a; right
    W r2b c2b; right
    W c1a 1d; right=2
    M1 1d 1g 1s nmos; right
    W 1s 2d1; down
    W 2d1 2d; down
    M2 2d 2g 2s nmos; right
    W c2b 2s; right=2
    TF1 t1 t2 t3 t4 tapcore _t5 _t6; right, size=1.5
    W 2d1 t3; right
    Cb c1b t4; right=2
    ''')
    # First body diode
    cct.add('''
    W 1d dm1b; right=0.5
    W 1s dm1a; right=0.5
    DM1 dm1a dm1b; up
    ''')
    # Second body diode
    cct.add('''
    W 2d dm2b; right=0.5
    W 2s dm2a; right=0.5
    DM2 dm2a dm2b; up
    ''')
    # The output (transformer-rectifier-inductor-capacitor)
    cct.add('''
    W _t6 t6; right=2
    W t6 go1; down
    W go1 go2; right 
    D1 t1 d1b; right
    D2 t2 d2b; right
    W d2b d1b; up
    Lo d1b coa; right
    Co coa go2; down
    P2 p2a p2b; down, v=V_{out}
    W coa p2a; right
    W go2 p2b; right
    ''')
    image_path = static_image_file('half-bridge.png')
    if not os.path.exists(image_path) or not cache:
        cct.draw(image_path,
                 style='american',
                 draw_nodes='connections',
                 label_nodes=False,
                 scale=0.3)
    return image_path


def emi_filter_image(cache: bool = True):
    cct = Circuit()
    # the input (switches-capacitors-transformer)
    cct.add('''
    P1 vin_1 vin_2; down,l=V_{AC}
    W vin_1 cpi_up; rignt
    W vin_2 cpi_down; right
    Cpi cpi_up cpi_down; down
    W cpi_up t1_1up; right
    W cpi_down t1_2down; right
    TF1 t1_1 t1_2 t1_3 t1_4 core; right, rotate=90
    W t1_1up t1_1; down
    W t1_2down t1_3; up
    W t1_2 t1_2up; up
    W t1_4 t1_4down; down
    W t1_2up cy_up; right
    W t1_4down cy_down; right
    CY1 cy_up cy_j; down=1.5
    CY2 cy_j cy_down; down=1.5
    W cy_j cy_earth; right=0.8
    W cy_earth cy_earth2; down=2
    W cy_earth2 cy_earth3; left=4.8,l=earth
    W cy_up ldm1_1; right
    W cy_down ldm2_1; right
    LDM1 ldm1_1 ldm1_2; right,l=1/2L_{DM}
    LDM2 ldm2_1 ldm2_2; right,l=1/2L_{DM}
    CX ldm1_2 ldm2_2; down
    W ldm1_2 damp_up; right
    W ldm2_2 damp_down; right
    Cd damp_up damp_j; down, dashed
    Rd damp_j damp_down; down, dashed
    W damp_up vout_1; right
    W damp_down vout_2; right
    P vout_1 vout_2; l=Vout
    ''')
    image_path = static_image_file('input-emi-filter.png')
    if not os.path.exists(image_path) or not cache:
        cct.draw(image_path,
                 style='american',
                 draw_nodes='connections',
                 label_nodes=False,
                 scale=0.3)
    return image_path

