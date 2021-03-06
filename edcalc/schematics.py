import os.path

from typing import Callable
from lcapy import Circuit


def static_file_path(*path_parts: str) -> str:
    dir_name = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(dir_name, *path_parts)


def static_image_file(image_file_name: str) -> str:
    return static_file_path('static/images', image_file_name)


def static_circuit_image(image_name: str, label_nodes=False, help_lines=False,
                         scale=0.3):
    def decorator(circuit_generator: Callable[[], Circuit]) -> Callable[
            [bool], str]:

        def image_generator(cache: bool = True) -> str:
            image_path = static_image_file(image_name)
            if not os.path.exists(image_path) or not cache:
                cct = circuit_generator()
                cct.draw(image_path,
                         style='american',
                         draw_nodes='connections',
                         label_nodes=label_nodes,
                         help_lines=help_lines,
                         scale=scale)
            return image_path

        return image_generator
    return decorator


@static_circuit_image('half-bridge.png')
def half_bridge_topology_image():
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
    return cct


@static_circuit_image('input-emi-filter.png')
def emi_filter_image():
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
    return cct


@static_circuit_image('push-pull-topology.png')
def push_pull_topology_image():
    cct = Circuit()
    # the input (switches-capacitors-transformer)
    cct.add('''
    P1 p1a p1b; down=3,v=V_{in}
    W p1a c1a; right
    W p1b c1b; right
    C1 c1a c1b; down
    W c1a tf1t5; right=3
    TF1 t1 t2 t3 t4 tapcore _t5 _t6; right, size=1.5
    W tf1t5 tf1t52; down
    W tf1t52 _t5; right
    M2 2d 2g 2s nmos; right
    W 2d t4; up
    M1 1d 1g 1s nmos; right
    W 1d tf1t3a; up=2.5
    W tf1t3a t3; right=2.5
    ''')
    cct.add('''
    W 1s 2s; right=2
    Rs 1s rsb; down
    W rsb c1b; left=3
    W 1s isense; left,l=$I_{sense}$
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
    return cct


@static_circuit_image('push-pull-rcd-snubber.png')
def push_pull_rcd_snubber_image():
    cct = Circuit()
    cct.add('''
    P1 p1a p1b; down=3,v=V_{in}
    W p1a c1a; right
    W p1b c1b; right
    C1 c1a c1b; down
    W c1a tf1t5; right=6
    TF1 t1 t2 t3 t4 tapcore _t5 _t6; right, size=1.5
    W tf1t5 tf1t52; down
    W tf1t52 _t5; right
    M2 2d 2g 2s nmos; right
    Lleakage2 2d t4; up,color=red
    M1 1d 1g 1s nmos; right
    Lleakage1 1d lleak1b; up,color=red
    W lleak1b tf1t3a; up=1.5
    W tf1t3a t3; right=5
    ''')
    # Second snubber
    cct.add('''
    Csnubber2 2d cc2b; right,color=blue
    Dsnubber2 cc2b dc2b; down,color=blue
    W cc2b rc2a; right,color=blue
    Rsnubber2 rc2a rc2b; down,color=blue
    W rc2b dc2b; left,color=blue
    ''')
    # First snubber
    cct.add('''
    Csnubber1 1d cc1b; right,color=blue
    Dsnubber1 cc1b dc1b; down,color=blue
    W cc1b rc1a; right,color=blue
    Rsnubber1 rc1a rc1b; down,color=blue
    W rc1b dc1b; left,color=blue
    ''')
    cct.add('''
    W 2s 1s; left=4.5
    Rs 1s rsb; down
    W c1b rsb; right
    W dc1b rsb; left
    W dc2b rc1b; left=3.5
    ''')
    return cct


@static_circuit_image('push-pull-lossless-clamp.png')
def push_pull_lossless_clamp_image():
    cct = Circuit()
    cct.add('''
    P1 p1a p1b; down=3,v=V_{in}
    W p1a c1a; right
    W p1b c1b; right
    C1 c1a c1b; down
    W c1a clampret1; right
    W clampret1 clampret2; right=4
    W clampret2 tf1t5; right
    TF1 t1 t2 t3 t4 tapcore _t5 _t6; right
    W tf1t5 tf1t52; down
    W tf1t52 _t5; right
    M2 2d 2g 2s nmos; right
    W 2d 2d2; up,fixed
    Lleakage2 2d2 t4; up,color=red
    M1 1d 1g 1s nmos; right
    W 1d 1d2; up,fixed
    Lleakage1 1d2 lleak1b; up,color=red
    W lleak1b tf1t3a; up
    W tf1t3a t3; right=4
    ''')
    # First clamp
    cct.add('''
    CS1 1d2 cc1b; left=2,color=blue
    DS1_A ds1a ls1a; up, color=blue
    LS1 ls1a cc1b; up, color=blue
    DS1_B cc1b clampret1; up, color=blue 
    W ds1a ds1a1; down, color=blue
    ''')
    # Second clamp
    cct.add('''
    CS2 2d2 cc2b; left=2,color=blue
    DS2_A ds2a ls2a; up, color=blue
    LS2 ls2a cc2b; up, color=blue
    DS2_B cc2b clampret2; up, color=blue 
    W ds2a ds2a1; down, color=blue
    ''')
    cct.add('''
    Rs 1s rsb; down
    W ds1a1 rsb; right=2
    W ds2a1 rsb; left=2
    W 1s 2s; right=2
    W c1b ds1a1; right
    ''')
    return cct


@static_circuit_image('push-pull-coupled-lossless-clamp.png', scale=0.25)
def push_pull_coupled_lossless_clamp_image():
    cct = Circuit()
    cct.add('''
    P1 p1a p1b; down=10,v=V_{in}
    W p1a c1a; right
    W p1b c1b; right
    C1 c1a c1b; down
    W c1a clampret1; right=1
    W clampret1 clampret2; right=4
    W clampret2 t1c1; right=3
    W t1c1 t1c2; down=2
    W t1c2 t1c3; right
    LW1 l1a t1c3; down
    LW2 t1c3 l2b; down
    W l1a lleak21; left=2
    W lleak21 lleak2a; down=3
    Lleak2 lleak2a lleak2b; down,color=red
    CS2 lleak2b cs2b; left=2,color=blue,v=$V_{CS}$
    W cs2b ds2a; up=2,color=blue
    DS2 ds2a ds2b; up=2.5,color=blue
    W ds2b clampret2; up,color=blue
    W ds2a lw3a1; right=3,color=blue
    W lw3a1 lw3a2; down=1,color=blue
    LW3 lw3a2 lw3b; down,color=blue
    LW4 lw3b lw4b; down,color=blue,fixed
    W l2b lleak1a1; left=6
    W lleak1a1 lleak1a; down=1.5
    Lleak1 lleak1a lleak1b; down,color=red
    CS1 lleak1b cs1b; left=2,color=blue,v=$V_{CS}$
    W cs1b ds1a; up=2,color=blue
    DS1 ds1a ds1b; up=2.5,color=blue
    W ds1b clampret1; up,color=blue
    W lleak1b m1d; down=0.5
    W lleak2b m2d; down=0.5
    M1 m1d m1g m1s nmos; right,size=1
    M2 m2d m2g m2s nmos; right,size=1
    W cs1b lw4b1; down=2.5,color=blue
    W lw4b1 lw4b2; right=8,color=blue
    W lw4b2 lw4b; up,color=blue
    W m1s m2s; right=4
    RS m1s rsb; down=2
    W c1b rsb; right
    W lw3b lw3b1; left=1,color=blue
    W lw3b1 lw3b2; down=4,color=blue
    W lw3b2 rsb; left=5,color=blue
    A1 lw3a2; l={\ \ \ \ \ •}
    A2 lw3b;  l={\ \ \ \ \ •}
    A3 l1a;   l={\ \ \ \ \ •}
    A4 t1c3;  l={\ \ \ \ \ •}
    ''')
    return cct


@static_circuit_image('forward-clamped-topology.png')
def forward_topology_image():
    cct = Circuit()
    # the input (switches-capacitors-transformer)
    cct.add('''
    P1 p1a p1b; down=5.5,v=V_{in}
    W p1a c1a; right
    W p1b c1b; right
    CIN c1a c1b; down
    W c1a dsa; right=2
    DS dsb dsa; up=2
    W dsb l1a; down
    LS l1a l1b; down
    W l1b l2b2; down=1.5
    W dsa l2a1; right=2
    W l2a1 l2b; down=0.5
    TF1 l3a l3b l2b l2a tapcore _t5 _t6; right, size=1.5
    CS l2a dsb; left,color=gray,v=$V_{Cs}$
    A1 l1b; l={\ \ \ \ \ •}
    W l2a m1d; down=0.5
    M1 m1d m1g m1s nmos; right,size=1
    W m1s rsa; down=0.5
    RSENSE rsa rsb; down
    W rsa rsa1; right,l=$I_{sense}$
    W c1b l2b2; right
    W l2b2 rsb; right
    W l3b l3b1; down
    D1 l3a d1b; right=1.5,kind=schottky
    W l3b1 d2a; right=1.5
    D2 d2a d1b; up,kind=schottky
    W d1b loa; right
    LO loa lob; right=1.5
    W lob lob2; right
    CO lob2 cob; down
    W lob voa; right
    W cob vob; right
    P2 voa vob; down,v=V_O
    DX loa dxb; down,color=gray
    CX dxb cxb; down,color=gray
    RX dxb rxb; right,color=gray
    W rxb lob; up,color=gray
    W d2a cxb; right
    W cxb cob; right
    ''')
    return cct


@static_circuit_image('buck-rectifier-clamp.png')
def buck_rectifier_clamp_image():
    cct = Circuit()
    cct.add('''
    P1 d1b d2a; l=from the rectifier, v=$V_{IN}$
    W d1b loa; right
    LO loa lob; right=1.5
    W lob lob2; right
    CO lob2 cob; down
    W lob voa; right
    W cob vob; right
    P2 voa vob; down,v=V_O
    DX loa dxb; down,color=blue
    CX dxb cxb; down,color=blue
    RX dxb rxb; right,color=blue
    W rxb lob; up,color=blue
    W d2a cxb; right
    W cxb cob; right
    ''')
    return cct



@static_circuit_image('rectifier-rc-snubber.png')
def rectifier_rc_snubber_image():
    cct = Circuit()
    cct.add('''
    TF1 l3a l3b l2b l2a tapcore _t5 _t6; right, size=1, l=T
    W _t6 t6g; right=0.1, sground
    D1 l3a d1b; right=2
    D2 l3b d2b; right=2
    W l3a cs1a; up,color=blue
    W l3b cs2a; down, color=blue
    CS1 cs1a cs1b; right, l=$C_S$, color=blue
    CS2 cs2a cs2b; right, l=$C_S$, color=blue
    RS1 cs1b rs1b; right, l=$R_S$, color=blue
    RS2 cs2b rs2b; right, l=$R_S$, color=blue
    W d1b rs1b; up, color=blue
    W d2b rs2b; down, color=blue
    W d2b d1b; up
    LO d1b lob; right 
    CO lob cob; down
    W cob cob1; down=0.1, sground 
    W lob lob2; right, l=Vout
    ''')
    return cct
