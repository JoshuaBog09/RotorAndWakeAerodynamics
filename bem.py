import numpy as np
import math

def bem_procedure(U0: float, segment_c: float, r:float, R:float, tsr:float, segment_twist: float,
                   blade_pitch: float, RHO: float, polar_sheet: np.ndarray, BLADES: int, MU_ROOT: float,
                   dr: float):

    a, a_prime = [0.3], [0]
    phi_list = []
    beta_list = []
    iterating = True
    iteration = 0

    while iterating:
        
        iteration += 1
        print(f"{iteration}: a = {a[-1]}, a'= {a_prime[-1]}")
        
        V_axial = U0 * (1 - a[-1])
        omega = tsr * U0 / R
        V_tan   = omega * r * (1 + a_prime[-1])

        V_p = np.sqrt(V_axial**2 + V_tan**2)
        Phi = np.arctan2(V_axial, V_tan) # inflow angle [rad]
        beta = np.radians(segment_twist) + np.radians(blade_pitch)
        alpha = Phi - beta
        
        # Interpolate polar data to find cl and cd
        cl = np.interp(np.degrees(alpha), polar_sheet[0,:], polar_sheet[1,:])
        cd = np.interp(np.degrees(alpha), polar_sheet[0,:], polar_sheet[2,:])

        f_azi, f_axi = force_azi_axi(V_p, segment_c, Phi, RHO, cl, cd) 
        a_new, a_prime_new = induction(f_azi, f_axi, BLADES, U0, RHO, R, r, dr, tsr, MU_ROOT)

        a.append(0.25 * a_new + 0.75 * a[-1])
        a_prime.append(0.25 * a_prime_new + 0.75 * a_prime[-1])

        if abs(a[-1] - a[-2]) <=  0.0005 and abs(a_prime[-1] - a_prime[-2]) <= 0.0005:
            iterating = False
    
    return a[-1], a_prime[-1]


def force_azi_axi(V_p: float, segment_c: float, Phi:float, RHO: float,
                   cl: float, cd: float):
    lift = 0.5 * RHO * V_p ** 2 * segment_c * cl
    drag = 0.5 * RHO * V_p ** 2 * segment_c * cd
    f_azi = lift * np.sin(Phi) - drag * np.cos(Phi)
    f_axi = lift * np.cos(Phi) + drag * np.sin(Phi)
    return f_azi, f_axi


def induction(f_azi: float, f_axi: float, BLADES: int, U0: float, RHO: float, R: float, r: float, dr: float, 
              tsr: float, MU_ROOT: float) -> tuple:
    annuli_area = 2 * np.pi * r * dr
    ct = f_axi * BLADES * dr / (0.5 * RHO * U0 ** 2 * annuli_area)

    # Review order of corrections
    a = glauert(ct)
    a_prime = f_azi * BLADES / (2 * RHO * (2 * np.pi * r) * U0 ** 2 * (1 - a) * tsr * r / R)

    correction_prandtl = prandlt(a, BLADES, r, R, tsr, MU_ROOT)

    a /= correction_prandtl
    a_prime /= correction_prandtl

    return a, a_prime


def glauert(ct: float) -> float:
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    
    if ct < CT2:
        a = 0.5 * (1 - np.sqrt(1 - ct))
    elif ct >= CT2:
        a = 1 + (ct - CT1) / (4 * np.sqrt(CT1) - 4)
    else:
       raise BaseException(f"Unexptected behaviour in Glauert correction")
    
    return a


def prandlt(a: float, BLADES: int, r: float, R: float, tsr: float, MU_ROOT: float) -> float:
    mu = r / R  # Normalised annuli position
    f_tip = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (1 - mu) / mu * np.sqrt(1 + (tsr ** 2 * mu ** 2) / (1 - a) ** 2)))
    f_root = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (mu - MU_ROOT) / mu * np.sqrt(1 + (tsr ** 2 * mu ** 2) / (1 - a) ** 2)))
    f_total = f_tip * f_root
    return f_total

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import design

    r = np.linspace(design.start, design.end, 1000, endpoint=True) * design.R

    f_cor = prandlt(0.8, design.BLADES, r, design.R, design.TSR[0], design.start)

    plt.plot(r / design.R, f_cor)
    plt.grid()
    plt.show()

    cts = np.linspace(-4, 4, 1000, endpoint=True)
    a = []
    for ct in cts:
        a.append(glauert(ct))

    plt.plot(a, cts)
    plt.xlim([0,1])
    plt.ylim([0,3])
    plt.grid()
    plt.show()