import numpy as np
import math

def bem_procedure(U0: float, segment_c: float, r:float, R:float, tsr:float, segment_twist: float,
                   blade_pitch: float, RHO: float, polar_sheet: np.ndarray, BLADES: int, MU_ROOT: float,
                   dr: float):

    a, a_prime = [0.3], [0]
    
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
        beta = segment_twist + blade_pitch
        alpha = Phi - beta
        
        # Interpolate polar data to find cl and cd
        cl = np.interp(alpha, polar_sheet[0,:], polar_sheet[1,:])
        cd = np.interp(alpha, polar_sheet[0,:], polar_sheet[2,:])

        f_azi, f_axi = force_azi_axi(V_p, segment_c, Phi, RHO, cl, cd) 
        a_new, a_prime_new = induction(f_azi, f_axi, BLADES, U0, RHO, R, r, dr, tsr, MU_ROOT)

        a.append(0.25 * a_new + 0.75 * a[-1])
        a_prime.append(0.25 * a_prime_new + 0.75 * a_prime[-1])

        if abs(a[-1] - a[-2]) <=  0.0005 and abs(a_prime[-1] - a_prime[-2]) <= 0.0005:
            iterating = False
    
    return None


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
    a = 0.5 * (1 - np.sqrt(1 - ct))
    a_prime = f_azi * BLADES / (2 * RHO * (2 * np.pi * r) * U0 ** 2 * (1 - a) * tsr * r / R)

    a = glauert(a, ct)
    a, a_prime = prandlt(a, a_prime, BLADES, r, R, tsr, MU_ROOT)

    return a, a_prime


def glauert(a: float, ct: float) -> float:
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    if ct >= CT2:
        a = 1 + (ct - CT1) / (4 * np.sqrt(CT1) - 4)
    return a


def prandlt(a: float, a_prime: float, BLADES: int, r: float, R: float, tsr: float, MU_ROOT: float) -> tuple:
    mu = r / R  # Normalised annuli position
    f_tip = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (1 - mu) / mu * np.sqrt(1 + tsr ** 2 * mu ** 2 / (1 - a) ** 2)))
    f_root = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (mu - MU_ROOT) / mu * np.sqrt(1 + tsr ** 2 * mu ** 2 / (1 - a) ** 2)))
    f_total = f_tip * f_root
    a = a / f_total
    a_prime = a_prime / f_total
    return a, a_prime

