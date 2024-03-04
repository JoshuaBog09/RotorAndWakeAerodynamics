import numpy as np


def force_azi_axi(a: float, U0: float, c: float, RHO: float, cl: float, cd: float):
    vp = (1 - a) * U0
    lift = 0.5 * RHO * vp ** 2 * c * cl
    drag = 0.5 * RHO * vp ** 2 * c * cd
    f_azi = lift * np.sin(twist) - drag * np.cos(twist)
    f_axi = lift * np.cos(twist) + drag * np.sin(twist)
    return f_azi, f_axi


def induction(f_azi: float, f_axi: float, BLADES: int, U0: float, RHO: float, R: float, r: float, dr: float, tsr: float,
              MU_ROOT: float) -> tuple:
    annuli_area = 2 * np.pi * r * dr
    ct = f_axi * BLADES * dr / (0.5 * RHO * U0 ** 2 * annuli_area)
    a = 0.5 * (1 - np.sqrt(1 - ct))
    a_prime = f_azi * BLADES / (2 * RHO * (2 * np.pi * r) * U0 ** 2 * (1 - a) * tsr * r / R)

    a = glauert(a, ct)
    a, a_prime = prandlt(a, a_prime, BLADES, r, R)

    return a, a_prime


def glauert(a: float, ct: float) -> float:
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    if ct >= CT2:
        a = 1 + (ct - CT1) / (4 * np.sqrt(CT1) - 4)
    return a


def prandlt(a: float, a_prime: float, BLADES: int, r: float, R: float) -> tuple:
    mu = r / R  # Normalised annuli position
    f_tip = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (1 - mu) / mu * np.sqrt(1 + tsr ** 2 * mu ** 2 / (1 - a) ** 2)))
    f_root = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (mu - MU_ROOT) / mu * np.sqrt(1 + tsr ** 2 * mu ** 2 / (1 - a) ** 2)))
    f_total = f_tip * f_root
    a = a / f_total
    a_prime = a_prime / f_total
    return a, a_prime

