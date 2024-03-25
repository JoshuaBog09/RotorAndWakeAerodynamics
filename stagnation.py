import numpy as np
import matplotlib.pyplot as plt

import design

def compute_stagnation_pressure(U_inf: float, a: list[float], r_inner: list[float], r_outer: list[float], rho_inf: float=1.225,
                                p_inf:float = 101325) -> None:

    # Check to ensure right size of input arrays
    assert len(a) == len(r_inner) == len(r_outer), f"The length of a, r_inner and r_outer should be equivalent"

    # Initialise arrays to numpy arrays
    a = np.array(a)
    r_inner = np.array(r_inner)
    r_outer = np.array(r_outer)

    # Compute Velocity conditions at 3 positions based on induction factor
    U_upwind    = np.repeat(U_inf, len(a))
    U_downwind  = U_inf * (1 - 2 * a)
    U_rotor     = U_inf * (1 - a)
    
    # Evaluate array of annuli at 3 positions based on induction factor based on continuity
    A_downwind  = np.pi * (r_outer**2 - r_inner**2) * (1 - a) / (1 - 2*a)
    A_rotor     = np.pi * (r_outer**2 - r_inner**2)
    A_upwind    = np.pi * (r_outer**2 - r_inner**2) * (1 - a)

    # Compute stagnation pressures at 4 locations -->> Upwind, downwind, before and after rotor
    pstag_upwind = np.repeat(p_inf, len(a))
    pstag_rotor1 = p_inf + 0.5 * rho_inf * U_upwind**2 - 0.5 * rho_inf * U_rotor ** 2   # bernoulli

    F_rw = - (rho_inf * A_upwind * U_upwind**2 - rho_inf * A_downwind * U_downwind**2)  # Momentum
    pstag_rotor2 = pstag_rotor1 + (F_rw / A_rotor)

    pstag_downwind = pstag_rotor2 + 0.5 * rho_inf * U_rotor**2 - 0.5 * rho_inf * U_downwind ** 2    # bernoulli

    ## Plotting
    colors = ['#d7191c','#fdae61','#abdda4','#2b83ba']

    plt.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_upwind, label="Upwind", color=colors[0], linestyle="dotted", marker='o',markevery=13)
    plt.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_rotor1, label="Before rotor", color=colors[1], linestyle="solid")
    plt.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_rotor2, label="After rotor", color=colors[2], linestyle="solid")
    plt.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_downwind, label="Downwind", color=colors[3], linestyle="dotted", marker='d',markevery=11)
    plt.xlabel("Normalized spanwise position [-]")
    plt.ylabel("Stagnation (static) pressure [Pa]")
    plt.legend()
    plt.grid()
    plt.show()