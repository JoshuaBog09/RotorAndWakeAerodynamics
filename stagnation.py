import numpy as np
import matplotlib.pyplot as plt

import design

def compute_stagstat_pressure(U_inf: float, a: list[float], r_inner: list[float], r_outer: list[float], 
                              rho_inf: float=1.225, p_inf:float = 101325) -> None:
    """
    Generate two plots, the static pressure and stagnation (total) pressure as a function of spanwise 
    position on the rotor. The plot is created using the final data found after the running the bem
    procedure for all segmenents. Then for all segement the static and stagnation pressure are computed
    and displayed on two unique plots. The evaluation of the spanwise variation is performed for 4 unique
    locations. The upwind condition, the downwind condition, and just before and after the rotor.
    
    Parameters
    ----------

    U_inf : float
        Upstream wind velocity in `m/s` 

    a : list[float]
        List containing the final induction factor for each annuli segment in `-`

    r_inner : list[float]
        List containing the inner spanwise position (or inner radius) of the annuli of all segments in `m`

    rho_inf : list[float]
        Default : 1.225 kg/m^3
        The air density at upwind condition in `kg/m^3`
    
    p_inf : list[float]
        Default : 101325 Pa
        The air pressure at upwind condition in `pa`
    
    Returns
    -------

    None
        Plots are generated upon calling the function

    Raises
    ------

    asserts
        Checks to ensure the list of a, r_inner and r_outer are equal in size
    """
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
    A_upwind    = np.pi * (r_outer**2 - r_inner**2) * (1 - a)
    A_downwind  = np.pi * (r_outer**2 - r_inner**2) * (1 - a) / (1 - 2*a)
    A_rotor     = np.pi * (r_outer**2 - r_inner**2)

    # Compute stagnation pressures at 4 locations -->> Upwind, downwind, before and after rotor
    pstat_upwind = np.repeat(p_inf, len(a))
    pstat_rotor1 = p_inf + 0.5 * rho_inf * U_upwind**2 - 0.5 * rho_inf * U_rotor ** 2   # bernoulli

    F_rw = - (rho_inf * A_upwind * U_upwind**2 - rho_inf * A_downwind * U_downwind**2)  # Momentum
    pstat_rotor2 = pstat_rotor1 + (F_rw / A_rotor)

    pstat_downwind = pstat_rotor2 + 0.5 * rho_inf * U_rotor**2 - 0.5 * rho_inf * U_downwind ** 2    # bernoulli

    # Convert all to stagnation pressure
    pstag_upwind    = pstat_upwind + 0.5*rho_inf*U_upwind**2
    pstag_downwind  = pstat_downwind + 0.5*rho_inf*U_downwind**2
    pstag_rotor1    = pstat_rotor1 + 0.5*rho_inf*U_rotor**2
    pstag_rotor2    = pstat_rotor2 + 0.5*rho_inf*U_rotor**2
    
    ## Plotting
    colors = ['#d7191c','#fdae61','#abdda4','#2b83ba']

    fig_stat, ax_stat = plt.subplots(nrows=1, ncols=1)
    ax_stat.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstat_upwind, label="Upwind", color=colors[0], linestyle="dotted", marker='o',markevery=13)
    ax_stat.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstat_rotor1, label="Before rotor", color=colors[1], linestyle="solid")
    ax_stat.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstat_rotor2, label="After rotor", color=colors[2], linestyle="solid")
    ax_stat.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstat_downwind, label="Downwind", color=colors[3], linestyle="dotted", marker='d',markevery=11)
    ax_stat.set_xlabel("Normalized spanwise position [-]")
    ax_stat.set_ylabel("Static pressure [Pa]")
    ax_stat.legend()
    ax_stat.grid()

    fig_stag, ax_stag = plt.subplots(nrows=1, ncols=1)  
    ax_stag.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_upwind, label="Upwind", color=colors[0], linestyle="dotted", marker='o',markevery=13)
    ax_stag.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_rotor1, label="Before rotor", color=colors[1], linestyle="solid")
    ax_stag.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_rotor2, label="After rotor", color=colors[2], linestyle="solid")
    ax_stag.plot((r_inner + 0.5 * (r_outer - r_inner)) / design.R, pstag_downwind, label="Downwind", color=colors[3], linestyle="dotted", marker='d',markevery=11)
    ax_stag.set_xlabel("Normalized spanwise position [-]")
    ax_stag.set_ylabel("Stagnation pressure [Pa]")
    ax_stag.legend()
    ax_stag.grid()