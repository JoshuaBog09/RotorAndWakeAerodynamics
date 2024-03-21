import numpy as np

def bem_procedure(U0: float, segment_c: float, r: float, R: float, tsr: float, segment_twist: float,
                   blade_pitch: float, RHO: float, polar_sheet: np.ndarray, BLADES: int, MU_ROOT: float,
                   dr: float, tolerance: float) -> tuple[float, float, float, float, float, float]:
    """
    Main `BEM procedure`, which can be called on a given turbine design. The bem model utelises a singular anuli 
    segment as input for its computational procedure. And should therfore be called for each anuli (segment) 
    of the given turbine, to obtain the full spectrum of results. The bem model will iterate over a range of 
    induction factors, to balance out the loads, untill the specified tolarenaces are met. Upon conversion of the
    program will return the last found induction factors (axial and azimuthal) which lead to satisfaction of the 
    specified covergence tolerance requirement.

    The iteration methode which is utelised is based on the previous and current value in evaluation, and can 
    be desribed as follows: `xᵢ = 0.25 * xᵢ + 0.75 * xᵢ₋₁` where x is the parameter, which is altered to obtain 
    convergence.

    Parameters
    ----------

    U0 : float
        Upstream wind velocity in `m/s` 

    segment_c : float
        Chord length of the airfoil segment in the annuli being analysed expressed in `m`

    r : float 
        Current blade length wise position of the evaluated segment expressed in `m`

    R : float
        Blade length of the turbine expressed in `m` 

    tsr : float
        Tip speed ratio in `-`

    segment_twist : float
        Twist angle of the segment expressed in `deg`
    
    blade_pitch : float
        Pitch angle expressed in `deg`

    RHO : float
        Air density expressed in `kg/m^3`

    polar_sheet : np.ndarray
        Numpy array containing cl and cd as a function of the angle of attack alpha

    BLADES : int
        The amount of blades expressed as an integer `-`

    MU_ROOT : float
        Normalised start location of the airfoil of the wind turbine expressed in `-` 
    
    dr : float
        Length of the blade segment expressed in `m`

    tolerance : float
        Exit criteria discribing the relation between the current and previous iteration value in `-`.  

    Returns
    -------

    tuple
        The axial induction factor of the segment at convergence : float
        The azimuthal induction factor of the segment at convergence : float
        The inflow angle of the segment : float
        The angle of attack of the segment : float
        The azimuthal force (f_tan) of the segment : float
        The axial force (f_norm) of the segment : float

    Raises
    ------

    StopIteration 
        Raised when either of the induction factors is out of bounds 
    """
    a, a_prime = [0.3], [0.]    # Initial inductions
    # phi_list = []
    # beta_list = []
    iterating = True
    iteration = 0

    while iterating:
        
        iteration += 1
        #print(f"{iteration}: a = {a[-1]}, a'= {a_prime[-1]}")
        
        V_axial = U0 * (1 - a[-1])  # [m/s]
        omega = tsr * U0 / R    # [rad/s]
        V_tan   = omega * r * (1 + a_prime[-1]) # [m/s]

        V_p = np.sqrt(V_axial**2 + V_tan**2)    # [m/s]
        Phi = np.arctan2(V_axial, V_tan) # inflow angle [rad]
        beta = np.radians(segment_twist) + np.radians(blade_pitch) # [rad]
        alpha = Phi - beta # [rad]
        
        # Interpolate polar data to find cl and cd
        cl = float(np.interp(np.degrees(alpha), polar_sheet[0,:], polar_sheet[1,:]))
        cd = float(np.interp(np.degrees(alpha), polar_sheet[0,:], polar_sheet[2,:]))

        f_azi, f_axi = force_azi_axi(V_p, segment_c, Phi, RHO, cl, cd) 
        a_new, correction_prandtl = induction(f_azi, f_axi, BLADES, U0, RHO, R, r, dr, tsr, MU_ROOT)

        a.append(0.25 * a_new + 0.75 * a[-1])

        a_prime_new = f_azi * BLADES / (2 * RHO * (2 * np.pi * r) * U0 ** 2 * (1 - a[-1]) * tsr * r / R)
        a_prime_new /= correction_prandtl

        a_prime.append(0.25 * a_prime_new + 0.75 * a_prime[-1])

        if abs(a[-1] - a[-2]) <= tolerance and abs(a_prime[-1] - a_prime[-2]) <= tolerance:
            iterating = False

        if not 0 < a[-1] < 1 or not 0 < a_prime[-1] < 1:
            raise StopIteration(f"Itreration paramter a or a' out of bounds [0, 1]. \
                                Value at failure a={a[-1]:.5f}, a'={a_prime[-1]:.5f}")
    
    return a[-1], a_prime[-1], Phi, alpha, f_azi, f_axi


def force_azi_axi(V_p: float, segment_c: float, Phi:float, RHO: float,
                   cl: float, cd: float) -> tuple[float, float]:
    lift = 0.5 * RHO * segment_c * cl * V_p ** 2
    drag = 0.5 * RHO * segment_c * cd * V_p ** 2
    f_azi = lift * np.sin(Phi) - drag * np.cos(Phi)     # f_tan
    f_axi = lift * np.cos(Phi) + drag * np.sin(Phi)     # f_norm
    return f_azi, f_axi


def induction(f_azi: float, f_axi: float, BLADES: int, U0: float, RHO: float, R: float, r: float, dr: float, 
              tsr: float, MU_ROOT: float) -> tuple[float, float]:
    annuli_area = 2 * np.pi * r * dr
    ct = f_axi * BLADES * dr / (0.5 * RHO * U0 ** 2 * annuli_area)

    # Review order of corrections
    a = glauert(ct)
    correction_prandtl = prandlt(a, BLADES, r, R, tsr, MU_ROOT)
    a /= correction_prandtl

    # a_prime = f_azi * BLADES / (2 * RHO * (2 * np.pi * r) * U0 ** 2 * (1 - a) * tsr * r / R)

    # a_prime /= correction_prandtl

    return a, correction_prandtl


def glauert(ct: float) -> float:
    """
    Perform glauert (heavily loaded streamtubes) correction based on the thrust coefficient. 
    
    Parameters
    ----------

    ct : float
        Thrust coefficient expressed in `-`
    
    Returns
    -------

    float
        Corrected axial induction (a) factor `-`
    
    Raises
    ------

    None
    
    """
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
    """
    Provided the prandtl (correction for finite number of blades) corrective term based on the position on the blade. 
    Correction can be applied by deviding both the axial and azimuthal induction factor utelising the returned corective 
    factor.
    
    Parameters
    ----------

    a : float
        Axial induction factor expressed in `-`

    BLADES : int
        The amount of blades expressed as an integer `-`

    r : float 
        Current blade length wise position of the evaluated segment expressed in `m`

    R : float
        Blade length of the turbine expressed in `m`
    
    tsr : float
        Tip speed ratio in `-`

    MU_ROOT : float
        Normalised start location of the airfoil of the wind turbine expressed in `-`

    Returns
    -------

    float
        Corrected factor (f_cor) for both the axial and azimuthal induction factor `-`
    
    Raises
    ------

    None
    
    """
    mu = r / R  # Normalised annuli position
    f_tip = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (1 - mu) / mu * np.sqrt(1 + (tsr ** 2 * mu ** 2) / (1 - a) ** 2)))
    f_root = 2 / np.pi * np.arccos(
        np.exp(-0.5 * BLADES * (mu - MU_ROOT) / mu * np.sqrt(1 + (tsr ** 2 * mu ** 2) / (1 - a) ** 2)))
    f_total = f_tip * f_root
    return f_total

def prandtl_carlos(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction):
    temp1 = -NBlades / 2 * (tipradius_R - r_R) / r_R * np.sqrt(1 + (TSR * r_R)**2 / (1 - axial_induction)**2)
    ftip = 2/np.pi * np.arccos(np.exp(temp1))
    for i in range(len(ftip)):
        if (np.isnan(ftip[i])):
            ftip[i]=0
    temp1 = NBlades / 2 * (rootradius_R - r_R) / r_R * np.sqrt(
        1 + (TSR * r_R)**2 / (1 - axial_induction)**2)
    froot = 2/np.pi * np.arccos(np.exp(temp1))
    for i in range(len(froot)):
        if (np.isnan(froot[i])):
            froot[i]=0
    ftotal = froot*ftip
    return ftotal

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import design

    r = np.linspace(design.start, design.end, 10000, endpoint=True) * design.R

    f_cor = prandlt(0.4, design.BLADES, r, design.R, design.TSR[0], design.start)
    f_cor_carlos = prandtl_carlos(r/design.R, design.start, design.end, design.TSR[0], design.BLADES, 0.4)
    fig_prandtl, ax_prandtl = plt.subplots(nrows=1, ncols=1)
    ax_prandtl.plot(r / design.R, f_cor), ax_prandtl.set_title("Prandtl's correction")
    ax_prandtl.plot(r / design.R, f_cor_carlos)
    ax_prandtl.set_ylabel("Prandtl's correction factor"), ax_prandtl.set_xlabel('r/R')
    ax_prandtl.grid()
    # plt.show()

    cts = np.linspace(-4, 4, 10000, endpoint=True)
    a = []
    for ct in cts:
        a.append(glauert(ct))

    fig_glauert, ax_glauert = plt.subplots(nrows=1, ncols=1)
    ax_glauert.plot(a, cts), plt.title(r"$C_T(a)$ including Glauert's correction")
    ax_glauert.set_ylabel(r'$C_T$'), plt.xlabel('a')
    ax_glauert.set_xlim([0,1]), ax_glauert.set_xticks(np.arange(0, 1.1, 0.1))
    ax_glauert.set_ylim([0,2])
    ax_glauert.grid()
    plt.show()