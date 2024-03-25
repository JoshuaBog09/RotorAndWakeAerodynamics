import numpy as np
import matplotlib.pyplot as plt

import math

def compute_stagnation_pressure(U_inf: float, a: np.ndarray, r_inner: np.ndarray, r_outer: np.ndarray):

    assert len(a) == len(r_inner) == len(r_outer), f"The length of a, r_inner and r_outer should be equivalent"

    A_upwind    = np.zeros(len(a)) 
    A_downwind  = np.zeros(len(a)) 
    A_rotor     = np.zeros(len(a))

    U_upwind    = np.zeros(len(a)) 
    U_downwind  = np.zeros(len(a)) 
    U_rotor     = np.zeros(len(a))

    for i in range(len(a)):

        U_upwind[i] = U_inf
        U_downwind[i] = U_inf * (1 - 2 * a[i])
        U_rotor[i] = U_inf * (1 - a[i])

        A_downwind[i] = math.pi * (r_outer[i]**2 - r_inner[i]**2) * (1 - a[i]) / (1 - 2*a[i])
        A_rotor[i] = math.pi * (r_outer[i]**2 - r_inner[i]**2)
        A_upwind[i] = math.pi * (r_outer[i]**2 - r_inner[i]**2) * (1 - a[i])

    print(U_upwind, U_downwind, U_rotor)
    print(sum(A_upwind), sum(A_downwind), sum(A_rotor))










