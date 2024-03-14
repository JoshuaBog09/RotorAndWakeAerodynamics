import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import math
import statistics

import design
import bem


file_path = 'polar_DU95W180.xlsx'

# Read the Excel file into a pandas DataFrame
df = pd.read_excel(file_path)

# Extract columns into NumPy arrays
alfa_array = df['Alfa'].to_numpy()
cl_array = df['Cl'].to_numpy()
cd_array = df['Cd'].to_numpy()
cm_array = df['Cm'].to_numpy()

polar_sheet = np.array([alfa_array, cl_array, cd_array, cm_array])

# Constant
RHO = 1.225 # [kg/m^3]

# Data over alpha plots

# plt.plot(alfa_array, cl_array, label="cl")
# plt.plot(alfa_array, cd_array, label="cd")
# plt.xlabel("Angle of attack [deg]")
# plt.ylabel("Coefficient [-]")
# plt.legend()
# plt.show()


# plt.plot(cd_array, cl_array, label="cd - cl")
# plt.xlabel("cd [-]")
# plt.ylabel("cl [-]")
# plt.legend()
# plt.show()

## Define blade elemenets 

# resolution = 942
resolution = 80
r = np.linspace(design.start, design.end, resolution, endpoint=True) * design.R

# Loop over all segments and take mean conditions for further evaluation (ASSUMPTION)
# n point of evaluation -->> leading to n-1 segments
segment = 0
a_list = []
a_prime_list = []
r_loc = []
for i in np.arange(0, len(r)-1):
    segment += 1
    print(f"Segment number = {segment}")
    segment_start   = r[i]
    segment_end     = r[i + 1]
    segment_dr      = segment_end - segment_start

    segment_mean = statistics.mean([segment_start, segment_end])

    segment_chord = design.chord(segment_mean / design.R)
    segment_twist = design.twist(segment_mean / design.R)

    # For each segment solve the Blade element momentum theory model
    a, a_prime = bem.bem_procedure(design.U0, segment_chord, segment_mean, design.R, design.TSR[1], segment_twist, design.pitch,
                       RHO, polar_sheet, design.BLADES, design.start, segment_dr)
    a_list.append(a)
    a_prime_list.append(a_prime)
    r_loc.append(segment_mean)

plt.plot(r_loc, a_list, label='a')
plt.plot(r_loc, a_prime_list, label='a_prime')
plt.xlabel("Position of blade [m]")
plt.ylabel("Induction factor [-]")
plt.legend()
plt.show()


    