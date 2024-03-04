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

plt.plot(alfa_array, cl_array, label="cl")
plt.plot(alfa_array, cd_array, label="cd")
plt.xlabel("Angle of attack [deg]")
plt.ylabel("Coefficient [-]")
plt.legend()
plt.show()


plt.plot(cd_array, cl_array, label="cd - cl")
plt.xlabel("cd [-]")
plt.ylabel("cl [-]")
plt.legend()
plt.show()

## Define blade elemenets 

resolution = 1000
r = np.linspace(design.start, design.end, resolution, endpoint=True) * design.R

# Loop over all segments and take mean conditions for further evaluation (ASSUMPTION)
# n point of evaluation -->> leading to n-1 segments
for i in [0]:

    segment_start   = r[i]
    segment_end     = r[i + 1]
    segment_dr      = segment_end - segment_start

    segment_mean = statistics.mean([segment_start, segment_end])

    segment_chord = design.chord(segment_mean)
    segment_twist = design.twist(segment_mean)

    # For each segment solve the Blade element momentum theory model
    bem.bem_procedure(design.U0, segment_chord, segment_mean, design.R, design.TSR[1], segment_twist, design.twist,
                       RHO, polar_sheet, design.BLADES, design.start, segment_dr)
    