import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import math
import statistics

import design


file_path = 'polar_DU95W180.xlsx'

# Read the Excel file into a pandas DataFrame
df = pd.read_excel(file_path)

# Extract columns into NumPy arrays
alfa_array = df['Alfa'].to_numpy()
cl_array = df['Cl'].to_numpy()
cd_array = df['Cd'].to_numpy()
cm_array = df['Cm'].to_numpy()

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
r_R = np.linspace(design.start, design.end, resolution,endpoint=True)

# Loop over all segments and take mean conditions
for i in range(len(r_R) - 1):

    segment_start   = r_R[i]
    segment_end     = r_R[i + 1]

    segment_mean = statistics.mean([segment_start, segment_end])

    segment_chord = design.chord(segment_mean)
    segment_twist = design.twist(segment_mean)

    # For each segment solve the Blade element momentum theory model

