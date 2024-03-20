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
RHO = 1.225  # [kg/m^3]

# Data over alpha plots

# plt.plot(alfa_array, cl_array, label="$C_L$")
# plt.plot(alfa_array, cd_array, label="$C_D$")
# plt.xlim(-30, 30)
# plt.xlabel("Angle of attack [deg]")
# plt.ylabel("Coefficient [-]")
# plt.legend()
# plt.grid()
# plt.show()
#
#
# plt.plot(cd_array, cl_array, label="$C_D$ vs. $C_L$")
# plt.xlim(0, 0.1)
# plt.xlabel(r"$C_D$ [-]")
# plt.ylabel(r"$C_L$ [-]")
# plt.legend()
# plt.grid()
# plt.show()

## Define blade elemenets 

# resolution = 942
resolution = 100
r = np.linspace(design.start, design.end, resolution, endpoint=True) * design.R

# Loop over each TSR value
all_results = []
for tsr_value in design.TSR:
    print(f"Calculations for TSR: {tsr_value}")

    # Loop over all segments and take mean conditions for further evaluation (ASSUMPTION)
    # n point of evaluation -->> leading to n-1 segments
    segment = 0
    a_list = []
    a_prime_list = []
    phi_list = []
    alpha_list = []
    r_loc = []
    for i in np.arange(0, len(r) - 1):
        segment += 1
        print(f"Segment number = {segment}")
        segment_start = r[i]
        segment_end = r[i + 1]
        segment_dr = segment_end - segment_start

        segment_mean = statistics.mean([segment_start, segment_end])

        segment_chord = design.chord(segment_mean / design.R)
        segment_twist = design.twist(segment_mean / design.R)

        # For each segment solve the Blade element momentum theory model
        a, a_prime, phi, alpha = bem.bem_procedure(design.U0, segment_chord, segment_mean, design.R, tsr_value,
                                                   segment_twist, design.pitch,
                                                   RHO, polar_sheet, design.BLADES, design.start, segment_dr, 0.0001,
                                                   np.radians(design.YAW[2]))
        a_list.append(a)
        a_prime_list.append(a_prime)
        phi_list.append(phi)
        alpha_list.append(alpha)
        r_loc.append(segment_mean / design.R)

    all_results.append({
        'tsr': tsr_value,
        'a_list': a_list,
        'a_prime_list': a_prime_list,
        'phi_list': phi_list,
        'alpha_list': alpha_list,
        'r_loc': r_loc
    })

blue_color = ['#00CCFF', '#005FFF', '#0000FF']
red_color = ['#FFCCCC', '#FF5F5F', '#FF0000']

fig_induction, ax_induction = plt.subplots(nrows=1, ncols=1)
for i, result in enumerate(all_results):
    ax_induction.plot(result['r_loc'], result['a_list'], label=f"TSR: {result['tsr']}, a", color=blue_color[i])
    ax_induction.plot(result['r_loc'], result['a_prime_list'], label=f"TSR: {result['tsr']}, a'", color=red_color[i])
ax_induction.set_xlabel("Normalized position of blade (r/R)")
ax_induction.set_ylabel("Induction factor [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_induction.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
ax_induction.set_title(f'Induction factor for TSR: 6,8,10')
ax_induction.set_ylim([0, 1])
ax_induction.grid()
# plt.show()

fig_phi, ax_phi = plt.subplots(nrows=1, ncols=1)
for i, result in enumerate(all_results):
    ax_phi.plot(result['r_loc'], np.degrees(result['phi_list']), label=f"TSR: {result['tsr']}, $\\phi$", color=blue_color[i])
    ax_phi.plot(result['r_loc'], np.degrees(result['alpha_list']), label=f"TSR: {result['tsr']}, $\\alpha$", color=red_color[i])
ax_phi.set_xlabel("Normalized position of blade (r/R)")
ax_phi.set_ylabel(r"$\phi$, $\alpha$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_phi.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
ax_phi.grid()
plt.show()

