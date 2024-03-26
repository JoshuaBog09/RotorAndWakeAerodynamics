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
tot_thrust_list, tot_power_list, tot_torque_list,  = [], [], []
C_T_list, C_P_list, C_Q_list = [], [], []
for tsr_value in design.TSR:
    print(f"Calculations for TSR: {tsr_value}")

    # Loop over all segments and take mean conditions for further evaluation (ASSUMPTION)
    # n point of evaluation -->> leading to n-1 segments
    segment = 0
    a_list = []
    a_prime_list = []
    phi_list = []
    alpha_list = []
    r_loc_list = []
    f_azi_list, f_axi_list = [], []
    for i in np.arange(0, len(r) - 1):
        segment += 1
        # print(f"Segment number = {segment}")
        segment_start = r[i]
        segment_end = r[i + 1]
        segment_dr = segment_end - segment_start

        segment_mean = statistics.mean([segment_start, segment_end])

        segment_chord = design.chord(segment_mean / design.R)
        segment_twist = design.twist(segment_mean / design.R)

        # For each segment solve the Blade element momentum theory model
        a, a_prime, phi, alpha, f_azi, f_axi = bem.bem_procedure(design.U0, segment_chord, segment_mean, design.R, tsr_value,
                                                   segment_twist, design.pitch,
                                                   RHO, polar_sheet, design.BLADES, design.start, segment_dr, 0.0001,
                                                   np.radians(design.YAW[0]))
        a_list.append(a)
        a_prime_list.append(a_prime)
        phi_list.append(phi)
        alpha_list.append(alpha)
        r_loc_list.append(segment_mean / design.R)
        f_azi_list.append(f_azi)  # f_tan
        f_axi_list.append(f_axi)  # f_norm

    # Calculating the total thrust, power and torque
    total_thrust = np.sum(design.BLADES * np.array(f_axi_list) * segment_dr)
    total_power = np.sum(segment_dr * design.BLADES * np.array(f_azi_list) * np.array(r_loc_list) * tsr_value/design.U0)
    total_torque = np.sum(segment_dr * design.BLADES * np.array(f_azi_list) * (np.array(r_loc_list) * design.R))
    C_T = total_thrust / (0.5 * RHO * design.U0 ** 2 * np.pi * design.R ** 2)
    C_P = total_power / (0.5 * RHO * design.U0 ** 3 * np.pi * design.R ** 2)
    C_Q = total_torque / (0.5 * RHO * design.U0 ** 2 * np.pi *design.R ** 3)
    tot_thrust_list.append(total_thrust), tot_power_list.append(total_power), tot_torque_list.append(total_torque)
    C_T_list.append(C_T), C_P_list.append(C_P), C_Q_list.append(C_Q)

    all_results.append({
        'tsr': tsr_value,
        'a': a_list,
        'a_prime': a_prime_list,
        'phi': phi_list,
        'alpha': alpha_list,
        'r_loc': r_loc_list,
        'f_azi': np.array(f_azi_list),
        'f_axi': np.array(f_axi_list),
    })

blue_color = ['#00CCFF', '#005FFF', '#0000FF']
red_color = ['#FFCCCC', '#FF5F5F', '#FF0000']

fig_induction, ax_induction = plt.subplots(nrows=1, ncols=1)
for i, result in enumerate(all_results):
    ax_induction.plot(result['r_loc'], result['a'], label=f"TSR: {result['tsr']}, a", color=blue_color[i])
    ax_induction.plot(result['r_loc'], result['a_prime'], label=f"TSR: {result['tsr']}, a'", color=red_color[i])
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
    ax_phi.plot(result['r_loc'], np.degrees(result['phi']), label=f"TSR: {result['tsr']}, $\\phi$", color=blue_color[i])
    ax_phi.plot(result['r_loc'], np.degrees(result['alpha']), label=f"TSR: {result['tsr']}, $\\alpha$", color=red_color[i])
ax_phi.set_xlabel("Normalized position of blade (r/R)")
ax_phi.set_ylabel(r"$\phi$, $\alpha$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_phi.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
ax_phi.grid()
# plt.show()

fig_force, ax_force = plt.subplots(nrows=1, ncols=1)
n_dim_force = 0.5*RHO*design.U0**2*design.R
for i, result in enumerate(all_results):
    ax_force.plot(result['r_loc'], result['f_azi']/n_dim_force, label=f"TSR: {result['tsr']}, $F_{{tan}}$", color=blue_color[i])
    ax_force.plot(result['r_loc'], result['f_axi']/n_dim_force, label=f"TSR: {result['tsr']}, $F_{{norm}}$", color=red_color[i])
ax_force.set_xlabel("Normalized position of blade (r/R)")
ax_force.set_ylabel(r"$C_t$, $C_n$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_force.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
ax_force.grid()
plt.show()

fig_tot_force, (ax_thrust, ax_torque) = plt.subplots(nrows=1, ncols=2)
ax_thrust.scatter(design.TSR, tot_thrust_list, label=r"$T_{total}$", color='blue')
ax_thrust.set_xlabel("Tip speed ratio (TSR)")
ax_thrust.set_ylabel(r"$T_{tot}$ [N]")
ax_thrust.grid()
ax_thrust.set_title("Total Thrust")
ax_torque.scatter(design.TSR, tot_torque_list, label=r"$Q_{total}$", color='red')
ax_torque.set_xlabel("Tip speed ratio (TSR)")
ax_torque.set_ylabel(r"$Q_{tot}$ [Nm]")
ax_torque.grid()
ax_torque.set_title("Total Torque")
plt.tight_layout()

fig_C_force, (ax_C_thrust, ax_C_torque) = plt.subplots(nrows=1, ncols=2)
ax_C_thrust.scatter(design.TSR, C_T_list, label=r"$T_{total}$", color='blue')
ax_C_thrust.set_xlabel("Tip speed ratio (TSR)")
ax_C_thrust.set_ylabel(r"$T_{tot}$ [N]")
ax_C_thrust.grid()
ax_C_thrust.set_title("Total Thrust")
ax_C_torque.scatter(design.TSR, C_Q_list, label=r"$Q_{total}$", color='red')
ax_C_torque.set_xlabel("Tip speed ratio (TSR)")
ax_C_torque.set_ylabel(r"$Q_{tot}$ [Nm]")
ax_C_torque.grid()
ax_C_torque.set_title("Total Torque")
plt.tight_layout()
plt.show()