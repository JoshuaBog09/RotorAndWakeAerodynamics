import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

import math
import statistics

import design
import bem
import stagnation
import seaborn as sns

sns.set_theme(rc={'grid.color': 'white'})
matplotlib.rcParams.update({'font.size': 15, 'figure.figsize': (6, 5), 'figure.dpi': 150})
# matplotlib.rcParams['figure.subplot.top'] = 0.95
# matplotlib.rcParams['figure.subplot.right'] = 0.98

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
yaw_angle = design.YAW[1]

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
r = np.linspace(design.start, design.end, resolution, endpoint=True) * design.R  # Equal spacing
# r = (np.sin(np.linspace(0, np.pi/2, resolution, endpoint=True)) * (design.end - design.start) + design.start)* design.R
#          # High density @ tip, low density @ root
# r = (((np.cos(np.linspace(np.pi/2, np.pi, resolution, endpoint=True)) + design.end) * (design.end - design.start)) + design.start)* design.R
#          # High density @ root, low density @ tip
# r = ((((np.cos(np.linspace(0, np.pi, resolution, endpoint=True)) + design.end) / 2 ) * (design.end - design.start)) + design.start)* design.R
#          # High @ both endpoints
# r = (((((np.linspace(-1, 1, resolution, endpoint=True)**3) + design.end) / 2) * (design.end - design.start)) + design.start)* design.R
         # High @ central points

# Loop over each TSR value
all_results = []
tot_thrust_list, tot_power_list, tot_torque_list = [], [], []
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
    CT_list = []
    CP_list = []
    cl_list = []
    cd_list = []
    r_loc_list = []
    f_azi_list, f_axi_list = [], []
    r_inner_list = []
    r_outer_list = []

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
        # a, a_prime, phi, alpha, f_azi, f_axi, ct, cp
        output = bem.bem_procedure(design.U0, segment_chord, segment_mean,
                                 design.R, tsr_value,
                                 segment_twist, design.pitch,
                                 RHO, polar_sheet, design.BLADES, design.start,
                                 segment_dr, 0.0001,
                                 np.radians(yaw_angle))
        a_list.append(output['a_new'])
        a_prime_list.append(output['a_prime_new'])
        phi_list.append(output['Phi'])
        alpha_list.append(output['alpha'])
        r_loc_list.append(segment_mean / design.R)
        f_azi_list.append(output['f_azi'])  # f_tan
        f_axi_list.append(output['f_axi'])  # f_norm
        cl_list.append(output['cl'])
        cd_list.append(output['cd'])
        CT_list.append(output['CT'])
        CP_list.append(output['CP'])
        
        r_inner_list.append(segment_start)
        r_outer_list.append(segment_end)

    fig_stat, fig_stag = stagnation.compute_stagstat_pressure(design.U0, a_list, r_inner_list, r_outer_list)

    # Calculating the total thrust, power and torque
    total_thrust = np.sum(design.BLADES * np.array(f_axi_list) * segment_dr)
    total_power = np.sum(
        segment_dr * design.BLADES * np.array(f_azi_list) * np.array(r_loc_list) * tsr_value / design.U0)
    total_torque = np.sum(
        segment_dr * design.BLADES * np.array(f_azi_list) * (np.array(r_loc_list) * design.R))
    C_T = total_thrust / (0.5 * RHO * design.U0 ** 2 * np.pi * design.R ** 2)
    C_P = total_power / (0.5 * RHO * design.U0 ** 3 * np.pi * design.R ** 2)
    C_Q = total_torque / (0.5 * RHO * design.U0 ** 2 * np.pi * design.R ** 3)
    tot_thrust_list.append(total_thrust), tot_power_list.append(total_torque), tot_torque_list.append(total_torque)
    C_T_list.append(C_T), C_P_list.append(C_P), C_Q_list.append(C_Q)
    print(C_P)

    all_results.append({
        'tsr': tsr_value,
        'a': a_list,
        'a_prime': a_prime_list,
        'phi': phi_list,
        'alpha': alpha_list,
        'r/R_loc': r_loc_list,
        'f_azi': np.array(f_azi_list),
        'f_axi': np.array(f_axi_list),
        'cl': np.array(cl_list),
        'cd': np.array(cd_list),
        'CP': np.array(CP_list),
    })

blue_color = ['#00CCFF', '#005FFF', '#0000FF']
red_color = ['#FFCCCC', '#FF5F5F', '#FF0000']

fig_induction, ax_induction = plt.subplots(nrows=1, ncols=1)
for i, result in enumerate(all_results):
    ax_induction.plot(result['r/R_loc'], result['a'], label=f"TSR: {result['tsr']}, a", color=blue_color[i])
    ax_induction.plot(result['r/R_loc'], result['a_prime'], label=f"TSR: {result['tsr']}, a'", color=red_color[i])
ax_induction.set_xlabel("r/R")
ax_induction.set_ylabel("Induction factor [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_induction.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
ax_induction.set_title(f'Induction factor for TSR: 6,8,10')
ax_induction.set_ylim([0, 1])
# ax_induction.grid()
# fig_induction.tight_layout()
# plt.show()

fig_phi, ax_phi = plt.subplots(nrows=1, ncols=1)
for i, result in enumerate(all_results):
    ax_phi.plot(result['r/R_loc'], np.degrees(result['phi']), label=f"TSR: {result['tsr']}, $\\phi$", color=blue_color[i])
    ax_phi.plot(result['r/R_loc'], np.degrees(result['alpha']), label=f"TSR: {result['tsr']}, $\\alpha$",
                color=red_color[i])
ax_phi.set_xlabel("r/R")
ax_phi.set_ylabel(r"$\phi$, $\alpha$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_phi.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
# ax_phi.grid()
# fig_phi.tight_layout()
# plt.show()

fig_force, ax_force = plt.subplots(nrows=1, ncols=1)
n_dim_force = 0.5 * RHO * design.U0 ** 2 * (2*np.pi*np.array(result['r/R_loc'])*design.R)*segment_dr
for i, result in enumerate(all_results):
    ax_force.plot(result['r/R_loc'], result['f_azi'] / n_dim_force, label=f"TSR: {result['tsr']}, $C_{{t}}$",
                  color=blue_color[i])
    ax_force.plot(result['r/R_loc'], result['f_axi'] / n_dim_force, label=f"TSR: {result['tsr']}, $C_{{n}}$",
                  color=red_color[i])
ax_force.set_xlabel("r/R")
ax_force.set_ylabel(r"$C_t$, $C_n$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_force.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
# ax_force.grid()
# fig_force.tight_layout()
# plt.show()

fig_thrust, ax_thrust = plt.subplots(nrows=1, ncols=1)
fig_torque, ax_torque = plt.subplots(nrows=1, ncols=1)
ax_thrust.scatter(design.TSR, C_T_list, label=r"$C_T$", color='blue')
ax_thrust.scatter(design.TSR, C_P_list, label=r"$C_P$", color='green')
ax_thrust.set_xlabel("TSR [-])")
ax_thrust.set_ylabel(r"$C_T$, $C_P$ [-]")
ax_thrust.legend()
# ax_thrust.grid()
# ax_thrust.set_title("Total Thrust")
ax_torque.scatter(design.TSR, C_Q_list, label=r"$C_Q$", color='red')
ax_torque.set_xlabel("TSR [-]")
ax_torque.set_ylabel(r"$C_Q$ [-]")
# ax_torque.grid()
# ax_torque.set_title("Total Torque")
# fig_tot_force.tight_layout()

fig_lift, ax_lift = plt.subplots(nrows=1, ncols=1)
fig_drag, ax_drag = plt.subplots(nrows=1, ncols=1)
n_dim_force = 0.5 * RHO * design.U0 ** 2 * (2*np.pi*np.array(result['r/R_loc'])*design.R)*segment_dr
for i, result in enumerate(all_results):
    ax_lift.plot(result['r/R_loc'], result['cl'], label=f"TSR: {result['tsr']}, $C_{{l}}$",
                  color=blue_color[i])
    ax_drag.plot(result['r/R_loc'], result['cd'], label=f"TSR: {result['tsr']}, $C_{{d}}$",
                  color=red_color[i])
ax_lift.set_xlabel("r/R [-]")
ax_lift.set_ylabel(r"$C_l$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_lift.legend()
ax_lift.set_ylim(0, 1.3)
# ax_lift.grid()
# fig_lift.tight_layout()

ax_drag.set_xlabel("r/R [-]")
ax_drag.set_ylabel(r"$C_d$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_drag.legend()
ax_drag.set_ylim(0, 0.3)
# ax_drag.grid()
# fig_drag.tight_layout()

fig_CP, ax_CP = plt.subplots(nrows=1, ncols=1)
for i, result in enumerate(all_results):
    ax_CP.plot(result['r/R_loc'], result['CP'], label=f"TSR: {result['tsr']}, $C_{{P}}$",
                  color=blue_color[i])
ax_CP.set_xlabel("r/R [-]")
ax_CP.set_ylabel(r"$C_P$ [-]")
# Rearrange legend
handles, labels = plt.gca().get_legend_handles_labels()
ax_CP.legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2])
ax_CP.set_ylim(0, 0.6)
# ax_CP.grid()
# fig_CP.tight_layout()

# Data over alpha plots

fig_lift_drag_polar, ax_lift_drag_polar = plt.subplots(nrows=1, ncols=1)
ax_lift_drag_polar.plot(alfa_array, cl_array, label="$C_l$")
ax_lift_drag_polar.plot(alfa_array, cd_array, label="$C_d$")
# ax_lift_drag_polar.set_xlim(-30, 30)
ax_lift_drag_polar.set_xlabel(r'$\alpha$ [deg]')
ax_lift_drag_polar.set_ylabel('$C_{l}$, $C_{d}$ [-]')
# ax_lift_drag_polar.legend()
# ax_lift_drag_polar.grid()
# fig_lift_drag_polar.tight_layout()

fig_lift_drag_bucket, ax_lift_drag_bucket = plt.subplots(nrows=1, ncols=1)
ax_lift_drag_bucket.plot(cd_array, cl_array)
ax_lift_drag_bucket.set_xlim(0, 0.15)
ax_lift_drag_bucket.set_xlabel(r"$C_d$ [-]")
ax_lift_drag_bucket.set_ylabel(r"$C_l$ [-]")
# ax_lift_drag_bucket.legend()
# ax_lift_drag_bucket.grid()
# fig_lift_drag_bucket.tight_layout()

# path = r'C:\Users\roelv\Documents\School\TU_Delft\MSC_1\Physics\AE4135_Rotor_Wake_Aerodynamics\RWA_figures'
# fig_CP.savefig(rf'{path}\CP_rR_yaw={yaw_angle}.pdf')
# fig_induction.savefig(rf'{path}\a_rR_yaw={yaw_angle}.pdf')
# fig_phi.savefig(rf'{path}\phi_alpha_rR_yaw={yaw_angle}.pdf')
# fig_drag.savefig(rf'{path}\Cd_rR_yaw={yaw_angle}.pdf')
# fig_lift.savefig(rf'{path}\Cl_rR_yaw={yaw_angle}.pdf')
# fig_lift_drag_polar.savefig(rf'{path}\Cl_Cd_polar.pdf')
# fig_lift_drag_bucket.savefig(rf'{path}\Cl_Cd_bucket.pdf')
# fig_force.savefig(rf'{path}\Ct_Cn_rR_yaw={yaw_angle}.pdf')
# fig_thrust.savefig(rf'{path}\Trust_TSR_yaw={yaw_angle}.pdf')
# fig_torque.savefig(rf'{path}\Torque_rR_yaw={yaw_angle}.pdf')
# fig_stat.savefig(rf'{path}\p_stat_distr_yaw={yaw_angle}.pdf')
# fig_stag.savefig(rf'{path}\p_stag_distr_yaw={yaw_angle}.pdf')


plt.show()

# plt.show()