Airfoil = "DU 95-W-180"

R = 50          # Blade length [m]
BLADES = 3      # Numbder of blades [-]

start = 0.2     # start position of the actual airfoil section [r_over_R]
end = 1         # End postion of the blades [r_over_R]

twist = lambda r_over_R : 14 * (1 - r_over_R)       # [deg]
chord = lambda r_over_R : (3 * (1 - r_over_R) + 1)  # [m]
pitch = - 2     # [deg]

U0 = 10             # [m/s]
TSR = [6, 8, 10]
YAW = [0, 15, 30]   # Yaw expressed in degress 

if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plt

    x = np.linspace(start, end, 1000)
    y = twist(x)

    plt.plot(x, y)
    plt.show()
