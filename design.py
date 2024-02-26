Airfoil = "DU 95-W-180"

R = 50
BLADES = 3

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

    x = []
    y = []
    
    for i in np.arange(start, end, 0.01):
        x.append(i)
        y.append(twist(i))

    plt.plot(x, y)
    plt.show()
