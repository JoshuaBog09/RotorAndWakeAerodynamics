import numpy as np
import matplotlib.pyplot as plt

import design 

segments = 25
segmentation1 = np.linspace(design.start, design.end, segments, endpoint=True)   # Equal spacing
segmentation2 = np.sin(np.linspace(0, np.pi/2, segments, endpoint=True)) * (design.end - design.start) + design.start
         # High density @ tip, low density @ root
segmentation3 = ((np.cos(np.linspace(np.pi/2, np.pi, segments, endpoint=True)) + design.end) * (design.end - design.start)) + design.start
         # High density @ root, low density @ tip
segmentation4 = (((np.cos(np.linspace(0, np.pi, segments, endpoint=True)) + design.end) / 2 ) * (design.end - design.start)) + design.start
         # High @ both endpoints
segmentation5 = ((((np.linspace(-1, 1, segments, endpoint=True)**3) + design.end) / 2) * (design.end - design.start)) + design.start
         # High @ central points

plt.scatter(segmentation1, np.zeros(len(segmentation1)) + 0)
plt.scatter(segmentation2, np.zeros(len(segmentation2)) + 1)
plt.scatter(segmentation3, np.zeros(len(segmentation3)) + 2)
plt.scatter(segmentation4, np.zeros(len(segmentation4)) + 3)
plt.scatter(segmentation5, np.zeros(len(segmentation5)) + 4)
plt.show()

