from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
import numpy as np

# Fixing random state for reproducibility


arr = np.loadtxt("field.txt")
for i in range(120):
    for j in range(120):
        if arr[i,j] > 100:
            arr[i,j] = 100;



fig = plt.figure()
ax = fig.add_subplot(111)

plt.imshow(arr)


plt.show()