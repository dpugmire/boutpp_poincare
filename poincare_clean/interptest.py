import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt

# Original data
x = np.arange(0, 6)  # x-coordinates
y = np.arange(0, 6)  # y-coordinates
X, Y = np.meshgrid(x, y)
Z = np.sin(X) * np.cos(Y)  # Some sample 2D data

# New points for interpolation
xq = np.linspace(0, 5, 50)  # Query x-coordinates
yq = np.linspace(0, 5, 50)  # Query y-coordinates

# Spline interpolation
spline = RectBivariateSpline(x, y, Z)
Zq = spline(xq, yq)

# Plot original and interpolated data
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
axes[0].contourf(X, Y, Z, cmap='viridis')
axes[0].set_title('Original Data')
axes[1].contourf(xq, yq, Zq, cmap='viridis')
axes[1].set_title('Interpolated Data')
plt.show()
