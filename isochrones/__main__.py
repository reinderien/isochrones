from matplotlib import pyplot as plt

from isochrones.files import load_home
from isochrones.plot import animate_spherical_realtime, animate_spherical_fast, plot_spherical


# fig, anim = animate_spherical_realtime(home_coord=load_home())
plot_spherical(home_coord=load_home())
plt.show()
