from matplotlib import pyplot as plt

from isochrones.files import load_home, load_blue_marble
from isochrones.plot import animate_spherical_fast, animate_spherical_realtime, plot_spherical

load_blue_marble()

# fig, anim = animate_spherical_realtime(home_coord=load_home())
plot_spherical(home_coord=load_home())
plt.show()
