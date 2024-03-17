print('Loading libraries...')
from matplotlib import pyplot as plt
from .files import load_home, load_blue_marble
from .plot import (
    animate_spherical_fast, animate_spherical_realtime, plot_spherical,
)

print('Loading files...')
load_blue_marble()
home_coord = load_home()

print('Setting up artists...')
# fig, anim = animate_spherical_realtime(home_coord=load_home())

# The update step has negligible time, so we don't split these out to log them separately
fig = plot_spherical(home_coord=home_coord)

print('Plotting...')
plt.show()
