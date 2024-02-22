from datetime import datetime, timezone

from matplotlib import pyplot as plt

from isochrones.files import load_home
from isochrones.plot import plot_spherical


utcnow = datetime.now().astimezone(timezone.utc)
home_coord = load_home()
plot_spherical(home_coord=home_coord, utcnow=utcnow)
plt.show()
