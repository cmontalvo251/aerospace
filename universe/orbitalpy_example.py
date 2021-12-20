from numpy import radians
from scipy.constants import kilo
import matplotlib.pyplot as plt

from orbital import earth, KeplerianElements, plot
#Here, several convenience methods for defining orbits are shown.

# Create circular orbit with 90 minute period
orbit1 = KeplerianElements.with_period(90 * 60, body=earth)
plot(orbit1, title='Orbit 1')
#plt.plot(orbit1)
