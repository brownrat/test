"""
peaked or "band-pass" spectra. Blue, green and purple objects usually have this type
of reflectance, including green leaves and some flowers. Modeled with a normal distribution.
See https://mathworld.wolfram.com/NormalDistribution.html
The standard deviation is set to 50, and the maximum value is set to 1 rather than the scalar
value included in the normal distribution function. What I really want is to have them be 100 nm
wide at half-maximum, but standard deviation doesn't quite do that.
"""

import central as c
import numpy as np
import math
import matplotlib.pyplot as plt

# normal distribution
def normal(mu, std_dev, x):
	return (1 / std_dev*math.sqrt(2*math.pi)) * math.exp((-1/2)*((x - mu)/std_dev)**2)

spectra = []

# 350-nm peak
peak350 = np.empty(401)
for i in range(401):
        peak350[i] = normal(350, 50, i+300)
spectra.append(peak350 / max(*peak350))

# 400-nm peak
peak400 = np.empty(401)
for i in range(401):
        peak400[i] = normal(400, 50, i+300)
spectra.append(peak400 / max(*peak400))

# 450-nm peak
peak450 = np.empty(401)
for i in range(401):
        peak450[i] = normal(450, 50, i+300)
spectra.append(peak450 / max(*peak450))

# 500-nm peak
peak500 = np.empty(401)
for i in range(401):
        peak500[i] = normal(500, 50, i+300)
spectra.append(peak500 / max(*peak500))

# 550-nm peak
peak550 = np.empty(401)
for i in range(401):
        peak550[i] = normal(550, 50, i+300)
spectra.append(peak550 / max(*peak550))

# 600-nm peak
peak600 = np.empty(401)
for i in range(401):
        peak600[i] = normal(600, 50, i+300)
spectra.append(peak600 / max(*peak600))

# 650-nm peak
peak650 = np.empty(401)
for i in range(401):
        peak650[i] = normal(650, 50, i+300)
spectra.append(peak650 / max(*peak650))

# plots
colors = []
text = []
for i in range(len(spectra)):
	colors.append(c.spec2rgb(spectra[i]))
	text.append(i*50 + 350)

for i in range(len(spectra)):
	plt.plot(c.x_1nm, spectra[i], color=colors[i])
plt.show()

c.triangle(
	spectra=spectra,
	colors=colors,
	gamut=True,
	gamutcolor="0.7",
	gamutedge='',
	text=text,
	legend=True
	)

# narrower peaks
spectra = []

# 325-nm peak
peak325 = np.empty(401)
for i in range(401):
        peak325[i] = normal(325, 25, i+300)
spectra.append(peak325 / max(*peak325))

# 350-nm peak
peak350 = np.empty(401)
for i in range(401):
        peak350[i] = normal(350, 25, i+300)
spectra.append(peak350 / max(*peak350))

# 375-nm peak
peak375 = np.empty(401)
for i in range(401):
        peak375[i] = normal(375, 25, i+300)
spectra.append(peak375 / max(*peak375))

# 400-nm peak
peak400 = np.empty(401)
for i in range(401):
        peak400[i] = normal(400, 25, i+300)
spectra.append(peak400 / max(*peak400))

# 425-nm peak
peak425 = np.empty(401)
for i in range(401):
        peak425[i] = normal(425, 25, i+300)
spectra.append(peak425 / max(*peak425))

# 450-nm peak
peak450 = np.empty(401)
for i in range(401):
        peak450[i] = normal(450, 25, i+300)
spectra.append(peak450 / max(*peak450))

# 475-nm peak
peak475 = np.empty(401)
for i in range(401):
        peak475[i] = normal(475, 25, i+300)
spectra.append(peak475 / max(*peak475))

# 500-nm peak
peak500 = np.empty(401)
for i in range(401):
        peak500[i] = normal(500, 25, i+300)
spectra.append(peak500 / max(*peak500))

# 525-nm peak
peak525 = np.empty(401)
for i in range(401):
        peak525[i] = normal(525, 25, i+300)
spectra.append(peak525 / max(*peak525))

# 550-nm peak
peak550 = np.empty(401)
for i in range(401):
        peak550[i] = normal(550, 25, i+300)
spectra.append(peak550 / max(*peak550))

# 575-nm peak
peak575 = np.empty(401)
for i in range(401):
        peak575[i] = normal(575, 25, i+300)
spectra.append(peak575 / max(*peak575))

# 600-nm peak
peak600 = np.empty(401)
for i in range(401):
        peak600[i] = normal(600, 25, i+300)
spectra.append(peak600 / max(*peak600))

# 625-nm peak
peak625 = np.empty(401)
for i in range(401):
        peak625[i] = normal(625, 25, i+300)
spectra.append(peak625 / max(*peak625))

# 650-nm peak
peak650 = np.empty(401)
for i in range(401):
        peak650[i] = normal(650, 25, i+300)
spectra.append(peak650 / max(*peak650))

# plots
colors = []
text = []
for i in range(len(spectra)):
	colors.append(c.spec2rgb(spectra[i]))
	text.append(i*25 + 325)

for i in range(len(spectra)):
	plt.plot(c.x_1nm, spectra[i], color=colors[i])
plt.show()

c.triangle(
	spectra=spectra,
	colors=colors,
	#gamut=True,
	gamutcolor="0.7",
	gamutedge='',
	text=text,
	legend=True
	)
