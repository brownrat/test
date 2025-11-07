"""
Peaked or "band-pass" spectra. Blue, green and purple objects usually have this
type of reflectance, including green leaves and some flowers. See also:
https://mathworld.wolfram.com/GaussianFunction.html
"""

import central as c
import numpy as np
import math
import matplotlib.pyplot as plt

spectra = []
for i in range(17):
        peak = c.gaussian(i*25+300, 100)
        spectra.append(peak)

# plots
colors = []
text = []
for i in range(len(spectra)):
	colors.append(c.spec2rgb(spectra[i]))
	text.append(i*25 + 300)

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
for i in range(1,20):
        peak = c.gaussian(i*20+300, 50)
        spectra.append(peak*2)

# plots
colors = []
text = []
for i in range(len(spectra)):
	colors.append(c.spec2rgb(spectra[i]))
	text.append(i*20 + 320)

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

# even narrower
spectra = []
for i in range(41):
        peak = c.gaussian(i*10+300, 20)
        spectra.append(peak*10)

# plots
colors = []
text = []
for i in range(len(spectra)):
	colors.append(c.spec2rgb(spectra[i]))
	text.append(i*10 + 300)

for i in range(len(spectra)):
	plt.plot(c.x_1nm, spectra[i], color=colors[i])
plt.show()

c.triangle(
	spectra=spectra,
	colors=colors,
	#gamut=True,
	gamutcolor="0.7",
	gamutedge='',
	#text=text,
	#legend=True
	)
