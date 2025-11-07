"""
Three types of optimal colors: see https://en.wikipedia.org/wiki/Gamut#Optimal_colors

* long-pass step functions (0 up to a certain wavelength and then 1)
* short-pass step functions (1 then 0)
* monochromatic single wavelengths (set to 100 because otherwise they come out
really dark)
"""

import central as c
import numpy as np
import matplotlib.pyplot as plt

# "long-pass" step functions: 0 for wavelengths less than the step and 1 for
# longer wavelengths.
# Many white, yellow, orange and red flowers have reflectance spectra like
# these, and they're also found in paints and camera filters.
print("Long-pass step functions")
longpass = []
for i in range(1,int(c.range10/2)):
        step = c.optimal(start=i*20)
        longpass.append(step)

# plots
colors = []
#text = []
for i in range(len(longpass)):
	colors.append(c.spec2rgb(longpass[i]))
	#text.append(i*20 + 20 + c.start)

for i in range(len(longpass)):
	plt.plot(c.x_1nm, longpass[i], color=colors[i])
plt.show()

c.triangle(
	spectra=longpass,
	colors=colors,
	gamut=True,
	gamutcolor="0.7",
	gamutedge='',
	#text=text,
	#legend=True
	)

# "short-pass" step functions, the reverse of the above. These form the other
# edge of the optimal color solid (see Wikipedia) but don't seem to occur in
# nature or have any applications. Blue/green surfaces usually are the
# "band-pass" type with one or two peaks.
print("Short-pass step functions")
shortpass = []
for i in range(1,int(c.range10/2)):
        step = c.optimal(end=i*20)
        shortpass.append(step)

# plots
colors = []
#text = []
for i in range(len(shortpass)):
	colors.append(c.spec2rgb(shortpass[i]))
	#text.append(i*20 + 20 + c.start)

for i in range(len(shortpass)):
	plt.plot(c.x_1nm, shortpass[i], color=colors[i])
plt.show()

c.triangle(
	spectra=shortpass,
	colors=colors,
	gamut=True,
	gamutcolor="0.7",
	gamutedge='',
##	text=text,
##	legend=True
	)

# single wavelengths
# These are bumped up to 50x intensity so they don't look too dark. Actual
# light with only a single precise wavelength would of course be invisible.
print("Single wavelengths")

mono = []
for i in range(c.range1):
        wl = c.optimal(i, i+1)
        mono.append(wl*50)

# plots
colors = []
text = []
for i in range(len(mono)):
	colors.append(c.spec2rgb(mono[i], average=True))

for i in range(len(mono)):
	plt.plot(c.x_1nm, mono[i], color=colors[i])
plt.show()

c.triangle(
	spectra=mono,
	colors=colors,
        mec='#0000'
	)
