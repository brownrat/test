"""
sky and daylight
The units on these are arbitrary, so they're set to a maximum value of 100.

See also twilight and night skies:
* https://www.researchgate.net/profile/S-Hengst/publication/252090838/figure/fig5/AS:670394587566112@1536845933805/Typical-spectrum-of-the-twilight-sky-from-Dome-A-with-approximate-flux-calibration-using.png
* https://live.staticflickr.com/4642/25663810348_322daa5ced_b.jpg
* https://www.researchgate.net/figure/Relative-and-absolute-downwelling-illumination-at-night-twilight-and-day-a_fig1_303852485

The D65 illuminant peaks in the blue wavelengths (461 nm), but this
represents direct sunlight (sunlight is ~5800 K blackbody radiation) rather
than the color of a blue sky produced by Rayleigh scattering. Blue sky peaks
in or near the UV. The twilight spectra above with similar peaks to D65 may
represent sunlight that passes through the atmosphere without being scattered,
since the blueness of twilight has another cause:
https://en.wikipedia.org/wiki/Chappuis_absorption
"""

import central as c
import numpy as np
import matplotlib.pyplot as plt

daylight_table = np.empty(41)
for i in range(41):
        w = i*10 + 300
        # normalize to 100% at 550
        daylight_table[i] = c.blackbody(w, 5800)
daylight_table = 100 * daylight_table / max(*daylight_table)

sky_table = np.empty(41)
for i in range(41):
        w = i*10 + 300
        sky_table[i] = 100*c.blackbody(w, 5800) * w**-4
sky_table = 100 * sky_table / max(*sky_table)

sky1_table = np.empty(41)
for i in range(41):
        w = i*10 + 300
        sky1_table[i] = 100*c.d65_10nm[i] * w**-4
sky1_table = 100 * sky1_table / max(*sky1_table)

plt.plot(c.x_10nm, c.d65_10nm, 'o-k', label='D65')
plt.plot(c.x_10nm, c.wp_10nm, 'o-w', mec='k', label='current illuminant')
plt.plot(c.x_10nm, daylight_table, 's-', mec='k', color=c.spec2rgb(daylight_table, reflect=False), label="Blackbody daylight")
plt.plot(c.x_10nm, sky_table, 'D-', mec='k', color=c.spec2rgb(sky_table, reflect=False), label="Sky (blackbody)")
plt.plot(c.x_10nm, sky1_table, 'v-', mec='k', color=c.spec2rgb(sky1_table, reflect=False), label="Sky (D65)")
plt.legend()
plt.show()

# triangle
spectra=[daylight_table, sky_table, sky1_table]
colors = []
for i in range(3): colors.append(c.spec2rgb(spectra[i], reflect=False, mode="cmf"))
c.triangle(
	spectra=spectra,
	colors=colors,
	markers=['s', 'D', 'v'],
	reflect=False,
	text=['Blackbody (daylight)', 'Sky (blackbody)', 'Sky (D65)'],
	legend=True,
	gamut=True,
	gamutcolor="0.7",
	gamutedge=''
	)
