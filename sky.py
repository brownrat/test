"""
sky and daylight
These are both normalized to "1" at 550 nm. If I try to make "sky" relative to "daylight",
the values are way too low and probably not meaningful.

See also twilight and night skies:
* https://www.researchgate.net/profile/S-Hengst/publication/252090838/figure/fig5/AS:670394587566112@1536845933805/Typical-spectrum-of-the-twilight-sky-from-Dome-A-with-approximate-flux-calibration-using.png
* https://live.staticflickr.com/4642/25663810348_322daa5ced_b.jpg
* https://www.researchgate.net/figure/Relative-and-absolute-downwelling-illumination-at-night-twilight-and-day-a_fig1_303852485

The D65 illuminant peaks in the blue wavelengths (about 480 nm), but this
represents direct sunlight (sunlight is ~5800 K blackbody radiation) rather
than the color of a blue sky produced by Rayleigh scattering. Blue sky peaks
in or near the UV. The twilight spectra above that peak near 450 nm presumably
represent sunlight that passes through the atmosphere without being scattered,
since the blueness of twilight has another cause:
https://en.wikipedia.org/wiki/Chappuis_absorption
"""

import central as c
import numpy as np

daylight_table = np.empty(41)
for i in range(41):
        w = i*10 + 300
        # normalize to 100% at 550
        daylight_table[i] = 100*c.blackbody(w, 5800) / c.blackbody(550, 5800)

print("Black body spectrum approximating daylight")
c.spectral_rendering(daylight_table, reflect=False)

sky_table = np.empty(41)
for i in range(0, 41):
        w = i*10 + 300
        sky_table[i] = 100*c.blackbody(w, 5800) * w**-4 / (c.blackbody(550, 5800) * 550**-4)
        
print("Black body spectrum with Rayleigh scattering approximating a blue sky")
c.spectral_rendering(sky_table, reflect=False)

sky1_table = np.empty(41)
for i in range(0, 41):
        w = i*10 + 300
        sky1_table[i] = 100*c.d65[i] * w**-4 * 5e8

print("D65 with Rayleigh scattering")
c.spectral_rendering(sky1_table, reflect=False)
