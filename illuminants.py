# example illuminants

import central as c
import numpy as np
import colormath
from colormath import spectral_constants
e = c.e_10nm
d65 = c.d65_10nm
a = np.zeros(41)
for i in range(4, 41):
	a[i] = spectral_constants.REF_ILLUM_TABLE["a"][i-4]
incandescent = c.incandescent_10nm

c.triangle(
	[e, d65, a, incandescent],
	markers=['o', 's', 'D', 'v'],
	colors=[c.spec2rgb(e), c.spec2rgb(d65), c.spec2rgb(a), c.spec2rgb(incandescent)],
	text=['E', 'D65', 'A', 'incandescent'],
	legend=True,
	reflect=False)

# color/brightness contrast test
print("Color contrast between E and E: " + str(c.color_contrast(e, e)))
print("Color contrast between E and D65: " + str(c.color_contrast(e, d65)))
print("Color contrast between E and A: " + str(c.color_contrast(e, a)))
print("Color contrast between E and incandescent: " + str(c.color_contrast(e, incandescent)))
print("Color contrast between D65 and D65: " + str(c.color_contrast(d65, d65)))
print("Color contrast between D65 and A: " + str(c.color_contrast(d65, a)))
print("Color contrast between D65 and incandescent: " + str(c.color_contrast(d65, incandescent)))
print("Color contrast between A and A: " + str(c.color_contrast(a, a)))
print("Color contrast between A and incandescent: " + str(c.color_contrast(a, incandescent)))
print("Color contrast between incandescent and incandescent: " + str(c.color_contrast(incandescent, incandescent)))
print("Brightness contrast between E and E: " + str(c.brightness_contrast(e, e)))
print("Brightness contrast between E and D65: " + str(c.brightness_contrast(e, d65)))
print("Brightness contrast between E and A: " + str(c.brightness_contrast(e, a)))
print("Brightness contrast between E and incandescent: " + str(c.brightness_contrast(e, incandescent)))
print("Brightness contrast between D65 and D65: " + str(c.brightness_contrast(d65, d65)))
print("Brightness contrast between D65 and A: " + str(c.brightness_contrast(d65, a)))
print("Brightness contrast between D65 and incandescent: " + str(c.brightness_contrast(d65, incandescent)))
print("Brightness contrast between A and A: " + str(c.brightness_contrast(a, a)))
print("Brightness contrast between A and incandescent: " + str(c.brightness_contrast(a, incandescent)))
print("Brightness contrast between incandescent and incandescent: " + str(c.brightness_contrast(incandescent, incandescent)))
