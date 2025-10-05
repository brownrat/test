# example illuminants

import central as c
e = c.e
d65 = c.d65
a = c.a
incandescent = c.incandescent

print("White (E, equal-energy)")
c.spectral_rendering(e, reflect=False)

print("White (D65)")
c.spectral_rendering(d65, reflect=False)

print("Incandescent lighting (A)")
c.spectral_rendering(a, reflect=False)

print("Incandescent lighting (approximated with 2856 K blackbody spectrum)")
c.spectral_rendering(incandescent, reflect=False)

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
