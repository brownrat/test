"""
10/10/2025 -- bringing the 24 May 2024 version back as a separate file so OpenCV doesn't have to share
space with matplotlib.

We skip the step of finding a dominant wavelength and directly transform the source chromaticity into
the target chromaticity, as described here:

* https://www.empiricalimaging.com/knowledge-base/spectral-sensitivity-based-cone-catch-models/
* https://doi.org/10.1111/2041-210X.12439

Interactions can be set to either one-way (none; default) or two-way. I couldn't get two-way fits to
produce good results for anything other than pure hues like the monochromatic wavelengths we fit
the curves on, so like CVS we separate hue from saturation before converting the color. The final
saturation is altered if the smallest output rgb value is not equal to 0 (>0 = less, <0 = more). I
also convert the curves to chromaticity before fitting them because two-way interactions won't come
out right otherwise. This avoids the problem CVS seems to have with dichromatic simulations where
setting the L/M cone close to the human L cone removes the green channel and turns blue-green areas
pure blue (earlier versions of my code also did this, as well as the reverse where setting it close
to or below the M cone removes the red channel).

Luminance is held constant by comparing the original XYZ Y value with the final value. This can
possibly be changed to be relative to the target's luminosity function.
"""

import math
import numpy as np
import time
import argparse
import cv2
import central as c
args = c.args

# execution time
start_time = time.time()

# get arrays we generated in central.py
red = c.sr1
green = c.sr2
blue = c.sr3
uv = np.zeros(401) # not yet implemented

wptr1 = c.wptr1
wptr2 = c.wptr2
wptr3 = c.wptr3

# image processing
imagename = args.image
# Convert image from BGR to HLS. We use BGR because this is how OpenCV reads images.
# If we use RGB, the output appears normal with settings approximating human vision,
# but shifting the cones produces the opposite of the expected result.
img = cv2.imread(imagename)
img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

# transform hues for each pixel
for x in range(img.shape[0]):
	for y in range(img.shape[1]):
		# report current pixel (for testing)
		#print("pixel: " + str(x) + ", " + str(y))

		# get pixel
		pixel = img_rgb[x][y]
		
		# 0-255 -> 0-1
		rgb_r = pixel[0] / 255
		rgb_g = pixel[1] / 255
		rgb_b = pixel[2] / 255
		
		# undo gamma
		rgb_r = c.linearize(rgb_r)
		rgb_g = c.linearize(rgb_g)
		rgb_b = c.linearize(rgb_b)
		
		# luminance
		lum_in = rgb_r*0.2126729 + rgb_g*0.7151522 + rgb_b*0.0721750
		
		# saturate
		max_rgb = max(rgb_r,rgb_g,rgb_b)
		min_rgb = min(rgb_r,rgb_g,rgb_b)
		if (max_rgb > 0): sat_in = min_rgb / max_rgb
		else: sat_in = 0
		rgb_r -= min_rgb
		rgb_g -= min_rgb
		rgb_b -= min_rgb
		max_rgb_diff = max(rgb_r,rgb_g,rgb_b)
		
		# chromaticity
		total = rgb_r + rgb_g + rgb_b
		if (total > 0):
			rgb_r /= total
			rgb_g /= total
			rgb_b /= total
		
		# combine
		rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_r*rgb_g, rgb_r*rgb_b, rgb_g*rgb_b]
		rgb_r1 = sum(rgb_terms * c.s2t_r1)
		rgb_r2 = sum(rgb_terms * c.s2t_r2)
		rgb_r3 = sum(rgb_terms * c.s2t_r3)
		
		# remove junk channels for lower dimensions
		if (c.dimension == 2):
			rgb_r3 = rgb_r2
			# standard tritanopia colors
			if (args.receptors[1] < 490): rgb_r2 = rgb_r1
		if (c.dimension == 1):
			rgb_r2 = rgb_r1
			rgb_r3 = rgb_r1
		
		# desaturate
		max_rgb = max(rgb_r1,rgb_r2,rgb_r3)
		if (max_rgb > 0):
			rgb_r1 *= max_rgb_diff / max_rgb
			rgb_r2 *= max_rgb_diff / max_rgb
			rgb_r3 *= max_rgb_diff / max_rgb
		rgb_r1 += min_rgb
		rgb_r2 += min_rgb
		rgb_r3 += min_rgb
		
		# cut off
		# This is done after the previous step because negative numbers correspond to
		# increased saturation.
		rgb_r1 = max(rgb_r1,0)
		rgb_r2 = max(rgb_r2,0)
		rgb_r3 = max(rgb_r3,0)
		
		# output luminance
		if (max_rgb > 0):
			lum_out = rgb_r1*0.2126729 + rgb_r2*0.7151522 + rgb_r3*0.0721750
			lum_diff = lum_in / lum_out
			rgb_r1 *= lum_diff
			rgb_r2 *= lum_diff
			rgb_r3 *= lum_diff
		
		# redo gamma
		rgb_r1 = c.gamma(rgb_r1)
		rgb_r2 = c.gamma(rgb_r2)
		rgb_r3 = c.gamma(rgb_r3)

		# XYZ luminosity
		# source: http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
		# The values have to be both clipped and rounded off or white parts of the image
		# get mangled.
		img_rgb[x][y][0] = min(255, round(rgb_r1 * 255))
		img_rgb[x][y][1] = min(255, round(rgb_r2 * 255))
		img_rgb[x][y][2] = min(255, round(rgb_r3 * 255))
		if (canary):
			print((rgb_r1*255,rgb_r2*255,rgb_r3*255))
			print(img_rgb[x][y])

# convert back to BGR
img_result = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2BGR)

# print execution time
print("%s seconds" % (time.time() - start_time))

# display result
# Now actually just displays it. Saving is optional.
cv2.imshow("output", img_result)
cv2.waitKey(0)
cv2.destroyAllWindows()
