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

Adding a UV photo can be done with --uvimage. Only the red channel is used, so if you want the other
channels they should be merged/equalized first. If a UV image is not specified, white/gray parts of
the non-UV image will stay that way even if the target has a UV receptor.

By default, luminance is held constant by comparing the original XYZ Y value with the final value.
The argument --changelum allows for changing the luminance to match the target in an analogous way
using their sensitivity to typical sRGB primaries and (if applicable) the selected UV channel. This
feature is more experimental as I'm not aware of an "official" way to do this (CVS and micaToolbox
don't bother).
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
uv = c.sr4

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
if (args.uvimage != "none"):
	uvimagename = args.uvimage
	uv_img = cv2.imread(uvimagename)

# transform hues for each pixel
for x in range(img.shape[0]):
	for y in range(img.shape[1]):
		# report current pixel (for testing)
		#print("pixel: " + str(x) + ", " + str(y))

		# 0-255 -> 0-1
		rgb_r = img_rgb[x][y][0] / 255
		rgb_g = img_rgb[x][y][1] / 255
		rgb_b = img_rgb[x][y][2] / 255
		
		rgb_uv = 0
		if (args.uvimage != "none" and y < uv_img.shape[1]):
			rgb_uv = uv_img[x][y][0] / 255

		# gamma decode
		rgb_r = c.gamma_decode(rgb_r)
		rgb_g = c.gamma_decode(rgb_g)
		rgb_b = c.gamma_decode(rgb_b)
		rgb_uv = c.gamma_decode(rgb_uv)

		# default input luminance
		if (args.changelum):
			lum_in = (rgb_r*c.lum_r1 + rgb_g*c.lum_g1 + rgb_b*c.lum_b1 + rgb_uv*c.lum_uv) / c.lum_max
		else:
			lum_in = rgb_r*c.lum_r + rgb_g*c.lum_g + rgb_b*c.lum_b

		# saturate
		if (args.uvimage == "none"):
			min_rgb = min(rgb_r,rgb_g,rgb_b)
			rgb_r -= min_rgb
			rgb_g -= min_rgb
			rgb_b -= min_rgb
			max_rgb_diff = max(rgb_r,rgb_g,rgb_b)
		else:
			min_rgb = min(rgb_r,rgb_g,rgb_b,rgb_uv)
			rgb_r -= min_rgb
			rgb_g -= min_rgb
			rgb_b -= min_rgb
			rgb_uv -= min_rgb
			max_rgb_diff = max(rgb_r,rgb_g,rgb_b,rgb_uv)
		
		# chromaticity
		total = rgb_r + rgb_g + rgb_b + rgb_uv
		if (total > 0):
			rgb_r /= total
			rgb_g /= total
			rgb_b /= total
			rgb_uv /= total

		# combine
		rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_uv, rgb_r*rgb_g, rgb_r*rgb_b, rgb_r*rgb_uv, rgb_g*rgb_b, rgb_g*rgb_uv, rgb_b*rgb_uv]
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
			lum_out = rgb_r1*c.lum_r + rgb_r2*c.lum_g + rgb_r3*c.lum_b
			lum_diff = lum_in / lum_out
			if (lum_diff > 0):
				rgb_r1 *= lum_diff
				rgb_r2 *= lum_diff
				rgb_r3 *= lum_diff

		# gamma encode
		rgb_r1 = c.gamma_encode(rgb_r1)
		rgb_r2 = c.gamma_encode(rgb_r2)
		rgb_r3 = c.gamma_encode(rgb_r3)

		# 0-1 -> 0-255
		# The values have to be both clipped and rounded off or white parts of the image
		# get mangled.
		img_rgb[x][y][0] = min(255, round(rgb_r1 * 255))
		img_rgb[x][y][1] = min(255, round(rgb_r2 * 255))
		img_rgb[x][y][2] = min(255, round(rgb_r3 * 255))

# convert back to BGR
img_result = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2BGR)

# print execution time
print("%s seconds" % (time.time() - start_time))

# display result
# Now actually just displays it. Saving is optional.
cv2.imshow("output", img_result)
cv2.waitKey(0)
cv2.destroyAllWindows()
