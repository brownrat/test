"""
10/10/2025 -- bringing the 24 May 2024 version back as a separate file so OpenCV doesn't have to share
space with matplotlib.

This has two modes depending on --s2tmode:

* chromaticity (c, default): similar to Color Vision Simulator and gives comparable results when
provided with the same images and parameters. An equation is found for transforming the source
chromaticity for wavelengths 300-700 into the target chromaticity, and this is used on a
fully saturated version of the pixel's color. The final saturation differs if the transformed
color has more or less than 100% saturation (smallest RGB value >0 = less, <0 = more).
* sensitivity (s): similar to micaToolbox. Each channel as a whole is transformed, and the
fit is based on the sensitivity curves for source and target rather than chromaticity.
Before finding the equation, all the curves are multiplied by 2147483647.

For how the transformations work, see here:

* https://www.empiricalimaging.com/knowledge-base/spectral-sensitivity-based-cone-catch-models/
* https://doi.org/10.1111/2041-210X.12439

Interactions between channels can be set to either none (default), two-way or three-way.

Adding a UV photo can be done with --uvimage. In fact the code always behaves as though there's a
fourth channel, so when a UV image isn't given it tends to produce warnings about "Covariance of the
parameters could not be estimated". Still works though.

How luminance is handled depends on the argument --luminance:

* constant (default for chromaticity): held constant based on the XYZ Y value
* raw (default for sensitivity): transformed RGB values used as is
* adjust: adjusted to the target's luminance vision (may be inaccurate)

Since there's no built-in recognition of whether something is an "L cone" etc., the target's
receptors are assigned as r1=R, r2=G, r3=B, so you can choose which is which just by changing the
order. If there are more than three, the fourth and up are ignored. Dichromatic images are assigned
based on the peak sensitivity of receptor 2: <490 means r1=RG, r2=B; >=490 means r1=R, r2=GB. This
follows the usual convention for representations of human color deficiency where
protanopia/deuteranopia is yellow/blue and tritanopia is red/cyan.
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
# OpenCV reads images in BGR, so we convert to RGB
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
		rgb_r = img_rgb[x][y][0]
		rgb_g = img_rgb[x][y][1]
		rgb_b = img_rgb[x][y][2]
		rgb_r *= 1/255
		rgb_g *= 1/255
		rgb_b *= 1/255
		
		rgb_uv = 0
		if (args.uvimage != "none" and y < uv_img.shape[1]):
			rgb_uv = uv_img[x][y][0]
			rgb_uv *= 1/255

		# gamma decode
		rgb_r = c.gamma_decode(rgb_r)
		rgb_g = c.gamma_decode(rgb_g)
		rgb_b = c.gamma_decode(rgb_b)
		rgb_uv = c.gamma_decode(rgb_uv)

		# input luminance
		lum_in = rgb_r*c.lum_r + rgb_g*c.lum_g + rgb_b*c.lum_b
		lum_in1 = (rgb_r*c.lum_r1 + rgb_g*c.lum_g1 + rgb_b*c.lum_b1 + rgb_uv*c.lum_uv) / c.lum_max

		if (args.s2tmode == "c"):
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
		
		# 0-1 -> 0-2147483647
		# Removed. It seems like the RGB values should be on the same scale as
		# the curves we fit, but scaling them all the way up just causes white
		# areas to turn into a pure saturated color while everything else
		# is the same. I don't think that's supposed to happen.
		#if (args.s2tmode == "s"):
		#	rgb_r *= 2147483647
		#	rgb_g *= 2147483647
		#	rgb_b *= 2147483647
		#	rgb_uv *= 2147483647

		# combine
		if (args.interactions == 1): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_uv]
		elif (args.interactions == 2): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_uv, rgb_r*rgb_g, rgb_r*rgb_b, rgb_r*rgb_uv, rgb_g*rgb_b, rgb_g*rgb_uv, rgb_b*rgb_uv]
		elif (args.interactions == 3): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_uv, rgb_r*rgb_g, rgb_r*rgb_b, rgb_r*rgb_uv, rgb_g*rgb_b, rgb_g*rgb_uv, rgb_b*rgb_uv, rgb_r*rgb_g*rgb_b, rgb_r*rgb_g*rgb_uv, rgb_r*rgb_b*rgb_uv, rgb_g*rgb_b*rgb_uv, rgb_r**2, rgb_g**2, rgb_b**2, rgb_uv**2]
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

		max_rgb = max(rgb_r1,rgb_r2,rgb_r3)
		if (args.s2tmode == "c"):
			# desaturate
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
			if (args.luminance == "constant"
			or (args.luminance == "default" and args.s2tmode == "c")):
				lum_diff = lum_in / lum_out
			elif (args.luminance == "adjust"):
				lum_diff = lum_in1 / lum_out
			elif (args.luminance == "raw"
			or (args.luminance == "default" and args.s2tmode == "s")):
				lum_diff = 1
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
