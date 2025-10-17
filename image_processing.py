"""
10/10/2025 -- bringing the 24 May 2024 version back as a separate file so OpenCV doesn't have to
share space with matplotlib.

This has two modes depending on --s2tmode:

* chromaticity (c, default): similar to Color Vision Simulator and gives comparable results when
provided with the same images and parameters. An equation is found for transforming the source
chromaticity for wavelengths 300-700 into the target chromaticity, and this is used on a
fully saturated version of the pixel's color. The final saturation differs if the transformed
color has more or less than 100% saturation (smallest RGB value >0 = less, <0 = more).
* sensitivity (s): similar to micaToolbox. Each channel as a whole is transformed, and the
fit is based on the sensitivity curves for source and target rather than chromaticity.
Before finding the equation, all the curves are multiplied by 65536.

For how the transformations work, see here:

* https://www.empiricalimaging.com/knowledge-base/spectral-sensitivity-based-cone-catch-models/
* https://doi.org/10.1111/2041-210X.12439

Interactions between channels can be set to either none (default), two-way or three-way.

Adding a UV photo can be done with --uvimage.

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

I fixed the OpenCV/matplotlib problem by installing opencv-headless and going back to saving
the output image with imwrite() instead of displaying it with imshow(). I think I'll keep this
as a separate file to reduce clutter.
"""

import math
import numpy as np
import time
import argparse
import cv2
import central as c
import matplotlib.pyplot as plt
import scipy
args = c.args

# execution time
start_time = time.time()

# set source receptors
srlist = []
if (args.source[0] == "cie2"): srlist = c.cie2
# CIE 10-deg fundamentals
elif (args.source[0] == "cie10"): srlist = c.cie10
# UVS bird
elif (args.source[0] == "uvbird"): srlist = c.uvbird
# VS bird
elif (args.source[0] == "vbird"): srlist = c.vbird
else:
	for i in range(len(args.source)):
		srlist.append(csv2spec(args.source[i]))

sr1 = srlist[0]
sr2 = srlist[1]
sr3 = srlist[2]
sr4 = np.zeros(401)
if ((len(srlist)) == 4): sr4 = srlist[3]

# luminance adjustment
# sRGB -> XYZ Y: http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
lum_r = 0.2126729
lum_g = 0.7151522
lum_b = 0.0721750

# contribution of each channel to the target's luminance vision
lum_r1 = sum(sr1 * c.luminosity)
lum_g1 = sum(sr2 * c.luminosity)
lum_b1 = sum(sr3 * c.luminosity)
lum_uv = sum(sr4 * c.luminosity)
lum_max = lum_r1 + lum_g1 + lum_b1 + lum_uv

# adapt to illuminant
sr1 *= c.wp_1nm
sr2 *= c.wp_1nm
sr3 *= c.wp_1nm
sr4 *= c.wp_1nm
wpsr1 = sum(sr1)
wpsr2 = sum(sr2)
wpsr3 = sum(sr3)
wpsr4 = sum(sr4)
sr1 /= wpsr1
sr2 /= wpsr2
sr3 /= wpsr3
if (wpsr4 > 0): sr4 /= wpsr4

# declare variables
tr1 = np.zeros(401)
tr2 = np.zeros(401)
tr3 = np.zeros(401)
wptr1 = 1
wptr2 = 1
wptr3 = 1

# use CMFs instead of raw receptors
# No von Kries transform because they're already adapted.
if (args.fcmode == 'cmf'):
	matrix_cmf = cmf()
	tr1 = matrix_cmf[0]
	if (dimension > 1): tr2 = matrix_cmf[1]
	if (dimension > 2): tr3 = matrix_cmf[2]
	
	# cut off negative numbers -- they mess with chromaticity and represent something
	# that can't be displayed anyway
	for i in range(401):
		if (tr1[i] < 0): tr1[i] = 0
		if (tr2[i] < 0): tr2[i] = 0
		if (tr3[i] < 0): tr3[i] = 0
else:
	tr1 = c.r1_1nm * c.media_1nm * c.wp_1nm
	tr2 = c.r2_1nm * c.media_1nm * c.wp_1nm
	tr3 = c.r3_1nm * c.media_1nm * c.wp_1nm

	# von Kries transform
	wptr1 = sum(tr1)
	wptr2 = sum(tr2)
	wptr3 = sum(tr3)
	tr1 /= wptr1
	if (wptr2 > 0): tr2 /= wptr2
	if (wptr3 > 0): tr3 /= wptr3

if (args.s2tmode == "c"):
	# convert to chromaticity
	sr1_c = np.zeros(401)
	sr2_c = np.zeros(401)
	sr3_c = np.zeros(401)
	sr4_c = np.zeros(401)
	tr1_c = np.zeros(401)
	tr2_c = np.zeros(401)
	tr3_c = np.zeros(401)
	for i in range(401):
		stotal = sr1[i] + sr2[i] + sr3[i] + sr4[i]
		ttotal = tr1[i] + tr2[i] + tr3[i]
		if (stotal > 0):
			sr1_c[i] = sr1[i] / stotal
			sr2_c[i] = sr2[i] / stotal
			sr3_c[i] = sr3[i] / stotal
			sr4_c[i] = sr4[i] / stotal
		if (ttotal > 0):
			tr1_c[i] = tr1[i] / ttotal
			tr2_c[i] = tr2[i] / ttotal
			tr3_c[i] = tr3[i] / ttotal
elif (args.s2tmode == "s"):
	# 0 - 65535
	smax = max(*sr1,*sr2,*sr3,*sr4)
	tmax = max(*tr1,*tr2,*tr3)
	sr1_c = sr1 * 65535 / smax
	sr2_c = sr2 * 65535 / smax
	sr3_c = sr3 * 65535 / smax
	sr4_c = sr4 * 65535 / smax
	tr1_c = tr1 * 65535 / smax
	tr2_c = tr2 * 65535 / smax
	tr3_c = tr3 * 65535 / smax

# find best fit to target
if (args.uvimage == "none"):
	def source2target(xdata, r, g, b, rg=0, rb=0, gb=0, rgb=0, r2=0, g2=0, b2=0):
		first_order = r*sr1_c + g*sr2_c + b*sr3_c
		second_order = rg*(sr1_c*sr2_c) + rb*(sr1_c*sr3_c) + gb*(sr2_c*sr3_c)
		third_order = rgb*(sr1_c*sr2_c*sr3_c) + r2*(sr1_c**2) + g2*(sr2_c**2) + b2*(sr3_c**2)
		if (args.interactions == 3): ydata = first_order + second_order + third_order
		elif (args.interactions == 2): ydata = first_order + second_order
		else: ydata = first_order
		return(ydata)

	if (args.interactions == 1): p0 = [1,1,1]
	elif (args.interactions == 2): p0 = [1,1,1,0,0,0]
	elif (args.interactions == 3): p0 = [1,1,1,0,0,0,0,0,0,0]
else:
	def source2target(xdata, r, g, b, u, rg=0, rb=0, ru=0, gb=0, gu=0, bu=0, rgb=0, rgu=0, rbu=0, gbu=0, r2=0, g2=0, b2=0, u2=0):
		first_order = r*sr1_c + g*sr2_c + b*sr3_c + u*sr4_c
		second_order = rg*(sr1_c*sr2_c) + rb*(sr1_c*sr3_c) + ru*(sr1_c*sr4_c) + gb*(sr2_c*sr3_c) + gu*(sr2_c*sr4_c) + bu*(sr3_c*sr4_c)
		third_order = rgb*(sr1_c*sr2_c*sr3_c) + rgu*(sr1_c*sr2_c*sr4_c) + rbu*(sr1_c*sr3_c*sr4_c) + gbu*(sr2_c*sr3_c*sr4_c) + r2*(sr1_c**2) + g2*(sr2_c**2) + b2*(sr3_c**2) + u2*(sr4_c**2)
		if (args.interactions == 3): ydata = first_order + second_order + third_order
		elif (args.interactions == 2): ydata = first_order + second_order
		else: ydata = first_order
		return(ydata)

	if (args.interactions == 1): p0 = [1,1,1,1]
	elif (args.interactions == 2): p0 = [1,1,1,1,0,0,0,0,0,0]
	elif (args.interactions == 3): p0 = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

popt, pcov = scipy.optimize.curve_fit(source2target, c.x_1nm, tr1_c, p0=p0, bounds=[args.bound, np.inf])
if (args.verbose): print(popt)
s2t_r1 = popt
popt, pcov = scipy.optimize.curve_fit(source2target, c.x_1nm, tr2_c, p0=p0, bounds=[args.bound, np.inf])
if (args.verbose): print(popt)
s2t_r2 = popt
popt, pcov = scipy.optimize.curve_fit(source2target, c.x_1nm, tr3_c, p0=p0, bounds=[args.bound, np.inf])
if (args.verbose): print(popt)
s2t_r3 = popt

if (args.s2tplot):
	plt.plot(c.x_1nm, tr1_c, color='k')
	plt.plot(c.x_1nm, source2target(c.x_1nm, *s2t_r1), '--', color='k')
	if (c.dimension > 1):
		plt.plot(c.x_1nm, tr2_c, color='k')
		plt.plot(c.x_1nm, source2target(c.x_1nm, *s2t_r2), '--', color='k')
	if (c.dimension > 2):
		plt.plot(c.x_1nm, tr3_c, color='k')
		plt.plot(c.x_1nm, source2target(c.x_1nm, *s2t_r3), '--', color='k')
	plt.show()

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
		lum_in = rgb_r*lum_r + rgb_g*lum_g + rgb_b*lum_b
		lum_in1 = (rgb_r*lum_r1 + rgb_g*lum_g1 + rgb_b*lum_b1 + rgb_uv*lum_uv) / lum_max

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
		if (args.uvimage == "none"):
			if (args.interactions == 1): rgb_terms = [rgb_r, rgb_g, rgb_b]
			elif (args.interactions == 2): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_r*rgb_g, rgb_r*rgb_b, rgb_g*rgb_b]
			elif (args.interactions == 3): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_r*rgb_g, rgb_r*rgb_b, rgb_g*rgb_b, rgb_r*rgb_g*rgb_b, rgb_r**2, rgb_g**2, rgb_b**2]
		else:
			if (args.interactions == 1): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_uv]
			elif (args.interactions == 2): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_uv, rgb_r*rgb_g, rgb_r*rgb_b, rgb_r*rgb_uv, rgb_g*rgb_b, rgb_g*rgb_uv, rgb_b*rgb_uv]
			elif (args.interactions == 3): rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_uv, rgb_r*rgb_g, rgb_r*rgb_b, rgb_r*rgb_uv, rgb_g*rgb_b, rgb_g*rgb_uv, rgb_b*rgb_uv, rgb_r*rgb_g*rgb_b, rgb_r*rgb_g*rgb_uv, rgb_r*rgb_b*rgb_uv, rgb_g*rgb_b*rgb_uv, rgb_r**2, rgb_g**2, rgb_b**2, rgb_uv**2]

		rgb_r1 = sum(rgb_terms * s2t_r1)
		rgb_r2 = sum(rgb_terms * s2t_r2)
		rgb_r3 = sum(rgb_terms * s2t_r3)
		
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
			lum_out = rgb_r1*lum_r + rgb_r2*lum_g + rgb_r3*lum_b
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

# save result at specified file name
# Can't use imshow() anymore because we require opencv-headless.
cv2.imwrite(args.outimage, img_result)
