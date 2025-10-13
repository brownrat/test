"""
10/10/2025 -- bringing the 24 May 2024 version back as a separate file so OpenCV doesn't have to share
space with matplotlib.

I abandoned the approach of finding a "dominant wavelength" because I couldn't find a satisfactory way
to do that. Human color appearance is not derived from LMS in the way CVS seems to assume.
Instead I decided to skip that step and directly transform the source chromaticity into the target
chromaticity by fitting the source sensitivity curves to the target, somewhat similar to a method
available in micaToolbox (https://www.empiricalimaging.com/knowledge-base/spectral-sensitivity-based-cone-catch-models/).

The RGB values are scaled to the desired intensity using HSV. This is obviously not as good as
converting to a color space designed to be more perceptually uniform like Lab or the LCH spaces,
but any kind of conversion between sRGB and another color space seems to distort the hue and/or
saturation.

micaToolbox uses not just a linear sum of the RGB channels but also up to 3-way interactions between
them. I tried including 2-way interactions. The fitted curves look great on a plot, but the image
looks like garbage. Using the monochromatic colors from optimal.py as a color chart, it seems the fit is
only valid for monochromatic or fully saturated colors whereas anything closer to white is increasingly
wrong. It's also strongly affected by the scale of both the input curves and input RGB values because
they're multiplied. I get the closest to reasonable results by converting the curves and input to
chromaticity and correcting the white balance of the output, but the original paper (https://doi.org/10.1111/2041-210X.12439)
doesn't mention doing either of these. I also find the results are less likely to turn bright magenta
and such if the lower bound is set to 0.
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
img_hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)

# white balance
gray = [1/3,1/3,1/3,1/9,1/9,1/9]
gray_r1 = sum(gray * c.s2t_r1)
gray_r2 = sum(gray * c.s2t_r2)
gray_r3 = sum(gray * c.s2t_r3)

# transform hues for each pixel
for x in range(img.shape[0]):
	for y in range(img.shape[1]):
		# report current pixel to determine whether there's an infinite loop
		#print("pixel: " + str(x) + ", " + str(y))

		# get pixel
		pixel = img_rgb[x][y]
		
		# 0-255 -> 0-1
		rgb_r = pixel[0] / 255
		rgb_g = pixel[1] / 255
		rgb_b = pixel[2] / 255
		#print(rgb_r)
		#print(type(rgb_r))
		
		# undo gamma
		rgb_r = c.linearize(rgb_r)
		rgb_g = c.linearize(rgb_g)
		rgb_b = c.linearize(rgb_b)
		#print(type(rgb_r))
		
		# chromaticity
		total = rgb_r + rgb_g + rgb_b
		if (total > 0):
			rgb_r /= total
			rgb_g /= total
			rgb_b /= total
		#print(type(rgb_r))
		
		# combine
		rgb_terms = [rgb_r, rgb_g, rgb_b, rgb_r*rgb_g, rgb_r*rgb_b, rgb_g*rgb_b]
		##print(rgb_terms)
		#print(rgb_terms * c.s2t_r1)
		#print(rgb_terms * c.s2t_r2)
		#print(rgb_terms * c.s2t_r3)
		rgb_r1 = sum(rgb_terms * c.s2t_r1)
		rgb_r2 = sum(rgb_terms * c.s2t_r2)
		rgb_r3 = sum(rgb_terms * c.s2t_r3)
		#print((rgb_r1,rgb_r2,rgb_r3))
		#print(type(rgb_r1))
		
		if (args.whitebalance == "on"):
			# white balance
			if (gray_r1 > 0): rgb_r1 /= gray_r1
			if (gray_r2 > 0): rgb_r2 /= gray_r2
			if (gray_r3 > 0): rgb_r3 /= gray_r3
		
		# cut off
		#rgb_r1 = min(rgb_r1,1)
		#rgb_r2 = min(rgb_r2,1)
		#rgb_r3 = min(rgb_r3,1)
		rgb_r1 = max(rgb_r1,0)
		rgb_r2 = max(rgb_r2,0)
		rgb_r3 = max(rgb_r3,0)
		
		# redo gamma
		rgb_r1 = c.gamma(rgb_r1)
		rgb_r2 = c.gamma(rgb_r2)
		rgb_r3 = c.gamma(rgb_r3)
		#print(type(rgb_r1))
		
		# remove junk channels for lower dimensions
		if (c.dimension == 2):
			rgb_r3 = rgb_r2
			# standard tritanopia colors
			if (args.receptors[1] < 490): rgb_r2 = rgb_r1
		if (c.dimension == 1):
			rgb_r2 = rgb_r1
			rgb_r3 = rgb_r1

		# HSV lightness
		c_max = max(rgb_r1,rgb_r2,rgb_r3)
		
		if (c_max > 0):
			img_rgb[x][y][0] = rgb_r1*img_hsv[x][y][2]/c_max
			img_rgb[x][y][1] = rgb_r2*img_hsv[x][y][2]/c_max
			img_rgb[x][y][2] = rgb_r3*img_hsv[x][y][2]/c_max
		else: img_rgb[x][y] = [0,0,0]

# convert back to BGR
img_result = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2BGR)

# print execution time
print("%s seconds" % (time.time() - start_time))

# display result
# Now actually just displays it. Saving is optional.
cv2.imshow("output", img_result)
cv2.waitKey(0)
cv2.destroyAllWindows()
