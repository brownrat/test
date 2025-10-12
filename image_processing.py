"""
10/10/2025 -- bringing the 24 May 2024 version back as a separate file so OpenCV doesn't have to share
space with matplotlib.

I abandoned the approach of finding a "dominant wavelength" because I couldn't find a satisfactory way
to do that. Human color appearance is not derived from LMS in the way CVS seems to assume.
Instead I decided to skip that step and directly transform the source chromaticity into the target
chromaticity by fitting the source sensitivity curves to the target, somewhat similar to a method
available in micaToolbox (https://www.empiricalimaging.com/knowledge-base/spectral-sensitivity-based-cone-catch-models/).
If the "source" and "target" are specified as described (CIE 10-deg fundamentals -> standard templates
+ human ocular media), the output closely resembles CVS images except when the target includes an M
cone with a peak wavelength less than the human M cone (about 540 nm). This prevents the red channel
from being mapped on to it at all. CVS has a similar flaw where dichromatic simulations with an L cone
close to the human L cone cause the green channel to be removed and blue-green regions to turn pure
blue, which my code also replicates.

The RGB values are scaled to the desired intensity using HSV. This is obviously not as good as
converting to a color space designed to be more perceptually uniform like Lab or the LCH spaces,
but any kind of conversion between sRGB and another color space seems to distort the hue and/or
saturation.
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

red_r1 = c.red_r1
green_r1 = c.green_r1
blue_r1 = c.blue_r1

red_r2 = c.red_r2
green_r2 = c.green_r2
blue_r2 = c.blue_r2

red_r3 = c.red_r3
green_r3 = c.green_r3
blue_r3 = c.blue_r3

# gamma correction
def linear(v):
            if v <= 0.04045: return v / 12.92
            else: return math.pow((v + 0.055) / 1.055, 2.4)
            #return math.pow(v, 2.2)

# image processing
imagename = args.image
# Convert image from BGR to HLS. We use BGR because this is how OpenCV reads images.
# If we use RGB, the output appears normal with settings approximating human vision,
# but shifting the cones produces the opposite of the expected result.
img = cv2.imread(imagename)
img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
img_hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)

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
		
		# undo gamma
		rgb_r = linear(rgb_r)
		rgb_g = linear(rgb_g)
		rgb_b = linear(rgb_b)
		
		# combine
		rgb_r1 = red_r1*rgb_r + green_r1*rgb_g + blue_r1*rgb_b
		rgb_r2 = red_r2*rgb_r + green_r2*rgb_g + blue_r2*rgb_b
		rgb_r3 = red_r3*rgb_r + green_r3*rgb_g + blue_r3*rgb_b
		
		# redo gamma
		rgb_r1 = c.gamma(rgb_r1)
		rgb_r2 = c.gamma(rgb_r2)
		rgb_r3 = c.gamma(rgb_r3)
		
		# remove junk channels for lower dimensions
		if (c.dimension == 2):
			rgb_r3 = rgb_r2
			# standard tritanopia colors
			if (args.receptors[1] < 490): rgb_r2 = rgb_r1
		if (c.dimension == 1):
			rgb_r2 = rgb_r1
			rgb_r3 = rgb_r1

		# HSV lightness. This makes the chromaticity conversion redundant.
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
