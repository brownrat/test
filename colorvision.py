# This is an attempt to create something similar to Color Vision Simulator (Melin et al. 2013).
# Unlike CVS, it only alters the hue and doesn't account for saturation. I don't know (yet)
# how to determine the distance of a color from the white point of a color space.
# Like CVS, it also doesn't account for differences in lightness/intensity. This could be
# done by finding the global photopic sensitivity curve and dividing by the human version. It
# would also be nice to use LCh instead of HSL so the lightness result would make sense.

import cv2
import sys
import colorsys
import math
import numpy as np

# image
imagename = sys.argv[1]

# cone peak sensitivities
lcone = int(sys.argv[2])
mcone = int(sys.argv[3])
scone = int(sys.argv[4])

# method to use
version = 2
if (len(sys.argv) > 5):
	version = int(sys.argv[5])

# convert hue to closest wavelength match (linear)
#def hue_to_wavelength(hue):
#	# spectral colors
#	if (hue <= 240):
#		w = 625 - hue * 175/240 # using wavelengths 450-625
#		return int(w) # round off
#	# extra-spectral colors (purple/magenta/pink)
#	else:
#		r = hue/120 - 2
#		b = -hue/120 + 3
#		return [r, b]

# non-linear (piecewise) version based on sRGB primaries
def hue_to_wavelength(hue):
	#print(hue)
	# red-yellow: 0-59 -> 620-581
	if (0 <= hue < 60):
		w = 620 - hue * (620 - 581)/60
	# yellow-green: 60-119 -> 580-551
	elif (60 <= hue < 120):
		w = 580 - (hue - 60) * (580 - 551)/60
	# green-cyan: 120-179 -> 550-501
	elif (120 <= hue < 180):
		w = 550 - (hue - 120) * (550 - 501)/60
	# cyan-blue: 180-239 -> 500-471
	elif (180 <= hue < 240):
		w = 500 - (hue - 180) * (500 - 471)/60
	# blue-violet: 240-269 -> 470-400
	elif (240 <= hue < 270):
		w = 470 - (hue - 240) * (470 - 420)/30
	# purple-magenta: 270-360
	else:
		r = hue/120 - 2
		b = -hue/120 + 3
		return [r, b]
#	print(w)
	return w

# find sensitivity of a given cone type to a wavelength
def wl_to_s(wl, peak):
	value = (math.exp(69.7*(0.8795+0.0459*math.exp(-(peak-300)**2/11940)-(peak/wl)))+math.exp(28*(0.922-(peak/wl)))+math.exp(-14.9*(1.104-(peak/wl)))+0.674)**-1 + 0.26*math.exp(-((wl-(189+0.315*peak))/(-40.5+0.195*peak))**2)
	return value

# convert wavelength to RGB values for target color space

# simple version: use unaltered cone fundamentals
# This gives something in the right ballpark but is less and less accurate the farther away
# the primaries are from the cone peaks. To produce a result close to the original (i.e. human
# vision), you have to use the values 620, 550 and 470 that I chose for the reference sRGB
# primaries.
def wl_to_rgb_0(wl, l, m, s):
	red = wl_to_s(wl, l)
	green = wl_to_s(wl, m)
	blue = wl_to_s(wl, s)
	return [red, green, blue]

# complex version: use color matching functions
# In theory this should be the most accurate, but it does something weird where red and green
# seem to be switched around. It's also extremely slow due to iterating over two 250x3 matrices
# for every pixel.
def wl_to_rgb_1(wl, l, m, s):
	# primaries
	red = l + 40 - (535 - m) # roughly shift red primary depending on M and L peaks: less distance -> "red" is farther away from L peak
	green = m
	blue = s

	# sensitivity of cones to primaries
	matrix_a = np.array([
		[wl_to_s(red, l), wl_to_s(green, l), wl_to_s(blue, l)],
		[wl_to_s(red, m), wl_to_s(green, m), wl_to_s(blue, m)],
		[wl_to_s(red, s), wl_to_s(green, s), wl_to_s(blue, s)]
	])
	print(matrix_a)

	# sensitivity of cones to each wavelength
	matrix_c = np.empty((251, 3))
	for i in range(0, 250):
		matrix_c[i][0] = wl_to_s(i + 400, l)
		matrix_c[i][1] = wl_to_s(i + 400, m)
		matrix_c[i][2] = wl_to_s(i + 400, s)

	# initial color match matrix
	matrix_m = np.matmul(matrix_c, np.linalg.inv(matrix_a)) # A x M = C, so M = C x A^-1 (if I remember right)

	# white matching

	# initialize (do you need to do this?)
	rw = 0.0
	gw = 0.0
	bw = 0.0

	# there's probably a better way to find these sums, but whatever
	for i in range (0, 250):
		rw += matrix_m[i][0]
		gw += matrix_m[i][1]
		bw += matrix_m[i][2]
	#print(rw)

	# final color match matrix
	matrix_cmf = np.empty((251, 3))
	for i in range(0, 250):
		matrix_cmf[i][0] = matrix_m[i][0] / rw
		matrix_cmf[i][1] = matrix_m[i][1] / gw
		matrix_cmf[i][2] = matrix_m[i][2] / bw

	# find match given a wavelength
	#rgb = [matrix_cmf[wl - 400][0], matrix_cmf[wl - 400][1], matrix_cmf[wl - 400][2]]
	#print(rgb)
	rgb = [matrix_m[wl - 400][0], matrix_m[wl - 400][1], matrix_m[wl - 400][2]]
	return rgb

# less accurate version: cone response ratios
# we find the wavelengths that produce ratios matching those for the sRGB primaries and
# secondaries in human vision, then use those to convert the wavelength to a hue as the
# reverse of the piecewise hue_to_wavelength. To account for differences in the range of
# ratios depending on LM overlap, we find the ratios at the L and M peaks for human vision
# and the target and use them to scale the "red", "yellow" and "green" ratios to the
# target's range.
# When used with the images from studies that used CVS (see also "Experimental evidence that
# primate trichromacy is well suited for detecting primate social colour signals", Hiramatsu
# et al. 2017) this version gives results very similar to the actual CVS output shown in these
# papers besides not accounting for saturation. Melin et al. says CVS uses a table of
# "relative sensitivities", so my method is probably roughly what it does besides using
# fewer points and more linear interpolation.
# Note that I've used the Stockton and Sharpe values of 440-540-565 instead of
# 420-530-560 as used in the studies, so using these values produces a red-shifted image.
# (Using the Stockton and Sharpe values produces a nearly 1:1 transformation.)
# A cone spacing of less than 10 nm usually produces a math overflow error. Fortunately there
# are not many species with this kind of arrangement (but there are some).

def find_primaries(l, m, s):
	# find ratios to match
	rr = wl_to_s(620, 565) / wl_to_s(620, 540)
	ry = wl_to_s(580, 565) / wl_to_s(580, 540)
	rg = wl_to_s(550, 565) / wl_to_s(550, 540)
	rc = wl_to_s(500, 540) / wl_to_s(500, 440)
	rb = wl_to_s(470, 540) / wl_to_s(470, 440)
	rv = wl_to_s(420, 565) / wl_to_s(420, 540)
	
	# scale ratios
	rr1 = wl_to_s(565, 565) / wl_to_s(540, 565) # LM ratio at L peak
	rg1 = wl_to_s(540, 540) / wl_to_s(565, 540) # ML ratio at M peak
	rr2 = wl_to_s(l, l) / wl_to_s(m, l) # LM ratio at target L peak
	rg2 = wl_to_s(m, m) / wl_to_s(l, m) # ML ratio at target M peak
	
	# find "red" L primary
	red = 650
	ratio = wl_to_s(red, l) / wl_to_s(red, m)
	while (ratio > rr * rr2 / rr1):
		red -= 1
		ratio = wl_to_s(red, l) / wl_to_s(red, m)
	
	# find "yellow" LM secondary
	yellow = red
	ratio = wl_to_s(yellow, l) / wl_to_s(yellow, m)
	while (ratio > ry * rr2 / rr1):
		yellow -= 1
		ratio = wl_to_s(yellow, l) / wl_to_s(yellow, m)
	
	# find "green" M primary
	green = yellow
	ratio = wl_to_s(green, l) / wl_to_s(green, m)
	while (ratio > 1): # count back to the crossover point
		green -= 1
		ratio = wl_to_s(green, l) / wl_to_s(green, m)
	while (ratio > rg * rg2 / rg1): # continue on the other side
		green -= 1
		ratio = wl_to_s(green, l) / wl_to_s(green, m)
	
	# find "cyan" MS secondary
	cyan = m
	ratio = wl_to_s(cyan, m) / wl_to_s(cyan, s)
	while (ratio > rc):
		cyan -= 1
		ratio = wl_to_s(cyan, m) / wl_to_s(cyan, s)
	
	# find "blue" S primary
	blue = cyan
	ratio = wl_to_s(blue, m) / wl_to_s(blue, s)
	while (ratio > rb):
		blue -= 1
		ratio = wl_to_s(blue, m) / wl_to_s(blue, s)
	
	# find "violet" LS secondary
	violet = blue
	ratio = wl_to_s(violet, l) / wl_to_s(violet, m)
	while (ratio < rv):
		violet -= 1
		ratio = wl_to_s(violet, l) / wl_to_s(violet, m)
	print([red, yellow, green, cyan, blue, violet])
	return [red, yellow, green, cyan, blue, violet]
	
# second less accurate version: nearby cone sensitivities
# Instead of using ratios, we find the points where a type of cone "falls off" and gives way
# to another color. This blue-shifts the "yellow point" and "green point" more significantly
# when L and M approach each other and fixes the "cyan point" at the S cutoff, which may be
# more realistic. However, the red-green difference is overestimated relative to actual CVS
# output.

def find_primaries_1(l, m, s):
	# find sensitivities to match
	sr = wl_to_s(620, 540)
	sy = wl_to_s(580, 540)
	sg = wl_to_s(550, 565)
	sc = wl_to_s(500, 440)
	sb = wl_to_s(470, 540)
	sv = wl_to_s(420, 565)
	
	# find "red" L primary
	red = 650
	r = wl_to_s(red, m)
	while (r < sr):
		red -= 5
		r = wl_to_s(red, m)
	
	# find "yellow" LM secondary
	yellow = red
	y = wl_to_s(yellow, m)
	while (y < sy):
		yellow -= 5
		y = wl_to_s(yellow, m)
	
	# find "green" M primary
	green = l # count backwards from L
	g = wl_to_s(green, l)
	while (g > sg):
		green -= 5
		g = wl_to_s(green, l)
		print(green)
	
	# if yellow and green are in impossible places, switch L and M
	if (yellow - green < 5):
		sy = wl_to_s(580, 565)
		sg = wl_to_s(550, 540)
		
		y = wl_to_s(yellow, l)
		while (y < sy):
			yellow += 5
			y = wl_to_s(yellow, l)
		
		g = wl_to_s(green, m)
		while (g < sg):
			green -= 5
			g = wl_to_s(green, m)
	
	# find "cyan" MS secondary
	cyan = s # count forwards from S
	c = wl_to_s(cyan, s)
	while (c > sc):
		cyan += 5
		c = wl_to_s(cyan, s)
	
	# find "blue" S primary
	blue = cyan
	b = wl_to_s(blue, m)
	while (b > sb):
		blue -= 5
		b = wl_to_s(blue, m)
	
	# find "violet" LS secondary
	violet = blue
	v = wl_to_s(violet, l)
	while (v < sv):
		violet -= 5
		v = wl_to_s(violet, l)
	
	return [red, yellow, green, cyan, blue, violet]

# third less accurate version: cone response differences (doesn't work)
def find_primaries_2(l, m, s):
	# find ratios to match
	dr = wl_to_s(620, 565) - wl_to_s(620, 540)
	dy = wl_to_s(580, 565) - wl_to_s(580, 540)
	dg = wl_to_s(550, 565) - wl_to_s(550, 540)
	dc = wl_to_s(500, 540) - wl_to_s(500, 440)
	db = wl_to_s(470, 540) - wl_to_s(470, 440)
	dv = wl_to_s(420, 565) - wl_to_s(420, 540)
	
	# find "red" L primary
	red = 650
	d = wl_to_s(red, l) - wl_to_s(red, m)
	print(d)
	print(dr)
	while (d < dr):
		red -= 5
		d = wl_to_s(red, l) - wl_to_s(red, m)
		print(red)
	
	# find "yellow" LM secondary
	yellow = red
	d = wl_to_s(yellow, l) - wl_to_s(yellow, m)
	while (d > dy):
		yellow -= 5
		d = wl_to_s(yellow, l) - wl_to_s(yellow, m)
	
	# find "green" M primary
	green = yellow
	d = wl_to_s(green, l) - wl_to_s(green, m)
	while (d > dg):
		green -= 5
		d = wl_to_s(green, l) - wl_to_s(green, m)
	
	# find "cyan" MS secondary
	cyan = m
	d = wl_to_s(cyan, m) - wl_to_s(cyan, s)
	while (d > dc):
		cyan -= 5
		d = wl_to_s(cyan, m) - wl_to_s(cyan, s)
	
	# find "blue" S primary
	blue = cyan
	d = wl_to_s(blue, m) - wl_to_s(blue, s)
	while (d > db):
		blue -= 5
		d = wl_to_s(blue, m) - wl_to_s(blue, s)
	
	# find "violet" LS secondary
	violet = blue
	d = wl_to_s(violet, l) - wl_to_s(violet, m)
	while (d < dv):
		violet -= 5
		d = wl_to_s(violet, l) - wl_to_s(violet, m)
	
	return [red, yellow, green, cyan, blue, violet]

def wl_to_rgb_2(wl, l, m, s, p):
	# piecewise hue conversion
	
	# red-yellow: 0-59
	if (p[0] >= wl > p[1]):
		hue = 0 + (p[0] - wl) * 60/(p[0] - p[1])
	# yellow-green: 60-119
	elif (p[1] >= wl > p[2]):
		hue = 60 + (p[1] - wl) * 60/(p[1] - p[2])
	# green-cyan: 120-179
	elif (p[2] >= wl > p[3]):
		hue = 120 + (p[2] - wl) * 60/(p[2] - p[3])
	# cyan-blue: 180-239
	elif (p[3] >= wl > p[4]):
		hue = 180 + (p[3] - wl) * 60/(p[3] - p[4])
	# blue-violet: 240-269
	elif (wl <= p[4]):
		# if S is at a shorter wavelength, "violet" probably isn't seen
		if (p[5] >= s):
			hue = 240
		else:
			hue = 240 + (p[4] - wl) * 30/(p[4] - p[5])
	else:
		hue = 0
	
	return colorsys.hls_to_rgb(hue / 360, 0.5, 1)

# Convert image from BGR to HLS. We use BGR because this is how OpenCV reads images.
# If we use RGB, the output appears normal with settings approximating human vision,
# but shifting the cones produces the opposite of the expected result.
img = cv2.imread(imagename)
img_hls = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)

# switch method used based on input -- this function calls the desired version of wl_to_rgb
# default is 2 (cone response ratios)
if (version == 2):
	primaries = find_primaries(lcone, mcone, scone) # save some time
elif (version == 3):
	primaries = find_primaries_1(lcone, mcone, scone)
elif (version == 4):
	primaries = find_primaries_2(lcone, mcone, scone)

def wl_to_rgb(w, l, m, s):
	if (version == 1):
		return wl_to_rgb_1(w, l, m, s)
	elif (version == 2 or version == 3 or version == 4):
		return wl_to_rgb_2(w, l, m, s, primaries)
	else:
		return wl_to_rgb_0(w, l, m, s)

# transform hues for each pixel
for x in range(img.shape[0]):
	for y in range(img.shape[1]):
		#print("pixel: " + str(x) + ", " + str(y))
		# find pixel
		hls = img_hls[x][y]
		
		# convert hue from 0-180 (OpenCV format) to 0-360
		hue = hls[0]*2
		
		# convert hue to predominant wavelength(s)
		wl = hue_to_wavelength(hue)
		
		if (type(wl) != list):
			hue_target = wl_to_rgb(wl, lcone, mcone, scone)
		else:
			hue_r = wl_to_rgb(620, lcone, mcone, scone)
			hue_b = wl_to_rgb(470, lcone, mcone, scone)
			# sum the amounts of "red" and "blue"
			hue_target = [hue_r[0] * wl[0] + hue_b[0] * wl[1], hue_r[1] * wl[0] + hue_b[1] * wl[1], hue_r[2] * wl[0] + hue_b[2] * wl[1]]
		
		# convert predominant wavelengths back into a hue
		hls_target = colorsys.rgb_to_hls(hue_target[0], hue_target[1], hue_target[2])
		#print(hls_target)
		# shift hue in our pixel. Colorsys uses 0-1, so we have to convert back to
		# OpenCV format.
		img_hls[x][y] = [hls_target[0]*180, hls[1], hls[2]]

# convert back to BGR
img_result = cv2.cvtColor(img_hls, cv2.COLOR_HLS2BGR)

# fix brightness
img_lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
img_result_lab = cv2.cvtColor(img_result, cv2.COLOR_BGR2LAB)
for x in range(img.shape[0]):
	for y in range(img.shape[1]):
		img_result_lab[x][y][0] = img_lab[x][y][0]

img_result = cv2.cvtColor(img_result_lab, cv2.COLOR_LAB2BGR)

# display result
cv2.imwrite("colorvisionpy-result.png", img_result)

# testing
#foo = wl_to_rgb(650, 565, 535, 450)
#print(foo)
#print(colorsys.rgb_to_hls(foo[0], foo[1], foo[2]))
#bar = wl_to_rgb(530, 565, 535, 450)
#print(bar)
#print(colorsys.rgb_to_hls(bar[0], bar[1], bar[2]))
#baz = wl_to_rgb(450, 565, 535, 450)
#print(baz)
#print(colorsys.rgb_to_hls(baz[0], baz[1], baz[2]))
