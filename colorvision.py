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
	# purple-magenta: 270-360 -- just convert it to RGB and use the R and B values
	else:
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)
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
# reverse of the piecewise hue_to_wavelength. The current version uses the sensitivity
# of each cone relative to all three and finds the closest match for these three values
# between 300 and 700 nm.
# When used with the images from studies that used CVS (see also "Experimental evidence that
# primate trichromacy is well suited for detecting primate social colour signals", Hiramatsu
# et al. 2017) this version gives results similar to the actual CVS output shown in these
# papers besides not accounting for saturation. Melin et al. says CVS uses a table of
# "relative sensitivities", so my method should be similar to what it does besides using
# fewer points and more linear interpolation. However, the location of "yellow" is far more
# red-shifted, resulting in a green/blue shift where CVS output shows a red shift.
# Note that I've used the Stockton and Sharpe values of 440-540-565 whereas the studies
# assume 420-530-560 is the unaltered default. Inputting these values produces a red-shifted
# image. (Using the Stockton and Sharpe values produces a 1:1 transformation, as it
# should.)
# In earlier versions, a cone spacing of less than about 8-10 nm (as in uakari vision,
# 420-550-556) usually produced a math overflow error. Now the program will run with this
# type of input but produces a counterintuitive result due to "yellow" being so red-shifted
# that it crosses "red".
# The results are probably not accurate for marsupial-type vision with "yellow-green",
# "cyan" and UV cones (leaving aside the obvious UV issue). Arrese et al. give the primary
# wavelengths for the fat-tailed dunnart as 350, 450 and 620, but inputting its sensitivities
# (363-509-535) gives "red" as 580.5 and "green" as 520.

# match wavelength to sensitivity ratios
def find_match(w, l, m, s):
	lp = wl_to_s(w, 565) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
	mp = wl_to_s(w, 540) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
	sp = wl_to_s(w, 440) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
	
	p0 = 700
	p1 = 300
	
	while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
		+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
		+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > 0.001
		and p0 > 300):
		p0 -= 1
	while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
		+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
		+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > 0.001
		and p1 < 700):
		p1 += 1
	# lower precision
	if (p0 == 300 or p1 == 700):
		print("0.005")
		p0 = 700
		p1 = 300
		while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > 0.005
			and p0 > 300):
			p0 -= 1
		while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > 0.005
			and p1 < 700):
			p1 += 1
	# lower precision
	if (p0 == 300 or p1 == 700):
		print("0.01")
		p0 = 700
		p1 = 300
		while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > 0.01
			and p0 > 300):
			p0 -= 1
		while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > 0.01
			and p1 < 700):
			p1 += 1
	# lower precision
	if (p0 == 300 or p1 == 700):
		print("0.05")
		p0 = 700
		p1 = 300
		while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > 0.05
			and p0 > 300):
			p0 -= 1
		while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > 0.05
			and p1 < 700):
			p1 += 1
	# lower precision
	if (p0 == 300 or p1 == 700):
		print("0.1")
		p0 = 700
		p1 = 300
		while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > 0.1
			and p0 > 300):
			p0 -= 1
		while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > 0.1
			and p1 < 700):
			p1 += 1
	# lower precision
	if (p0 == 300 or p1 == 700):
		print("0.5")
		p0 = 700
		p1 = 300
		while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
			+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > 0.5
			and p0 > 300):
			p0 -= 1
		while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
			+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > 0.5
			and p1 < 700):
			p1 += 1
	
	print(p0)
	print(p1)
	p = (p0 + p1) / 2
	
	return p

def find_primaries(l, m, s):
	# find ratios to match
	
	# red: L/M
	#rr = wl_to_s(620, 565) / wl_to_s(620, 540)
	ry = wl_to_s(580, 565) / wl_to_s(580, 540)
	#rg = wl_to_s(550, 565) / wl_to_s(550, 540)
	
	# yellow: L/(M+S)
	#ry = wl_to_s(580, 565) / (wl_to_s(580, 540) + wl_to_s(580, 440))
	
	# green: (L+M)/S
	rg = (wl_to_s(550, 565) + wl_to_s(550, 540)) / wl_to_s(550, 440)
	
	# cyan and blue: M/S
	rc = wl_to_s(500, 540) / wl_to_s(500, 440)
	rb = wl_to_s(470, 540) / wl_to_s(470, 440)
	
	# violet: L/S
	rv = wl_to_s(420, 565) / wl_to_s(420, 440)
	
	# scale ratios
	#rr1 = wl_to_s(565, 565) / wl_to_s(565, 540) # LM ratio at L peak
	#rg1 = wl_to_s(540, 540) / wl_to_s(540, 565) # ML ratio at M peak
	#rr2 = wl_to_s(l, l) / wl_to_s(l, m) # LM ratio at target L peak
	#rg2 = wl_to_s(m, m) / wl_to_s(m, l) # ML ratio at target M peak
	
	# find "red" L primary
	#red = 650
	#while (wl_to_s(red, l) / wl_to_s(red, m) > rr * rr2 / rr1):
	#while (wl_to_s(red, l) / wl_to_s(red, m) > rr):
	red = find_match(620, l, m, s)
	
	# find "yellow" LM secondary
	#yellow = red
	#while (wl_to_s(yellow, l) / wl_to_s(yellow, m) > ry * rr2 / rr1):
	#while (wl_to_s(yellow, l) / wl_to_s(yellow, m) > ry):
	yellow = find_match(580, l, m, s)
	
	# find "green" M primary
	#green = yellow
	#while (wl_to_s(green, l) / wl_to_s(green, m) > rg * rg2 / rg1):
	#while ((wl_to_s(green, l) + wl_to_s(green, m)) / wl_to_s(green, s) > rg):
	green = find_match(550, l, m, s)
	
	# find "cyan" MS secondary
	#cyan = m
	#while (wl_to_s(cyan, m) / wl_to_s(cyan, s) > rc):
	cyan = find_match(500, l, m, s)
	
	# find "blue" S primary
	#blue = cyan
	#while (wl_to_s(blue, m) / wl_to_s(blue, s) > rb):
	blue = find_match(470, l, m, s)
	
	# find "violet" LS secondary
	#violet = s
	#print(rv)
	#while (wl_to_s(violet, l) / wl_to_s(violet, s) < rv):
	violet = find_match(420, l, m, s)
	
	# print out colors
	print("red: " + str(red))
	print("yellow: " + str(yellow))
	print("green: " + str(green))
	print("cyan: " + str(cyan))
	print("blue: " + str(blue))
	print("violet: " + str(violet))
	
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

def wl_to_rgb(w, l, m, s):
	if (version == 1):
		return wl_to_rgb_1(w, l, m, s)
	elif (version == 2):
		return wl_to_rgb_2(w, l, m, s, primaries)
	else:
		return wl_to_rgb_0(w, l, m, s)

# transform hues for each pixel
for x in range(img.shape[0]):
	for y in range(img.shape[1]):
		# find pixel
		hls = img_hls[x][y]
		
		# convert hue from 0-180 (OpenCV format) to 0-360
		hue = hls[0]*2
		
		# convert hue to predominant wavelength(s)
		wl = hue_to_wavelength(hue)
		
		if (type(wl) != tuple):
			hue_target = wl_to_rgb(wl, lcone, mcone, scone)
		else:
			hue_r = wl_to_rgb(620, lcone, mcone, scone)
			hue_b = wl_to_rgb(470, lcone, mcone, scone)
			# sum the amounts of "red" and "blue"
			hue_target = [hue_r[0] * wl[0] + hue_b[0] * wl[2], hue_r[1] * wl[0] + hue_b[1] * wl[2], hue_r[2] * wl[0] + hue_b[2] * wl[2]]
		
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
