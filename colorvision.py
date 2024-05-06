# This is an attempt to create something similar to Color Vision Simulator (Melin et al. 2013).
# Unlike CVS, it only alters the hue and doesn't account for saturation. I don't know (yet)
# how to determine the distance of a color from the white point of a color space.
# Like CVS, it also doesn't account for differences in lightness/intensity. This could be
# done by finding the global photopic sensitivity curve and dividing by the human version. It
# would also be nice to use LCh instead of HSL so the lightness result would make sense.
#
# References:
# A.D. Melin, D.W. Kline, C.M. Hickey, L.M. Fedigan, Food search through the eyes of a monkey: A functional substitution approach for assessing the ecology of primate color vision, Vision Research, Volume 86, 2013, Pages 87-96, ISSN 0042-6989, https://doi.org/10.1016/j.visres.2013.04.013. (https://www.sciencedirect.com/science/article/pii/S0042698913001119)
# Hogan JD, Fedigan LM, Hiramatsu C, Kawamura S, Melin AD. Trichromatic perception of flower colour improves resource detection among New World monkeys. Sci Rep. 2018 Jul 18;8(1):10883. doi: 10.1038/s41598-018-28997-4. PMID: 30022096; PMCID: PMC6052032.
# Melin AD, Chiou KL, Walco ER, Bergstrom ML, Kawamura S, Fedigan LM. Trichromacy increases fruit intake rates of wild capuchins (Cebus capucinus imitator). Proc Natl Acad Sci U S A. 2017 Sep 26;114(39):10402-10407. doi: 10.1073/pnas.1705957114. Epub 2017 Sep 11. PMID: 28894009; PMCID: PMC5625910.
# Fedigan LM, Melin AD, Addicott JF, Kawamura S. The heterozygote superiority hypothesis for polymorphic color vision is not supported by long-term fitness data from wild neotropical monkeys. PLoS One. 2014 Jan 3;9(1):e84872. doi: 10.1371/journal.pone.0084872. PMID: 24404195; PMCID: PMC3880319.
# Melin AD, Kline DW, Hiramatsu C, Caro T (2016) Zebra Stripes through the Eyes of Their Predators, Zebras, and Humans. PLoS ONE 11(1): e0145679. https://doi.org/10.1371/journal.pone.0145679
# Corso J, Bowler M, Heymann EW, Roos C, Mundy NI. Highly polymorphic colour vision in a New World monkey with red facial skin, the bald uakari (Cacajao calvus). Proc Biol Sci. 2016 Apr 13;283(1828):20160067. doi: 10.1098/rspb.2016.0067. PMID: 27053753; PMCID: PMC4843651.
# Veilleux CC, Scarry CJ, Di Fiore A, Kirk EC, Bolnick DA, Lewis RJ. Group benefit associated with polymorphic trichromacy in a Malagasy primate (Propithecus verreauxi). Sci Rep. 2016 Dec 2;6:38418. doi: 10.1038/srep38418. PMID: 27910919; PMCID: PMC5133583.
# Hiramatsu C, Melin AD, Allen WL, Dubuc C, Higham JP. Experimental evidence that primate trichromacy is well suited for detecting primate social colour signals. Proc Biol Sci. 2017 Jun 14;284(1856):20162458. doi: 10.1098/rspb.2016.2458. PMID: 28615496; PMCID: PMC5474062.
# Stavenga, D.G. On visual pigment templates and the spectral shape of invertebrate rhodopsins and metarhodopsins. J Comp Physiol A 196, 869–878 (2010). https://doi.org/10.1007/s00359-010-0568-7
# Hempel de Ibarra N, Vorobyev M, Menzel R. Mechanisms, functions and ecology of colour vision in the honeybee. J Comp Physiol A Neuroethol Sens Neural Behav Physiol. 2014 Jun;200(6):411-33. doi: 10.1007/s00359-014-0915-1. Epub 2014 May 15. PMID: 24828676; PMCID: PMC4035557.
# CMF From Cone Fundamentals. Horizon Lab @ UCRS. https://horizon-lab.org/colorvis/cone2cmf.html
# Petroc Sumner, Catherine A. Arrese, Julian C. Partridge; The ecology of visual pigment tuning in an Australian marsupial: the honey possum Tarsipes rostratus. J Exp Biol 15 May 2005; 208 (10): 1803–1815. doi: https://doi.org/10.1242/jeb.01610
# Arrese, A. C., Beazley, L. D. and Neumeyer, C. Behavioural evidence for marsupial trichromacy. doi:10.1016/j.cub.2006.02.036
#
# Further reading: https://xkcd.com/1926/

import cv2
import sys
import colorsys
import math
import numpy as np
import time
import argparse
import colormath
from colormath.color_objects import HSLColor, LabColor, sRGBColor
from colormath.color_conversions import convert_color

# execution time
start_time = time.time()

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--lw", type=int, help="longwave cone sensitivity")
parser.add_argument("-m", "--mw", type=int, help="mediumwave cone sensitivity")
parser.add_argument("-s", "--sw", type=int, help="shortwave cone sensitivity")
parser.add_argument("--slw", type=int, default=560, help="source longwave cone sensitivity")
parser.add_argument("--smw", type=int, default=530, help="source mediumwave cone sensitivity")
parser.add_argument("--ssw", type=int, default=420, help="source shortwave cone sensitivity")
parser.add_argument("-i", "--image", help="image name")
parser.add_argument("-v", "--version", type=int, default=2, help="color conversion method")
args = parser.parse_args()

# image
imagename = args.image

# cone peak sensitivities
lcone = args.lw
mcone = args.mw
scone = args.sw

# default reference cone peaks
l0 = args.slw
m0 = args.smw
s0 = args.ssw

# default reference primary/secondary wavelengths
r0 = 650
y0 = 570
g0 = 545
c0 = 490
b0 = 460
v0 = 420

# method to use
version = args.version

# convert hue to wavelength based on location of primary and secondary colors
# Originally I used primary colors based on an sRGB monitor display and my observations
# of white light through a spectroscope (620, 580, 550, 500, 470, 420). I've changed it
# to more closely match the CIE color space so these colors stand in for "perceptual red",
# etc. rather than the limitations of a monitor. This produces a result closer to both CVS
# output and the fat-tailed dunnart color space (Arrese et al.).
def hue_to_wavelength_0(hue):
	#print(hue)
	# red-yellow: 0-59
	if (0 <= hue < 60):
		w = r0 - hue * (r0 - (y0 + 1))/60
	# yellow-green: 60-119
	elif (60 <= hue < 120):
		w = y0 - (hue - 60) * (y0 - (g0 + 1))/60
	# green-cyan: 120-179
	elif (120 <= hue < 180):
		w = g0 - (hue - 120) * (g0 - (c0 + 1))/60
	# cyan-blue: 180-239
	elif (180 <= hue < 240):
		w = c0 - (hue - 180) * (c0 - (b0 + 1))/60
	# blue-violet: 240-269
	elif (240 <= hue < 270):
		w = b0 - (hue - 240) * (b0 - v0)/30
	# purple-magenta: 270-360 -- just convert it to RGB and use the R and B values
	else:
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)
#	print(w)
	return w

# find sensitivity of a given cone type to a wavelength
# This does not use the same function as CVS. The function comes from a paper describing
# templates for visual pigment absorbance that are meant to fit both vertebrates and
# invertebrates (Stavenga 2010), whereas CVS uses shifted versions of the 10° human cone
# fundamentals. This is probably close enough and better for non-primates.
def wl_to_s(wl, peak):
	value = (math.exp(69.7*(0.8795+0.0459*math.exp(-(peak-300)**2/11940)-(peak/wl)))+math.exp(28*(0.922-(peak/wl)))+math.exp(-14.9*(1.104-(peak/wl)))+0.674)**-1 + 0.26*math.exp(-((wl-(189+0.315*peak))/(-40.5+0.195*peak))**2)
	return value

# convert wavelength to RGB values for target color space
# How CVS does this is described as follows:
# "Using the idealized cone fundamentals, and based on standard daylight (D65) illumination, the program generates a color space table for the source and target phenotypes that lists the relative sensitivities of the component photopigments to each wavelength in 1 nm increments – intermediate values are derived from the table via linear interpolation. For each image in the trial, the starting RGB values for each pixel are converted to chromaticity and intensity (luminance) values. The chromaticity values are located in the source color space table to find the predominant hue (wavelength); the saturation value is determined by the relative distance of the chromaticity value from the white point. The hue value is then located in the target color space table, and the target saturation value is projected from the white point, to find the modified chromaticity values. An inverse transformation to the new RGB color space generates the task stimulus. The luminance values are held constant during this process. In the case of monochromatic simulations, only the luminance values of each pixel are used to generate gray-scale images."

# version 0: use unaltered cone fundamentals
# This just takes the sensitivity values for all three cones and uses them as the R, G and B
# values. For wide spacing (>40 nm) as in birds, insects and some marsupials, where the
# primary wavelengths are received almost exclusively by one type of cone, this is probably
# pretty close and may be better than the other versions. For narrower spacing as in primates,
# it's way off because the interactions become more complex. There are no values that produce
# a 1:1 transformation because the hue-wavelength conversion is done completely differently
# using the linear hue_to_wavelength_0().
def wl_to_rgb_0(w, l, m, s):
	red = wl_to_s(w, l)
	green = wl_to_s(w, m)
	blue = wl_to_s(w, s)
	return [red, green, blue]

# version 1: use color matching functions
# This is an attempt to create a proper color space, based on the equations on the
# Horizon Lab page. Horizon Lab says primaries can be arbitrary, and I don't know a way to
# find them "directly", so we use the method in find_match/find_primaries (matching ratios).
# Since we can directly convert between hues and wavelengths once we have the primaries, we
# don't use a linear approximation (this doesn't work properly).
# With 420-530-560 as the source LMS values, this produces similar results to CVS but with
# increased red-green contrast. This may be more realistic (see Webster et al. 2011) but
# indicates CVS uses a different method to construct color spaces. Possibly it uses L-M
# and S-LM contrast similarly to CIELAB.
# For highly blue-shifted target M and S values as in marsupials or the 420-490-560 example,
# the results are blatantly wrong, with much of what should be the yellow/orange region
# turning red due to negative or zero values. I don't know if this is because of poorly
# chosen primaries or another reason. Note the monkey face image appears correct because
# it only contains yellow/red hues; any image containing green or blue makes the problem
# obvious. The results are also wrong for very narrow spacing -- the uakari image with
# 420-550-556 turns cyan.
def cmf(l, m, s, r, g, b):

	# sensitivity of cones to primaries
	matrix_a = np.array([
		[wl_to_s(r, l), wl_to_s(g, l), wl_to_s(b, l)],
		[wl_to_s(r, m), wl_to_s(g, m), wl_to_s(b, m)],
		[wl_to_s(r, s), wl_to_s(g, s), wl_to_s(b, s)]
	])

	# sensitivity of cones to each wavelength from 300 to 800 nm in 1 nm increments
	matrix_c = np.empty((3, 501)) # 500x3 matrix -- this is height x width, not width x height
	
	# red row
	for i in range(0, 500):
		matrix_c[0][i] = wl_to_s(i + 300, l)
	# green row
	for i in range(0, 500):
		matrix_c[1][i] = wl_to_s(i + 300, m)
	# blue row
	for i in range(0, 500):
		matrix_c[2][i] = wl_to_s(i + 300, s)

	# initial color match matrix
	matrix_m = np.matmul(np.linalg.inv(matrix_a), matrix_c) # A x M = C, so M = A^-1 x C
	
	# sum every row of the matrix to find the response to equal energy white (I don't know
	# how to find D50)
	rw = 0.0
	gw = 0.0
	bw = 0.0
	for i in range (0, 500):
		rw += matrix_m[0][i]
		gw += matrix_m[1][i]
		bw += matrix_m[2][i]

	# final color match matrix
	matrix_cmf = np.empty((3, 501))
	for i in range(0, 500):
		matrix_cmf[0][i] = matrix_m[0][i] / rw
		matrix_cmf[1][i] = matrix_m[1][i] / gw
		matrix_cmf[2][i] = matrix_m[2][i] / bw
	
	return matrix_cmf

# tables
def cs_tables_1(l, m, s):
	cs_source = np.empty((501, 3))
	cs_target = np.empty((501, 3))
	
	# primaries
	r = find_match(r0, l, m, s)
	g = find_match(g0, l, m, s)
	b = find_match(b0, l, m, s)
	
	r = (r[0] + r[1]) / 2
	g = (g[0] + g[1]) / 2
	b = (b[0] + b[1]) / 2
	
	# find raw CMF tables
	cmf_source = cmf(l0, m0, s0, r0, g0, b0)
	cmf_target = cmf(l, m, s, r, g, b)
	
	for i in range(0, 500):
		rgb_source = [cmf_source[0][i], cmf_source[1][i], cmf_source[2][i]]
		rgb_target = [cmf_target[0][i], cmf_target[1][i], cmf_target[2][i]]
		
		# adjust negative values
		if (rgb_source[0] < 0):
			rgb_source[1] += abs(rgb_source[0])
			rgb_source[2] += abs(rgb_source[0])
			rgb_source[0] = 0
		if (rgb_source[1] < 0):
			rgb_source[0] += abs(rgb_source[1])
			rgb_source[2] += abs(rgb_source[1])
			rgb_source[1] = 0
		if (rgb_source[2] < 0):
			rgb_source[0] += abs(rgb_source[2])
			rgb_source[1] += abs(rgb_source[2])
			rgb_source[2] = 0
		
		if (rgb_target[0] < 0):
			rgb_target[1] += abs(rgb_target[0])
			rgb_target[2] += abs(rgb_target[0])
			rgb_target[0] = 0
		if (rgb_target[1] < 0):
			rgb_target[0] += abs(rgb_target[1])
			rgb_target[2] += abs(rgb_target[1])
			rgb_target[1] = 0
		if (rgb_target[2] < 0):
			rgb_target[0] += abs(rgb_target[2])
			rgb_target[1] += abs(rgb_target[2])
			rgb_target[2] = 0
		
		cs_source[i][0] = rgb_source[0] / (rgb_source[0] + rgb_source[1] + rgb_source[2])
		cs_source[i][1] = rgb_source[1] / (rgb_source[0] + rgb_source[1] + rgb_source[2])
		cs_source[i][2] = rgb_source[2] / (rgb_source[0] + rgb_source[1] + rgb_source[2])
		
		cs_target[i][0] = rgb_target[0] / (rgb_target[0] + rgb_target[1] + rgb_target[2])
		cs_target[i][1] = rgb_target[1] / (rgb_target[0] + rgb_target[1] + rgb_target[2])
		cs_target[i][2] = rgb_target[2] / (rgb_target[0] + rgb_target[1] + rgb_target[2])
		
	return [cs_source, cs_target]

def hue_to_wavelength_1(h, source):
	h = round(h)
	
	# spectral colors
	if (0 <= h <= 270):
		i = 0
		
		while (i < 500):
			h_source = round(colorsys.rgb_to_hls(source[i][0], source[i][1], source[i][2])[0] * 360)
			h_prev = round(colorsys.rgb_to_hls(source[i - 1][0], source[i - 1][1],
			source[i - 1][2])[0] * 360)
			# first exact match. I tried averaging all the exact matches, but this seems to just
			# increase the execution time.
			if (h_source == h):
				return i + 300
			# if no exact match, use an intermediate value
			elif (h_prev > h and h_source < h):
				return i - (h - h_source) / (h_prev - h_source) + 300
			else:
				i += 1
		# otherwise
		return i + 300
	
	# non-spectral colors -- we don't handle them, we just do what the other hue_to_wavelength does
	else:
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)

def wl_to_rgb_1(w, target):
	if (w == round(w)):
		return target[round(w) - 300]
	else: # if non-integer value, use a weighted average
		d = w - int(w)
		return target[math.floor(w) - 300] * (1 - d) + target[math.ceil(w) - 300] * d

# version 2: cone response ratios
# we find the wavelengths that produce ratios matching those for the primaries and
# secondaries in human vision, then use those to convert the wavelength to a hue as the
# reverse of the piecewise hue_to_wavelength. The current version uses the sensitivity
# of each cone relative to all three and finds the closest match for these three values
# between 300 and 800 nm.
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
# type of input but produces a counterintuitive result due to the green/blue shift.
# The problem is that while red/green colors approach "yellow" as the L and M curves
# approach each other, perceptual "yellow" to a normal human has an L:M ratio much greater
# than 1:1, meaning a transformation based on this will cause a shift toward green rather
# than yellow.
# The results are also not accurate for marsupial-type vision with "yellow-green",
# "cyan" and UV cones (leaving aside the obvious UV issue). Arrese et al. give the primary
# wavelengths for the fat-tailed dunnart as 350, 450 and 620, but inputting its sensitivities
# (363-509-535) gives "red" and "green" values that are significantly closer together.
# Due to the above issues, the results are probably better for more widely spaced L/M
# cones as in the honey possum (350-505-557).

# match wavelength to sensitivity ratios
def find_match(w, l, m, s, d=0.001):
	lp = wl_to_s(w, l0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
	mp = wl_to_s(w, m0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
	sp = wl_to_s(w, s0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
	
	p0 = 700
	p1 = 300
	
	#print(d)
	while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
		+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
		+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > d
		and p0 > 300):
		p0 -= 1
	while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
		+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
		+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > d
		and p1 < 700):
		p1 += 1
	# lower precision if not found
	if (p0 == 300 and p1 == 700):
		# recursion
		p2 = find_match(w, l, m, s, d*2) # avoid running this twice
		p0 = p2[0]
		p1 = p2[1]
	
	#print(p0)
	#print(p1)
	
	return([p0, p1])

def find_primaries(l, m, s):
	# find "red" L primary
	red = (find_match(r0, l, m, s)[0] + find_match(r0, l, m, s)[1]) / 2
	
	# find "yellow" LM secondary
	yellow = (find_match(y0, l, m, s)[0] + find_match(y0, l, m, s)[1]) / 2
	
	# find "green" M primary
	green = (find_match(g0, l, m, s)[0] + find_match(g0, l, m, s)[1]) / 2
	
	# find "cyan" MS secondary
	cyan = (find_match(c0, l, m, s)[0] + find_match(c0, l, m, s)[1]) / 2
	
	# find "blue" S primary
	blue = (find_match(b0, l, m, s)[0] + find_match(b0, l, m, s)[1]) / 2
	
	# find "violet" LS secondary
	violet = (find_match(v0, l, m, s)[0] + find_match(v0, l, m, s)[1]) / 2
	
	# patch results if they're out of order
	if (red < yellow):
		red = 800
	if (cyan < s):
		cyan = (green + s) / 2
	if (blue < violet):
		violet = 300
	
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

# version 3
# These functions try to match every wavelength to the closest set of ratios as CVS does
# instead of finding specific "primary colors" and filling in the gaps linearly. It doesn't
# fully solve the yellow/green issue, but for very narrow spacing the results are significantly
# better because it doesn't force a "green point". For wider spacing, on the other hand, the
# shifting of red/yellow (in either direction) seems to be underestimated. Running this version
# takes >100 times longer than version 2 because it has to search a table with 500 entries for
# every pixel.
def cs_tables(l, m, s):
	# source color space table
	cs_source = np.empty([501, 4])
	for i in range(cs_source.shape[0]):
		w = i + 300
		rgb = wl_to_rgb_2(w, l0, m0, s0, [r0, y0, g0, c0, b0, v0])
		cs_source[i][0] = wl_to_s(w, l0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
		cs_source[i][1] = wl_to_s(w, m0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
		cs_source[i][2] = wl_to_s(w, s0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
		cs_source[i][3] = colorsys.rgb_to_hls(rgb[0], rgb[1], rgb[2])[0] * 360
	
	# target color space table (doesn't work, I don't know why)
	#cs_target = np.empty([501, 3])
	#for i in range(cs_target.shape[0]):
	#	w = i + 300
	#	cs_target[i][0] = wl_to_s(w, l) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
	#	cs_target[i][1] = wl_to_s(w, m) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
	#	cs_target[i][2] = wl_to_s(w, s) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
	
	#return [cs_source, cs_target]
	return cs_source

# Note this now returns a position in the source table rather than a hue.
def find_match_1(w, source, l, m, s, d=0.001):
	# target sensitivities
	#lw = target[int(w) - 300][0]
	#mw = target[int(w) - 300][1]
	#sw = target[int(w) - 300][2]
	lw = wl_to_s(w, l) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
	mw = wl_to_s(w, m) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
	sw = wl_to_s(w, s) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
	#print(l)
	#print(m)
	#print(s)
	#print(w)
	#print(lw)
	#print(mw)
	#print(sw)
	
	# test
	#print("source table values:")
	#print(source[0]) # red
	#print(source[60]) # yellow
	#print(source[120]) # green
	#print(source[180]) # cyan
	#print(source[240]) # blue
	#print(source[269]) # violet
	
	# find nearest match
	i = 0
	j = 500
	
	# forward
	while (abs(source[i][0] - lw)
		+ abs(source[i][1] - mw)
		+ abs(source[i][2] - sw) > d
		and i < 500):
		i += 1
		#print(source[i])
		#print(abs(source[i][0] - l)
		#+ abs(source[i][1] - m)
		#+ abs(source[i][2] - s))
		#print(i)
	# back
	while (abs(source[j][0] - lw)
		+ abs(source[j][1] - mw)
		+ abs(source[j][2] - sw) > d
		and j > 0):
		j -= 1
		#print(j)
	#print("i")
	#print(i)
	#print("j")
	#print(j)
	
	# lower precision if not found
	if (i == 500 and j == 0):
		#print("recursion")
		# recursion
		match = find_match_1(w, source, l, m, s, d*2) # avoid running this twice
		i = match[0]
		j = match[1]
	
	return([i, j])

# version 4
# Another attempt at constructing a color space. This is my best guess at what CVS does:
# a color space based on the L-M and (L+M)-S differences, using these as X and Y coordinates.
# The results are surprisingly close for red/green hues, but blue/violet isn't handled properly
# and converts into unrelated colors.
def rg(w, l, m, s):
	return (wl_to_s(w, l) - wl_to_s(w, m)) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))

def yb(w, l, m, s):
	return (wl_to_s(w, l) + wl_to_s(w, m) - wl_to_s(w, s)) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))

def cs_tables_2(l, m, s):
	cs_source = np.empty([501, 3])
	cs_target = np.empty([501, 3])
	
	# scaling factors
	red_rg = rg(r0, l0, m0, s0)
	green_rg = rg(g0, l0, m0, s0)
	yellow_yb = yb(y0, l0, m0, s0)
	blue_yb = yb(b0, l0, m0, s0)
	
	# source white point
	white_rg = 0
	white_yb = 0
	for i in range(0, 500):
		white_rg += rg(i + 300, l0, m0, s0)
		white_yb += yb(i + 300, l0, m0, s0)
	
	white_rg = white_rg / 500
	white_yb = white_yb / 500
	
	# target white point
	white_rg1 = 0
	white_yb1 = 0
	for i in range(0, 500):
		white_rg1 += rg(i + 300, l, m, s)
		white_yb1 += yb(i + 300, l, m, s)
	
	white_rg1 = white_rg1 / 500
	white_yb1 = white_yb1 / 500
	
	# coordinates of "primary colors"
	red = ((rg(r0, l0, m0, s0) - white_rg) / abs(red_rg - white_rg), (yb(r0, l0, m0, s0) - white_yb) / abs(yellow_yb - white_yb))
	green = ((rg(g0, l0, m0, s0) - white_rg) / abs(red_rg - white_rg), (yb(g0, l0, m0, s0) - white_yb) / abs(yellow_yb - white_yb))
	yellow = ((rg(y0, l0, m0, s0) - white_rg) / abs(red_rg - white_rg), (yb(y0, l0, m0, s0) - white_yb) / abs(yellow_yb - white_yb))
	blue = ((rg(b0, l0, m0, s0) - white_rg) / abs(red_rg - white_rg), (yb(b0, l0, m0, s0) - white_yb) / abs(yellow_yb - white_yb))
	cyan = ((rg(c0, l0, m0, s0) - white_rg) / abs(red_rg - white_rg), (yb(c0, l0, m0, s0) - white_yb) / abs(yellow_yb - white_yb))
	
	# find and scale RG and YB coordinates for wavelengths 300-800
	for i in range(0, 500):
		# populate source table
		rg_source = rg(i + 300, l0, m0, s0) - white_rg
		#if (rg_source > 0):
		rg_source = rg_source / abs(red_rg - white_rg)
		#elif (rg_source < 0):
		#	rg_source = rg_source / abs(green_rg - white_rg)
		cs_source[i][0] = rg_source
		
		yb_source = yb(i + 300, l0, m0, s0) - white_yb
		#if (yb_source > 0):
		yb_source = yb_source / abs(yellow_yb - white_yb)
		#elif (rg_source < 0):
		#	yb_source = yb_source / abs(blue_yb - white_yb)
		cs_source[i][1] = yb_source
		
		# populate target table
		rg_target = rg(i + 300, l, m, s) - white_rg1
		#if (rg_target > 0):
		rg_target = rg_target / abs(red_rg - white_rg)
		#elif (rg_target < 0):
		#	rg_target = rg_target / abs(green_rg - white_rg)
		cs_target[i][0] = rg_target
		
		yb_target = yb(i + 300, l, m, s) - white_yb1
		#if (yb_target > 0):
		yb_target = yb_target / abs(yellow_yb - white_yb)
		#elif (yb_target < 0):
		#	yb_target = yb_target / abs(blue_yb - white_yb)
		cs_target[i][1] = yb_target
		
		# convert from RGYB to RGB
		rgb_source = rgyb_to_rgb(cs_source[i][0], cs_source[i][1], red, green, blue, (white_rg, white_yb))
		rgb_target = rgyb_to_rgb(cs_target[i][0], cs_target[i][1], red, green, blue, (white_rg, white_yb))
		#print(rgb_source)
		#print(rgb_target)
		cs_source[i][2] = colorsys.rgb_to_hls(rgb_source[0], rgb_source[1], rgb_source[2])[0] * 360
		cs_target[i][2] = colorsys.rgb_to_hls(rgb_target[0], rgb_target[1], rgb_target[2])[0] * 360
		#print(i + 300)
		#print(cs_source[i])
		#print(cs_target[i])
	
	return [cs_source, cs_target]

def rgyb_to_rgb(rg, yb, red, green, blue, white):
	#print(values)
	# clamp
	if (rg < blue[0]):
		rg = blue[0]
	if (yb < blue[1]):
		yb = blue[1]
	if (rg > red[0]):
		rg = red[0]
	if (yb > red[1]):
		yb = red[1]
	
	r = 1
	g = 1
	b = 1
	
	# This doesn't work that way. "Blue" is way off to the left of "green".
	# first quadrant
	#if (rg > 0 and yb > 0):
	if (yb > 0):
		r = abs(1 + rg)
		g = abs(1 - rg)
		b = abs(1 - yb)
	#print([r, g, b])
	# second quadrant
	#elif (rg > 0 and yb < 0):
	#	r = math.dist(blue, (rg, yb)) # distance from blue
	#	g = 1 - math.dist((0, 0), (rg, yb)) # distance from white
	#	b = 1 - math.dist(blue, (rg, yb)) / 2*math.sqrt(2) # negative distance from blue
	# third quadrant
	#elif (rg < 0 and yb < 0):
	#	r = 1 - math.dist((0, 0), (rg, yb)) # distance from white
	#	g = math.dist(blue, (rg, yb)) # distance from blue
	#	b = 1 - math.dist(blue, (rg, yb)) / 2*math.sqrt(2) # negative distance from blue
	# fourth quadrant
	#elif (rg < 0 and yb > 0):
	#	r = 1 + rg
	#	g = 1
	#	b = 1 - yb
	elif (yb < 0):
		r = abs(1 + rg + yb)
		g = abs(1 - rg + yb)
		b = abs(1 - yb)
	if (r > 1):
		#g = g / r
		#b = b / r
		r = 1
	if (g > 1):
		#r = r / g
		#b = b / g
		g = 1
	if (b > 1):
		r = r / b
		g = g / b
		b = 1

	# distance to lines defining the color space (this probably doesn't work that way either
	# because non-linear, but it should be closer)
	#dist_r = abs((blue[0] - green[0]) * (red[1] - green[1]) - (red[0] - green[0]) * (blue[1] - green[1])) / math.sqrt((blue[0] - green[0])**2 + (blue[1] - green[1])**2)
	#dist_g = abs((blue[0] - red[0]) * (green[1] - red[1]) - (green[0] - red[0]) * (blue[1] - red[1])) / math.sqrt((blue[0] - red[0])**2 + (blue[1] - red[1])**2)
	#dist_b = abs((red[0] - green[0]) * (blue[1] - green[1]) - (blue[0] - green[0]) * (red[1] - green[1])) / math.sqrt((red[0] - green[0])**2 + (red[1] - green[1])**2)
	
	#r = (abs((blue[0] - green[0]) * (yb - green[1]) - (rg - green[0]) * (blue[1] - green[1])) / math.sqrt((blue[0] - green[0])**2 + (blue[1] - green[1])**2)) / dist_r
	#g = (abs((blue[0] - red[0]) * (yb - red[1]) - (rg - red[0]) * (blue[1] - red[1])) / math.sqrt((blue[0] - red[0])**2 + (blue[1] - red[1])**2)) / dist_g
	#b = (abs((red[0] - green[0]) * (yb - green[1]) - (rg - green[0]) * (red[1] - green[1])) / math.sqrt((red[0] - green[0])**2 + (red[1] - green[1])**2)) / dist_b
	
	return [r, g, b]

def hue_to_wavelength_2(h, source):
	#hsl = HSLColor(h, 0.5, 0.5)
	#lab = convert_color(hsl, LabColor).get_value_tuple()
	#a = round(lab[1])
	#b = round(lab[2])
	
	for i in range(0, 500):
		h_cur = source[i][2]
		h_next = source[i + 1][2]
		#print(rgb_cur)
		
		
		# exact match
		if (round(h) == round(h_cur)):
			#print(i)
			return i + 300
		# intermediate match
		elif (h_next <= h <= h_cur or h_cur <= h <= h_next):
			#print(i)
			return i + 0.5 + 300
	
	# not found
	#print(i)
	return 800

def wl_to_rgb_4(w, target):
	#print(w)
	if (w == round(w)):
		return colorsys.hls_to_rgb(target[round(w) - 300][2] / 360, 0.5, 1)
	else:
		return colorsys.hls_to_rgb(((target[round(w) - 301][2] + target[round(w) - 300][2]) / 2) / 360, 0.5, 1)
	
	#return convert_color(LabColor(50, ab[0], ab[1]), sRGBColor).get_value_tuple()

# Convert image from BGR to HLS. We use BGR because this is how OpenCV reads images.
# If we use RGB, the output appears normal with settings approximating human vision,
# but shifting the cones produces the opposite of the expected result.
img = cv2.imread(imagename)
img_hls = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)

# switch method used based on input -- this function calls the desired version of wl_to_rgb
# default is 2 (cone response ratios)

# save some time by taking these out of wl_to_rgb so we don't call them for every pixel
if (version == 1):
	tables = cs_tables_1(lcone, mcone, scone)
elif (version == 2):
	primaries = find_primaries(lcone, mcone, scone)
elif (version == 3):
	tables = cs_tables(lcone, mcone, scone)
elif (version == 4):
	tables = cs_tables_2(lcone, mcone, scone)

# wrapper functions
def wl_to_rgb(w, l, m, s):
	if (version == 1):
		#return wl_to_rgb_1(w, l, m, s, primaries)
		return wl_to_rgb_1(w, tables[1])
	elif (version == 2):
		return wl_to_rgb_2(w, l, m, s, primaries)
	elif (version == 3):
		#match = find_match_1(w, tables[0], l, m, s)
		match = find_match_1(w, tables, l, m, s)
		#hue = (match[0] + match[1]) / 2
		hue = (tables[match[0]][3] + tables[match[1]][3]) / 2
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)
	elif (version == 4):
		return wl_to_rgb_4(w, tables[1])
	else:
		return wl_to_rgb_0(w, l, m, s)

def hue_to_wavelength(h):
	if (version == 1):
		return hue_to_wavelength_1(h, tables[0])
	elif (version == 4):
		return hue_to_wavelength_2(h, tables[0])
	else:
		return hue_to_wavelength_0(h)

# transform hues for each pixel
for x in range(img.shape[0]):
	for y in range(img.shape[1]):
		# report current pixel to determine whether there's an infinite loop
		#print("pixel: " + str(x) + ", " + str(y))
		# find pixel
		hls = img_hls[x][y]
		
		# convert hue from 0-180 (OpenCV format) to 0-360
		hue = hls[0]*2
		
		# convert hue to predominant wavelength(s)
		wl = hue_to_wavelength(hue)
		#print(hue)
		
		if (type(wl) != tuple):
			hue_target = wl_to_rgb(wl, lcone, mcone, scone)
		else:
			hue_r = wl_to_rgb(r0, lcone, mcone, scone)
			hue_b = wl_to_rgb(b0, lcone, mcone, scone)
			# sum the amounts of "red" and "blue"
			hue_target = [hue_r[0] * wl[0] + hue_b[0] * wl[2], hue_r[1] * wl[0] + hue_b[1] * wl[2], hue_r[2] * wl[0] + hue_b[2] * wl[2]]
		
		# convert predominant wavelengths back into a hue
		#print(hue_target)
		hls_target = colorsys.rgb_to_hls(hue_target[0], hue_target[1], hue_target[2])
		#print(hls_target)
		# shift hue in our pixel. Colorsys uses 0-1, so we have to convert back to
		# OpenCV format.
		img_hls[x][y] = [hls_target[0]*180, hls[1], hls[2]]

# convert back to BGR
img_result = cv2.cvtColor(img_hls, cv2.COLOR_HLS2BGR)

# fix brightness
#img_lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
#img_result_lab = cv2.cvtColor(img_result, cv2.COLOR_BGR2LAB)
#for x in range(img.shape[0]):
#	for y in range(img.shape[1]):
#		img_result_lab[x][y][0] = img_lab[x][y][0]

#img_result = cv2.cvtColor(img_result_lab, cv2.COLOR_LAB2BGR)

# display result
cv2.imwrite("colorvisionpy-result.png", img_result)

# print execution time
print("%s seconds" % (time.time() - start_time))
