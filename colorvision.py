# This is an attempt to create something similar to Color Vision Simulator (Melin et al. 2013).
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
from colormath.color_objects import LabColor, LCHabColor, sRGBColor, XYZColor, HSLColor
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
parser.add_argument("--cutoff", type=int, default=300, help="lower limit of lens transmission")
parser.add_argument("-r", "--red", type=int, default=650, help="reference red")
parser.add_argument("-y", "--yellow", type=int, default=570, help="reference yellow")
parser.add_argument("-g", "--green", type=int, default=545, help="reference green")
parser.add_argument("-c", "--cyan", type=int, default=490, help="reference cyan")
parser.add_argument("-b", "--blue", type=int, default=460, help="reference blue")
parser.add_argument("--violet", type=int, default=420, help="reference violet")
parser.add_argument("-i", "--image", help="image name")
parser.add_argument("-v", "--version", type=int, default=2, help="color conversion method")
args = parser.parse_args()

# cone peak sensitivities
l1 = args.lw
m1 = args.mw
s1 = args.sw

# default reference cone peaks
l0 = args.slw
m0 = args.smw
s0 = args.ssw

# default reference primary/secondary wavelengths
r0 = args.red
y0 = args.yellow
g0 = args.green
c0 = args.cyan
b0 = args.blue
v0 = args.violet

# method to use
version = args.version

# relative wavelength sensitivity
# This is not derived from any specific source. Instead we add the L, M and S functions weighted
# by the proportion of L, M and S cones, assumed to be 10:5:1, and multiply this by a function
# resembling a graph of typical lens transmission adjusted to reach 1 at 800 nm. The result
# for human-like values is similar (but not identical) to the CIE luminous efficiency function.
def sensitivity(w, l, m, s, c=300):
	value = (wl_to_s(w, l) + wl_to_s(w, m) / 2 + wl_to_s(w, s) / 10) * np.log(w - c) / np.log(800 - c)
	if (value > 0): # a log function has the wrong shape near the cutoff point
		return value
	return 0

# black body
# This doesn't work because Python doesn't like small numbers and the denominator is rounded
# to 1 - 1 = 0. This is probably also why Gnuplot refuses to plot a black body curve.
#def blackbody(w, t):
#	h = 6.62607015e-34
#	c = 299792458
#	k = 1.380649e-23
	#print(math.exp(h*c / w*k*t))
#	return w**-4 / (math.exp(0.1) - 1)

# convert hue to wavelength based on location of primary and secondary colors
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
# This is a conversion to "LMS color space": it just takes the sensitivity values for all three
# cones and uses them as the R, G and B values. For wide spacing (>40 nm) as in birds, insects
# and some marsupials, where the primary wavelengths are received almost exclusively by one type
# of cone, this is probably pretty close and may be better than the other versions. For narrower
# spacing as in primates, it's way off because the interactions become more complex. There are no
# values that produce a 1:1 transformation because the hue-wavelength conversion is done completely
# differently using the linear hue_to_wavelength_0().
# This produces a passable simulation of dichromacy or monochromacy if two or all of the LMS
# values are identical. Dichromacy is slightly wrong because of the "beta-band" causing blue
# hues to appear desaturated. Monochromacy should be identical to CVS because it only uses
# the original lightness values, but the method for determining lightness/intensity is
# evidently slightly different and produces darker shadows.
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
# papers besides not accounting for saturation. However, the location of "yellow" is far more
# red-shifted, resulting in a green/blue shift where CVS output shows a red shift.
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
	#print("red: " + str(red))
	#print("yellow: " + str(yellow))
	#print("green: " + str(green))
	#print("cyan: " + str(cyan))
	#print("blue: " + str(blue))
	#print("violet: " + str(violet))
	
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
# The origin is set as the white point based on equal energy white, and the values are scaled
# based on arbitrary "red" and "yellow" primary wavelengths. The transformation to RGB from
# this "RGYB" space is assumed to be linear.
# Results with similar input are noticeably different but surprisingly close and better than
# any of the other versions, particularly with very narrow spacing; this is the only version
# that gets uakari trichromacy right. It can also produce a passable simulation of dichromacy
# if the target L and M values are identical. However, with other arrangements such as the
# 420-490-560 example, the blue/green region disintegrates and red becomes pink. Some
# blue/violet hues can't be processed properly regardless of input and are assigned red or
# cyan.
# "Dichromatic" output is much less blue than what CVS produces. Based on other simulations
# of dichromacy, I think CVS heavily overestimates the amount of blue, possibly because
# it cuts the "yellow" number in half by reducing the number of cone types and then tries
# to use the resulting (L+M) - S value the same way. Giving the L/M cone double weight by
# treating it as "identical L and M cones" is probably closer to the truth.
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
	blue = ((rg(b0, l0, m0, s0) - white_rg) / abs(red_rg - white_rg), (yb(b0, l0, m0, s0) - white_yb) / abs(yellow_yb - white_yb))
	
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
		rgb_source = rgyb_to_rgb(cs_source[i][0], cs_source[i][1], red, blue)
		rgb_target = rgyb_to_rgb(cs_target[i][0], cs_target[i][1], red, blue)
		#print(rgb_source)
		#print(rgb_target)
		cs_source[i][2] = colorsys.rgb_to_hls(rgb_source[0], rgb_source[1], rgb_source[2])[0] * 360
		cs_target[i][2] = colorsys.rgb_to_hls(rgb_target[0], rgb_target[1], rgb_target[2])[0] * 360
		#print(i + 300)
		#print(cs_source[i])
		#print(cs_target[i])
	
	return [cs_source, cs_target]

def rgyb_to_rgb(rg, yb, red, blue):
	#print(values)
	# clamp
	#if (rg < blue[0]):
	#	rg = blue[0]
	#if (yb < blue[1]):
	#	yb = blue[1]
	#if (rg > red[0]):
	#	rg = red[0]
	#if (yb > red[1]):
	#	yb = red[1]
	
	r = 1
	g = 1
	b = 1
	
	if (yb > 0):
		r = abs(1 + rg)
		g = abs(1 - rg)
		b = abs(1 - yb)
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
	
	return [r, g, b]

def hue_to_wavelength_2(h, source):
	if (0 <= h <= 270):
		for i in range(0, 500):
			h_cur = source[i][2]
			h_next = source[i + 1][2]
			
			# exact match
			if (round(h) == round(h_cur)):
				return i + 300
			# intermediate match
			elif (h_next <= h <= h_cur or h_cur <= h <= h_next):
				return i + 0.5 + 300
		
		return 800
	else:
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)

def wl_to_rgb_4(w, target):
	if (w == round(w)):
		return colorsys.hls_to_rgb(target[round(w) - 300][2] / 360, 0.5, 1)
	else:
		return colorsys.hls_to_rgb(((target[round(w) - 301][2] + target[round(w) - 300][2]) / 2) / 360, 0.5, 1)

# It's not clear how exactly CVS alters the saturation. What I've done here is compare the
# distances from the white point for a wavelength's entries in the source and target tables.
# The results of this are a bit patchy but broadly resemble CVS output.
def saturation(w, source, target):
	if (type(w) != tuple):
		if (w == round(w)):
			sx = source[round(w) - 300][0]
			sy = source[round(w) - 300][1]
			tx = target[round(w) - 300][0]
			ty = target[round(w) - 300][1]
		else:
			sx = (source[round(w) - 301][0] + source[round(w) - 300][0]) / 2
			sy = (source[round(w) - 301][1] + target[round(w) - 300][1]) / 2
			tx = (source[round(w) - 301][0] + source[round(w) - 300][0]) / 2
			ty = (source[round(w) - 301][1] + target[round(w) - 300][1]) / 2
	else:
		sx = source[r0 - 300][0] * w[0] + source[b0 - 300][0] * w[1]
		sy = source[r0 - 300][1] * w[0] + source[b0 - 300][1] * w[1]
		tx = target[r0 - 300][0] * w[0] + target[b0 - 300][0] * w[1]
		ty = target[r0 - 300][1] * w[0] + target[b0 - 300][1] * w[1]
	
	sat_s = math.dist((0, 0), (sx, sy))
	sat_t = math.dist((0, 0), (tx, ty))
	if (sat_s == 0):
		return 0
	return sat_t / sat_s

# version 5
# This converts LMS values to CIE XYZ and then to HLS. See https://en.wikipedia.org/wiki/LMS_color_space for the conversion matrix.
# The results are very similar to CVS and to version 4 aside from a slight green shift for
# narrow spacing, probably for the same reason as in version 2. The main problem is the
# range of hues provided for integer wavelengths is full of holes you could drive a truck
# through. These gaps can be filled in convincingly with some careful calculations (probably
# the reason for the long execution time) for primate-like values, but for other arrangements
# it doesn't work as well.
# XYZ conversion can also simulate dichromacy in the same way as v0 and v4. Tritanopia is a
# bit odd-looking as red becomes purple for some reason.
def cs_tables_3(l, m, s):
	cs_source = np.empty([501, 4])
	cs_target = np.empty([501, 4])
	lms_to_xyz = np.array([
		[1.91020, -1.11212, 0.20191],
		[0.37095, 0.62905, 0],
		[0, 0, 1]
	])
	
	for i in range(0, 500):
		cs_source[i][0] = wl_to_s(i + 300, l0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		cs_source[i][1] = wl_to_s(i + 300, m0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		cs_source[i][2] = wl_to_s(i + 300, s0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		lms_array = np.empty([3, 1])
		lms_array[0] = cs_source[i][0]
		lms_array[1] = cs_source[i][1]
		lms_array[2] = cs_source[i][2]
		xyz_array = np.matmul(lms_to_xyz, lms_array)
		xyz = XYZColor(xyz_array[0], xyz_array[1], xyz_array[2])
		hsl = convert_color(xyz, HSLColor)
		cs_source[i][3] = hsl.get_value_tuple()[0]
		
		cs_target[i][0] = wl_to_s(i + 300, l) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][1] = wl_to_s(i + 300, m) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][2] = wl_to_s(i + 300, s) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		lms_array1 = np.empty([3, 1])
		lms_array1[0] = cs_target[i][0]
		lms_array1[1] = cs_target[i][1]
		lms_array1[2] = cs_target[i][2]
		xyz_array1 = np.matmul(lms_to_xyz, lms_array1)
		xyz1 = XYZColor(xyz_array1[0], xyz_array1[1], xyz_array1[2])
		hsl1 = convert_color(xyz1, HSLColor)
		cs_target[i][3] = hsl1.get_value_tuple()[0]
	
	return [cs_source, cs_target]

# As with version 1, this one does some intense gap-filling to ensure the result looks sensible.
def hue_to_wavelength_3(h, source):
	if (0 <= h <= 270):
		i = 0
		while (i < 500):
			h_cur = round(source[i][3])
			h_next = round(source[i + 1][3])
			
			# exact match
			if (round(h) == round(h_cur)):
				return i + 300
			# intermediate match
			elif (h_next <= h <= h_cur):
				return i + 1 - (h - h_next) / (h_cur - h_next) + 300
			else:
				i += 1
		
		return 800
	else:
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)

def wl_to_rgb_5(w, target):
	if (w == round(w)):
		return colorsys.hls_to_rgb(target[round(w) - 300][3] / 360, 0.5, 1)
	else: # if non-integer value, use a weighted average
		d = w - int(w)
		return colorsys.hls_to_rgb((target[math.floor(w) - 300][3] * (1 - d) + target[math.ceil(w) - 300][3] * d) / 360, 0.5, 1)

# An XYZ version of the saturation function. The shift is similar to CVS and smoother than
# version 4 but less noticeable. Averaging the "white" LMS values or turning them into
# ratios produces a much more extreme shift. The issue may be E vs. D65.
def saturation_1(w, source, target, l, m, s):
	if (type(w) != tuple):
		if (w == round(w)):
			sx = source[round(w) - 300][0]
			sy = source[round(w) - 300][1]
			sz = source[round(w) - 300][2]
			tx = target[round(w) - 300][0]
			ty = target[round(w) - 300][1]
			tz = target[round(w) - 300][2]
		else:
			sx = (source[round(w) - 301][0] + source[round(w) - 300][0]) / 2
			sy = (source[round(w) - 301][1] + source[round(w) - 300][1]) / 2
			sz = (source[round(w) - 301][2] + source[round(w) - 300][2]) / 2
			tx = (source[round(w) - 301][0] + target[round(w) - 300][0]) / 2
			ty = (source[round(w) - 301][1] + target[round(w) - 300][1]) / 2
			tz = (source[round(w) - 301][2] + target[round(w) - 300][2]) / 2
	else:
		sx = source[r0 - 300][0] * w[0] + source[b0 - 300][0] * w[1]
		sy = source[r0 - 300][1] * w[0] + source[b0 - 300][1] * w[1]
		sz = source[r0 - 300][2] * w[0] + source[b0 - 300][2] * w[1]
		tx = target[r0 - 300][0] * w[0] + target[b0 - 300][0] * w[1]
		ty = target[r0 - 300][1] * w[0] + target[b0 - 300][1] * w[1]
		tz = target[r0 - 300][2] * w[0] + target[b0 - 300][2] * w[1]
	
	# source white point
	white_l = 0
	white_m = 0
	white_s = 0
	for i in range(300, 800):
		white_l += wl_to_s(i, l0)
		white_m += wl_to_s(i, m0)
		white_s += wl_to_s(i, s0)
	
	#white_l = white_l / 500
	#white_m = white_m / 500
	#white_s = white_s / 500
	#white_total = white_l + white_m + white_s
	#white_l = white_l / white_total
	#white_m = white_m / white_total
	#white_s = white_s / white_total
	
	# target white point
	white_l1 = 0
	white_m1 = 0
	white_s1 = 0
	for i in range(300, 800):
		white_l1 += wl_to_s(i, l)
		white_m1 += wl_to_s(i, m)
		white_s1 += wl_to_s(i, s)
	
	#white_l1 = white_l1 / 500
	#white_m1 = white_m1 / 500
	#white_s1 = white_s1 / 500
	#white_total1 = white_l1 + white_m1 + white_s1
	#white_l1 = white_l1 / white_total1
	#white_m1 = white_m1 / white_total1
	#white_s1 = white_s1 / white_total1
	
	# convert to XYZ
	lms_to_xyz = np.array([
		[1.91020, -1.11212, 0.20191],
		[0.37095, 0.62905, 0],
		[0, 0, 1]
	])
	white_lms = np.array([
		[white_l],
		[white_m],
		[white_s]
	])
	white_lms1 = np.array([
		[white_l1],
		[white_m1],
		[white_s1]
	])
	s = np.array([
		[sx],
		[sy],
		[sz]
	])
	t = np.array([
		[tx],
		[ty],
		[tz]
	])
	white_xyz = np.matmul(lms_to_xyz, white_lms)
	white_xyz1 = np.matmul(lms_to_xyz, white_lms1)
	s_xyz = np.matmul(lms_to_xyz, s)
	t_xyz = np.matmul(lms_to_xyz, t)
	
	sat_s = math.dist(white_xyz, s_xyz)
	sat_t = math.dist(white_xyz1, t_xyz)
	if (sat_s == 0):
		return 0
	return sat_t / sat_s

# switch method used based on input -- this function calls the desired version of wl_to_rgb
# default is 2 (cone response ratios)

# save some time by taking these out of wl_to_rgb so we don't call them for every pixel
if (version == 1):
	tables = cs_tables_1(l1, m1, s1)
elif (version == 2):
	primaries = find_primaries(l1, m1, s1)
elif (version == 3):
	tables = cs_tables(l1, m1, s1)
elif (version == 4):
	tables = cs_tables_2(l1, m1, s1)
elif (version == 5):
	tables = cs_tables_3(l1, m1, s1)

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
	elif (version == 5):
		return wl_to_rgb_5(w, tables[1])
	else:
		return wl_to_rgb_0(w, l, m, s)

def hue_to_wavelength(h):
	if (version == 1):
		return hue_to_wavelength_1(h, tables[0])
	elif (version == 4):
		return hue_to_wavelength_2(h, tables[0])
	elif (version == 5):
		return hue_to_wavelength_3(h, tables[0])
	else:
		return hue_to_wavelength_0(h)

# image processing
if (args.image):
	imagename = args.image
	# Convert image from BGR to HLS. We use BGR because this is how OpenCV reads images.
	# If we use RGB, the output appears normal with settings approximating human vision,
	# but shifting the cones produces the opposite of the expected result.
	img = cv2.imread(imagename)
	img_hls = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
	img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

	# transform hues for each pixel
	for x in range(img.shape[0]):
		for y in range(img.shape[1]):
			# report current pixel to determine whether there's an infinite loop
			#print("pixel: " + str(x) + ", " + str(y))
			# find pixel
			hls = img_hls[x][y]
			#lch = convert_color(sRGBColor(img_rgb[x][y][0], img_rgb[x][y][1], img_rgb[x][y][2]), LCHabColor).get_value_tuple()
			
			# convert hue from 0-180 (OpenCV format) to 0-360
			hue = hls[0]*2
			#hue = lch[2]
			#print("hue")
			#print(hue)
			#print(img_hls[x][y][0]*2)
			#
			# convert hue to predominant wavelength(s)
			wl = hue_to_wavelength(hue)
			
			if (type(wl) != tuple):
				hue_target = wl_to_rgb(wl, l1, m1, s1)
			else:
				hue_r = wl_to_rgb(r0, l1, m1, s1)
				hue_b = wl_to_rgb(b0, l1, m1, s1)
				# sum the amounts of "red" and "blue"
				hue_target = [hue_r[0] * wl[0] + hue_b[0] * wl[2], hue_r[1] * wl[0] + hue_b[1] * wl[2], hue_r[2] * wl[0] + hue_b[2] * wl[2]]
				#print(hue_target)
			
			# convert predominant wavelengths back into a hue
			hls_target = colorsys.rgb_to_hls(hue_target[0], hue_target[1], hue_target[2])
			#lch_target = convert_color(sRGBColor(hue_target[0], hue_target[1], hue_target[2]), LCHabColor).get_value_tuple()
			
			if (version == 4):
				sat_diff = saturation(wl, tables[0], tables[1])
				#sat_diff = 1
			elif (version == 5):
				sat_diff = saturation_1(wl, tables[0], tables[1], l1, m1, s1)
			else:
				sat_diff = hls_target[2]
				#sat_diff = lch_target[1]
			
			# shift hue in our pixel. Colorsys uses 0-1, so we have to convert back to
			# OpenCV format.
			img_hls[x][y] = [hls_target[0]*180, hls[1], hls[2] * sat_diff]
			#lch_final = LCHabColor(lch[0], lch[1] * sat_diff, lch_target[2])
			#img_rgb[x][y] = convert_color(lch_final, sRGBColor).get_value_tuple()

	# convert back to BGR
	img_result = cv2.cvtColor(img_hls, cv2.COLOR_HLS2BGR)
	#img_result = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2BGR)

	# fix brightness
	# I tried using LCh instead of HLS, but it doesn't work the way I expected. This extra step
	# doesn't add much time.
	img_lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
	img_result_lab = cv2.cvtColor(img_result, cv2.COLOR_BGR2LAB)
	for x in range(img.shape[0]):
		for y in range(img.shape[1]):
			img_result_lab[x][y][0] = img_lab[x][y][0]

	img_result = cv2.cvtColor(img_result_lab, cv2.COLOR_LAB2BGR)

	# display result
	cv2.imwrite("colorvisionpy-result.png", img_result)

# print general information
else:
	print("L: " + str(l1) + ", M: " + str(m1) + ", S: " + str(s1))
	print("")
	
	# primary/secondary colors
	primaries = find_primaries(l1, m1, s1)
	print("Estimated primary and secondary wavelengths based on cone response ratios:")
	print("red: " + str(primaries[0]))
	print("yellow: " + str(primaries[1]))
	print("green: " + str(primaries[2]))
	print("cyan: " + str(primaries[3]))
	print("blue: " + str(primaries[4]))
	print("violet: " + str(primaries[5]))
	print("")
	
	# white
	white_l = 0
	white_m = 0
	white_s = 0
	for i in range(300, 800):
		white_l += wl_to_s(i, l1) * sensitivity(i, l1, m1, s1, args.cutoff)
		white_m += wl_to_s(i, m1) * sensitivity(i, l1, m1, s1, args.cutoff)
		white_s += wl_to_s(i, s1) * sensitivity(i, l1, m1, s1, args.cutoff)
	
	white_l = white_l / (white_l + white_m + white_s)
	white_m = white_m / (white_l + white_m + white_s)
	white_s = white_s / (white_l + white_m + white_s)
	
	print("LMS response to white (equal-energy or \"E\"): l=" + str(white_l)
		+ ", m=" + str(white_m) + ", s=" + str(white_s))
	
	white_lms = np.array([
		[white_l],
		[white_m],
		[white_s]
	])
	lms_to_xyz = np.array([
		[1.91020, -1.11212, 0.20191],
		[0.37095, 0.62905, 0],
		[0, 0, 1]
	])
	white_xyz = np.matmul(lms_to_xyz, white_lms)
	xyz = XYZColor(white_xyz[0], white_xyz[1], white_xyz[2], illuminant="e")
	print("Color space coordinates of white (with E illuminant):")
	print("CIE XYZ: " + str(xyz.get_value_tuple()))
	print("sRGB: " + str(convert_color(xyz, sRGBColor).get_value_tuple()))
	print("HSL: " + str(convert_color(xyz, HSLColor).get_value_tuple()))
	print("CIE LAB: " + str(convert_color(xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(xyz, LCHabColor).get_value_tuple()))
	print("")
	
	# non-UV white
	white_l = 0
	white_m = 0
	white_s = 0
	for i in range(400, 800):
		white_l += wl_to_s(i, l1) * sensitivity(i, l1, m1, s1, args.cutoff)
		white_m += wl_to_s(i, m1) * sensitivity(i, l1, m1, s1, args.cutoff)
		white_s += wl_to_s(i, s1) * sensitivity(i, l1, m1, s1, args.cutoff)
	
	#white_l = white_l / 500
	#white_m = white_m / 500
	#white_s = white_s / 500
	
	print("LMS response to non-UV white (cut off at 400 nm): l=" + str(white_l)
		+ ", m=" + str(white_m) + ", s=" + str(white_s))
	
	white_lms = np.array([
		[white_l],
		[white_m],
		[white_s]
	])
	white_matrix = np.matmul(lms_to_xyz, white_lms)
	white_xyz = XYZColor(white_matrix[0], white_matrix[1], white_matrix[2], illuminant="e")
	print("Color space coordinates of non-UV white:")
	print("CIE XYZ: " + str(white_xyz.get_value_tuple()))
	print("sRGB: " + str(convert_color(white_xyz, sRGBColor).get_value_tuple()))
	print("HSL: " + str(convert_color(white_xyz, HSLColor).get_value_tuple()))
	print("CIE LAB: " + str(convert_color(white_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(white_xyz, LCHabColor).get_value_tuple()))
	print("")
	
	# leaves
	
	# reflectance spectrum of green maple leaves estimated from this graph:
	# https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgKwT9GhKlvdbOROD_JtEEym7Ovzq8GPfSCORiVcKklEQI9DSTuNjoaIJMWMdpJpc4ijq0T1m_PXF2seWczauKLz-4VPIY9TSXQqXdp1B80vu4w5O4lWaAF0k2kaA5ThrJrjCzlck8Ez1fF/s1600/LeafSpectra.jpg
	# The graph has information for 300-400 but appears to be roughly 0 there.
	leaf_l = 0
	leaf_m = 0
	leaf_s = 0
	leaf_l += wl_to_s(400, l1) * sensitivity(400, l1, m1, s1, args.cutoff) * 0.01
	leaf_m += wl_to_s(400, m1) * sensitivity(400, l1, m1, s1, args.cutoff) * 0.01
	leaf_s += wl_to_s(400, s1) * sensitivity(400, l1, m1, s1, args.cutoff) * 0.01
	leaf_l += wl_to_s(410, l1) * sensitivity(410, l1, m1, s1, args.cutoff) * 0.01
	leaf_m += wl_to_s(410, m1) * sensitivity(410, l1, m1, s1, args.cutoff) * 0.01
	leaf_s += wl_to_s(410, s1) * sensitivity(410, l1, m1, s1, args.cutoff) * 0.01
	leaf_l += wl_to_s(420, l1) * sensitivity(420, l1, m1, s1, args.cutoff) * 0.01
	leaf_m += wl_to_s(420, m1) * sensitivity(420, l1, m1, s1, args.cutoff) * 0.01
	leaf_s += wl_to_s(420, s1) * sensitivity(420, l1, m1, s1, args.cutoff) * 0.01
	leaf_l += wl_to_s(430, l1) * sensitivity(430, l1, m1, s1, args.cutoff) * 0.02
	leaf_m += wl_to_s(430, m1) * sensitivity(430, l1, m1, s1, args.cutoff) * 0.02
	leaf_s += wl_to_s(430, s1) * sensitivity(430, l1, m1, s1, args.cutoff) * 0.02
	leaf_l += wl_to_s(440, l1) * sensitivity(440, l1, m1, s1, args.cutoff) * 0.02
	leaf_m += wl_to_s(440, m1) * sensitivity(440, l1, m1, s1, args.cutoff) * 0.02
	leaf_s += wl_to_s(440, s1) * sensitivity(440, l1, m1, s1, args.cutoff) * 0.02
	leaf_l += wl_to_s(450, l1) * sensitivity(450, l1, m1, s1, args.cutoff) * 0.02
	leaf_m += wl_to_s(450, m1) * sensitivity(450, l1, m1, s1, args.cutoff) * 0.02
	leaf_s += wl_to_s(450, s1) * sensitivity(450, l1, m1, s1, args.cutoff) * 0.02
	leaf_l += wl_to_s(460, l1) * sensitivity(460, l1, m1, s1, args.cutoff) * 0.02
	leaf_m += wl_to_s(460, m1) * sensitivity(460, l1, m1, s1, args.cutoff) * 0.02
	leaf_s += wl_to_s(460, s1) * sensitivity(460, l1, m1, s1, args.cutoff) * 0.02
	leaf_l += wl_to_s(470, l1) * sensitivity(470, l1, m1, s1, args.cutoff) * 0.03
	leaf_m += wl_to_s(470, m1) * sensitivity(470, l1, m1, s1, args.cutoff) * 0.03
	leaf_s += wl_to_s(470, s1) * sensitivity(470, l1, m1, s1, args.cutoff) * 0.03
	leaf_l += wl_to_s(480, l1) * sensitivity(480, l1, m1, s1, args.cutoff) * 0.04
	leaf_m += wl_to_s(480, m1) * sensitivity(480, l1, m1, s1, args.cutoff) * 0.04
	leaf_s += wl_to_s(480, s1) * sensitivity(480, l1, m1, s1, args.cutoff) * 0.04
	leaf_l += wl_to_s(490, l1) * sensitivity(490, l1, m1, s1, args.cutoff) * 0.04
	leaf_m += wl_to_s(490, m1) * sensitivity(490, l1, m1, s1, args.cutoff) * 0.04
	leaf_s += wl_to_s(490, s1) * sensitivity(490, l1, m1, s1, args.cutoff) * 0.04
	leaf_l += wl_to_s(500, l1) * sensitivity(500, l1, m1, s1, args.cutoff) * 0.05
	leaf_m += wl_to_s(500, m1) * sensitivity(500, l1, m1, s1, args.cutoff) * 0.05
	leaf_s += wl_to_s(500, s1) * sensitivity(500, l1, m1, s1, args.cutoff) * 0.05
	leaf_l += wl_to_s(510, l1) * sensitivity(510, l1, m1, s1, args.cutoff) * 0.08
	leaf_m += wl_to_s(510, m1) * sensitivity(510, l1, m1, s1, args.cutoff) * 0.08
	leaf_s += wl_to_s(510, s1) * sensitivity(510, l1, m1, s1, args.cutoff) * 0.08
	leaf_l += wl_to_s(520, l1) * sensitivity(520, l1, m1, s1, args.cutoff) * 0.15
	leaf_m += wl_to_s(520, m1) * sensitivity(520, l1, m1, s1, args.cutoff) * 0.15
	leaf_s += wl_to_s(520, s1) * sensitivity(520, l1, m1, s1, args.cutoff) * 0.15
	leaf_l += wl_to_s(530, l1) * sensitivity(530, l1, m1, s1, args.cutoff) * 0.2
	leaf_m += wl_to_s(530, m1) * sensitivity(530, l1, m1, s1, args.cutoff) * 0.2
	leaf_s += wl_to_s(530, s1) * sensitivity(530, l1, m1, s1, args.cutoff) * 0.2
	leaf_l += wl_to_s(540, l1) * sensitivity(540, l1, m1, s1, args.cutoff) * 0.22
	leaf_m += wl_to_s(540, m1) * sensitivity(540, l1, m1, s1, args.cutoff) * 0.22
	leaf_s += wl_to_s(540, s1) * sensitivity(540, l1, m1, s1, args.cutoff) * 0.22
	leaf_l += wl_to_s(550, l1) * sensitivity(550, l1, m1, s1, args.cutoff) * 0.26
	leaf_m += wl_to_s(550, m1) * sensitivity(550, l1, m1, s1, args.cutoff) * 0.26
	leaf_s += wl_to_s(550, s1) * sensitivity(550, l1, m1, s1, args.cutoff) * 0.26
	leaf_l += wl_to_s(560, l1) * sensitivity(560, l1, m1, s1, args.cutoff) * 0.25
	leaf_m += wl_to_s(560, m1) * sensitivity(560, l1, m1, s1, args.cutoff) * 0.25
	leaf_s += wl_to_s(560, s1) * sensitivity(560, l1, m1, s1, args.cutoff) * 0.25
	leaf_l += wl_to_s(570, l1) * sensitivity(570, l1, m1, s1, args.cutoff) * 0.2
	leaf_m += wl_to_s(570, m1) * sensitivity(570, l1, m1, s1, args.cutoff) * 0.2
	leaf_s += wl_to_s(570, s1) * sensitivity(570, l1, m1, s1, args.cutoff) * 0.2
	leaf_l += wl_to_s(580, l1) * sensitivity(580, l1, m1, s1, args.cutoff) * 0.19
	leaf_m += wl_to_s(580, m1) * sensitivity(580, l1, m1, s1, args.cutoff) * 0.19
	leaf_s += wl_to_s(580, s1) * sensitivity(580, l1, m1, s1, args.cutoff) * 0.19
	leaf_l += wl_to_s(590, l1) * sensitivity(590, l1, m1, s1, args.cutoff) * 0.18
	leaf_m += wl_to_s(590, m1) * sensitivity(590, l1, m1, s1, args.cutoff) * 0.18
	leaf_s += wl_to_s(590, s1) * sensitivity(590, l1, m1, s1, args.cutoff) * 0.18
	leaf_l += wl_to_s(600, l1) * sensitivity(600, l1, m1, s1, args.cutoff) * 0.17
	leaf_m += wl_to_s(600, m1) * sensitivity(600, l1, m1, s1, args.cutoff) * 0.17
	leaf_s += wl_to_s(600, s1) * sensitivity(600, l1, m1, s1, args.cutoff) * 0.17
	leaf_l += wl_to_s(610, l1) * sensitivity(610, l1, m1, s1, args.cutoff) * 0.15
	leaf_m += wl_to_s(610, m1) * sensitivity(610, l1, m1, s1, args.cutoff) * 0.15
	leaf_s += wl_to_s(610, s1) * sensitivity(610, l1, m1, s1, args.cutoff) * 0.15
	leaf_l += wl_to_s(620, l1) * sensitivity(620, l1, m1, s1, args.cutoff) * 0.15
	leaf_m += wl_to_s(620, m1) * sensitivity(620, l1, m1, s1, args.cutoff) * 0.15
	leaf_s += wl_to_s(620, s1) * sensitivity(620, l1, m1, s1, args.cutoff) * 0.15
	leaf_l += wl_to_s(630, l1) * sensitivity(630, l1, m1, s1, args.cutoff) * 0.14
	leaf_m += wl_to_s(630, m1) * sensitivity(630, l1, m1, s1, args.cutoff) * 0.14
	leaf_s += wl_to_s(630, s1) * sensitivity(630, l1, m1, s1, args.cutoff) * 0.14
	leaf_l += wl_to_s(640, l1) * sensitivity(640, l1, m1, s1, args.cutoff) * 0.13
	leaf_m += wl_to_s(640, m1) * sensitivity(640, l1, m1, s1, args.cutoff) * 0.13
	leaf_s += wl_to_s(640, s1) * sensitivity(640, l1, m1, s1, args.cutoff) * 0.13
	leaf_l += wl_to_s(650, l1) * sensitivity(650, l1, m1, s1, args.cutoff) * 0.09
	leaf_m += wl_to_s(650, m1) * sensitivity(650, l1, m1, s1, args.cutoff) * 0.09
	leaf_s += wl_to_s(650, s1) * sensitivity(650, l1, m1, s1, args.cutoff) * 0.09
	leaf_l += wl_to_s(660, l1) * sensitivity(660, l1, m1, s1, args.cutoff) * 0.05
	leaf_m += wl_to_s(660, m1) * sensitivity(660, l1, m1, s1, args.cutoff) * 0.05
	leaf_s += wl_to_s(660, s1) * sensitivity(660, l1, m1, s1, args.cutoff) * 0.05
	leaf_l += wl_to_s(670, l1) * sensitivity(670, l1, m1, s1, args.cutoff) * 0.04
	leaf_m += wl_to_s(670, m1) * sensitivity(670, l1, m1, s1, args.cutoff) * 0.04
	leaf_s += wl_to_s(670, s1) * sensitivity(670, l1, m1, s1, args.cutoff) * 0.04
	leaf_l += wl_to_s(680, l1) * sensitivity(680, l1, m1, s1, args.cutoff) * 0.05
	leaf_m += wl_to_s(680, m1) * sensitivity(680, l1, m1, s1, args.cutoff) * 0.05
	leaf_s += wl_to_s(680, s1) * sensitivity(680, l1, m1, s1, args.cutoff) * 0.05
	leaf_l += wl_to_s(690, l1) * sensitivity(690, l1, m1, s1, args.cutoff) * 0.09
	leaf_m += wl_to_s(690, m1) * sensitivity(690, l1, m1, s1, args.cutoff) * 0.09
	leaf_s += wl_to_s(690, s1) * sensitivity(690, l1, m1, s1, args.cutoff) * 0.09
	leaf_l += wl_to_s(700, l1) * sensitivity(700, l1, m1, s1, args.cutoff) * 0.27
	leaf_m += wl_to_s(700, m1) * sensitivity(700, l1, m1, s1, args.cutoff) * 0.27
	leaf_s += wl_to_s(700, s1) * sensitivity(700, l1, m1, s1, args.cutoff) * 0.27
	leaf_l += wl_to_s(710, l1) * sensitivity(710, l1, m1, s1, args.cutoff) * 0.35
	leaf_m += wl_to_s(710, m1) * sensitivity(710, l1, m1, s1, args.cutoff) * 0.35
	leaf_s += wl_to_s(710, s1) * sensitivity(710, l1, m1, s1, args.cutoff) * 0.35
	leaf_l += wl_to_s(720, l1) * sensitivity(720, l1, m1, s1, args.cutoff) * 0.4
	leaf_m += wl_to_s(720, m1) * sensitivity(720, l1, m1, s1, args.cutoff) * 0.4
	leaf_s += wl_to_s(720, s1) * sensitivity(720, l1, m1, s1, args.cutoff) * 0.4
	leaf_l += wl_to_s(730, l1) * sensitivity(730, l1, m1, s1, args.cutoff) * 0.5
	leaf_m += wl_to_s(730, m1) * sensitivity(730, l1, m1, s1, args.cutoff) * 0.5
	leaf_s += wl_to_s(730, s1) * sensitivity(730, l1, m1, s1, args.cutoff) * 0.5
	leaf_l += wl_to_s(740, l1) * sensitivity(740, l1, m1, s1, args.cutoff) * 0.55
	leaf_m += wl_to_s(740, m1) * sensitivity(740, l1, m1, s1, args.cutoff) * 0.55
	leaf_s += wl_to_s(740, s1) * sensitivity(740, l1, m1, s1, args.cutoff) * 0.55
	leaf_l += wl_to_s(750, l1) * sensitivity(750, l1, m1, s1, args.cutoff) * 0.6
	leaf_m += wl_to_s(750, m1) * sensitivity(750, l1, m1, s1, args.cutoff) * 0.6
	leaf_s += wl_to_s(750, s1) * sensitivity(750, l1, m1, s1, args.cutoff) * 0.6
	
	print("LMS response to green maple leaf: l=" + str(leaf_l) +", m=" + str(leaf_m) + ", s=" + str(leaf_s))
	leaf_lms = np.array([
		[leaf_l],
		[leaf_m],
		[leaf_s]
	])
	leaf_matrix = np.matmul(lms_to_xyz, leaf_lms)
	leaf_xyz = XYZColor(leaf_matrix[0], leaf_matrix[1], leaf_matrix[2], illuminant="d65")
	print("Color coordinates of green maple leaf (with D65):")
	print("CIE XYZ: " + str(leaf_xyz.get_value_tuple()))
	print("sRGB: " + str(convert_color(leaf_xyz, sRGBColor).get_value_tuple()))
	print("HSL: " + str(convert_color(leaf_xyz, HSLColor).get_value_tuple()))
	print("CIE LAB: " + str(convert_color(leaf_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(leaf_xyz, LCHabColor).get_value_tuple()))
	print("")

# print execution time
print("%s seconds" % (time.time() - start_time))
