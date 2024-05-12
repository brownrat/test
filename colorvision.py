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
from colormath.color_objects import LabColor, LCHabColor, LCHuvColor, sRGBColor, XYZColor, HSLColor, SpectralColor, IPTColor
from colormath.color_conversions import convert_color
from colormath import spectral_constants

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
def lens_filter(w, c=300):
	value = np.log(w - c) / np.log(800 - c)
	if (value > 0): # a log function has the wrong shape near the cutoff point
		return value
	return 0

def sensitivity(w, l, m, s, c=300):
	value = (wl_to_s(w, l) + wl_to_s(w, m) / 2 + wl_to_s(w, s) / 10) * lens_filter(w, c)
	return value

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

def rgyb_to_rgb(rg, yb, red=(0,0), blue=(0,0)):
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
		i = 0
		while (i < 500):
			h_cur = round(source[i][2])
			h_next = round(source[i + 1][2])
			
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

def wl_to_rgb_4(w, target):
	if (w == round(w)):
		return colorsys.hls_to_rgb(target[round(w) - 300][2] / 360, 0.5, 1)
	else: # if non-integer value, use a weighted average
		d = w - int(w)
		return colorsys.hls_to_rgb((target[math.floor(w) - 300][2] * (1 - d) + target[math.ceil(w) - 300][2] * d) / 360, 0.5, 1)

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
# This converts LMS values to CIE XYZ and then to HLS. See https://en.wikipedia.org/wiki/CIE_1931_color_space#Meaning_of_X,_Y_and_Z for the conversion matrix.
# The results are very similar to CVS and to version 4 aside from a slight green shift for
# narrow spacing, probably for the same reason as in version 2. The main problem is the
# range of hues provided for integer wavelengths is full of holes you could drive a truck
# through. These gaps can be filled in convincingly with some careful calculations (probably
# the reason for the long execution time) for primate-like values, but for other arrangements
# it doesn't work as well.
# XYZ conversion can also simulate dichromacy in the same way as v0 and v4. Tritanopia is a
# bit odd-looking as red becomes purple for some reason.
# Looking at what's actually in the tables this generates reveals that the conversion from
# XYZ to sRGB/HSL is extremely wrong and it only looks good because it goes back the same
# way it came in. I don't know why.
lms_to_xyz = np.array([
	[1.94735469, -1.41445123, 0.36476327],
	[0.68990272, 0.34832189, 0],
	[0, 0, 1.93485343]
])

def cs_tables_3(l, m, s):
	cs_source = np.empty([501, 4])
	cs_target = np.empty([501, 4])
#	lms_to_xyz = np.array([
#		[1.91020, -1.11212, 0.20191],
#		[0.37095, 0.62905, 0],
#		[0, 0, 1]
#	])
	
	for i in range(0, 500):
		cs_source[i][0] = wl_to_s(i + 300, l0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		cs_source[i][1] = wl_to_s(i + 300, m0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		cs_source[i][2] = wl_to_s(i + 300, s0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		lms_array = np.empty([3, 1])
		lms_array[0] = cs_source[i][0]
		lms_array[1] = cs_source[i][1]
		lms_array[2] = cs_source[i][2]
		xyz_array = np.matmul(lms_to_xyz, lms_array)
		xyz = XYZColor(*xyz_array)
		hsl = convert_color(xyz, HSLColor)
		cs_source[i][3] = hsl.get_value_tuple()[0]
		#lch = convert_color(xyz, LCHabColor)
		#cs_source[i][3] = lch.get_value_tuple()[2]
		
		cs_target[i][0] = wl_to_s(i + 300, l) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][1] = wl_to_s(i + 300, m) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][2] = wl_to_s(i + 300, s) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		lms_array1 = np.empty([3, 1])
		lms_array1[0] = cs_target[i][0]
		lms_array1[1] = cs_target[i][1]
		lms_array1[2] = cs_target[i][2]
		xyz_array1 = np.matmul(lms_to_xyz, lms_array1)
		xyz1 = XYZColor(*xyz_array1)
		hsl1 = convert_color(xyz1, HSLColor)
		cs_target[i][3] = hsl1.get_value_tuple()[0]
		#lch1 = convert_color(xyz1, LCHabColor)
		#cs_target[i][3] = lch1.get_value_tuple()[2]
	
	return [cs_source, cs_target]

# An LCHab version. This produces more reasonable hue values, but image transformations
# look worse. (See below for what's going on with this)
def cs_tables_4(l, m, s):
	cs_source = np.empty([501, 4])
	cs_target = np.empty([501, 4])
#	lms_to_xyz = np.array([
#		[1.91020, -1.11212, 0.20191],
#		[0.37095, 0.62905, 0],
#		[0, 0, 1]
#	])
	
	for i in range(0, 500):
		cs_source[i][0] = wl_to_s(i + 300, l0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		cs_source[i][1] = wl_to_s(i + 300, m0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		cs_source[i][2] = wl_to_s(i + 300, s0) / (wl_to_s(i + 300, l0) + wl_to_s(i + 300, m0) + wl_to_s(i + 300, s0))
		lms_array = np.empty([3, 1])
		lms_array[0] = cs_source[i][0]
		lms_array[1] = cs_source[i][1]
		lms_array[2] = cs_source[i][2]
		xyz_array = np.matmul(lms_to_xyz, lms_array)
		xyz = XYZColor(*xyz_array)
		#hsl = convert_color(xyz, HSLColor)
		#cs_source[i][3] = hsl.get_value_tuple()[0]
		lch = convert_color(xyz, LCHabColor)
		cs_source[i][3] = lch.get_value_tuple()[2]
		
		cs_target[i][0] = wl_to_s(i + 300, l) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][1] = wl_to_s(i + 300, m) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][2] = wl_to_s(i + 300, s) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		lms_array1 = np.empty([3, 1])
		lms_array1[0] = cs_target[i][0]
		lms_array1[1] = cs_target[i][1]
		lms_array1[2] = cs_target[i][2]
		xyz_array1 = np.matmul(lms_to_xyz, lms_array1)
		xyz1 = XYZColor(*xyz_array1)
		#hsl1 = convert_color(xyz1, HSLColor)
		#cs_target[i][3] = hsl1.get_value_tuple()[0]
		lch1 = convert_color(xyz1, LCHabColor)
		cs_target[i][3] = lch1.get_value_tuple()[2]
	
	return [cs_source, cs_target]

# LCHuv version
def cs_tables_5(l, m, s):
	cs_source = np.empty([501, 4])
	cs_target = np.empty([501, 4])
#	lms_to_xyz = np.array([
#		[1.91020, -1.11212, 0.20191],
#		[0.37095, 0.62905, 0],
#		[0, 0, 1]
#	])
	
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
		#hsl = convert_color(xyz, HSLColor)
		#cs_source[i][3] = hsl.get_value_tuple()[0]
		lch = convert_color(xyz, LCHuvColor)
		cs_source[i][3] = lch.get_value_tuple()[2]
		
		cs_target[i][0] = wl_to_s(i + 300, l) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][1] = wl_to_s(i + 300, m) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		cs_target[i][2] = wl_to_s(i + 300, s) / (wl_to_s(i + 300, l) + wl_to_s(i + 300, m) + wl_to_s(i + 300, s))
		lms_array1 = np.empty([3, 1])
		lms_array1[0] = cs_target[i][0]
		lms_array1[1] = cs_target[i][1]
		lms_array1[2] = cs_target[i][2]
		xyz_array1 = np.matmul(lms_to_xyz, lms_array1)
		xyz1 = XYZColor(xyz_array1[0], xyz_array1[1], xyz_array1[2])
		#hsl1 = convert_color(xyz1, HSLColor)
		#cs_target[i][3] = hsl1.get_value_tuple()[0]
		lch1 = convert_color(xyz1, LCHabColor)
		cs_target[i][3] = lch1.get_value_tuple()[2]
	
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
	
	# maximum sensitivity
	ms = 300
	for i in range(300, 800):
		if (sensitivity(i, l1, m1, s1, args.cutoff) > (sensitivity(ms, l1, m1, s1, args.cutoff))):
			ms = i
	print("Maximum sensitivity: " + str(ms))
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
	
	#table = cs_tables_2(l1, m1, s1)[1]
	#print("Estimated primary and secondary wavelengths based on \"RGYB\" (L-M and LM-S differences):")
	#print("red: " + str(hue_to_wavelength_2(0, table)))
	#print("yellow: " + str(hue_to_wavelength_2(60, table)))
	#print("green: " + str(hue_to_wavelength_2(120, table)))
	#print("cyan: " + str(hue_to_wavelength_2(180, table)))
	#print("blue: " + str(hue_to_wavelength_2(240, table)))
	#print("violet: " + str(hue_to_wavelength_2(270, table)))
	#print("")
	
	table = cs_tables_3(l1, m1, s1)[1]
	
	# red and violet are cutoff points
	i = 0
	while (i < 500):
		h_cur = round(table[i][3])
		h_next = round(table[i + 1][3])
		
		# exact match
		if (0 == round(h_cur)):
			break
		# intermediate match
		elif (h_next <= 0 <= h_cur):
			i = i + 1 - (0 - h_next) / (h_cur - h_next)
			break
		else:
			i += 1
	
	r = i + 300
	
	j = 500
	while (j > 0):
		h_cur = round(table[j][3])
		h_next = round(table[j - 1][3])
		
		# exact match
		if (270 == round(h_cur)):
			break
		# intermediate match
		elif (h_next >= 270 >= h_cur):
			j = j - 1 + (270 - h_next) / (h_cur - h_next)
			break
		else:
			j -= 1
	
	v = j + 300
	
	# yellow through blue may be a range of hues, so we try to find the middle
	y1 = 0
	y2 = 0
	i = 0
	while (i < 500):
		h_cur = round(table[i][3])
		h_next = round(table[i + 1][3])
		
		# exact match
		if (60 == round(h_cur)):
			y1 = i + 300
			break
		# intermediate match
		elif (h_next <= 60 <= h_cur):
			y1 = i + 1 - (60 - h_next) / (h_cur - h_next) + 300
			break
		else:
			i += 1
	
	j = 500
	while (j > 0):
		h_cur = round(table[j][3])
		h_next = round(table[j - 1][3])
		
		# exact match
		if (60 == round(h_cur)):
			y2 = j + 300
			break
		# intermediate match
		elif (h_next >= 60 >= h_cur):
			y2 = j - 1 + (60 - h_next) / (h_cur - h_next) + 300
			break
		else:
			j -= 1
	
	y = (y1 + y2) / 2
	
	g1 = 0
	g2 = 0
	i = 0
	while (i < 500):
		h_cur = round(table[i][3])
		h_next = round(table[i + 1][3])
		
		# exact match
		if (120 == round(h_cur)):
			g1 = i + 300
			break
		# intermediate match
		elif (h_next <= 120 <= h_cur):
			g1 = i + 1 - (120 - h_next) / (h_cur - h_next) + 300
			break
		else:
			i += 1
	
	j = 500
	while (j > 0):
		h_cur = round(table[j][3])
		h_next = round(table[j - 1][3])
		
		# exact match
		if (120 == round(h_cur)):
			g2 = j + 300
			break
		# intermediate match
		elif (h_next >= 120 >= h_cur):
			g2 = j - 1 + (120 - h_next) / (h_cur - h_next) + 300
			break
		else:
			j -= 1
	
	g = (g1 + g2) / 2
	
	c1 = 0
	c2 = 0
	i = 0
	while (i < 500):
		h_cur = round(table[i][3])
		h_next = round(table[i + 1][3])
		
		# exact match
		if (180 == round(h_cur)):
			c1 = i + 300
			break
		# intermediate match
		elif (h_next <= 180 <= h_cur):
			c1 = i + 1 - (180 - h_next) / (h_cur - h_next) + 300
			break
		else:
			i += 1
	
	j = 500
	while (j > 0):
		h_cur = round(table[j][3])
		h_next = round(table[j - 1][3])
		
		# exact match
		if (180 == round(h_cur)):
			c2 = j + 300
			break
		# intermediate match
		elif (h_next >= 180 >= h_cur):
			c2 = j - 1 + (180 - h_next) / (h_cur - h_next) + 300
			break
		else:
			j -= 1
	
	c = (c1 + c2) / 2
	
	b1 = 0
	b2 = 0
	i = 0
	while (i < 500):
		h_cur = round(table[i][3])
		h_next = round(table[i + 1][3])
		
		# exact match
		if (240 == round(h_cur)):
			b1 = i + 300
			break
		# intermediate match
		elif (h_next <= 240 <= h_cur):
			b1 = i + 1 - (240 - h_next) / (h_cur - h_next) + 300
			break
		else:
			i += 1
	
	j = 500
	while (j > 0):
		h_cur = round(table[j][3])
		h_next = round(table[j - 1][3])
		
		# exact match
		if (240 == round(h_cur)):
			b2 = j + 300
			break
		# intermediate match
		elif (h_next >= 240 >= h_cur):
			b2 = j - 1 + (240 - h_next) / (h_cur - h_next) + 300
			break
		else:
			j -= 1
	
	b = (b1 + b2) / 2
	
	print("Estimated primary and secondary wavelengths based on CIE XYZ (converted through sRGB):")
	print("red: " + str(r))
	print("yellow: " + str(y) + " (" + str(y1) + "-" + str(y2) + ")")
	print("green: " + str(g) + " (" + str(g1) + "-" + str(g2) + ")")
	print("cyan: " + str(c) + " (" + str(c1) + "-" + str(c2) + ")")
	print("blue: " + str(b) + " (" + str(b1) + "-" + str(b2) + ")")
	print("violet: " + str(v))
	print("")
	
	#table1 = cs_tables_4(l1, m1, s1)[1]
	#print("Estimated primary and secondary wavelengths based on CIE XYZ (converted through LCHab):")
	#print("red: " + str(hue_to_wavelength_3(0, table1)))
	#print("yellow: " + str(hue_to_wavelength_3(60, table1)))
	#print("green: " + str(hue_to_wavelength_3(120, table1)))
	#print("cyan: " + str(hue_to_wavelength_3(180, table1)))
	#print("blue: " + str(hue_to_wavelength_3(240, table1)))
	#print("violet: " + str(hue_to_wavelength_3(270, table1)))
	#print("")
	
	#table2 = cs_tables_5(l1, m1, s1)[1]
	#print("Estimated primary and secondary wavelengths based on CIE XYZ (converted through LCHuv):")
	#print("red: " + str(hue_to_wavelength_3(0, table2)))
	#print("yellow: " + str(hue_to_wavelength_3(60, table2)))
	#print("green: " + str(hue_to_wavelength_3(120, table2)))
	#print("cyan: " + str(hue_to_wavelength_3(180, table2)))
	#print("blue: " + str(hue_to_wavelength_3(240, table2)))
	#print("violet: " + str(hue_to_wavelength_3(270, table2)))
	#print("")
	
	# estimated chromatic distances based on formula in Arrese et al. (2005) without the
	# noise terms. The idea was to produce something similar to wavelength discrimination
	# functions, which show secondary colors fairly well, but the results are not what
	# one would expect.
	#table1 = np.empty((501, 1))
	#for i in range(0, 500):
	#	la = wl_to_s(i + 299, l1) / (wl_to_s(i + 299, l1) + wl_to_s(i + 299, m1) + wl_to_s(i + 299, s1))
	#	ma = wl_to_s(i + 299, m1) / (wl_to_s(i + 299, l1) + wl_to_s(i + 299, m1) + wl_to_s(i + 299, s1))
	#	sa = wl_to_s(i + 299, s1) / (wl_to_s(i + 299, l1) + wl_to_s(i + 299, m1) + wl_to_s(i + 299, s1))
	#	lb = wl_to_s(i + 300, l1) / (wl_to_s(i + 300, l1) + wl_to_s(i + 300, m1) + wl_to_s(i + 300, s1))
	#	mb = wl_to_s(i + 300, m1) / (wl_to_s(i + 300, l1) + wl_to_s(i + 300, m1) + wl_to_s(i + 300, s1))
	#	sb = wl_to_s(i + 300, s1) / (wl_to_s(i + 300, l1) + wl_to_s(i + 300, m1) + wl_to_s(i + 300, s1))
	#	lc = wl_to_s(i + 301, l1) / (wl_to_s(i + 301, l1) + wl_to_s(i + 301, m1) + wl_to_s(i + 301, s1))
	#	mc = wl_to_s(i + 301, m1) / (wl_to_s(i + 301, l1) + wl_to_s(i + 301, m1) + wl_to_s(i + 301, s1))
	#	sc = wl_to_s(i + 301, s1) / (wl_to_s(i + 301, l1) + wl_to_s(i + 301, m1) + wl_to_s(i + 301, s1))
	#	dl0 = abs(lb - la)
	#	dm0 = abs(mb - ma)
	#	ds0 = abs(sb - sa)
	#	d0 = (dl0 - dm0)**2 + (dl0 - ds0)**2 + (ds0 - dm0)**2
	#	#d0 = math.dist((la, ma, sa), (lb, mb, sb))
	#	dl1 = abs(lc - lb)
	#	dm1 = abs(mc - mb)
	#	ds1 = abs(sc - sb)
	#	d1 = (dl1 - dm1)**2 + (dl1 - ds1)**2 + (ds1 - dm1)**2
	#	#d1 = math.dist((lb, mb, sb), (lc, mc, sc))
	#	table1[i][0] = d0 + d1
		#print(i + 300)
		#print(d0)
		#print(d1)
	
	#d0 = table[0][0]
	#w0 = 300
#	for i in range(s1, m1):
#		if (table[i - 300][0] > d0):
#			d0 = table[i - 300][0]
#			w0 = i
#	
#	d1 = table[m1 - 300][0]
#	w1 = m1 - 300
#	for i in range(m1, 800):
#		if (table[i - 300][0] > d1):
#			d1 = table[i - 300][0]
#			w1 = i
	
	#print("Wavelengths with maximum chromatic separation:")
	#print("L-M: " + str(w1))
	#print("M-S: " + str(w0))
	#print("")
	
	# crossover points
	# These seem like they should correspond to secondary colors but don't really, as
	# discussed earlier. The reason for this is probably that combinations of primary
	# and/or complementary wavelengths have to appear "white" and white light does not
	# excite all three cones equally. Thus "yellow" and "cyan" have L:M and M:S ratios
	# much higher than 1:1. The location of yellow/green also seems to be influenced
	# by overlap with "blue".
	# Crossover points can only tell us this much:
	# * The actual secondary wavelengths are longer.
	# * The less overlap, the closer the secondary wavelengths are to the crossover points.
	
	# "cyan"
	i = s1
	while (i < m1 and wl_to_s(i, s1) / wl_to_s(i, m1) > 1):
		i += 1
	j = m1
	while (j > s1 and wl_to_s(j, s1) / wl_to_s(j, m1) < 1):
		j -= 1
	c = (i + j) / 2
	
	# "yellow"
	i = m1
	while (i < 800 and wl_to_s(i, m1) / wl_to_s(i, l1) > 1):
		i += 1
	j = 800
	while (j > m1 and wl_to_s(j, m1) / wl_to_s(j, l1) < 1):
		j -= 1
	y = (i + j) / 2
	
	# blue/violet
	i = 300
	while (i < 800 and wl_to_s(i, m1) / wl_to_s(i, l1) > 1):
		i += 1
	j = s1
	while (j > m1 and wl_to_s(j, m1) / wl_to_s(j, l1) < 1):
		j -= 1
	b = (i + j) / 2
	
	print("Crossover points:")
	print("L-M 1: " + str(y))
	print("L-M 2: " + str(b))
	print("M-S: " + str(c))
	print("")
	
	# white
	# HSL comes up with a negative saturation. Either Python expects something different
	# from XYZ input or it just doesn't handle color spaces properly. Also note that for
	# both the "white" and "maple leaf" colors, the XYZ, LAB and LCh values specify an
	# impossibly high lightness/intensity (>100) whereas sRGB and HSL (derived from sRGB)
	# specify a very low value. (No, it's out of 1, not 100)
	white_l = 0
	white_m = 0
	white_s = 0
	for i in range(300, 800):
		white_l += wl_to_s(i, l1)
		white_m += wl_to_s(i, m1)
		white_s += wl_to_s(i, s1)
	
	print("LMS response to white (equal-energy or \"E\"): l=" + str(white_l)
		+ ", m=" + str(white_m) + ", s=" + str(white_s))
	
	white_lms = np.array([
		[white_l],
		[white_m],
		[white_s]
	])
	# LMS->XYZ matrix (D65?) taken from Wikipedia
	#lms_to_xyz = np.array([
	#	[1.91020, -1.11212, 0.20191],
	#	[0.37095, 0.62905, 0],
	#	[0, 0, 1]
	#])
	
	# Python's matrix from IPTColor. Very similar to the above, and may be wrong as
	# it expects D65.
	#lms_to_xyz = np.linalg.inv(IPTColor.conversion_matrices['xyz_to_lms'])
	#print(IPTColor.conversion_matrices['xyz_to_lms'])
	#print(np.linalg.inv(IPTColor.conversion_matrices['xyz_to_lms']))
	
	# inverse of XYZ->LMS (E) from Wikipedia
	# This is even more wrong. On inspection, it multiplies L by a larger number and
	# M and S by smaller (magnitude) numbers.
	#lms_to_xyz = np.linalg.inv([
	#	[0.38971, 0.68898, -0.07868],
	#	[-0.22981, 1.18340, 0.04641],
	#	[0, 0, 1]
	#])
	#print(np.linalg.inv([
	#	[0.38971, 0.68898, -0.07868],
	#	[-0.22981, 1.18340, 0.04641],
	#	[0, 0, 1]
	#]))
	# another Wikipedia matrix
	lms_to_xyz = np.array([
		[1.94735469, -1.41445123, 0.3647327],
		[0.68990272, 0.34832189, 0],
		[0, 0, 1.93485343]
	])
	
	#xyz_to_rgb = (1 / 3400850) * np.array([
	#	[8041697, -3049000, -1591847],
	#	[-1752003, 4851000, 301853],
	#	[17697, -49000, 3432153]
	#])
	white_matrix = np.matmul(lms_to_xyz, white_lms)
	#white_matrix1 = np.matmul(xyz_to_rgb, white_matrix)
	white_xyz = XYZColor(*white_matrix, illuminant="e")
	#white_rgb = sRGBColor(white_matrix1[0], white_matrix1[1], white_matrix1[2])
	print("Color space coordinates of white (with E illuminant):")
	print("CIE XYZ: " + str(white_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(white_rgb.get_value_tuple()))
	#print("sRGB: " + str(convert_color(white_xyz, sRGBColor).get_value_tuple()))
	#print("HSL (converted through sRGB): " + str(convert_color(white_xyz, HSLColor).get_value_tuple()))
	#print("HSL (converted through CIE RGB): " + str(convert_color(white_rgb, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(white_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(white_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(white_xyz, LCHuvColor).get_value_tuple()))
	print("")
	
	# non-UV white
	white_l = 0
	white_m = 0
	white_s = 0
	for i in range(400, 800):
		white_l += wl_to_s(i, l1)
		white_m += wl_to_s(i, m1)
		white_s += wl_to_s(i, s1)
	
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
	#white_matrix1 = np.matmul(xyz_to_rgb, white_matrix)
	white_xyz = XYZColor(*white_matrix, illuminant="e")
	#white_rgb = sRGBColor(white_matrix1[0], white_matrix1[1], white_matrix1[2])
	print("Color space coordinates of non-UV white:")
	print("CIE XYZ: " + str(white_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(white_rgb.get_value_tuple()))
	#print("sRGB: " + str(convert_color(white_xyz, sRGBColor).get_value_tuple()))
	#print("HSL (converted through sRGB): " + str(convert_color(white_xyz, HSLColor).get_value_tuple()))
	#print("HSL (converted through CIE RGB): " + str(convert_color(white_rgb, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(white_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(white_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(white_xyz, LCHuvColor).get_value_tuple()))
	print("")
	
	# leaves
	# Again, something is wrong with the XYZ conversions. The LCHab hues kind of make sense
	# if you interpret them as HSL/HSV hues, but they're not supposed to work that way --
	# the hue varies with lightness. This means the three leaf spectra are interpreted as
	# "yellow". This is probably why using LCH rather than HSL for image transformations
	# totally mangles the colors.
	# The output is somewhat better when using Python's built-in XYZ<->LMS matrix, E instead
	# of D65 and 440-540-565 instead of 420-530-560, but it's still wrong.
	# Using the "physiological" matrix from Wikipedia produces a much better result as it
	# increases the contribution from M and S, though now somewhat green-shifted. I'm unsure
	# if the luminous efficiency function should be included. Is the absorption ratio
	# affected by the relative sensitivity and density of cone types, or are they
	# adjusted relative to each other?
	
	# I finally found D65.
	d65 = spectral_constants.REF_ILLUM_TABLE["d65"]
	
	# reflectance spectrum of maple leaves estimated from this graph:
	# https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgKwT9GhKlvdbOROD_JtEEym7Ovzq8GPfSCORiVcKklEQI9DSTuNjoaIJMWMdpJpc4ijq0T1m_PXF2seWczauKLz-4VPIY9TSXQqXdp1B80vu4w5O4lWaAF0k2kaA5ThrJrjCzlck8Ez1fF/s1600/LeafSpectra.jpg
	# The graph has information for 300-400 but appears to be roughly 0 there.
	leaf_table = np.array([
		[340, 0.0],
		[350, 0.0],
		[360, 0.0],
		[370, 0.0],
		[380, 0.0],
		[390, 0.0],
		[400, 0.01],
		[410, 0.01],
		[420, 0.02],
		[430, 0.02],
		[440, 0.02],
		[450, 0.03],
		[460, 0.03],
		[470, 0.03],
		[480, 0.04],
		[490, 0.04],
		[500, 0.05],
		[510, 0.08],
		[520, 0.19],
		[530, 0.23],
		[540, 0.25],
		[550, 0.27],
		[560, 0.25],
		[570, 0.23],
		[580, 0.19],
		[590, 0.17],
		[600, 0.17],
		[610, 0.15],
		[620, 0.14],
		[630, 0.14],
		[640, 0.1],
		[650, 0.07],
		[660, 0.05],
		[670, 0.04],
		[680, 0.1],
		[690, 0.3],
		[700, 0.35],
		[710, 0.45],
		[720, 0.55],
		[730, 0.6],
		[740, 0.63],
		[750, 0.65]
	])
	leaf_l = 0
	leaf_m = 0
	leaf_s = 0
	for i in range(0, leaf_table.shape[0]):
		f = lens_filter(leaf_table[i][0], args.cutoff)
		leaf_l += wl_to_s(leaf_table[i][0], l1) * leaf_table[i][1] * f
		leaf_m += wl_to_s(leaf_table[i][0], m1) * leaf_table[i][1] * f
		leaf_s += wl_to_s(leaf_table[i][0], s1) * leaf_table[i][1] * f
	
	# ratio comparison
	similar = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, l0) / wl_to_s(i, m0), 2)
		ratio1 = round(leaf_l / leaf_m, 2)
		ratio2 = round(wl_to_s(i + 1, l0) / wl_to_s(i + 1, m0), 2)
		if (ratio0 == ratio1):
			similar = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar = i + 0.5
	
	similar1 = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, m0) / wl_to_s(i, s0), 2)
		ratio1 = round(leaf_m / leaf_s, 2)
		ratio2 = round(wl_to_s(i + 1, m0) / wl_to_s(i + 1, s0), 2)
		if (ratio0 == ratio1):
			similar1 = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar1 = i + 0.5
	
	# normalize
	total = leaf_l + leaf_m + leaf_s
	leaf_l = leaf_l / total
	leaf_m = leaf_m / total
	leaf_s = leaf_s / total
	
	print("LMS response to green maple leaf: l=" + str(leaf_l) +", m=" + str(leaf_m) + ", s=" + str(leaf_s))
	print("L:M ratio: " + str(leaf_l / leaf_m) + " (looks like " + str(similar) + " nm)")
	print("M:S ratio: " + str(leaf_m / leaf_s) + " (looks like " + str(similar1) + " nm)")
	
	leaf_lms = np.array([
		[leaf_l],
		[leaf_m],
		[leaf_s]
	])
	
	leaf_matrix = np.matmul(lms_to_xyz, leaf_lms)
	#leaf_matrix1 = np.matmul(xyz_to_rgb, leaf_matrix)
	leaf_xyz = XYZColor(*leaf_matrix, illuminant="e")
	#leaf_rgb = sRGBColor(leaf_matrix1[0], leaf_matrix1[1], leaf_matrix1[2])
	print("Color coordinates of green maple leaf:")
	print("CIE XYZ: " + str(leaf_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(leaf_rgb.get_value_tuple()))
	#print("sRGB: " + str(convert_color(leaf_xyz, sRGBColor).get_value_tuple()))
	#print("RGB (converted through RGYB): " + str(rgyb_to_rgb((leaf_l - leaf_m) / (leaf_l + leaf_m + leaf_s), (leaf_l + leaf_m - leaf_s) / (leaf_l + leaf_m + leaf_s))))
	print("HSL: " + str(convert_color(leaf_xyz, HSLColor).get_value_tuple()))
	#print("HSL (converted through CIE RGB): " + str(convert_color(leaf_rgb, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(leaf_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(leaf_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(leaf_xyz, LCHuvColor).get_value_tuple()))
	
	# "SpectralColor" object
	leaf_spectral = SpectralColor(spec_340nm=leaf_table[0][1],
	spec_350nm=leaf_table[1][1],
	spec_360nm=leaf_table[2][1],
	spec_370nm=leaf_table[3][1],
	spec_380nm=leaf_table[4][1],
	spec_390nm=leaf_table[5][1],
	spec_400nm=leaf_table[6][1],
	spec_410nm=leaf_table[7][1],
	spec_420nm=leaf_table[8][1],
	spec_430nm=leaf_table[9][1],
	spec_440nm=leaf_table[10][1],
	spec_450nm=leaf_table[11][1],
	spec_460nm=leaf_table[12][1],
	spec_470nm=leaf_table[13][1],
	spec_480nm=leaf_table[14][1],
	spec_490nm=leaf_table[15][1],
	spec_500nm=leaf_table[16][1],
	spec_510nm=leaf_table[17][1],
	spec_520nm=leaf_table[18][1],
	spec_530nm=leaf_table[19][1],
	spec_540nm=leaf_table[20][1],
	spec_550nm=leaf_table[21][1],
	spec_560nm=leaf_table[22][1],
	spec_570nm=leaf_table[23][1],
	spec_580nm=leaf_table[24][1],
	spec_590nm=leaf_table[25][1],
	spec_600nm=leaf_table[26][1],
	spec_610nm=leaf_table[27][1],
	spec_620nm=leaf_table[28][1],
	spec_630nm=leaf_table[29][1],
	spec_640nm=leaf_table[30][1],
	spec_650nm=leaf_table[31][1],
	spec_660nm=leaf_table[32][1],
	spec_670nm=leaf_table[33][1],
	spec_680nm=leaf_table[34][1],
	spec_690nm=leaf_table[35][1],
	spec_700nm=leaf_table[36][1],
	spec_710nm=leaf_table[37][1],
	spec_720nm=leaf_table[38][1],
	spec_730nm=leaf_table[39][1],
	spec_740nm=leaf_table[40][1],
	spec_750nm=leaf_table[41][1], illuminant='e')
	
	print("Python spectral color conversion:")
	print(convert_color(leaf_spectral, XYZColor))
	print(convert_color(leaf_spectral, sRGBColor))
	print(convert_color(leaf_spectral, HSLColor))
	print(convert_color(leaf_spectral, LabColor))
	print(convert_color(leaf_spectral, LCHabColor))
	print(convert_color(leaf_spectral, LCHuvColor))
	
	# I found Python's canonical conversion matrix. Maybe.
	#leaf_xyz1 = np.dot(np.linalg.inv(IPTColor.conversion_matrices['xyz_to_lms']), leaf_lms)
	#leaf_xyz2 = XYZColor(*leaf_xyz1, illuminant='d65')
	#print(leaf_xyz2)
	#print(convert_color(leaf_xyz2, sRGBColor))
	#print(convert_color(leaf_xyz2, HSLColor))
	#print(convert_color(leaf_xyz2, LabColor))
	#print(convert_color(leaf_xyz2, LCHabColor))
	#print(convert_color(leaf_xyz2, LCHuvColor))
	
	#leaf_xyz3 = np.array([*convert_color(leaf_spectral, XYZColor).get_value_tuple()])
	#print(np.dot(IPTColor.conversion_matrices['xyz_to_lms'], leaf_xyz3))
	#print(np.dot(IPTColor.conversion_matrices['xyz_to_lms'], leaf_matrix))
	
	print("")
	
	red_leaf_table = np.array([
		[340, 0.0],
		[350, 0.0],
		[360, 0.0],
		[370, 0.0],
		[380, 0.0],
		[390, 0.0],
		[400, 0.0],
		[410, 0.0],
		[420, 0.0],
		[430, 0.0],
		[440, 0.0],
		[450, 0.0],
		[460, 0.0],
		[470, 0.0],
		[480, 0.0],
		[490, 0.0],
		[500, 0.0],
		[510, 0.0],
		[520, 0.0],
		[530, 0.0],
		[540, 0.0],
		[550, 0.01],
		[560, 0.01],
		[570, 0.02],
		[580, 0.02],
		[590, 0.03],
		[600, 0.05],
		[610, 0.06],
		[620, 0.08],
		[630, 0.1],
		[640, 0.1],
		[650, 0.09],
		[660, 0.08],
		[670, 0.05],
		[680, 0.1],
		[690, 0.3],
		[700, 0.35],
		[710, 0.45],
		[720, 0.55],
		[730, 0.6],
		[740, 0.63],
		[750, 0.65]
	])
	red_leaf_l = 0
	red_leaf_m = 0
	red_leaf_s = 0
	for i in range(0, red_leaf_table.shape[0]):
		f = lens_filter(red_leaf_table[i][0], args.cutoff)
		red_leaf_l += wl_to_s(red_leaf_table[i][0], l1) * red_leaf_table[i][1] * f
		red_leaf_m += wl_to_s(red_leaf_table[i][0], m1) * red_leaf_table[i][1] * f
		red_leaf_s += wl_to_s(red_leaf_table[i][0], s1) * red_leaf_table[i][1] * f
	
	total = red_leaf_l + red_leaf_m + red_leaf_s
	red_leaf_l = red_leaf_l / total
	red_leaf_m = red_leaf_m / total
	red_leaf_s = red_leaf_s / total
	
	# ratio comparison
	similar = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, l0) / wl_to_s(i, m0), 2)
		ratio1 = round(red_leaf_l / red_leaf_m, 2)
		ratio2 = round(wl_to_s(i + 1, l0) / wl_to_s(i + 1, m0), 2)
		if (ratio0 == ratio1):
			similar = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar = i + 0.5
	
	similar1 = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, m0) / wl_to_s(i, s0), 2)
		ratio1 = round(red_leaf_m / red_leaf_s, 2)
		ratio2 = round(wl_to_s(i + 1, m0) / wl_to_s(i + 1, s0), 2)
		if (ratio0 == ratio1):
			similar1 = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar1 = i + 0.5
	
	print("LMS response to red maple leaf: l=" + str(red_leaf_l) +", m=" + str(red_leaf_m) + ", s=" + str(red_leaf_s))
	red_leaf_lms = np.array([
		[red_leaf_l],
		[red_leaf_m],
		[red_leaf_s]
	])
	print("L:M ratio: " + str(red_leaf_l / red_leaf_m) + " (looks like " + str(similar) + " nm)")
	print("M:S ratio: " + str(red_leaf_m / red_leaf_s) + " (looks like " + str(similar1) + " nm)")
	
	red_leaf_matrix = np.matmul(lms_to_xyz, red_leaf_lms)
	#red_leaf_matrix1 = np.matmul(xyz_to_rgb, red_leaf_matrix)
	red_leaf_xyz = XYZColor(*red_leaf_matrix, illuminant="e")
	#red_leaf_rgb = XYZColor(red_leaf_matrix[0], red_leaf_matrix[1], red_leaf_matrix[2])
	print("Color coordinates of red maple leaf:")
	print("CIE XYZ: " + str(red_leaf_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(red_leaf_rgb.get_value_tuple()))
	print("sRGB: " + str(convert_color(red_leaf_xyz, sRGBColor).get_value_tuple()))
	#print("RGB (converted through RGYB): " + str(rgyb_to_rgb((red_leaf_l - red_leaf_m) / (red_leaf_l + red_leaf_m + red_leaf_s), (red_leaf_l + red_leaf_m - red_leaf_s) / (red_leaf_l + red_leaf_m + red_leaf_s))))
	print("HSL: " + str(convert_color(red_leaf_xyz, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(red_leaf_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(red_leaf_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(red_leaf_xyz, LCHuvColor).get_value_tuple()))
	
	# "SpectralColor" object
	red_leaf_spectral = SpectralColor(spec_340nm=red_leaf_table[0][1],
	spec_350nm=red_leaf_table[1][1],
	spec_360nm=red_leaf_table[2][1],
	spec_370nm=red_leaf_table[3][1],
	spec_380nm=red_leaf_table[4][1],
	spec_390nm=red_leaf_table[5][1],
	spec_400nm=red_leaf_table[6][1],
	spec_410nm=red_leaf_table[7][1],
	spec_420nm=red_leaf_table[8][1],
	spec_430nm=red_leaf_table[9][1],
	spec_440nm=red_leaf_table[10][1],
	spec_450nm=red_leaf_table[11][1],
	spec_460nm=red_leaf_table[12][1],
	spec_470nm=red_leaf_table[13][1],
	spec_480nm=red_leaf_table[14][1],
	spec_490nm=red_leaf_table[15][1],
	spec_500nm=red_leaf_table[16][1],
	spec_510nm=red_leaf_table[17][1],
	spec_520nm=red_leaf_table[18][1],
	spec_530nm=red_leaf_table[19][1],
	spec_540nm=red_leaf_table[20][1],
	spec_550nm=red_leaf_table[21][1],
	spec_560nm=red_leaf_table[22][1],
	spec_570nm=red_leaf_table[23][1],
	spec_580nm=red_leaf_table[24][1],
	spec_590nm=red_leaf_table[25][1],
	spec_600nm=red_leaf_table[26][1],
	spec_610nm=red_leaf_table[27][1],
	spec_620nm=red_leaf_table[28][1],
	spec_630nm=red_leaf_table[29][1],
	spec_640nm=red_leaf_table[30][1],
	spec_650nm=red_leaf_table[31][1],
	spec_660nm=red_leaf_table[32][1],
	spec_670nm=red_leaf_table[33][1],
	spec_680nm=red_leaf_table[34][1],
	spec_690nm=red_leaf_table[35][1],
	spec_700nm=red_leaf_table[36][1],
	spec_710nm=red_leaf_table[37][1],
	spec_720nm=red_leaf_table[38][1],
	spec_730nm=red_leaf_table[39][1],
	spec_740nm=red_leaf_table[40][1],
	spec_750nm=red_leaf_table[41][1], illuminant='e')
	
	print("Python spectral color conversion:")
	print(convert_color(red_leaf_spectral, XYZColor))
	print(convert_color(red_leaf_spectral, sRGBColor))
	print(convert_color(red_leaf_spectral, HSLColor))
	print(convert_color(red_leaf_spectral, LabColor))
	print(convert_color(red_leaf_spectral, LCHabColor))
	print(convert_color(red_leaf_spectral, LCHuvColor))
	print("")
	
	# corn
	# estimated from https://www.yorku.ca/planters/photosynthesis/2014_08_15_lab_manual_static_html/images/Corn_leaf_reflectance.png. Values below 400 are estimated to be the same as 400.
	corn_table = np.array([
		#[300, 0.45],
		#[310, 0.45],
		#[320, 0.45],
		#[330, 0.45],
		[340, 0.45],
		[350, 0.45],
		[360, 0.45],
		[370, 0.45],
		[380, 0.45],
		[390, 0.45],
		[400, 0.45],
		[410, 0.4],
		[420, 0.35],
		[430, 0.35],
		[440, 0.35],
		[450, 0.35],
		[460, 0.37],
		[470, 0.37],
		[480, 0.37],
		[490, 0.37],
		[500, 0.43],
		[510, 0.5],
		[520, 0.6],
		[530, 0.77],
		[540, 0.8],
		[550, 0.8],
		[560, 0.77],
		[570, 0.75],
		[580, 0.65],
		[590, 0.63],
		[600, 0.62],
		[610, 0.6],
		[620, 0.55],
		[630, 0.55],
		[640, 0.55],
		[650, 0.45],
		[660, 0.35],
		[670, 0.35],
		[680, 0.4],
		[690, 0.6],
		[700, 0.65],
		[710, 0.75],
		[720, 0.8],
		[730, 0.83],
		[740, 0.83],
		[750, 0.83],
		[760, 0.83],
		[770, 0.8],
		[780, 0.8],
		[790, 0.8],
		[800, 0.8]
	])
	corn_l = 0
	corn_m = 0
	corn_s = 0
	for i in range(0, corn_table.shape[0]):
		f = lens_filter(corn_table[i][0], args.cutoff)
		corn_l += wl_to_s(corn_table[i][0], l1) * corn_table[i][1] * f
		corn_m += wl_to_s(corn_table[i][0], m1) * corn_table[i][1] * f
		corn_s += wl_to_s(corn_table[i][0], s1) * corn_table[i][1] * f
	
	total = corn_l + corn_m + corn_s
	corn_l = corn_l / total
	corn_m = corn_m / total
	corn_s = corn_s / total
	
	# LM ratio comparison
	similar = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, l0) / wl_to_s(i, m0), 2)
		ratio1 = round(corn_l / corn_m, 2)
		ratio2 = round(wl_to_s(i + 1, l0) / wl_to_s(i + 1, m0), 2)
		if (ratio0 == ratio1):
			similar = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar = i + 0.5
	
	similar1 = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, m0) / wl_to_s(i, s0), 2)
		ratio1 = round(corn_m / corn_s, 2)
		ratio2 = round(wl_to_s(i + 1, m0) / wl_to_s(i + 1, s0), 2)
		if (ratio0 == ratio1):
			similar1 = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar1 = i + 0.5
	
	print("LMS response to corn leaf: l=" + str(corn_l) +", m=" + str(corn_m) + ", s=" + str(corn_s))
	corn_lms = np.array([
		[corn_l],
		[corn_m],
		[corn_s]
	])
	print("L:M ratio: " + str(corn_l / corn_m) + " (looks like " + str(similar) + " nm)")
	print("M:S ratio: " + str(corn_m / corn_s) + " (looks like " + str(similar1) + " nm)")
	
	corn_matrix = np.matmul(lms_to_xyz, corn_lms)
	#corn_matrix1 = np.matmul(xyz_to_rgb, corn_matrix)
	corn_xyz = XYZColor(*corn_matrix, illuminant="e")
	#corn_rgb = XYZColor(corn_matrix1[0], corn_matrix1[1], corn_matrix1[2])
	print("Color coordinates of corn leaf:")
	print("CIE XYZ: " + str(corn_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(corn_rgb.get_value_tuple()))
	print("sRGB: " + str(convert_color(corn_xyz, sRGBColor).get_value_tuple()))
	#print("RGB (converted through RGYB): " + str(rgyb_to_rgb((corn_l - corn_m) / (corn_l + corn_m + corn_s), (corn_l + corn_m - corn_s) / (corn_l + corn_m + corn_s))))
	print("HSL: " + str(convert_color(corn_xyz, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(corn_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(corn_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(corn_xyz, LCHuvColor).get_value_tuple()))
	
	# "SpectralColor" object
	corn_spectral = SpectralColor(spec_340nm=corn_table[0][1],
	spec_350nm=corn_table[1][1],
	spec_360nm=corn_table[2][1],
	spec_370nm=corn_table[3][1],
	spec_380nm=corn_table[4][1],
	spec_390nm=corn_table[5][1],
	spec_400nm=corn_table[6][1],
	spec_410nm=corn_table[7][1],
	spec_420nm=corn_table[8][1],
	spec_430nm=corn_table[9][1],
	spec_440nm=corn_table[10][1],
	spec_450nm=corn_table[11][1],
	spec_460nm=corn_table[12][1],
	spec_470nm=corn_table[13][1],
	spec_480nm=corn_table[14][1],
	spec_490nm=corn_table[15][1],
	spec_500nm=corn_table[16][1],
	spec_510nm=corn_table[17][1],
	spec_520nm=corn_table[18][1],
	spec_530nm=corn_table[19][1],
	spec_540nm=corn_table[20][1],
	spec_550nm=corn_table[21][1],
	spec_560nm=corn_table[22][1],
	spec_570nm=corn_table[23][1],
	spec_580nm=corn_table[24][1],
	spec_590nm=corn_table[25][1],
	spec_600nm=corn_table[26][1],
	spec_610nm=corn_table[27][1],
	spec_620nm=corn_table[28][1],
	spec_630nm=corn_table[29][1],
	spec_640nm=corn_table[30][1],
	spec_650nm=corn_table[31][1],
	spec_660nm=corn_table[32][1],
	spec_670nm=corn_table[33][1],
	spec_680nm=corn_table[34][1],
	spec_690nm=corn_table[35][1],
	spec_700nm=corn_table[36][1],
	spec_710nm=corn_table[37][1],
	spec_720nm=corn_table[38][1],
	spec_730nm=corn_table[39][1],
	spec_740nm=corn_table[40][1],
	spec_750nm=corn_table[41][1],
	spec_760nm=corn_table[42][1],
	spec_770nm=corn_table[43][1],
	spec_780nm=corn_table[44][1],
	spec_790nm=corn_table[45][1],
	spec_800nm=corn_table[46][1], illuminant='e')
	
	spec_750nm=corn_table[41][1],
	print("Python spectral color conversion:")
	print(convert_color(corn_spectral, XYZColor))
	print(convert_color(corn_spectral, sRGBColor))
	print(convert_color(corn_spectral, HSLColor))
	print(convert_color(corn_spectral, LabColor))
	print(convert_color(corn_spectral, LCHabColor))
	print(convert_color(corn_spectral, LCHuvColor))
	print("")
	
	# hydrangea
	# estimated from https://spectralevolution.com/wp-content/uploads/2024/04/RT_hydrang_ref.jpg
	hydrangea_table = np.array([
		[340, 0.2],
		[350, 0.2],
		[360, 0.175],
		[370, 0.1375],
		[380, 0.125],
		[390, 0.1],
		[400, 0.0875],
		[410, 0.075],
		[420, 0.075],
		[430, 0.0625],
		[440, 0.0625],
		[450, 0.0625],
		[460, 0.0625],
		[470, 0.0625],
		[480, 0.0625],
		[490, 0.0625],
		[500, 0.0625],
		[510, 0.075],
		[520, 0.075],
		[530, 0.0875],
		[540, 0.1],
		[550, 0.1],
		[560, 0.0875],
		[570, 0.0875],
		[580, 0.075],
		[590, 0.075],
		[600, 0.075],
		[610, 0.0625],
		[620, 0.0625],
		[630, 0.0625],
		[640, 0.0625],
		[650, 0.0625],
		[660, 0.0625],
		[670, 0.0625],
		[680, 0.0625],
		[690, 0.0875],
		[700, 0.125],
		[710, 0.1875],
		[720, 0.3],
		[730, 0.375],
		[740, 0.45],
		[750, 0.4875],
		[760, 0.5125],
		[770, 0.5125],
		[780, 0.525],
		[790, 0.525],
		[800, 0.525],
		[810, 0.525],
		[820, 0.525],
		[830, 0.525]
	])
	hydrangea_l = 0
	hydrangea_m = 0
	hydrangea_s = 0
	for i in range(0, hydrangea_table.shape[0]):
		f = lens_filter(hydrangea_table[i][0], args.cutoff)
		hydrangea_l += wl_to_s(hydrangea_table[i][0], l1) * hydrangea_table[i][1] * f
		hydrangea_m += wl_to_s(hydrangea_table[i][0], m1) * hydrangea_table[i][1] * f
		hydrangea_s += wl_to_s(hydrangea_table[i][0], s1) * hydrangea_table[i][1] * f
	
	total = hydrangea_l + hydrangea_m + hydrangea_s
	hydrangea_l = hydrangea_l / total
	hydrangea_m = hydrangea_m / total
	hydrangea_s = hydrangea_s / total
	
	# LM ratio comparison
	similar = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, l0) / wl_to_s(i, m0), 2)
		ratio1 = round(hydrangea_l / hydrangea_m, 2)
		ratio2 = round(wl_to_s(i + 1, l0) / wl_to_s(i + 1, m0), 2)
		if (ratio0 == ratio1):
			similar = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar = i + 0.5
	
	similar1 = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, m0) / wl_to_s(i, s0), 2)
		ratio1 = round(hydrangea_m / hydrangea_s, 2)
		ratio2 = round(wl_to_s(i + 1, m0) / wl_to_s(i + 1, s0), 2)
		if (ratio0 == ratio1):
			similar1 = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar1 = i + 0.5
	
	print("LMS response to hydrangea leaf: l=" + str(hydrangea_l) +", m=" + str(hydrangea_m) + ", s=" + str(hydrangea_s))
	hydrangea_lms = np.array([
		[hydrangea_l],
		[hydrangea_m],
		[hydrangea_s]
	])
	print("L:M ratio: " + str(hydrangea_l / hydrangea_m) + " (looks like " + str(similar) + " nm)")
	print("M:S ratio: " + str(hydrangea_m / hydrangea_s) + " (looks like " + str(similar1) + " nm)")
	
	hydrangea_matrix = np.matmul(lms_to_xyz, hydrangea_lms)
	#hydrangea_matrix1 = np.matmul(xyz_to_rgb, hydrangea_matrix)
	hydrangea_xyz = XYZColor(*hydrangea_matrix, illuminant="e")
	#hydrangea_rgb = XYZColor(*hydrangea_matrix1)
	print("Color coordinates of hydrangea leaf:")
	print("CIE XYZ: " + str(hydrangea_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(hydrangea_rgb.get_value_tuple()))
	print("sRGB: " + str(convert_color(hydrangea_xyz, sRGBColor).get_value_tuple()))
	print("HSL: " + str(convert_color(hydrangea_xyz, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(hydrangea_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(hydrangea_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(hydrangea_xyz, LCHuvColor).get_value_tuple()))
	
	# "SpectralColor" object
	hydrangea_spectral = SpectralColor(spec_340nm=hydrangea_table[0][1],
	spec_350nm=hydrangea_table[1][1],
	spec_360nm=hydrangea_table[2][1],
	spec_370nm=hydrangea_table[3][1],
	spec_380nm=hydrangea_table[4][1],
	spec_390nm=hydrangea_table[5][1],
	spec_400nm=hydrangea_table[6][1],
	spec_410nm=hydrangea_table[7][1],
	spec_420nm=hydrangea_table[8][1],
	spec_430nm=hydrangea_table[9][1],
	spec_440nm=hydrangea_table[10][1],
	spec_450nm=hydrangea_table[11][1],
	spec_460nm=hydrangea_table[12][1],
	spec_470nm=hydrangea_table[13][1],
	spec_480nm=hydrangea_table[14][1],
	spec_490nm=hydrangea_table[15][1],
	spec_500nm=hydrangea_table[16][1],
	spec_510nm=hydrangea_table[17][1],
	spec_520nm=hydrangea_table[18][1],
	spec_530nm=hydrangea_table[19][1],
	spec_540nm=hydrangea_table[20][1],
	spec_550nm=hydrangea_table[21][1],
	spec_560nm=hydrangea_table[22][1],
	spec_570nm=hydrangea_table[23][1],
	spec_580nm=hydrangea_table[24][1],
	spec_590nm=hydrangea_table[25][1],
	spec_600nm=hydrangea_table[26][1],
	spec_610nm=hydrangea_table[27][1],
	spec_620nm=hydrangea_table[28][1],
	spec_630nm=hydrangea_table[29][1],
	spec_640nm=hydrangea_table[30][1],
	spec_650nm=hydrangea_table[31][1],
	spec_660nm=hydrangea_table[32][1],
	spec_670nm=hydrangea_table[33][1],
	spec_680nm=hydrangea_table[34][1],
	spec_690nm=hydrangea_table[35][1],
	spec_700nm=hydrangea_table[36][1],
	spec_710nm=hydrangea_table[37][1],
	spec_720nm=hydrangea_table[38][1],
	spec_730nm=hydrangea_table[39][1],
	spec_740nm=hydrangea_table[40][1],
	spec_750nm=hydrangea_table[41][1],
	spec_760nm=hydrangea_table[42][1],
	spec_770nm=hydrangea_table[43][1],
	spec_780nm=hydrangea_table[44][1],
	spec_790nm=hydrangea_table[45][1],
	spec_800nm=hydrangea_table[46][1],
	spec_810nm=hydrangea_table[47][1],
	spec_820nm=hydrangea_table[48][1],
	spec_830nm=hydrangea_table[49][1], illuminant='e')
	
	print("Python spectral color conversion:")
	print(convert_color(hydrangea_spectral, XYZColor))
	print(convert_color(hydrangea_spectral, sRGBColor))
	print(convert_color(hydrangea_spectral, HSLColor))
	print(convert_color(hydrangea_spectral, LabColor))
	print(convert_color(hydrangea_spectral, LCHabColor))
	print(convert_color(hydrangea_spectral, LCHuvColor))
	print("")
	
	# flowers
	# estimated from https://www.researchgate.net/profile/Robert-Gegear/publication/222035355/figure/fig1/AS:632345413578753@1527774303378/The-spectral-reflectance-curves-of-coloured-flowers-and-the-green-array-background-used.png
	flower0_table = np.array([
		#[300, 0.05],
		#[310, 0.08],
		#[320, 0.08],
		#[330, 0.08],
		[340, 0.08],
		[350, 0.08],
		[360, 0.08],
		[370, 0.1],
		[380, 0.13],
		[390, 0.2],
		[400, 0.35],
		[410, 0.55],
		[420, 0.65],
		[430, 0.73],
		[440, 0.83],
		[450, 0.85],
		[460, 0.88],
		[470, 0.9],
		[480, 0.9],
		[490, 0.88],
		[500, 0.85],
		[510, 0.8],
		[520, 0.75],
		[530, 0.65],
		[540, 0.55],
		[550, 0.4],
		[560, 0.25],
		[570, 0.15],
		[580, 0.1],
		[590, 0.08],
		[600, 0.08],
		[610, 0.05],
		[620, 0.05],
		[630, 0.05],
		[640, 0.05],
		[650, 0.05],
		[660, 0.08],
		[670, 0.08],
		[680, 0.08],
		[690, 0.08],
		[700, 0.08]
	])
	flower0_l = 0
	flower0_m = 0
	flower0_s = 0
	for i in range(0, flower0_table.shape[0]):
		f = lens_filter(flower0_table[i][0], args.cutoff)
		flower0_l += wl_to_s(flower0_table[i][0], l1) * flower0_table[i][1] * f
		flower0_m += wl_to_s(flower0_table[i][0], m1) * flower0_table[i][1] * f
		flower0_s += wl_to_s(flower0_table[i][0], s1) * flower0_table[i][1] * f
	
	total = flower0_l + flower0_m + flower0_s
	flower0_l = flower0_l / total
	flower0_m = flower0_m / total
	flower0_s = flower0_s / total
	
	# LM ratio comparison
	similar = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, l0) / wl_to_s(i, m0), 2)
		ratio1 = round(flower0_l / flower0_m, 2)
		ratio2 = round(wl_to_s(i + 1, l0) / wl_to_s(i + 1, m0), 2)
		if (ratio0 == ratio1):
			similar = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar = i + 0.5
	
	similar1 = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, m0) / wl_to_s(i, s0), 2)
		ratio1 = round(flower0_m / flower0_s, 2)
		ratio2 = round(wl_to_s(i + 1, m0) / wl_to_s(i + 1, s0), 2)
		if (ratio0 == ratio1):
			similar1 = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar1 = i + 0.5
	
	print("LMS response to blue flower: l=" + str(flower0_l) +", m=" + str(flower0_m) + ", s=" + str(flower0_s))
	flower0_lms = np.array([
		[flower0_l],
		[flower0_m],
		[flower0_s]
	])
	print("L:M ratio: " + str(flower0_l / flower0_m) + " (looks like " + str(similar) + " nm)")
	print("M:S ratio: " + str(flower0_m / flower0_s) + " (looks like " + str(similar1) + " nm)")
	
	flower0_matrix = np.matmul(lms_to_xyz, flower0_lms)
	#flower0_matrix1 = np.matmul(xyz_to_rgb, flower0_matrix)
	flower0_xyz = XYZColor(*flower0_matrix, illuminant="e")
	#flower0_rgb = XYZColor(*flower0_matrix1)
	print("Color coordinates of blue flower:")
	print("CIE XYZ: " + str(flower0_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(flower0_rgb.get_value_tuple()))
	print("sRGB: " + str(convert_color(flower0_xyz, sRGBColor).get_value_tuple()))
	print("HSL: " + str(convert_color(flower0_xyz, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(flower0_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(flower0_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(flower0_xyz, LCHuvColor).get_value_tuple()))
	flower0_spectral = SpectralColor(spec_340nm=flower0_table[0][1],
	spec_350nm=flower0_table[1][1],
	spec_360nm=flower0_table[2][1],
	spec_370nm=flower0_table[3][1],
	spec_380nm=flower0_table[4][1],
	spec_390nm=flower0_table[5][1],
	spec_400nm=flower0_table[6][1],
	spec_410nm=flower0_table[7][1],
	spec_420nm=flower0_table[8][1],
	spec_430nm=flower0_table[9][1],
	spec_440nm=flower0_table[10][1],
	spec_450nm=flower0_table[11][1],
	spec_460nm=flower0_table[12][1],
	spec_470nm=flower0_table[13][1],
	spec_480nm=flower0_table[14][1],
	spec_490nm=flower0_table[15][1],
	spec_500nm=flower0_table[16][1],
	spec_510nm=flower0_table[17][1],
	spec_520nm=flower0_table[18][1],
	spec_530nm=flower0_table[19][1],
	spec_540nm=flower0_table[20][1],
	spec_550nm=flower0_table[21][1],
	spec_560nm=flower0_table[22][1],
	spec_570nm=flower0_table[23][1],
	spec_580nm=flower0_table[24][1],
	spec_590nm=flower0_table[25][1],
	spec_600nm=flower0_table[26][1],
	spec_610nm=flower0_table[27][1],
	spec_620nm=flower0_table[28][1],
	spec_630nm=flower0_table[29][1],
	spec_640nm=flower0_table[30][1],
	spec_650nm=flower0_table[31][1],
	spec_660nm=flower0_table[32][1],
	spec_670nm=flower0_table[33][1],
	spec_680nm=flower0_table[34][1],
	spec_690nm=flower0_table[35][1],
	spec_700nm=flower0_table[36][1], illuminant='e')
	
	print("Python spectral color conversion:")
	print(convert_color(flower0_spectral, XYZColor))
	print(convert_color(flower0_spectral, sRGBColor))
	print(convert_color(flower0_spectral, HSLColor))
	print(convert_color(flower0_spectral, LabColor))
	print(convert_color(flower0_spectral, LCHabColor))
	print(convert_color(flower0_spectral, LCHuvColor))
	print("")
	
	flower1_table = np.array([
		#[310, 0.08],
		#[320, 0.08],
		#[330, 0.08],
		[340, 0.08],
		[350, 0.1],
		[360, 0.1],
		[370, 0.1],
		[380, 0.1],
		[390, 0.08],
		[400, 0.05],
		[410, 0.05],
		[420, 0.05],
		[430, 0.05],
		[440, 0.05],
		[450, 0.05],
		[460, 0.05],
		[470, 0.05],
		[480, 0.05],
		[490, 0.08],
		[500, 0.1],
		[510, 0.2],
		[520, 0.55],
		[530, 0.8],
		[540, 0.95],
		[550, 0.95],
		[560, 0.98],
		[570, 0.98],
		[580, 0.98],
		[590, 0.98],
		[600, 0.98],
		[610, 0.1],
		[620, 0.1],
		[630, 0.1],
		[640, 0.1],
		[650, 0.1],
		[660, 0.1],
		[670, 0.1],
		[680, 0.1],
		[690, 0.1],
		[700, 0.1]
	])
	flower1_l = 0
	flower1_m = 0
	flower1_s = 0
	for i in range(0, flower1_table.shape[0]):
		f = lens_filter(flower1_table[i][0], args.cutoff)
		flower1_l += wl_to_s(flower1_table[i][0], l1) * flower1_table[i][1] * f
		flower1_m += wl_to_s(flower1_table[i][0], m1) * flower1_table[i][1] * f
		flower1_s += wl_to_s(flower1_table[i][0], s1) * flower1_table[i][1] * f
	
	total = flower1_l + flower1_m + flower1_s
	flower1_l = flower1_l / total
	flower1_m = flower1_m / total
	flower1_s = flower1_s / total
	
	# LM ratio comparison
	similar = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, l0) / wl_to_s(i, m0), 2)
		ratio1 = round(flower1_l / flower1_m, 2)
		ratio2 = round(wl_to_s(i + 1, l0) / wl_to_s(i + 1, m0), 2)
		if (ratio0 == ratio1):
			similar = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar = i + 0.5
	
	similar1 = 0
	for i in range(300, 800):
		ratio0 = round(wl_to_s(i, m0) / wl_to_s(i, s0), 2)
		ratio1 = round(flower1_m / flower1_s, 2)
		ratio2 = round(wl_to_s(i + 1, m0) / wl_to_s(i + 1, s0), 2)
		if (ratio0 == ratio1):
			similar1 = i
		elif (ratio0 <= ratio1 <= ratio2
			or ratio0 >= ratio1 >= ratio2):
			similar1 = i + 0.5
	
	print("LMS response to yellow flower: l=" + str(flower1_l) +", m=" + str(flower1_m) + ", s=" + str(flower1_s))
	flower1_lms = np.array([
		[flower1_l],
		[flower1_m],
		[flower1_s]
	])
	print("L:M ratio: " + str(flower1_l / flower1_m) + " (looks like " + str(similar) + " nm)")
	print("M:S ratio: " + str(flower1_m / flower1_s) + " (looks like " + str(similar1) + " nm)")
	
	flower1_matrix = np.matmul(lms_to_xyz, flower1_lms)
	#flower1_matrix1 = np.matmul(xyz_to_rgb, flower1_matrix)
	flower1_xyz = XYZColor(*flower1_matrix, illuminant="e")
	#flower1_rgb = XYZColor(*flower1_matrix1)
	print("Color coordinates of yellow flower:")
	print("CIE XYZ: " + str(flower1_xyz.get_value_tuple()))
	#print("CIE RGB: " + str(flower1_rgb.get_value_tuple()))
	print("sRGB: " + str(convert_color(flower1_xyz, sRGBColor).get_value_tuple()))
	print("HSL: " + str(convert_color(flower1_xyz, HSLColor).get_value_tuple()))
	#print("CIE LAB: " + str(convert_color(flower1_xyz, LabColor).get_value_tuple()))
	print("LCh (LCHab): " + str(convert_color(flower1_xyz, LCHabColor).get_value_tuple()))
	#print("LCh (LCHuv): " + str(convert_color(flower1_xyz, LCHuvColor).get_value_tuple()))
	flower1_spectral = SpectralColor(spec_340nm=flower1_table[0][1],
	spec_350nm=flower1_table[1][1],
	spec_360nm=flower1_table[2][1],
	spec_370nm=flower1_table[3][1],
	spec_380nm=flower1_table[4][1],
	spec_390nm=flower1_table[5][1],
	spec_400nm=flower1_table[6][1],
	spec_410nm=flower1_table[7][1],
	spec_420nm=flower1_table[8][1],
	spec_430nm=flower1_table[9][1],
	spec_440nm=flower1_table[10][1],
	spec_450nm=flower1_table[11][1],
	spec_460nm=flower1_table[12][1],
	spec_470nm=flower1_table[13][1],
	spec_480nm=flower1_table[14][1],
	spec_490nm=flower1_table[15][1],
	spec_500nm=flower1_table[16][1],
	spec_510nm=flower1_table[17][1],
	spec_520nm=flower1_table[18][1],
	spec_530nm=flower1_table[19][1],
	spec_540nm=flower1_table[20][1],
	spec_550nm=flower1_table[21][1],
	spec_560nm=flower1_table[22][1],
	spec_570nm=flower1_table[23][1],
	spec_580nm=flower1_table[24][1],
	spec_590nm=flower1_table[25][1],
	spec_600nm=flower1_table[26][1],
	spec_610nm=flower1_table[27][1],
	spec_620nm=flower1_table[28][1],
	spec_630nm=flower1_table[29][1],
	spec_640nm=flower1_table[30][1],
	spec_650nm=flower1_table[31][1],
	spec_660nm=flower1_table[32][1],
	spec_670nm=flower1_table[33][1],
	spec_680nm=flower1_table[34][1],
	spec_690nm=flower1_table[35][1],
	spec_700nm=flower1_table[36][1], illuminant='e')
	
	print("Python spectral color conversion:")
	print(convert_color(flower1_spectral, XYZColor))
	print(convert_color(flower1_spectral, sRGBColor))
	print(convert_color(flower1_spectral, HSLColor))
	print(convert_color(flower1_spectral, LabColor))
	print(convert_color(flower1_spectral, LCHabColor))
	print(convert_color(flower1_spectral, LCHuvColor))
	print("")

# print execution time
print("%s seconds" % (time.time() - start_time))
