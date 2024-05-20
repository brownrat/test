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
# Ortín-Martínez A, Nadal-Nicolás FM, Jiménez-López M, Alburquerque-Béjar JJ, Nieto-López L, García-Ayuso D, Villegas-Pérez MP, Vidal-Sanz M, Agudo-Barriuso M. Number and distribution of mouse retinal cone photoreceptors: differences between an albino (Swiss) and a pigmented (C57/BL6) strain. PLoS One. 2014 Jul 16;9(7):e102392. doi: 10.1371/journal.pone.0102392. PMID: 25029531; PMCID: PMC4100816.
# Jacobs, G. H., Fenwick, J. A. and Williams, G. A. Cone-based vision of rats for ultraviolet and visible lights. Journal of Experimental Biology, 204(14), 15 July 2001. https://doi.org/10.1242/jeb.204.14.2439
# CIE 2-deg CMFs. http://cvrl.ucl.ac.uk/database/text/cienewxyz/cie2012xyz2.htm
#
# Further reading: https://xkcd.com/1926/

import cv2
import sys
import colorsys
import math
import numpy as np
from numpy.polynomial.polynomial import Polynomial
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
parser.add_argument("-l", "--lw", type=int, default=0, help="longwave cone sensitivity")
parser.add_argument("-m", "--mw", type=int, default=0, help="mediumwave cone sensitivity")
parser.add_argument("-s", "--sw", type=int, default=0, help="shortwave cone sensitivity")
parser.add_argument("--slw", type=int, default=565, help="source longwave cone sensitivity")
parser.add_argument("--smw", type=int, default=540, help="source mediumwave cone sensitivity")
parser.add_argument("--ssw", type=int, default=440, help="source shortwave cone sensitivity")
parser.add_argument("--cutoff", type=int, default=300, help="lower limit of lens transmission")
parser.add_argument("-r", "--red", type=int, default=650, help="reference red")
parser.add_argument("-y", "--yellow", type=int, default=570, help="reference yellow")
parser.add_argument("-g", "--green", type=int, default=540, help="reference green")
parser.add_argument("-c", "--cyan", type=int, default=490, help="reference cyan")
parser.add_argument("-b", "--blue", type=int, default=470, help="reference blue")
parser.add_argument("--violet", type=int, default=420, help="reference violet")
parser.add_argument("-i", "--image", help="image name")
parser.add_argument("-v", "--version", type=int, default=2, help="color conversion method")
parser.add_argument("--peak", help="show wavelength with maximum sensitivity", action="store_true")
parser.add_argument("--primaries", help="show primary and secondary wavelengths", action="store_true")
parser.add_argument("--white", default="e", help="reference white")
parser.add_argument("--lighting", help="show standard light sources", action="store_true")
parser.add_argument("--leaves", help="show several types of leaves", action="store_true")
parser.add_argument("--flowers", help="show several types of flowers", action="store_true")
parser.add_argument("--sky", help="show estimates for daylight and a clear sky", action="store_true")
parser.add_argument("--check", help="show colormath spectral rendering for comparison", action="store_true")
parser.add_argument("--lb", type=float, default=0.68990272, help="contribution of L cones to perception of brightness")
parser.add_argument("--mb", type=float, default=0.34832189, help="contribution of M cones to perception of brightness")
parser.add_argument("--sb", type=float, default=0.0, help="contribution of S cones to perception of brightness")
parser.add_argument("--triangle", help="print color triangle coordinates", action="store_true")
args = parser.parse_args()

# cone peak sensitivities
l1 = args.lw
m1 = args.mw
s1 = args.sw
if (l1 == 0 and m1 != 0):
	l1 = m1
if (m1 == 0 and l1 != 0):
	m1 = l1
if (s1 == 0 and m1 != 0):
	s1 = m1

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

# LMS to XYZ
lms_to_xyz = np.array([
	[1.94735469, -1.41445123, 0.36476327],
	[0.68990272, 0.34832189, 0],
	[0, 0, 1.93485343]
])

# relative wavelength sensitivity
# By default this does the same thing as the second row of the LMS->XYZ matrix, so the value
# it returns is the same as XYZ Y and so should be similar to the CIE luminous efficiency
# function. (The maximum sensitivity given for l=565 and m=540 is 555, which is the same.)
# This is equivalent to assuming brightness perception is mediated entirely by the L and
# (if present) M cones and the contribution from L cones is about twice that of M cones.
# There used to be an extra "lens filtering" factor, but I think this is misleading outside
# of removing wavelengths beyond a specified point. Note the relationship between XYZ Y,
# HSL L and Lab/LCH L is non-linear.
# The percentages of L, M and S cones are about 51-76%, 20-44% and 2% in humans (Wikipedia)
# and 68-73%, 20-25% and 7% in the fat-tailed dunnart ("Diversity of Color Vision [sic]: Not All
# Australian Marsupials are Trichromatic"). In Norway rats, the percentage of S cones is 11-12%
# ("Cone-based vision of rats for ultraviolet and visible lights") and in the house mouse the
# percentage of exclusive non-coexpressing S cones is 26% ("Number and Distribution of Mouse
# Retinal Cone Photoreceptors").
def lens_filter(w, c=300):
	# remove filter when cutoff is set to 0
	if (c == 0):
		return 1
	elif (w < c):
		return 0
	#else:
		#value = np.log(w - c) / np.log(800 - c)
		#if (value > 0): # a log function has the wrong shape near the cutoff point
			#return value
	#return 0
	return 1

def sensitivity(w, l=l1, m=m1, s=s1, c=300):
	value = (args.lb*wl_to_s(w, l) + args.mb*wl_to_s(w, m) + args.sb*wl_to_s(w, m)) * lens_filter(w, c)
	return value

# black body
# Note this expects a wavelength in meters.
def blackbody(w, t):
	c1 = 3.74183e-16
	c2 = 1.4388e-2
	return c1*w**-5 / (math.exp(c2/(w*t)) - 1)
	
# light sources
d65 = spectral_constants.REF_ILLUM_TABLE["d65"]
e = spectral_constants.REF_ILLUM_TABLE["e"]
incandescent = np.empty(50)
for i in range(0, 50):
	w = i*10 + 340
	# normalize to 1% at 550 nm
	incandescent[i] = blackbody(w / 1000000000, 2000) / blackbody(550 / 1000000000, 2000)

# white point
if (args.white == "d65"):
	wp = d65
elif (args.white == "i"):
	wp = incandescent
else:
	wp = e

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
# This function now prints a warning instead of causing a crash when it produces an overflow
# exception, but this usually seems to be a sign that something went wrong, so that doesn't
# help much.
def wl_to_s(w, lmax):
	try:
		value = 1 / (math.exp(69.7*(0.8795 + 0.0459*math.exp(-(lmax - 300)**2 / 11940) - lmax/w)) + math.exp(28*(0.922 - lmax/w)) + math.exp(-14.9*(1.104 - lmax/w)) + 0.674) + 0.26*math.exp(-((w - (189 + 0.315*lmax))/(-40.5 + 0.195*lmax))**2)
	except OverflowError:
		print("Warning: math overflow, clipping to 2.2250738585072014e-308")
		return 2.2250738585072014e-308
	
	return value

def wl_to_s1(w, lmax):
	try:
		value = 1 / (math.exp(69.7*(0.8795 - lmax/w)) + math.exp(28*(0.922 - lmax/w)) + math.exp(-14.9*(1.104 - lmax/w)) + 0.674)
	except OverflowError:
		print("Warning: math overflow, clipping to 2.2250738585072014e-308")
		return 2.2250738585072014e-308
	
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

# slightly clean up blue/green hues by normalizing the LMS values to equal at the white point
if (version == 0):
	white_l = 0
	white_m = 0
	white_s = 0
	for i in range(wp.shape[0]):
		w = i*10 + 340
		if (w >= args.cutoff):
			white_l += wl_to_s(w, l1) * wp[i]
			white_m += wl_to_s(w, m1) * wp[i]
			white_s += wl_to_s(w, s1) * wp[i]

def wl_to_rgb_0(w, l=l1, m=m1, s=s1):
	red = wl_to_s(w, l)
	green = wl_to_s(w, m)
	blue = wl_to_s(w, s)
	
	red = red * white_s / white_l
	green = green * white_s / white_m
	
	total = red+green+blue
	return [red/total, green/total, blue/total]

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
def find_match(w, l=l1, m=m1, s=s1, d=0.001):
	lp = wl_to_s(w, l0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
	mp = wl_to_s(w, m0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
	sp = wl_to_s(w, s0) / (wl_to_s(w, l0) + wl_to_s(w, m0) + wl_to_s(w, s0))
	
	p0 = 700
	p1 = 300
	
	#print(d)
	#while (abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
	#	+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
	#	+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s))) > d
	#	and p0 > 300):
	#	p0 -= 1
	#while (abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
	#	+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
	#	+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s))) > d
	#	and p1 < 700):
	#	p1 += 1
	# lower precision if not found
	#if (p0 == 300 and p1 == 700):
	#	# recursion
	#	p2 = find_match(w, l, m, s, d*2) # avoid running this twice
	#	p0 = p2[0]
	#	p1 = p2[1]

	i = 300
	while (i < 800):
		diff0 = abs(lp - wl_to_s(p0, l) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
		+ abs(mp - wl_to_s(p0, m) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
		+ abs(sp - wl_to_s(p0, s) / (wl_to_s(p0, l) + wl_to_s(p0, m) + wl_to_s(p0, s)))
		diff1 = abs(lp - wl_to_s(i, l) / (wl_to_s(i, l) + wl_to_s(i, m) + wl_to_s(i, s)))
		+ abs(mp - wl_to_s(i, m) / (wl_to_s(i, l) + wl_to_s(i, m) + wl_to_s(i, s)))
		+ abs(sp - wl_to_s(i, s) / (wl_to_s(i, l) + wl_to_s(i, m) + wl_to_s(i, s)))
		if (diff1 < diff0):
			p0 = i
		i += 1
	
	j = 800
	while (j > 300):
		diff0 = abs(lp - wl_to_s(p1, l) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
		+ abs(mp - wl_to_s(p1, m) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
		+ abs(sp - wl_to_s(p1, s) / (wl_to_s(p1, l) + wl_to_s(p1, m) + wl_to_s(p1, s)))
		diff1 = abs(lp - wl_to_s(j, l) / (wl_to_s(j, l) + wl_to_s(j, m) + wl_to_s(j, s)))
		+ abs(mp - wl_to_s(j, m) / (wl_to_s(j, l) + wl_to_s(j, m) + wl_to_s(j, s)))
		+ abs(sp - wl_to_s(j, s) / (wl_to_s(j, l) + wl_to_s(j, m) + wl_to_s(j, s)))
		if (diff1 < diff0):
			p1 = j
		j -= 1
	
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
# This converts LMS values to CIE XYZ and then to HLS.
# The results are very similar to CVS and to version 4 aside from a slight green shift for
# narrow spacing, probably for the same reason as in version 2. The main problem is the
# range of hues provided for integer wavelengths is full of holes you could drive a truck
# through. These gaps can be filled in convincingly with some careful calculations (probably
# the reason for the long execution time) for primate-like values, but for other arrangements
# it doesn't work as well.
# XYZ conversion can also simulate dichromacy in the same way as v0 and v4. Tritanopia is a
# bit odd-looking as red becomes purple for some reason.
# I've fixed some of the problems, but it's still not appropriate for non-primate vision.
# Almost everything ends up as pure red/green. This is probably because (a) the "threshold"
# for pure hues is based on a system with strong overlap, (b) a hue may be represented by
# a range of wavelengths in the source but one or none in the target and the inevitable
# averaging leads to posterization.

def cs_tables_3(l, m, s):
	cs_source = np.empty([501, 4])
	cs_target = np.empty([501, 4])
	
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
		#print(i+300)
		#print(hsl)
		
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
		#print(i+300)
		#print(hsl1)
	
	return [cs_source, cs_target]

# As with version 1, this one does some intense gap-filling to ensure the result looks sensible.
def hue_to_wavelength_3(h, source):
	if (0 <= h < 260): # XYZ does not produce any hues greater than this, so we cut off here
		w1 = 300
		w2 = 800
		i = 300
		while (i < 800):
			diff0 = abs(h - source[w1 - 300][3])
			diff1 = abs(h - source[i - 300][3])
			if (diff1 < diff0):
				w1 = i
			i += 1
		
		j = 800
		while (j > 300):
			diff0 = abs(h - source[w2 - 300][3])
			diff1 = abs(h - source[j - 300][3])
			if (diff1 < diff0):
				w2 = j
			j -= 1
		
		# red is a cutoff point and blue has some weird behavior
		if (h == 0):
			return w1
		elif (200 <= h < 260):
			return w2
		else:
			return (w1 + w2) / 2
	else:
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)

def wl_to_rgb_5(w, target):
	if (w == round(w)):
		return colorsys.hls_to_rgb(target[round(w) - 300][3] / 360, 0.5, 1)
	else:
		return colorsys.hls_to_rgb((target[int(w) - 300][3] + target[int(w) - 299][3]) / 720, 0.5, 1)

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

# version 6: dichromacy
# 0, 4 and 5 can produce something like a dichromatic simulation but aren't really built
# for it. This uses a similar method to v2. The results are more realistic, but with S
# values < 430 the neutral point is placed at an unexpectedly short wavelength. This has
# been improved by replacing E with D65.
def find_primaries_dc(l=l1, m=m1, s=s1, white=d65):

	# yellow/blue
	if (l == m):
		# find white/neutral point based on D65
		white_l = 0
		white_s = 0
		for i in range(white.shape[0]):
			w = i*10 + 340
			if (w >= args.cutoff):
				white_l += wl_to_s(w, l) * white[i]
				white_s += wl_to_s(w, s) * white[i]
		
		#white_ratio = white_l / white_s
		white_ls = white_l / (white_l + white_s)
		white_ls = (white_l - white_s) / (white_l + white_s)
		
		n = 300
		for i in range(300, 800):
			diff0 = abs(white_ls - (wl_to_s(n, l) - wl_to_s(n, s)) / (wl_to_s(n, l) + wl_to_s(n, s)))
			diff1 = abs(white_ls - (wl_to_s(i, l) - wl_to_s(i, s)) / (wl_to_s(i, l) + wl_to_s(i, s)))
			if (diff1 < diff0):
				n = i
		
		# find "yellow" and "blue" points based on edges of color space
		
		y = 300
		for i in range(300, 800):
			diff0 = abs(1 - round((wl_to_s(y, l) - wl_to_s(y, s)) / (wl_to_s(y, l) + wl_to_s(y, s)), 2))
			diff1 = abs(1 - round((wl_to_s(i, l) - wl_to_s(i, s)) / (wl_to_s(i, l) + wl_to_s(i, s)), 2))
			if (diff1 < diff0):
				y = i
		
		b = 800
		i = 800
		while (i > 300):
			diff0 = abs(-1 - round((wl_to_s(b, l) - wl_to_s(b, s)) / (wl_to_s(b, l) + wl_to_s(b, s)), 2))
			diff1 = abs(-1 - round((wl_to_s(i, l) - wl_to_s(i, s)) / (wl_to_s(i, l) + wl_to_s(i, s)), 2))
			if (diff1 < diff0):
				b = i
			i -= 1
		
		return [n, y, b]
	
	# red/cyan (tritanopia)
	elif (m == s):
		# find white/neutral point
		white_l = 0
		white_m = 0
		for i in range(300, 800):
			white_l += wl_to_s(i, l) * lens_filter(i, args.cutoff)
			white_m += wl_to_s(i, m) * lens_filter(i, args.cutoff)
		
		white_ratio = white_l / white_m
		
		n = m
		for i in range(m, 800):
			diff0 = abs(white_ratio - wl_to_s(n, l) / wl_to_s(n, m))
			diff1 = abs(white_ratio - wl_to_s(i, l) / wl_to_s(i, m))
			if (diff1 < diff0):
				n = i
		
		# find "red" and "cyan" points based on reference red and green
		r_ratio = wl_to_s(r0, l0) / (wl_to_s(r0, m0) + wl_to_s(r0, s0)) 
		c_ratio = wl_to_s(g0, l0) / (wl_to_s(g0, m0) + wl_to_s(g0, s0))
		
		r = n
		for i in range(n, 800):
			diff0 = abs(r_ratio - wl_to_s(r, l) / wl_to_s(r, s))
			diff1 = abs(r_ratio - wl_to_s(i, l) / wl_to_s(i, s))
			if (diff1 < diff0):
				r = i
		
		c = m
		for i in range(m, n):
			diff0 = abs(c_ratio - wl_to_s(c, l) / wl_to_s(c, s))
			diff1 = abs(c_ratio - wl_to_s(i, l) / wl_to_s(i, s))
			if (diff1 < diff0):
				c = i
		
		return [n, r, c]
	
	else:
		print("Please specify a valid dichromatic phenotype.")
		return False

def hue_saturation_dc(w, p):
	if (l1 == m1):
		n = p[0]
		y = p[1]
		b = p[2]
		
		if (w > n): # yellow
			if (w <= y):
				return [60, (w - n) / (y - n)] # desaturate
			else:
				return [60, 1]
		elif (w == n): # neutral
			return [0, 0]
		else: # blue
			if (w >= b):
				return [240, (w - n) / (b - n)]
			else:
				return [240, 1]
		
	elif (m1 == s1):
		n = p[0]
		r = p[1]
		c = p[2]
	
		if (w > n): # red
			if (w <= r):
				return [0, (w - n) / (r - n)] # desaturate
			else:
				return [0, 1]
		elif (w == n): # neutral
			return [0, 0]
		else: # cyan
			if (w >= c):
				return [180, (w - n) / (c - n)]
			else:
				return [180, 1]
	
	else:
		print("Please specify a valid dichromatic phenotype.")
		return False

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
elif (version == 6):
	primaries = find_primaries_dc(l1, m1, s1)

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

# estimate hue, saturation and lightness for a spectral power distribution
# The "brightness" value is basically the same thing as XYZ Y and the CIE luminous efficiency
# function (which are basically the same thing), so it should be a reasonable estimate for
# a species that derives brightness from cone signals in a similar way viewing an object
# under daylight (photopic) conditions. Otherwise it might not be reasonable.
def spectral_rendering(table, light_source=wp/100):
	table_l = 0
	table_m = 0
	table_s = 0
	brightness = 0
	for i in range(0, table.shape[0]):
		w = i*10 + 340 # wavelength
		if (w >= args.cutoff): # we do a "hard" cutoff because sensitivity is
			# probably adjusted to compensate for filtering
			table_l += wl_to_s(w, l1) * table[i] * light_source[i]
			# remove either M or S for dichromacy
			if (m1 != l1):
				table_m += wl_to_s(w, m1) * table[i] * light_source[i]
			if (s1 != m1):
				table_s += wl_to_s(w, s1) * table[i] * light_source[i]
			
			# brightness
			brightness += sensitivity(w, c=args.cutoff) * table[i] * light_source[i]
	
	# normalize according to provided light source: the total sensitivity of the
	# target phenotype to this light spectrum is "1"
	n = 0
	for i in range(0, light_source.shape[0]):
		w = i*10 + 340
		n += sensitivity(w, c=args.cutoff) * light_source[i]
	
	table_l = table_l / n
	table_m = table_m / n
	table_s = table_s / n
	
	# express brightness of reflected/emitted light proportional to the light source
	brightness = brightness / n
	
	# estimate color
	# For trichromacy we convert the LMS values to other color spaces to roughly visualize
	# the hue and saturation. This assumes a given ratio of cone absorption produces the
	# same perception as it does in humans and can be adequately modeled by an LMS<->XYZ
	# transformation matrix I found on Wikipedia. For dichromacy and monochromacy this
	# doesn't really work, so we just report the estimated brightness and (for dichromacy)
	# the L:S ratio.
	if (l1 != m1 != s1):
		print("Cone response: l=" + str(table_l) + ", m=" + str(table_m) + ", s=" + str(table_s))
		print("L:M ratio: " + str(table_l / table_m))
		print("M:S ratio: " + str(table_m / table_s))
		
		lms = np.array([
			[table_l],
			[table_m],
			[table_s]
		])
		
		matrix = np.matmul(lms_to_xyz, lms)
		xyz = XYZColor(*matrix, illuminant="e")
		print("Color coordinates:")
		print("CIE XYZ: " + str(xyz.get_value_tuple()))
		print("sRGB: " + str(convert_color(xyz, sRGBColor).get_value_tuple()))
		print("CIE LAB: " + str(convert_color(xyz, LabColor).get_value_tuple()))
	elif (l1 == m1 != s1):
		# typical dichromacy
		print("Cone response: l/m=" + str(table_l) + ", s=" + str(table_s))
		print("L:S ratio: " + str(table_l / table_s))
		
		# find matching wavelength
		p = find_primaries_dc(l1, m1, s1, white=wp)
		neutral = p[0]
		yellow = p[1]
		blue = p[2]
		
		match = args.cutoff
		for i in range(args.cutoff, 800):
			diff0 = abs(table_l / table_s - wl_to_s(match, l1) / wl_to_s(match, s1))
			diff1 = abs(table_l / table_s - wl_to_s(i, l1) / wl_to_s(i, s1))
			if (diff1 < diff0):
				match = i
		print("Matching wavelength: " + str(match))
		
		n_ratio = wl_to_s(neutral, l1) / wl_to_s(neutral, s1)
		y_ratio = wl_to_s(yellow, l1) / wl_to_s(yellow, s1)
		b_ratio = wl_to_s(blue, l1) / wl_to_s(blue, s1)
		match_ratio = table_l / table_s
		if (match == neutral):
			print("Color type: neutral")
		elif (match_ratio >= y_ratio):
			print("Color type: L")
		elif (n_ratio < match_ratio < y_ratio):
			print("Color type: L>S")
		elif (b_ratio < match_ratio < n_ratio):
			print("Color type: S>L")
		else:
			print("Color type: S")
	elif (l1 != m1 == s1):
		# tritanopia
		print("Cone response: l=" + str(table_l) + ", m=" + str(table_m))
	
	# estimate brightness
	print("Relative brightness: " + str(round(100*brightness, 2)) + "%")
	
	# "SpectralColor" conversion for comparison to check if our estimates are reasonable
	if (args.check):
		spectral = SpectralColor(spec_340nm=table[0],
		spec_350nm=table[1],
		spec_360nm=table[2],
		spec_370nm=table[3],
		spec_380nm=table[4],
		spec_390nm=table[5],
		spec_400nm=table[6],
		spec_410nm=table[7],
		spec_420nm=table[8],
		spec_430nm=table[9],
		spec_440nm=table[10],
		spec_450nm=table[11],
		spec_460nm=table[12],
		spec_470nm=table[13],
		spec_480nm=table[14],
		spec_490nm=table[15],
		spec_500nm=table[16],
		spec_510nm=table[17],
		spec_520nm=table[18],
		spec_530nm=table[19],
		spec_540nm=table[20],
		spec_550nm=table[21],
		spec_560nm=table[22],
		spec_570nm=table[23],
		spec_580nm=table[24],
		spec_590nm=table[25],
		spec_600nm=table[26],
		spec_610nm=table[27],
		spec_620nm=table[28],
		spec_630nm=table[29],
		spec_640nm=table[30],
		spec_650nm=table[31],
		spec_660nm=table[32],
		spec_670nm=table[33],
		spec_680nm=table[34],
		spec_690nm=table[35],
		spec_700nm=table[36],
		spec_710nm=table[37],
		spec_720nm=table[38],
		spec_730nm=table[39],
		spec_740nm=table[40],
		spec_750nm=table[41],
		spec_760nm=table[42],
		spec_770nm=table[43],
		spec_780nm=table[44],
		spec_790nm=table[45],
		spec_800nm=table[46],
		spec_810nm=table[47],
		spec_820nm=table[48],
		spec_830nm=table[49], illuminant='e')
		
		print("colormath conversion:")
		print(convert_color(spectral, XYZColor))
		print(convert_color(spectral, sRGBColor))
		print(convert_color(spectral, LabColor))
	
	# break up text
	print("")

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
			
			# convert hue from 0-180 (OpenCV format) to 0-360
			hue = hls[0]*2
			
			# dichromacy
			if (version == 6):
				if (hue >= 270):
					rb = colorsys.hls_to_rgb(hue / 360, 0.5, 1)
					hue_sat_r = hue_saturation_dc(r0, primaries)
					hue_sat_b = hue_saturation_dc(b0, primaries)
					rgb_r = colorsys.hls_to_rgb(hue_sat_r[0]/360, 0.5, hue_sat_r[1])
					rgb_b = colorsys.hls_to_rgb(hue_sat_b[0]/360, 0.5, hue_sat_b[1])
					
					if (l1 == m1):
						# "regular" dichromacy (protanopia, deuteranopia and typical
						# placental mammal vision)
						rg = (rb[0]*rgb_r[0] + rb[2]*rgb_b[0]) / 2
						b = rb[0]*rgb_r[2] + rb[2]*rgb_b[2]
						total = rg+b
						rgb = [rg/total, rg/total, b/total]
					else:
						# tritanopia
						r = rb[0]*rgb_r[0] + rb[2]*rgb_b[0]
						gb = (rb[0]*rgb_r[2] + rb[2]*rgb_b[2]) / 2
						rgb = [r, gb, gb]
					hue_sat0 = colorsys.rgb_to_hls(*rgb)
					hue_sat = [hue_sat0[0]*360, hue_sat0[1]]
				else:
					w = hue_to_wavelength(hue)
					hue_sat = hue_saturation_dc(w, primaries)
				
				# Hue is slightly shifted for a more natural appearance. The default
				# "dark yellow" pea-soup color is really unpleasant looking and has
				# more to do with how an sRGB display works than what real dichromacy
				# looks like. CVS doesn't do this, but some simulations do. An LCh
				# hue shift would look even better but may not be possible.
				if (l1 == m1):
					hue1 = hue_sat[0] - 10
				else:
					hue1 = hue_sat[0] + 10
				img_hls[x][y] = [hue1 / 2, hls[1], abs(hls[2] * hue_sat[1])]
			# trichromacy
			else:
				# convert hue to predominant wavelength(s)
				wl = hue_to_wavelength(hue)
				
				if (type(wl) != tuple):
					hue_target = wl_to_rgb(wl, l1, m1, s1)
				else:
					hue_r = wl_to_rgb(r0, l1, m1, s1)
					hue_b = wl_to_rgb(b0, l1, m1, s1)
					# sum the amounts of "red" and "blue"
					hue_target = [hue_r[0] * wl[0] + hue_b[0] * wl[2], hue_r[1] * wl[0] + hue_b[1] * wl[2], hue_r[2] * wl[0] + hue_b[2] * wl[2]]
				
				# convert predominant wavelengths back into a hue
				hls_target = colorsys.rgb_to_hls(*hue_target)
				
				if (version == 4):
					sat_diff = saturation(wl, *tables)
					#sat_diff = 1
				elif (version == 5):
					sat_diff = saturation_1(wl, *tables, l1, m1, s1)
				else:
					sat_diff = hls_target[2]
				
				# shift hue in our pixel. Colorsys uses 0-1, so we have to convert back to
				# OpenCV format.
				img_hls[x][y] = [hls_target[0]*180, hls[1], hls[2] * sat_diff]

	# convert back to BGR
	img_result = cv2.cvtColor(img_hls, cv2.COLOR_HLS2BGR)

	# fix brightness
	# I tried using LCh instead of HLS, but it doesn't work the way I expected. This
	# extra step doesn't add much time. However, "fixing" the lightness may be just as
	# misleading as keeping the HSL values.
	#img_lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
	#img_result_lab = cv2.cvtColor(img_result, cv2.COLOR_BGR2LAB)
	#for x in range(img.shape[0]):
	#	for y in range(img.shape[1]):
	#		img_result_lab[x][y][0] = img_lab[x][y][0]

	#img_result = cv2.cvtColor(img_result_lab, cv2.COLOR_LAB2BGR)

	# display result
	cv2.imwrite("colorvisionpy-result.png", img_result)

# if no image given, print general information
else:
	print("L: " + str(l1) + ", M: " + str(m1) + ", S: " + str(s1))
	print("")
	
	# maximum sensitivity
	if (args.peak):
		ms = 300
		for i in range(300, 800):
			if (sensitivity(i, c=args.cutoff) > (sensitivity(ms, c=args.cutoff))):
				ms = i
		print("Maximum sensitivity: " + str(ms))
		print("")
	
	# primary/secondary colors
	if (args.primaries):
		if (l1 == m1 or m1 == s1): # dichromacy
			primaries = find_primaries_dc(white=wp)
			print("Estimated primary and secondary wavelengths:")
			if (l1 == m1):
				print("L/M: " + str(primaries[1]))
				print("Neutral point: " + str(primaries[0]))
				print("S: " + str(primaries[2]))
			else:
				print("L: " + str(primaries[1]))
				print("Neutral point: " + str(primaries[0]))
				print("M: " + str(primaries[2]))
			print("")
		else:
			r = find_match(r0)
			y = find_match(y0)
			g = find_match(g0)
			c = find_match(c0)
			b = find_match(b0)
			v = find_match(v0)
			print("Estimated primary and secondary wavelengths based on cone response ratios:")
			print("red: " + str((r[0] + r[1]) / 2) + " (" + str(r[0]) + "-" + str(r[1]) + ")")
			print("yellow: " + str((y[0] + y[1]) / 2) + " (" + str(y[0]) + "-" + str(y[1]) + ")")
			print("green: " + str((g[0] + g[1]) / 2) + " (" + str(g[0]) + "-" + str(g[1]) + ")")
			print("cyan: " + str((c[0] + c[1]) / 2) + " (" + str(c[0]) + "-" + str(c[1]) + ")")
			print("blue: " + str((b[0] + b[1]) / 2) + " (" + str(b[0]) + "-" + str(b[1]) + ")")
			print("violet: " + str((v[0] + v[1]) / 2) + " (" + str(v[0]) + "-" + str(v[1]) + ")")
			print("")
			
			table = cs_tables_3(l1, m1, s1)[1]
			
			# red and violet are cutoff points
			r = 300
			v = 800
			i = 300
			while (i < 800):
				diff0 = abs(0 - table[r - 300][3])
				diff1 = abs(0 - table[i - 300][3])
				if (diff1 < diff0):
					r = i
				i += 1
			
			j = 800
			while (j > 300):
				diff0 = abs(270 - table[v - 300][3])
				diff1 = abs(270 - table[j - 300][3])
				if (diff1 < diff0):
					v = j
				j -= 1
			
			# yellow through blue may be a range of hues, so we try to find the middle
			y1 = 300
			y2 = 800
			i = 300
			while (i < 800):
				diff0 = abs(60 - table[y1 - 300][3])
				diff1 = abs(60 - table[i - 300][3])
				if (diff1 < diff0):
					y1 = i
				i += 1
			
			j = 800
			while (j > 300):
				diff0 = abs(60 - table[y2 - 300][3])
				diff1 = abs(60 - table[j - 300][3])
				if (diff1 < diff0):
					y2 = j
				j -= 1
			
			y = (y1 + y2) / 2
			
			g1 = 300
			g2 = 800
			i = 300
			while (i < 800):
				diff0 = abs(120 - table[g1 - 300][3])
				diff1 = abs(120 - table[i - 300][3])
				if (diff1 < diff0):
					g1 = i
				i += 1
			
			j = 800
			while (j > 300):
				diff0 = abs(120 - table[g2 - 300][3])
				diff1 = abs(120 - table[j - 300][3])
				if (diff1 < diff0):
					g2 = j
				j -= 1
			
			g = (g1 + g2) / 2
			
			c1 = 300
			c2 = 800
			i = 300
			while (i < 800):
				diff0 = abs(180 - table[c1 - 300][3])
				diff1 = abs(180 - table[i - 300][3])
				if (diff1 < diff0):
					c1 = i
				i += 1
			
			j = 800
			while (j > 300):
				diff0 = abs(180 - table[c2 - 300][3])
				diff1 = abs(180 - table[j - 300][3])
				if (diff1 < diff0):
					c2 = j
				j -= 1
			
			c = (c1 + c2) / 2
			
			b1 = 300
			b2 = 800
			i = 300
			while (i < 800):
				diff0 = abs(230 - table[b1 - 300][3])
				diff1 = abs(230 - table[i - 300][3])
				if (diff1 < diff0):
					b1 = i
				i += 1
			
			j = 800
			while (j > 300):
				diff0 = abs(220 - table[b2 - 300][3])
				diff1 = abs(220 - table[j - 300][3])
				if (diff1 < diff0):
					b2 = j
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
			
			# color triangle
			
			# primaries
			r = 300
			for i in range(300, 800):
				diff0 = abs(math.dist((1, -1), (round((wl_to_s(r, l1) - wl_to_s(r, s1)) / (wl_to_s(r, l1) + (wl_to_s(r, m1) + wl_to_s(r, s1))), 1), round((wl_to_s(r, m1) - wl_to_s(r, l1) - wl_to_s(r, s1)) / (wl_to_s(r, l1) + wl_to_s(r, m1) + wl_to_s(r, s1)), 1))))
				diff1 = abs(math.dist((1, -1), (round((wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + (wl_to_s(i, m1) + wl_to_s(i, s1))), 1), round((wl_to_s(i, m1) - wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + wl_to_s(i, m1) + wl_to_s(i, s1)), 1))))
				if (diff1 < diff0):
					r = i
			
			g = 300
			for i in range(300, 800):
				diff0 = abs(math.dist((0, 1), ((wl_to_s(g, l1) - wl_to_s(g, s1)) / (wl_to_s(g, l1) + (wl_to_s(g, m1) + wl_to_s(g, s1))), (wl_to_s(g, m1) - wl_to_s(g, l1) - wl_to_s(g, s1)) / (wl_to_s(g, l1) + wl_to_s(g, m1) + wl_to_s(g, s1)))))
				diff1 = abs(math.dist((0, 1), ((wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + (wl_to_s(i, m1) + wl_to_s(i, s1))), (wl_to_s(i, m1) - wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + wl_to_s(i, m1) + wl_to_s(i, s1)))))
				if (diff1 < diff0):
					g = i
			
			b = 800
			i = 800
			while (i > 300):
				diff0 = abs(math.dist((-1, -1), (round((wl_to_s(b, l1) - wl_to_s(b, s1)) / (wl_to_s(b, l1) + (wl_to_s(b, m1) + wl_to_s(b, s1))), 1), round((wl_to_s(b, m1) - wl_to_s(b, l1) - wl_to_s(b, s1)) / (wl_to_s(b, l1) + wl_to_s(b, m1) + wl_to_s(b, s1)), 1))))
				diff1 = abs(math.dist((-1, -1), (round((wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + (wl_to_s(i, m1) + wl_to_s(i, s1))), 1), round((wl_to_s(i, m1) - wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + wl_to_s(i, m1) + wl_to_s(i, s1)), 1))))
				if (diff1 < diff0):
					b = i
				i -= 1
			
			# secondaries
			white_l = 0
			white_m = 0
			white_s = 0
			for i in range(0, 50):
				w = i*10 + 340
				if (w >= args.cutoff):
					white_l += wl_to_s(w, l1) * wp[i]
					white_m += wl_to_s(w, m1) * wp[i]
					white_s += wl_to_s(w, s1) * wp[i]
			
			c = s1
			x1 = (wl_to_s(r, l1) - wl_to_s(r, s1)) / (wl_to_s(r, l1) + wl_to_s(r, m1) + wl_to_s(r, s1))
			y1 = (wl_to_s(r, m1) - (wl_to_s(r, l1) - wl_to_s(r, s1)) / wl_to_s(r, l1) + wl_to_s(r, m1) + wl_to_s(r, s1))
			x2 = (white_l - white_s) / (white_l + white_m + white_s)
			y2 = (white_m - white_l - white_s) / (white_l + white_m + white_s)
			
			for i in range(s1, m1):
				x3 = (wl_to_s(c, l1) - wl_to_s(c, s1)) / (wl_to_s(c, l1) + wl_to_s(c, m1) + wl_to_s(c, s1))
				y3 = (wl_to_s(c, m1) - (wl_to_s(c, l1) - wl_to_s(c, s1)) / wl_to_s(c, l1) + wl_to_s(c, m1) + wl_to_s(c, s1))
				x4 = (wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + wl_to_s(i, m1) + wl_to_s(i, s1))
				y4 = (wl_to_s(i, m1) - (wl_to_s(i, l1) - wl_to_s(i, s1)) / wl_to_s(i, l1) + wl_to_s(i, m1) + wl_to_s(i, s1))
				
				#d3 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x3, y3)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))
				#d4 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x4, y4)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))
				line = np.polyfit((x1, y1), (x2, y2), 1)
				A = line[0]
				B = line[1]
				C = A*x1 + B*x2
				
				d3 = abs(A*x3 + B*y3 + C) / math.sqrt(A**2 + B**2)
				d4 = abs(A*x4 + B*y4 + C) / math.sqrt(A**2 + B**2)
				
				if (d4 < d3):
					c = i
			
			y = m1
			x1 = (wl_to_s(b, l1) - wl_to_s(b, s1)) / (wl_to_s(b, l1) + wl_to_s(b, m1) + wl_to_s(b, s1))
			y1 = (wl_to_s(b, m1) - (wl_to_s(b, l1) - wl_to_s(b, s1)) / wl_to_s(b, l1) + wl_to_s(b, m1) + wl_to_s(b, s1))
			
			for i in range(m1, l1):
				x3 = (wl_to_s(y, l1) - wl_to_s(y, s1)) / (wl_to_s(y, l1) + wl_to_s(y, m1) + wl_to_s(y, s1))
				y3 = (wl_to_s(y, m1) - (wl_to_s(y, l1) - wl_to_s(y, s1)) / wl_to_s(y, l1) + wl_to_s(y, m1) + wl_to_s(y, s1))
				x4 = (wl_to_s(i, l1) - wl_to_s(i, s1)) / (wl_to_s(i, l1) + wl_to_s(i, m1) + wl_to_s(i, s1))
				y4 = (wl_to_s(i, m1) - (wl_to_s(i, l1) - wl_to_s(i, s1)) / wl_to_s(i, l1) + wl_to_s(i, m1) + wl_to_s(i, s1))
				
				#d3 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x3, y3)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))
				#d4 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x4, y4)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))
				
				line = np.polyfit((x1, y1), (x2, y2), 1)
				A = line[0]
				B = line[1]
				C = A*x1 + B*x2
				
				d3 = abs(A*x3 + B*y3 + C) / math.sqrt(A**2 + B**2)
				d4 = abs(A*x4 + B*y4 + C) / math.sqrt(A**2 + B**2)
				if (d4 < d3):
					y = i
			
			print("Estimated primary colors based on color triangle:")
			print("red: " + str(r))
			print("green: " + str(g))
			print("blue: " + str(b))
			print("cyan: " + str(c))
			print("yellow: " + str(y))
	
	# white
	# HSL comes up with a negative saturation. Either Python expects something different
	# from XYZ input or it just doesn't handle color spaces properly. Also note that for
	# both the "white" and "maple leaf" colors, the XYZ, LAB and LCh values specify an
	# impossibly high lightness/intensity (>100) whereas sRGB and HSL (derived from sRGB)
	# specify a very low value. (No, it's out of 1, not 100 or 255, so the RGB and HSL values
	# are also impossibly high.)
	# The negative saturation seems to happen when the input values are too high and the
	# brightness is more than 100%. I've scaled everything so the LMS values always add
	# up to 1, so we don't see this anymore. This means the lightness should be ignored.
	# We still have >100% lightness sometimes because XYZ Y is determined by the L and M
	# values and wider cone spacing means many visible wavelengths have a sky-high
	# L:M/S or M:L/S ratio that would be impossible for a human.
	
	# define these outside the if block so we can refer to them elsewhere
	
	if (args.lighting):
		print("White (E, equal-energy)")
		spectral_rendering(e/100)
		
		print("White (D65)")
		spectral_rendering(d65/100)
		
		print("Black body spectrum approximating 2000 K incandescent light bulb")
		spectral_rendering(incandescent)
	
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
	
	if (args.leaves):
		# reflectance spectrum of maple leaves estimated from this graph:
		# https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgKwT9GhKlvdbOROD_JtEEym7Ovzq8GPfSCORiVcKklEQI9DSTuNjoaIJMWMdpJpc4ijq0T1m_PXF2seWczauKLz-4VPIY9TSXQqXdp1B80vu4w5O4lWaAF0k2kaA5ThrJrjCzlck8Ez1fF/s1600/LeafSpectra.jpg
		# The graph has information for 300-400 but appears to be roughly 0 there.
		leaf_table = np.array([
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.01,
			0.01,
			0.02,
			0.02,
			0.02,
			0.03,
			0.03,
			0.03,
			0.04,
			0.04,
			0.05,
			0.08,
			0.19,
			0.23,
			0.25,
			0.27,
			0.25,
			0.23,
			0.19,
			0.17,
			0.17,
			0.15,
			0.14,
			0.14,
			0.10,
			0.07,
			0.05,
			0.04,
			0.10,
			0.30,
			0.35,
			0.45,
			0.55,
			0.60,
			0.63,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65
		])
		print("Green maple leaf")
		spectral_rendering(leaf_table)
		
		red_leaf_table = np.array([
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.01,
			0.01,
			0.02,
			0.02,
			0.03,
			0.05,
			0.06,
			0.08,
			0.10,
			0.10,
			0.09,
			0.08,
			0.05,
			0.10,
			0.30,
			0.35,
			0.45,
			0.55,
			0.60,
			0.63,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65,
			0.65
		])
		print("Red maple leaf")
		spectral_rendering(red_leaf_table)
		
		# corn
		# estimated from https://www.yorku.ca/planters/photosynthesis/2014_08_15_lab_manual_static_html/images/Corn_leaf_reflectance.png. Values below 400 are estimated to be the same as 400.
		corn_table = np.array([
			0.45,
			0.45,
			0.45,
			0.45,
			0.45,
			0.45,
			0.45,
			0.40,
			0.35,
			0.35,
			0.35,
			0.35,
			0.37,
			0.37,
			0.37,
			0.37,
			0.43,
			0.50,
			0.60,
			0.77,
			0.80,
			0.80,
			0.77,
			0.75,
			0.65,
			0.63,
			0.62,
			0.60,
			0.55,
			0.55,
			0.55,
			0.45,
			0.35,
			0.35,
			0.40,
			0.60,
			0.65,
			0.75,
			0.80,
			0.83,
			0.83,
			0.83,
			0.83,
			0.80,
			0.80,
			0.80,
			0.80,
			0.80,
			0.80,
			0.80
		])
		print("Corn leaf")
		spectral_rendering(corn_table)
		
		# hydrangea
		# estimated from https://spectralevolution.com/wp-content/uploads/2024/04/RT_hydrang_ref.jpg
		hydrangea_table = np.array([
			0.2,
			0.2,
			0.175,
			0.1375,
			0.125,
			0.1,
			0.0875,
			0.075,
			0.075,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.075,
			0.075,
			0.0875,
			0.1,
			0.1,
			0.0875,
			0.0875,
			0.075,
			0.075,
			0.075,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0625,
			0.0875,
			0.125,
			0.1875,
			0.3,
			0.375,
			0.45,
			0.4875,
			0.5125,
			0.5125,
			0.525,
			0.525,
			0.525,
			0.525,
			0.525,
			0.525
		])
		print("Hydrangea leaf")
		spectral_rendering(hydrangea_table)
	
	# flowers
	# estimated from https://www.researchgate.net/profile/Robert-Gegear/publication/222035355/figure/fig1/AS:632345413578753@1527774303378/The-spectral-reflectance-curves-of-coloured-flowers-and-the-green-array-background-used.png
	
	if (args.flowers):
		flower0_table = np.array([
			#0.05,
			#0.08,
			#0.08,
			#0.08,
			0.08,
			0.08,
			0.08,
			0.1,
			0.13,
			0.2,
			0.35,
			0.55,
			0.65,
			0.73,
			0.83,
			0.85,
			0.88,
			0.9,
			0.9,
			0.88,
			0.85,
			0.8,
			0.75,
			0.65,
			0.55,
			0.4,
			0.25,
			0.15,
			0.1,
			0.08,
			0.08,
			0.05,
			0.05,
			0.05,
			0.05,
			0.05,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08,
			0.08
		])
		print("Blue flower")
		spectral_rendering(flower0_table)
		
		flower1_table = np.array([
			#0.08,
			#0.08,
			#0.08,
			0.08,
			0.1,
			0.1,
			0.1,
			0.1,
			0.08,
			0.05,
			0.05,
			0.05,
			0.05,
			0.05,
			0.05,
			0.05,
			0.05,
			0.05,
			0.08,
			0.1,
			0.2,
			0.55,
			0.8,
			0.95,
			0.95,
			0.98,
			0.98,
			0.98,
			0.98,
			0.98,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1,
			0.1
		])
		print("Yellow flower")
		spectral_rendering(flower1_table)
		
		# Banksia attenuata flowers
		# estimated from fig. 5 in Arrese et al. 2005
		
		# immature
		banksia0_table = np.array([
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0625,
			0.125,
			0.125,
			0.125,
			0.125,
			0.125,
			0.25,
			0.375,
			0.5,
			0.625,
			0.75,
			0.75,
			0.6875,
			0.625,
			0.5625,
			0.625,
			0.5625,
			0.5,
			0.5,
			0.5,
			0.375,
			0.3125,
			0.1875,
			0.1875,
			0.625,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375,
			0.9375
		])
		print("Immature Banksia attenuata flower")
		spectral_rendering(banksia0_table)
		
		# mature
		banksia1_table = np.array([
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0625,
			0.125,
			0.125,
			0.1875,
			0.25,
			0.25,
			0.25,
			0.375,
			0.375,
			0.5,
			0.5,
			0.625,
			0.625,
			0.625,
			0.625,
			0.5625,
			0.625,
			0.625,
			0.625,
			0.625,
			0.625,
			0.625,
			0.625,
			0.5625,
			0.5625,
			0.625,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75,
			0.75
		])
		print("Mature Banksia attenuata flower")
		spectral_rendering(banksia1_table)
	
	# sky and daylight
	# These are both normalized to "1" at 550 nm. If I try to make "sky" relative to "daylight",
	# the values are way too low and probably not meaningful.
	if (args.sky):
		daylight_table = np.empty(50)
		for i in range(0, 50):
			w = i*10 + 340
			# normalize to 1 at 550
			daylight_table[i] = blackbody(w / 1000000000, 5800) / blackbody(550 / 1000000000, 5800)
		
		print("Black body spectrum approximating daylight")
		spectral_rendering(daylight_table)
		
		sky_table = np.empty(50)
		for i in range(0, 50):
			w = i*10 + 340
			sky_table[i] = blackbody(w / 1000000000, 5800) * w**-4 / (blackbody(550 / 1000000000, 5800) * 550**-4)
		
		print("Black body spectrum with Rayleigh scattering approximating a blue sky")
		spectral_rendering(sky_table)
		
	if (args.triangle):
		el = 0
		em = 0
		es = 0
		d65l = 0
		d65m = 0
		d65s = 0
		for i in range(0, 37):
			w = i*10 + 340
			el += wl_to_s(w, l1) * e[i]
			em += wl_to_s(w, m1) * e[i]
			es += wl_to_s(w, s1) * e[i]
			d65l += wl_to_s(w, l1) * d65[i]
			d65m += wl_to_s(w, m1) * d65[i]
			d65s += wl_to_s(w, s1) * d65[i]
			l = wl_to_s(w, l1) / (wl_to_s(w, l1) + wl_to_s(w, m1) + wl_to_s(w, s1))
			m = wl_to_s(w, m1) / (wl_to_s(w, l1) + wl_to_s(w, m1) + wl_to_s(w, s1))
			s = wl_to_s(w, s1) / (wl_to_s(w, l1) + wl_to_s(w, m1) + wl_to_s(w, s1))
			print(str(w) + " " + str(l) + " " + str(m) + " " + str(s))
		print("E " + str(el/(el+em+es)) + " " + str(em/(el+em+es)) + " " + str(es/(el+em+es)))
		print("D65 " + str(d65l/(d65l+d65m+d65s)) + " " + str(d65m/(d65l+d65m+d65s)) + " " + str(d65s/(d65l+d65m+d65s)))

# print execution time
print("%s seconds" % (time.time() - start_time))
