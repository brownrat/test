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

# convert hue to wavelength based on location of primary and secondary colors
# Originally I used primary colors based on an sRGB monitor display and my observations
# of white light through a spectroscope (620, 580, 550, 500, 470, 420). I've changed it
# to more closely match the CIE color space so these colors stand in for "perceptual red",
# etc. rather than the limitations of a monitor. This produces a result closer to both CVS
# output and the fat-tailed dunnart color space (Arrese et al.).
def hue_to_wavelength(hue):
	#print(hue)
	# red-yellow: 0-59 -> 650-571
	if (0 <= hue < 60):
		w = 650 - hue * (650 - 571)/60
	# yellow-green: 60-119 -> 570-551
	elif (60 <= hue < 120):
		w = 570 - (hue - 60) * (570 - 526)/60
	# green-cyan: 120-179 -> 525-491
	elif (120 <= hue < 180):
		w = 525 - (hue - 120) * (525 - 491)/60
	# cyan-blue: 180-239 -> 490-471
	elif (180 <= hue < 240):
		w = 490 - (hue - 180) * (490 - 471)/60
	# blue-violet: 240-269 -> 470-420
	elif (240 <= hue < 270):
		w = 470 - (hue - 240) * (470 - 420)/30
	# purple-magenta: 270-360 -- just convert it to RGB and use the R and B values
	else:
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)
#	print(w)
	return w

# find sensitivity of a given cone type to a wavelength
# This does not use the same function as CVS. The function comes from a paper describing
# templates for visual pigment absorbance that are meant to fit both vertebrates and
# invertebrates (Stavenga 2010), whereas CVS uses shifted versions of the 10° human cone
# fundamentals.
def wl_to_s(wl, peak):
	value = (math.exp(69.7*(0.8795+0.0459*math.exp(-(peak-300)**2/11940)-(peak/wl)))+math.exp(28*(0.922-(peak/wl)))+math.exp(-14.9*(1.104-(peak/wl)))+0.674)**-1 + 0.26*math.exp(-((wl-(189+0.315*peak))/(-40.5+0.195*peak))**2)
	return value

# convert wavelength to RGB values for target color space
# How CVS does this is described as follows:
# "Using the idealized cone fundamentals, and based on standard daylight (D65) illumination, the program generates a color space table for the source and target phenotypes that lists the relative sensitivities of the component photopigments to each wavelength in 1 nm increments – intermediate values are derived from the table via linear interpolation. For each image in the trial, the starting RGB values for each pixel are converted to chromaticity and intensity (luminance) values. The chromaticity values are located in the source color space table to find the predominant hue (wavelength); the saturation value is determined by the relative distance of the chromaticity value from the white point. The hue value is then located in the target color space table, and the target saturation value is projected from the white point, to find the modified chromaticity values. An inverse transformation to the new RGB color space generates the task stimulus. The luminance values are held constant during this process. In the case of monochromatic simulations, only the luminance values of each pixel are used to generate gray-scale images."
# So far I haven't been able to figure out how to construct a color space and derive anything
# from it. One method is described here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3056573/

# version 0: use unaltered cone fundamentals
# This just takes the sensitivity values for all three cones and uses them as the R, G and B
# values. For wide spacing (> 40 nm) as in birds, insects and some marsupials, where the
# primary wavelengths are received almost exclusively by one type of cone, this is probably
# pretty close and may be better than the other versions. For narrower spacing as in primates,
# it's way off because the interactions become more complex.
def wl_to_rgb_0(wl, l, m, s):
	red = wl_to_s(wl, l)
	green = wl_to_s(wl, m)
	blue = wl_to_s(wl, s)
	return [red, green, blue]

# version 1: use color matching functions
# This was the first attempt to create a proper color space, based on the equations on the
# Horizon Lab page. It doesn't work for reasons I don't entirely understand. The most obvious
# problem is it has no way to find the primary wavelengths.
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
	lp = wl_to_s(w, 565) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
	mp = wl_to_s(w, 540) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
	sp = wl_to_s(w, 440) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
	
	p0 = 800
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
		and p1 < 800):
		p1 += 1
	# lower precision if not found
	if (p0 == 300 and p1 == 800):
		# recursion
		p2 = find_match(w, l, m, s, d*2) # avoid running this twice
		p0 = p2[0]
		p1 = p2[1]
	
	#print(p0)
	#print(p1)
	
	return([p0, p1])

def find_primaries(l, m, s):
	# find "red" L primary
	red = (find_match(650, l, m, s)[0] + find_match(650, l, m, s)[1]) / 2
	
	# find "yellow" LM secondary
	yellow = (find_match(570, l, m, s)[0] + find_match(570, l, m, s)[1]) / 2
	
	# find "green" M primary
	green = (find_match(525, l, m, s)[0] + find_match(525, l, m, s)[1]) / 2
	
	# find "cyan" MS secondary
	cyan = (find_match(490, l, m, s)[0] + find_match(490, l, m, s)[1]) / 2
	
	# find "blue" S primary
	blue = (find_match(470, l, m, s)[0] + find_match(470, l, m, s)[1]) / 2
	
	# find "violet" LS secondary
	violet = (find_match(420, l, m, s)[0] + find_match(420, l, m, s)[1]) / 2
	
	# patch results if they're out of order
	if (red < yellow):
		red = 800
	if (cyan < blue):
		cyan = (green + blue) / 2
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
# shifting of red/yellow (in either direction) seems to be underestimated.
def cs_tables(l, m, s):
	# source color space table
	cs_source = np.empty([270, 3])
	for i in range(cs_source.shape[0]):
		w = hue_to_wavelength(i)
		cs_source[i][0] = wl_to_s(w, 565) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
		cs_source[i][1] = wl_to_s(w, 540) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
		cs_source[i][2] = wl_to_s(w, 440) / (wl_to_s(w, 565) + wl_to_s(w, 540) + wl_to_s(w, 440))
	
	# target color space table (doesn't work, I don't know why)
	cs_target = np.empty([501, 3])
	for i in range(cs_target.shape[0]):
		w = i + 300
		cs_target[i][0] = wl_to_s(w, l) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
		cs_target[i][1] = wl_to_s(w, m) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
		cs_target[i][2] = wl_to_s(w, s) / (wl_to_s(w, l) + wl_to_s(w, m) + wl_to_s(w, s))
	
	return [cs_source, cs_target]
	
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
	j = 269
	
	# forward
	while (abs(source[i][0] - lw)
		+ abs(source[i][1] - mw)
		+ abs(source[i][2] - sw) > d
		and i < 269):
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
	if (i == 269 and j == 0):
		#print("recursion")
		# recursion
		match = find_match_1(w, source, l, m, s, d*2) # avoid running this twice
		i = match[0]
		j = match[1]
	
	return([i, j])

# Convert image from BGR to HLS. We use BGR because this is how OpenCV reads images.
# If we use RGB, the output appears normal with settings approximating human vision,
# but shifting the cones produces the opposite of the expected result.
img = cv2.imread(imagename)
img_hls = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)

# switch method used based on input -- this function calls the desired version of wl_to_rgb
# default is 2 (cone response ratios)

# save some time by taking these out of wl_to_rgb so we don't call them for every pixel
if (version == 2):
	primaries = find_primaries(lcone, mcone, scone)
elif (version == 3):
	tables = cs_tables(lcone, mcone, scone)

# wrapper function
def wl_to_rgb(w, l, m, s):
	if (version == 1):
		return wl_to_rgb_1(w, l, m, s)
	elif (version == 2):
		return wl_to_rgb_2(w, l, m, s, primaries)
	elif (version == 3):
		match = find_match_1(w, tables[0], l, m, s)
		hue = (match[0] + match[1]) / 2
		return colorsys.hls_to_rgb(hue / 360, 0.5, 1)
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
		
		if (type(wl) != tuple):
			hue_target = wl_to_rgb(wl, lcone, mcone, scone)
		else:
			hue_r = wl_to_rgb(650, lcone, mcone, scone)
			hue_b = wl_to_rgb(470, lcone, mcone, scone)
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
