# This is an attempt to create something similar to Color Vision Simulator (Melin et al. 2013).
#
# References:
# A.D. Melin, D.W. Kline, C.M. Hickey, L.M. Fedigan, Food search through the eyes of a monkey: A functional substitution approach for assessing the ecology of primate color vision, Vision Research, Volume 86, 2013, Pages 87-96, ISSN 0042-6989, https://doi.org/10.1016/j.visres.2013.04.013. (https://www.sciencedirect.com/science/article/pii/S0042698913001119)
# Hogan JD, Fedigan LM, Hiramatsu C, Kawamura S, Melin AD. Trichromatic perception of flower colour improves resource detection among New World monkeys. Sci Rep. 2018 Jul 18;8(1):10883. doi: 10.1038/s41598-018-28997-4. PMID: 30022096; PMCID: PMC6052032. https://www.nature.com/articles/s41598-018-28997-4
# Melin AD, Chiou KL, Walco ER, Bergstrom ML, Kawamura S, Fedigan LM. Trichromacy increases fruit intake rates of wild capuchins (Cebus capucinus imitator). Proc Natl Acad Sci U S A. 2017 Sep 26;114(39):10402-10407. doi: 10.1073/pnas.1705957114. Epub 2017 Sep 11. PMID: 28894009; PMCID: PMC5625910. https://www.pnas.org/doi/10.1073/pnas.1705957114
# Fedigan LM, Melin AD, Addicott JF, Kawamura S. The heterozygote superiority hypothesis for polymorphic color vision is not supported by long-term fitness data from wild neotropical monkeys. PLoS One. 2014 Jan 3;9(1):e84872. doi: 10.1371/journal.pone.0084872. PMID: 24404195; PMCID: PMC3880319.
# Melin AD, Kline DW, Hiramatsu C, Caro T (2016) Zebra Stripes through the Eyes of Their Predators, Zebras, and Humans. PLoS ONE 11(1): e0145679. https://doi.org/10.1371/journal.pone.0145679
# Corso J, Bowler M, Heymann EW, Roos C, Mundy NI. Highly polymorphic colour vision in a New World monkey with red facial skin, the bald uakari (Cacajao calvus). Proc Biol Sci. 2016 Apr 13;283(1828):20160067. doi: 10.1098/rspb.2016.0067. PMID: 27053753; PMCID: PMC4843651. https://doi.org/10.1098/rspb.2016.0067
# Veilleux CC, Scarry CJ, Di Fiore A, Kirk EC, Bolnick DA, Lewis RJ. Group benefit associated with polymorphic trichromacy in a Malagasy primate (Propithecus verreauxi). Sci Rep. 2016 Dec 2;6:38418. doi: 10.1038/srep38418. PMID: 27910919; PMCID: PMC5133583. https://www.nature.com/articles/srep38418
# Hiramatsu C, Melin AD, Allen WL, Dubuc C, Higham JP. Experimental evidence that primate trichromacy is well suited for detecting primate social colour signals. Proc Biol Sci. 2017 Jun 14;284(1856):20162458. doi: 10.1098/rspb.2016.2458. PMID: 28615496; PMCID: PMC5474062. https://doi.org/10.1098/rspb.2016.2458
# Stavenga, D.G. On visual pigment templates and the spectral shape of invertebrate rhodopsins and metarhodopsins. J Comp Physiol A 196, 869–878 (2010). https://doi.org/10.1007/s00359-010-0568-7
# Hempel de Ibarra N, Vorobyev M, Menzel R. Mechanisms, functions and ecology of colour vision in the honeybee. J Comp Physiol A Neuroethol Sens Neural Behav Physiol. 2014 Jun;200(6):411-33. doi: 10.1007/s00359-014-0915-1. Epub 2014 May 15. PMID: 24828676; PMCID: PMC4035557.
# CMF From Cone Fundamentals. Horizon Lab @ UCRS. https://horizon-lab.org/colorvis/cone2cmf.html
# Petroc Sumner, Catherine A. Arrese, Julian C. Partridge; The ecology of visual pigment tuning in an Australian marsupial: the honey possum Tarsipes rostratus. J Exp Biol 15 May 2005; 208 (10): 1803–1815. doi: https://doi.org/10.1242/jeb.01610
# Arrese, A. C., Beazley, L. D. and Neumeyer, C. Behavioural evidence for marsupial trichromacy. doi:10.1016/j.cub.2006.02.036
# Ortín-Martínez A, Nadal-Nicolás FM, Jiménez-López M, Alburquerque-Béjar JJ, Nieto-López L, García-Ayuso D, Villegas-Pérez MP, Vidal-Sanz M, Agudo-Barriuso M. Number and distribution of mouse retinal cone photoreceptors: differences between an albino (Swiss) and a pigmented (C57/BL6) strain. PLoS One. 2014 Jul 16;9(7):e102392. doi: 10.1371/journal.pone.0102392. PMID: 25029531; PMCID: PMC4100816.
# Jacobs, G. H., Fenwick, J. A. and Williams, G. A. Cone-based vision of rats for ultraviolet and visible lights. Journal of Experimental Biology, 204(14), 15 July 2001. https://doi.org/10.1242/jeb.204.14.2439
# CIE 2-deg CMFs. http://cvrl.ucl.ac.uk/database/text/cienewxyz/cie2012xyz2.htm
# Pridmore, R.W. A new transformation of cone responses to opponent color responses. Atten Percept Psychophys 83, 1797–1803 (2021). https://doi.org/10.3758/s13414-020-02216-7
# https://chittkalab.sbcs.qmul.ac.uk/1992/Chittka%201992%20J%20Comp%20Physiol_R.pdf
#
# Further reading: https://xkcd.com/1926/

"""
10/10/2025 -- bringing the 24 May 2024 version back as a separate file so OpenCV doesn't have to share
space with matplotlib.

The methods I implemented are way messier than they need to be. If you really want to make a
false-color animal-vision image out of an ordinary photo, it should go something like this:

* find the spectral sensitivities of the camera, or at least an average/typical camera (PlotDigitizer to
the rescue)
* find the coefficients for those spectral sensitivities that provide the closest match to your
target phenotype
* undo the gamma function
* set the new r, g and b values according to the coefficients you found earlier
* redo gamma correction
* ???
* profit

micaToolbox can do something like this: https://www.empiricalimaging.com/knowledge-base/spectral-sensitivity-based-cone-catch-models/ (but much more complicated)
CVS doesn't, though. The parts I don't understand are (a) how they're getting the new saturation
if they completely separate it from the hue, (b) how they avoid doing hundreds of operations
for every pixel if they look up hues in a table containing "each wavelength in 1 nm increments".
Hash table?

The current version produces the desired intensity manually using the original HSV value. This is
obviously not as good as converting to a color space designed to be more perceptually uniform like
Lab or the LCH spaces, but those don't have linear hue so we lose the chromaticity. Converting between
RGB and HSL is better but distorts the saturation. Setting the CIE 2-deg or 10-deg fundamentals as the
source rather than digital camera channels allows for more similar output to real CVS images.
"""

import math
import numpy as np
import time
import argparse
import cv2
import central as c
import scipy
args = c.args

# execution time
start_time = time.time()

# I can't be bothered reading all that so I'll just put this down here

red = np.zeros(401)
green = np.zeros(401)
blue = np.zeros(401)
uv = np.zeros(401)

if (args.source == "camera"):
	# data for typical camera channels (Kolláth et al. 2020 https://doi.org/10.1016/j.jqsrt.2020.107162)
	# digitized with WebPlotDigitizer
	# I tried to pick a single curve that runs close to the middle, but they overlap a lot. Oh well.
	red_xy = c.csv2spec('red.csv')
	green_xy = c.csv2spec('green2.csv')
	blue_xy = c.csv2spec('blue.csv')
	red = np.interp(c.x_1nm, *red_xy)
	green = np.interp(c.x_1nm, *green_xy)
	blue = np.interp(c.x_1nm, *blue_xy)

	# remove extrapolation below ~400 nm
	for i in range(99):
		red[i] = 0
		green[i] = 0
		blue[i] = 0
elif (args.source == "cie2"):
	for i in range(401):
		w = i+300
		if (w < 390):
			red[i] = 0
			green[i] = 0
			blue[i] = 0
		else:
			red[i] = 10**(c.hcf2deg[int(round(w))-390][1])
			green[i] = 10**(c.hcf2deg[int(round(w))-390][2])
			blue[i] = 10**(c.hcf2deg[int(round(w))-390][3])
elif (args.source == "cie10"):
	for i in range(401):
		w = i+300
		if (w < 390):
			red[i] = 0
			green[i] = 0
			blue[i] = 0
		else:
			red[i] = 10**(c.hcf10deg[int(round(w))-390][1])
			green[i] = 10**(c.hcf10deg[int(round(w))-390][2])
			blue[i] = 10**(c.hcf10deg[int(round(w))-390][3])

# find best fit to target
def camera2target(xdata, rscale, gscale, bscale):
	ydata = np.empty(xdata.shape[0])
	for i in range(xdata.shape[0]):
		ydata[i] = rscale*red[i] + gscale*green[i] + bscale*blue[i]
	return(ydata)

wpred = sum(red * c.wp_1nm)
wpgreen = sum(green * c.wp_1nm)
wpblue = sum(blue * c.wp_1nm)
wpr1 = 1
wpr2 = 1
wpr3 = 1

r1 = np.zeros(401)
r2 = np.zeros(401)
r3 = np.zeros(401)
if (args.fcmode == 'cmf'):
	cmf = c.cmf()
	r1 = cmf[0]*c.wp_1nm
	if (c.dimension > 1): r2 = cmf[1]*c.wp_1nm
	if (c.dimension > 2): r3 = cmf[2]*c.wp_1nm
else:
	r1 = c.r1_1nm*c.media_1nm
	r2 = c.r2_1nm*c.media_1nm
	r3 = c.r3_1nm*c.media_1nm

	wpr1 = sum(r1 * c.wp_1nm)
	wpr2 = sum(r2 * c.wp_1nm)
	wpr3 = sum(r3 * c.wp_1nm)

# set lower bound to 0 because weird stuff happens if I let it choose negative numbers
popt, pcov = scipy.optimize.curve_fit(camera2target, c.x_1nm, r1, p0=[1, 1, 1], bounds=(0,np.inf))
print(popt)
red_r1 = popt[0]
green_r1 = popt[1]
blue_r1 = popt[2]
popt, pcov = scipy.optimize.curve_fit(camera2target, c.x_1nm, r2, p0=[1, 1, 1], bounds=(0,np.inf))
print(popt)
red_r2 = popt[0]
green_r2 = popt[1]
blue_r2 = popt[2]
popt, pcov = scipy.optimize.curve_fit(camera2target, c.x_1nm, r3, p0=[1, 1, 1], bounds=(0,np.inf))
print(popt)
red_r3 = popt[0]
green_r3 = popt[1]
blue_r3 = popt[2]
"""
plt.plot(c.x_1nm, r1, 'r')
plt.plot(c.x_1nm, r2, 'g')
plt.plot(c.x_1nm, r3, 'b')
plt.plot(c.x_1nm, red_r1*red + green_r1*green + blue_r1*blue, '--r')
plt.plot(c.x_1nm, red_r2*red + green_r2*green + blue_r2*blue, '--g')
plt.plot(c.x_1nm, red_r3*red + green_r3*green + blue_r3*blue, '--b')
plt.show()
"""
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

		# reconstruct original response based on illuminant
		rgb_r *= wpred
		rgb_g *= wpgreen
		rgb_b *= wpblue
		
		# combine
		rgb_r1 = red_r1*rgb_r + green_r1*rgb_g + blue_r1*rgb_b
		rgb_r2 = red_r2*rgb_r + green_r2*rgb_g + blue_r2*rgb_b
		rgb_r3 = red_r3*rgb_r + green_r3*rgb_g + blue_r3*rgb_b
		
		# compensate for illuminant
		if (wpr1 > 0): rgb_r1 /= wpr1
		if (wpr2 > 0): rgb_r2 /= wpr2
		if (wpr3 > 0): rgb_r3 /= wpr3
		
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
		
		# convert to chromaticity
		# This looks better if it's done after the gamma correction, probably
		# because the chromaticity is being used by sRGB.
		#total = rgb_r1 + rgb_r2 + rgb_r3
		#if (total > 0):
		#	rgb_r1 = rgb_r1 / total
		#	rgb_r2 = rgb_r2 / total
		#	rgb_r3 = rgb_r3 / total

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
