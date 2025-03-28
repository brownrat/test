# This began as an attempt to create something similar to Color Vision Simulator (Melin et al.
# 2013) but now tests general properties of a specified color vision system and perception
# of various spectral power distributions. The image processing functions were no longer
# usable because OpenCV is not compatible with Matplotlib. See https://forum.qt.io/topic/119109/using-pyqt5-with-opencv-python-cv2-causes-error-could-not-load-qt-platform-plugin-xcb-even-though-it-was-found
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
# Pridmore, R.W. A new transformation of cone responses to opponent color responses. Atten Percept Psychophys 83, 1797–1803 (2021). https://doi.org/10.3758/s13414-020-02216-7
# https://chittkalab.sbcs.qmul.ac.uk/1992/Chittka%201992%20J%20Comp%20Physiol_R.pdf
# GOVARDOVSKII VI, FYHRQUIST N, REUTER T, KUZMIN DG, DONNER K. In search of the visual pigment template. Visual Neuroscience. 2000;17(4):509-528. doi:10.1017/S0952523800174036
# https://journals.biologists.com/jeb/article/207/14/2471/14748/Interspecific-and-intraspecific-views-of-color
# T.D. Lamb, Photoreceptor spectral sensitivities: Common shape in the long-wavelength region, Vision Research, Volume 35, Issue 22, 1995, Pages 3083-3091, ISSN 0042-6989, https://doi.org/10.1016/0042-6989(95)00114-F. (https://www.sciencedirect.com/science/article/pii/004269899500114F)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3652215/
# https://royalsocietypublishing.org/doi/10.1098/rspb.1985.0045
# https://labs.mcdb.ucsb.edu/fisher/steven/pubs/Fisher,Pfeffer&Anderson-1983.pdf
# Nikonov SS, Kholodenko R, Lem J, Pugh EN Jr. Physiological features of the S- and M-cone photoreceptors of wild-type mice from single-cell recordings. J Gen Physiol. 2006 Apr;127(4):359-74. doi: 10.1085/jgp.200609490. PMID: 16567464; PMCID: PMC2151510.
#
# Further reading: https://xkcd.com/1926/

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
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binomtest
import statistics
from scipy.integrate import quad

# execution time
start_time = time.time()

# arguments
# 07/03/2025 -- replaced some default values:
# * weberb: 0.14 -> 0.11. 0.11 is the value determined for humans in a study of manatees (https://pubmed.ncbi.nlm.nih.gov/9202447/).
# 0.14 was found in a previous study (Cornsweet & Pinsker), but the number seems to only appear in citations
# and not in the actual paper. All I can find in the paper is graphs.
# * lp/mp: 9/4.5 -> 10/5. The ratio of L to S cones in Didelphis aurita and other marsupials is closer to 10:1 than
# 9:1. The L:M ratio is kept at 2:1.
# * lb/mb: 0.68990272/0.34832189 -> 1/0. The former values are for humans. Since we're now modeling human luminance vision with the CIE
# luminosity function, we don't need this.
# Earlier I also changed the chromatic Weber fraction from the human value of 0.018 to the semi-standard estimate of 0.05
# because we're not trying to model human color vision with RNL and we're not set up properly for that anyway.
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--lw", type=int, default=0, help="longwave cone sensitivity")
parser.add_argument("-m", "--mw", type=int, default=0, help="mediumwave cone sensitivity")
parser.add_argument("-s", "--sw", type=int, default=0, help="shortwave cone sensitivity")
parser.add_argument("--rod", type=int, default=0, help="rod sensitivity")
parser.add_argument("--weber", type=float, default=0.05, help="Weber fraction for L cones") # formerly 0.018
parser.add_argument("--weberm", type=float, default=0, help="override Weber fraction for M cones")
parser.add_argument("--webers", type=float, default=0, help="override Weber fraction for S cones")
parser.add_argument("--weberb", type=float, default=0.11, help="Weber fraction for brightness")
parser.add_argument("--lp", type=float, default=10, help="number of L cones per receptive field")
parser.add_argument("--mp", type=float, default=5, help="number of M cones per receptive field")
parser.add_argument("--sp", type=float, default=1, help="number of S cones per receptive field")
parser.add_argument("--filter", type=str, default="none", help="type of lens filtering")
parser.add_argument("--qn", help="use quantum noise in color differences", action="store_true")
parser.add_argument("--lb", type=float, default=1, help="contribution of L cones to perception of brightness") # formerly 0.68990272
parser.add_argument("--mb", type=float, default=0, help="contribution of M cones to perception of brightness")
parser.add_argument("--sb", type=float, default=0, help="contribution of S cones to perception of brightness")
parser.add_argument("--rb", type=float, default=0, help="contribution of rods to perception of brightness")
parser.add_argument("--novk", help="disable von Kries transform", action="store_true")
parser.add_argument("--luminosity", help="show wavelength with maximum sensitivity", action="store_true")
parser.add_argument("--primaries", help="show primary and secondary wavelengths", action="store_true")
parser.add_argument("--white", default="e", help="reference white")
parser.add_argument("--lighting", help="standard light sources", action="store_true")
parser.add_argument("--leaves", help="several types of leaves", action="store_true")
parser.add_argument("--flowers", help="several types of flowers", action="store_true")
parser.add_argument("--step", help="long-pass step-function spectra", action="store_true")
parser.add_argument("--steps", help="short-pass step-function spectra", action="store_true")
parser.add_argument("--peak", help="band-pass spectra", action="store_true")
parser.add_argument("--ramp", help="ramp-function spectra", action="store_true")
parser.add_argument("--sky", help="sky-like and sun-like spectra", action="store_true")
parser.add_argument("--kodak", help="Kodak Wratten camera filters", action="store_true")
parser.add_argument("--kv2", help="use tabulated Kodak Wratten data from Kodak Photographic Filters Handbook (1990)", action="store_true")
parser.add_argument("--kcheck", help="check camera filter brightness matches", action="store_true")
parser.add_argument("--render", help="show colormath spectral rendering for comparison", action="store_true")
parser.add_argument("--qcheck", help="show absolute quantum catches", action="store_true")
parser.add_argument("--vpt", help="plot visual pigment templates", action="store_true")
parser.add_argument("--triangle", help="plot color triangle", action="store_true")
parser.add_argument("--triangle1", help="plot logarithmic noise-scaled color triangle", action="store_true")
parser.add_argument("--hexagon", help="plot color hexagon", action="store_true")
parser.add_argument("--wd", help="plot wavelength discrimination function", action="store_true")
parser.add_argument("--blackbody", type=int, default=0, help="plot black body curve for a given temperature in Kelvin")
parser.add_argument("--munsell", help="Munsell color cards", action="store_true")
parser.add_argument("--mv2", help="use alternative Munsell spectra", action="store_true")
parser.add_argument("--mv3", help="use alternative Munsell spectra for <420nm only", action="store_true")
parser.add_argument("--mv4", help="remove wavelengths <400nm from Munsell spectra", action="store_true")
parser.add_argument("--mopt1", help="find optimal visual pigments for Munsell cards (dichromacy)", action="store_true")
parser.add_argument("--mopt2", help="find optimal visual pigments for Munsell cards (trichromacy)", action="store_true")
parser.add_argument("--mopt3", help="find L-M chromatic spread/overlap of Munsell cards", action="store_true")
parser.add_argument("--erg", help="show example of ERG (electroretinography) graph adjustment", action="store_true")
parser.add_argument("--convergence", type=float, default=1, help="cone-to-ganglion cell convergence ratio")
parser.add_argument("--ff", type=float, default=22.4, help="flicker fusion frequency (Hz)")
parser.add_argument("--osd", type=float, default=1.2, help="outer segment diameter (um)")
parser.add_argument("--fl", type=float, default=6.81, help="focal length (mm)")
parser.add_argument("--pupil", type=float, default=6, help="pupil diameter (mm)")
parser.add_argument("--od", type=float, default=0.015, help="optical density of pigment (um^-1)")
parser.add_argument("--osl", type=float, default=30, help="outer segment length (um)")
parser.add_argument("--dist", type=float, default=5, help="distance from light source (cm)")
parser.add_argument("--width", type=float, default=3.8, help="width of light/object (cm)")
parser.add_argument("--height", type=float, default=0, help="height of light/object (cm)")
parser.add_argument("--radius", type=float, default=1, help="radius of light source (cm)")
parser.add_argument("--watts", type=int, default=12, help="number of watts produced by light source")
parser.add_argument("--lux", type=int, default=115, help="illuminance in lux")
parser.add_argument("--ct", type=int, default=2856, help="color temperature of incandescent/blackbody illuminant")
parser.add_argument("--selfscreen", help="enable pigment self-screening (varying optical density)", action="store_true")
parser.add_argument("--selfscreen1", help="enable pigment self-screening (simple version)", action="store_true")
parser.add_argument("--warnings", help="show warnings", action="store_true")
parser.add_argument("--cmf", help="color-matching functions", action="store_true")
parser.add_argument("-r", "--red", type=float, default=700, help="'red' primary")
parser.add_argument("-g", "--green", type=float, default=546.1, help="'green' primary")
parser.add_argument("-b", "--blue", type=float, default=435.8, help="'blue' primary")
parser.add_argument("--q10deg", help="use 10-degree instead of 2-degree cone fundamentals", action="store_true")
parser.add_argument("--e2deg", help="use 2-degree energy cone fundamentals", action="store_true")
args = parser.parse_args()

# cone and rod peak sensitivities
l1 = args.lw
m1 = args.mw
s1 = args.sw
rod = args.rod
if (l1 == 0 and m1 != 0):
	l1 = m1
if (m1 == 0 and l1 != 0):
	m1 = l1
if (s1 == 0 and m1 != 0):
	s1 = m1

# Weber fractions
wl = args.weber
wm = wl * math.sqrt(args.lp) / math.sqrt(args.mp)
ws = wl * math.sqrt(args.lp) / math.sqrt(args.sp)
if (args.weberm > 0):
	wm = args.weberm
if (args.webers > 0):
	ws = args.webers

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
# It can only tell us the direction of a difference in brightness, not how large or
# perceptible it is -- for that you need Weber fractions.
# The percentages of L, M and S cones are about 51-76%, 20-44% and 2% in humans (Wikipedia)
# and 68-73%, 20-25% and 7% in the fat-tailed dunnart ("Diversity of Color* Vision: Not All
# Australian Marsupials are Trichromatic"). In Norway rats, the percentage of S cones is 11-12%
# ("Cone-based vision of rats for ultraviolet and visible lights") and in the house mouse the
# percentage of exclusive non-coexpressing S cones is 26% ("Number and Distribution of Mouse
# Retinal Cone Photoreceptors").
# * normally spelled "colour" in Australia, but this is the correct title.

# function approximating human lens filtering, for real this time (Lamb 1995)
# 21/03/2025 -- optical density uses 10 not e, again.
def template_filter(w, a, b, c, d, e, f):
	try:
		density = a*math.exp((b - w) / c) + d*math.exp((e - w) / f)
	except OverflowError:
		if (args.warnings): print("Warning (template_filter): math overflow, returning 0")
		return 0
	if (density < 0):
		return 0
	return 10**(-density)
def human_filter(w):
	return template_filter(w, 1.1, 400, 15, 0.11, 500, 80)

# array version for curve fitting
def filter_fit(xdata, a, b, c, d, e, f):
	ydata = np.empty(xdata.shape[0])
	for i in range(xdata.shape[0]):
		ydata[i] = template_filter(xdata[i], a, b, c, d, e, f)
	return(ydata)

# mouse (Jacobs & Williams 2007)
# Transmission below 310 nm is not provided and has been set to 0. This means we're effectively
# only considering wavelengths >=310, but it's easier to keep track of the numbers if they begin
# with a multiple of 100.
mouse_lens = np.array([
	0, # 300
	0, #309
	10.3, # 310
	33.0, # 320
	44.9, # 330
	51.9, # 340
	58.1, # 350
	63.8, # 360
	68.0, # 370
	71.9, # 380
	75.2, # 390
	78.1, # 400
	80.3, # 410
	82.0, # 420
	83.5, # 430
	85.1, # 440
	86.4, # 450
	87.6, # 460
	88.7, # 470
	89.6, # 480
	90.2, # 490
	91.3, # 500
	92.1, # 510
	92.8, # 520
	93.4, # 530
	93.9, # 540
	94.5, # 550
	95.0, # 560
	95.5, # 570
	96.1, # 580
	96.6, # 590
	97.1, # 600
	97.6, # 610
	97.9, # 620
	98.1, # 630
	98.5, # 640
	98.8, # 650
	99.2, # 660
	99.4, # 670
	99.6, # 680
	99.8, # 690
	100.0 # 700
])

# interpolate
x_10nm = np.empty(42)
x_10nm[0] = 300
x_10nm[1] = 309
for i in range(2, 42):
	x_10nm[i] = (i-1)*10 + 300
x_1nm = np.empty(401)
for i in range(401): x_1nm[i] = i + 300

mouse_lens_1nm = np.interp(x_1nm, x_10nm, mouse_lens/100)
mouse_lens_1nm_squared = np.interp(x_1nm, x_10nm, (mouse_lens/100)**2)

def lens_filter(w):
	if (args.filter == "human"):
		return human_filter(w)
	elif (args.filter == "mouse"):
		w = round(w)
		return mouse_lens_1nm[w-300]
	elif (args.filter == "mousesquared"):
		w = round(w)
		return mouse_lens_1nm_squared[w-300]
	return 1

def sensitivity(w, l=l1, m=m1, s=s1):
	if (args.filter == "human2"):
		return luminosity[w-300]
	return (args.lb*vpt(w, l) + args.mb*vpt(w, m) + args.sb*vpt(w, s) + args.rb*vpt(w, rod)) * lens_filter(w)

# black body -- SI units, returns energy in joules/m^2/nm/sec/steradian = watts/m^2/nm/steradian
# The "1e-9" converts from meters to nanometers, as this is not cubic meters but "per
# wavelength".
# https://yceo.yale.edu/sites/default/files/files/ComputingThePlanckFunction.pdf
h = 6.626068e-34 # Planck's constant (joule sec)
k = 1.38066e-23 # Boltzmann's constant (joule/deg)
c = 2.997925e+8 # speed of light (m/s)

def blackbody(w, t):
	l = w / 1000000000 # wavelength (m)
	try:
		value = 1e-9*2*h*c**2*l**-5 / (math.exp((h*c) / (k*l*t)) - 1)
	except OverflowError:
		if (args.warnings): print("Warning (blackbody): math overflow, returning 1")
		return 1
	return value

# normal distribution
def normal(mu, std_dev, x):
	return (1 / std_dev*math.sqrt(2*math.pi)) * math.exp((-1/2)*((x - mu)/std_dev)**2)

# illuminants

# D65 300-830 (CIE)
d65_1nm = np.array([
	0.0341,
	0.36014,
	0.68618,
	1.01222,
	1.33826,
	1.6643,
	1.99034,
	2.31638,
	2.64242,
	2.96846,
	3.2945,
	4.98865,
	6.6828,
	8.37695,
	10.0711,
	11.7652,
	13.4594,
	15.1535,
	16.8477,
	18.5418,
	20.236,
	21.9177,
	23.5995,
	25.2812,
	26.963,
	28.6447,
	30.3265,
	32.0082,
	33.69,
	35.3717,
	37.0535,
	37.343,
	37.6326,
	37.9221,
	38.2116,
	38.5011,
	38.7907,
	39.0802,
	39.3697,
	39.6593,
	39.9488,
	40.4451,
	40.9414,
	41.4377,
	41.934,
	42.4302,
	42.9265,
	43.4228,
	43.9191,
	44.4154,
	44.9117,
	45.0844,
	45.257,
	45.4297,
	45.6023,
	45.775,
	45.9477,
	46.1203,
	46.293,
	46.4656,
	46.6383,
	47.1834,
	47.7285,
	48.2735,
	48.8186,
	49.3637,
	49.9088,
	50.4539,
	50.9989,
	51.544,
	52.0891,
	51.8777,
	51.6664,
	51.455,
	51.2437,
	51.0323,
	50.8209,
	50.6096,
	50.3982,
	50.1869,
	49.9755,
	50.4428,
	50.91,
	51.3773,
	51.8446,
	52.3118,
	52.7791,
	53.2464,
	53.7137,
	54.1809,
	54.6482,
	57.4589,
	60.2695,
	63.0802,
	65.8909,
	68.7015,
	71.5122,
	74.3229,
	77.1336,
	79.9442,
	82.7549,
	83.628,
	84.5011,
	85.3742,
	86.2473,
	87.1204,
	87.9936,
	88.8667,
	89.7398,
	90.6129,
	91.486,
	91.6806,
	91.8752,
	92.0697,
	92.2643,
	92.4589,
	92.6535,
	92.8481,
	93.0426,
	93.2372,
	93.4318,
	92.7568,
	92.0819,
	91.4069,
	90.732,
	90.057,
	89.3821,
	88.7071,
	88.0322,
	87.3572,
	86.6823,
	88.5006,
	90.3188,
	92.1371,
	93.9554,
	95.7736,
	97.5919,
	99.4102,
	101.228,
	103.047,
	104.865,
	106.079,
	107.294,
	108.508,
	109.722,
	110.936,
	112.151,
	113.365,
	114.579,
	115.794,
	117.008,
	117.088,
	117.169,
	117.249,
	117.33,
	117.41,
	117.49,
	117.571,
	117.651,
	117.732,
	117.812,
	117.517,
	117.222,
	116.927,
	116.632,
	116.336,
	116.041,
	115.746,
	115.451,
	115.156,
	114.861,
	114.967,
	115.073,
	115.18,
	115.286,
	115.392,
	115.498,
	115.604,
	115.711,
	115.817,
	115.923,
	115.212,
	114.501,
	113.789,
	113.078,
	112.367,
	111.656,
	110.945,
	110.233,
	109.522,
	108.811,
	108.865,
	108.92,
	108.974,
	109.028,
	109.082,
	109.137,
	109.191,
	109.245,
	109.3,
	109.354,
	109.199,
	109.044,
	108.888,
	108.733,
	108.578,
	108.423,
	108.268,
	108.112,
	107.957,
	107.802,
	107.501,
	107.2,
	106.898,
	106.597,
	106.296,
	105.995,
	105.694,
	105.392,
	105.091,
	104.79,
	105.08,
	105.37,
	105.66,
	105.95,
	106.239,
	106.529,
	106.819,
	107.109,
	107.399,
	107.689,
	107.361,
	107.032,
	106.704,
	106.375,
	106.047,
	105.719,
	105.39,
	105.062,
	104.733,
	104.405,
	104.369,
	104.333,
	104.297,
	104.261,
	104.225,
	104.19,
	104.154,
	104.118,
	104.082,
	104.046,
	103.641,
	103.237,
	102.832,
	102.428,
	102.023,
	101.618,
	101.214,
	100.809,
	100.405,
	100,
	99.6334,
	99.2668,
	98.9003,
	98.5337,
	98.1671,
	97.8005,
	97.4339,
	97.0674,
	96.7008,
	96.3342,
	96.2796,
	96.225,
	96.1703,
	96.1157,
	96.0611,
	96.0065,
	95.9519,
	95.8972,
	95.8426,
	95.788,
	95.0778,
	94.3675,
	93.6573,
	92.947,
	92.2368,
	91.5266,
	90.8163,
	90.1061,
	89.3958,
	88.6856,
	88.8177,
	88.9497,
	89.0818,
	89.2138,
	89.3459,
	89.478,
	89.61,
	89.7421,
	89.8741,
	90.0062,
	89.9655,
	89.9248,
	89.8841,
	89.8434,
	89.8026,
	89.7619,
	89.7212,
	89.6805,
	89.6398,
	89.5991,
	89.4091,
	89.219,
	89.029,
	88.8389,
	88.6489,
	88.4589,
	88.2688,
	88.0788,
	87.8887,
	87.6987,
	87.2577,
	86.8167,
	86.3757,
	85.9347,
	85.4936,
	85.0526,
	84.6116,
	84.1706,
	83.7296,
	83.2886,
	83.3297,
	83.3707,
	83.4118,
	83.4528,
	83.4939,
	83.535,
	83.576,
	83.6171,
	83.6581,
	83.6992,
	83.332,
	82.9647,
	82.5975,
	82.2302,
	81.863,
	81.4958,
	81.1285,
	80.7613,
	80.394,
	80.0268,
	80.0456,
	80.0644,
	80.0831,
	80.1019,
	80.1207,
	80.1395,
	80.1583,
	80.177,
	80.1958,
	80.2146,
	80.4209,
	80.6272,
	80.8336,
	81.0399,
	81.2462,
	81.4525,
	81.6588,
	81.8652,
	82.0715,
	82.2778,
	81.8784,
	81.4791,
	81.0797,
	80.6804,
	80.281,
	79.8816,
	79.4823,
	79.0829,
	78.6836,
	78.2842,
	77.4279,
	76.5716,
	75.7153,
	74.859,
	74.0027,
	73.1465,
	72.2902,
	71.4339,
	70.5776,
	69.7213,
	69.9101,
	70.0989,
	70.2876,
	70.4764,
	70.6652,
	70.854,
	71.0428,
	71.2315,
	71.4203,
	71.6091,
	71.8831,
	72.1571,
	72.4311,
	72.7051,
	72.979,
	73.253,
	73.527,
	73.801,
	74.075,
	74.349,
	73.0745,
	71.8,
	70.5255,
	69.251,
	67.9765,
	66.702,
	65.4275,
	64.153,
	62.8785,
	61.604,
	62.4322,
	63.2603,
	64.0885,
	64.9166,
	65.7448,
	66.573,
	67.4011,
	68.2293,
	69.0574,
	69.8856,
	70.4057,
	70.9259,
	71.446,
	71.9662,
	72.4863,
	73.0064,
	73.5266,
	74.0467,
	74.5669,
	75.087,
	73.9376,
	72.7881,
	71.6387,
	70.4893,
	69.3398,
	68.1904,
	67.041,
	65.8916,
	64.7421,
	63.5927,
	61.8752,
	60.1578,
	58.4403,
	56.7229,
	55.0054,
	53.288,
	51.5705,
	49.8531,
	48.1356,
	46.4182,
	48.4569,
	50.4956,
	52.5344,
	54.5731,
	56.6118,
	58.6505,
	60.6892,
	62.728,
	64.7667,
	66.8054,
	66.4631,
	66.1209,
	65.7786,
	65.4364,
	65.0941,
	64.7518,
	64.4096,
	64.0673,
	63.7251,
	63.3828,
	63.4749,
	63.567,
	63.6592,
	63.7513,
	63.8434,
	63.9355,
	64.0276,
	64.1198,
	64.2119,
	64.304,
	63.8188,
	63.3336,
	62.8484,
	62.3632,
	61.8779,
	61.3927,
	60.9075,
	60.4223,
	59.9371,
	59.4519,
	58.7026,
	57.9533,
	57.204,
	56.4547,
	55.7054,
	54.9562,
	54.2069,
	53.4576,
	52.7083,
	51.959,
	52.5072,
	53.0553,
	53.6035,
	54.1516,
	54.6998,
	55.248,
	55.7961,
	56.3443,
	56.8924,
	57.4406,
	57.7278,
	58.015,
	58.3022,
	58.5894,
	58.8765,
	59.1637,
	59.4509,
	59.7381,
	60.0253,
	60.3125
])

d65 = np.zeros(41)
for i in range(0, 41):
	d65[i] = d65_1nm[i*10]
	#print(d65[i])

e = np.empty(41)
for i in range(0, 41):
	e[i] = 100
e_1nm = np.empty(401)
for i in range(0, 401):
	e_1nm[i] = 100

a = np.zeros(41)
for i in range(4, 41):
	a[i] = spectral_constants.REF_ILLUM_TABLE["a"][i-4]

# incandescent lighting with specified color temperature
incandescent = np.empty(41)
for i in range(0, 41):
	w = i*10 + 300
	# normalize to 100% at 560 nm
	incandescent[i] = 100*blackbody(w, args.ct) / blackbody(560, args.ct)
# 1-nm resolution
incandescent_1nm = np.empty(401)
for i in range(401):
	w = i + 300
	# normalize to 100% at 560 nm
	incandescent_1nm[i] = 100*blackbody(w, args.ct) / blackbody(560, args.ct)

# absolute number of photons
# Since the output of blackbody() is joules/m^3/sec/sr, this should be photons/m^3/sec/sr.
# I think what we really want is square meters. (No, the meters cancel, see the quantum
# noise calculations)
# Changing this to 1-nm intervals because 10-nm isn't very useful.
ia = np.empty(501)
# scale it down to what amount of energy would be produced by the specified number of
# watts

# find total energy produced using Stefan-Boltzmann constant
# We divide by pi steradians to get the radiance and make the units equivalent to blackbody()
# (not 4pi, this gives the wrong value -- see p. 61 here: https://eodg.atm.ox.ac.uk/user/grainger/research/book/protected/Chapter3.pdf)
sigma = 5.670374419e-8
energy = sigma * args.ct**4 #/ math.pi
#power = energy * 4*math.pi*args.dist**2 / 10000

# then find the radiance we "should" have in terms of the sphere corresponding to
# the distance we chose
sphere_area = 4*math.pi*args.dist**2 # sphere surface area
if (args.height > 0):
	area = args.width * args.height
else:
	area = math.pi*(args.width/2)**2
# steradians -- not divided by the whole sphere surface area, just the square of the radius
sr = area / args.dist**2
watts_cm2 = args.watts / (sphere_area)
watts_m2 = watts_cm2 * 100**2
scale = watts_m2 / energy # scaling factor -- units should cancel (W/m^2)
#scale = args.watts / power
#print("W/m^2 scaling factor: " + str(scale))

for i in range(501):
	w = i + 300 # nanometers
	ia[i] = scale*1e-9*w*blackbody(w, args.ct) / (h*c) # meters * joules/m^2/nm/sec/sr / joule*sec * m/s
	#print(ia[i])

# 2-deg quantal human cone fundamentals -- http://www.cvrl.org/cones.htm
# Negative infinities (log(0)) were inserted for consistent array length.
hcf2deg = np.array([
	[390,  -3.21863,  -3.29076,  -1.96598],
	[391,  -3.13652,  -3.20684,  -1.88705],
	[392,  -3.05545,  -3.12367,  -1.80825],
	[393,  -2.97558,  -3.04148,  -1.72974],
	[394,  -2.89711,  -2.96047,  -1.65172],
	[395,  -2.82023,  -2.88089,  -1.57436],
	[396,  -2.74510,  -2.80295,  -1.49783],
	[397,  -2.67193,  -2.72686,  -1.42232],
	[398,  -2.60089,  -2.65287,  -1.34800],
	[399,  -2.53217,  -2.58117,  -1.27505],
	[400,  -2.46595,  -2.51201,  -1.20366],
	[401,  -2.40234,  -2.44551,  -1.13397],
	[402,  -2.34112,  -2.38148,  -1.06610],
	[403,  -2.28199,  -2.31962,  -1.00014],
	[404,  -2.22466,  -2.25966,  -0.93617],
	[405,  -2.16882,  -2.20131,  -0.87429],
	[406,  -2.11430,  -2.14436,  -0.81458],
	[407,  -2.06146,  -2.08897,  -0.75719],
	[408,  -2.01078,  -2.03539,  -0.70223],
	[409,  -1.96274,  -1.98384,  -0.64986],
	[410,  -1.91782,  -1.93457,  -0.60021],
	[411,  -1.87633,  -1.88775,  -0.55336],
	[412,  -1.83798,  -1.84325,  -0.50921],
	[413,  -1.80230,  -1.80089,  -0.46762],
	[414,  -1.76883,  -1.76047,  -0.42843],
	[415,  -1.73708,  -1.72181,  -0.39150],
	[416,  -1.70674,  -1.68475,  -0.35674],
	[417,  -1.67790,  -1.64931,  -0.32428],
	[418,  -1.65082,  -1.61558,  -0.29433],
	[419,  -1.62575,  -1.58360,  -0.26709],
	[420,  -1.60291,  -1.55345,  -0.24275],
	[421,  -1.58242,  -1.52513,  -0.22138],
	[422,  -1.56383,  -1.49833,  -0.20247],
	[423,  -1.54655,  -1.47270,  -0.18540],
	[424,  -1.52999,  -1.44787,  -0.16953],
	[425,  -1.51356,  -1.42347,  -0.15422],
	[426,  -1.49681,  -1.39922,  -0.13899],
	[427,  -1.47980,  -1.37505,  -0.12401],
	[428,  -1.46270,  -1.35100,  -0.10960],
	[429,  -1.44571,  -1.32708,  -0.09609],
	[430,  -1.42902,  -1.30331,  -0.08379],
	[431,  -1.41277,  -1.27973,  -0.07291],
	[432,  -1.39694,  -1.25647,  -0.06320],
	[433,  -1.38146,  -1.23366,  -0.05429],
	[434,  -1.36626,  -1.21144,  -0.04578],
	[435,  -1.35128,  -1.18995,  -0.03732],
	[436,  -1.33650,  -1.16933,  -0.02868],
	[437,  -1.32214,  -1.14970,  -0.02028],
	[438,  -1.30846,  -1.13119,  -0.01268],
	[439,  -1.29574,  -1.11391,  -0.00645],
	[440,  -1.28423,  -1.09798,  -0.00216],
	[441,  -1.27411,  -1.08346,  -0.00020],
	[442,  -1.26512,  -1.07012,  -0.00022],
	[443,  -1.25691,  -1.05766,  -0.00169],
	[444,  -1.24911,  -1.04578,  -0.00408],
	[445,  -1.24136,  -1.03418,  -0.00685],
	[446,  -1.23339,  -1.02265,  -0.00965],
	[447,  -1.22525,  -1.01126,  -0.01272],
	[448,  -1.21706,  -1.00015,  -0.01649],
	[449,  -1.20895,  -0.98948,  -0.02138],
	[450,  -1.20104,  -0.97939,  -0.02782],
	[451,  -1.19340,  -0.96994,  -0.03607],
	[452,  -1.18581,  -0.96085,  -0.04576],
	[453,  -1.17799,  -0.95173,  -0.05636],
	[454,  -1.16968,  -0.94220,  -0.06735],
	[455,  -1.16057,  -0.93189,  -0.07819],
	[456,  -1.15043,  -0.92047,  -0.08843],
	[457,  -1.13914,  -0.90788,  -0.09796],
	[458,  -1.12659,  -0.89415,  -0.10671],
	[459,  -1.11269,  -0.87926,  -0.11466],
	[460,  -1.09736,  -0.86324,  -0.12174],
	[461,  -1.08057,  -0.84614,  -0.12801],
	[462,  -1.06266,  -0.82826,  -0.13387],
	[463,  -1.04403,  -0.80994,  -0.13983],
	[464,  -1.02509,  -0.79155,  -0.14638],
	[465,  -1.00625,  -0.77342,  -0.15403],
	[466,  -0.98786,  -0.75586,  -0.16321],
	[467,  -0.97000,  -0.73897,  -0.17401],
	[468,  -0.95272,  -0.72280,  -0.18646],
	[469,  -0.93604,  -0.70740,  -0.20059],
	[470,  -0.91999,  -0.69281,  -0.21642],
	[471,  -0.90459,  -0.67906,  -0.23395],
	[472,  -0.88975,  -0.66604,  -0.25306],
	[473,  -0.87536,  -0.65362,  -0.27362],
	[474,  -0.86130,  -0.64167,  -0.29547],
	[475,  -0.84748,  -0.63006,  -0.31849],
	[476,  -0.83380,  -0.61866,  -0.34252],
	[477,  -0.82024,  -0.60745,  -0.36735],
	[478,  -0.80680,  -0.59641,  -0.39278],
	[479,  -0.79348,  -0.58551,  -0.41860],
	[480,  -0.78028,  -0.57474,  -0.44460],
	[481,  -0.76720,  -0.56409,  -0.47063],
	[482,  -0.75426,  -0.55360,  -0.49678],
	[483,  -0.74151,  -0.54330,  -0.52323],
	[484,  -0.72896,  -0.53324,  -0.55014],
	[485,  -0.71664,  -0.52346,  -0.57766],
	[486,  -0.70454,  -0.51395,  -0.60586],
	[487,  -0.69243,  -0.50450,  -0.63447],
	[488,  -0.68007,  -0.49484,  -0.66310],
	[489,  -0.66718,  -0.48470,  -0.69138],
	[490,  -0.65351,  -0.47381,  -0.71893],
	[491,  -0.63886,  -0.46198,  -0.74546],
	[492,  -0.62333,  -0.44931,  -0.77105],
	[493,  -0.60706,  -0.43594,  -0.79586],
	[494,  -0.59023,  -0.42206,  -0.82008],
	[495,  -0.57300,  -0.40782,  -0.84385],
	[496,  -0.55550,  -0.39336,  -0.86738],
	[497,  -0.53779,  -0.37871,  -0.89092],
	[498,  -0.51991,  -0.36389,  -0.91475],
	[499,  -0.50188,  -0.34889,  -0.93916],
	[500,  -0.48375,  -0.33372,  -0.96443],
	[501,  -0.46554,  -0.31841,  -0.99080],
	[502,  -0.44731,  -0.30301,  -1.01838],
	[503,  -0.42910,  -0.28759,  -1.04726],
	[504,  -0.41096,  -0.27221,  -1.07751],
	[505,  -0.39294,  -0.25695,  -1.10920],
	[506,  -0.37509,  -0.24186,  -1.14231],
	[507,  -0.35745,  -0.22701,  -1.17636],
	[508,  -0.34005,  -0.21243,  -1.21077],
	[509,  -0.32292,  -0.19817,  -1.24494],
	[510,  -0.30611,  -0.18429,  -1.27831],
	[511,  -0.28964,  -0.17081,  -1.31046],
	[512,  -0.27356,  -0.15775,  -1.34171],
	[513,  -0.25788,  -0.14509,  -1.37256],
	[514,  -0.24263,  -0.13282,  -1.40350],
	[515,  -0.22786,  -0.12094,  -1.43503],
	[516,  -0.21358,  -0.10947,  -1.46755],
	[517,  -0.19990,  -0.09850,  -1.50099],
	[518,  -0.18689,  -0.08817,  -1.53521],
	[519,  -0.17465,  -0.07861,  -1.57006],
	[520,  -0.16327,  -0.06995,  -1.60537],
	[521,  -0.15280,  -0.06226,  -1.64102],
	[522,  -0.14314,  -0.05546,  -1.67701],
	[523,  -0.13417,  -0.04940,  -1.71337],
	[524,  -0.12575,  -0.04393,  -1.75014],
	[525,  -0.11776,  -0.03890,  -1.78735],
	[526,  -0.11007,  -0.03421,  -1.82501],
	[527,  -0.10272,  -0.02986,  -1.86305],
	[528,  -0.09574,  -0.02586,  -1.90139],
	[529,  -0.08916,  -0.02227,  -1.93996],
	[530,  -0.08304,  -0.01910,  -1.97866],
	[531,  -0.07738,  -0.01637,  -2.01743],
	[532,  -0.07209,  -0.01400,  -2.05636],
	[533,  -0.06704,  -0.01189,  -2.09553],
	[534,  -0.06210,  -0.00994,  -2.13503],
	[535,  -0.05714,  -0.00805,  -2.17496],
	[536,  -0.05208,  -0.00615,  -2.21538],
	[537,  -0.04701,  -0.00433,  -2.25624],
	[538,  -0.04206,  -0.00269,  -2.29745],
	[539,  -0.03735,  -0.00136,  -2.33891],
	[540,  -0.03303,  -0.00043,  -2.38056],
	[541,  -0.02918,  -0.00002,  -2.42230],
	[542,  -0.02583,  -0.00013,  -2.46413],
	[543,  -0.02296,  -0.00077,  -2.50607],
	[544,  -0.02058,  -0.00191,  -2.54812],
	[545,  -0.01866,  -0.00358,  -2.59028],
	[546,  -0.01718,  -0.00572,  -2.63256],
	[547,  -0.01601,  -0.00822,  -2.67498],
	[548,  -0.01498,  -0.01093,  -2.71754],
	[549,  -0.01395,  -0.01368,  -2.76026],
	[550,  -0.01276,  -0.01634,  -2.80313],
	[551,  -0.01130,  -0.01880,  -2.84617],
	[552,  -0.00966,  -0.02118,  -2.88932],
	[553,  -0.00798,  -0.02365,  -2.93253],
	[554,  -0.00640,  -0.02637,  -2.97574],
	[555,  -0.00505,  -0.02953,  -3.01890],
	[556,  -0.00403,  -0.03323,  -3.06197],
	[557,  -0.00330,  -0.03741,  -3.10495],
	[558,  -0.00276,  -0.04192,  -3.14788],
	[559,  -0.00233,  -0.04663,  -3.19077],
	[560,  -0.00193,  -0.05142,  -3.23365],
	[561,  -0.00150,  -0.05619,  -3.27653],
	[562,  -0.00105,  -0.06102,  -3.31942],
	[563,  -0.00063,  -0.06601,  -3.36228],
	[564,  -0.00028,  -0.07126,  -3.40512],
	[565,  -0.00006,  -0.07689,  -3.44790],
	[566,   0.00000,  -0.08298,  -3.49061],
	[567,  -0.00011,  -0.08952,  -3.53326],
	[568,  -0.00039,  -0.09648,  -3.57583],
	[569,  -0.00085,  -0.10380,  -3.61832],
	[570,  -0.00149,  -0.11147,  -3.66073],
	[571,  -0.00232,  -0.11946,  -3.70305],
	[572,  -0.00339,  -0.12784,  -3.74527],
	[573,  -0.00476,  -0.13670,  -3.78740],
	[574,  -0.00649,  -0.14613,  -3.82942],
	[575,  -0.00863,  -0.15622,  -3.87133],
	[576,  -0.01121,  -0.16702,  -3.91312],
	[577,  -0.01408,  -0.17840,  -3.95478],
	[578,  -0.01705,  -0.19018,  -3.99632],
	[579,  -0.01993,  -0.20220,  -4.03772],
	[580,  -0.02252,  -0.21429,  -4.07899],
	[581,  -0.02472,  -0.22632,  -4.12012],
	[582,  -0.02663,  -0.23834,  -4.16109],
	[583,  -0.02843,  -0.25044,  -4.20192],
	[584,  -0.03033,  -0.26272,  -4.24259],
	[585,  -0.03249,  -0.27526,  -4.28311],
	[586,  -0.03507,  -0.28817,  -4.32346],
	[587,  -0.03806,  -0.30148,  -4.36364],
	[588,  -0.04142,  -0.31525,  -4.40365],
	[589,  -0.04510,  -0.32951,  -4.44349],
	[590,  -0.04907,  -0.34432,  -4.48316],
	[591,  -0.05329,  -0.35969,  -4.52264],
	[592,  -0.05775,  -0.37562,  -4.56194],
	[593,  -0.06247,  -0.39207,  -4.60105],
	[594,  -0.06746,  -0.40900,  -4.63998],
	[595,  -0.07271,  -0.42637,  -4.67871],
	[596,  -0.07823,  -0.44417,  -4.71725],
	[597,  -0.08401,  -0.46240,  -4.75559],
	[598,  -0.09001,  -0.48107,  -4.79373],
	[599,  -0.09622,  -0.50020,  -4.83167],
	[600,  -0.10261,  -0.51981,  -4.86941],
	[601,  -0.10916,  -0.53990,  -4.90694],
	[602,  -0.11593,  -0.56045,  -4.94426],
	[603,  -0.12295,  -0.58145,  -4.98137],
	[604,  -0.13029,  -0.60287,  -5.01827],
	[605,  -0.13801,  -0.62468,  -5.05496],
	[606,  -0.14613,  -0.64686,  -5.09143],
	[607,  -0.15465,  -0.66939,  -5.12769],
	[608,  -0.16354,  -0.69227,  -5.16373],
	[609,  -0.17277,  -0.71546,  -5.19955],
	[610,  -0.18231,  -0.73896,  -5.23514],
	[611,  -0.19215,  -0.76275,  -5.27052],
	[612,  -0.20228,  -0.78684,  -5.30568],
	[613,  -0.21273,  -0.81124,  -5.34061],
	[614,  -0.22350,  -0.83597,  -5.37532],
	[615,  -0.23461,  -0.86103,  -5.40980],
	[616,  -0.24607,  -0.88645,  -float('inf')],
	[617,  -0.25783,  -0.91220,  -float('inf')],
	[618,  -0.26982,  -0.93830,  -float('inf')],
	[619,  -0.28201,  -0.96473,  -float('inf')],
	[620,  -0.29432,  -0.99150,  -float('inf')],
	[621,  -0.30674,  -1.01858,  -float('inf')],
	[622,  -0.31939,  -1.04596,  -float('inf')],
	[623,  -0.33243,  -1.07358,  -float('inf')],
	[624,  -0.34601,  -1.10142,  -float('inf')],
	[625,  -0.36030,  -1.12943,  -float('inf')],
	[626,  -0.37542,  -1.15759,  -float('inf')],
	[627,  -0.39127,  -1.18590,  -float('inf')],
	[628,  -0.40776,  -1.21441,  -float('inf')],
	[629,  -0.42475,  -1.24314,  -float('inf')],
	[630,  -0.44212,  -1.27212,  -float('inf')],
	[631,  -0.45978,  -1.30137,  -float('inf')],
	[632,  -0.47769,  -1.33086,  -float('inf')],
	[633,  -0.49582,  -1.36056,  -float('inf')],
	[634,  -0.51416,  -1.39045,  -float('inf')],
	[635,  -0.53269,  -1.42049,  -float('inf')],
	[636,  -0.55139,  -1.45068,  -float('inf')],
	[637,  -0.57023,  -1.48110,  -float('inf')],
	[638,  -0.58918,  -1.51185,  -float('inf')],
	[639,  -0.60820,  -1.54306,  -float('inf')],
	[640,  -0.62726,  -1.57482,  -float('inf')],
	[641,  -0.64637,  -1.60717,  -float('inf')],
	[642,  -0.66566,  -1.63988,  -float('inf')],
	[643,  -0.68529,  -1.67262,  -float('inf')],
	[644,  -0.70542,  -1.70510,  -float('inf')],
	[645,  -0.72621,  -1.73699,  -float('inf')],
	[646,  -0.74780,  -1.76809,  -float('inf')],
	[647,  -0.77013,  -1.79861,  -float('inf')],
	[648,  -0.79311,  -1.82889,  -float('inf')],
	[649,  -0.81665,  -1.85925,  -float('inf')],
	[650,  -0.84067,  -1.89000,  -float('inf')],
	[651,  -0.86509,  -1.92141,  -float('inf')],
	[652,  -0.88985,  -1.95344,  -float('inf')],
	[653,  -0.91491,  -1.98600,  -float('inf')],
	[654,  -0.94025,  -2.01900,  -float('inf')],
	[655,  -0.96582,  -2.05233,  -float('inf')],
	[656,  -0.99158,  -2.08591,  -float('inf')],
	[657,  -1.01754,  -2.11971,  -float('inf')],
	[658,  -1.04369,  -2.15368,  -float('inf')],
	[659,  -1.07006,  -2.18780,  -float('inf')],
	[660,  -1.09663,  -2.22203,  -float('inf')],
	[661,  -1.12342,  -2.25633,  -float('inf')],
	[662,  -1.15043,  -2.29061,  -float('inf')],
	[663,  -1.17764,  -2.32476,  -float('inf')],
	[664,  -1.20506,  -2.35868,  -float('inf')],
	[665,  -1.23268,  -2.39227,  -float('inf')],
	[666,  -1.26050,  -2.42545,  -float('inf')],
	[667,  -1.28853,  -2.45830,  -float('inf')],
	[668,  -1.31676,  -2.49092,  -float('inf')],
	[669,  -1.34521,  -2.52341,  -float('inf')],
	[670,  -1.37388,  -2.55589,  -float('inf')],
	[671,  -1.40278,  -2.58845,  -float('inf')],
	[672,  -1.43191,  -2.62108,  -float('inf')],
	[673,  -1.46127,  -2.65378,  -float('inf')],
	[674,  -1.49089,  -2.68655,  -float('inf')],
	[675,  -1.52075,  -2.71938,  -float('inf')],
	[676,  -1.55087,  -2.75226,  -float('inf')],
	[677,  -1.58124,  -2.78520,  -float('inf')],
	[678,  -1.61183,  -2.81818,  -float('inf')],
	[679,  -1.64264,  -2.85123,  -float('inf')],
	[680,  -1.67365,  -2.88433,  -float('inf')],
	[681,  -1.70485,  -2.91750,  -float('inf')],
	[682,  -1.73631,  -2.95079,  -float('inf')],
	[683,  -1.76808,  -2.98425,  -float('inf')],
	[684,  -1.80024,  -3.01795,  -float('inf')],
	[685,  -1.83284,  -3.05194,  -float('inf')],
	[686,  -1.86590,  -3.08623,  -float('inf')],
	[687,  -1.89927,  -3.12069,  -float('inf')],
	[688,  -1.93274,  -3.15516,  -float('inf')],
	[689,  -1.96611,  -3.18945,  -float('inf')],
	[690,  -1.99917,  -3.22339,  -float('inf')],
	[691,  -2.03177,  -3.25686,  -float('inf')],
	[692,  -2.06399,  -3.28991,  -float('inf')],
	[693,  -2.09595,  -3.32262,  -float('inf')],
	[694,  -2.12779,  -3.35508,  -float('inf')],
	[695,  -2.15963,  -3.38741,  -float('inf')],
	[696,  -2.19157,  -3.41966,  -float('inf')],
	[697,  -2.22362,  -3.45187,  -float('inf')],
	[698,  -2.25573,  -3.48405,  -float('inf')],
	[699,  -2.28788,  -3.51622,  -float('inf')],
	[700,  -2.32004,  -3.54839,  -float('inf')],
	[701,  -2.35220,  -3.58057,  -float('inf')],
	[702,  -2.38440,  -3.61281,  -float('inf')],
	[703,  -2.41671,  -3.64515,  -float('inf')],
	[704,  -2.44920,  -3.67763,  -float('inf')],
	[705,  -2.48194,  -3.71031,  -float('inf')],
	[706,  -2.51497,  -3.74321,  -float('inf')],
	[707,  -2.54826,  -3.77626,  -float('inf')],
	[708,  -2.58174,  -3.80941,  -float('inf')],
	[709,  -2.61534,  -3.84257,  -float('inf')],
	[710,  -2.64902,  -3.87567,  -float('inf')],
	[711,  -2.68271,  -3.90864,  -float('inf')],
	[712,  -2.71635,  -3.94146,  -float('inf')],
	[713,  -2.74989,  -3.97412,  -float('inf')],
	[714,  -2.78329,  -4.00660,  -float('inf')],
	[715,  -2.81649,  -4.03889,  -float('inf')],
	[716,  -2.84945,  -4.07097,  -float('inf')],
	[717,  -2.88223,  -4.10289,  -float('inf')],
	[718,  -2.91489,  -4.13469,  -float('inf')],
	[719,  -2.94750,  -4.16642,  -float('inf')],
	[720,  -2.98012,  -4.19812,  -float('inf')],
	[721,  -3.01280,  -4.22983,  -float('inf')],
	[722,  -3.04551,  -4.26151,  -float('inf')],
	[723,  -3.07819,  -4.29312,  -float('inf')],
	[724,  -3.11077,  -4.32461,  -float('inf')],
	[725,  -3.14321,  -4.35592,  -float('inf')],
	[726,  -3.17547,  -4.38704,  -float('inf')],
	[727,  -3.20755,  -4.41798,  -float('inf')],
	[728,  -3.23951,  -4.44878,  -float('inf')],
	[729,  -3.27137,  -4.47948,  -float('inf')],
	[730,  -3.30318,  -4.51014,  -float('inf')],
	[731,  -3.33497,  -4.54077,  -float('inf')],
	[732,  -3.36677,  -4.57140,  -float('inf')],
	[733,  -3.39861,  -4.60204,  -float('inf')],
	[734,  -3.43052,  -4.63269,  -float('inf')],
	[735,  -3.46253,  -4.66338,  -float('inf')],
	[736,  -3.49464,  -4.69408,  -float('inf')],
	[737,  -3.52678,  -4.72473,  -float('inf')],
	[738,  -3.55883,  -4.75522,  -float('inf')],
	[739,  -3.59069,  -4.78546,  -float('inf')],
	[740,  -3.62227,  -4.81534,  -float('inf')],
	[741,  -3.65348,  -4.84480,  -float('inf')],
	[742,  -3.68442,  -4.87394,  -float('inf')],
	[743,  -3.71519,  -4.90288,  -float('inf')],
	[744,  -3.74590,  -4.93176,  -float('inf')],
	[745,  -3.77667,  -4.96070,  -float('inf')],
	[746,  -3.80758,  -4.98981,  -float('inf')],
	[747,  -3.83860,  -5.01904,  -float('inf')],
	[748,  -3.86965,  -5.04833,  -float('inf')],
	[749,  -3.90068,  -5.07760,  -float('inf')],
	[750,  -3.93161,  -5.10679,  -float('inf')],
	[751,  -3.96239,  -5.13584,  -float('inf')],
	[752,  -3.99302,  -5.16474,  -float('inf')],
	[753,  -4.02351,  -5.19349,  -float('inf')],
	[754,  -4.05387,  -5.22209,  -float('inf')],
	[755,  -4.08412,  -5.25055,  -float('inf')],
	[756,  -4.11426,  -5.27886,  -float('inf')],
	[757,  -4.14433,  -5.30708,  -float('inf')],
	[758,  -4.17435,  -5.33524,  -float('inf')],
	[759,  -4.20436,  -5.36340,  -float('inf')],
	[760,  -4.23438,  -5.39159,  -float('inf')],
	[761,  -4.26443,  -5.41985,  -float('inf')],
	[762,  -4.29447,  -5.44814,  -float('inf')],
	[763,  -4.32445,  -5.47640,  -float('inf')],
	[764,  -4.35431,  -5.50458,  -float('inf')],
	[765,  -4.38399,  -5.53262,  -float('inf')],
	[766,  -4.41348,  -5.56049,  -float('inf')],
	[767,  -4.44282,  -5.58824,  -float('inf')],
	[768,  -4.47210,  -5.61594,  -float('inf')],
	[769,  -4.50141,  -5.64365,  -float('inf')],
	[770,  -4.53081,  -5.67146,  -float('inf')],
	[771,  -4.56037,  -5.69940,  -float('inf')],
	[772,  -4.59000,  -5.72740,  -float('inf')],
	[773,  -4.61960,  -5.75534,  -float('inf')],
	[774,  -4.64907,  -5.78314,  -float('inf')],
	[775,  -4.67830,  -5.81067,  -float('inf')],
	[776,  -4.70721,  -5.83787,  -float('inf')],
	[777,  -4.73587,  -5.86479,  -float('inf')],
	[778,  -4.76433,  -5.89151,  -float('inf')],
	[779,  -4.79271,  -5.91813,  -float('inf')],
	[780,  -4.82106,  -5.94471,  -float('inf')],
	[781,  -4.84947,  -5.97134,  -float('inf')],
	[782,  -4.87791,  -5.99800,  -float('inf')],
	[783,  -4.90636,  -6.02466,  -float('inf')],
	[784,  -4.93478,  -6.05128,  -float('inf')],
	[785,  -4.96314,  -6.07785,  -float('inf')],
	[786,  -4.99142,  -6.10433,  -float('inf')],
	[787,  -5.01964,  -6.13074,  -float('inf')],
	[788,  -5.04780,  -6.15710,  -float('inf')],
	[789,  -5.07592,  -6.18343,  -float('inf')],
	[790,  -5.10404,  -6.20974,  -float('inf')],
	[791,  -5.13215,  -6.23606,  -float('inf')],
	[792,  -5.16024,  -6.26233,  -float('inf')],
	[793,  -5.18829,  -6.28853,  -float('inf')],
	[794,  -5.21625,  -6.31461,  -float('inf')],
	[795,  -5.24411,  -6.34051,  -float('inf')],
	[796,  -5.27185,  -6.36622,  -float('inf')],
	[797,  -5.29945,  -6.39173,  -float('inf')],
	[798,  -5.32693,  -6.41707,  -float('inf')],
	[799,  -5.35428,  -6.44226,  -float('inf')],
	[800,  -5.38150,  -6.46732,  -float('inf')],
	[801,  -5.40862,  -6.49227,  -float('inf')],
	[802,  -5.43566,  -6.51716,  -float('inf')],
	[803,  -5.46268,  -6.54203,  -float('inf')],
	[804,  -5.48972,  -6.56693,  -float('inf')],
	[805,  -5.51683,  -6.59192,  -float('inf')],
	[806,  -5.54404,  -6.61701,  -float('inf')],
	[807,  -5.57132,  -6.64218,  -float('inf')],
	[808,  -5.59861,  -6.66736,  -float('inf')],
	[809,  -5.62586,  -6.69251,  -float('inf')],
	[810,  -5.65303,  -6.71758,  -float('inf')],
	[811,  -5.68006,  -6.74251,  -float('inf')],
	[812,  -5.70694,  -6.76731,  -float('inf')],
	[813,  -5.73368,  -6.79197,  -float('inf')],
	[814,  -5.76026,  -6.81650,  -float('inf')],
	[815,  -5.78669,  -6.84090,  -float('inf')],
	[816,  -5.81296,  -6.86518,  -float('inf')],
	[817,  -5.83911,  -6.88937,  -float('inf')],
	[818,  -5.86518,  -6.91353,  -float('inf')],
	[819,  -5.89120,  -6.93768,  -float('inf')],
	[820,  -5.91722,  -6.96187,  -float('inf')],
	[821,  -5.94325,  -6.98614,  -float('inf')],
	[822,  -5.96929,  -7.01048,  -float('inf')],
	[823,  -5.99533,  -7.03487,  -float('inf')],
	[824,  -6.02135,  -7.05929,  -float('inf')],
	[825,  -6.04732,  -7.08374,  -float('inf')],
	[826,  -6.07324,  -7.10819,  -float('inf')],
	[827,  -6.09908,  -7.13264,  -float('inf')],
	[828,  -6.12483,  -7.15706,  -float('inf')],
	[829,  -6.15048,  -7.18143,  -float('inf')],
	[830,  -6.17600,  -7.20576,  -float('inf')]
])

# 10-deg
hcf10deg = np.array([
	[390,  -3.22745,  -3.30405,  -2.15199],
	[391,  -3.14240,  -3.21720,  -2.06972],
	[392,  -3.05838,  -3.13108,  -1.98752],
	[393,  -2.97555,  -3.04592,  -1.90557],
	[394,  -2.89410,  -2.96192,  -1.82403],
	[395,  -2.81418,  -2.87930,  -1.74308],
	[396,  -2.73598,  -2.79825,  -1.66289],
	[397,  -2.65966,  -2.71900,  -1.58364],
	[398,  -2.58540,  -2.64174,  -1.50551],
	[399,  -2.51337,  -2.56670,  -1.42867],
	[400,  -2.44375,  -2.49408,  -1.35331],
	[401,  -2.37662,  -2.42401,  -1.27959],
	[402,  -2.31180,  -2.35633,  -1.20764],
	[403,  -2.24903,  -2.29078,  -1.13759],
	[404,  -2.18804,  -2.22711,  -1.06955],
	[405,  -2.12856,  -2.16507,  -1.00365],
	[406,  -2.07045,  -2.10449,  -0.93999],
	[407,  -2.01407,  -2.04551,  -0.87871],
	[408,  -1.95988,  -1.98836,  -0.81992],
	[409,  -1.90834,  -1.93326,  -0.76376],
	[410,  -1.85993,  -1.88044,  -0.71035],
	[411,  -1.81495,  -1.83002,  -0.65976],
	[412,  -1.77303,  -1.78185,  -0.61182],
	[413,  -1.73362,  -1.73565,  -0.56631],
	[414,  -1.69620,  -1.69118,  -0.52301],
	[415,  -1.66022,  -1.64817,  -0.48169],
	[416,  -1.62529,  -1.60642,  -0.44222],
	[417,  -1.59158,  -1.56601,  -0.40477],
	[418,  -1.55938,  -1.52705,  -0.36962],
	[419,  -1.52901,  -1.48966,  -0.33702],
	[420,  -1.50076,  -1.45398,  -0.30724],
	[421,  -1.47479,  -1.42005,  -0.28040],
	[422,  -1.45075,  -1.38766,  -0.25609],
	[423,  -1.42809,  -1.35651,  -0.23373],
	[424,  -1.40633,  -1.32633,  -0.21277],
	[425,  -1.38493,  -1.29683,  -0.19265],
	[426,  -1.36350,  -1.26776,  -0.17294],
	[427,  -1.34214,  -1.23911,  -0.15386],
	[428,  -1.32104,  -1.21092,  -0.13573],
	[429,  -1.30040,  -1.18321,  -0.11890],
	[430,  -1.28044,  -1.15602,  -0.10372],
	[431,  -1.26129,  -1.12940,  -0.09039],
	[432,  -1.24291,  -1.10342,  -0.07861],
	[433,  -1.22519,  -1.07819,  -0.06795],
	[434,  -1.20800,  -1.05382,  -0.05797],
	[435,  -1.19125,  -1.03040,  -0.04826],
	[436,  -1.17488,  -1.00804,  -0.03852],
	[437,  -1.15899,  -0.98674,  -0.02908],
	[438,  -1.14375,  -0.96654,  -0.02042],
	[439,  -1.12931,  -0.94742,  -0.01301],
	[440,  -1.11586,  -0.92941,  -0.00732],
	[441,  -1.10345,  -0.91247,  -0.00367],
	[442,  -1.09184,  -0.89636,  -0.00169],
	[443,  -1.08068,  -0.88080,  -0.00086],
	[444,  -1.06962,  -0.86551,  -0.00065],
	[445,  -1.05834,  -0.85023,  -0.00054],
	[446,  -1.04658,  -0.83475,  -0.00019],
	[447,  -1.03452,  -0.81929,  -0.00001],
	[448,  -1.02241,  -0.80411,  -0.00057],
	[449,  -1.01053,  -0.78952,  -0.00246],
	[450,  -0.99914,  -0.77581,  -0.00627],
	[451,  -0.98840,  -0.76313,  -0.01238],
	[452,  -0.97811,  -0.75120,  -0.02038],
	[453,  -0.96794,  -0.73959,  -0.02968],
	[454,  -0.95758,  -0.72787,  -0.03969],
	[455,  -0.94671,  -0.71563,  -0.04981],
	[456,  -0.93508,  -0.70254,  -0.05958],
	[457,  -0.92267,  -0.68866,  -0.06900],
	[458,  -0.90955,  -0.67415,  -0.07818],
	[459,  -0.89577,  -0.65918,  -0.08723],
	[460,  -0.88139,  -0.64390,  -0.09626],
	[461,  -0.86648,  -0.62847,  -0.10540],
	[462,  -0.85124,  -0.61307,  -0.11496],
	[463,  -0.83589,  -0.59785,  -0.12524],
	[464,  -0.82065,  -0.58297,  -0.13656],
	[465,  -0.80575,  -0.56860,  -0.14924],
	[466,  -0.79135,  -0.55486,  -0.16352],
	[467,  -0.77744,  -0.54173,  -0.17940],
	[468,  -0.76394,  -0.52918,  -0.19678],
	[469,  -0.75078,  -0.51714,  -0.21558],
	[470,  -0.73790,  -0.50557,  -0.23572],
	[471,  -0.72522,  -0.49442,  -0.25712],
	[472,  -0.71267,  -0.48358,  -0.27966],
	[473,  -0.70018,  -0.47296,  -0.30324],
	[474,  -0.68768,  -0.46246,  -0.32773],
	[475,  -0.67508,  -0.45196,  -0.35305],
	[476,  -0.66233,  -0.44139,  -0.37906],
	[477,  -0.64942,  -0.43071,  -0.40559],
	[478,  -0.63638,  -0.41994,  -0.43244],
	[479,  -0.62322,  -0.40906,  -0.45942],
	[480,  -0.60995,  -0.39808,  -0.48632],
	[481,  -0.59660,  -0.38702,  -0.51304],
	[482,  -0.58330,  -0.37600,  -0.53976],
	[483,  -0.57016,  -0.36517,  -0.56674],
	[484,  -0.55732,  -0.35469,  -0.59426],
	[485,  -0.54491,  -0.34470,  -0.62257],
	[486,  -0.53301,  -0.33528,  -0.65182],
	[487,  -0.52139,  -0.32622,  -0.68173],
	[488,  -0.50977,  -0.31723,  -0.71193],
	[489,  -0.49788,  -0.30800,  -0.74201],
	[490,  -0.48543,  -0.29825,  -0.77161],
	[491,  -0.47224,  -0.28777,  -0.80043],
	[492,  -0.45847,  -0.27674,  -0.82866],
	[493,  -0.44440,  -0.26545,  -0.85655],
	[494,  -0.43031,  -0.25420,  -0.88438],
	[495,  -0.41648,  -0.24328,  -0.91242],
	[496,  -0.40312,  -0.23290,  -0.94092],
	[497,  -0.39018,  -0.22299,  -0.97003],
	[498,  -0.37753,  -0.21340,  -0.99988],
	[499,  -0.36507,  -0.20396,  -1.03062],
	[500,  -0.35265,  -0.19453,  -1.06238],
	[501,  -0.34020,  -0.18498,  -1.09528],
	[502,  -0.32773,  -0.17535,  -1.12940],
	[503,  -0.31530,  -0.16571,  -1.16483],
	[504,  -0.30295,  -0.15613,  -1.20165],
	[505,  -0.29073,  -0.14670,  -1.23992],
	[506,  -0.27870,  -0.13748,  -1.27961],
	[507,  -0.26687,  -0.12851,  -1.32023],
	[508,  -0.25527,  -0.11981,  -1.36117],
	[509,  -0.24390,  -0.11142,  -1.40183],
	[510,  -0.23278,  -0.10337,  -1.44160],
	[511,  -0.22193,  -0.09567,  -1.48005],
	[512,  -0.21130,  -0.08826,  -1.51744],
	[513,  -0.20087,  -0.08104,  -1.55421],
	[514,  -0.19058,  -0.07393,  -1.59080],
	[515,  -0.18040,  -0.06684,  -1.62764],
	[516,  -0.17030,  -0.05972,  -1.66507],
	[517,  -0.16040,  -0.05271,  -1.70307],
	[518,  -0.15082,  -0.04599,  -1.74148],
	[519,  -0.14169,  -0.03974,  -1.78018],
	[520,  -0.13314,  -0.03415,  -1.81903],
	[521,  -0.12528,  -0.02935,  -1.85794],
	[522,  -0.11805,  -0.02530,  -1.89695],
	[523,  -0.11136,  -0.02188,  -1.93614],
	[524,  -0.10513,  -0.01898,  -1.97559],
	[525,  -0.09927,  -0.01651,  -2.01539],
	[526,  -0.09370,  -0.01435,  -2.05559],
	[527,  -0.08839,  -0.01248,  -2.09608],
	[528,  -0.08333,  -0.01086,  -2.13673],
	[529,  -0.07851,  -0.00947,  -2.17741],
	[530,  -0.07390,  -0.00827,  -2.21798],
	[531,  -0.06948,  -0.00723,  -2.25837],
	[532,  -0.06517,  -0.00630,  -2.29865],
	[533,  -0.06089,  -0.00543,  -2.33898],
	[534,  -0.05655,  -0.00456,  -2.37948],
	[535,  -0.05206,  -0.00362,  -2.42030],
	[536,  -0.04737,  -0.00260,  -2.46154],
	[537,  -0.04259,  -0.00160,  -2.50315],
	[538,  -0.03787,  -0.00074,  -2.54504],
	[539,  -0.03337,  -0.00017,  -2.58714],
	[540,  -0.02922,  -0.00001,  -2.62935],
	[541,  -0.02557,  -0.00037,  -2.67161],
	[542,  -0.02244,  -0.00128,  -2.71391],
	[543,  -0.01982,  -0.00274,  -2.75627],
	[544,  -0.01771,  -0.00475,  -2.79870],
	[545,  -0.01612,  -0.00730,  -2.84122],
	[546,  -0.01501,  -0.01038,  -2.88383],
	[547,  -0.01421,  -0.01381,  -2.92655],
	[548,  -0.01354,  -0.01741,  -2.96936],
	[549,  -0.01280,  -0.02100,  -3.01228],
	[550,  -0.01179,  -0.02439,  -3.05529],
	[551,  -0.01040,  -0.02747,  -3.09840],
	[552,  -0.00873,  -0.03038,  -3.14157],
	[553,  -0.00696,  -0.03334,  -3.18476],
	[554,  -0.00527,  -0.03656,  -3.22794],
	[555,  -0.00385,  -0.04025,  -3.27107],
	[556,  -0.00282,  -0.04455,  -3.31412],
	[557,  -0.00214,  -0.04937,  -3.35710],
	[558,  -0.00169,  -0.05457,  -3.40003],
	[559,  -0.00138,  -0.05999,  -3.44293],
	[560,  -0.00112,  -0.06548,  -3.48582],
	[561,  -0.00081,  -0.07095,  -3.52871],
	[562,  -0.00050,  -0.07645,  -3.57160],
	[563,  -0.00023,  -0.08212,  -3.61446],
	[564,  -0.00005,  -0.08807,  -3.65729],
	[565,  -0.00001,  -0.09442,  -3.70007],
	[566,  -0.00015,  -0.10125,  -3.74278],
	[567,  -0.00049,  -0.10856,  -3.78543],
	[568,  -0.00103,  -0.11629,  -3.82800],
	[569,  -0.00177,  -0.12441,  -3.87049],
	[570,  -0.00272,  -0.13287,  -3.91290],
	[571,  -0.00388,  -0.14165,  -3.95522],
	[572,  -0.00532,  -0.15084,  -3.99745],
	[573,  -0.00711,  -0.16051,  -4.03957],
	[574,  -0.00930,  -0.17077,  -4.08159],
	[575,  -0.01197,  -0.18171,  -4.12350],
	[576,  -0.01513,  -0.19338,  -4.16529],
	[577,  -0.01861,  -0.20562,  -4.20696],
	[578,  -0.02219,  -0.21827,  -4.24849],
	[579,  -0.02565,  -0.23113,  -4.28990],
	[580,  -0.02878,  -0.24402,  -4.33116],
	[581,  -0.03142,  -0.25681,  -4.37229],
	[582,  -0.03374,  -0.26956,  -4.41327],
	[583,  -0.03592,  -0.28236,  -4.45409],
	[584,  -0.03820,  -0.29532,  -4.49477],
	[585,  -0.04077,  -0.30854,  -4.53528],
	[586,  -0.04381,  -0.32210,  -4.57563],
	[587,  -0.04731,  -0.33606,  -4.61581],
	[588,  -0.05121,  -0.35047,  -4.65583],
	[589,  -0.05547,  -0.36537,  -4.69567],
	[590,  -0.06004,  -0.38081,  -4.73533],
	[591,  -0.06487,  -0.39681,  -4.77481],
	[592,  -0.06997,  -0.41336,  -4.81411],
	[593,  -0.07533,  -0.43041,  -4.85323],
	[594,  -0.08096,  -0.44792,  -4.89215],
	[595,  -0.08686,  -0.46587,  -4.93088],
	[596,  -0.09303,  -0.48422,  -4.96942],
	[597,  -0.09944,  -0.50297,  -5.00776],
	[598,  -0.10609,  -0.52216,  -5.04590],
	[599,  -0.11295,  -0.54179,  -5.08384],
	[600,  -0.11999,  -0.56189,  -5.12158],
	[601,  -0.12720,  -0.58245,  -5.15911],
	[602,  -0.13464,  -0.60347,  -5.19643],
	[603,  -0.14234,  -0.62492,  -5.23354],
	[604,  -0.15037,  -0.64676,  -5.27045],
	[605,  -0.15877,  -0.66899,  -5.30713],
	[606,  -0.16758,  -0.69156,  -5.34361],
	[607,  -0.17680,  -0.71447,  -5.37986],
	[608,  -0.18639,  -0.73770,  -5.41590],
	[609,  -0.19631,  -0.76124,  -5.45172],
	[610,  -0.20655,  -0.78506,  -5.48732],
	[611,  -0.21708,  -0.80917,  -5.52270],
	[612,  -0.22791,  -0.83355,  -5.55785],
	[613,  -0.23904,  -0.85824,  -5.59278],
	[614,  -0.25050,  -0.88324,  -5.62749],
	[615,  -0.26228,  -0.90856,  -5.66197],
	[616,  -0.27440,  -0.93422,  -float('inf')],
	[617,  -0.28681,  -0.96020,  -float('inf')],
	[618,  -0.29945,  -0.98652,  -float('inf')],
	[619,  -0.31225,  -1.01316,  -float('inf')],
	[620,  -0.32517,  -1.04012,  -float('inf')],
	[621,  -0.33817,  -1.06739,  -float('inf')],
	[622,  -0.35140,  -1.09495,  -float('inf')],
	[623,  -0.36501,  -1.12274,  -float('inf')],
	[624,  -0.37917,  -1.15073,  -float('inf')],
	[625,  -0.39404,  -1.17889,  -float('inf')],
	[626,  -0.40974,  -1.20719,  -float('inf')],
	[627,  -0.42619,  -1.23564,  -float('inf')],
	[628,  -0.44326,  -1.26427,  -float('inf')],
	[629,  -0.46083,  -1.29312,  -float('inf')],
	[630,  -0.47876,  -1.32220,  -float('inf')],
	[631,  -0.49696,  -1.35155,  -float('inf')],
	[632,  -0.51538,  -1.38114,  -float('inf')],
	[633,  -0.53402,  -1.41093,  -float('inf')],
	[634,  -0.55284,  -1.44091,  -float('inf')],
	[635,  -0.57184,  -1.47103,  -float('inf')],
	[636,  -0.59098,  -1.50129,  -float('inf')],
	[637,  -0.61025,  -1.53178,  -float('inf')],
	[638,  -0.62960,  -1.56260,  -float('inf')],
	[639,  -0.64901,  -1.59387,  -float('inf')],
	[640,  -0.66844,  -1.62568,  -float('inf')],
	[641,  -0.68790,  -1.65809,  -float('inf')],
	[642,  -0.70752,  -1.69085,  -float('inf')],
	[643,  -0.72746,  -1.72364,  -float('inf')],
	[644,  -0.74790,  -1.75616,  -float('inf')],
	[645,  -0.76901,  -1.78809,  -float('inf')],
	[646,  -0.79089,  -1.81922,  -float('inf')],
	[647,  -0.81351,  -1.84978,  -float('inf')],
	[648,  -0.83678,  -1.88009,  -float('inf')],
	[649,  -0.86060,  -1.91047,  -float('inf')],
	[650,  -0.88489,  -1.94126,  -float('inf')],
	[651,  -0.90956,  -1.97269,  -float('inf')],
	[652,  -0.93457,  -2.00474,  -float('inf')],
	[653,  -0.95987,  -2.03733,  -float('inf')],
	[654,  -0.98543,  -2.07035,  -float('inf')],
	[655,  -1.01120,  -2.10370,  -float('inf')],
	[656,  -1.03716,  -2.13730,  -float('inf')],
	[657,  -1.06331,  -2.17111,  -float('inf')],
	[658,  -1.08964,  -2.20510,  -float('inf')],
	[659,  -1.11617,  -2.23923,  -float('inf')],
	[660,  -1.14290,  -2.27348,  -float('inf')],
	[661,  -1.16984,  -2.30779,  -float('inf')],
	[662,  -1.19699,  -2.34208,  -float('inf')],
	[663,  -1.22433,  -2.37624,  -float('inf')],
	[664,  -1.25188,  -2.41017,  -float('inf')],
	[665,  -1.27961,  -2.44377,  -float('inf')],
	[666,  -1.30755,  -2.47696,  -float('inf')],
	[667,  -1.33568,  -2.50981,  -float('inf')],
	[668,  -1.36401,  -2.54244,  -float('inf')],
	[669,  -1.39255,  -2.57495,  -float('inf')],
	[670,  -1.42131,  -2.60743,  -float('inf')],
	[671,  -1.45029,  -2.63999,  -float('inf')],
	[672,  -1.47950,  -2.67263,  -float('inf')],
	[673,  -1.50894,  -2.70534,  -float('inf')],
	[674,  -1.53862,  -2.73811,  -float('inf')],
	[675,  -1.56855,  -2.77095,  -float('inf')],
	[676,  -1.59873,  -2.80383,  -float('inf')],
	[677,  -1.62915,  -2.83677,  -float('inf')],
	[678,  -1.65980,  -2.86976,  -float('inf')],
	[679,  -1.69066,  -2.90280,  -float('inf')],
	[680,  -1.72172,  -2.93591,  -float('inf')],
	[681,  -1.75297,  -2.96908,  -float('inf')],
	[682,  -1.78447,  -3.00237,  -float('inf')],
	[683,  -1.81628,  -3.03584,  -float('inf')],
	[684,  -1.84847,  -3.06954,  -float('inf')],
	[685,  -1.88110,  -3.10352,  -float('inf')],
	[686,  -1.91420,  -3.13782,  -float('inf')],
	[687,  -1.94760,  -3.17229,  -float('inf')],
	[688,  -1.98110,  -3.20675,  -float('inf')],
	[689,  -2.01449,  -3.24104,  -float('inf')],
	[690,  -2.04757,  -3.27499,  -float('inf')],
	[691,  -2.08020,  -3.30846,  -float('inf')],
	[692,  -2.11244,  -3.34151,  -float('inf')],
	[693,  -2.14442,  -3.37422,  -float('inf')],
	[694,  -2.17627,  -3.40669,  -float('inf')],
	[695,  -2.20813,  -3.43901,  -float('inf')],
	[696,  -2.24009,  -3.47126,  -float('inf')],
	[697,  -2.27214,  -3.50347,  -float('inf')],
	[698,  -2.30427,  -3.53566,  -float('inf')],
	[699,  -2.33644,  -3.56783,  -float('inf')],
	[700,  -2.36861,  -3.59999,  -float('inf')],
	[701,  -2.40077,  -3.63217,  -float('inf')],
	[702,  -2.43298,  -3.66441,  -float('inf')],
	[703,  -2.46530,  -3.69675,  -float('inf')],
	[704,  -2.49780,  -3.72924,  -float('inf')],
	[705,  -2.53055,  -3.76192,  -float('inf')],
	[706,  -2.56359,  -3.79481,  -float('inf')],
	[707,  -2.59688,  -3.82787,  -float('inf')],
	[708,  -2.63037,  -3.86102,  -float('inf')],
	[709,  -2.66398,  -3.89418,  -float('inf')],
	[710,  -2.69766,  -3.92728,  -float('inf')],
	[711,  -2.73135,  -3.96025,  -float('inf')],
	[712,  -2.76500,  -3.99308,  -float('inf')],
	[713,  -2.79855,  -4.02574,  -float('inf')],
	[714,  -2.83195,  -4.05822,  -float('inf')],
	[715,  -2.86515,  -4.09050,  -float('inf')],
	[716,  -2.89812,  -4.12258,  -float('inf')],
	[717,  -2.93090,  -4.15450,  -float('inf')],
	[718,  -2.96356,  -4.18630,  -float('inf')],
	[719,  -2.99617,  -4.21803,  -float('inf')],
	[720,  -3.02880,  -4.24974,  -float('inf')],
	[721,  -3.06148,  -4.28145,  -float('inf')],
	[722,  -3.09419,  -4.31313,  -float('inf')],
	[723,  -3.12687,  -4.34473,  -float('inf')],
	[724,  -3.15946,  -4.37622,  -float('inf')],
	[725,  -3.19190,  -4.40754,  -float('inf')],
	[726,  -3.22415,  -4.43865,  -float('inf')],
	[727,  -3.25624,  -4.46959,  -float('inf')],
	[728,  -3.28820,  -4.50039,  -float('inf')],
	[729,  -3.32006,  -4.53110,  -float('inf')],
	[730,  -3.35187,  -4.56175,  -float('inf')],
	[731,  -3.38367,  -4.59238,  -float('inf')],
	[732,  -3.41547,  -4.62301,  -float('inf')],
	[733,  -3.44731,  -4.65365,  -float('inf')],
	[734,  -3.47922,  -4.68430,  -float('inf')],
	[735,  -3.51123,  -4.71499,  -float('inf')],
	[736,  -3.54334,  -4.74570,  -float('inf')],
	[737,  -3.57547,  -4.77634,  -float('inf')],
	[738,  -3.60753,  -4.80684,  -float('inf')],
	[739,  -3.63939,  -4.83707,  -float('inf')],
	[740,  -3.67097,  -4.86695,  -float('inf')],
	[741,  -3.70219,  -4.89641,  -float('inf')],
	[742,  -3.73312,  -4.92555,  -float('inf')],
	[743,  -3.76389,  -4.95449,  -float('inf')],
	[744,  -3.79460,  -4.98337,  -float('inf')],
	[745,  -3.82538,  -5.01231,  -float('inf')],
	[746,  -3.85629,  -5.04142,  -float('inf')],
	[747,  -3.88730,  -5.07065,  -float('inf')],
	[748,  -3.91836,  -5.09994,  -float('inf')],
	[749,  -3.94938,  -5.12921,  -float('inf')],
	[750,  -3.98031,  -5.15841,  -float('inf')],
	[751,  -4.01109,  -5.18746,  -float('inf')],
	[752,  -4.04172,  -5.21636,  -float('inf')],
	[753,  -4.07221,  -5.24511,  -float('inf')],
	[754,  -4.10258,  -5.27371,  -float('inf')],
	[755,  -4.13282,  -5.30216,  -float('inf')],
	[756,  -4.16297,  -5.33048,  -float('inf')],
	[757,  -4.19304,  -5.35869,  -float('inf')],
	[758,  -4.22306,  -5.38686,  -float('inf')],
	[759,  -4.25307,  -5.41501,  -float('inf')],
	[760,  -4.28309,  -5.44320,  -float('inf')],
	[761,  -4.31314,  -5.47146,  -float('inf')],
	[762,  -4.34318,  -5.49975,  -float('inf')],
	[763,  -4.37315,  -5.52801,  -float('inf')],
	[764,  -4.40301,  -5.55619,  -float('inf')],
	[765,  -4.43270,  -5.58423,  -float('inf')],
	[766,  -4.46219,  -5.61211,  -float('inf')],
	[767,  -4.49153,  -5.63985,  -float('inf')],
	[768,  -4.52081,  -5.66755,  -float('inf')],
	[769,  -4.55011,  -5.69527,  -float('inf')],
	[770,  -4.57952,  -5.72308,  -float('inf')],
	[771,  -4.60907,  -5.75102,  -float('inf')],
	[772,  -4.63871,  -5.77901,  -float('inf')],
	[773,  -4.66831,  -5.80696,  -float('inf')],
	[774,  -4.69778,  -5.83475,  -float('inf')],
	[775,  -4.72701,  -5.86228,  -float('inf')],
	[776,  -4.75592,  -5.88948,  -float('inf')],
	[777,  -4.78457,  -5.91640,  -float('inf')],
	[778,  -4.81304,  -5.94313,  -float('inf')],
	[779,  -4.84141,  -5.96974,  -float('inf')],
	[780,  -4.86977,  -5.99633,  -float('inf')],
	[781,  -4.89818,  -6.02296,  -float('inf')],
	[782,  -4.92662,  -6.04961,  -float('inf')],
	[783,  -4.95507,  -6.07627,  -float('inf')],
	[784,  -4.98349,  -6.10290,  -float('inf')],
	[785,  -5.01185,  -6.12946,  -float('inf')],
	[786,  -5.04013,  -6.15594,  -float('inf')],
	[787,  -5.06834,  -6.18236,  -float('inf')],
	[788,  -5.09650,  -6.20871,  -float('inf')],
	[789,  -5.12463,  -6.23504,  -float('inf')],
	[790,  -5.15275,  -6.26136,  -float('inf')],
	[791,  -5.18086,  -6.28767,  -float('inf')],
	[792,  -5.20895,  -6.31395,  -float('inf')],
	[793,  -5.23700,  -6.34015,  -float('inf')],
	[794,  -5.26496,  -6.36622,  -float('inf')],
	[795,  -5.29282,  -6.39213,  -float('inf')],
	[796,  -5.32055,  -6.41783,  -float('inf')],
	[797,  -5.34816,  -6.44334,  -float('inf')],
	[798,  -5.37563,  -6.46869,  -float('inf')],
	[799,  -5.40298,  -6.49387,  -float('inf')],
	[800,  -5.43021,  -6.51893,  -float('inf')],
	[801,  -5.45733,  -6.54388,  -float('inf')],
	[802,  -5.48437,  -6.56877,  -float('inf')],
	[803,  -5.51139,  -6.59364,  -float('inf')],
	[804,  -5.53843,  -6.61855,  -float('inf')],
	[805,  -5.56554,  -6.64353,  -float('inf')],
	[806,  -5.59275,  -6.66862,  -float('inf')],
	[807,  -5.62003,  -6.69379,  -float('inf')],
	[808,  -5.64732,  -6.71898,  -float('inf')],
	[809,  -5.67457,  -6.74413,  -float('inf')],
	[810,  -5.70174,  -6.76919,  -float('inf')],
	[811,  -5.72876,  -6.79413,  -float('inf')],
	[812,  -5.75565,  -6.81892,  -float('inf')],
	[813,  -5.78238,  -6.84358,  -float('inf')],
	[814,  -5.80897,  -6.86811,  -float('inf')],
	[815,  -5.83539,  -6.89251,  -float('inf')],
	[816,  -5.86167,  -6.91679,  -float('inf')],
	[817,  -5.88782,  -6.94099,  -float('inf')],
	[818,  -5.91389,  -6.96514,  -float('inf')],
	[819,  -5.93991,  -6.98929,  -float('inf')],
	[820,  -5.96593,  -7.01349,  -float('inf')],
	[821,  -5.99196,  -7.03775,  -float('inf')],
	[822,  -6.01800,  -7.06209,  -float('inf')],
	[823,  -6.04404,  -7.08648,  -float('inf')],
	[824,  -6.07005,  -7.11091,  -float('inf')],
	[825,  -6.09603,  -7.13535,  -float('inf')],
	[826,  -6.12195,  -7.15981,  -float('inf')],
	[827,  -6.14779,  -7.18425,  -float('inf')],
	[828,  -6.17354,  -7.20867,  -float('inf')],
	[829,  -6.19919,  -7.23305,  -float('inf')],
	[830,  -6.22470,  -7.25737,  -float('inf')]
])

# 2-deg energy
hcfe2deg = np.array([
	[390,  -3.38195,  -3.43374,  -2.02012],
	[391,  -3.29873,  -3.34871,  -1.94008],
	[392,  -3.21655,  -3.26443,  -1.86017],
	[393,  -3.13558,  -3.18113,  -1.78056],
	[394,  -3.05601,  -3.09902,  -1.70143],
	[395,  -2.97802,  -3.01834,  -1.62297],
	[396,  -2.90179,  -2.93929,  -1.54534],
	[397,  -2.82753,  -2.86212,  -1.46874],
	[398,  -2.75539,  -2.78703,  -1.39332],
	[399,  -2.68558,  -2.71424,  -1.31929],
	[400,  -2.61828,  -2.64399,  -1.24680],
	[401,  -2.55358,  -2.57641,  -1.17603],
	[402,  -2.49128,  -2.51129,  -1.10708],
	[403,  -2.43107,  -2.44836,  -1.04004],
	[404,  -2.37266,  -2.38732,  -0.97500],
	[405,  -2.31575,  -2.32789,  -0.91204],
	[406,  -2.26016,  -2.26988,  -0.85126],
	[407,  -2.20626,  -2.21342,  -0.79280],
	[408,  -2.15451,  -2.15877,  -0.73678],
	[409,  -2.10541,  -2.10616,  -0.68334],
	[410,  -2.05942,  -2.05583,  -0.63263],
	[411,  -2.01688,  -2.00795,  -0.58472],
	[412,  -1.97747,  -1.96240,  -0.53952],
	[413,  -1.94074,  -1.91898,  -0.49687],
	[414,  -1.90621,  -1.87751,  -0.45664],
	[415,  -1.87342,  -1.83780,  -0.41866],
	[416,  -1.84203,  -1.79970,  -0.38285],
	[417,  -1.81215,  -1.76322,  -0.34935],
	[418,  -1.78403,  -1.72844,  -0.31836],
	[419,  -1.75792,  -1.69543,  -0.29008],
	[420,  -1.73405,  -1.66424,  -0.26471],
	[421,  -1.71253,  -1.63489,  -0.24230],
	[422,  -1.69291,  -1.60706,  -0.22237],
	[423,  -1.67460,  -1.58040,  -0.20427],
	[424,  -1.65701,  -1.55454,  -0.18737],
	[425,  -1.63956,  -1.52913,  -0.17103],
	[426,  -1.62179,  -1.50385,  -0.15479],
	[427,  -1.60375,  -1.47866,  -0.13879],
	[428,  -1.58564,  -1.45359,  -0.12336],
	[429,  -1.56764,  -1.42866,  -0.10883],
	[430,  -1.54994,  -1.40388,  -0.09553],
	[431,  -1.53268,  -1.37930,  -0.08364],
	[432,  -1.51585,  -1.35503,  -0.07293],
	[433,  -1.49936,  -1.33121,  -0.06300],
	[434,  -1.48316,  -1.30799,  -0.05350],
	[435,  -1.46718,  -1.28550,  -0.04404],
	[436,  -1.45140,  -1.26389,  -0.03440],
	[437,  -1.43604,  -1.24326,  -0.02500],
	[438,  -1.42137,  -1.22376,  -0.01641],
	[439,  -1.40766,  -1.20549,  -0.00919],
	[440,  -1.39517,  -1.18857,  -0.00392],
	[441,  -1.38406,  -1.17306,  -0.00097],
	[442,  -1.37409,  -1.15874,  -0.00001],
	[443,  -1.36489,  -1.14529,  -0.00049],
	[444,  -1.35611,  -1.13244,  -0.00190],
	[445,  -1.34739,  -1.11987,  -0.00370],
	[446,  -1.33844,  -1.10736,  -0.00552],
	[447,  -1.32933,  -1.09499,  -0.00761],
	[448,  -1.32017,  -1.08292,  -0.01042],
	[449,  -1.31109,  -1.07127,  -0.01434],
	[450,  -1.30221,  -1.06022,  -0.01982],
	[451,  -1.29361,  -1.04981,  -0.02710],
	[452,  -1.28506,  -1.03976,  -0.03583],
	[453,  -1.27628,  -1.02968,  -0.04547],
	[454,  -1.26701,  -1.01919,  -0.05550],
	[455,  -1.25695,  -1.00792,  -0.06538],
	[456,  -1.24585,  -0.99554,  -0.07467],
	[457,  -1.23361,  -0.98201,  -0.08325],
	[458,  -1.22011,  -0.96732,  -0.09105],
	[459,  -1.20527,  -0.95149,  -0.09805],
	[460,  -1.18899,  -0.93452,  -0.10419],
	[461,  -1.17126,  -0.91648,  -0.10951],
	[462,  -1.15240,  -0.89766,  -0.11443],
	[463,  -1.13283,  -0.87841,  -0.11945],
	[464,  -1.11296,  -0.85907,  -0.12507],
	[465,  -1.09318,  -0.84001,  -0.13179],
	[466,  -1.07386,  -0.82152,  -0.14003],
	[467,  -1.05507,  -0.80370,  -0.14990],
	[468,  -1.03686,  -0.78660,  -0.16142],
	[469,  -1.01925,  -0.77026,  -0.17462],
	[470,  -1.00228,  -0.75475,  -0.18953],
	[471,  -0.98596,  -0.74008,  -0.20613],
	[472,  -0.97020,  -0.72614,  -0.22432],
	[473,  -0.95488,  -0.71280,  -0.24396],
	[474,  -0.93991,  -0.69993,  -0.26490],
	[475,  -0.92518,  -0.68740,  -0.28700],
	[476,  -0.91058,  -0.67510,  -0.31012],
	[477,  -0.89611,  -0.66297,  -0.33404],
	[478,  -0.88176,  -0.65102,  -0.35856],
	[479,  -0.86753,  -0.63921,  -0.38347],
	[480,  -0.85342,  -0.62754,  -0.40856],
	[481,  -0.83944,  -0.61599,  -0.43369],
	[482,  -0.82560,  -0.60459,  -0.45894],
	[483,  -0.81195,  -0.59340,  -0.48449],
	[484,  -0.79850,  -0.58244,  -0.51050],
	[485,  -0.78528,  -0.57176,  -0.53712],
	[486,  -0.77229,  -0.56136,  -0.56443],
	[487,  -0.75929,  -0.55102,  -0.59215],
	[488,  -0.74604,  -0.54046,  -0.61989],
	[489,  -0.73226,  -0.52943,  -0.64728],
	[490,  -0.71770,  -0.51766,  -0.67394],
	[491,  -0.70217,  -0.50494,  -0.69958],
	[492,  -0.68575,  -0.49138,  -0.72429],
	[493,  -0.66860,  -0.47714,  -0.74822],
	[494,  -0.65089,  -0.46238,  -0.77156],
	[495,  -0.63278,  -0.44726,  -0.79446],
	[496,  -0.61441,  -0.43192,  -0.81711],
	[497,  -0.59582,  -0.41640,  -0.83977],
	[498,  -0.57707,  -0.40070,  -0.86273],
	[499,  -0.55817,  -0.38483,  -0.88627],
	[500,  -0.53916,  -0.36880,  -0.91066],
	[501,  -0.52009,  -0.35262,  -0.93617],
	[502,  -0.50099,  -0.33635,  -0.96288],
	[503,  -0.48192,  -0.32006,  -0.99090],
	[504,  -0.46292,  -0.30382,  -1.02028],
	[505,  -0.44404,  -0.28770,  -1.05112],
	[506,  -0.42533,  -0.27175,  -1.08337],
	[507,  -0.40683,  -0.25604,  -1.11656],
	[508,  -0.38857,  -0.24060,  -1.15011],
	[509,  -0.37059,  -0.22549,  -1.18343],
	[510,  -0.35293,  -0.21076,  -1.21595],
	[511,  -0.33561,  -0.19643,  -1.24725],
	[512,  -0.31867,  -0.18252,  -1.27765],
	[513,  -0.30215,  -0.16901,  -1.30765],
	[514,  -0.28606,  -0.15590,  -1.33774],
	[515,  -0.27044,  -0.14318,  -1.36843],
	[516,  -0.25532,  -0.13086,  -1.40010],
	[517,  -0.24080,  -0.11905,  -1.43271],
	[518,  -0.22695,  -0.10788,  -1.46609],
	[519,  -0.21387,  -0.09749,  -1.50010],
	[520,  -0.20165,  -0.08799,  -1.53457],
	[521,  -0.19035,  -0.07947,  -1.56938],
	[522,  -0.17986,  -0.07183,  -1.60454],
	[523,  -0.17006,  -0.06494,  -1.64008],
	[524,  -0.16081,  -0.05864,  -1.67602],
	[525,  -0.15198,  -0.05279,  -1.71240],
	[526,  -0.14348,  -0.04727,  -1.74923],
	[527,  -0.13530,  -0.04209,  -1.78644],
	[528,  -0.12749,  -0.03727,  -1.82397],
	[529,  -0.12010,  -0.03286,  -1.86171],
	[530,  -0.11315,  -0.02887,  -1.89959],
	[531,  -0.10668,  -0.02532,  -1.93755],
	[532,  -0.10057,  -0.02214,  -1.97565],
	[533,  -0.09470,  -0.01921,  -2.01400],
	[534,  -0.08894,  -0.01644,  -2.05269],
	[535,  -0.08317,  -0.01374,  -2.09181],
	[536,  -0.07730,  -0.01103,  -2.13142],
	[537,  -0.07143,  -0.00840,  -2.17147],
	[538,  -0.06567,  -0.00595,  -2.21187],
	[539,  -0.06015,  -0.00381,  -2.25253],
	[540,  -0.05502,  -0.00208,  -2.29337],
	[541,  -0.05037,  -0.00087,  -2.33431],
	[542,  -0.04622,  -0.00018,  -2.37534],
	[543,  -0.04255,  -0.00001,  -2.41648],
	[544,  -0.03937,  -0.00036,  -2.45772],
	[545,  -0.03665,  -0.00122,  -2.49909],
	[546,  -0.03438,  -0.00257,  -2.54058],
	[547,  -0.03241,  -0.00428,  -2.58220],
	[548,  -0.03059,  -0.00619,  -2.62397],
	[549,  -0.02876,  -0.00815,  -2.66589],
	[550,  -0.02678,  -0.01002,  -2.70798],
	[551,  -0.02453,  -0.01169,  -2.75022],
	[552,  -0.02211,  -0.01328,  -2.79259],
	[553,  -0.01964,  -0.01496,  -2.83501],
	[554,  -0.01728,  -0.01690,  -2.87744],
	[555,  -0.01514,  -0.01928,  -2.91982],
	[556,  -0.01334,  -0.02220,  -2.96210],
	[557,  -0.01183,  -0.02559,  -3.00431],
	[558,  -0.01051,  -0.02932,  -3.04645],
	[559,  -0.00931,  -0.03326,  -3.08856],
	[560,  -0.00813,  -0.03728,  -3.13067],
	[561,  -0.00692,  -0.04127,  -3.17278],
	[562,  -0.00570,  -0.04533,  -3.21489],
	[563,  -0.00451,  -0.04954,  -3.25698],
	[564,  -0.00339,  -0.05402,  -3.29904],
	[565,  -0.00240,  -0.05888,  -3.34105],
	[566,  -0.00157,  -0.06421,  -3.38300],
	[567,  -0.00092,  -0.06998,  -3.42488],
	[568,  -0.00043,  -0.07617,  -3.46669],
	[569,  -0.00013,  -0.08273,  -3.50841],
	[570,   0.00000,  -0.08964,  -3.55006],
	[571,  -0.00007,  -0.09687,  -3.59162],
	[572,  -0.00038,  -0.10449,  -3.63308],
	[573,  -0.00099,  -0.11259,  -3.67445],
	[574,  -0.00196,  -0.12126,  -3.71571],
	[575,  -0.00335,  -0.13060,  -3.75687],
	[576,  -0.00518,  -0.14064,  -3.79790],
	[577,  -0.00729,  -0.15126,  -3.83881],
	[578,  -0.00951,  -0.16230,  -3.87960],
	[579,  -0.01164,  -0.17357,  -3.92025],
	[580,  -0.01348,  -0.18490,  -3.96077],
	[581,  -0.01493,  -0.19619,  -4.00115],
	[582,  -0.01609,  -0.20746,  -4.04138],
	[583,  -0.01715,  -0.21881,  -4.08146],
	[584,  -0.01830,  -0.23035,  -4.12139],
	[585,  -0.01972,  -0.24215,  -4.16116],
	[586,  -0.02156,  -0.25431,  -4.20077],
	[587,  -0.02381,  -0.26689,  -4.24021],
	[588,  -0.02643,  -0.27991,  -4.27948],
	[589,  -0.02937,  -0.29344,  -4.31859],
	[590,  -0.03261,  -0.30751,  -4.35751],
	[591,  -0.03609,  -0.32215,  -4.39626],
	[592,  -0.03982,  -0.33734,  -4.43482],
	[593,  -0.04381,  -0.35306,  -4.47320],
	[594,  -0.04806,  -0.36925,  -4.51140],
	[595,  -0.05258,  -0.38590,  -4.54940],
	[596,  -0.05738,  -0.40297,  -4.58721],
	[597,  -0.06242,  -0.42047,  -4.62482],
	[598,  -0.06770,  -0.43841,  -4.66224],
	[599,  -0.07318,  -0.45682,  -4.69945],
	[600,  -0.07884,  -0.47570,  -4.73646],
	[601,  -0.08467,  -0.49506,  -4.77327],
	[602,  -0.09072,  -0.51490,  -4.80987],
	[603,  -0.09702,  -0.53518,  -4.84626],
	[604,  -0.10365,  -0.55587,  -4.88244],
	[605,  -0.11064,  -0.57696,  -4.91841],
	[606,  -0.11805,  -0.59843,  -4.95417],
	[607,  -0.12585,  -0.62025,  -4.98970],
	[608,  -0.13403,  -0.64240,  -5.02503],
	[609,  -0.14254,  -0.66488,  -5.06013],
	[610,  -0.15137,  -0.68767,  -5.09502],
	[611,  -0.16049,  -0.71075,  -5.12969],
	[612,  -0.16992,  -0.73413,  -5.16413],
	[613,  -0.17965,  -0.75782,  -5.19835],
	[614,  -0.18972,  -0.78184,  -5.23235],
	[615,  -0.20013,  -0.80620,  -5.26613],
	[616,  -0.21088,  -0.83091,  -float('inf')],
	[617,  -0.22193,  -0.85596,  -float('inf')],
	[618,  -0.23322,  -0.88135,  -float('inf')],
	[619,  -0.24470,  -0.90708,  -float('inf')],
	[620,  -0.25631,  -0.93315,  -float('inf')],
	[621,  -0.26803,  -0.95953,  -float('inf')],
	[622,  -0.27998,  -0.98621,  -float('inf')],
	[623,  -0.29232,  -1.01314,  -float('inf')],
	[624,  -0.30521,  -1.04028,  -float('inf')],
	[625,  -0.31881,  -1.06759,  -float('inf')],
	[626,  -0.33323,  -1.09505,  -float('inf')],
	[627,  -0.34839,  -1.12268,  -float('inf')],
	[628,  -0.36419,  -1.15050,  -float('inf')],
	[629,  -0.38048,  -1.17854,  -float('inf')],
	[630,  -0.39717,  -1.20682,  -float('inf')],
	[631,  -0.41414,  -1.23538,  -float('inf')],
	[632,  -0.43136,  -1.26418,  -float('inf')],
	[633,  -0.44880,  -1.29320,  -float('inf')],
	[634,  -0.46646,  -1.32240,  -float('inf')],
	[635,  -0.48431,  -1.35176,  -float('inf')],
	[636,  -0.50232,  -1.38126,  -float('inf')],
	[637,  -0.52048,  -1.41100,  -float('inf')],
	[638,  -0.53874,  -1.44107,  -float('inf')],
	[639,  -0.55708,  -1.47160,  -float('inf')],
	[640,  -0.57547,  -1.50268,  -float('inf')],
	[641,  -0.59390,  -1.53436,  -float('inf')],
	[642,  -0.61251,  -1.56638,  -float('inf')],
	[643,  -0.63146,  -1.59845,  -float('inf')],
	[644,  -0.65092,  -1.63025,  -float('inf')],
	[645,  -0.67104,  -1.66147,  -float('inf')],
	[646,  -0.69196,  -1.69190,  -float('inf')],
	[647,  -0.71361,  -1.72175,  -float('inf')],
	[648,  -0.73592,  -1.75136,  -float('inf')],
	[649,  -0.75879,  -1.78104,  -float('inf')],
	[650,  -0.78215,  -1.81113,  -float('inf')],
	[651,  -0.80589,  -1.84187,  -float('inf')],
	[652,  -0.82999,  -1.87323,  -float('inf')],
	[653,  -0.85439,  -1.90513,  -float('inf')],
	[654,  -0.87906,  -1.93746,  -float('inf')],
	[655,  -0.90396,  -1.97013,  -float('inf')],
	[656,  -0.92906,  -2.00305,  -float('inf')],
	[657,  -0.95436,  -2.03619,  -float('inf')],
	[658,  -0.97986,  -2.06950,  -float('inf')],
	[659,  -1.00556,  -2.10295,  -float('inf')],
	[660,  -1.03148,  -2.13653,  -float('inf')],
	[661,  -1.05761,  -2.17017,  -float('inf')],
	[662,  -1.08396,  -2.20379,  -float('inf')],
	[663,  -1.11051,  -2.23729,  -float('inf')],
	[664,  -1.13728,  -2.27055,  -float('inf')],
	[665,  -1.16425,  -2.30349,  -float('inf')],
	[666,  -1.19142,  -2.33602,  -float('inf')],
	[667,  -1.21879,  -2.36821,  -float('inf')],
	[668,  -1.24637,  -2.40018,  -float('inf')],
	[669,  -1.27417,  -2.43203,  -float('inf')],
	[670,  -1.30219,  -2.46386,  -float('inf')],
	[671,  -1.33044,  -2.49577,  -float('inf')],
	[672,  -1.35892,  -2.52775,  -float('inf')],
	[673,  -1.38765,  -2.55981,  -float('inf')],
	[674,  -1.41661,  -2.59194,  -float('inf')],
	[675,  -1.44584,  -2.62412,  -float('inf')],
	[676,  -1.47531,  -2.65636,  -float('inf')],
	[677,  -1.50503,  -2.68865,  -float('inf')],
	[678,  -1.53499,  -2.72100,  -float('inf')],
	[679,  -1.56516,  -2.75340,  -float('inf')],
	[680,  -1.59552,  -2.78586,  -float('inf')],
	[681,  -1.62609,  -2.81840,  -float('inf')],
	[682,  -1.65691,  -2.85105,  -float('inf')],
	[683,  -1.68805,  -2.88388,  -float('inf')],
	[684,  -1.71957,  -2.91694,  -float('inf')],
	[685,  -1.75153,  -2.95029,  -float('inf')],
	[686,  -1.78396,  -2.98395,  -float('inf')],
	[687,  -1.81670,  -3.01778,  -float('inf')],
	[688,  -1.84954,  -3.05161,  -float('inf')],
	[689,  -1.88228,  -3.08527,  -float('inf')],
	[690,  -1.91471,  -3.11859,  -float('inf')],
	[691,  -1.94668,  -3.15143,  -float('inf')],
	[692,  -1.97827,  -3.18384,  -float('inf')],
	[693,  -2.00960,  -3.21592,  -float('inf')],
	[694,  -2.04082,  -3.24777,  -float('inf')],
	[695,  -2.07203,  -3.27946,  -float('inf')],
	[696,  -2.10335,  -3.31109,  -float('inf')],
	[697,  -2.13477,  -3.34268,  -float('inf')],
	[698,  -2.16626,  -3.37424,  -float('inf')],
	[699,  -2.19779,  -3.40578,  -float('inf')],
	[700,  -2.22933,  -3.43733,  -float('inf')],
	[701,  -2.26087,  -3.46889,  -float('inf')],
	[702,  -2.29245,  -3.50051,  -float('inf')],
	[703,  -2.32414,  -3.53223,  -float('inf')],
	[704,  -2.35602,  -3.56410,  -float('inf')],
	[705,  -2.38814,  -3.59616,  -float('inf')],
	[706,  -2.42056,  -3.62844,  -float('inf')],
	[707,  -2.45323,  -3.66088,  -float('inf')],
	[708,  -2.48609,  -3.69342,  -float('inf')],
	[709,  -2.51909,  -3.72596,  -float('inf')],
	[710,  -2.55215,  -3.75845,  -float('inf')],
	[711,  -2.58522,  -3.79081,  -float('inf')],
	[712,  -2.61826,  -3.82303,  -float('inf')],
	[713,  -2.65119,  -3.85508,  -float('inf')],
	[714,  -2.68398,  -3.88695,  -float('inf')],
	[715,  -2.71657,  -3.91862,  -float('inf')],
	[716,  -2.74892,  -3.95010,  -float('inf')],
	[717,  -2.78110,  -3.98141,  -float('inf')],
	[718,  -2.81315,  -4.01261,  -float('inf')],
	[719,  -2.84515,  -4.04373,  -float('inf')],
	[720,  -2.87717,  -4.07483,  -float('inf')],
	[721,  -2.90925,  -4.10594,  -float('inf')],
	[722,  -2.94136,  -4.13702,  -float('inf')],
	[723,  -2.97344,  -4.16803,  -float('inf')],
	[724,  -3.00542,  -4.19891,  -float('inf')],
	[725,  -3.03726,  -4.22963,  -float('inf')],
	[726,  -3.06892,  -4.26014,  -float('inf')],
	[727,  -3.10040,  -4.29048,  -float('inf')],
	[728,  -3.13176,  -4.32069,  -float('inf')],
	[729,  -3.16303,  -4.35080,  -float('inf')],
	[730,  -3.19425,  -4.38086,  -float('inf')],
	[731,  -3.22544,  -4.41089,  -float('inf')],
	[732,  -3.25665,  -4.44093,  -float('inf')],
	[733,  -3.28790,  -4.47097,  -float('inf')],
	[734,  -3.31922,  -4.50104,  -float('inf')],
	[735,  -3.35063,  -4.53113,  -float('inf')],
	[736,  -3.38215,  -4.56125,  -float('inf')],
	[737,  -3.41370,  -4.59131,  -float('inf')],
	[738,  -3.44516,  -4.62121,  -float('inf')],
	[739,  -3.47643,  -4.65085,  -float('inf')],
	[740,  -3.50742,  -4.68015,  -float('inf')],
	[741,  -3.53805,  -4.70902,  -float('inf')],
	[742,  -3.56840,  -4.73757,  -float('inf')],
	[743,  -3.59859,  -4.76593,  -float('inf')],
	[744,  -3.62872,  -4.79423,  -float('inf')],
	[745,  -3.65890,  -4.82259,  -float('inf')],
	[746,  -3.68923,  -4.85111,  -float('inf')],
	[747,  -3.71967,  -4.87976,  -float('inf')],
	[748,  -3.75014,  -4.90847,  -float('inf')],
	[749,  -3.78058,  -4.93716,  -float('inf')],
	[750,  -3.81093,  -4.96577,  -float('inf')],
	[751,  -3.84113,  -4.99424,  -float('inf')],
	[752,  -3.87119,  -5.02257,  -float('inf')],
	[753,  -3.90110,  -5.05074,  -float('inf')],
	[754,  -3.93089,  -5.07876,  -float('inf')],
	[755,  -3.96056,  -5.10664,  -float('inf')],
	[756,  -3.99013,  -5.13438,  -float('inf')],
	[757,  -4.01962,  -5.16203,  -float('inf')],
	[758,  -4.04907,  -5.18961,  -float('inf')],
	[759,  -4.07851,  -5.21720,  -float('inf')],
	[760,  -4.10795,  -5.24482,  -float('inf')],
	[761,  -4.13743,  -5.27251,  -float('inf')],
	[762,  -4.16690,  -5.30022,  -float('inf')],
	[763,  -4.19631,  -5.32792,  -float('inf')],
	[764,  -4.22560,  -5.35553,  -float('inf')],
	[765,  -4.25472,  -5.38300,  -float('inf')],
	[766,  -4.28364,  -5.41030,  -float('inf')],
	[767,  -4.31241,  -5.43749,  -float('inf')],
	[768,  -4.34113,  -5.46462,  -float('inf')],
	[769,  -4.36987,  -5.49177,  -float('inf')],
	[770,  -4.39871,  -5.51901,  -float('inf')],
	[771,  -4.42770,  -5.54639,  -float('inf')],
	[772,  -4.45677,  -5.57382,  -float('inf')],
	[773,  -4.48581,  -5.60121,  -float('inf')],
	[774,  -4.51472,  -5.62844,  -float('inf')],
	[775,  -4.54339,  -5.65541,  -float('inf')],
	[776,  -4.57174,  -5.68205,  -float('inf')],
	[777,  -4.59983,  -5.70841,  -float('inf')],
	[778,  -4.62774,  -5.73457,  -float('inf')],
	[779,  -4.65556,  -5.76063,  -float('inf')],
	[780,  -4.68336,  -5.78666,  -float('inf')],
	[781,  -4.71121,  -5.81273,  -float('inf')],
	[782,  -4.73909,  -5.83884,  -float('inf')],
	[783,  -4.76698,  -5.86494,  -float('inf')],
	[784,  -4.79485,  -5.89101,  -float('inf')],
	[785,  -4.82266,  -5.91702,  -float('inf')],
	[786,  -4.85039,  -5.94295,  -float('inf')],
	[787,  -4.87805,  -5.96881,  -float('inf')],
	[788,  -4.90566,  -5.99462,  -float('inf')],
	[789,  -4.93323,  -6.02039,  -float('inf')],
	[790,  -4.96080,  -6.04616,  -float('inf')],
	[791,  -4.98836,  -6.07192,  -float('inf')],
	[792,  -5.01591,  -6.09765,  -float('inf')],
	[793,  -5.04340,  -6.12330,  -float('inf')],
	[794,  -5.07082,  -6.14883,  -float('inf')],
	[795,  -5.09813,  -6.17419,  -float('inf')],
	[796,  -5.12532,  -6.19935,  -float('inf')],
	[797,  -5.15238,  -6.22431,  -float('inf')],
	[798,  -5.17931,  -6.24911,  -float('inf')],
	[799,  -5.20612,  -6.27376,  -float('inf')],
	[800,  -5.23280,  -6.29827,  -float('inf')],
	[801,  -5.25937,  -6.32268,  -float('inf')],
	[802,  -5.28587,  -6.34703,  -float('inf')],
	[803,  -5.31235,  -6.37136,  -float('inf')],
	[804,  -5.33885,  -6.39572,  -float('inf')],
	[805,  -5.36542,  -6.42016,  -float('inf')],
	[806,  -5.39209,  -6.44472,  -float('inf')],
	[807,  -5.41883,  -6.46934,  -float('inf')],
	[808,  -5.44559,  -6.49399,  -float('inf')],
	[809,  -5.47230,  -6.51861,  -float('inf')],
	[810,  -5.49893,  -6.54314,  -float('inf')],
	[811,  -5.52542,  -6.56753,  -float('inf')],
	[812,  -5.55177,  -6.59179,  -float('inf')],
	[813,  -5.57797,  -6.61592,  -float('inf')],
	[814,  -5.60402,  -6.63992,  -float('inf')],
	[815,  -5.62992,  -6.66378,  -float('inf')],
	[816,  -5.65566,  -6.68753,  -float('inf')],
	[817,  -5.68128,  -6.71119,  -float('inf')],
	[818,  -5.70682,  -6.73482,  -float('inf')],
	[819,  -5.73231,  -6.75844,  -float('inf')],
	[820,  -5.75779,  -6.78210,  -float('inf')],
	[821,  -5.78329,  -6.80584,  -float('inf')],
	[822,  -5.80881,  -6.82965,  -float('inf')],
	[823,  -5.83432,  -6.85351,  -float('inf')],
	[824,  -5.85981,  -6.87741,  -float('inf')],
	[825,  -5.88525,  -6.90133,  -float('inf')],
	[826,  -5.91064,  -6.92525,  -float('inf')],
	[827,  -5.93596,  -6.94917,  -float('inf')],
	[828,  -5.96119,  -6.97307,  -float('inf')],
	[829,  -5.98631,  -6.99692,  -float('inf')],
	[830,  -6.01130,  -7.02072,  -float('inf')]
])

hcf = hcf2deg
if (args.q10deg): hcf = hcf10deg
if (args.e2deg): hcf = hcfe2deg

# find sensitivity of a given photoreceptor type to a wavelength using visual pigment templates
# from Govardovskii et al. (2000). l is short for lambda. See also https://pmc.ncbi.nlm.nih.gov/articles/PMC2962788/
# This produces the human cone fundamentals instead if the "lmax" (peak sensitivity) argument
# is set to -1 (L cone), -2 (M cone) or -3 (S cone), since negative values are not meaningful
# here. This option should not be used with ocular media filtering because it's already
# built in. The results are similar to -l 560 -m 530 -s 420 --filter human, but not identical
# because they include the macular pigment. So far I haven't found any tabulated data for
# this.
# Palacios et al. 2010 use a different equation for the beta-band peak (123 + 0.429 * alpha
# peak).
def vpt(w, lmax):
	# use human cone fundamentals for testing
	if (lmax == -1): # L cone
		if (w < 390): return 0
		else: return 10**(hcf[int(round(w))-390][1])
	if (lmax == -2): # M
		if (w < 390): return 0
		else: return 10**(hcf[int(round(w))-390][2])
	if (lmax == -3): # S
		if (w < 390): return 0
		else: return 10**(hcf[int(round(w))-390][3])
	
	# coefficients
	A = 69.7
	a = 0.8795 + 0.0459*math.exp(-(lmax - 300)**2 / 11940)
	B = 28
	b = 0.922
	C = -14.9
	c = 1.104
	D = 0.674
	Abeta = 0.26
	lmbeta = 189 + 0.315 * lmax
	b1 = -40.5 + 0.195 * lmax
	try:
		alpha = 1 / (math.exp(A*(a - lmax/w)) + math.exp(B*(b - lmax/w)) + math.exp(C*(c - lmax/w)) + D)
		beta = Abeta * math.exp(-((w - lmbeta) / b1)**2)
	except OverflowError:
		if (args.warnings): print("Warning (vpt): math overflow, clipping to 2.2250738585072014e-308")
		return 2.2250738585072014e-308
	
	# self-screening (see The Optics of Life equation 4.12 and https://pmc.ncbi.nlm.nih.gov/articles/PMC3269789/)
	# We have a problem: self-screening is described as broadening spectral sensitivity by
	# reducing sensitivity at the peak, but this formula actually does the opposite. This
	# happens because A(lambda) is defined as the result of the above equation and the larger
	# values near the peak produce a smaller value for the exponential function (due to the
	# negative sign), which is then subtracted from 1. It's correct that the highest optical
	# density is found at the peak of visual pigment absorption, but higher density leads to
	# less light being transmitted whereas this equation has less light transmitted at
	# lower densities.
	# This comes back to the issue of exactly what "optical density" represents and what unit
	# it's measured in. For the Wratten filters, converting optical density to transmission
	# is done with 10^-density. The self-screening equation is similar with e^-density but then
	# subtracts that from 1, giving us an inverse relationship. On the other hand, increasing
	# the peak density causes a broader curve, which is correct.
	# Lamb (1995) uses a more complex formula: log10[1 - S(1 - 10^-d)]/(-D), where S is the
	# pigment sensitivity/absorption. The key here seems to be dividing by the optical density
	# as well as subtracting it. When we add this feature, we see the expected relationship.
	value = alpha + beta
	if (args.selfscreen):
		return ((value * (1 - math.exp(-args.od*value*args.osl))) / (args.od*value*args.osl))
	elif (args.selfscreen1):
		return value * (1 - math.exp(-args.od*args.osl))
	
	return value

# human photopic luminosity function (2019 CIE standard)
# No values are provided below 360 nm, so we assume 0.
luminosity = np.array([
	0, # 300
	0,
	0,
	0,
	0,
	0, # 305
	0,
	0,
	0,
	0,
	0, # 310
	0,
	0,
	0,
	0,
	0, # 315
	0,
	0,
	0,
	0,
	0, # 320
	0,
	0,
	0,
	0,
	0, # 325
	0,
	0,
	0,
	0,
	0, # 330
	0,
	0,
	0,
	0,
	0, # 335
	0,
	0,
	0,
	0,
	0, # 340
	0,
	0,
	0,
	0,
	0, # 345
	0,
	0,
	0,
	0,
	0, # 350
	0,
	0,
	0,
	0,
	0, # 355
	0,
	0,
	0,
	0,
	0.000003917, # 360
	0.000004393581,
	0.000004929604,
	0.000005532136,
	0.000006208245,
	0.000006965, # 365
	0.000007813219,
	0.000008767336,
	0.000009839844,
	0.00001104323,
	0.00001239, # 370
	0.00001388641,
	0.00001555728,
	0.00001744296,
	0.00001958375,
	0.00002202, # 375
	0.00002483965,
	0.00002804126,
	0.00003153104,
	0.00003521521,
	0.000039, # 380
	0.0000428264,
	0.0000469146,
	0.0000515896,
	0.0000571764,
	0.000064, # 385
	0.00007234421,
	0.00008221224,
	0.00009350816,
	0.0001061361,
	0.00012, # 390
	0.000134984,
	0.000151492,
	0.000170208,
	0.000191816,
	0.000217, # 395
	0.0002469067,
	0.00028124,
	0.00031852,
	0.0003572667,
	0.000396, # 400
	0.0004337147,
	0.000473024,
	0.000517876,
	0.0005722187,
	0.00064, # 405
	0.00072456,
	0.0008255,
	0.00094116,
	0.00106988,
	0.00121, # 410
	0.001362091,
	0.001530752,
	0.001720368,
	0.001935323,
	0.00218, # 415
	0.0024548,
	0.002764,
	0.0031178,
	0.0035264,
	0.004, # 420
	0.00454624,
	0.00515932,
	0.00582928,
	0.00654616,
	0.0073, # 425
	0.008086507,
	0.00890872,
	0.00976768,
	0.01066443,
	0.0116, # 430
	0.01257317,
	0.01358272,
	0.01462968,
	0.01571509,
	0.01684, # 435
	0.01800736,
	0.01921448,
	0.02045392,
	0.02171824,
	0.023, # 440
	0.02429461,
	0.02561024,
	0.02695857,
	0.02835125,
	0.0298, # 445
	0.03131083,
	0.03288368,
	0.03452112,
	0.03622571,
	0.038, # 450
	0.03984667,
	0.041768,
	0.043766,
	0.04584267,
	0.048, # 455
	0.05024368,
	0.05257304,
	0.05498056,
	0.05745872,
	0.06, # 460
	0.06260197,
	0.06527752,
	0.06804208,
	0.07091109,
	0.0739, # 465
	0.077016,
	0.0802664,
	0.0836668,
	0.0872328,
	0.09098, # 470
	0.09491755,
	0.09904584,
	0.1033674,
	0.1078846,
	0.1126, # 475
	0.117532,
	0.1226744,
	0.1279928,
	0.1334528,
	0.13902, # 480
	0.1446764,
	0.1504693,
	0.1564619,
	0.1627177,
	0.1693, # 485
	0.1762431,
	0.1835581,
	0.1912735,
	0.199418,
	0.20802, # 490
	0.2171199,
	0.2267345,
	0.2368571,
	0.2474812,
	0.2586, # 495
	0.2701849,
	0.2822939,
	0.2950505,
	0.308578,
	0.323, # 500
	0.3384021,
	0.3546858,
	0.3716986,
	0.3892875,
	0.4073, # 505
	0.4256299,
	0.4443096,
	0.4633944,
	0.4829395,
	0.503, # 510
	0.5235693,
	0.544512,
	0.56569,
	0.5869653,
	0.6082, # 515
	0.6293456,
	0.6503068,
	0.6708752,
	0.6908424,
	0.71, # 520
	0.7281852,
	0.7454636,
	0.7619694,
	0.7778368,
	0.7932, # 525
	0.8081104,
	0.8224962,
	0.8363068,
	0.8494916,
	0.862, # 530
	0.8738108,
	0.8849624,
	0.8954936,
	0.9054432,
	0.9148501, # 535
	0.9237348,
	0.9320924,
	0.9399226,
	0.9472252,
	0.954, # 540
	0.9602561,
	0.9660074,
	0.9712606,
	0.9760225,
	0.9803, # 545
	0.9840924,
	0.9874182,
	0.9903128,
	0.9928116,
	0.9949501, # 550
	0.9967108,
	0.9980983,
	0.999112,
	0.9997482,
	1, # 555
	0.9998567,
	0.9993046,
	0.9983255,
	0.9968987,
	0.995, # 560
	0.9926005,
	0.9897426,
	0.9864444,
	0.9827241,
	0.9786, # 565
	0.9740837,
	0.9691712,
	0.9638568,
	0.9581349,
	0.952, # 570
	0.9454504,
	0.9384992,
	0.9311628,
	0.9234576,
	0.9154, # 575
	0.9070064,
	0.8982772,
	0.8892048,
	0.8797816,
	0.87, # 580
	0.8598613,
	0.849392,
	0.838622,
	0.8275813,
	0.8163, # 585
	0.8047947,
	0.793082,
	0.781192,
	0.7691547,
	0.757, # 590
	0.7447541,
	0.7324224,
	0.7200036,
	0.7074965,
	0.6949, # 595
	0.6822192,
	0.6694716,
	0.6566744,
	0.6438448,
	0.631, # 600
	0.6181555,
	0.6053144,
	0.5924756,
	0.5796379,
	0.5668, # 605
	0.5539611,
	0.5411372,
	0.5283528,
	0.5156323,
	0.503, # 610
	0.4904688,
	0.4780304,
	0.4656776,
	0.4534032,
	0.4412, # 615
	0.42908,
	0.417036,
	0.405032,
	0.393032,
	0.381, # 620
	0.3689184,
	0.3568272,
	0.3447768,
	0.3328176,
	0.321, # 625
	0.3093381,
	0.2978504,
	0.2865936,
	0.2756245,
	0.265, # 630
	0.2547632,
	0.2448896,
	0.2353344,
	0.2260528,
	0.217, # 635
	0.2081616,
	0.1995488,
	0.1911552,
	0.1829744,
	0.175, # 640
	0.1672235,
	0.1596464,
	0.1522776,
	0.1451259,
	0.1382, # 645
	0.1315003,
	0.1250248,
	0.1187792,
	0.1127691,
	0.107, # 650
	0.1014762,
	0.09618864,
	0.09112296,
	0.08626485,
	0.0816, # 655
	0.07712064,
	0.07282552,
	0.06871008,
	0.06476976,
	0.061, # 660
	0.05739621,
	0.05395504,
	0.05067376,
	0.04754965,
	0.04458, # 665
	0.04175872,
	0.03908496,
	0.03656384,
	0.03420048,
	0.032, # 670
	0.02996261,
	0.02807664,
	0.02632936,
	0.02470805,
	0.0232, # 675
	0.02180077,
	0.02050112,
	0.01928108,
	0.01812069,
	0.017, # 680
	0.01590379,
	0.01483718,
	0.01381068,
	0.01283478,
	0.01192, # 685
	0.01106831,
	0.01027339,
	0.009533311,
	0.008846157,
	0.00821, # 690
	0.007623781,
	0.007085424,
	0.006591476,
	0.006138485,
	0.005723, # 695
	0.005343059,
	0.004995796,
	0.004676404,
	0.004380075,
	0.004102, # 700
	0.003838453,
	0.003589099,
	0.003354219,
	0.003134093,
	0.002929, # 705
	0.002738139,
	0.002559876,
	0.002393244,
	0.002237275,
	0.002091, # 710
	0.001953587,
	0.00182458,
	0.00170358,
	0.001590187,
	0.001484, # 715
	0.001384496,
	0.001291268,
	0.001204092,
	0.001122744,
	0.001047, # 720
	0.0009765896,
	0.0009111088,
	0.0008501332,
	0.0007932384,
	0.00074, # 725
	0.0006900827,
	0.00064331,
	0.000599496,
	0.0005584547,
	0.00052, # 730
	0.0004839136,
	0.0004500528,
	0.0004183452,
	0.0003887184,
	0.0003611, # 735
	0.0003353835,
	0.0003114404,
	0.0002891656,
	0.0002684539,
	0.0002492, # 740
	0.0002313019,
	0.0002146856,
	0.0001992884,
	0.0001850475,
	0.0001719, # 745
	0.0001597781,
	0.0001486044,
	0.0001383016,
	0.0001287925,
	0.00012, # 750
	0.0001118595,
	0.0001043224,
	0.0000973356,
	0.00009084587,
	0.0000848, # 755
	0.00007914667,
	0.000073858,
	0.000068916,
	0.00006430267,
	0.00006, # 760
	0.00005598187,
	0.0000522256,
	0.0000487184,
	0.00004544747,
	0.0000424, # 765
	0.00003956104,
	0.00003691512,
	0.00003444868,
	0.00003214816,
	0.00003, # 770
	0.00002799125,
	0.00002611356,
	0.00002436024,
	0.00002272461,
	0.0000212, # 775
	0.00001977855,
	0.00001845285,
	0.00001721687,
	0.00001606459,
	0.00001499, # 780
	0.00001398728,
	0.00001305155,
	0.00001217818,
	0.00001136254,
	0.0000106, # 785
	0.000009885877,
	0.000009217304,
	0.000008592362,
	0.000008009133,
	0.0000074657, # 790
	0.000006959567,
	0.000006487995,
	0.000006048699,
	0.000005639396,
	0.0000052578, # 795
	0.000004901771,
	0.00000456972,
	0.000004260194,
	0.000003971739,
	0.0000037029, # 800
	0.000003452163,
	0.000003218302,
	0.0000030003,
	0.000002797139,
	0.0000026078, # 805
	0.00000243122,
	0.000002266531,
	0.000002113013,
	0.000001969943,
	0.0000018366, # 810
	0.00000171223,
	0.000001596228,
	0.00000148809,
	0.000001387314,
	0.0000012934, # 815
	0.00000120582,
	0.000001124143,
	0.000001048009,
	0.0000009770578,
	0.00000091093, # 820
	0.0000008492513,
	0.0000007917212,
	0.0000007380904,
	0.0000006881098,
	0.00000064153, # 825
	0.0000005980895,
	0.0000005575746,
	0.000000519808,
	0.0000004846123,
	0.00000045181 # 830
])

# absolute number of photons, part 2: scaling D65 to a specified value in lux
# For values not included in Python's D65, we approximate this with a black body with
# temperature 6504 K. (Never mind, I found the real thing)

# integral of D65 x human photopic luminosity function (we hardcode this instead of using
# sensitivity() because that varies with settings for the specified species)
# The units we want at the end are W/m^2, so we multiply this by the lux-to-watts scaling
# factor for 555 nm (683.002 lm/W). This will end up on the bottom.
integral = 0
for i in range(531):
	integral += d65_1nm[i] * luminosity[i] * 683.002

# lux scaling factor (units = lx / lx/nm = nm)
lxscale = args.lux / integral
#print(lxscale)

# scaled version (W/m^2 * m / joule*sec * m/s * nm = joule/sec/m^2 * m / joule/sec * m/s * nm * sr)
# Per nm is implied. The values we just found do not include steradians (probably), so we divide
# by the sr value we found earlier. (No we don't)
# Also changing this to 1 nm.
d65a = np.empty(531)
for i in range(d65a.shape[0]):
	w = i + 300
	d65a[i] = lxscale*1e-9*w*d65_1nm[i] / (h*c*math.pi)
	#print(d65a[i])

e_abs = np.empty(531)
for i in range(e_abs.shape[0]):
	w = i + 300
	e_abs[i] = lxscale*1e-9*w*100 / (h*c*math.pi)

# white point
if (args.white == "d65"):
	wp = d65
	wp_1nm = d65_1nm
elif (args.white == "a"):
	wp = a
elif (args.white == "i"):
	wp = incandescent
	wp_1nm = incandescent_1nm
else:
	wp = e
	wp_1nm = e_1nm

# Reviving the CMF calculation to create a perceptual color space. I kept the names r(ed),
# g(reen) and b(lue) for convenience even though they have nothing to do with how the chosen
# wavelengths appear to humans.
# Note that two of the CIE RGB primary wavelengths (700, 546.1 and 435.8) are not integers.
# Non-integer values are supported, but when using the human cone fundamentals they're
# rounded off because we look up the sensitivity in an array rather than calculating it directly.
# This probably doesn't matter much.
# This part prevents the script from running if we don't specify all of L, M and S as separate
# values, so it's been made conditional on this. Determining a set of dichromatic CMFs would
# follow a different process.
if (l1 != m1 != s1):
	# primaries
	r = args.red
	g = args.green
	b = args.blue

	# sensitivity of cones to primaries
	matrix_a = np.array([
		[vpt(r, l1)*lens_filter(r), vpt(g, l1)*lens_filter(g), vpt(b, l1)*lens_filter(b)],
		[vpt(r, m1)*lens_filter(r), vpt(g, m1)*lens_filter(g), vpt(b, m1)*lens_filter(b)],
		[vpt(r, s1)*lens_filter(r), vpt(g, s1)*lens_filter(g), vpt(b, s1)*lens_filter(b)]
	])

	# sensitivity of cones to each wavelength from 300 to 800 nm in 1 nm increments
	matrix_c = np.empty((3, 401)) # 500x3 matrix -- this is height x width, not width x height

	# red row
	for i in range(401):
		matrix_c[0][i] = vpt(i + 300, l1)*lens_filter(i+300)
	# green row
		matrix_c[1][i] = vpt(i + 300, m1)*lens_filter(i+300)
	# blue row
		matrix_c[2][i] = vpt(i + 300, s1)*lens_filter(i+300)

	# initial color match matrix
	matrix_m = np.matmul(np.linalg.inv(matrix_a), matrix_c) # A x M = C, so M = A^-1 x C

	# adjust to reference white
	#watts = 0
	#for i in range(401): watts += wp_1nm[i] # scale to 1 watt
	rw = 0.0
	gw = 0.0
	bw = 0.0
	for i in range (401):
		rw += matrix_m[0][i] * wp_1nm[i]#/watts
		gw += matrix_m[1][i] * wp_1nm[i]#/watts
		bw += matrix_m[2][i] * wp_1nm[i]#/watts

	# final color match matrix
	matrix_cmf = np.empty((3, 401))
	matrix_cmf[0] = matrix_m[0] / rw
	matrix_cmf[1] = matrix_m[1] / gw
	matrix_cmf[2] = matrix_m[2] / bw

if (args.cmf):
	print("White units: r=" + str(rw) + ", g=" + str(gw) + ", b=" + str(bw))
	# maxima, minima, intersections and negative numbers
	maxrx = 0
	maxgx = 0
	maxbx = 0
	maxry = 0
	maxgy = 0
	maxby = 0
	minrx = 0
	mingx = 0
	minbx = 0
	minry = 0
	mingy = 0
	minby = 0
	intrg = 0
	intgb = 0
	rneg = 0
	gneg = 0
	bneg = 0
	for i in range(401):
		# find maxima
		if (matrix_cmf[0][i] > maxry):
			maxry = matrix_cmf[0][i]
			maxrx = i+300
		if (matrix_cmf[1][i] > maxgy):
			maxgy = matrix_cmf[1][i]
			maxgx = i+300
		if (matrix_cmf[2][i] > maxby):
			maxby = matrix_cmf[2][i]
			maxbx = i+300
		# find minima
		if (matrix_cmf[0][i] < minry):
			minry = matrix_cmf[0][i]
			minrx = i+300
		if (matrix_cmf[1][i] < mingy):
			mingy = matrix_cmf[1][i]
			mingx = i+300
		if (matrix_cmf[2][i] < minby):
			minby = matrix_cmf[2][i]
			minbx = i+300
		# find intersections
		if ((matrix_cmf[0][i] - matrix_cmf[1][i]) >= 0 and (matrix_cmf[0][i-1] - matrix_cmf[1][i-1]) <= 0): intrg = i+300
		if ((matrix_cmf[1][i] - matrix_cmf[2][i]) >= 0 and (matrix_cmf[1][i-1] - matrix_cmf[2][i-1]) <= 0): intgb = i+300
		# find negatives
		if (matrix_cmf[0][i] < 0): rneg += matrix_cmf[0][i]
		if (matrix_cmf[1][i] < 0): gneg += matrix_cmf[1][i]
		if (matrix_cmf[2][i] < 0): bneg += matrix_cmf[2][i]
	print("Peaks: r=" + str(maxrx) + ", g=" + str(maxgx) + ", b=" + str(maxbx))
	print("Troughs: r=" + str(minrx) + ", g=" + str(mingx) + ", b=" + str(minbx))
	print("Intersections: rg=" + str(intrg) + ", gb=" + str(intgb))
	print("Total negative values: r=" + str(rneg) + ", g=" + str(gneg) + ", b=" + str(bneg) + ", total=" + str(rneg+gneg+bneg))
	
	# plot curves
	xvalues = np.empty(401)
	for i in range(401): xvalues[i] = i+300
	plt.plot(xvalues, matrix_m[0], 'r')
	plt.plot(xvalues, matrix_m[1], 'g')
	plt.plot(xvalues, matrix_m[2], 'b')
	plt.show()
	plt.plot(xvalues, matrix_cmf[0], 'r')
	plt.plot(xvalues, matrix_cmf[1], 'g')
	plt.plot(xvalues, matrix_cmf[2], 'b')
	plt.show()

# gamma correction (see color_conversions.py from the colormath module)
def gamma(v):
	if v <= 0.0031308:
                result = v * 12.92
	else:
		result = 1.055 * math.pow(v, 1 / 2.4) - 0.055
	return result

# estimate hue, saturation and lightness for a spectral power distribution
# Lens filtering is now used for colors properly with a von Kries transform.
# 03/03/2025 -- added a "reflect" argument for specifying whether this is a
# reflectance/transmission spectrum and should be multiplied by our light source. If it's
# an emission spectrum as in the sky and illuminant examples, we don't want this.
def spectral_rendering(table, light_source=wp_1nm, output=True, lw=l1, mw=m1, sw=s1, reflect=True):
	table_l = 0
	table_m = 0
	table_s = 0
	table_r = 0
	table_g = 0
	table_b = 0
	brightness = 0
	
	# interpolate
	if(len(table) < 401):
		x10nm = np.empty(41)
		for i in range(41):
			x10nm[i] = i*10 + 300
		
		x1nm = np.empty(401)
		for i in range(401):
			x1nm[i] = i + 300
		
		table = np.interp(x1nm, x10nm, table)
	
	for i in range(0, table.shape[0]):
		w = i + 300
		if (table[i] > 0): # don't add zeroes or nans
			if (reflect == True): table_l += vpt(w, lw) * table[i] * light_source[i] * lens_filter(w)
			else: table_l += vpt(w, lw) * table[i] * lens_filter(w)
			# remove either M or S for dichromacy
			if (mw != lw):
				if (reflect == True): table_m += vpt(w, mw) * table[i] * light_source[i] * lens_filter(w)
				else: table_m += vpt(w, mw) * table[i] * lens_filter(w)
			if (sw != mw):
				if (reflect == True): table_s += vpt(w, sw) * table[i] * light_source[i] * lens_filter(w)
				else: table_s += vpt(w, sw) * table[i] * lens_filter(w)
			
			# brightness
			if (reflect == True): brightness += sensitivity(w) * table[i] * light_source[i]
			else: brightness += sensitivity(w) * table[i]
			
			# CMF
			if (lw != mw != sw):
				if (reflect == True): table_r += matrix_cmf[0][i] * table[i] * light_source[i]
				else: table_r += matrix_cmf[0][i] * table[i]
				if (reflect == True): table_g += matrix_cmf[1][i] * table[i] * light_source[i]
				else: table_g += matrix_cmf[1][i] * table[i]
				if (reflect == True): table_b += matrix_cmf[2][i] * table[i] * light_source[i]
				else: table_b += matrix_cmf[2][i] * table[i]
	
	# normalize according to provided light source
	n = 0
	wpl = 0
	wpm = 0
	wps = 0
	for i in range(0, table.shape[0]):
		w = i + 300
		n += sensitivity(w) * light_source[i]
		wpl += vpt(w, lw) * light_source[i] * lens_filter(w)
		wpm += vpt(w, mw) * light_source[i] * lens_filter(w)
		wps += vpt(w, sw) * light_source[i] * lens_filter(w)
	
	# von Kries transform
	if (args.novk == False):
		table_l = table_l / wpl
		table_m = table_m / wpm
		table_s = table_s / wps
	
	# gamma correction for CMF colors
	table_r = gamma(table_r)
	table_g = gamma(table_g)
	table_b = gamma(table_b)
	
	# express brightness of reflected/emitted light proportional to the light source
	brightness = brightness / n
	
	# L-M position
	lm = (table_l - table_m) / (table_l + table_m)
	
	if (output):
		# estimate color
		# For trichromacy we convert the LMS values to other color spaces to roughly visualize
		# the hue and saturation. This assumes a given ratio of cone absorption produces the
		# same perception as it does in humans and can be adequately modeled by an LMS<->XYZ
		# transformation matrix I found on Wikipedia. For dichromacy and monochromacy this
		# doesn't really work, so we just report the estimated brightness and (for dichromacy)
		# the L:S difference.
		# Now estimates trichromatic colors based on a color triangle. The wavelengths
		# reported agree reasonably well with what the hue should be. I've fixed this so
		# it uses proper Maxwell coordinates, lens filtering and a von Kries transform.
		if (lw != mw != sw):
			print("Cone response: l=" + str(table_l) + ", m=" + str(table_m) + ", s=" + str(table_s))
			print("Position on L-M axis: " + str(lm))
			
			# simple false color rendering
			scale_max = max(table_l,table_m,table_s)
			#scale_eq = (table_l+table_m+table_s)
			print("RGB=LMS false color: (" + str(table_l * 255) + ", " + str(table_m * 255) + ", " + str(table_s * 255) + ")")
			print("RGB=LMS false color (maximum brightness): (" + str(table_l * 255 / scale_max) + ", " + str(table_m * 255 / scale_max) + ", " + str(table_s * 255 / scale_max) + ")")
			#print("RGB=LMS false color (equal brightness): (" + str(table_l * 255 / scale_eq) + ", " + str(table_m * 255 / scale_eq) + ", " + str(table_s * 255 / scale_eq) + ")")
			hsl = convert_color(sRGBColor(table_l, table_m, table_s), HSLColor)
			print("LMS->HSL: " + str(hsl.get_value_tuple()))
			
			# CMF false color rendering. Out-of-gamut colors are scaled down if
			# too large and clipped toward the white point if negative.
			rgb_max = max(table_r,table_g,table_b)
			rgb_min = min(table_r,table_g,table_b,0)
			rgb_max2 = max(table_r-rgb_min,table_g-rgb_min,table_b-rgb_min,rgb_max,1)
			print("CMF false color: (" + str(table_r * 255) + ", " + str(table_g * 255) + ", " + str(table_b * 255) + ")")
			#print("CMF false color (gamma corrected): (" + str(gamma(table_r) * 255) + ", " + str(gamma(table_g) * 255) + ", " + str(gamma(table_b) * 255) + ")")
			print("CMF false color (maximum brightness): (" + str(table_r * 255 / rgb_max) + ", " + str(table_g * 255 / rgb_max) + ", " + str(table_b * 255 / rgb_max) + ")")
			print("CMF false color (within gamut): (" + str((table_r-rgb_min) * 255 / rgb_max2) + ", " + str((table_g-rgb_min) * 255 / rgb_max2) + ", " + str((table_b-rgb_min) * 255 / rgb_max2) + ")")
			hsl = convert_color(sRGBColor(table_r-rgb_min, table_g-rgb_min, table_b-rgb_min), HSLColor)
			print("CMF->HSL: " + str(hsl.get_value_tuple()))
				
			
			# This should be the nearest wavelength on a line from the white point.
			match = 300
			x0 = math.sqrt(1/2)*(table_l - table_s) / (table_l + table_m + table_s)
			y0 = math.sqrt(2/3)*(table_m - 0.5*(table_l + table_s)) / (table_l + table_m + table_s)
			wx = math.sqrt(1/2)*(wpl - wps) / (wpl + wpm + wps)
			wy = math.sqrt(2/3)*(wpm - 0.5*(wpl + wps)) / (wpl + wpm + wps)
			
			# adjust to chosen visible range
			while (lens_filter(match) == 0):
				match += 1
			
			for i in range(300, 701):
				matchl1 = vpt(match, lw)*lens_filter(match) / wpl
				matchm1 = vpt(match, mw)*lens_filter(match) / wpm
				matchs1 = vpt(match, sw)*lens_filter(match) / wps
				matchl2 = vpt(i, lw)*lens_filter(i) / wpl
				matchm2 = vpt(i, mw)*lens_filter(i) / wpm
				matchs2 = vpt(i, sw)*lens_filter(i) / wps
				x1 = math.sqrt(1/2)*(matchl1 - matchs1) / (matchl1 + matchm1 + matchs1)
				y1 = math.sqrt(2/3)*(matchm1 - 0.5*(matchl1 + matchs1)) / (matchl1 + matchm1 + matchs1)
				x2 = math.sqrt(1/2)*(matchl2 - matchs2) / (matchl2 + matchm2 + matchs2)
				y2 = math.sqrt(2/3)*(matchm2 - 0.5*(matchl2 + matchs2)) / (matchl2 + matchm2 + matchs2)
				d1 = np.linalg.norm(np.cross(np.asarray((wx, wy)) - np.asarray((x0, y0)), np.asarray((x0, y0)) - np.asarray((x1, y1)))) / np.linalg.norm(np.asarray((wx, wy)) - np.asarray((x0, y0)))
				d2 = np.linalg.norm(np.cross(np.asarray((wx, wy)) - np.asarray((x0, y0)), np.asarray((x0, y0)) - np.asarray((x2, y2)))) / np.linalg.norm(np.asarray((wx, wy)) - np.asarray((x0, y0)))
				
				# This should be on the same "side" of the white point, not just on the line.
				# Otherwise we have the complementary color.
				if (math.dist((x2, y2), (x0, y0)) < math.dist((x2, y2), (wx, wy)) and d2 < d1):
					match = i
			print("Closest wavelength: " + str(match))
			
			posx = math.sqrt(1/2)*(table_l - table_s) / (table_l + table_m + table_s)
			posy = math.sqrt(2/3)*(table_m - (table_l + table_s)/2) / (table_l + table_m + table_s)
			dist = math.dist((posx, posy), (wx, wy))
			print("Position in color triangle: " + str((posx, posy)))
			
			# These conversions turned out to be not very useful.
			#lms = np.array([
			#	[table_l],
			#	[table_m],
			#	[table_s]
			#])
			
			#if (args.render):
				#matrix = np.matmul(lms_to_xyz, lms)
				#xyz = XYZColor(*matrix, illuminant="e")
				#print("Color coordinates:")
				#print("CIE XYZ: " + str(xyz.get_value_tuple()))
				#print("sRGB: " + str(convert_color(xyz, sRGBColor).get_value_tuple()))
				#print("CIE LAB: " + str(convert_color(xyz, LabColor).get_value_tuple()))
		elif (lw == mw != sw):
			# typical dichromacy
			print("Cone response: l/m=" + str(table_l) + ", s=" + str(table_s))
			
			match = 300
			for i in range(300, 800):
				diff0 = abs(table_l / table_s - vpt(match, lw) / vpt(match, sw))
				diff1 = abs(table_l / table_s - vpt(i, lw) / vpt(i, sw))
				if (diff1 < diff0):
					match = i
			print("Matching wavelength: " + str(match))
			
			pos = (table_l - table_s) / (table_l + table_s)
			dist = pos - (wpl - wps) / (wpl + wps)
			print("Position on color line: " + str(pos))
		elif (lw != mw == sw):
			# tritanopia
			print("Cone response: l=" + str(table_l) + ", m=" + str(table_m))
		elif (lw == mw == sw):
			# monochromacy
			print("Cone response: " + str(table_l))
		
		# estimate brightness
		print("Brightness relative to light source: " + str(brightness))
		
		# "SpectralColor" conversion for comparison to check if our estimates are reasonable
		if (args.render):
			# colormath doesn't handle nans, so replace them with zeroes first
			for i in range(table.shape[0]):
				if (np.isnan(table[i])): table[i] = 0
			
			spectral = SpectralColor(spec_340nm=table[40],
			spec_350nm=table[50],
			spec_360nm=table[60],
			spec_370nm=table[70],
			spec_380nm=table[80],
			spec_390nm=table[90],
			spec_400nm=table[100],
			spec_410nm=table[110],
			spec_420nm=table[120],
			spec_430nm=table[130],
			spec_440nm=table[140],
			spec_450nm=table[150],
			spec_460nm=table[160],
			spec_470nm=table[170],
			spec_480nm=table[180],
			spec_490nm=table[190],
			spec_500nm=table[200],
			spec_510nm=table[210],
			spec_520nm=table[220],
			spec_530nm=table[230],
			spec_540nm=table[240],
			spec_550nm=table[250],
			spec_560nm=table[260],
			spec_570nm=table[270],
			spec_580nm=table[280],
			spec_590nm=table[290],
			spec_600nm=table[300],
			spec_610nm=table[310],
			spec_620nm=table[320],
			spec_630nm=table[330],
			spec_640nm=table[340],
			spec_650nm=table[350],
			spec_660nm=table[360],
			spec_670nm=table[370],
			spec_680nm=table[380],
			spec_690nm=table[390],
			spec_700nm=table[400],
			spec_710nm=0,
			spec_720nm=0,
			spec_730nm=0,
			spec_740nm=0,
			spec_750nm=0,
			spec_760nm=0,
			spec_770nm=0,
			spec_780nm=0,
			spec_790nm=0,
			spec_800nm=0,
			spec_810nm=0,
			spec_820nm=0,
			spec_830nm=0, illuminant='e')
			
			# If an illuminant is set, we'll also set it here. Note that colormath
			# only adapts the CMFs and not the SPD itself, so this may not do what
			# we want. D65 vs. E doesn't make much difference because they're very
			# similar for human vision.
			if (args.white == 'd65'): spectral.set_illuminant('d65')
			if (args.white == 'a' or args.white == 'i'): spectral.set_illuminant('a')
			
			srgb = convert_color(spectral, sRGBColor)
			srgb_tuple = srgb.get_value_tuple()
			scale = max(srgb_tuple[0], srgb_tuple[1], srgb_tuple[2]) / 255
			scale_eq = (srgb_tuple[0] + srgb_tuple[1] + srgb_tuple[2]) / 255
			print("colormath conversion:")
			print("XYZ: " + str(convert_color(spectral, XYZColor).get_value_tuple()))
			print("sRGB (0.0-1.0): " + str(srgb.get_value_tuple()))
			print("sRGB (0-255): " + str(srgb.get_upscaled_value_tuple()))
			#print("Lab: " + str(convert_color(spectral, LabColor).get_value_tuple()))
			print("HSL: " + str(convert_color(spectral, HSLColor).get_value_tuple()))
			if (scale > 0): # avoid dividing by 0
				print("sRGB (maximum brightness): (" + str(srgb_tuple[0] / scale) + ", " + str(srgb_tuple[1] / scale) + ", " + str(srgb_tuple[2] / scale) + ")")
				print("sRGB (equal brightness): (" + str(srgb_tuple[0] / scale_eq) + ", " + str(srgb_tuple[1] / scale_eq) + ", " + str(srgb_tuple[2] / scale_eq) + ")")
		
		# break up text
		print("")
	
		# return brightness and color for comparisons
		if (lw != mw):
			return(brightness, posx, posy)
		elif (mw != s1): # add zeroes for compatibility
			return(brightness, pos, 0)
		else:
			return(brightness, 0, 0)
	# if in no-output mode, just return the LM value
	return lm

# contrast sensitivity (Vorobyev and Osorio (1998); Vorobyev and Osorio (2001))
# This predicts that two very small values like 1e-5 and 1e-6 can produce a
# large amount of contrast despite probably neither of them being noticeable.
# Is this where quantum noise comes in? An alternative is to use the linear
# model rather than log-linear, which makes the contrast proportional to the
# size of the signals, but this reduces the contrast too much. For reference,
# human trichromacy should show all the Kodak Wratten colors as "different" and
# human-like dichromacy should show red, yellow and probably green as "the same".
# Neither of these models do both of these. (This has been fixed.)
# Absolute quantum catch integrals need values at 1-nm intervals, so we find them
# with linear interpolation if we don't already have them.
def color_contrast(table1, table2, quantum_noise=args.qn, lw = args.lw, mw = args.mw, sw = args.sw):
	# interpolate 1-nm intervals if we're provided with 10-nm
	if (len(table1) < 401):
		x10nm = np.empty(41)
		for i in range(41):
			x10nm[i] = i*10 + 300
		
		x1nm = np.empty(401)
		for i in range(401):
			x1nm[i] = i + 300
		
		table1 = np.interp(x1nm, x10nm, table1)
	if (len(table2) < 401):
		x10nm = np.empty(41)
		for i in range(41):
			x10nm[i] = i*10 + 300
		
		x1nm = np.empty(401)
		for i in range(401):
			x1nm[i] = i + 300
		
		table2 = np.interp(x1nm, x10nm, table2)
	
	# background light
	wpl = 0
	wpm = 0
	wps = 0
	
	ql1 = 0
	qm1 = 0
	qs1 = 0
	ql2 = 0
	qm2 = 0
	qs2 = 0
	for i in range(0, table1.shape[0]):
		w = i + 300
		if (table1[i] > 0): # zero/nan check
			ql1 += vpt(w, lw) * table1[i] * wp_1nm[i] * lens_filter(w)
			qm1 += vpt(w, mw) * table1[i] * wp_1nm[i] * lens_filter(w)
			qs1 += vpt(w, sw) * table1[i] * wp_1nm[i] * lens_filter(w)
		wpl += vpt(w, lw) * wp_1nm[i] * lens_filter(w)
		wpm += vpt(w, mw) * wp_1nm[i] * lens_filter(w)
		wps += vpt(w, sw) * wp_1nm[i] * lens_filter(w)
	for i in range(0, table2.shape[0]):
		w = i + 300
		if (table2[i] > 0):
			ql2 += vpt(w, lw) * table2[i] * wp_1nm[i] * lens_filter(w)
			qm2 += vpt(w, mw) * table2[i] * wp_1nm[i] * lens_filter(w)
			qs2 += vpt(w, sw) * table2[i] * wp_1nm[i] * lens_filter(w)
	
	# normalize
	if (args.novk == False):
		ql1 = ql1 / wpl
		qm1 = qm1 / wpm
		qs1 = qs1 / wps
		ql2 = ql2 / wpl
		qm2 = qm2 / wpm
		qs2 = qs2 / wps
	
	# differences
	dfl = math.log(ql1 / ql2)
	dfm = math.log(qm1 / qm2)
	dfs = math.log(qs1 / qs2)
	
	# quantum noise (needs work)
	# Based on these equations:
	# https://journals.biologists.com/jeb/article/218/2/184/14274/Bird-colour-vision-behavioural-thresholds-reveal
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5043323/
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1691374/pdf/12816638.pdf
	if (quantum_noise):
		# v: number of cones per receptive field (spatial summation)
		# d: receptor diameter (um)
		# f: focal length (um, cancels d)
		# D: pupil diameter (um)
		# O/Ai: transmittance of ocular media. "Ai" is for oil droplets, so we don't need it
		# because marsupial oil droplets are transparent.
		# T: integration time (ms), inverse of flicker fusion frequency (Hz)
		# ui/k: optical density of pigment (um^-1)
		# l: length of outer segment (um, cancels ui)
		# I tried to keep this as simple as possible to avoid making too many assumptions,
		# since almost none of the relevant information is available for any opossum
		# species, so the equation I use is equivalent to the "simple model" without
		# any terms for refractive indices or oil droplets. This means it probably
		# produces overestimates. Including oil droplets would be more complicated than
		# in birds or fish because in opossums there are separate classes of cones with
		# and without them.
		# I'm not sure how the units work. The amount of light going in should have units
		# of photons/m^3/s/sr, and the value we get should have just photons, so we need
		# to have units of m^3*s*sr for the product of these terms. I understand where
		# s and sr come from (integration time and size/distance), but what about meters?
		# All the terms with meters cancel out except for D, which is squared and gives us
		# m^2. Where's the other meter term? According to the Yale page and Vorobyev (2003),
		# the photons are not really expressed in cubic meters but in square meters "per
		# wavelength" with the wavelength part given as nanometers or micrometers, which
		# reduces the number by a factor of 1e-9 or 1e-6 and implies we have some extra
		# term that cancels this unit. Do we multiply it by the wavelength?
		# The answers to this StackExchange question says the missing term is the range of
		# wavelengths being integrated over: https://physics.stackexchange.com/questions/690117/what-does-spectral-flux-density-mean-per-wavelength We're doing this one wavelength at
		# a time, so this term would be 1 nm (I think) and we can then ignore it
		# aside from the units and order of magnitude.
		# According to The Optics of Life, the "per nanometer" unit is there because
		# the energy/photons are binned into 1-nm intervals. In this case we're using
		# 10 nm intervals, which would suggest we should multiply (or divide?) by 10 nm,
		# but this probably isn't needed because at the end I multiply the sum by 600 nm
		# (value of dλ, the whole range). Without the 600 nm the numbers look way too
		# small. (Never mind, I'm pretty sure this was wrong)
		
		# Estimates used:
		# Where necessary, I've converted centimeters, millimeters and micrometers to meters.
		# v = 3000 * 450: cone density for Didelphis aurita (Ahnelt et al. 2009), retinal
		# area for cats (Vaney 1985). This number is probably an overestimate because
		# 3000/mm^2 is the maximum density, not the average density. (Not needed, see
		# above)
		# * The retinal area in D. virginiana is 100-180 (Kolb & Wang 1985), so they
		# have smaller eyes than cats. Not sure how you derive the axial/focal length
		# from this. (never mind, we have it)
		# v = 3: cone to ganglion cell convergence in D. virginiana, finally the right
		# species for once (Kolb & Wang, 1985; cited in Vlahos et al. 2014) now how do
		# I turn this into radians? (You don't, you already have the angle from the
		# d/f term. It's basically the same as increasing d.) I don't know if this applies
		# to all cone classes equally or if some are more convergent than others. The
		# chicken study seems to assume the latter. With the greater sparsity of S cones,
		# surely connecting three of them to one ganglion cell and excluding L cones
		# from that connection would be more difficult.
		# dt = 1/22.4: 22.4 Hz = 45 ms, brushtail possum (Signal, Temple & Foster 2011)
		# d = 1.2 / 1000000: 1.2 um, mice (Nikonov et al. 2006, Fu & Yau 2007)
		# f = 6.81 / 1000: Didelphis aurita (Oswaldo-Cruz, Hokoc & Sousa 1979)
		# D = 6 / 1000: Didelphis aurita (")
		# O: mouse (see the lens filter functions)
		# k = 0.015: mice (Yin et al. 2013)
		# l = 30: Didelphis aurita (Ahnelt et al. 1995)
		# sphere = 4*math.pi*5**2: the filters were located 5 cm away from the light bulbs
		# sr = math.pi*3.8**2 / sphere: 3.8 cm diameter circle (oops, this is wrong, it's
		# diameter not radius -- divide by 2) (this was still wrong, a sphere has 4pi
		# steradians not 1) (also we don't need this at all)
		# K = 0.5 (see The Optics of Life)
		
		# Some values for other species for comparison:
		# area of retina: 1094 mm^2 (humans), 80 mm^2 (rats; Mayhew & Astle 1997)
		# total number of cones: ~6 million (humans)
		# outer segment length: 30 um (chickens), 6.2 um (cats; Fisher, Pfeffer & Anderson
		# 1983), 7 um (gray squirrels; "), 20 um (honey possums), 25.1-40 um (humans),
		# 15.2 um (goldfish red single cones, Harosi & Flamarique 2012), 5 um (tammar
		# wallaby and fat-tailed dunnart; Ebeling, Natoli & Hemmi 2010)
		# outer segment width: 1.73 um (chickens), 6.1 um (goldfish red single cones)
		# optical density: 0.015 (coral reef triggerfish, birds)
		# focal length: 8300 um = 8.3 mm (chickens), 2.6 mm (mice; Geng et al. 2011),
		# 22.3 mm (humans; "), 2.5 mm (triggerfish)
		# pupil size: 4900-3500 um = 4.9-3.5 mm (chickens), 2 mm (mice), 6 mm (humans)
		# transmittance of ocular media: 80% (general estimate from The Optics of Light)
		# integration time/flicker fusion: 50-12 ms = 20-83 Hz (chickens),
		# 70-80 Hz = 14-12 ms (cats and dogs), 18 Hz = 56 ms (rats; Gil, Valente & Shemesh 2024),
		# 14 Hz = 71 ms (mice; "), 60 Hz = 17 ms (humans), 25 Hz = 40 ms (some fish)
		
		v = args.convergence
		d = args.osd / 1000000
		f = args.fl / 1000
		D = args.pupil / 1000
		K = 0.5
		dt = 1 / args.ff
		k = args.od
		l = args.osl
		# We've already calculated the steradians elsewhere, so we won't do it again.
		# Actually forget the steradians. (pi/4)(d/f)^2 is an angle, so we already have
		# it. The size of the visual stimulus isn't supposed to be in there.
		
		# lens filtering scaling factor
		#T = 0.8 / lens_filter(700)
		T = 1 / lens_filter(700) # set to 1 at 700 nm
		
		# select illuminant
		if (args.white == "i"):
			photons = ia
		elif (args.white == "d65"):
			photons = d65a
		else:
			photons = e_abs
	
		aql1 = 0
		aqm1 = 0
		aqs1 = 0
		aql2 = 0
		aqm2 = 0
		aqs2 = 0
		for i in range(table1.shape[0]):
			w = i + 300
			if (table1[i] > 0):
				aql1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, lw)*l)) * vpt(w, lw) * table1[i] * photons[i] * lens_filter(w)
				aqm1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, mw)*l)) * vpt(w, mw) * table1[i] * photons[i] * lens_filter(w) 
				aqs1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, sw)*l)) * vpt(w, sw) * table1[i] * photons[i] * lens_filter(w)
		for i in range(table2.shape[0]):
			w = i + 300
			if (table2[i] > 0):
				aql2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, lw)*l)) * vpt(w, lw) * table2[i] * photons[i] * lens_filter(w)
				aqm2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, mw)*l)) * vpt(w, mw) * table2[i] * photons[i] * lens_filter(w)
				aqs2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, sw)*l)) * vpt(w, sw) * table2[i] * photons[i] * lens_filter(w)
		
		# adjust size of integrals to compensate for sampling every 10 nm
		# The Optics of Life recommends multiplying by the size of the
		# wavelength interval, but I don't think that's what's called for
		# here. The size seems to be underestimated by a factor of about 10,
		# not a factor the size of the whole interval, which is what you would
		# intuitively expect. The reason that "looked right" is I was confusing
		# energy and radiance (per steradian) elsewhere.
		# We don't need this anymore -- see above
		#integral_full = 0
		#integral_10nm = 0
		#if (args.white == 'd65'): # use D65 table
		#	for i in range(d65_1nm.shape[0]):
		#		integral_full += d65_1nm[i]
		#	for i in range(0, int(d65_1nm.shape[0]/10)):
		#		integral_10nm += d65_1nm[i*10]
		#	iscale = integral_full / integral_10nm
		#else:
		#	integral_full = quad(blackbody, args.start, args.end, args=args.ct)
		#	for i in range(int(args.start/10), int((args.end+1)/10)):
		#		integral_10nm += blackbody(i*10, args.ct)
		#	iscale = integral_full[0] / integral_10nm
		#aql1 = aql1 * iscale
		#aqm1 = aqm1 * iscale
		#aqs1 = aqs1 * iscale
		#aql2 = aql2 * iscale
		#aqm2 = aqm2 * iscale
		#aqs2 = aqs2 * iscale
		
		# This defines "number of cones per receptive field" to mean the convergence
		# ratio for each type of cone. This can in fact be less than 1 (see Vlahos
		# et al. 2014). For example, if 90% of cones are L cones, 10% are S cones
		# and the cone-to-ganglion cell ratio is 3:1, there are 2.7 L cones and
		# 0.3 S cones for every ganglion cell. I don't know if this is really what
		# that means, though. Does 1 receptive field = 1 ganglion cell?
		# Redoing this again so the minimum is 1, which is how it was done for
		# chickens and honey possums. The default for v is now "1". Setting it to 0
		# disables this feature.
		if (v != 0):
			aql1 = aql1*v*args.lp
			aql2 = aql2*v*args.lp
			aqm1 = aqm1*v*args.mp
			aqm2 = aqm2*v*args.mp
			aqs1 = aqs1*v*args.sp
			aqs2 = aqs2*v*args.sp
		
		el = math.sqrt((1 / aql1 + 1 / aql2) + 2*wl**2)
		em = math.sqrt((1 / aqm1 + 1 / aqm2) + 2*wm**2)
		es = math.sqrt((1 / aqs1 + 1 / aqs2) + 2*ws**2)
		
		if (args.qcheck):
			print("L1: " + str(aql1) + ", M1: " + str(aqm1) + ", S1: " + str(aqs1))
			print("L2: " + str(aql2) + ", M2: " + str(aqm2) + ", S2: " + str(aqs2))
			print("el: " + str(el) + ", em: " + str(em) + ", es: " + str(es))
			print("wl: " + str(wl) + ", wm: " + str(wm) + ", ws: " + str(ws))
	else:
		el = wl
		em = wm
		es = ws
	
	if (l1 != m1):
		delta_s = (es**2*(dfl - dfm)**2 + em**2*(dfl - dfs)**2 + el**2*(dfs - dfm)**2) / ((el*em)**2 + (el*es)**2 + (em*es)**2)
	else:
		delta_s = (dfl - dfs)**2 / (el**2 + es**2)
	return math.sqrt(delta_s)

# brightness contrast based on https://journals.biologists.com/jeb/article/207/14/2471/14748/Interspecific-and-intraspecific-views-of-color
def brightness_contrast(table1, table2):
	# interpolate 1-nm intervals if we're provided with 10-nm
	if (len(table1) < 401):
		x10nm = np.empty(41)
		for i in range(41):
			x10nm[i] = i*10 + 300
		
		x1nm = np.empty(401)
		for i in range(401):
			x1nm[i] = i + 300
		
		table1 = np.interp(x1nm, x10nm, table1)
	if (len(table2) < 401):
		x10nm = np.empty(41)
		for i in range(41):
			x10nm[i] = i*10 + 300
		
		x1nm = np.empty(401)
		for i in range(401):
			x1nm[i] = i + 300
		
		table2 = np.interp(x1nm, x10nm, table2)
	
	# background light
	# This doesn't seem to be used for anything? I don't think a von Kries transform is
	# relevant here. It gets divided out anyway.
	wpt = 0
	
	#ql1 = 0 # L cones
	#qm1 = 0 # M cones
	#qr1 = 0 # rods
	#ql2 = 0
	#qm2 = 0
	#qr2 = 0
	# We're redoing this to use sensitivity() so we can use the photopic luminosity
	# function for humans.
	q1 = 0
	q2 = 0
	for i in range(table1.shape[0]):
		w = i + 300
		wpt += sensitivity(w, l1) * wp_1nm[i]
		#ql1 += vpt(w, l1) * table1[i] * wp[i] * lens_filter(w)
		#qm1 += vpt(w, m1) * table1[i] * wp[i] * lens_filter(w)
		#qr1 += vpt(w, rod) * table1[i] * wp[i] * lens_filter(w)
		q1 += sensitivity(w) * table1[i] * wp_1nm[i]
	for i in range(table2.shape[0]):
		w = i + 300
		#ql2 += vpt(w, l1) * table2[i] * wp[i] * lens_filter(w)
		#qm2 += vpt(w, m1) * table2[i] * wp[i] * lens_filter(w)
		#qr2 += vpt(w, rod) * table2[i] * wp[i] * lens_filter(w)
		q2 += sensitivity(w) * table2[i] * wp_1nm[i]
	
	#if (l1 != m1):
	#	df = math.log((args.lb*ql1 + args.mb*qm1 + args.rb*qr1) / (args.lb*ql2 + args.mb*qm2 + args.rb*qr2))
	#else:
	#	df = math.log((args.lb*ql1 + args.rb*qr1) / (args.lb*ql2 + args.rb*qr2))
	df = math.log(q1 / q2)
	delta_s = abs(df / args.weberb)
	return delta_s

# print specified information
print("L: " + str(l1) + ", M: " + str(m1) + ", S: " + str(s1) + ", rod: " + str(rod))
print("")

# plot visual pigment templates
# I changed the line colors back to red, green and blue. I've been trying to avoid this kind
# of color coding, but in this case it's generally understood which visual pigment template
# is which based on its position on the graph, and the shades-of-gray "color" scheme was too
# pale for my liking.
if (args.vpt):
	xvalues = np.empty(400)
	lvalues = np.empty(400)
	mvalues = np.empty(400)
	svalues = np.empty(400)
	rvalues = np.empty(400)
	lmax = 0
	mmax = 0
	smax = 0
	rmax = 0
	for i in range(0, 400):
		xvalues[i] = i + 300
		if (args.filter != "none"):
			lvalues[i] = vpt(i+300, l1) * lens_filter(i+300)
			mvalues[i] = vpt(i+300, m1) * lens_filter(i+300)
			svalues[i] = vpt(i+300, s1) * lens_filter(i+300)
			rvalues[i] = vpt(i+300, rod) * lens_filter(i+300)
			lmax = max(lmax, lvalues[i])
			mmax = max(mmax, mvalues[i])
			smax = max(smax, svalues[i])
			rmax = max(rmax, rvalues[i])
		else:
			lvalues[i] = vpt(i+300, l1)
			mvalues[i] = vpt(i+300, m1)
			svalues[i] = vpt(i+300, s1)
			rvalues[i] = vpt(i+300, rod)
			lmax = 1
			mmax = 1
			smax = 1
			rmax = 1

	plt.plot(xvalues, lvalues/lmax, 'r', label="L (" + str(args.lw) + ")")
	# don't plot redundant curves
	if (l1 != m1):
		plt.plot(xvalues, mvalues/mmax, 'g', label="M (" + str(args.mw) + ")")
	if (m1 != s1):
		plt.plot(xvalues, svalues/smax, 'b', label="S (" + str(args.sw) + ")")
	if (args.rod != 0):
		plt.plot(xvalues, rvalues/rmax, ':k', label="rod (" + str(args.rod) + ")")
	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.legend()
	plt.show()

# show luminous efficiency function and lens filtering
if (args.luminosity):
	xvalues = np.empty(401)
	yvalues = np.empty(401)
	yvalues1 = np.empty(401)
	ms = 300
	for i in range(300, 701):
		if (sensitivity(i) > (sensitivity(ms))):
			ms = i
		xvalues[i-300] = i
		yvalues[i-300] = sensitivity(i)
		yvalues1[i-300] = lens_filter(i)
	print("Maximum sensitivity: " + str(ms))
	print("")
	
	plt.plot(xvalues, yvalues/sensitivity(ms), 'k')
	if (args.filter != "none"):
		plt.plot(xvalues, yvalues1, ':k')
	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.show()

# primary/secondary colors
if (args.primaries):
		
	# Maxwell triangle
		
	# primaries
	# For red and blue I've rounded off the values because wavelengths at the
	# ends of the spectrum get progressively closer to the vertices but the
	# differences stop being noticeable well before we get to the "closest"
	# point. For red we also look at the difference from the previous wavelength.
	# For blue I don't do this because with UV cones the points at that end
	# are bunched together and the "canonical" dunnart triangle suggests they
	# should be more spread out.
	r = 300
	for i in range(30, 80):
		rx = (vpt(r, l1) - vpt(r, s1)) / (vpt(r, l1) + vpt(r, m1) + vpt(r, s1))
		ry = (vpt(r, m1) - vpt(r, l1) - vpt(r, s1)) / (vpt(r, l1) + vpt(r, m1) + vpt(r, s1))
		ix = (vpt(i*10, l1) - vpt(i*10, s1)) / (vpt(i*10, l1) + vpt(i*10, m1) + vpt(i*10, s1))
		iy = (vpt(i*10, m1) - vpt(i*10, l1) - vpt(i*10, s1)) / (vpt(i*10, l1) + vpt(i*10, m1) + vpt(i*10, s1))
		ix1 = (vpt((i-1)*10, l1) - vpt((i-1)*10, s1)) / (vpt((i-1)*10, l1) + vpt((i-1)*10, m1) + vpt((i-1)*10, s1))
		iy1 = (vpt((i-1)*10, m1) - vpt((i-1)*10, l1) - vpt((i-1)*10, s1)) / (vpt((i-1)*10, l1) + vpt((i-1)*10, m1) + vpt((i-1)*10, s1))
		ix2 = (vpt((i+1)*10, l1) - vpt((i+1)*10, s1)) / (vpt((i+1)*10, l1) + vpt((i+1)*10, m1) + vpt((i+1)*10, s1))
		iy2 = (vpt((i+1)*10, m1) - vpt((i+1)*10, l1) - vpt((i+1)*10, s1)) / (vpt((i+1)*10, l1) + vpt((i+1)*10, m1) + vpt((i+1)*10, s1))
		diff0 = round(abs(math.dist((1, -1), (rx, ry))), 1) # distance to triangle vertex
		diff1 = round(abs(math.dist((1, -1), (ix, iy))), 1)
		diff2 = round(abs(math.dist((ix1, iy1), (ix, iy))), 1) # distance to previous
		diff3 = round(abs(math.dist((ix2, iy2), (ix, iy))), 1) # distance to next
		if (diff1 < diff0 and diff3 < diff2):
			r = i*10
	
	g = 300
	for i in range(30, 80):
		gx = (vpt(g, l1) - vpt(g, s1)) / (vpt(g, l1) + vpt(g, m1) + vpt(g, s1))
		gy = (vpt(g, m1) - vpt(g, l1) - vpt(g, s1)) / (vpt(g, l1) + vpt(g, m1) + vpt(g, s1))
		ix = (vpt(i*10, l1) - vpt(i*10, s1)) / (vpt(i*10, l1) + vpt(i*10, m1) + vpt(i*10, s1))
		iy = (vpt(i*10, m1) - vpt(i*10, l1) - vpt(i*10, s1)) / (vpt(i*10, l1) + vpt(i*10, m1) + vpt(i*10, s1))
		diff0 = abs(math.dist((0, 1), (gx, gy)))
		diff1 = abs(math.dist((0, 1), (ix, iy)))
		# distance to LM line: "green" is not only the closest point to the
		# M vertex but seems to be found on the LM side of the triangle.
		# For humans, the closest point is around 510-520 nm but is not on
		# the LM line, whereas the closest point that is on the line is
		# well within the range recognized as "green". For marsupial color
		# spaces with UV cones, the closest point also lies on the line
		# and is thus a primary color.
		diffrg = round(np.linalg.norm(np.cross((1, -1) - np.asarray((0, 1)), np.asarray((0, 1)) - np.asarray((ix, iy)))) / np.linalg.norm(np.asarray((1, -1)) - np.asarray((0, 1))), 2)
		if (diff1 < diff0 and diffrg == 0):
			g = i*10
	
	# The "blue" we get for human/primate values is really a violet. The
	# issue may be similar to green.
	b = 800
	for i in range(30, 80):
		bx = (vpt(b, l1) - vpt(b, s1)) / (vpt(b, l1) + vpt(b, m1) + vpt(b, s1))
		by = (vpt(b, m1) - vpt(b, l1) - vpt(b, s1)) / (vpt(b, l1) + vpt(b, m1) + vpt(b, s1))
		ix = (vpt(i*10, l1) - vpt(i*10, s1)) / (vpt(i*10, l1) + vpt(i*10, m1) + vpt(i*10, s1))
		iy = (vpt(i*10, m1) - vpt(i*10, l1) - vpt(i*10, s1)) / (vpt(i*10, l1) + vpt(i*10, m1) + vpt(i*10, s1))
		diff0 = round(abs(math.dist((-1, -1), (bx, by))), 2)
		diff1 = round(abs(math.dist((-1, -1), (ix, iy))), 2)
		if (diff1 < diff0):
			b = i*10
	
	# secondaries
	# As in the study, these are found by drawing a line from the corresponding
	# primary through the white point and finding the nearest wavelength to
	# the other end of the line. This is easy to do by hand but surprisingly
	# difficult to do in Python. The np.polyfit/Polynomial.fit method does not
	# work. How hard is it to find a line passing through two points? I don't
	# get it.
	white_l = 0
	white_m = 0
	white_s = 0
	for i in range(0, 50):
		w = i*10 + 340
		if (w >= 300):
			white_l += vpt(w, l1) * wp[i]
			white_m += vpt(w, m1) * wp[i]
			white_s += vpt(w, s1) * wp[i]
	
	c = s1
	x1 = (vpt(r, l1) - vpt(r, s1)) / (vpt(r, l1) + vpt(r, m1) + vpt(r, s1))
	y1 = (vpt(r, m1) - vpt(r, l1) - vpt(r, s1)) / (vpt(r, l1) + vpt(r, m1) + vpt(r, s1))
	x2 = (white_l - white_s) / (white_l + white_m + white_s)
	y2 = (white_m - white_l - white_s) / (white_l + white_m + white_s)
	xy1 = (x1, y1)
	xy2 = (x2, y2)
	
	for i in range(int(s1/5), int(m1/5)):
		x3 = (vpt(c, l1) - vpt(c, s1)) / (vpt(c, l1) + vpt(c, m1) + vpt(c, s1))
		y3 = (vpt(c, m1) - vpt(c, l1) - vpt(c, s1)) / (vpt(c, l1) + vpt(c, m1) + vpt(c, s1))
		x4 = (vpt(i*5, l1) - vpt(i*5, s1)) / (vpt(i*5, l1) + vpt(i*5, m1) + vpt(i*5, s1))
		y4 = (vpt(i*5, m1) - vpt(i*5, l1) - vpt(i*5, s1)) / (vpt(i*5, l1) + vpt(i*5, m1) + vpt(i*5, s1))
		
		d3 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x3, y3)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))
		d4 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x4, y4)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))

		if (d4 < d3):
			c = i*5
	
	y = m1
	x1 = (vpt(b, l1) - vpt(b, s1)) / (vpt(b, l1) + vpt(b, m1) + vpt(b, s1))
	y1 = (vpt(b, m1) - vpt(b, l1) - vpt(b, s1)) / (vpt(b, l1) + vpt(b, m1) + vpt(b, s1))
		
	for i in range(int(m1/5), 160):
		x3 = (vpt(y, l1) - vpt(y, s1)) / (vpt(y, l1) + vpt(y, m1) + vpt(y, s1))
		y3 = (vpt(y, m1) - vpt(y, l1) - vpt(y, s1)) / (vpt(y, l1) + vpt(y, m1) + vpt(y, s1))
		x4 = (vpt(i*5, l1) - vpt(i*5, s1)) / (vpt(i*5, l1) + vpt(i*5, m1) + vpt(i*5, s1))
		y4 = (vpt(i*5, m1) - vpt(i*5, l1) - vpt(i*5, s1)) / (vpt(i*5, l1) + vpt(i*5, m1) + vpt(i*5, s1))
			
		d3 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x3, y3)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))
		d4 = np.linalg.norm(np.cross((x2, y2) - np.asarray((x1, y1)), np.asarray((x1, y1)) - np.asarray((x4, y4)))) / np.linalg.norm(np.asarray((x2, y2)) - np.asarray((x1, y1)))
			
		if (d4 < d3):
			y = i*5
		
	print("Estimated primary colors based on color triangle:")
	print("red: " + str(r))
	print("yellow: " + str(y))
	print("green: " + str(g))
	print("cyan: " + str(c))
	print("blue: " + str(b))
	
# plot color triangle
# Now a real Maxwell triangle with Cartesian coordinates. Note these triangles often (more commonly?)
# have M on the left and S on top rather than the other way around as I do here. I chose this
# orientation because I find it more intuitive.
if (args.triangle):
	el = 0
	em = 0
	es = 0
	d65l = 0
	d65m = 0
	d65s = 0
	wl = 0
	wm = 0
	ws = 0
	xvalues = np.empty(41)
	yvalues = np.empty(41)
	labels = np.empty(41)
	
	# white
	for i in range(0, 41):
		w = i*10 + 300
		wl += vpt(w, l1) * wp[i] * lens_filter(w)
		wm += vpt(w, m1) * wp[i] * lens_filter(w)
		ws += vpt(w, s1) * wp[i] * lens_filter(w)
		el += vpt(w, l1) * e[i] * lens_filter(w)
		em += vpt(w, m1) * e[i] * lens_filter(w)
		es += vpt(w, s1) * e[i] * lens_filter(w)
		d65l += vpt(w, l1) * d65[i] * lens_filter(w)
		d65m += vpt(w, m1) * d65[i] * lens_filter(w)
		d65s += vpt(w, s1) * d65[i] * lens_filter(w)
	
	for i in range(0, 41):
		w = i*10 + 300
		labels[i] = w
		l = vpt(w, l1) * lens_filter(w)
		m = vpt(w, m1) * lens_filter(w)
		s = vpt(w, s1) * lens_filter(w)
		if (args.novk == False):
			l = l / wl
			m = m / wm
			s = s / ws
		total = l + m + s
		l = l / total
		m = m / total
		s = s / total
		xvalues[i] = math.sqrt(1/2)*(l - s)
		yvalues[i] = math.sqrt(2/3)*(m - (l + s)/2)
	plt.plot(xvalues, yvalues, '-k')
	for i in range(0, 41):
		plt.plot(xvalues[i], yvalues[i], 'ok')
		plt.text(xvalues[i], yvalues[i], labels[i])
	
	# E and D65
	if (args.novk == False):
		el = el / wl
		em = em / wm
		es = es / ws
		d65l = d65l / wl
		d65m = d65m / wm
		d65s = d65s / ws
	
	ex = math.sqrt(1/2)*(el - es) / (el + em + es)
	ey = math.sqrt(2/3)*(em - (el + es)/2) / (el + em + es)
	d65x = math.sqrt(1/2)*(d65l - d65s) / (d65l + d65m + d65s)
	d65y = math.sqrt(2/3)*(d65m - (d65l + d65s)/2) / (d65l + d65m + d65s)
	plt.plot(ex, ey, 'ok')
	plt.plot(d65x, d65y, 'ok')
	plt.text(ex, ey, 'E')
	plt.text(d65x, d65y, 'D65')
	
	xborder = np.array([-math.sqrt(1/2), 0, math.sqrt(1/2), -math.sqrt(1/2)])
	yborder = np.array([-math.sqrt(2/3)/2, math.sqrt(2/3), -math.sqrt(2/3)/2, -math.sqrt(2/3)/2])
	plt.plot(xborder, yborder, '-k')
	plt.text(-math.sqrt(1/2) - 0.05, -math.sqrt(2/3)/2 - 0.025, 'S')
	plt.text(0 - 0.025, math.sqrt(2/3) + 0.0125, 'M')
	plt.text(math.sqrt(1/2) + 0.0125, -math.sqrt(2/3)/2 - 0.025, 'L')
	plt.show()

# receptor noise-scaled version
if (args.triangle1):
	el = 0
	em = 0
	es = 0
	d65l = 0
	d65m = 0
	d65s = 0
	xvalues = np.empty(37)
	yvalues = np.empty(37)
	labels = np.empty(37)
	
	a = ws**2 / (ws**2 + wl**2)
	b = wl**2 / (ws**2 + wl**2)
	A = math.sqrt(1 / (ws**2 + wl**2))
	B = math.sqrt((ws**2 + wl**2) / ((wm*ws)**2 + (wm*wl)**2 + (ws*wl)**2))
		
	for i in range(0, 37):
		w = i*10 + 340
		labels[i] = w
		el += vpt(w, l1) * e[i]
		em += vpt(w, m1) * e[i]
		es += vpt(w, s1) * e[i]
		d65l += vpt(w, l1) * d65[i]
		d65m += vpt(w, m1) * d65[i]
		d65s += vpt(w, s1) * d65[i]
		fl = math.log(vpt(w, l1) / (vpt(w, l1) + vpt(w, m1) + vpt(w, s1)))
		fm = math.log(vpt(w, m1) / (vpt(w, l1) + vpt(w, m1) + vpt(w, s1)))
		fs = math.log(vpt(w, s1) / (vpt(w, l1) + vpt(w, m1) + vpt(w, s1)))
		xvalues[i] = A*(fl - fs)
		yvalues[i] = B*(fm - (a*fl + b*fs))
	plt.plot(xvalues, yvalues, '-k')
	for i in range(0, 17):
		plt.plot(xvalues[i*2+1], yvalues[i*2+1], 'ok')
		plt.text(xvalues[i*2+1], yvalues[i*2+1], labels[i*2+1])
	efl = math.log(el)
	efm = math.log(em)
	efs = math.log(es)
	d65fl = math.log(d65l)
	d65fm = math.log(d65m)
	d65fs = math.log(d65s)
	ex = A*(efl - efs) / (efl + efm + efs)
	ey = B*(efm - (a*efl + b*efs)) / (efl + efm + efs)
	d65x = A*(d65fl - d65fs) / (d65fl + d65fm + d65fs)
	d65y = B*(d65fm - (a*d65fl + b*d65fs)) / (d65fl + d65fm + d65fs)
	plt.plot(ex, ey, 'ok')
	plt.plot(d65x, d65y, 'ok')
	plt.text(ex, ey, 'E')
	plt.text(d65x, d65y, 'D65')
	
	#xborder = np.array([-math.sqrt(1/2), 0, math.sqrt(1/2), -math.sqrt(1/2)])
	#yborder = np.array([-math.sqrt(2/3)/2, math.sqrt(2/3), -math.sqrt(2/3)/2, -math.sqrt(2/3)/2])
	#plt.plot(xborder, yborder, '-k')
	#plt.text(-math.sqrt(1/2) - 0.05, -math.sqrt(2/3)/2 - 0.025, 'S')
	#plt.text(0 - 0.025, math.sqrt(2/3) + 0.0125, 'M')
	#plt.text(math.sqrt(1/2) + 0.0125, -math.sqrt(2/3)/2 - 0.025, 'L')
	plt.show()
	
# plot color hexagon (Chittka 1992)
# This doesn't quite work as advertised. The shape is supposed to be more of
# a circle, and the white point should probably be in the middle.
if (args.hexagon):
	# convert quantum catch to excitation
	wpl = 0
	wpm = 0
	wps = 0
	if (args.white == 'i'):
		start = 300
		end = 900
	else:
		start = 340
		end = 830
	t = 0
	for i in range(wp.shape[0]):
		w = i*10 + start
		wpl += vpt(w, l1) * wp[i]
		wpm += vpt(w, m1) * wp[i]
		wps += vpt(w, s1) * wp[i]
		t += wp[i]
	wpl = wpl/t
	wpm = wpm/t
	wps = wps/t
	wle = (wpl * wpl) / (wpl * wpl + 1)
	wme = (wpm * wpm) / (wpm * wpm + 1)
	wse = (wps * wps) / (wps * wps + 1)
	
	lvmax = 0
	mvmax = 0
	svmax = 0
	
	# find maximum values
	for i in range(0, 41):
		lvmax = max((wpl * vpt(i*10 + 300, l1) / (wpl * vpt(i*10 + 300, l1) + 1)), lvmax)
		mvmax = max((wpm * vpt(i*10 + 300, m1) / (wpm * vpt(i*10 + 300, m1) + 1)), mvmax)
		svmax = max((wps * vpt(i*10 + 300, s1) / (wps * vpt(i*10 + 300, s1) + 1)), svmax)
		
	# find coordinates
	xvalues = np.empty(41)
	yvalues = np.empty(41)
	labels = np.empty(41)
	for i in range(0, 41):
		le = (wpl * vpt(i*10 + 300, l1) / (wpl * vpt(i*10 + 300, l1) + 1)) / lvmax
		me = (wpm * vpt(i*10 + 300, m1) / (wpm * vpt(i*10 + 300, m1) + 1)) / mvmax
		se = (wps * vpt(i*10 + 300, s1) / (wps * vpt(i*10 + 300, s1) + 1)) / svmax
		yvalues[i] = (me - 0.5*(le + se))
		xvalues[i] = (math.sqrt(3)/2) * (le - se)
		labels[i] = i*10 + 300
	
	plt.plot(xvalues, yvalues, '-k')
	for i in range(0, 8): # every 50nm for less messiness
		plt.plot(xvalues[i*5], yvalues[i*5], 'ok')
		plt.text(xvalues[i*5], yvalues[i*5], labels[i*5])
	
	wle = wle / lvmax
	wme = wme / mvmax
	wse = wse / svmax
	wx = (math.sqrt(3)/2) * (wle - wse)
	wy = wme - 0.5*(wle + wse)
	plt.plot(wx, wy, 'ok')
	plt.text(wx, wy, args.white)
	
	xborder = np.array([0, -math.sqrt(3)/2, -math.sqrt(3)/2, 0, math.sqrt(3)/2, math.sqrt(3)/2, 0])
	yborder = np.array([-1, -0.5, 0.5, 1, 0.5, -0.5, -1])
	plt.plot(xborder, yborder, '-k')
	plt.text(-math.sqrt(3)/2 - 0.05, -0.5 - 0.025, 'S')
	plt.text(0 - 0.025, 1 + 0.0125, 'M')
	plt.text(math.sqrt(3)/2 + 0.0125, -0.5 - 0.025, 'L')
	plt.show()
	
# wavelength discrimination take three (http://www.eeza.csic.es/Pollination_ecology/Publicaciones/Telles%20et%20al%202016%20JEB.pdf)
# This is mostly right but for honeybees and butterflies gives an extra minimum
# where a maximum should be (really an extra peak in the first trough) and way too
# high values on either side of it. I cut the values off at 500 because otherwise
# the graph is unreadable.
# It looks much better if I include chromatic adaptation and use the "i" light. Also
# better if the white values are divided by 100?
# No matter what I do to this thing, the graph it produces doesn't make sense. I think we
# have to give up for now.
# With quantum noise added it more resembles what we see in the literature and what we
# would expect from Maxwell triangle distances, but only if the Weber fraction is set
# very low (0.001-0.01 -- even lower values don't make much difference). With values
# higher than 0.01, the minima are strongly shifted toward longer wavelengths. I don't
# know why this happens. Also, for humans the values make much more sense if L, M and S
# are set to 565, 540 and 440 (the peaks of the standard observer functions) instead of
# the "real" values of 560, 530 and 420 because humans have particularly strong
# short-wavelength filters that shift the apparent peak sensitivity to much longer
# wavelengths. I can't use the setting --filter human because it produces too small
# numbers somewhere and the program crashes.
# Some other examples:
# * https://www.researchgate.net/figure/Wavelength-discrimination-A-Discrimination-functions-from-3-participants-One-curve-is_fig2_24218078
if (args.wd):
	# background light
	wpl = 0
	wpm = 0
	wps = 0
	
	for i in range(wp_1nm.shape[0]):
		w = i + 300
		wpl += vpt(w, l1)*lens_filter(w) * wp_1nm[i]
		wpm += vpt(w, m1)*lens_filter(w) * wp_1nm[i]
		wps += vpt(w, s1)*lens_filter(w) * wp_1nm[i]
	kl = 1
	km = wpl / wpm
	ks = wpl / wps
	
	# minima
	min1 = 300
	min1y = 100
	min2 = 300
	min2y = 100
	max1 = 650
	max1y = 0
	
	xvalues = np.empty(300)
	yvalues = np.empty(300)
	for i in range(350, 650):
		xvalues[i-350] = i
		
		# derivatives -- pretty sure the lens filtering cancels
		#scale1 = (vpt(i+1, l1) + vpt(i+1, m1) + vpt(i+1, s1))
#		scale1 = 1
#		dl1 = math.log(1 + kl*vpt(i+1, l1)) / scale1
#		dm1 = math.log(1 + km*vpt(i+1, m1)) / scale1
#		ds1 = math.log(1 + ks*vpt(i+1, s1)) / scale1
#		#scale2 = (vpt(i-1, l1) + vpt(i-1, m1) + vpt(i-1, s1))
#		scale2 = 1
#		dl2 = math.log(1 + kl*vpt(i-1, l1)) / scale2
#		dm2 = math.log(1 + km*vpt(i-1, m1)) / scale2
#		ds2 = math.log(1 + ks*vpt(i-1, s1)) / scale2
#		# Not totally sure if this is supposed to be logarithmic or linear. The
#		# contrast functions are logarithmic.
#		dll = (dl1 - dl2) / 2
#		dlm = (dm1 - dm2) / 2
#		dls = (ds1 - ds2) / 2
#		dfl = (kl / (1 + kl*vpt(i, l1))) * dll
#		dfm = (km / (1 + km*vpt(i, m1))) * dlm
#		dfs = (ks / (1 + ks*vpt(i, s1))) * dls
		
#		v = math.sqrt(((wl*wm)**2 + (wl*ws)**2 + (wm*ws)**2) / (wl**2*(dfm - dfs)**2 + wm**2*(dfl - dfs)**2 + ws**2*(dfm - dfl)**2))
		# intensity?
		# No.
		#intensity = 1
		#v = 0
		#while (intensity < 100):
		#	dl1i = intensity*dl1
		#	dm1i = intensity*dm1
		#	ds1i = intensity*ds1
		#	dll = (dl1 - dl2) / 2
		#	dlm = (dm1 - dm2) / 2
		#	dls = (ds1 - ds2) / 2
		#	dfl = (kl / (1 + kl*vpt(i, l1))) * dll
		#	dfm = (km / (1 + km*vpt(i, m1))) * dlm
		#	dfs = (ks / (1 + ks*vpt(i, s1))) * dls
		#	v2 = math.sqrt(((wl*wm)**2 + (wl*ws)**2 + (wm*ws)**2) / (wl**2*(dfm - dfs)**2 + wm**2*(dfl - dfs)**2 + ws**2*(dfm - dfl)**2))
		#	print(v)
		#	v = max(v, v2)
		#	print(intensity)
		#	intensity += 1
		
		# outsourcing this to color_contrast()
		minus1 = np.zeros(501)
		plus1 = np.zeros(501)
		minus1[i-301] = 100
		plus1[i-299] = 100
		v = color_contrast(minus1, plus1)
		
		#if (v > 100):
		#	yvalues[i-350] = float('inf')
		#else:
		yvalues[i-350] = v
		
		# set minima
		if (i < args.mw and 2/v < min1y and i < max1):
			min1 = i
			min1y = 2/v
		if (i > (args.mw + args.sw) / 2 and i < args.lw and 2/v > max1y):
			max1 = i
			max1y = 2/v
		if (i > args.mw and 2/v < min2y):
			min2 = i
			min2y = 2/v
	
	print("First minimum: " + str(min1) + " nm (" + str(min1y) + ")")
	print("Maximum: " + str(max1) + " nm (" + str(max1y) + ")")
	print("Second minimum: " + str(min2) + " nm (" + str(min2y) + ")")
	plt.plot(xvalues, yvalues/2, 'k')
	plt.xlabel("Wavelength (nm)")
	plt.ylabel('Contrast (JND/nm)')
	plt.show()
	plt.plot(xvalues, 2/yvalues, 'k')
	plt.xlabel("Wavelength (nm)")
	plt.ylabel('Δλ (nm)')
	plt.ylim(0,50) # we can do this, it's pretty cool
	plt.show()
	
# example illuminants
# I tried forcing the "light source" as E to prevent multiplying these by another illuminant like they're
# reflectance spectra, but this doesn't let you add a von Kries transform based on an illuminant
# other than E. Just set --white e --novk.
# 11/03 -- We can now get the desired behavior with the option reflect=False.

if (args.lighting):
	print("White (E, equal-energy)")
	spectral_rendering(e, reflect=False)
	
	print("White (D65)")
	spectral_rendering(d65, reflect=False)
	
	print("Incandescent lighting (A)")
	spectral_rendering(a, reflect=False)
	
	print("Incandescent lighting (approximated with 2856 K blackbody spectrum)")
	spectral_rendering(incandescent, reflect=False)
	
	# color/brightness contrast test
	print("Color contrast between E and E: " + str(color_contrast(e, e)))
	print("Color contrast between E and D65: " + str(color_contrast(e, d65)))
	print("Color contrast between E and A: " + str(color_contrast(e, a)))
	print("Color contrast between E and incandescent: " + str(color_contrast(e, incandescent)))
	print("Color contrast between D65 and D65: " + str(color_contrast(d65, d65)))
	print("Color contrast between D65 and A: " + str(color_contrast(d65, a)))
	print("Color contrast between D65 and incandescent: " + str(color_contrast(d65, incandescent)))
	print("Color contrast between A and A: " + str(color_contrast(a, a)))
	print("Color contrast between A and incandescent: " + str(color_contrast(a, incandescent)))
	print("Color contrast between incandescent and incandescent: " + str(color_contrast(incandescent, incandescent)))
	print("Brightness contrast between E and E: " + str(brightness_contrast(e, e)))
	print("Brightness contrast between E and D65: " + str(brightness_contrast(e, d65)))
	print("Brightness contrast between E and A: " + str(brightness_contrast(e, a)))
	print("Brightness contrast between E and incandescent: " + str(brightness_contrast(e, incandescent)))
	print("Brightness contrast between D65 and D65: " + str(brightness_contrast(d65, d65)))
	print("Brightness contrast between D65 and A: " + str(brightness_contrast(d65, a)))
	print("Brightness contrast between D65 and incandescent: " + str(brightness_contrast(d65, incandescent)))
	print("Brightness contrast between A and A: " + str(brightness_contrast(a, a)))
	print("Brightness contrast between A and incandescent: " + str(brightness_contrast(a, incandescent)))
	print("Brightness contrast between incandescent and incandescent: " + str(brightness_contrast(incandescent, incandescent)))

# "long-pass" step functions: 0 for wavelengths less than the step and 1 for longer wavelengths.
# Many white, yellow, orange and red flowers have reflectance spectra like these, and they're
# also found in paints and camera filters.
if (args.step):
	# 350-nm step
	step350 = np.empty(401)
	for i in range(49): step350[i] = 0
	for i in range(50, 401): step350[i] = 1
	print("350-nm step function")
	spectral_rendering(step350)
	
	# 400-nm step
	step400 = np.empty(401)
	for i in range(99): step400[i] = 0
	for i in range(100, 401): step400[i] = 1
	print("400-nm step function")
	spectral_rendering(step400)
	
	# 450-nm step
	step450 = np.empty(401)
	for i in range(149): step450[i] = 0
	for i in range(150, 401): step450[i] = 1
	print("450-nm step function")
	spectral_rendering(step450)
	
	# 500-nm step
	step500 = np.empty(401)
	for i in range(199): step500[i] = 0
	for i in range(200, 401): step500[i] = 1
	print("500-nm step function")
	spectral_rendering(step500)
	
	# 550-nm step
	step550 = np.empty(401)
	for i in range(249): step550[i] = 0
	for i in range(250, 401): step550[i] = 1
	print("550-nm step function")
	spectral_rendering(step550)
	
	# 600-nm step
	step600 = np.empty(401)
	for i in range(299): step600[i] = 0
	for i in range(299, 401): step600[i] = 1
	print("600-nm step function")
	spectral_rendering(step600)
	
	# 650-nm step
	step650 = np.empty(401)
	for i in range(349): step650[i] = 0
	for i in range(350, 401): step650[i] = 1
	print("650-nm step function")
	spectral_rendering(step650)

# "short-pass" step functions, the reverse of the above. These form the other edge of the
# optimal color solid (see Wikipedia) but don't seem to occur in nature or have any
# applications. Blue/green surfaces usually are the "band-pass" type with one or two
# peaks.
if (args.steps):
	# 350-nm step
	step350 = np.empty(401)
	for i in range(49): step350[i] = 1
	for i in range(50, 401): step350[i] = 0
	print("350-nm step function")
	spectral_rendering(step350)
	
	# 400-nm step
	step400 = np.empty(401)
	for i in range(99): step400[i] = 1
	for i in range(100, 401): step400[i] = 0
	print("400-nm step function")
	spectral_rendering(step400)
	
	# 450-nm step
	step450 = np.empty(401)
	for i in range(149): step450[i] = 1
	for i in range(150, 401): step450[i] = 0
	print("450-nm step function")
	spectral_rendering(step450)
	
	# 500-nm step
	step500 = np.empty(401)
	for i in range(199): step500[i] = 1
	for i in range(200, 401): step500[i] = 0
	print("500-nm step function")
	spectral_rendering(step500)
	
	# 550-nm step
	step550 = np.empty(401)
	for i in range(249): step550[i] = 1
	for i in range(250, 401): step550[i] = 0
	print("550-nm step function")
	spectral_rendering(step550)
	
	# 600-nm step
	step600 = np.empty(401)
	for i in range(299): step600[i] = 1
	for i in range(300, 401): step600[i] = 0
	print("600-nm step function")
	spectral_rendering(step600)
	
	# 650-nm step
	step650 = np.empty(401)
	for i in range(349): step650[i] = 1
	for i in range(350, 401): step650[i] = 0
	print("650-nm step function")
	spectral_rendering(step650)

# peaked or "band-pass" spectra. Blue, green and purple objects usually have this type
# of reflectance, including green leaves and some flowers. Modeled with a normal distribution.
# See https://mathworld.wolfram.com/NormalDistribution.html
# The peak region is usually set to about 100 nm wide by choosing 50 as the standard
# deviation. The width at half-maximum is slightly narrower.
if (args.peak):
	# 350-nm peak
	peak350 = np.empty(401)
	for i in range(401):
		peak350[i] = math.exp(-(i+300-350)**2/(2*25)**2)
	print("350-nm peak")
	spectral_rendering(peak350)
	
	# 400-nm peak
	peak400 = np.empty(401)
	for i in range(401):
		peak400[i] = math.exp(-(i+300-400)**2/(2*25)**2)
	print("400-nm peak")
	spectral_rendering(peak400)
	
	# 450-nm peak
	peak450 = np.empty(401)
	for i in range(401):
		peak450[i] = math.exp(-(i+300-450)**2/(2*25)**2)
	print("450-nm peak")
	spectral_rendering(peak450)
	
	# 500-nm peak
	peak500 = np.empty(401)
	for i in range(401):
		peak500[i] = math.exp(-(i+300-500)**2/(2*25)**2)
	print("500-nm peak")
	spectral_rendering(peak500)
	
	# 550-nm peak. This is similar to the typical reflectance peak of green
	# vegetation created by the combination of chlorophyll (absorbs red and blue)
	# and carotenoids (absorbs blue and blue-green), except the long-wavelength side of
	# the slope is slightly shallower and is more so for lighter/yellower leaves.
	# This color would be a pure, idealized "leaf-green".
	peak550 = np.empty(401)
	for i in range(401):
		peak550[i] = math.exp(-(i+300-550)**2/(2*25)**2)
	print("550-nm peak")
	spectral_rendering(peak550)
	
	# 600-nm peak
	peak600 = np.empty(401)
	for i in range(401):
		peak600[i] = math.exp(-(i+300-600)**2/(2*25)**2)
	print("600-nm peak")
	spectral_rendering(peak600)
	
	# 650-nm peak
	peak650 = np.empty(401)
	for i in range(401):
		peak650[i] = math.exp(-(i+300-650)**2/(2*25)**2)
	print("650-nm peak")
	spectral_rendering(peak650)
	
	# carotenoid-like: green
	# See https://doi.org/10.1098/rspb.2013.2209. The text says the green peak is 520,
	# but Figure 6 shows it at more like 550. The color of green budgerigar feathers
	# appears yellow-green in photographs, so I'm going with the 550 version. This
	# type of green is a subtractive color formed from structural blue (minus-red,
	# sort of) and carotenoid yellow (minus-blue).
	c_green = np.empty(401)
	for i in range(401):
		c_green[i] = math.exp(-(i+300-320)**2/(2*20)**2) + math.exp(-(i+300-550)**2/(2*25)**2)
	print("Carotenoid-like green")
	spectral_rendering(c_green)
	
	# carotenoid-like: yellow (320 + 500-nm step)
	# For the pure carotenoid colors, the UV peak is about half the height of
	# the long-wavelength step. See http://dx.doi.org/10.1093/jmammal/gyz035.
	c_yellow = np.empty(401)
	for i in range(200):
		c_yellow[i] = 0.5*math.exp(-(i+300-320)**2/(2*20)**2)
	for i in range(200, 401): c_yellow[i] = 1
	print("Carotenoid-like yellow")
	spectral_rendering(c_yellow)
	
	# carotenoid-like: red (320 + 600-nm step)
	c_red = np.empty(401)
	for i in range(300):
		c_red[i] = 0.5*math.exp(-(i+300-320)**2/(2*20)**2)
	for i in range(300, 401): c_red[i] = 1
	print("Carotenoid-like red")
	spectral_rendering(c_red)
	
	# flower-like: purple (450 + 600-nm step)
	purple = np.empty(401)
	for i in range(300):
		purple[i] = math.exp(-(i+300-450)**2/(2*25)**2)
	for i in range(300, 401): purple[i] = 1
	print("Flower-like purple")
	spectral_rendering(purple)
	
	# leaf-like: 550-nm peak + 1/4-height flat reflectance
	leaf1 = np.empty(401)
	for i in range(401):
		leaf1[i] = 0.4*math.exp(-(i+300-550)**2/(2*25)**2) + 0.1
	print("Leaf-like 1 (diffuse/transmitted)")
	spectral_rendering(leaf1)
	
	# leaf with equal-height UV peak
	leaf2 = np.empty(401)
	for i in range(401):
		leaf2[i] = 0.4*math.exp(-(i+300-550)**2/(2*25)**2) + 0.4*math.exp(-(i+300-350)**2/(2*25)**2) + 0.1
	print("Leaf-like 2 (specular)")
	spectral_rendering(leaf2)
	
	# leaf with half-height UV peak
	leaf3 = np.empty(401)
	for i in range(401):
		leaf3[i] = 0.4*math.exp(-(i+300-550)**2/(2*25)**2) + 0.2*math.exp(-(i+300-350)**2/(2*25)**2) + 0.1
	print("Leaf-like 3 (specular, less UV)")
	spectral_rendering(leaf3)
	
	# wider 550-nm peak
	peak550w = np.empty(401)
	for i in range(401):
		peak550w[i] = 0.5*math.exp(-(i+300-550)**2/(2*50)**2)
	print("Leaf-like 4 (wider peak)")
	spectral_rendering(peak550w)
	
	#x = np.empty(401)
	#for i in range(401): x[i] = i+300
	#plt.plot(x, c_green, 'g')
	#plt.plot(x, c_yellow, 'y')
	#plt.plot(x, c_red, 'r')
	#plt.plot(x, leaf1, 'k')
	#plt.plot(x, leaf2, 'k')
	#plt.show()

# ramp functions
# Soil, dead/dying vegetation and the fur of some Arctic mammals has ramp-like
# reflectance with an increasing slope, causing a brown/yellow appearance. Snow
# is a negative, shallower ramp.
# Examples:
# * soil and vegetation: https://www.researchgate.net/profile/Yingying-Yang-5/publication/340281560/figure/fig5/AS:878006243430404@1586344409333/a-A-comparison-of-spectral-reflectance-curves-between-soil-and-withered-vegetation.png
# * different soil types: https://www.researchgate.net/publication/47426815/figure/fig1/AS:11431281212212810@1702556618056/Typical-spectral-reflectance-of-vegetation-soils-water-and-snow-The-original-data.tif
# * fur and snow: https://www.mdpi.com/2072-4292/8/4/273
# * red vs. green leaves: http://2.bp.blogspot.com/-DfLsYNxVoMo/UbwHSrq_KxI/AAAAAAAAAYM/OZyOnBPRxf4/s1600/LeafSpectra.jpg The color of red/yellow leaves is mainly due to
# carotenoids, and carotenoid spectra seem to be somewhere between a step and a ramp.
# Some mammals with white fur (e.g. polar bear) absorb UV whereas others (e.g. Arctic
# hare) do not. An example of non-absorbing white fur is the Honduran white bat (see
# carotenoid example).
if (args.ramp):
	# 300-nm ramp
	ramp300 = np.empty(401)
	for i in range(401): ramp300[i] = i/401
	print("300-nm ramp function")
	spectral_rendering(ramp300)
	
	# 350-nm ramp
	ramp350 = np.empty(401)
	for i in range(49): ramp350[i] = 0
	for i in range(50, 401): ramp350[i] = (i-50)/351
	print("350-nm ramp function")
	spectral_rendering(ramp350)
	
	# 400-nm ramp
	ramp400 = np.empty(401)
	for i in range(99): ramp400[i] = 0
	for i in range(100, 401): ramp400[i] = (i-100)/301
	print("400-nm ramp function")
	spectral_rendering(ramp400)
	
	# 450-nm ramp
	ramp450 = np.empty(401)
	for i in range(149): ramp450[i] = 0
	for i in range(150, 401): ramp450[i] = (i-150)/251
	print("450-nm ramp function")
	spectral_rendering(ramp450)
	
	# 500-nm ramp
	ramp500 = np.empty(401)
	for i in range(199): ramp500[i] = 0
	for i in range(200, 401): ramp500[i] = (i-200)/201
	print("500-nm ramp function")
	spectral_rendering(ramp500)
	
	# 550-nm ramp
	ramp550 = np.empty(401)
	for i in range(249): ramp550[i] = 0
	for i in range(250, 401): ramp550[i] = (i-250)/151
	print("550-nm ramp function")
	spectral_rendering(ramp550)
	
	# 600-nm ramp
	ramp600 = np.empty(401)
	for i in range(299): ramp600[i] = 0
	for i in range(300, 401): ramp600[i] = (i-300)/101
	print("600-nm ramp function")
	spectral_rendering(ramp600)
	
	# snow-like: 300-nm reversed ramp + flat reflectance
	ramp_snow = np.empty(401)
	for i in range(401): ramp_snow[i] = 1 - 0.0004*i
	print("Snow-like ramp function")
	spectral_rendering(ramp_snow)
	
	# 300-nm ramp with custom slope (beige, fur-like)
	ramp300half = np.empty(401)
	for i in range(401): ramp300half[i] = 0.0015*i
	print("Fur-like ramp function")
	spectral_rendering(ramp300half)
	
	# 300-nm quarter-height ramp (dull brown, soil-like)
	ramp300q = np.empty(401)
	for i in range(401): ramp300q[i] = (i/401)/4
	print("300-nm ramp function (quarter height)")
	spectral_rendering(ramp300q)
	
	# 500-nm quarter-height ramp (brown-orange, clay soil-like)
	ramp500q = np.empty(401)
	for i in range(199): ramp500q[i] = 0
	for i in range(200, 401): ramp500q[i] = ((i-200)/201)/4
	print("500-nm ramp function (quarter height)")
	spectral_rendering(ramp500q)

# sky and daylight
# These are both normalized to "1" at 550 nm. If I try to make "sky" relative to "daylight",
# the values are way too low and probably not meaningful.
if (args.sky):
	daylight_table = np.empty(41)
	for i in range(41):
		w = i*10 + 300
		# normalize to 100% at 550
		daylight_table[i] = 100*blackbody(w, 5800) / blackbody(550, 5800)
	
	print("Black body spectrum approximating daylight")
	spectral_rendering(daylight_table, reflect=False)
	
	sky_table = np.empty(41)
	for i in range(0, 41):
		w = i*10 + 300
		sky_table[i] = 100*blackbody(w, 5800) * w**-4 / (blackbody(550, 5800) * 550**-4)
		
	print("Black body spectrum with Rayleigh scattering approximating a blue sky")
	spectral_rendering(sky_table, reflect=False)
	
	sky1_table = np.empty(41)
	for i in range(0, 41):
		w = i*10 + 300
		sky1_table[i] = 100*d65[i] * w**-4 * 5e8
	
	print("D65 with Rayleigh scattering")
	spectral_rendering(sky1_table, reflect=False)

# Kodak Wratten camera filters
# These were obtained from the graphs at https://www.kodak.com/en/motion/page/wratten-2-filters/
# using Plot Digitizer (https://plotdigitizer.sourceforge.net/). The wavelength intervals are
# roughly 5 nm (getting exact integer X values usually isn't possible).

# yellow 15
yellow15_5nm = np.array([[
	300,
	301.510,
	305.402,
	310.343,
	314.208,
	317.014,
	320.513,
	325.750,
	329.928,
	335.146,
	336.536,
	337,
	501,
	502.211,
	505.413,
	510.383,
	515.688,
	519.566,
	525.542,
	530.806,
	535.015,
	540.273,
	545.178,
	550.433,
	555.688,
	560.241,
	565.146,
	570.399,
	575.303,
	579.857,
	585.111,
	590.365,
	595.619,
	599.822,
	605.426,
	609.279,
	615.233,
	619.437,
	625.040,
	629.944,
	635.549,
	639.752,
	644.656,
	650.259,
	655.513,
	659.367,
	665.671,
	669.875,
	675.128,
	680.032,
	685.286,
	690.190,
	695.793,
	700.347,
	705.251,
	710.504,
	715.409,
	720.312,
	725.216,
	729.768,
	734.673,
	739.927,
	744.830,
	750.785,
	756.039,
	760.592,
	765.495,
	770.750,
	774.603,
	779.857,
	786.512,
	791.416,
	795.619,
	799.822,
	805.776,
	810.680,
	815.934,
	820.137,
	825.041,
	830.295,
	835.199,
	840.102,
	844.305,
	849.559,
	854.463,
	859.717,
	864.971,
	870.224,
	874.428,
	880.032,
	884.936,
	890.189,
	895.443,
	900.347,
],[
	0,
	0.00100584,
	0.00311012,
	0.00896471,
	0.0128078,
	0.0141464,
	0.0127323,
	0.00792732,
	0.00386065,
	0.00137094,
	0.00101140,
	0,
	0,
	0.00101569,
	0.00425698,
	0.0279919,
	0.122221,
	0.245141,
	0.463741,
	0.621248,
	0.727500,
	0.803516,
	0.856872,
	0.877112,
	0.897829,
	0.897789,
	0.908310,
	0.908263,
	0.908220,
	0.913507,
	0.918819,
	0.913413,
	0.913366,
	0.918687,
	0.918637,
	0.923991,
	0.913192,
	0.923899,
	0.907779,
	0.918417,
	0.918367,
	0.918329,
	0.923672,
	0.907556,
	0.912833,
	0.923540,
	0.907420,
	0.923445,
	0.918013,
	0.907293,
	0.923307,
	0.917878,
	0.907153,
	0.923171,
	0.923127,
	0.907023,
	0.923035,
	0.917608,
	0.912213,
	0.896306,
	0.922862,
	0.922815,
	0.912039,
	0.922717,
	0.922670,
	0.917248,
	0.901250,
	0.917157,
	0.922503,
	0.911728,
	0.917016,
	0.927763,
	0.922314,
	0.916897,
	0.916844,
	0.922179,
	0.922131,
	0.911370,
	0.911326,
	0.916625,
	0.921958,
	0.921914,
	0.905841,
	0.905795,
	0.916409,
	0.916362,
	0.916315,
	0.905612,
	0.916231,
	0.921555,
	0.916137,
	0.910747,
	0.905390,
	0.915999,
]])

# red 25
red25_5nm = np.array([[
	299.592,
	304.516,
	310.143,
	314.722,
	320.361,
	325.306,
	330.239,
	334.810,
	339.713,
	344.954,
	349.841,
	354.723,
	359.945,
	365.152,
	370.356,
	374.853,
	376.932,
	377,
	573,
	573.600,
	575.761,
	580.171,
	584.921,
	589.604,
	594.587,
	600.227,
	605.150,
	610.415,
	615.324,
	620.234,
	625.493,
	630.751,
	635.310,
	640.219,
	645.476,
	650.385,
	654.592,
	660.201,
	664.758,
	669.667,
	674.576,
	680.534,
	684.741,
	689.650,
	694.908,
	700.517,
	705.424,
	710.685,
	715.241,
	720.851,
	725.407,
	729.615,
	735.224,
	740.132,
	745.740,
	750.298,
	755.207,
	760.465,
	765.022,
	769.929,
	774.837,
	780.448,
	785.005,
	789.913,
	795.521,
	800.430,
	805.338,
	810.598,
	815.154,
	820.412,
	825.320,
	830.229,
	835.137,
	840.396,
	845.302,
	850.210,
	855.118,
	860.027,
	864.584,
	869.843,
	875.802,
	880.360,
	885.267,
	890.526,
	895.435,
	900.342,
],[
	0.00372910,
	0.00437114,
	0.00518432,
	0.00637012,
	0.00844965,
	0.0121000,
	0.0154022,
	0.0175304,
	0.0166230,
	0.0140939,
	0.0115349,
	0.00895313,
	0.00632429,
	0.00387854,
	0.00229604,
	0.00128908,
	0.00101250,
	0,
	0,
	0.00100074,
	0.00174062,
	0.0121518,
	0.0767534,
	0.256647,
	0.526375,
	0.710655,
	0.818422,
	0.867931,
	0.878087,
	0.893609,
	0.904056,
	0.898609,
	0.909133,
	0.914371,
	0.903525,
	0.914097,
	0.919382,
	0.913833,
	0.908345,
	0.918974,
	0.924268,
	0.913286,
	0.913173,
	0.923858,
	0.923715,
	0.918140,
	0.912617,
	0.928739,
	0.923162,
	0.923010,
	0.912080,
	0.922772,
	0.928069,
	0.922486,
	0.911534,
	0.922210,
	0.933001,
	0.927379,
	0.916398,
	0.910885,
	0.910753,
	0.926834,
	0.926709,
	0.921135,
	0.915575,
	0.926288,
	0.926155,
	0.936982,
	0.925887,
	0.920308,
	0.920175,
	0.925476,
	0.925342,
	0.930664,
	0.919634,
	0.919501,
	0.919368,
	0.924664,
	0.924540,
	0.924397,
	0.924235,
	0.929569,
	0.918552,
	0.918410,
	0.929157,
	0.918145,
]])

# blue 47
blue47_5nm = np.array([[
	299.664,
	304.912,
	310.505,
	315.063,
	319.985,
	324.907,
	329.819,
	335.070,
	340.319,
	345.925,
	350.488,
	355.408,
	359.630,
	365.257,
	370.180,
	375.103,
	380.376,
	384.946,
	389.515,
	395.140,
	400.416,
	405.340,
	410.261,
	415.181,
	420.445,
	425.003,
	430.261,
	435.165,
	440.419,
	444.971,
	449.872,
	455.473,
	460.373,
	465.272,
	470.519,
	475.766,
	480.312,
	485.555,
	489.746,
	494.636,
	499.870,
	505.101,
	509.979,
	514.850,
	520.415,
	525.274,
	529.783,
	531.164,
	532,
	689,
	689.833,
	694.770,
	700.408,
	704.986,
	709.915,
	715.190,
	720.468,
	725.392,
	730.319,
	734.540,
	740.520,
	745.440,
	750.362,
	754.929,
	760.197,
	764.759,
	770.022,
	774.582,
	780.192,
	785.449,
	790.354,
	794.910,
	800.165,
	805.071,
	810.676,
	815.579,
	820.483,
	825.387,
	830.292,
	835.197,
	840.100,
	845.003,
	850.257,
	855.512,
	860.065,
	864.969,
	869.873,
	875.477,
	880.380,
	885.284,
	890.188,
	895.443,
	900.347,
],[
	0.00155619,
	0.00140765,
	0.00113164,
	0.00123621,
	0.00175044,
	0.00250799,
	0.00290613,
	0.00273951,
	0.00249268,
	0.00259751,
	0.00311826,
	0.00423686,
	0.00614255,
	0.00944611,
	0.0138572,
	0.0199717,
	0.0287840,
	0.0393410,
	0.0534539,
	0.0793447,
	0.122018,
	0.181121,
	0.253456,
	0.344377,
	0.415860,
	0.459676,
	0.490441,
	0.499158,
	0.493268,
	0.484586,
	0.459507,
	0.435722,
	0.403541,
	0.369355,
	0.322486,
	0.279910,
	0.245839,
	0.198805,
	0.157025,
	0.119712,
	0.0815928,
	0.0524272,
	0.0321349,
	0.0173005,
	0.00813278,
	0.00341796,
	0.00147076,
	0.000990734,
	0,
	0,
	0.000994265,
	0.00190178,
	0.00363760,
	0.00576142,
	0.00939824,
	0.0141159,
	0.0224895,
	0.0331866,
	0.0522534,
	0.0726932,
	0.119985,
	0.164960,
	0.233580,
	0.304539,
	0.399399,
	0.471065,
	0.552318,
	0.625082,
	0.707422,
	0.745920,
	0.768184,
	0.805235,
	0.829264,
	0.859066,
	0.869185,
	0.869122,
	0.869059,
	0.879305,
	0.894933,
	0.910839,
	0.910773,
	0.900031,
	0.899961,
	0.910567,
	0.915890,
	0.910440,
	0.915758,
	0.910299,
	0.910233,
	0.910167,
	0.915483,
	0.926272,
	0.931682,
]])

# green 58
green58_5nm = np.array([[
	300,
	462,
	462.347,
	465.499,
	469.702,
	480.210,
	485.464,
	489.667,
	495.447,
	499.650,
	505.429,
	510.158,
	515.412,
	520.140,
	524.869,
	530.123,
	535.026,
	540.280,
	544.834,
	550.438,
	555.692,
	560.245,
	565.149,
	570.403,
	575.306,
	580.210,
	585.814,
	590.018,
	594.571,
	600.525,
	605.429,
	609.982,
	614.536,
	617.338,
	618,
	702,
	702.102,
	705.254,
	710.158,
	715.061,
	719.965,
	725.219,
	730.823,
	735.026,
	739.930,
	745.884,
	750.088,
	754.991,
	760.245,
	765.149,
	770.403,
	774.956,
	780.210,
	785.814,
	790.368,
	795.271,
	800.175,
	805.779,
	810.683,
	815.236,
	820.490,
	825.744,
	829.947,
	835.201,
	840.105,
	845.359,
	850.963,
	855.166,
	860.070,
	864.623,
	870.928,
	875.832,
	880.385,
	884.939,
	890.543,
	895.447,
	900,
],[
	0,
	0,
	0.00100884,
	0.00236877,
	0.00634670,
	0.0331909,
	0.0614519,
	0.110812,
	0.199819,
	0.294300,
	0.407560,
	0.485987,
	0.540112,
	0.549702,
	0.544886,
	0.521432,
	0.494615,
	0.455616,
	0.409958,
	0.360320,
	0.309347,
	0.256400,
	0.207586,
	0.155724,
	0.112120,
	0.0774769,
	0.0436004,
	0.0280807,
	0.0167573,
	0.00843557,
	0.00469178,
	0.00278347,
	0.00145993,
	0.00101477,
	0,
	0,
	0.00101477,
	0.00187882,
	0.00494615,
	0.0125708,
	0.0287474,
	0.0591527,
	0.108879,
	0.166105,
	0.248988,
	0.339790,
	0.414796,
	0.494615,
	0.576115,
	0.636534,
	0.678965,
	0.719987,
	0.754583,
	0.781616,
	0.800175,
	0.814382,
	0.828841,
	0.838623,
	0.843557,
	0.853513,
	0.853513,
	0.858535,
	0.863586,
	0.863586,
	0.863586,
	0.868667,
	0.868667,
	0.868667,
	0.878919,
	0.878919,
	0.878919,
	0.878919,
	0.878919,
	0.878919,
	0.884091,
	0.873778,
	0.878919,
]])

# neutral density 0.1
nd01_5nm = np.array([[
	299.973,
	304.881,
	310.490,
	315.398,
	319.956,
	324.864,
	330.122,
	335.029,
	339.586,
	344.494,
	350.102,
	355.009,
	359.916,
	365.174,
	370.081,
	375.338,
	380.245,
	385.502,
	390.760,
	395.316,
	400.223,
	404.779,
	410.036,
	415.644,
	420.200,
	424.757,
	430.715,
	434.920,
	440.178,
	445.435,
	450.692,
	455.949,
	459.804,
	465.061,
	469.968,
	475.225,
	480.482,
	485.389,
	490.295,
	495.202,
	499.758,
	505.365,
	509.922,
	514.828,
	520.436,
	525.342,
	530.599,
	535.155,
	540.412,
	545.319,
	550.226,
	555.483,
	560.039,
	565.296,
	570.553,
	575.459,
	580.016,
	585.623,
	590.529,
	595.436,
	599.992,
	605.600,
	609.455,
	614.712,
	619.618,
	625.226,
	630.483,
	635.389,
	639.945,
	645.202,
	650.459,
	655.716,
	660.272,
	665.880,
	670.436,
	675.342,
	680.249,
	685.506,
	690.763,
	695.319,
	699.525,
	705.132,
	710.740,
	715.297,
	720.554,
	725.811,
	730.367,
	735.624,
	740.181,
	745.788,
	750.344,
	755.251,
	760.157,
	765.414,
	770.321,
	774.877,
	780.485,
	785.742,
	790.298,
	794.854,
	800.462,
	805.368,
	810.625,
	815.181,
	820.438,
	825.695,
	830.251,
	834.807,
	840.064,
	844.971,
	850.228,
	855.485,
	860.392,
	864.948,
	870.205,
	875.111,
	880.368,
	885.625,
	890.532,
	895.088,
	900.345,
],[
	0.462957,
	0.479427,
	0.505246,
	0.532459,
	0.561140,
	0.577722,
	0.594793,
	0.608808,
	0.623155,
	0.641571,
	0.652863,
	0.668246,
	0.680013,
	0.687959,
	0.700073,
	0.704133,
	0.712363,
	0.716495,
	0.724867,
	0.733343,
	0.737598,
	0.737566,
	0.741844,
	0.750510,
	0.750476,
	0.754834,
	0.759207,
	0.759176,
	0.768047,
	0.772501,
	0.767968,
	0.776942,
	0.776913,
	0.776873,
	0.781382,
	0.781342,
	0.781302,
	0.781265,
	0.781227,
	0.785761,
	0.785726,
	0.781113,
	0.790246,
	0.790208,
	0.790165,
	0.790127,
	0.790087,
	0.790052,
	0.794634,
	0.794596,
	0.794558,
	0.794518,
	0.799131,
	0.794442,
	0.794401,
	0.799011,
	0.798976,
	0.798932,
	0.794247,
	0.794209,
	0.798821,
	0.798777,
	0.798747,
	0.798706,
	0.794022,
	0.793979,
	0.793938,
	0.793900,
	0.793865,
	0.793825,
	0.793784,
	0.793744,
	0.793708,
	0.793665,
	0.793630,
	0.793592,
	0.793554,
	0.793514,
	0.793473,
	0.793438,
	0.798048,
	0.802674,
	0.807326,
	0.816766,
	0.816724,
	0.821461,
	0.821425,
	0.826189,
	0.835849,
	0.835803,
	0.835766,
	0.835726,
	0.840576,
	0.835643,
	0.835604,
	0.845374,
	0.845327,
	0.845284,
	0.845247,
	0.845209,
	0.850109,
	0.850068,
	0.850025,
	0.849987,
	0.849943,
	0.849900,
	0.849862,
	0.844881,
	0.844838,
	0.844798,
	0.849697,
	0.849654,
	0.854585,
	0.849576,
	0.859503,
	0.854462,
	0.854419,
	0.849405,
	0.854334,
	0.849327,
	0.849283,
]])

# 0.2
nd02_5nm = np.array([[
	299.825,
	305.439,
	310.351,
	314.561,
	319.474,
	325.439,
	330,
	335.263,
	339.825,
	345.088,
	350,
	355.263,
	359.825,
	365.088,
	370.351,
	374.561,
	380.175,
	385.439,
	390,
	395.263,
	399.825,
	404.386,
	410.351,
	415.263,
	419.825,
	425.439,
	430.351,
	435.263,
	439.825,
	445.439,
	450,
	455.263,
	460.175,
	465.088,
	470,
	474.912,
	479.825,
	485.088,
	490,
	494.912,
	499.825,
	505.088,
	510,
	515.263,
	519.825,
	525.088,
	530,
	535.263,
	539.825,
	545.088,
	550,
	554.912,
	560.175,
	565.088,
	570,
	574.912,
	580.175,
	585.088,
	590,
	594.912,
	600.175,
	605.088,
	610,
	615.263,
	619.825,
	625.088,
	630,
	634.912,
	640.175,
	645.088,
	649.649,
	654.912,
	660.175,
	664.737,
	670,
	674.912,
	680.175,
	685.088,
	690,
	694.912,
	700.175,
	705.088,
	710,
	714.912,
	720.175,
	725.088,
	730,
	734.912,
	740.175,
	745.088,
	750,
	755.263,
	759.825,
	765.088,
	770,
	774.912,
	780.175,
	785.088,
	790.351,
	794.912,
	800.175,
	805.088,
	810,
	814.912,
	820.175,
	825.088,
	830,
	835.263,
	839.825,
	845.088,
	850,
	855.614,
	860.175,
	864.737,
	870.351,
	875.263,
	880.175,
	885.088,
	890.351,
	894.912,
	900.175,
],[
	0.304441,
	0.319012,
	0.344193,
	0.373538,
	0.407761,
	0.432302,
	0.447729,
	0.458319,
	0.469159,
	0.483070,
	0.500309,
	0.515144,
	0.527328,
	0.536655,
	0.546147,
	0.552568,
	0.555807,
	0.562341,
	0.568953,
	0.575642,
	0.579016,
	0.582409,
	0.592711,
	0.596185,
	0.603194,
	0.606729,
	0.610285,
	0.613862,
	0.621080,
	0.617460,
	0.624720,
	0.628381,
	0.628381,
	0.635769,
	0.635769,
	0.639496,
	0.643244,
	0.643244,
	0.643244,
	0.643244,
	0.643244,
	0.643244,
	0.643244,
	0.647014,
	0.643244,
	0.647014,
	0.643244,
	0.647014,
	0.647014,
	0.647014,
	0.650806,
	0.647014,
	0.647014,
	0.647014,
	0.643244,
	0.647014,
	0.647014,
	0.650806,
	0.654621,
	0.654621,
	0.654621,
	0.658458,
	0.650806,
	0.647014,
	0.647014,
	0.643244,
	0.643244,
	0.643244,
	0.639496,
	0.639496,
	0.639496,
	0.639496,
	0.639496,
	0.635769,
	0.635769,
	0.643244,
	0.635769,
	0.635769,
	0.639496,
	0.647014,
	0.650806,
	0.650806,
	0.662317,
	0.674031,
	0.685953,
	0.694018,
	0.702177,
	0.710433,
	0.722998,
	0.727236,
	0.735786,
	0.740098,
	0.744436,
	0.748799,
	0.748799,
	0.748799,
	0.753188,
	0.766510,
	0.757603,
	0.762043,
	0.753188,
	0.762043,
	0.762043,
	0.771003,
	0.766510,
	0.762043,
	0.762043,
	0.762043,
	0.771003,
	0.771003,
	0.775522,
	0.771003,
	0.766510,
	0.771003,
	0.775522,
	0.775522,
	0.780067,
	0.775522,
	0.775522,
	0.775522,
	0.775522,
]])

# 0.3
nd03_5nm = np.array([[
	299.719,
	304.620,
	309.870,
	314.768,
	320.017,
	325.267,
	329.468,
	335.071,
	339.974,
	344.526,
	349.778,
	355.029,
	359.931,
	364.834,
	370.437,
	375.691,
	380.243,
	385.147,
	390.750,
	394.953,
	400.206,
	404.759,
	409.662,
	414.565,
	420.169,
	425.072,
	429.625,
	434.528,
	440.132,
	445.036,
	449.939,
	455.193,
	460.447,
	464.650,
	470.604,
	475.508,
	480.411,
	485.665,
	490.919,
	494.772,
	499.675,
	505.279,
	510.183,
	515.087,
	520.341,
	525.244,
	529.798,
	534.701,
	539.605,
	545.209,
	550.463,
	555.717,
	559.920,
	564.824,
	570.078,
	574.631,
	579.535,
	584.438,
	590.042,
	594.945,
	599.849,
	605.103,
	610.358,
	614.911,
	619.465,
	625.069,
	630.323,
	635.578,
	639.781,
	645.735,
	650.639,
	655.192,
	659.395,
	664.299,
	670.254,
	675.508,
	680.061,
	684.615,
	690.219,
	695.122,
	700.025,
	704.928,
	710.181,
	715.084,
	719.987,
	724.540,
	730.143,
	735.396,
	740.299,
	745.202,
	750.455,
	755.008,
	759.561,
	764.115,
	769.368,
	774.972,
	780.226,
	784.779,
	790.383,
	795.987,
	800.190,
	805.794,
	809.997,
	815.251,
	820.155,
	825.409,
	830.663,
	834.866,
	839.769,
	844.673,
	850.277,
	855.881,
	860.084,
	865.338,
	870.242,
	875.846,
	880.399,
	885.303,
	890.557,
	895.461,
	900.014,
],[
	0.140340,
	0.148742,
	0.170015,
	0.197747,
	0.228670,
	0.255371,
	0.267531,
	0.278650,
	0.288548,
	0.300538,
	0.318532,
	0.341549,
	0.355742,
	0.368379,
	0.374876,
	0.381486,
	0.388211,
	0.395055,
	0.399694,
	0.406738,
	0.411512,
	0.418766,
	0.423680,
	0.433662,
	0.436212,
	0.443903,
	0.449111,
	0.457029,
	0.459716,
	0.459738,
	0.465133,
	0.467867,
	0.470617,
	0.473379,
	0.476164,
	0.478962,
	0.481775,
	0.484607,
	0.484632,
	0.484650,
	0.484673,
	0.484699,
	0.484722,
	0.484745,
	0.484770,
	0.484793,
	0.484814,
	0.484837,
	0.487685,
	0.490553,
	0.490578,
	0.484936,
	0.484956,
	0.484979,
	0.482194,
	0.485025,
	0.490717,
	0.493599,
	0.499395,
	0.505256,
	0.499443,
	0.499468,
	0.493723,
	0.493745,
	0.488063,
	0.485262,
	0.482475,
	0.476926,
	0.476945,
	0.479752,
	0.476995,
	0.477016,
	0.479815,
	0.477058,
	0.477086,
	0.477110,
	0.477131,
	0.477152,
	0.477178,
	0.479981,
	0.488443,
	0.497055,
	0.508767,
	0.517737,
	0.533022,
	0.545578,
	0.564962,
	0.578274,
	0.591898,
	0.598844,
	0.612955,
	0.620146,
	0.623787,
	0.627449,
	0.631137,
	0.638548,
	0.642301,
	0.642329,
	0.642364,
	0.649907,
	0.649933,
	0.653756,
	0.653782,
	0.653816,
	0.653847,
	0.657690,
	0.661556,
	0.661582,
	0.665469,
	0.665500,
	0.665537,
	0.669451,
	0.665600,
	0.665634,
	0.669544,
	0.669580,
	0.669610,
	0.673543,
	0.673577,
	0.673609,
	0.673639,
]])

# 0.4
nd04_5nm = np.array([[
	300,
	305.254,
	310.158,
	315.412,
	319.615,
	324.518,
	329.422,
	334.676,
	339.930,
	344.834,
	349.387,
	354.291,
	359.895,
	364.799,
	369.002,
	375.306,
	380.210,
	385.814,
	390.718,
	395.271,
	400.525,
	405.079,
	409.632,
	415.587,
	419.790,
	424.694,
	430.648,
	435.552,
	440.455,
	445.709,
	450.263,
	455.517,
	460.070,
	464.623,
	469.877,
	475.482,
	479.685,
	484.588,
	490.543,
	494.746,
	500,
	505.254,
	510.508,
	515.762,
	519.965,
	525.569,
	530.473,
	535.377,
	540.280,
	545.534,
	550.788,
	555.692,
	559.545,
	565.149,
	570.053,
	574.606,
	580.210,
	585.464,
	589.667,
	594.921,
	599.825,
	604.729,
	610.333,
	615.937,
	620.140,
	625.394,
	630.298,
	635.201,
	640.455,
	645.359,
	650.613,
	655.517,
	659.720,
	664.623,
	669.527,
	674.431,
	680.035,
	685.289,
	690.893,
	695.447,
	700.350,
	705.954,
	710.508,
	714.711,
	719.965,
	724.869,
	729.422,
	735.026,
	740.280,
	745.184,
	749.387,
	754.291,
	759.895,
	765.149,
	770.403,
	775.657,
	780.210,
	785.114,
	790.718,
	795.271,
	799.825,
	805.079,
	810.333,
	815.587,
	820.140,
	825.044,
	829.597,
	834.501,
	840.105,
	844.659,
	849.562,
	854.816,
	860.070,
	864.623,
	869.877,
	875.131,
	875.131,
	880.035,
	884.588,
	889.142,
	894.746,
	900,
],[
	0.0936174,
	0.100389,
	0.114766,
	0.135077,
	0.154423,
	0.175515,
	0.186032,
	0.193767,
	0.203001,
	0.213917,
	0.230728,
	0.245981,
	0.262242,
	0.271561,
	0.279579,
	0.287833,
	0.292902,
	0.298061,
	0.305081,
	0.310454,
	0.317765,
	0.325249,
	0.332909,
	0.340749,
	0.346750,
	0.350810,
	0.359072,
	0.361167,
	0.367528,
	0.371831,
	0.374001,
	0.378379,
	0.380588,
	0.382809,
	0.387291,
	0.387291,
	0.387291,
	0.391825,
	0.391825,
	0.394112,
	0.394112,
	0.394112,
	0.394112,
	0.394112,
	0.394112,
	0.391825,
	0.394112,
	0.394112,
	0.394112,
	0.396412,
	0.398725,
	0.396412,
	0.394112,
	0.394112,
	0.394112,
	0.394112,
	0.398725,
	0.405748,
	0.408116,
	0.408116,
	0.410498,
	0.408116,
	0.405748,
	0.401053,
	0.398725,
	0.394112,
	0.391825,
	0.387291,
	0.387291,
	0.387291,
	0.387291,
	0.387291,
	0.389551,
	0.389551,
	0.389551,
	0.387291,
	0.389551,
	0.389551,
	0.387291,
	0.391825,
	0.394112,
	0.405748,
	0.417727,
	0.430061,
	0.442758,
	0.455831,
	0.472028,
	0.491654,
	0.503233,
	0.518091,
	0.524156,
	0.536501,
	0.545950,
	0.549136,
	0.555565,
	0.562069,
	0.565349,
	0.571968,
	0.571968,
	0.575306,
	0.575306,
	0.582041,
	0.585438,
	0.585438,
	0.585438,
	0.585438,
	0.592292,
	0.588855,
	0.592292,
	0.595749,
	0.595749,
	0.595749,
	0.599226,
	0.602723,
	0.599226,
	0.599226,
	0.599226,
	0.599226,
	0.602723,
	0.606241,
	0.606241,
	0.609779,
]])

# 0.5
nd05_5nm = np.array([[
	299.688,
	305.302,
	310.565,
	314.778,
	319.691,
	324.953,
	329.163,
	334.073,
	339.336,
	345.299,
	350.212,
	354.422,
	359.685,
	364.244,
	370.557,
	374.765,
	380.376,
	385.636,
	389.143,
	394.052,
	400.364,
	404.572,
	409.131,
	414.040,
	420.001,
	424.910,
	430.169,
	435.078,
	439.637,
	444.195,
	449.454,
	455.064,
	460.323,
	465.932,
	470.841,
	475.398,
	480.307,
	484.514,
	489.773,
	495.382,
	499.589,
	505.199,
	510.457,
	515.716,
	520.274,
	525.532,
	530.440,
	535.699,
	539.555,
	544.814,
	550.073,
	555.332,
	559.889,
	565.148,
	570.406,
	575.314,
	580.222,
	585.131,
	589.689,
	595.298,
	600.206,
	605.815,
	610.723,
	615.980,
	619.836,
	625.095,
	630.353,
	635.611,
	639.818,
	645.076,
	650.686,
	655.944,
	660.151,
	664.358,
	669.616,
	675.225,
	680.484,
	685.041,
	690.651,
	695.559,
	700.468,
	704.675,
	710.637,
	715.898,
	720.106,
	724.665,
	730.276,
	734.835,
	740.096,
	745.005,
	750.966,
	755.525,
	760.433,
	765.342,
	770.601,
	775.860,
	780.418,
	785.326,
	790.936,
	795.844,
	800.752,
	805.660,
	810.919,
	815.126,
	820.385,
	825.994,
	830.902,
	835.810,
	840.018,
	844.926,
	850.886,
	855.093,
	860.001,
	865.260,
	870.167,
	874.725,
	880.335,
	885.944,
	890.852,
	895.760,
	900.668,
],[
	0.0735785,
	0.0812574,
	0.0892152,
	0.0991054,
	0.109449,
	0.117392,
	0.123724,
	0.130396,
	0.139858,
	0.150885,
	0.163737,
	0.175621,
	0.189469,
	0.197368,
	0.206795,
	0.211673,
	0.219209,
	0.225692,
	0.229672,
	0.235087,
	0.243455,
	0.249197,
	0.255074,
	0.261089,
	0.268807,
	0.273542,
	0.278358,
	0.283262,
	0.288252,
	0.291622,
	0.295028,
	0.300221,
	0.301958,
	0.305483,
	0.309052,
	0.309032,
	0.310821,
	0.312623,
	0.316274,
	0.316248,
	0.318082,
	0.319921,
	0.321771,
	0.321746,
	0.323611,
	0.323586,
	0.323563,
	0.325434,
	0.325416,
	0.325391,
	0.327273,
	0.329166,
	0.327226,
	0.329118,
	0.329093,
	0.329070,
	0.330975,
	0.334842,
	0.336782,
	0.338728,
	0.340689,
	0.340661,
	0.336678,
	0.330802,
	0.328856,
	0.328830,
	0.324984,
	0.324959,
	0.321163,
	0.321139,
	0.322995,
	0.322970,
	0.322950,
	0.324823,
	0.322905,
	0.322879,
	0.320973,
	0.322832,
	0.324698,
	0.328492,
	0.330394,
	0.336217,
	0.350226,
	0.362696,
	0.375617,
	0.386728,
	0.400497,
	0.414762,
	0.429530,
	0.442235,
	0.452655,
	0.460630,
	0.466013,
	0.474221,
	0.479759,
	0.479723,
	0.488175,
	0.488140,
	0.493838,
	0.496697,
	0.496662,
	0.499537,
	0.502426,
	0.505340,
	0.505301,
	0.508221,
	0.511163,
	0.511127,
	0.514091,
	0.514054,
	0.517022,
	0.520020,
	0.526097,
	0.522991,
	0.522953,
	0.525983,
	0.529023,
	0.528980,
	0.532042,
	0.532004,
	0.528866,
]])

# 0.6
nd06_5nm = np.array([[
	300,
	305.954,
	310.858,
	314.711,
	320.315,
	324.168,
	329.422,
	335.026,
	340.280,
	344.483,
	349.387,
	354.291,
	360.245,
	365.499,
	370.753,
	374.956,
	380.560,
	384.764,
	389.317,
	394.221,
	400.525,
	405.079,
	410.683,
	415.937,
	420.490,
	424.343,
	429.247,
	434.851,
	440.806,
	445.359,
	450.613,
	455.517,
	460.420,
	465.324,
	470.928,
	475.832,
	480.736,
	485.990,
	490.543,
	495.096,
	500,
	504.553,
	510.858,
	515.412,
	520.666,
	525.919,
	530.823,
	535.377,
	540.630,
	545.884,
	550.788,
	555.692,
	560.595,
	565.849,
	570.053,
	575.306,
	579.860,
	585.114,
	590.718,
	595.271,
	600.175,
	604.378,
	609.632,
	615.587,
	619.790,
	624.343,
	629.597,
	635.201,
	640.455,
	645.359,
	650.963,
	655.867,
	660.070,
	664.273,
	670.228,
	674.431,
	680.385,
	685.639,
	690.893,
	695.797,
	700.701,
	705.954,
	710.858,
	715.412,
	720.315,
	724.869,
	730.823,
	735.026,
	740.630,
	745.884,
	750.788,
	755.341,
	760.595,
	765.149,
	769.702,
	775.306,
	780.560,
	785.464,
	790.018,
	794.571,
	800.525,
	804.028,
	809.282,
	815.587,
	820.490,
	825.744,
	830.998,
	835.552,
	840.105,
	845.709,
	850.263,
	854.816,
	860.420,
	865.674,
	870.578,
	875.832,
	880.385,
	885.289,
	890.193,
	895.096,
	900,
],[
	0.0469609,
	0.0509594,
	0.0572692,
	0.0643604,
	0.0744715,
	0.0817612,
	0.0882065,
	0.0924229,
	0.0974077,
	0.103867,
	0.116048,
	0.128904,
	0.142350,
	0.150906,
	0.157200,
	0.160913,
	0.166648,
	0.171583,
	0.175636,
	0.180837,
	0.187282,
	0.192828,
	0.199701,
	0.205615,
	0.211704,
	0.215443,
	0.220533,
	0.224428,
	0.229729,
	0.233788,
	0.237917,
	0.240711,
	0.244963,
	0.247839,
	0.249290,
	0.252217,
	0.252217,
	0.253693,
	0.255179,
	0.258175,
	0.258175,
	0.258175,
	0.258175,
	0.256672,
	0.256672,
	0.256672,
	0.256672,
	0.258175,
	0.259686,
	0.259686,
	0.259686,
	0.258175,
	0.256672,
	0.255179,
	0.255179,
	0.256672,
	0.261206,
	0.264273,
	0.268941,
	0.270516,
	0.273692,
	0.272099,
	0.268941,
	0.267376,
	0.262735,
	0.256672,
	0.252217,
	0.250749,
	0.247839,
	0.249290,
	0.247839,
	0.249290,
	0.247839,
	0.249290,
	0.250749,
	0.249290,
	0.247839,
	0.249290,
	0.250749,
	0.250749,
	0.255179,
	0.262735,
	0.273692,
	0.288454,
	0.302242,
	0.313014,
	0.333770,
	0.351772,
	0.366442,
	0.379502,
	0.393028,
	0.407036,
	0.416650,
	0.424010,
	0.428989,
	0.431500,
	0.441693,
	0.449495,
	0.449495,
	0.446879,
	0.452126,
	0.454773,
	0.460113,
	0.462806,
	0.462806,
	0.462806,
	0.465515,
	0.468240,
	0.473738,
	0.476512,
	0.476512,
	0.476512,
	0.479301,
	0.479301,
	0.482107,
	0.484929,
	0.484929,
	0.487768,
	0.487768,
	0.487768,
	0.484929,
]])

# 0.7
nd07_5nm = np.array([[
	300,
	304.904,
	309.807,
	315.762,
	319.965,
	324.168,
	330.123,
	334.676,
	339.930,
	344.483,
	349.387,
	354.641,
	360.245,
	365.849,
	370.753,
	374.956,
	380.560,
	384.413,
	390.368,
	395.622,
	400.525,
	404.378,
	409.282,
	414.886,
	420.490,
	425.744,
	429.947,
	435.201,
	440.105,
	445.009,
	450.613,
	455.867,
	460.420,
	464.974,
	470.928,
	475.482,
	480.035,
	485.639,
	490.543,
	495.797,
	500.350,
	505.954,
	510.508,
	515.762,
	519.965,
	524.518,
	530.123,
	535.727,
	540.630,
	545.534,
	550.088,
	554.991,
	560.245,
	565.149,
	570.403,
	574.256,
	580.560,
	585.114,
	589.317,
	594.921,
	600.175,
	605.429,
	610.683,
	615.937,
	620.490,
	625.744,
	630.998,
	635.902,
	640.806,
	645.709,
	649.562,
	654.816,
	660.070,
	665.324,
	670.228,
	675.832,
	680.736,
	684.939,
	690.193,
	694.746,
	700,
	705.254,
	710.158,
	715.061,
	720.666,
	724.869,
	729.422,
	735.377,
	740.630,
	745.884,
	750.788,
	755.341,
	760.245,
	765.849,
	770.753,
	775.657,
	780.560,
	785.114,
	790.368,
	795.622,
	800.175,
	804.729,
	809.632,
	814.536,
	820.140,
	825.394,
	830.998,
	835.552,
	840.455,
	844.659,
	849.212,
	854.466,
	859.720,
	865.324,
	870.578,
	875.832,
	880.736,
	885.289,
	890.543,
	895.797,
	900,
],[
	0.0285559,
	0.0308040,
	0.0348155,
	0.0414688,
	0.0476959,
	0.0520542,
	0.0568107,
	0.0602208,
	0.0638357,
	0.0692638,
	0.0778286,
	0.0889953,
	0.0988409,
	0.106003,
	0.111063,
	0.115016,
	0.119110,
	0.122633,
	0.128487,
	0.132288,
	0.136996,
	0.141048,
	0.146923,
	0.153042,
	0.158489,
	0.162228,
	0.167026,
	0.171966,
	0.176023,
	0.179128,
	0.182289,
	0.185504,
	0.187680,
	0.189881,
	0.193231,
	0.195497,
	0.196640,
	0.198946,
	0.200109,
	0.201279,
	0.203639,
	0.203639,
	0.203639,
	0.204829,
	0.204829,
	0.204829,
	0.204829,
	0.206027,
	0.208443,
	0.209662,
	0.209662,
	0.208443,
	0.207231,
	0.206027,
	0.206027,
	0.208443,
	0.210887,
	0.214608,
	0.218394,
	0.220955,
	0.223546,
	0.220955,
	0.218394,
	0.215863,
	0.212120,
	0.207231,
	0.206027,
	0.204829,
	0.202455,
	0.202455,
	0.202455,
	0.202455,
	0.202455,
	0.203639,
	0.204829,
	0.202455,
	0.201279,
	0.202455,
	0.203639,
	0.203639,
	0.208443,
	0.215863,
	0.224853,
	0.235587,
	0.248277,
	0.261650,
	0.278978,
	0.295724,
	0.309841,
	0.326531,
	0.340131,
	0.350190,
	0.358451,
	0.362655,
	0.371211,
	0.379968,
	0.382189,
	0.386671,
	0.388932,
	0.391206,
	0.395793,
	0.400435,
	0.402776,
	0.402776,
	0.405131,
	0.409881,
	0.412278,
	0.414688,
	0.414688,
	0.409881,
	0.414688,
	0.422004,
	0.424471,
	0.424471,
	0.426953,
	0.426953,
	0.429449,
	0.431960,
	0.434485,
	0.437025,
	0.437025,
]])

# 0.8
nd08_5nm = np.array([[
	300,
	305.954,
	310.508,
	315.412,
	320.315,
	324.869,
	330.123,
	334.326,
	340.280,
	344.834,
	350.788,
	354.641,
	359.895,
	365.499,
	370.753,
	374.606,
	380.210,
	384.063,
	389.667,
	395.622,
	400.525,
	405.429,
	410.333,
	415.236,
	420.841,
	425.394,
	430.298,
	435.902,
	440.806,
	445.709,
	450.613,
	455.166,
	460.420,
	464.974,
	470.578,
	475.482,
	480.736,
	485.639,
	489.842,
	495.447,
	500.350,
	504.904,
	510.158,
	515.061,
	520.315,
	525.219,
	530.473,
	535.377,
	540.280,
	545.184,
	549.387,
	554.641,
	559.895,
	565.149,
	570.053,
	575.657,
	580.560,
	585.464,
	590.718,
	595.622,
	600.525,
	605.779,
	610.333,
	615.937,
	620.490,
	625.044,
	629.247,
	635.201,
	640.105,
	645.359,
	649.912,
	655.166,
	659.720,
	665.324,
	669.877,
	675.482,
	680.385,
	684.939,
	689.842,
	695.447,
	700.350,
	705.254,
	710.508,
	714.711,
	720.315,
	724.518,
	729.422,
	735.026,
	740.280,
	745.184,
	749.737,
	754.641,
	760.245,
	765.149,
	770.053,
	774.956,
	780.210,
	785.464,
	790.368,
	794.921,
	800.175,
	805.079,
	810.333,
	815.236,
	820.490,
	825.044,
	830.298,
	835.201,
	840.455,
	845.359,
	850.263,
	854.816,
	859.720,
	864.623,
	870.228,
	875.131,
	880.035,
	884.939,
	890.193,
	894.746,
	900,
],[
	0.0182567,
	0.0201636,
	0.0227965,
	0.0270067,
	0.0314385,
	0.0353365,
	0.0388001,
	0.0408955,
	0.0436108,
	0.0476064,
	0.0570618,
	0.0633914,
	0.0712511,
	0.0773259,
	0.0819796,
	0.0844104,
	0.0884505,
	0.0916071,
	0.0954323,
	0.100000,
	0.103569,
	0.107893,
	0.113058,
	0.117778,
	0.122697,
	0.126335,
	0.130081,
	0.133938,
	0.137909,
	0.141171,
	0.143668,
	0.146209,
	0.148795,
	0.151427,
	0.153207,
	0.155917,
	0.157750,
	0.159605,
	0.160540,
	0.161481,
	0.162428,
	0.163380,
	0.164337,
	0.164337,
	0.164337,
	0.164337,
	0.164337,
	0.166269,
	0.168224,
	0.169210,
	0.170202,
	0.169210,
	0.167244,
	0.166269,
	0.167244,
	0.168224,
	0.170202,
	0.175249,
	0.178348,
	0.179394,
	0.181503,
	0.179394,
	0.177309,
	0.175249,
	0.173212,
	0.169210,
	0.167244,
	0.165300,
	0.163380,
	0.162428,
	0.163380,
	0.163380,
	0.164337,
	0.164337,
	0.165300,
	0.164337,
	0.164337,
	0.163380,
	0.165300,
	0.167244,
	0.169210,
	0.174228,
	0.183637,
	0.192426,
	0.206405,
	0.216284,
	0.231997,
	0.251777,
	0.266929,
	0.278075,
	0.289687,
	0.300025,
	0.310732,
	0.318081,
	0.323707,
	0.329433,
	0.337224,
	0.341189,
	0.345200,
	0.345200,
	0.347224,
	0.351306,
	0.355436,
	0.355436,
	0.357519,
	0.359615,
	0.361723,
	0.365975,
	0.370278,
	0.368120,
	0.370278,
	0.372448,
	0.374631,
	0.376827,
	0.379036,
	0.381257,
	0.383492,
	0.381257,
	0.383492,
	0.385740,
	0.388001,
]])

# 0.9
nd09_5nm = np.array([[
	300,
	304.904,
	310.158,
	315.412,
	319.965,
	325.219,
	330.473,
	335.026,
	339.930,
	344.483,
	349.387,
	354.641,
	360.245,
	365.499,
	370.753,
	375.306,
	380.560,
	385.814,
	390.018,
	394.921,
	399.825,
	405.429,
	409.982,
	415.236,
	420.490,
	425.394,
	429.597,
	434.501,
	440.105,
	444.659,
	450.613,
	455.867,
	460.420,
	464.623,
	469.877,
	475.832,
	480.035,
	485.289,
	489.842,
	495.096,
	500.350,
	505.254,
	510.158,
	515.412,
	519.965,
	524.869,
	530.123,
	535.026,
	540.280,
	545.534,
	550.438,
	554.991,
	559.895,
	565.149,
	570.403,
	574.606,
	579.860,
	585.114,
	589.667,
	594.571,
	600.175,
	605.079,
	610.333,
	615.236,
	620.140,
	625.044,
	630.298,
	635.201,
	639.755,
	645.009,
	650.263,
	655.166,
	660.070,
	664.974,
	669.877,
	675.482,
	680.035,
	684.939,
	690.193,
	695.096,
	700,
	705.254,
	710.158,
	715.061,
	719.965,
	724.869,
	730.123,
	735.026,
	740.280,
	745.184,
	749.387,
	754.991,
	760.245,
	765.149,
	770.403,
	774.956,
	780.210,
	785.114,
	790.018,
	795.271,
	799.825,
	804.729,
	810.333,
	815.236,
	819.790,
	825.044,
	829.947,
	835.201,
	840.455,
	844.659,
	849.912,
	855.166,
	860.420,
	865.324,
	870.228,
	874.781,
	879.685,
	884.939,
	890.193,
	895.096,
	900,
],[
	0.0128904,
	0.0145713,
	0.0171583,
	0.0205615,
	0.0239310,
	0.0265820,
	0.0295267,
	0.0313014,
	0.0349725,
	0.0386206,
	0.0436567,
	0.0508109,
	0.0584512,
	0.0634281,
	0.0680299,
	0.0708669,
	0.0742544,
	0.0778039,
	0.0810485,
	0.0839372,
	0.0874376,
	0.0921535,
	0.0959966,
	0.100000,
	0.104170,
	0.108515,
	0.110431,
	0.114367,
	0.117069,
	0.119834,
	0.123383,
	0.125562,
	0.127780,
	0.130037,
	0.131564,
	0.133109,
	0.135460,
	0.137051,
	0.138660,
	0.138660,
	0.140288,
	0.141110,
	0.141936,
	0.141936,
	0.141936,
	0.142766,
	0.142766,
	0.143602,
	0.145288,
	0.146139,
	0.146994,
	0.146139,
	0.146139,
	0.146139,
	0.146139,
	0.146994,
	0.148720,
	0.151347,
	0.153124,
	0.154922,
	0.154922,
	0.154922,
	0.153124,
	0.152233,
	0.149591,
	0.146994,
	0.145288,
	0.142766,
	0.141936,
	0.141110,
	0.141110,
	0.141936,
	0.141936,
	0.142766,
	0.142766,
	0.142766,
	0.141936,
	0.141936,
	0.141936,
	0.143602,
	0.147855,
	0.152233,
	0.158582,
	0.168113,
	0.180310,
	0.192267,
	0.206216,
	0.219890,
	0.231750,
	0.245679,
	0.254435,
	0.266597,
	0.276099,
	0.282621,
	0.287613,
	0.292694,
	0.297864,
	0.301361,
	0.304900,
	0.306685,
	0.310286,
	0.313929,
	0.315767,
	0.319474,
	0.319474,
	0.323226,
	0.327021,
	0.327021,
	0.328935,
	0.330861,
	0.330861,
	0.334746,
	0.336705,
	0.340659,
	0.340659,
	0.340659,
	0.342653,
	0.346676,
	0.348706,
	0.348706,
	0.350747,
]])

# 1.0
nd10_5nm = np.array([[
	299.914,
	304.829,
	309.395,
	314.664,
	319.580,
	324.846,
	330.111,
	335.023,
	339.585,
	344.501,
	350.119,
	355.387,
	359.950,
	364.862,
	370.827,
	375.036,
	379.947,
	385.208,
	390.820,
	395.380,
	399.939,
	405.551,
	410.462,
	414.671,
	419.932,
	425.192,
	430.453,
	435.362,
	439.921,
	445.181,
	450.091,
	455.701,
	459.909,
	465.168,
	470.778,
	475.336,
	479.894,
	484.802,
	490.412,
	495.321,
	499.878,
	504.787,
	510.396,
	515.304,
	520.212,
	525.120,
	530.028,
	535.637,
	540.195,
	545.104,
	550.012,
	555.972,
	559.828,
	565.436,
	569.994,
	575.253,
	579.811,
	584.720,
	590.331,
	595.239,
	599.796,
	605.405,
	610.313,
	615.220,
	620.127,
	625.034,
	629.941,
	634.848,
	640.106,
	644.663,
	649.571,
	654.829,
	659.738,
	664.295,
	670.255,
	675.514,
	680.421,
	685.329,
	690.237,
	695.497,
	700.056,
	704.966,
	710.228,
	715.139,
	720.401,
	724.963,
	730.225,
	735.488,
	740.049,
	745.310,
	750.221,
	754.780,
	760.041,
	765.301,
	769.859,
	774.768,
	779.677,
	785.287,
	790.196,
	795.455,
	800.363,
	805.272,
	810.180,
	814.738,
	820.348,
	825.256,
	830.164,
	835.073,
	839.981,
	844.539,
	850.148,
	855.407,
	859.965,
	865.224,
	870.132,
	874.690,
	880.300,
	884.857,
	890.116,
	895.024,
	900.283,
],[
	0.00685623,
	0.00788942,
	0.00923913,
	0.0113382,
	0.0132778,
	0.0152786,
	0.0171742,
	0.0186391,
	0.0204671,
	0.0236896,
	0.0282332,
	0.0336485,
	0.0369486,
	0.0405721,
	0.0445502,
	0.0466825,
	0.0489162,
	0.0518600,
	0.0546599,
	0.0569415,
	0.0589722,
	0.0625208,
	0.0658969,
	0.0686480,
	0.0715127,
	0.0744970,
	0.0771531,
	0.0794381,
	0.0817911,
	0.0842130,
	0.0867071,
	0.0882354,
	0.0903194,
	0.0919119,
	0.0935319,
	0.0946267,
	0.0957342,
	0.0957274,
	0.0974147,
	0.0979793,
	0.0985478,
	0.0991190,
	0.0996925,
	0.100270,
	0.0996782,
	0.0996710,
	0.0996639,
	0.100241,
	0.101414,
	0.102002,
	0.102593,
	0.102584,
	0.101980,
	0.100785,
	0.101370,
	0.101957,
	0.103755,
	0.106205,
	0.109349,
	0.110628,
	0.109976,
	0.109325,
	0.108045,
	0.106781,
	0.104915,
	0.102480,
	0.100690,
	0.0989306,
	0.0977720,
	0.0971951,
	0.0977585,
	0.0971806,
	0.0977439,
	0.0983110,
	0.0983024,
	0.0982949,
	0.0977144,
	0.0977074,
	0.0982737,
	0.100006,
	0.102969,
	0.106640,
	0.113722,
	0.121273,
	0.130085,
	0.141181,
	0.152328,
	0.165320,
	0.175271,
	0.184734,
	0.193574,
	0.201654,
	0.211302,
	0.216290,
	0.220105,
	0.225302,
	0.229275,
	0.233317,
	0.236046,
	0.237413,
	0.240191,
	0.244427,
	0.245844,
	0.247270,
	0.248701,
	0.251610,
	0.251592,
	0.254536,
	0.257514,
	0.257496,
	0.260506,
	0.262015,
	0.263535,
	0.266617,
	0.266598,
	0.269718,
	0.271279,
	0.271261,
	0.272832,
	0.276024,
	0.276002,
]])

# 2.0
nd20_5nm = np.array([[
	300,
	312,
	312.248,
	315.748,
	320.297,
	325.197,
	330.446,
	334.996,
	340.245,
	345.144,
	349.694,
	355.293,
	360.192,
	364.742,
	370.341,
	375.241,
	380.140,
	385.039,
	389.939,
	395.188,
	400.087,
	404.987,
	409.536,
	415.136,
	420.385,
	425.284,
	430.184,
	435.433,
	440.332,
	444.182,
	449.781,
	455.031,
	460.630,
	465.179,
	470.429,
	475.328,
	480.227,
	485.477,
	490.376,
	495.276,
	500.175,
	505.074,
	510.324,
	515.223,
	520.122,
	525.372,
	530.271,
	534.821,
	540.070,
	545.319,
	550.219,
	555.468,
	560.018,
	564.917,
	569.816,
	575.066,
	579.615,
	584.864,
	590.114,
	595.363,
	599.913,
	604.812,
	610.061,
	614.961,
	620.210,
	625.459,
	629.659,
	635.258,
	640.157,
	645.407,
	649.956,
	654.856,
	659.755,
	665.354,
	670.254,
	675.153,
	680.402,
	684.602,
	689.851,
	694.751,
	700.350,
	705.249,
	710.149,
	715.048,
	719.948,
	724.847,
	730.096,
	734.996,
	739.895,
	745.144,
	750.394,
	754.943,
	759.843,
	765.092,
	769.991,
	774.891,
	779.790,
	785.389,
	790.289,
	794.838,
	799.738,
	804.987,
	809.886,
	815.136,
	820.035,
	824.584,
	830.184,
	835.433,
	839.633,
	845.232,
	850.131,
	855.381,
	859.580,
	865.179,
	870.079,
	874.628,
	879.528,
	885.477,
	889.676,
	895.276,
	900.175,
],[
	0,
	0,
	9.92267e-05,
	0.000101565,
	0.000152072,
	0.000154451,
	0.000184642,
	0.000222454,
	0.000278616,
	0.000413943,
	0.000498713,
	0.000788415,
	0.00101861,
	0.00121772,
	0.00134701,
	0.00158551,
	0.00178130,
	0.00197044,
	0.00214609,
	0.00237396,
	0.00262603,
	0.00295032,
	0.00331465,
	0.00369518,
	0.00415150,
	0.00448660,
	0.00488654,
	0.00532213,
	0.00570724,
	0.00597935,
	0.00641203,
	0.00677009,
	0.00709287,
	0.00737357,
	0.00772513,
	0.00803086,
	0.00822007,
	0.00841375,
	0.00854539,
	0.00881489,
	0.00902258,
	0.00916375,
	0.00930713,
	0.00923516,
	0.00930713,
	0.00930713,
	0.00937966,
	0.00960066,
	0.00998061,
	0.0101368,
	0.0102158,
	0.0101368,
	0.0100584,
	0.00990344,
	0.00998061,
	0.0101368,
	0.0106201,
	0.0113886,
	0.0120245,
	0.0124037,
	0.0126960,
	0.0126960,
	0.0125004,
	0.0121182,
	0.0116569,
	0.0112131,
	0.0107862,
	0.0104565,
	0.0101368,
	0.0100584,
	0.0100584,
	0.0102158,
	0.0103756,
	0.0105379,
	0.0106201,
	0.0106201,
	0.0107028,
	0.0106201,
	0.0105379,
	0.0107862,
	0.0112131,
	0.0119315,
	0.0131984,
	0.0150603,
	0.0173188,
	0.0203851,
	0.0241814,
	0.0282427,
	0.0327311,
	0.0376394,
	0.0419605,
	0.0460570,
	0.0501625,
	0.0537924,
	0.0567962,
	0.0590440,
	0.0613807,
	0.0638099,
	0.0658223,
	0.0673732,
	0.0689606,
	0.0694980,
	0.0711354,
	0.0728115,
	0.0739507,
	0.0751078,
	0.0762829,
	0.0774765,
	0.0780802,
	0.0799199,
	0.0818029,
	0.0824404,
	0.0837303,
	0.0837303,
	0.0857030,
	0.0870440,
	0.0877223,
	0.0890948,
	0.0897891,
	0.0919047,
	0.0919047,
]])

# 3.0
nd30_5nm = np.array([[
	300,
	385,
	385.727,
	389.931,
	394.835,
	400.088,
	404.988,
	410.239,
	414.791,
	420.394,
	425.298,
	430.201,
	435.455,
	440.359,
	445.263,
	450.869,
	455.774,
	460.329,
	464.884,
	470.140,
	475.046,
	480.303,
	484.858,
	490.115,
	495.022,
	499.929,
	505.186,
	509.743,
	515.351,
	520.259,
	524.816,
	530.424,
	534.981,
	540.238,
	545.145,
	550.052,
	555.310,
	560.218,
	565.126,
	569.684,
	574.942,
	579.848,
	584.754,
	590.010,
	595.267,
	600.174,
	605.082,
	610.341,
	614.549,
	620.160,
	625.070,
	629.981,
	635.241,
	640.151,
	645.411,
	649.617,
	655.226,
	659.782,
	665.391,
	670.298,
	675.206,
	680.465,
	685.373,
	690.282,
	695.189,
	700.445,
	704.649,
	709.902,
	715.152,
	720.401,
	725.299,
	730.197,
	735.446,
	740.345,
	745.245,
	750.147,
	755.048,
	760.302,
	764.856,
	769.761,
	775.016,
	780.272,
	785.178,
	790.084,
	795.341,
	799.897,
	805.154,
	810.412,
	815.318,
	820.225,
	825.133,
	830.390,
	834.596,
	840.204,
	845.111,
	849.667,
	854.924,
	860.180,
	865.437,
	870.344,
	875.251,
	879.808,
	885.417,
	890.325,
	895.232,
	900.140,
],[
	0,
	0,
	0.000101343,
	0.000108705,
	0.000119357,
	0.000137319,
	0.000166827,
	0.000204261,
	0.000233174,
	0.000270361,
	0.000303866,
	0.000341524,
	0.000386849,
	0.000428076,
	0.000466380,
	0.000508119,
	0.000549294,
	0.000589198,
	0.000627103,
	0.000667454,
	0.000694004,
	0.000727252,
	0.000762083,
	0.000798593,
	0.000823922,
	0.000843465,
	0.000863478,
	0.000877100,
	0.000897916,
	0.000912089,
	0.000919296,
	0.000933818,
	0.000963430,
	0.000986288,
	0.00100968,
	0.00104171,
	0.00104995,
	0.00105005,
	0.00104201,
	0.00105025,
	0.00106683,
	0.00110927,
	0.00117148,
	0.00123720,
	0.00130660,
	0.00132722,
	0.00133772,
	0.00131719,
	0.00127692,
	0.00119998,
	0.00113647,
	0.00107632,
	0.00101936,
	0.000980547,
	0.000950591,
	0.000958096,
	0.000973231,
	0.000996308,
	0.00101205,
	0.00102802,
	0.00103615,
	0.00102025,
	0.00100459,
	0.000996900,
	0.00102055,
	0.00107779,
	0.00115608,
	0.00133006,
	0.00164123,
	0.00210556,
	0.00270123,
	0.00349251,
	0.00448059,
	0.00557200,
	0.00682225,
	0.00803424,
	0.00938820,
	0.0106342,
	0.0116762,
	0.0126223,
	0.0137519,
	0.0145233,
	0.0153379,
	0.0160726,
	0.0167120,
	0.0173766,
	0.0179279,
	0.0184967,
	0.0190833,
	0.0195360,
	0.0199993,
	0.0206338,
	0.0209592,
	0.0214566,
	0.0219656,
	0.0226621,
	0.0241203,
	0.0254733,
	0.0262815,
	0.0269048,
	0.0275430,
	0.0279775,
	0.0279806,
	0.0282019,
	0.0286470,
	0.0288736,
]]) 

# 4.0 
nd40_5nm = np.array([[
	300,
	406,
	406.398,
	409.904,
	414.461,
	420.070,
	424.628,
	430.587,
	434.794,
	440.053,
	444.961,
	449.518,
	455.478,
	459.684,
	464.943,
	469.851,
	474.759,
	480.368,
	484.926,
	490.184,
	495.092,
	499.649,
	504.908,
	510.167,
	515.425,
	519.982,
	525.241,
	530.149,
	535.057,
	539.614,
	544.873,
	550.482,
	555.390,
	560.298,
	564.505,
	569.763,
	575.022,
	579.930,
	584.838,
	590.096,
	595.355,
	599.912,
	604.820,
	609.728,
	614.987,
	620.245,
	624.803,
	629.711,
	634.969,
	639.877,
	645.486,
	650.394,
	655.302,
	659.860,
	665.118,
	669.676,
	675.285,
	680.193,
	685.101,
	689.658,
	695.267,
	699.825,
	704.733,
	709.991,
	714.899,
	720.158,
	724.715,
	729.974,
	734.882,
	740.491,
	745.048,
	750.307,
	754.514,
	760.473,
	765.381,
	770.289,
	774.847,
	780.105,
	785.364,
	790.272,
	795.180,
	800.438,
	804.996,
	809.553,
	815.162,
	820.421,
	824.628,
	829.886,
	835.145,
	840.053,
	845.311,
	849.167,
	855.127,
	860.386,
	865.294,
	870.202,
	874.759,
	880.368,
	884.925,
	889.833,
	895.443,
	900,
],[
	0,
	0,
	1.02059e-05,
	1.22959e-05,
	1.44729e-05,
	1.71349e-05,
	1.98195e-05,
	2.42990e-05,
	2.76194e-05,
	3.19466e-05,
	3.63120e-05,
	3.89396e-05,
	4.12740e-05,
	4.47790e-05,
	5.14941e-05,
	5.65213e-05,
	5.95619e-05,
	6.24016e-05,
	6.57585e-05,
	6.88937e-05,
	7.13428e-05,
	7.34501e-05,
	7.51807e-05,
	7.69520e-05,
	7.83079e-05,
	7.92250e-05,
	7.92250e-05,
	8.06209e-05,
	8.34869e-05,
	8.69596e-05,
	9.00509e-05,
	9.27108e-05,
	9.32521e-05,
	9.27108e-05,
	9.11056e-05,
	9.00509e-05,
	9.48951e-05,
	0.000101762,
	0.000112349,
	0.000124037,
	0.000133013,
	0.000136147,
	0.000136942,
	0.000135357,
	0.000127701,
	0.000117705,
	0.000109126,
	0.000101171,
	9.37966e-05,
	9.05767e-05,
	8.69596e-05,
	8.69596e-05,
	8.90084e-05,
	9.11056e-05,
	9.43443e-05,
	9.60066e-05,
	9.54492e-05,
	9.43443e-05,
	9.27108e-05,
	9.21726e-05,
	9.60066e-05,
	0.000103555,
	0.000117022,
	0.000143471,
	0.000186442,
	0.000255317,
	0.000347606,
	0.000498713,
	0.000682947,
	0.000940700,
	0.00116682,
	0.00151630,
	0.00179518,
	0.00215026,
	0.00242990,
	0.00269836,
	0.00299648,
	0.00326993,
	0.00348619,
	0.00365241,
	0.00384889,
	0.00405594,
	0.00420012,
	0.00440037,
	0.00461017,
	0.00477406,
	0.00488654,
	0.00500167,
	0.00520972,
	0.00542642,
	0.00568514,
	0.00581908,
	0.00592161,
	0.00609652,
	0.00635011,
	0.00673079,
	0.00713428,
	0.00730237,
	0.00725998,
	0.00730237,
	0.00806209,
	0.00830022,
]])

# 1-nm interpolation from 300 to 700 nm
x_1nm = np.empty(401)
for i in range(401): x_1nm[i] = i + 300
yellow15_1nm = np.interp(x_1nm, yellow15_5nm[0], yellow15_5nm[1])
red25_1nm = np.interp(x_1nm, red25_5nm[0], red25_5nm[1])
blue47_1nm = np.interp(x_1nm, blue47_5nm[0], blue47_5nm[1])
green58_1nm = np.interp(x_1nm, green58_5nm[0], green58_5nm[1])
nd01_1nm = np.interp(x_1nm, nd01_5nm[0], nd01_5nm[1])
nd02_1nm = np.interp(x_1nm, nd02_5nm[0], nd02_5nm[1])
nd03_1nm = np.interp(x_1nm, nd03_5nm[0], nd03_5nm[1])
nd04_1nm = np.interp(x_1nm, nd04_5nm[0], nd04_5nm[1])
nd05_1nm = np.interp(x_1nm, nd05_5nm[0], nd05_5nm[1])
nd06_1nm = np.interp(x_1nm, nd06_5nm[0], nd06_5nm[1])
nd07_1nm = np.interp(x_1nm, nd07_5nm[0], nd07_5nm[1])
nd08_1nm = np.interp(x_1nm, nd08_5nm[0], nd08_5nm[1])
nd09_1nm = np.interp(x_1nm, nd09_5nm[0], nd09_5nm[1])
nd10_1nm = np.interp(x_1nm, nd10_5nm[0], nd10_5nm[1])
nd20_1nm = np.interp(x_1nm, nd20_5nm[0], nd20_5nm[1])
nd30_1nm = np.interp(x_1nm, nd30_5nm[0], nd30_5nm[1])
nd40_1nm = np.interp(x_1nm, nd40_5nm[0], nd40_5nm[1])

# Alternate data for 400-700 nm from Kodak Photographic Filters Handbook (1990). There are
# also graphs going down to 300 that we could scan, but there's no information for neutral
# density filters other than 1.0. The biggest difference below 400 nm evident in these
# and earlier published graphs is that the Wratten 2 glass version of red 25 transmits a
# small amount of light in this range whereas all other sources show no transmission at all.
# The other four are consistently shown as either transmitting UV (yellow, blue, ND 1.0) or
# not (green).
if (args.kv2):
	# yellow 15 -- page 104
	yellow15_10nm = np.array([
		0, # 400
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0, # 500
		1.0,
		19.4,
		56.2,
		77.6,
		85.6,
		88.2,
		89.3,
		89.8,
		90.1,
		90.4, # 600
		90.5,
		90.6,
		90.7,
		90.9,
		91.0,
		91.1,
		91.1,
		91.1,
		91.1,
		91.1 # 700
	])

	# red 25 -- page 111
	red25_10nm = np.array([
		0, # 400
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0, # 500
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		12.6,
		50.0, # 600
		75.0,
		82.6,
		85.5,
		86.7,
		87.6,
		88.2,
		88.5,
		89.0,
		89.3,
		89.5 # 700
	])

	# blue 47 -- page 121
	blue47_10nm = np.array([
		9.7, # 400
		21.8,
		37.8,
		47.8,
		50.3,
		48.2,
		42.8,
		35.7,
		27.1,
		18.2,
		10.2, # 500
		4.3,
		1.2,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0, # 600
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0 # 700
	])

	# green 58 -- page 124
	green58_10nm = np.array([
		0, # 400
		0,
		0,
		0,
		0,
		0,
		0,
		0.23,
		1.38,
		4.90,
		17.7, # 500
		38.8,
		52.2,
		53.6,
		47.6,
		38.4,
		27.8,
		17.4,
		9.0,
		3.50,
		1.50, # 600
		0.41,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0.53 # 700
	])

	# neutral density 96 -- page 132
	nd10_10nm = np.array([
		4.28, # 400
		4.91,
		5.50,
		6.17,
		6.92,
		7.50,
		7.81,
		8.15,
		8.47,
		8.60,
		8.73, # 500
		8.85,
		8.90,
		9.01,
		9.07,
		9.20,
		9.30,
		9.20,
		9.19,
		9.54,
		9.64, # 600
		9.73,
		9.56,
		9.27,
		9.10,
		9.07,
		9.00,
		9.13,
		9.08,
		9.21,
		9.52 # 700
	])
	
	# interpolate
	x_400700_10nm = np.empty(31)
	for i in range(31): x_400700_10nm[i] = i*10 + 400
	x_400700_1nm = np.empty(301)
	for i in range(301): x_400700_1nm[i] = i + 400
	yellow15_400700 = np.interp(x_400700_1nm, x_400700_10nm, yellow15_10nm)
	red25_400700 = np.interp(x_400700_1nm, x_400700_10nm, red25_10nm)
	blue47_400700 = np.interp(x_400700_1nm, x_400700_10nm, blue47_10nm)
	green58_400700 = np.interp(x_400700_1nm, x_400700_10nm, green58_10nm)
	nd10_400700 = np.interp(x_400700_1nm, x_400700_10nm, nd10_10nm)
	
	# replace
	for i in range(100, 401):
		yellow15_1nm[i] = yellow15_400700[i-100]/100
		red25_1nm[i] = red25_400700[i-100]/100
		blue47_1nm[i] = blue47_400700[i-100]/100
		green58_1nm[i] = green58_400700[i-100]/100
		nd10_1nm[i] = nd10_400700[i-100]/100

# Find brightness matches for human vision. Since we begin with colors that are much
# brighter than blue 47, we choose the "brightest" ND filter or pair of filters that produces
# an achromatic contrast less than 1, but in practice there seems to be only one match for each
# color. This depends on the specified vision system and illuminant like everything else, so
# we have to set --filter human --white i (L, M and S don't matter). The current matches are:
# * red: 1.0 + 0.5
# * yellow: 2.0
# * green: 1.0 + 0.3
# * gray: 2.0 + 0.1
# When using the 1990 data, the matches are one or two steps brighter:
# * red: 1.0 + 0.3
# * yellow: 1.0 + 0.9
# * green: 1.0 + 0.2
# * gray: 2.0
if (args.kcheck):	
	# "dummy" array for gray as we assume no color filter
	placeholder = np.empty(401)
	for i in range(401):
		placeholder[i] = 1
	
	# brightness tests -- wrap this up in a function so we don't have to copy it a million
	# times
	def brightness_tests(f):
		f_01 = np.empty(401)
		for i in range(401):
			f_01[i] = f[i] * nd01_1nm[i]
		f_02 = np.empty(401)
		for i in range(401):
			f_02[i] = f[i] * nd02_1nm[i]
		f_03 = np.empty(401)
		for i in range(401):
			f_03[i] = f[i] * nd03_1nm[i]
		f_04 = np.empty(401)
		for i in range(401):
			f_04[i] = f[i] * nd04_1nm[i]
		f_05 = np.empty(401)
		for i in range(401):
			f_05[i] = f[i] * nd05_1nm[i]
		f_06 = np.empty(401)
		for i in range(401):
			f_06[i] = f[i] * nd06_1nm[i]
		f_07 = np.empty(401)
		for i in range(401):
			f_07[i] = f[i] * nd07_1nm[i]
		f_08 = np.empty(401)
		for i in range(401):
			f_08[i] = f[i] * nd08_1nm[i]
		f_09 = np.empty(401)
		for i in range(401):
			f_09[i] = f[i] * nd09_1nm[i]
		# 1.0 + others
		f_10 = np.empty(401)
		for i in range(401):
			f_10[i] = f[i] * nd10_1nm[i]
		f_1001 = np.empty(401)
		for i in range(401):
			f_1001[i] = f[i] * nd10_1nm[i] * nd01_1nm[i]
		f_1002 = np.empty(401)
		for i in range(401):
			f_1002[i] = f[i] * nd10_1nm[i] * nd02_1nm[i]
		f_1003 = np.empty(401)
		for i in range(401):
			f_1003[i] = f[i] * nd10_1nm[i] * nd03_1nm[i]
		f_1004 = np.empty(401)
		for i in range(401):
			f_1004[i] = f[i] * nd10_1nm[i] * nd04_1nm[i]
		f_1005 = np.empty(401)
		for i in range(401):
			f_1005[i] = f[i] * nd10_1nm[i] * nd05_1nm[i]
		f_1006 = np.empty(401)
		for i in range(401):
			f_1006[i] = f[i] * nd10_1nm[i] * nd06_1nm[i]
		f_1007 = np.empty(401)
		for i in range(401):
			f_1007[i] = f[i] * nd10_1nm[i] * nd07_1nm[i]
		f_1008 = np.empty(401)
		for i in range(401):
			f_1008[i] = f[i] * nd10_1nm[i] * nd08_1nm[i]
		f_1009 = np.empty(401)
		for i in range(401):
			f_1009[i] = f[i] * nd10_1nm[i] * nd09_1nm[i]
		# 2.0
		f_20 = np.empty(401)
		for i in range(401):
			f_20[i] = f[i] * nd20_1nm[i]
		f_2001 = np.empty(401)
		for i in range(401):
			f_2001[i] = f[i] * nd20_1nm[i] * nd01_1nm[i]
		f_2002 = np.empty(401)
		for i in range(401):
			f_2002[i] = f[i] * nd20_1nm[i] * nd02_1nm[i]
		f_2003 = np.empty(401)
		for i in range(401):
			f_2003[i] = f[i] * nd20_1nm[i] * nd03_1nm[i]
		f_2004 = np.empty(401)
		for i in range(401):
			f_2004[i] = f[i] * nd20_1nm[i] * nd04_1nm[i]
		f_2005 = np.empty(401)
		for i in range(401):
			f_2005[i] = f[i] * nd20_1nm[i] * nd05_1nm[i]
		f_2006 = np.empty(401)
		for i in range(401):
			f_2006[i] = f[i] * nd20_1nm[i] * nd06_1nm[i]
		f_2007 = np.empty(401)
		for i in range(401):
			f_2007[i] = f[i] * nd20_1nm[i] * nd07_1nm[i]
		f_2008 = np.empty(401)
		for i in range(401):
			f_2008[i] = f[i] * nd20_1nm[i] * nd08_1nm[i]
		f_2009 = np.empty(401)
		for i in range(401):
			f_2009[i] = f[i] * nd20_1nm[i] * nd09_1nm[i]
		# 3.0
		f_30 = np.empty(401)
		for i in range(401):
			f_30[i] = f[i] * nd30_1nm[i]
		f_3001 = np.empty(401)
		for i in range(401):
			f_3001[i] = f[i] * nd30_1nm[i] * nd01_1nm[i]
		f_3002 = np.empty(401)
		for i in range(401):
			f_3002[i] = f[i] * nd30_1nm[i] * nd02_1nm[i]
		f_3003 = np.empty(401)
		for i in range(401):
			f_3003[i] = f[i] * nd30_1nm[i] * nd03_1nm[i]
		f_3004 = np.empty(401)
		for i in range(401):
			f_3004[i] = f[i] * nd30_1nm[i] * nd04_1nm[i]
		f_3005 = np.empty(401)
		for i in range(401):
			f_3005[i] = f[i] * nd30_1nm[i] * nd05_1nm[i]
		f_3006 = np.empty(401)
		for i in range(401):
			f_3006[i] = f[i] * nd30_1nm[i] * nd06_1nm[i]
		f_3007 = np.empty(401)
		for i in range(401):
			f_3007[i] = f[i] * nd30_1nm[i] * nd07_1nm[i]
		f_3008 = np.empty(401)
		for i in range(401):
			f_3008[i] = f[i] * nd30_1nm[i] * nd08_1nm[i]
		f_3009 = np.empty(401)
		for i in range(401):
			f_3009[i] = f[i] * nd30_1nm[i] * nd09_1nm[i]
		# 4.0
		f_40 = np.empty(401)
		for i in range(401):
			f_40[i] = f[i] * nd40_1nm[i]
		f_4001 = np.empty(401)
		for i in range(401):
			f_4001[i] = f[i] * nd40_1nm[i] * nd01_1nm[i]
		f_4002 = np.empty(401)
		for i in range(401):
			f_4002[i] = f[i] * nd40_1nm[i] * nd02_1nm[i]
		f_4003 = np.empty(401)
		for i in range(401):
			f_4003[i] = f[i] * nd40_1nm[i] * nd03_1nm[i]
		f_4004 = np.empty(401)
		for i in range(401):
			f_4004[i] = f[i] * nd40_1nm[i] * nd04_1nm[i]
		f_4005 = np.empty(401)
		for i in range(401):
			f_4005[i] = f[i] * nd40_1nm[i] * nd05_1nm[i]
		f_4006 = np.empty(401)
		for i in range(401):
			f_4006[i] = f[i] * nd40_1nm[i] * nd06_1nm[i]
		f_4007 = np.empty(401)
		for i in range(401):
			f_4007[i] = f[i] * nd40_1nm[i] * nd07_1nm[i]
		f_4008 = np.empty(401)
		for i in range(401):
			f_4008[i] = f[i] * nd40_1nm[i] * nd08_1nm[i]
		f_4009 = np.empty(401)
		for i in range(401):
			f_4009[i] = f[i] * nd40_1nm[i] * nd09_1nm[i]
		
		print("Brightness contrast with blue:")
		print("0.1: " + str(brightness_contrast(f_01, blue47_1nm)))
		print("0.2: " + str(brightness_contrast(f_02, blue47_1nm)))
		print("0.3: " + str(brightness_contrast(f_03, blue47_1nm)))
		print("0.4: " + str(brightness_contrast(f_04, blue47_1nm)))
		print("0.5: " + str(brightness_contrast(f_05, blue47_1nm)))
		print("0.6: " + str(brightness_contrast(f_06, blue47_1nm)))
		print("0.7: " + str(brightness_contrast(f_07, blue47_1nm)))
		print("0.8: " + str(brightness_contrast(f_08, blue47_1nm)))
		print("0.9: " + str(brightness_contrast(f_09, blue47_1nm)))
		print("1.0: " + str(brightness_contrast(f_10, blue47_1nm)))
		print("1.0 + 0.1: " + str(brightness_contrast(f_1001, blue47_1nm)))
		print("1.0 + 0.2: " + str(brightness_contrast(f_1002, blue47_1nm)))
		print("1.0 + 0.3: " + str(brightness_contrast(f_1003, blue47_1nm)))
		print("1.0 + 0.4: " + str(brightness_contrast(f_1004, blue47_1nm)))
		print("1.0 + 0.5: " + str(brightness_contrast(f_1005, blue47_1nm)))
		print("1.0 + 0.6: " + str(brightness_contrast(f_1006, blue47_1nm)))
		print("1.0 + 0.7: " + str(brightness_contrast(f_1007, blue47_1nm)))
		print("1.0 + 0.8: " + str(brightness_contrast(f_1008, blue47_1nm)))
		print("1.0 + 0.9: " + str(brightness_contrast(f_1009, blue47_1nm)))
		print("2.0: " + str(brightness_contrast(f_20, blue47_1nm)))
		print("2.0 + 0.1: " + str(brightness_contrast(f_2001, blue47_1nm)))
		print("2.0 + 0.2: " + str(brightness_contrast(f_2002, blue47_1nm)))
		print("2.0 + 0.3: " + str(brightness_contrast(f_2003, blue47_1nm)))
		print("2.0 + 0.4: " + str(brightness_contrast(f_2004, blue47_1nm)))
		print("2.0 + 0.5: " + str(brightness_contrast(f_2005, blue47_1nm)))
		print("2.0 + 0.6: " + str(brightness_contrast(f_2006, blue47_1nm)))
		print("2.0 + 0.7: " + str(brightness_contrast(f_2007, blue47_1nm)))
		print("2.0 + 0.8: " + str(brightness_contrast(f_2008, blue47_1nm)))
		print("2.0 + 0.9: " + str(brightness_contrast(f_2009, blue47_1nm)))
		print("3.0: " + str(brightness_contrast(f_30, blue47_1nm)))
		print("3.0 + 0.1: " + str(brightness_contrast(f_3001, blue47_1nm)))
		print("3.0 + 0.2: " + str(brightness_contrast(f_3002, blue47_1nm)))
		print("3.0 + 0.3: " + str(brightness_contrast(f_3003, blue47_1nm)))
		print("3.0 + 0.4: " + str(brightness_contrast(f_3004, blue47_1nm)))
		print("3.0 + 0.5: " + str(brightness_contrast(f_3005, blue47_1nm)))
		print("3.0 + 0.6: " + str(brightness_contrast(f_3006, blue47_1nm)))
		print("3.0 + 0.7: " + str(brightness_contrast(f_3007, blue47_1nm)))
		print("3.0 + 0.8: " + str(brightness_contrast(f_3008, blue47_1nm)))
		print("3.0 + 0.9: " + str(brightness_contrast(f_3009, blue47_1nm)))
		print("4.0: " + str(brightness_contrast(f_40, blue47_1nm)))
		print("4.0 + 0.1: " + str(brightness_contrast(f_4001, blue47_1nm)))
		print("4.0 + 0.2: " + str(brightness_contrast(f_4002, blue47_1nm)))
		print("4.0 + 0.3: " + str(brightness_contrast(f_4003, blue47_1nm)))
		print("4.0 + 0.4: " + str(brightness_contrast(f_4004, blue47_1nm)))
		print("4.0 + 0.5: " + str(brightness_contrast(f_4005, blue47_1nm)))
		print("4.0 + 0.6: " + str(brightness_contrast(f_4006, blue47_1nm)))
		print("4.0 + 0.7: " + str(brightness_contrast(f_4007, blue47_1nm)))
		print("4.0 + 0.8: " + str(brightness_contrast(f_4008, blue47_1nm)))
		print("4.0 + 0.9: " + str(brightness_contrast(f_4009, blue47_1nm)))
	
	# red
	print("Red 25")
	brightness_tests(red25_1nm)
	
	# yellow
	print("")
	print("Yellow 15")
	brightness_tests(yellow15_1nm)
	
	# green
	print("")
	print("Green 58")
	brightness_tests(green58_1nm)
	
	# gray
	print("")
	print("Gray")
	brightness_tests(placeholder)

# These are used by both --blackbody and --kodak, so we keep them out here.
# red 25
red25_0 = np.empty(401)
for i in range(401):
	if (args.kv2): red25_0[i] = red25_1nm[i] * nd10_1nm[i] * nd03_1nm[i]
	else: red25_0[i] = red25_1nm[i] * nd10_1nm[i] * nd05_1nm[i]
red25_03 = np.empty(401)
for i in range(401):
	red25_03[i] = red25_0[i] * nd03_1nm[i]
red25_07 = np.empty(401)
for i in range(401):
	red25_07[i] = red25_0[i] * nd07_1nm[i]
red25_10 = np.empty(401)
for i in range(401):
	red25_10[i] = red25_0[i] * nd10_1nm[i]

# yellow 15
yellow15_0 = np.empty(401)
for i in range(401):
	if (args.kv2): yellow15_0[i] = yellow15_1nm[i] * nd10_1nm[i] * nd09_1nm[i]
	else: yellow15_0[i] = yellow15_1nm[i] * nd20_1nm[i]
yellow15_03 = np.empty(401)
for i in range(401):
	yellow15_03[i] = yellow15_0[i] * nd03_1nm[i]
yellow15_07 = np.empty(401)
for i in range(401):
	yellow15_07[i] = yellow15_0[i] * nd07_1nm[i]
yellow15_10 = np.empty(401)
for i in range(401):
	yellow15_10[i] = yellow15_0[i] * nd10_1nm[i]

# green 58
green58_0 = np.empty(401)
for i in range(401):
	if (args.kv2): green58_0[i] = green58_1nm[i] * nd10_1nm[i] * nd02_1nm[i]
	else: green58_0[i] = green58_1nm[i] * nd10_1nm[i] * nd03_1nm[i]
green58_03 = np.empty(401)
for i in range(401):
	green58_03[i] = green58_0[i] * nd03_1nm[i]
green58_07 = np.empty(401)
for i in range(401):
	green58_07[i] = green58_0[i] * nd07_1nm[i]
green58_10 = np.empty(401)
for i in range(401):
	green58_10[i] = green58_0[i] * nd10_1nm[i]

# blue 47
blue47_0 = blue47_1nm # keep this name for convenience
blue47_03 = np.empty(401)
for i in range(401):
	blue47_03[i] = blue47_0[i] * nd03_1nm[i]
blue47_07 = np.empty(401)
for i in range(401):
	blue47_07[i] = blue47_0[i] * nd07_1nm[i]
blue47_10 = np.empty(401)
for i in range(401):
	blue47_10[i] = blue47_0[i] * nd10_1nm[i]

# gray
gray_0 = np.empty(401)
for i in range(401):
	if (args.kv2): gray_0[i] = nd20_1nm[i]
	else: gray_0[i] = nd20_1nm[i] * nd01_1nm[i]
gray_03 = np.empty(401)
for i in range(401):
	gray_03[i] = gray_0[i] * nd03_1nm[i]
gray_07 = np.empty(401)
for i in range(401):
	gray_07[i] = gray_0[i] * nd07_1nm[i]
gray_10 = np.empty(401)
for i in range(401):
	gray_10[i] = gray_0[i] * nd10_1nm[i]

# used by both --kodak and --munsell
# This outputs a probability for both success (defined as at least the criterion) and failure (less
# than the criterion) because the behavioral results include both. We use different thresholds
# for "success" for the different sets of discriminations. In Friedman (1967) the criterion was
# 35 out of 40 for the color pairs and 68 of 80 for the color/gray pairs. In Gutierrez et al. (2011)
# significant performance is defined as more than 62.5% and the total number of trials for each
# pair was 64, so the criterion would be 41 of 64. This produces considerably different
# significance thresholds when comparing "success" or "failure" against expected values:
# * 35/40 and 68/80: success 0-12 (out of 16), failure 15-16
# * 41/64: success 0-8, failure 12-16
# For contrast where we don't consider overlap and expect chance (50%) performance for deltaS < 1,
# this means that under the first criterion "success" is defined as at least 10 of 16 distinguishable
# pairs and "failure" is defined as at most 13 of 16, whereas under the second criterion "success"
# is at least 1 and "failure" is at most 7. There are also several intermediate values where
# the probability of both success and failure exceeds 0.05. With the first two criteria, success
# becomes more likely than failure at 14 (12 distinguishable), and with the second this occurs at 10.5
# (5 distinguishable).
# This is now also used directly for optimization graphs, so all the parts relevant only to
# testing one set of pigments have been made conditional to reduce noise and computation.
def color_disc(f0, f1, f2, f3, s0, s1, s2, s3, correct=35, trials=40, order=0, alternative='greater', lw = args.lw, mw = args.mw, sw = args.sw, output=True):
	d = 0 # different
	s = 0 # same
	counter = 0
	
	for i in range(0, 4):
		if (i == 0):
			flevel = f0
		elif (i == 1):
			flevel = f1
		elif (i == 2):
			flevel = f2
		elif (i == 3):
			flevel = f3
		for j in range(0, 4):
			if (j == 0):
				slevel = s0
			elif (j == 1):
				slevel = s1
			elif (j == 2):
				slevel = s2
			elif (j == 3):
				slevel = s3
			contrast = color_contrast(flevel, slevel, lw = lw, mw = mw, sw = sw)
			box[counter][order] = contrast
			#print(box[counter][order])
			if (output):
				print(str(i) + " vs. " + str(j) + ": " + str(contrast))
				if (contrast < 1):
					s += 1
				else:
					d += 1
			counter += 1
	if (output):
		print("Different: " + str(d))
		print("Same: " + str(s))
		#if (l1 == m1):
		binomial = binomtest(correct, trials, p=(d + s/2)/16, alternative='greater')
		binomial2 = binomtest(correct-1, trials, p=(d + s/2)/16, alternative='less')
		#else:
		#	binomial = binomtest(correct, trials, p=(d + s/2)/16, alternative=alternative)
		print("P-value (success): " + str(binomial.pvalue))
		print("P-value (failure): " + str(binomial2.pvalue))
	
	# median contrast
	#print(box)
	#print([box[0][order], box[1][order], box[2][order], box[3][order], box[4][order], box[5][order], box[6][order], box[7][order], box[8][order], box[9][order], box[10][order], box[11][order], box[12][order], box[13][order], box[14][order], box[15][order]])
	mc = statistics.median([box[0][order], box[1][order], box[2][order], box[3][order], box[4][order], box[5][order], box[6][order], box[7][order], box[8][order], box[9][order], box[10][order], box[11][order], box[12][order], box[13][order], box[14][order], box[15][order]])
	minc = min([box[0][order], box[1][order], box[2][order], box[3][order], box[4][order], box[5][order], box[6][order], box[7][order], box[8][order], box[9][order], box[10][order], box[11][order], box[12][order], box[13][order], box[14][order], box[15][order]])
	maxc = max([box[0][order], box[1][order], box[2][order], box[3][order], box[4][order], box[5][order], box[6][order], box[7][order], box[8][order], box[9][order], box[10][order], box[11][order], box[12][order], box[13][order], box[14][order], box[15][order]])
	if (output):
		print("Median contrast: " + str(mc))
		print("")
		#return
	return [mc, minc, maxc]

if (args.kodak):
	
	# plot
	xvalues = np.empty(401)
	for i in range(401):
		xvalues[i] = i + 300
	plt.subplot(3, 2, 1)
	plt.plot(xvalues, red25_0*100, 'r')
	plt.plot(xvalues, red25_03*100, 'r')
	plt.plot(xvalues, red25_07*100, 'r')
	plt.plot(xvalues, red25_10*100, 'r')
	plt.title("Red 25")
	plt.subplot(3, 2, 2)
	plt.plot(xvalues, yellow15_0*100, 'y')
	plt.plot(xvalues, yellow15_03*100, 'y')
	plt.plot(xvalues, yellow15_07*100, 'y')
	plt.plot(xvalues, yellow15_10*100, 'y')
	plt.title("Yellow 15")
	plt.subplot(3, 2, 3)
	plt.plot(xvalues, green58_0*100, 'g')
	plt.plot(xvalues, green58_03*100, 'g')
	plt.plot(xvalues, green58_07*100, 'g')
	plt.plot(xvalues, green58_10*100, 'g')
	plt.title("Green 58")
	plt.subplot(3, 2, 4)
	plt.plot(xvalues, blue47_0*100, 'b')
	plt.plot(xvalues, blue47_03*100, 'b')
	plt.plot(xvalues, blue47_07*100, 'b')
	plt.plot(xvalues, blue47_10*100, 'b')
	plt.title("Blue 47")
	plt.subplot(3, 2, 5)
	plt.plot(xvalues, gray_0*100, color='gray')
	plt.plot(xvalues, gray_03*100, color='gray')
	plt.plot(xvalues, gray_07*100, color='gray')
	plt.plot(xvalues, gray_10*100, color='gray')
	plt.title("Gray")
	plt.show()
	
	# brightness levels and color space coordinates
	red25l = np.empty(4)
	red25cs = np.empty(4)
	red25cs1 = np.empty(4)
	print("Red 25 (+ 1.0 + 0.1)")
	red25data = spectral_rendering(red25_0)
	red25l[0] = red25data[0]
	red25cs[0] = red25data[1]
	red25cs1[0] = red25data[2]
	print("+ 0.3")
	red25data = spectral_rendering(red25_03)
	red25l[1] = red25data[0]
	red25cs[1] = red25data[1]
	red25cs1[1] = red25data[2]
	print("+ 0.7")
	red25data = spectral_rendering(red25_07)
	red25l[2] = red25data[0]
	red25cs[2] = red25data[1]
	red25cs1[2] = red25data[2]
	print("+ 1.0")
	red25data = spectral_rendering(red25_10)
	red25l[3] = red25data[0]
	red25cs[3] = red25data[1]
	red25cs1[3] = red25data[2]
	
	yellow15l = np.empty(4)
	yellow15cs = np.empty(4)
	yellow15cs1 = np.empty(4)
	print("Yellow 15 (+ 2.0 + 0.5)")
	yellow15data = spectral_rendering(yellow15_0)
	yellow15l[0] = yellow15data[0]
	yellow15cs[0] = yellow15data[1]
	yellow15cs1[0] = yellow15data[2]
	print("+ 0.3")
	yellow15data = spectral_rendering(yellow15_03)
	yellow15l[1] = yellow15data[0]
	yellow15cs[1] = yellow15data[1]
	yellow15cs1[1] = yellow15data[2]
	print("+ 0.7")
	yellow15data = spectral_rendering(yellow15_07)
	yellow15l[2] = yellow15data[0]
	yellow15cs[2] = yellow15data[1]
	yellow15cs1[2] = yellow15data[2]
	print("+ 1.0")
	yellow15data = spectral_rendering(yellow15_10)
	yellow15l[3] = yellow15data[0]
	yellow15cs[3] = yellow15data[1]
	yellow15cs1[3] = yellow15data[2]
	
	green58l = np.empty(4)
	green58cs = np.empty(4)
	green58cs1 = np.empty(4)
	print("Green 58 (+ 1.0 + 0.3)")
	green58data = spectral_rendering(green58_0)
	green58l[0] = green58data[0]
	green58cs[0] = green58data[1]
	green58cs1[0] = green58data[2]
	print("+ 0.3")
	green58data = spectral_rendering(green58_03)
	green58l[1] = green58data[0]
	green58cs[1] = green58data[1]
	green58cs1[1] = green58data[2]
	print("+ 0.7")
	green58data = spectral_rendering(green58_07)
	green58l[2] = green58data[0]
	green58cs[2] = green58data[1]
	green58cs1[2] = green58data[2]
	print("+ 1.0")
	green58data = spectral_rendering(green58_10)
	green58l[3] = green58data[0]
	green58cs[3] = green58data[1]
	green58cs1[3] = green58data[2]
	
	blue47l = np.empty(4)
	blue47cs = np.empty(4)
	blue47cs1 = np.empty(4)
	print("Blue 47")
	blue47data = spectral_rendering(blue47_0)
	blue47l[0] = blue47data[0]
	blue47cs[0] = blue47data[1]
	blue47cs1[0] = blue47data[2]
	print("+ 0.3")
	blue47data = spectral_rendering(blue47_03)
	blue47l[1] = blue47data[0]
	blue47cs[1] = blue47data[1]
	blue47cs1[1] = blue47data[2]
	print("+ 0.7")
	blue47data = spectral_rendering(blue47_07)
	blue47l[2] = blue47data[0]
	blue47cs[2] = blue47data[1]
	blue47cs1[2] = blue47data[2]
	print("+ 1.0")
	blue47data = spectral_rendering(blue47_10)
	blue47l[3] = blue47data[0]
	blue47cs[3] = blue47data[1]
	blue47cs1[3] = blue47data[2]
	
	grayl = np.empty(4)
	graycs = np.empty(4)
	graycs1 = np.empty(4)
	print("Gray")
	graydata = spectral_rendering(gray_0)
	grayl[0] = graydata[0]
	graycs[0] = graydata[1]
	graycs1[0] = graydata[2]
	print("+ 0.3")
	graydata = spectral_rendering(gray_03)
	grayl[1] = graydata[0]
	graycs[1] = graydata[1]
	graycs1[1] = graydata[2]
	print("+ 0.7")
	graydata = spectral_rendering(gray_07)
	grayl[2] = graydata[0]
	graycs[2] = graydata[1]
	graycs1[2] = graydata[2]
	print("+ 1.0")
	graydata = spectral_rendering(gray_10)
	grayl[3] = graydata[0]
	graycs[3] = graydata[1]
	graycs1[3] = graydata[2]
	
	# brightness differences between colors
	# We want to know both the direction and the strength of the difference
	# to judge whether it could be reliably used as a cue.
	
	# function
	def brightness_disc(f0, f1, f2, f3, fl, s0, s1, s2, s3, sl, correct=35, trials=40): 
		fb = 0 # first brighter
		sb = 0 # second brighter
		equal = 0
		# first set
		for i in range(0, 4):
			if (i == 0):
				flevel = f0
			elif (i == 1):
				flevel = f1
			elif (i == 2):
				flevel = f2
			elif (i == 3):
				flevel = f3
			# second set
			for j in range(0, 4):
				if (j == 0):
					slevel = s0
				elif (j == 1):
					slevel = s1
				elif (j == 2):
					slevel = s2
				elif (j == 3):
					slevel = s3
				diff = fl[i] - sl[j]
				contrast = brightness_contrast(flevel, slevel)
				if (contrast >= 1):
					if (diff > 0):
						fb += 1
					elif (diff < 0):
						sb += 1
				else:
					equal += 1
		print("First brighter: " + str(fb))
		print("Second brighter: " + str(sb))
		print("Equal: " + str(equal))
		print("Expected % correct: " + str(100*(max(fb, sb) + equal/2) / 16) + "% (" + str(100*max(fb, sb) / 16) + "-" + str(100*(max(fb, sb) + equal) / 16) + "%)")
		binomial = binomtest(correct, trials, (max(fb, sb) + equal/2) / 16, alternative='greater')
		binomial2 = binomtest(correct-1, trials, (max(fb, sb) + equal/2) / 16, alternative='less')
		print("P-value (success): " + str(binomial.pvalue))
		print("P-value (failure): " + str(binomial2.pvalue))
		print("")
	
	# red-yellow
	print("R-Y")
	brightness_disc(red25_0, red25_03, red25_07, red25_10, red25l, yellow15_0, yellow15_03, yellow15_07, yellow15_10, yellow15l)
	
	# red-green
	print("R-G")
	brightness_disc(red25_0, red25_03, red25_07, red25_10, red25l, green58_0, green58_03, green58_07, green58_10, green58l)
	
	# red-blue
	print("R-B")
	brightness_disc(red25_0, red25_03, red25_07, red25_10, red25l, blue47_0, blue47_03, blue47_07, blue47_10, blue47l)
	
	# yellow-green
	print("Y-G")
	brightness_disc(yellow15_0, yellow15_03, yellow15_07, yellow15_10, yellow15l, green58_0, green58_03, green58_07, green58_10, green58l)
	
	# yellow-blue
	print("Y-B")
	brightness_disc(yellow15_0, yellow15_03, yellow15_07, yellow15_10, yellow15l, blue47_0, blue47_03, blue47_07, blue47_10, blue47l)
	
	# green-blue
	print("G-B")
	brightness_disc(green58_0, green58_03, green58_07, green58_10, green58l, blue47_0, blue47_03, blue47_07, blue47_10, blue47l)
	
	# colors vs. gray
	
	# red
	print("Red vs. gray")
	brightness_disc(red25_0, red25_03, red25_07, red25_10, red25l, gray_0, gray_03, gray_07, gray_10, grayl, correct=68, trials=80)
	
	# yellow
	print("Yellow vs. gray")
	brightness_disc(yellow15_0, yellow15_03, yellow15_07, yellow15_10, yellow15l, gray_0, gray_03, gray_07, gray_10, grayl, correct=68, trials=80)
	
	# green
	print("Green vs. gray")
	brightness_disc(green58_0, green58_03, green58_07, green58_10, green58l, gray_0, gray_03, gray_07, gray_10, grayl, correct=68, trials=80)
	
	# blue
	print("Blue vs. gray")
	brightness_disc(blue47_0, blue47_03, blue47_07, blue47_10, blue47l, gray_0, gray_03, gray_07, gray_10, grayl, correct=68, trials=80)
	
	# plot
	x = np.array([0.0, 0.3, 0.7, 1.0])
	plt.plot(x, red25l, 'sr', mec='k', label="red 25")
	plt.plot(x, yellow15l, 'Dy', mec='k', label="yellow 15")
	plt.plot(x, green58l, '^g', mec='k', label="green 58")
	plt.plot(x, blue47l, 'ob', mec='k', label="blue 47")
	plt.plot(x, grayl, marker='v', linestyle='', color='gray', mec='k', label="gray")
	plt.xlabel("Filter optical density")
	plt.ylabel("Relative quantum catch")
	plt.yscale('log')
	plt.legend()
	plt.show()
	
	x = np.array(["R-Y", "R-G", "R-B", "Y-G", "Y-B", "G-B"])
	
	# color differences -- this is the hard part...
	
	# first we plot them in a two-dimensional color space with (L-S)/(L+S) on
	# the x-axis and brightness on the y-axis. Not sure how much this tells us
	# though.
	if (l1 == m1):
		plt.plot(red25cs, red25l, 'sr', mec='k', label="red 25")
		plt.plot(yellow15cs, yellow15l, 'Dy', mec='k', label="yellow 15")
		plt.plot(green58cs, green58l, '^g', mec='k', label="green 58")
		plt.plot(blue47cs, blue47l, 'ob', mec='k', label="blue 47")
		plt.plot(graycs, grayl, marker='v', linestyle='', color='gray', mec='k', label="gray")
		plt.xlabel("Chromaticity ((L-S)/(L+S))")
		plt.ylabel("Brightness (L)")
		plt.legend()
		plt.show()
	else:
		plt.plot(red25cs, red25cs1, 'sr', mec='k', label="red 25")
		plt.plot(yellow15cs, yellow15cs1, 'Dy', mec='k', label="yellow 15")
		plt.plot(green58cs, green58cs1, '^g', mec='k', label="green 58")
		plt.plot(blue47cs, blue47cs1, 'ob', mec='k', label="blue 47")
		plt.plot(graycs, graycs1, marker='v', linestyle='', color='gray', mec='k', label="gray")
		plt.xlabel("(L-S)/(L+M+S)")
		plt.ylabel("(M-0.5(L+S))/(L+M+S)")
		xborder = np.array([-math.sqrt(1/2), 0, math.sqrt(1/2), -math.sqrt(1/2)])
		yborder = np.array([-math.sqrt(2/3)/2, math.sqrt(2/3), -math.sqrt(2/3)/2, -math.sqrt(2/3)/2])
		plt.plot(xborder, yborder, '-k')
		plt.text(-math.sqrt(1/2) - 0.05, -math.sqrt(2/3)/2 - 0.025, 'S')
		plt.text(0 - 0.025, math.sqrt(2/3) + 0.0125, 'M')
		plt.text(math.sqrt(1/2) + 0.0125, -math.sqrt(2/3)/2 - 0.025, 'L')
		plt.legend()
		plt.show()
	
	# Next we try to assess whether they're distinguishable. As with brightness, we
	# check both the contrast and the direction.
	
	box = np.empty((16, 6))
	
	# R-Y
	print("R-Y")
	color_disc(red25_0, red25_03, red25_07, red25_10, yellow15_0, yellow15_03, yellow15_07, yellow15_10, order=0)
	
	# R-G
	print("R-G")
	color_disc(red25_0, red25_03, red25_07, red25_10, green58_0, green58_03, green58_07, green58_10, order=1)
	
	# R-B
	print("R-B")
	color_disc(red25_0, red25_03, red25_07, red25_10, blue47_0, blue47_03, blue47_07, blue47_10, order=2)
	
	# Y-G
	print("Y-G")
	color_disc(yellow15_0, yellow15_03, yellow15_07, yellow15_10, green58_0, green58_03, green58_07, green58_10, order=3)
	
	# Y-B
	print("Y-B")
	color_disc(yellow15_0, yellow15_03, yellow15_07, yellow15_10, blue47_0, blue47_03, blue47_07, blue47_10, order=4)
	
	# G-B
	print("G-B")
	color_disc(green58_0, green58_03, green58_07, green58_10, blue47_0, blue47_03, blue47_07, blue47_10, order=5)
	
	# medians
	mry = statistics.median([box[0][0], box[1][0], box[2][0], box[3][0], box[4][0], box[5][0], box[6][0], box[7][0], box[8][0], box[9][0], box[10][0], box[11][0], box[12][0], box[13][0], box[14][0], box[15][0]])
	mrg = statistics.median([box[0][1], box[1][1], box[2][1], box[3][1], box[4][1], box[5][1], box[6][1], box[7][1], box[8][1], box[9][1], box[10][1], box[11][1], box[12][1], box[13][1], box[14][1], box[15][1]])
	mrb = statistics.median([box[0][2], box[1][2], box[2][2], box[3][2], box[4][2], box[5][2], box[6][2], box[7][2], box[8][2], box[9][2], box[10][2], box[11][2], box[12][2], box[13][2], box[14][2], box[15][2]])
	myg = statistics.median([box[0][3], box[1][3], box[2][3], box[3][3], box[4][3], box[5][3], box[6][3], box[7][3], box[8][3], box[9][3], box[10][3], box[11][3], box[12][3], box[13][3], box[14][3], box[15][3]])
	myb = statistics.median([box[0][4], box[1][4], box[2][4], box[3][4], box[4][4], box[5][4], box[6][4], box[7][4], box[8][4], box[9][4], box[10][4], box[11][4], box[12][4], box[13][4], box[14][4], box[15][4]])
	mgb = statistics.median([box[0][5], box[1][5], box[2][5], box[3][5], box[4][5], box[5][5], box[6][5], box[7][5], box[8][5], box[9][5], box[10][5], box[11][5], box[12][5], box[13][5], box[14][5], box[15][5]])
	
	print("R-Y median contrast: " + str(mry))
	print("R-G median contrast: " + str(mrg))
	print("R-B median contrast: " + str(mrb))
	print("Y-G median contrast: " + str(myg))
	print("Y-B median contrast: " + str(myb))
	print("G-B median contrast: " + str(mgb))
	print("Lowest median contrast: " + str(min(mry, mrg, mrb, myg, myb, mgb)))
	print("Highest median contrast: " + str(max(mry, mrg, mrb, myg, myb, mgb)))
	print("")
	
	# plot contrast values
	plt.boxplot(x=box, tick_labels=x)
	plt.ylabel("ΔS (JND)")
	plt.show()
	
	# and again for colors vs gray. Is gray distinguishable from yellow and green? If
	# not, we have a major problem.
	box = np.empty((16, 4))
	labels = ["red", "yellow", "green", "blue"]
	
	# red
	print("Red vs. gray")
	color_disc(red25_0, red25_03, red25_07, red25_10, gray_0, gray_03, gray_07, gray_10, order=0, correct=68, trials=80)
	
	# yellow
	print("Yellow vs. gray")
	color_disc(yellow15_0, yellow15_03, yellow15_07, yellow15_10, gray_0, gray_03, gray_07, gray_10, correct=68, trials=80, order=1)
	
	# green
	print("Green vs. gray")
	color_disc(green58_0, green58_03, green58_07, green58_10, gray_0, gray_03, gray_07, gray_10, correct=68, trials=80, order=2)
	
	# blue
	print("Blue vs. gray")
	color_disc(blue47_0, blue47_03, blue47_07, blue47_10, gray_0, gray_03, gray_07, gray_10, correct=68, trials=80, order=3)
	
	plt.boxplot(x=box, tick_labels=labels)
	plt.ylabel("ΔS (JND)")
	plt.show()
	# Yes, it is! This came as a bit of a surprise considering the rods and cones apparently
	# have very similar responses to these.

if (args.blackbody != 0):
	energy = sigma * args.blackbody**4
	surface_area = 4*math.pi*args.radius**2
	print("Surface area of sphere: " + str(sphere_area))
	print("Energy produced by light bulb (W/cm^2): " + str(watts_cm2))
	print("Energy produced by light bulb (W/m^2): " + str(watts_m2))
	#print("Surface area of light bulb (cm^2): " + str(surface_area))
	print("Energy produced by blackbody model (W/m^2, 4pi sr): " + str(energy))
	print("Energy produced by blackbody model (W): " + str(energy * surface_area / 10000))
	print("W/m^2 scaling factor: " + str(scale))
	
	# how much energy is produced between 300-800 nm when paired with blue 47
	# Fix this -- should use the "full" integral with 1-nm steps. As written, the values
	# produced are about 1/10 what they should be. Also try to convert it to cd/m^2 so
	# we can match the value given in the study of 5.5 foot-lamberts ~= 19 cd/m^2 or
	# 59 lx (foot-lamberts are equivalent to both, and as above, the conversion factor
	# between radiance/luminance and irradiance/illuminance is pi sr). Candelas aren't
	# any less annoying than foot-lamberts just because they think they're an SI unit.
	
	# According to the new numbers, the bulbs would have to have a color temperature of
	# about 1700 K to actually produce the brightness reported in the study. I don't
	# think they make those. Alternatively, if the number corresponds to illuminance,
	# this would be scaled by steradians and the true value in cd/m^2 would be ~131,
	# which corresponds to a color temperature of between 2000 and 2100 K. This is still
	# much lower than a typical incandescent bulb. If I change the scaling factor back to
	# W/m^2 instead of W/m^2/sr, however, the illucd/mminance definition places the color
	# temperature at about 2650 K, which is close to the number 2700 I found
	# when searching for GE 656 bulbs (https://www.bulbs.com/product/656).
	
	# Now that I've fixed the density/transmission issues with the filters, the
	# filtered luminance is around 19 cd/m^2 for a color temperature near 2856 K.
	# I'm not sure I buy this because the Macbeth Illuminometer is supposed to measure
	# illuminance in foot-candles rather than luminance in "foot-lamberts". At least
	# the candela only has one definition. Also my quantum catch model now predicts
	# that neither humans nor opossums can reliably distinguish the R-Y pair unless
	# I assume a high degree of spatial summation. Since lux/foot-candles is measured
	# by light falling on some flat surface, I don't know what distinguishes light
	# falling only on that surface from light radiated in all directions from the
	# source. A "perfect" reflecting surface radiates 1/pi of the illuminance it
	# receives, but what does that say about the radiance/luminance of the light
	# source?
	
	# We want to scale the directional intensity. There is technically less irradiance
	# if spread over a larger area, but the radiance doesn't change. See the
	# Stefan-Boltzmann law -- "power per unit area" is power per the area of the
	# object doing the radiating, not power per how many square meters it's spread
	# over. We need an estimate of the light bulb size. I used to have that in here
	# but discarded it for some reason. The distance from the light source shouldn't
	# be involved for the same reason we don't use it in the lux calculation.
	# Actually I think this is wrong because a light bulb that's farther away from the
	# viewing point is basically the same as a bigger bulb. If we use the value of
	# 5 cm and define the number of cones per receptive field to include at least 1
	# of each type, we can nearly perfectly explain both the reported brightness
	# (assuming cd/m^2) and the opossums' behavior, though R-Y is apparently a very
	# difficult discrimination.
	
	radiance = quad(blackbody, 300, 700, args=args.blackbody)
	visible = radiance[0]
	x10nm = np.empty(41)
	for i in range(41):
		x10nm[i] = i*10 + 300
	x1nm = np.empty(501)
	for i in range(501):
		x1nm[i] = i + 300
	
	visible_r0 = 0
	for i in range(300, 701):
		visible_r0 += blackbody(i, args.blackbody) * red25_0[i-300]
	visible_r03 = 0
	for i in range(300, 701):
		visible_r03 += blackbody(i, args.blackbody) * red25_03[i-300]
	visible_r07 = 0
	for i in range(300, 701):
		visible_r07 += blackbody(i, args.blackbody) * red25_07[i-300]
	visible_r10 = 0
	for i in range(300, 701):
		visible_r10 += blackbody(i, args.blackbody) * red25_10[i-300]
	visible_y0 = 0
	for i in range(300, 701):
		visible_y0 += blackbody(i, args.blackbody) * yellow15_0[i-300]
	visible_y03 = 0
	for i in range(300, 701):
		visible_y03 += blackbody(i, args.blackbody) * yellow15_03[i-300]
	visible_y07 = 0
	for i in range(300, 701):
		visible_y07 += blackbody(i, args.blackbody) * yellow15_07[i-300]
	visible_y10 = 0
	for i in range(300, 701):
		visible_y10 += blackbody(i, args.blackbody) * yellow15_10[i-300]
	visible_g0 = 0
	for i in range(300, 701):
		visible_g0 += blackbody(i, args.blackbody) * green58_0[i-300]
	visible_g03 = 0
	for i in range(300, 701):
		visible_g03 += blackbody(i, args.blackbody) * green58_03[i-300]
	visible_g07 = 0
	for i in range(300, 701):
		visible_g07 += blackbody(i, args.blackbody) * green58_07[i-300]
	visible_g10 = 0
	for i in range(300, 701):
		visible_g10 += blackbody(i, args.blackbody) * green58_10[i-300]
	visible_b0 = 0
	for i in range(300, 701):
		visible_b0 += blackbody(i, args.blackbody) * blue47_0[i-300]
	visible_b03 = 0
	for i in range(300, 701):
		visible_b03 += blackbody(i, args.blackbody) * blue47_03[i-300]
	visible_b07 = 0
	for i in range(300, 701):
		visible_b07 += blackbody(i, args.blackbody) * blue47_07[i-300]
	visible_b10 = 0
	for i in range(300, 701):
		visible_b10 += blackbody(i, args.blackbody) * blue47_10[i-300]
	visible_gray0 = 0
	for i in range(300, 701):
		visible_gray0 += blackbody(i, args.blackbody) * gray_0[i-300]
	visible_gray03 = 0
	for i in range(300, 701):
		visible_gray03 += blackbody(i, args.blackbody) * gray_03[i-300]
	visible_gray07 = 0
	for i in range(300, 701):
		visible_gray07 += blackbody(i, args.blackbody) * gray_07[i-300]
	visible_gray10 = 0
	for i in range(300, 701):
		visible_gray10 += blackbody(i, args.blackbody) * gray_10[i-300]
	
	# lx/cd
	visible_bcd0 = 0
	for i in range(300, 701):
		visible_bcd0 += blackbody(i, args.blackbody) * blue47_0[i-300] *  luminosity[i-300] * 683.002
	visible_bcd03 = 0
	for i in range(300, 701):
		visible_bcd03 += blackbody(i, args.blackbody) * blue47_03[i-300] *  luminosity[i-300] * 683.002
	visible_bcd07 = 0
	for i in range(300, 701):
		visible_bcd07 += blackbody(i, args.blackbody) * blue47_07[i-300] *  luminosity[i-300] * 683.002
	visible_bcd10 = 0
	for i in range(300, 701):
		visible_bcd10 += blackbody(i, args.blackbody) * blue47_10[i-300] *  luminosity[i-300] * 683.002
	#visible_filtered_cd *= luminance_full / luminance_10nm
	visible_lx0 = visible_bcd0 * math.pi
	visible_lx03 = visible_bcd03 * math.pi
	visible_lx07 = visible_bcd07 * math.pi
	visible_lx10 = visible_bcd10 * math.pi
	# other colors
	# red
	visible_rcd0 = 0
	for i in range(300, 701):
		visible_rcd0 += blackbody(i, args.blackbody) * red25_0[i-300] *  luminosity[i-300] * 683.002
	visible_rcd03 = 0
	for i in range(300, 701):
		visible_rcd03 += blackbody(i, args.blackbody) * red25_03[i-300] *  luminosity[i-300] * 683.002
	visible_rcd07 = 0
	for i in range(300, 701):
		visible_rcd07 += blackbody(i, args.blackbody) * red25_07[i-300] *  luminosity[i-300] * 683.002
	visible_rcd10 = 0
	for i in range(300, 701):
		visible_rcd10 += blackbody(i, args.blackbody) * red25_10[i-300] *  luminosity[i-300] * 683.002
	# yellow
	visible_ycd0 = 0
	for i in range(300, 701):
		visible_ycd0 += blackbody(i, args.blackbody) * yellow15_0[i-300] *  luminosity[i-300] * 683.002
	visible_ycd03 = 0
	for i in range(300, 701):
		visible_ycd03 += blackbody(i, args.blackbody) * yellow15_03[i-300] *  luminosity[i-300] * 683.002
	visible_ycd07 = 0
	for i in range(300, 701):
		visible_ycd07 += blackbody(i, args.blackbody) * yellow15_07[i-300] *  luminosity[i-300] * 683.002
	visible_ycd10 = 0
	for i in range(300, 701):
		visible_ycd10 += blackbody(i, args.blackbody) * yellow15_10[i-300] *  luminosity[i-300] * 683.002
	# green
	visible_gcd0 = 0
	for i in range(300, 701):
		visible_gcd0 += blackbody(i, args.blackbody) * green58_0[i-300] *  luminosity[i-300] * 683.002
	visible_gcd03 = 0
	for i in range(300, 701):
		visible_gcd03 += blackbody(i, args.blackbody) * green58_03[i-300] *  luminosity[i-300] * 683.002
	visible_gcd07 = 0
	for i in range(300, 701):
		visible_gcd07 += blackbody(i, args.blackbody) * green58_07[i-300] *  luminosity[i-300] * 683.002
	visible_gcd10 = 0
	for i in range(300, 701):
		visible_gcd10 += blackbody(i, args.blackbody) * green58_10[i-300] *  luminosity[i-300] * 683.002
	# gray
	visible_graycd0 = 0
	for i in range(300, 701):
		visible_graycd0 += blackbody(i, args.blackbody) * gray_0[i-300] *  luminosity[i-300] * 683.002
	visible_graycd03 = 0
	for i in range(300, 701):
		visible_graycd03 += blackbody(i, args.blackbody) * gray_03[i-300] *  luminosity[i-300] * 683.002
	visible_graycd07 = 0
	for i in range(300, 701):
		visible_graycd07 += blackbody(i, args.blackbody) * gray_07[i-300] *  luminosity[i-300] * 683.002
	visible_graycd10 = 0
	for i in range(300, 701):
		visible_graycd10 += blackbody(i, args.blackbody) * gray_10[i-300] *  luminosity[i-300] * 683.002
	
	# candela equivalents weighted by a custom luminosity function in place of CIE
	# to compare brightness for another species
	# This is probably not very useful because we're already quantifying species-specific
	# luminosity in more standard ways, but see also https://pmc.ncbi.nlm.nih.gov/articles/PMC11562817/
	
	# blue
	visible_bcde0 = 0
	for i in range(300, 701):
		visible_bcde0 += blackbody(i, args.blackbody) * blue47_0[i-300] *  sensitivity(i) * 683.002
	visible_bcde03 = 0
	for i in range(300, 701):
		visible_bcde03 += blackbody(i, args.blackbody) * blue47_03[i-300] * sensitivity(i) * 683.002
	visible_bcde07 = 0
	for i in range(300, 701):
		visible_bcde07 += blackbody(i, args.blackbody) * blue47_07[i-300] * sensitivity(i) * 683.002
	visible_bcde10 = 0
	for i in range(300, 701):
		visible_bcde10 += blackbody(i, args.blackbody) * blue47_10[i-300] * sensitivity(i) * 683.002
	# red
	visible_rcde0 = 0
	for i in range(300, 701):
		visible_rcde0 += blackbody(i, args.blackbody) * red25_0[i-300] * sensitivity(i) * 683.002
	visible_rcde03 = 0
	for i in range(300, 701):
		visible_rcde03 += blackbody(i, args.blackbody) * red25_03[i-300] * sensitivity(i) * 683.002
	visible_rcde07 = 0
	for i in range(300, 701):
		visible_rcde07 += blackbody(i, args.blackbody) * red25_07[i-300] * sensitivity(i) * 683.002
	visible_rcde10 = 0
	for i in range(300, 701):
		visible_rcde10 += blackbody(i, args.blackbody) * red25_10[i-300] * sensitivity(i) * 683.002
	# yellow
	visible_ycde0 = 0
	for i in range(300, 701):
		visible_ycde0 += blackbody(i, args.blackbody) * yellow15_0[i-300] * sensitivity(i) * 683.002
	visible_ycde03 = 0
	for i in range(300, 701):
		visible_ycde03 += blackbody(i, args.blackbody) * yellow15_03[i-300] * sensitivity(i) * 683.002
	visible_ycde07 = 0
	for i in range(300, 701):
		visible_ycde07 += blackbody(i, args.blackbody) * yellow15_07[i-300] * sensitivity(i) * 683.002
	visible_ycde10 = 0
	for i in range(300, 701):
		visible_ycde10 += blackbody(i, args.blackbody) * yellow15_10[i-300] * sensitivity(i) * 683.002
	# green
	visible_gcde0 = 0
	for i in range(300, 701):
		visible_gcde0 += blackbody(i, args.blackbody) * green58_0[i-300] * sensitivity(i) * 683.002
	visible_gcde03 = 0
	for i in range(300, 701):
		visible_gcde03 += blackbody(i, args.blackbody) * green58_03[i-300] * sensitivity(i) * 683.002
	visible_gcde07 = 0
	for i in range(300, 701):
		visible_gcde07 += blackbody(i, args.blackbody) * green58_07[i-300] * sensitivity(i) * 683.002
	visible_gcde10 = 0
	for i in range(300, 701):
		visible_gcde10 += blackbody(i, args.blackbody) * green58_10[i-300] * sensitivity(i) * 683.002
	# gray
	visible_graycde0 = 0
	for i in range(300, 701):
		visible_graycde0 += blackbody(i, args.blackbody) * gray_0[i-300] * sensitivity(i) * 683.002
	visible_graycde03 = 0
	for i in range(300, 701):
		visible_graycde03 += blackbody(i, args.blackbody) * gray_03[i-300] * sensitivity(i) * 683.002
	visible_graycde07 = 0
	for i in range(300, 701):
		visible_graycde07 += blackbody(i, args.blackbody) * gray_07[i-300] * sensitivity(i) * 683.002
	visible_graycde10 = 0
	for i in range(300, 701):
		visible_graycde10 += blackbody(i, args.blackbody) * gray_10[i-300] * sensitivity(i) * 683.002
	
	print("Radiance from 300-700 nm (W/m^2/sr): " + str(visible))
	print("Irradiance from 300-700 nm (W/m^2): " + str(visible*sr))
	print("Irradiance from 300-700 nm (uW/cm^2): " + str(visible*sr*1000000/100**2))
	
	# scaled to specified watt number
	print("")
	print("Scaled radiance/irradiance:")
	print("Radiance from 300-700 nm (W/m^2/sr): " + str(scale*visible))
	radiance = 0
	for i in range(501):
		radiance += ia[i]
	print("Radiance from 300-700 nm (photons/sec/sr): " + str(radiance))
	print("Irradiance (W/m^2): " + str(scale*visible*sr))
	print("Irradiance (uW/cm^2): " + str(scale*visible*sr*1000000/100**2))
	# just multiply this by 100 to get uW/cm^2
	print("Irradiance filtered through red 25 (uW/cm^2):")
	print("0: " + str(scale*visible_r0*sr*100))
	print("0.3: " + str(scale*visible_r03*sr*100))
	print("0.7: " + str(scale*visible_r07*sr*100))
	print("1.0: " + str(scale*visible_r10*sr*100))
	print("Irradiance filtered through yellow 15 (uW/cm^2):")
	print("0: " + str(scale*visible_y0*sr*100))
	print("0.3: " + str(scale*visible_y03*sr*100))
	print("0.7: " + str(scale*visible_y07*sr*100))
	print("1.0: " + str(scale*visible_y10*sr*100))
	print("Irradiance filtered through green 58 (uW/cm^2):")
	print("0: " + str(scale*visible_g0*sr*100))
	print("0.3: " + str(scale*visible_g03*sr*100))
	print("0.7: " + str(scale*visible_g07*sr*100))
	print("1.0: " + str(scale*visible_g10*sr*100))
	print("Irradiance filtered through blue 47 (uW/cm^2):")
	print("0: " + str(scale*visible_b0*sr*100))
	print("0.3: " + str(scale*visible_b03*sr*100))
	print("0.7: " + str(scale*visible_b07*sr*100))
	print("1.0: " + str(scale*visible_b10*sr*100))
	print("Irradiance filtered through gray (uW/cm^2):")
	print("0: " + str(scale*visible_gray0*sr*100))
	print("0.3: " + str(scale*visible_gray03*sr*100))
	print("0.7: " + str(scale*visible_gray07*sr*100))
	print("1.0: " + str(scale*visible_gray10*sr*100))
	print("Luminance from 300-700 nm for blue 47 (cd/m^2):")
	print("0: " + str(scale*visible_bcd0))
	print("0.3: " + str(scale*visible_bcd03))
	print("0.7: " + str(scale*visible_bcd07))
	print("1.0: " + str(scale*visible_bcd10))
	print("Luminance from 300-700 nm for red 25 (cd/m^2):")
	print("0: " + str(scale*visible_rcd0))
	print("0.3: " + str(scale*visible_rcd03))
	print("0.7: " + str(scale*visible_rcd07))
	print("1.0: " + str(scale*visible_rcd10))
	print("Luminance from 300-700 nm for yellow 15 (cd/m^2):")
	print("0: " + str(scale*visible_ycd0))
	print("0.3: " + str(scale*visible_ycd03))
	print("0.7: " + str(scale*visible_ycd07))
	print("1.0: " + str(scale*visible_ycd10))
	print("Luminance from 300-700 nm for green 58 (cd/m^2):")
	print("0: " + str(scale*visible_gcd0))
	print("0.3: " + str(scale*visible_gcd03))
	print("0.7: " + str(scale*visible_gcd07))
	print("1.0: " + str(scale*visible_gcd10))
	print("Luminance from 300-700 nm for gray (cd/m^2):")
	print("0: " + str(scale*visible_graycd0))
	print("0.3: " + str(scale*visible_graycd03))
	print("0.7: " + str(scale*visible_graycd07))
	print("1.0: " + str(scale*visible_graycd10))
	print("Illuminance for blue 47 (lx):")
	print("0: " + str(scale*visible_lx0))
	print("0.3: " + str(scale*visible_lx03))
	print("0.7: " + str(scale*visible_lx07))
	print("1.0: " + str(scale*visible_lx10))
	# "candela equivalent" brightness
	print("Luminance in candela equivalents (cde/m^2)")
	print("Blue 47:")
	print("0: " + str(scale*visible_bcde0))
	print("0.3: " + str(scale*visible_bcde03))
	print("0.7: " + str(scale*visible_bcde07))
	print("1.0: " + str(scale*visible_bcde10))
	print("Red 25:")
	print("0: " + str(scale*visible_rcde0))
	print("0.3: " + str(scale*visible_rcde03))
	print("0.7: " + str(scale*visible_rcde07))
	print("1.0: " + str(scale*visible_rcde10))
	print("Yellow 15:")
	print("0: " + str(scale*visible_ycde0))
	print("0.3: " + str(scale*visible_ycde03))
	print("0.7: " + str(scale*visible_ycde07))
	print("1.0: " + str(scale*visible_ycde10))
	print("Green 58:")
	print("0: " + str(scale*visible_gcde0))
	print("0.3: " + str(scale*visible_gcde03))
	print("0.7: " + str(scale*visible_gcde07))
	print("1.0: " + str(scale*visible_gcde10))
	print("Gray:")
	print("0: " + str(scale*visible_graycde0))
	print("0.3: " + str(scale*visible_graycde03))
	print("0.7: " + str(scale*visible_graycde07))
	print("1.0: " + str(scale*visible_graycde10))
	
	# plot curve
	xvalues = np.empty(61)
	yvalues = np.empty(61)
	for i in range(61):
		w = i*10 + 300
		xvalues[i] = w
		yvalues[i] = blackbody(w, args.blackbody)
	plt.plot(xvalues, yvalues)
	plt.show()
	print("")
	
	# plot brightness
	plt.subplot(1, 2, 1)
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_rcd0, scale*visible_rcd03, scale*visible_rcd07, scale*visible_rcd10], 'sr', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_ycd0, scale*visible_ycd03, scale*visible_ycd07, scale*visible_ycd10], 'Dy', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_gcd0, scale*visible_gcd03, scale*visible_gcd07, scale*visible_gcd10], '^g', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_bcd0, scale*visible_bcd03, scale*visible_bcd07, scale*visible_bcd10], 'ob', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_graycd0, scale*visible_graycd03, scale*visible_graycd07, scale*visible_graycd10], marker='v', linestyle='', color='gray', mec='k')
	plt.yscale('log')
	plt.subplot(1, 2, 2)
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_rcde0, scale*visible_rcde03, scale*visible_rcde07, scale*visible_rcde10], 'sr', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_ycde0, scale*visible_ycde03, scale*visible_ycde07, scale*visible_ycde10], 'Dy', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_gcde0, scale*visible_gcde03, scale*visible_gcde07, scale*visible_gcde10], '^g', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_bcde0, scale*visible_bcde03, scale*visible_bcde07, scale*visible_bcde10], 'ob', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_graycde0, scale*visible_graycde03, scale*visible_graycde07, scale*visible_graycde10], marker='v', linestyle='', color='gray', mec='k')
	plt.yscale('log')
	plt.show()

""" Munsell color cards: 5B, 7.5B, 10YR, 5GY
Reflectance spectra from here: https://www.munsellcolourscienceforpainters.com/MunsellResources/SpectralReflectancesOf2007MunsellBookOfColorGlossy.txt
I can't find any Munsell reflectance spectra that go below 380 nm, so we'll just have to
make do with this. Yet more limitations of a system designed for human vision when applied to
anything else. I'm guessing none of these things reflect much UV, though.
Note Python variables can't begin with a digit, so the arrays begin with "m" for Munsell.
Lightness/chroma values used by Gutierrez et al. 2011: 10YR5-8/10, 5B4-7/6, 7.5B5-8/4, 5GY5-8/10
Lighting specifications for absolute quantum catches:
* card dimensions: 12.5 x 7.5 cm
* illuminant: D65 light bulb 1.5 m above the test box, filtered with baking paper to provide
an illuminance of 115 lx (I hate trying to interpret these human-scaled units, why can't you
guys just use watts?)
** I can't find any measurements of transmission for baking paper, so we'll have to assume
it transmits all wavelengths equally. It probably doesn't.
** Converting lux to watts: 1 W/m^2 = 683.002 lx at 555 nm. First multiply the D65 spectrum by
the photopic luminosity function we determined earlier, then find the scaling factor to
lux by dividing 115 by the integral of D65 x luminosity. To find the scaling factor to watts,
take the value at 555 nm and multiply it by the lux scaling factor we just found, then divide
it by 683.002. Except we don't have a value for 555 nm, so we have to pick a different wavelength
(say, 550) and multiply it by the luminosity value.
* As before, we assume the animals are arbitrarily close to the stimuli.
* wavelength interval: 730 - 380 = 350
"""
# 11/02/2025 -- note that a paper I found (https://www.researchgate.net/publication/6764034_Transforming_reflectance_spectra_into_Munsell_color_space_by_using_prime_colors#pf3) uses a different set of reflectance spectra from the same site, the ones for the matte colors. Unfortunately 10YR5/10 and 5GY5-7/10 aren't included.
if (args.munsell):
	# lux to watts test
	#integral = 0
	#for i in range(d65.shape[0]):
	#	integral += d65[i] * (0.68990272*vpt(i*10 + 300, 565) + 0.34832189*vpt(i*10 + 300, 535))*human_filter(i*10 + 300)
	#lxscale = 115 / integral
	#wscale = d65[25] * (0.68990272*vpt(550, 565) + 0.34832189*vpt(550, 535))*human_filter(550) * lxscale / 683.002
	#integral1 = 0
	#for i in range(d65.shape[0]):
	#	integral1 += d65[i]
	#watts = integral1 * wscale
	#print(watts)
	# 5B4/6
	m5b4 = np.array([
		# placeholders for compatibility with spectral_rendering() and color_contrast()
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.178745, # 380
		0.178802,
		0.178866,
		0.178937,
		0.179031,
		0.192057,
		0.204833,
		0.218286,
		0.231725,
		0.243877,
		0.249609,
		0.24821,
		0.239283,
		0.224377,
		0.20376,
		0.177949,
		0.149619,
		0.120234,
		0.093083,
		0.07434,
		0.06398,
		0.05786,
		0.052017,
		0.047748,
		0.047064,
		0.049879,
		0.054122,
		0.057921,
		0.060903,
		0.06179,
		0.059133,
		0.056417,
		0.056283
		#0.058632,
		#0.062006,
		#0.066158, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 5B5/6
	m5b5 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.25803, # 380
		0.258079,
		0.258134,
		0.258196,
		0.258278,
		0.27231,
		0.289218,
		0.306161,
		0.323487,
		0.343154,
		0.361538,
		0.371382,
		0.368246,
		0.350666,
		0.320172,
		0.280711,
		0.238678,
		0.196393,
		0.15764,
		0.131575,
		0.118184,
		0.11143,
		0.105163,
		0.100138,
		0.101005,
		0.10889,
		0.121262,
		0.134619,
		0.144545,
		0.1436,
		0.131814,
		0.120245,
		0.115894
		#0.119029,
		#0.127256,
		#0.144525, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 5B6/6
	m5b6 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.374389, # 380
		0.374422,
		0.374459,
		0.374501,
		0.374574,
		0.400569,
		0.432946,
		0.461152,
		0.482049,
		0.497958,
		0.507464,
		0.511782,
		0.508655,
		0.494051,
		0.465085,
		0.422165,
		0.371097,
		0.315289,
		0.260843,
		0.22209,
		0.200565,
		0.188663,
		0.177296,
		0.168512,
		0.16899,
		0.179156,
		0.19418,
		0.207963,
		0.217457,
		0.217895,
		0.207884,
		0.197167,
		0.196282
		#0.203308,
		#0.214796,
		#0.231492, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 5B7/6
	m5b7 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.56984, # 380
		0.570055,
		0.570298,
		0.570571,
		0.570887,
		0.588852,
		0.608936,
		0.624073,
		0.637839,
		0.651186,
		0.659546,
		0.660242,
		0.65223,
		0.634396,
		0.607741,
		0.57188,
		0.528807,
		0.476569,
		0.41871,
		0.369485,
		0.331277,
		0.304091,
		0.280917,
		0.26285,
		0.252911,
		0.247474,
		0.24495,
		0.246894,
		0.253137,
		0.255556,
		0.251243,
		0.24333,
		0.23233
		#0.223291,
		#0.227917,
		#0.250112, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 7.5B5/4
	m75b5 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.261784, # 380
		0.261876,
		0.26198,
		0.262097,
		0.262235,
		0.272347,
		0.281291,
		0.289641,
		0.296922,
		0.301383,
		0.302441,
		0.300763,
		0.295468,
		0.287397,
		0.274845,
		0.257209,
		0.236216,
		0.212185,
		0.187201,
		0.169182,
		0.158874,
		0.152215,
		0.143971,
		0.136691,
		0.135623,
		0.140887,
		0.149415,
		0.158116,
		0.166126,
		0.167674,
		0.160939,
		0.152913,
		0.148741
		#0.149692,
		#0.155841,
		#0.172085, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 7.5B6/4
	m75b6 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.400169, # 380
		0.400346,
		0.400546,
		0.400771,
		0.401027,
		0.410578,
		0.420516,
		0.426882,
		0.428536,
		0.425265,
		0.420981,
		0.416922,
		0.412783,
		0.407012,
		0.397587,
		0.381616,
		0.358683,
		0.328331,
		0.292359,
		0.26341,
		0.245639,
		0.236064,
		0.226997,
		0.220085,
		0.220996,
		0.23089,
		0.245082,
		0.256564,
		0.261856,
		0.256971,
		0.242915,
		0.229333,
		0.225293
		#0.229032,
		#0.236621,
		#0.246813, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 7.5B7/4
	m75b7 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.518759, # 380
		0.518908,
		0.519076,
		0.519265,
		0.519485,
		0.534988,
		0.546915,
		0.558654,
		0.571024,
		0.57923,
		0.582938,
		0.582082,
		0.57567,
		0.56471,
		0.549078,
		0.526794,
		0.498858,
		0.463611,
		0.42245,
		0.389295,
		0.367166,
		0.352971,
		0.338885,
		0.327905,
		0.326381,
		0.333958,
		0.346241,
		0.357724,
		0.366898,
		0.368351,
		0.361228,
		0.352629,
		0.350633
		#0.353927,
		#0.362951,
		#0.379432, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 7.5B8/4
	m75b8 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.704227, # 380
		0.704429,
		0.704657,
		0.704914,
		0.70521,
		0.723661,
		0.738402,
		0.756363,
		0.778921,
		0.79822,
		0.808638,
		0.806992,
		0.792298,
		0.770147,
		0.741298,
		0.705929,
		0.667284,
		0.622554,
		0.571451,
		0.530456,
		0.502657,
		0.483397,
		0.462854,
		0.445685,
		0.439752,
		0.445489,
		0.457057,
		0.466592,
		0.474052,
		0.476135,
		0.469816,
		0.46184,
		0.462779
		#0.469476,
		#0.477307,
		#0.489983, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 10YR5/10
	m10yr5 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.012332, # 380
		0.012043,
		0.011717,
		0.011352,
		0.010949,
		0.016199,
		0.017594,
		0.018379,
		0.018991,
		0.019719,
		0.020958,
		0.023837,
		0.032334,
		0.060861,
		0.106745,
		0.140234,
		0.156951,
		0.174581,
		0.208961,
		0.259478,
		0.302108,
		0.325517,
		0.333406,
		0.334332,
		0.332448,
		0.329418,
		0.326735,
		0.323567,
		0.320636,
		0.317198,
		0.314651,
		0.312537,
		0.310873,
		#0.309138,
		#0.305939,
		#0.302531, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 10YR6/10
	m10yr6 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.034696, # 380
		0.034299,
		0.03385,
		0.033346,
		0.032786,
		0.036301,
		0.038371,
		0.040094,
		0.041566,
		0.042744,
		0.044244,
		0.047367,
		0.059028,
		0.09866,
		0.166275,
		0.220733,
		0.25004,
		0.276238,
		0.318942,
		0.378976,
		0.429468,
		0.458221,
		0.469551,
		0.47277,
		0.472918,
		0.471793,
		0.471481,
		0.469715,
		0.468048,
		0.465578,
		0.463179,
		0.461716,
		0.462856
		#0.461125,
		#0.458641,
		#0.455757, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 10YR7/10
	m10yr7 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.066231, # 380
		0.065675,
		0.065047,
		0.064344,
		0.063561,
		0.06886,
		0.072356,
		0.075154,
		0.077795,
		0.079896,
		0.081973,
		0.08584,
		0.101576,
		0.155603,
		0.25041,
		0.333024,
		0.381987,
		0.421049,
		0.476274,
		0.549704,
		0.608034,
		0.640799,
		0.654189,
		0.659383,
		0.662338,
		0.664762,
		0.668311,
		0.670207,
		0.670854,
		0.669071,
		0.668089,
		0.668937,
		0.671756
		#0.673282,
		#0.672417,
		#0.670845, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 10YR8/10
	m10yr8 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.12431, # 380
		0.123607,
		0.122814,
		0.121924,
		0.120932,
		0.12709,
		0.130249,
		0.133506,
		0.136786,
		0.139491,
		0.142047,
		0.146205,
		0.164531,
		0.233954,
		0.361093,
		0.477832,
		0.548543,
		0.597321,
		0.658102,
		0.739037,
		0.805007,
		0.843502,
		0.86041,
		0.867769,
		0.871954,
		0.87487,
		0.879254,
		0.881513,
		0.882222,
		0.879091,
		0.878238,
		0.879804,
		0.884059
		#0.885425,
		#0.882944,
		#0.880923, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 5GY5/10
	m5gy5 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.012664, # 380
		0.012341,
		0.011977,
		0.011568,
		0.011115,
		0.015092,
		0.016432,
		0.017577,
		0.018949,
		0.021024,
		0.023912,
		0.027906,
		0.039837,
		0.084208,
		0.178723,
		0.282086,
		0.339094,
		0.342394,
		0.314772,
		0.278631,
		0.23888,
		0.197525,
		0.156496,
		0.124397,
		0.105102,
		0.094838,
		0.088166,
		0.082435,
		0.079019,
		0.079228,
		0.082376,
		0.087889,
		0.09505
		#0.100989,
		#0.100258,
		#0.098844, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 5GY6/10
	m5gy6 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.032315, # 380
		0.031884,
		0.031398,
		0.030852,
		0.030242,
		0.032329,
		0.03415,
		0.035964,
		0.037998,
		0.040435,
		0.043275,
		0.04794,
		0.061653,
		0.116183,
		0.23943,
		0.388322,
		0.485203,
		0.498968,
		0.458571,
		0.404354,
		0.346137,
		0.286732,
		0.227724,
		0.180557,
		0.151948,
		0.13625,
		0.126095,
		0.117313,
		0.112108,
		0.112244,
		0.117827,
		0.126935,
		0.138017,
		#0.145016,
		#0.144087,
		#0.140965, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 5GY7/10
	m5gy7 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.040501, # 380
		0.039827,
		0.039066,
		0.038213,
		0.037259,
		0.040034,
		0.044292,
		0.051027,
		0.062755,
		0.083806,
		0.120887,
		0.187442,
		0.289999,
		0.403185,
		0.487152,
		0.527593,
		0.542556,
		0.549095,
		0.551625,
		0.545624,
		0.517071,
		0.469334,
		0.408376,
		0.352003,
		0.31403,
		0.291419,
		0.2754,
		0.260754,
		0.251191,
		0.250301,
		0.257165,
		0.269653,
		0.286497
		#0.29664,
		#0.295497,
		#0.293521, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# 5GY8/10
	m5gy8 = np.array([
		np.nan, # 300
		np.nan, # 310
		np.nan, # 320
		np.nan, # 330
		np.nan, # 340
		np.nan, # 350
		np.nan, # 360
		np.nan, # 370
		0.065944, # 380
		0.065028,
		0.063995,
		0.062836,
		0.061541,
		0.067324,
		0.073985,
		0.084234,
		0.103077,
		0.138398,
		0.195627,
		0.279701,
		0.398815,
		0.543652,
		0.664439,
		0.725444,
		0.746586,
		0.745688,
		0.727935,
		0.704973,
		0.671825,
		0.628966,
		0.573402,
		0.519301,
		0.48158,
		0.459428,
		0.444344,
		0.429995,
		0.420806,
		0.421169,
		0.431991,
		0.451351,
		0.475862
		#0.491195,
		#0.491365,
		#0.487952, # 730
		#np.nan, # 740
		#np.nan, # 750
		#np.nan, # 760
		#np.nan, # 770
		#np.nan, # 780
		#np.nan, # 790
		#np.nan # 800
	])
	
	# linear interpolation
#	x10nm = np.empty(44)
#	for i in range(44):
#		x10nm[i] = i*10 + 300
#	
#	x1nm = np.empty(431)
#	for i in range(431):
#		x1nm[i] = i + 300
#	
#	# 5B
#	m5b4_1nm = np.interp(x1nm, x10nm, m5b4)
#	m5b5_1nm = np.interp(x1nm, x10nm, m5b5)
#	m5b6_1nm = np.interp(x1nm, x10nm, m5b6)
#	m5b7_1nm = np.interp(x1nm, x10nm, m5b7)
#	# 7.5B
#	m75b5_1nm = np.interp(x1nm, x10nm, m75b5)
#	m75b6_1nm = np.interp(x1nm, x10nm, m75b6)
#	m75b7_1nm = np.interp(x1nm, x10nm, m75b7)
#	m75b8_1nm = np.interp(x1nm, x10nm, m75b8)
#	# 10YR
#	m10yr5_1nm = np.interp(x1nm, x10nm, m10yr5)
#	m10yr6_1nm = np.interp(x1nm, x10nm, m10yr6)
#	m10yr7_1nm = np.interp(x1nm, x10nm, m10yr7)
#	m10yr8_1nm = np.interp(x1nm, x10nm, m10yr8)
#	# 5GY
#	m5gy5_1nm = np.interp(x1nm, x10nm, m5gy5)
#	m5gy6_1nm = np.interp(x1nm, x10nm, m5gy6)
#	m5gy7_1nm = np.interp(x1nm, x10nm, m5gy7)
#	m5gy8_1nm = np.interp(x1nm, x10nm, m5gy8)

	if (args.mv2 or args.mv3):
		# partial fix: use the correct matte spectra for 5B, 7.5B, 10YR6-8 and 5GY8
		# 5B4-7/6
		m5b4_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   9.0700000e-02,
		   8.7700000e-02,
		   8.6800000e-02,
		   8.8000000e-02,
		   8.9000000e-02,
		   9.3100000e-02,
		   9.4800000e-02,
		   9.7100000e-02,
		   1.0050000e-01,
		   1.0430000e-01,
		   1.0670000e-01,
		   1.1010000e-01,
		   1.1290000e-01,
		   1.1690000e-01,
		   1.1640000e-01,
		   1.1860000e-01,
		   1.2330000e-01,
		   1.2560000e-01,
		   1.2600000e-01,
		   1.2910000e-01,
		   1.3090000e-01,
		   1.3140000e-01,
		   1.3270000e-01,
		   1.3560000e-01,
		   1.3620000e-01,
		   1.3830000e-01,
		   1.3940000e-01,
		   1.4170000e-01,
		   1.4160000e-01,
		   1.4260000e-01,
		   1.4530000e-01,
		   1.4490000e-01,
		   1.4310000e-01,
		   1.4440000e-01,
		   1.4570000e-01,
		   1.4620000e-01,
		   1.4530000e-01,
		   1.4670000e-01,
		   1.4610000e-01,
		   1.4800000e-01,
		   1.4820000e-01,
		   1.4940000e-01,
		   1.4960000e-01,
		   1.5010000e-01,
		   1.5200000e-01,
		   1.5410000e-01,
		   1.5270000e-01,
		   1.5490000e-01,
		   1.5620000e-01,
		   1.5560000e-01,
		   1.5720000e-01,
		   1.6000000e-01,
		   1.6070000e-01,
		   1.6120000e-01,
		   1.6400000e-01,
		   1.6690000e-01,
		   1.6670000e-01,
		   1.6770000e-01,
		   1.7080000e-01,
		   1.7290000e-01,
		   1.7290000e-01,
		   1.7390000e-01,
		   1.7630000e-01,
		   1.7710000e-01,
		   1.7720000e-01,
		   1.7780000e-01,
		   1.8080000e-01,
		   1.8160000e-01,
		   1.8290000e-01,
		   1.8720000e-01,
		   1.8730000e-01,
		   1.8690000e-01,
		   1.8990000e-01,
		   1.9180000e-01,
		   1.9210000e-01,
		   1.9370000e-01,
		   1.9720000e-01,
		   1.9900000e-01,
		   1.9970000e-01,
		   2.0310000e-01,
		   2.0650000e-01,
		   2.0790000e-01,
		   2.0880000e-01,
		   2.1250000e-01,
		   2.1460000e-01,
		   2.1440000e-01,
		   2.1700000e-01,
		   2.1890000e-01,
		   2.1890000e-01,
		   2.2140000e-01,
		   2.2400000e-01,
		   2.2350000e-01,
		   2.2350000e-01,
		   2.2540000e-01,
		   2.2720000e-01,
		   2.2740000e-01,
		   2.2720000e-01,
		   2.2860000e-01,
		   2.2910000e-01,
		   2.2960000e-01,
		   2.2960000e-01,
		   2.3070000e-01,
		   2.2900000e-01,
		   2.2830000e-01,
		   2.2940000e-01,
		   2.2910000e-01,
		   2.2710000e-01,
		   2.2790000e-01,
		   2.2850000e-01,
		   2.2580000e-01,
		   2.2470000e-01,
		   2.2490000e-01,
		   2.2350000e-01,
		   2.2270000e-01,
		   2.2250000e-01,
		   2.2100000e-01,
		   2.1870000e-01,
		   2.1740000e-01,
		   2.1870000e-01,
		   2.1770000e-01,
		   2.1540000e-01,
		   2.1190000e-01,
		   2.1230000e-01,
		   2.1000000e-01,
		   2.0850000e-01,
		   2.0710000e-01,
		   2.0560000e-01,
		   2.0200000e-01,
		   2.0260000e-01,
		   2.0020000e-01,
		   1.9710000e-01,
		   1.9390000e-01,
		   1.9310000e-01,
		   1.9040000e-01,
		   1.8740000e-01,
		   1.8430000e-01,
		   1.8440000e-01,
		   1.8070000e-01,
		   1.7520000e-01,
		   1.7420000e-01,
		   1.7450000e-01,
		   1.6990000e-01,
		   1.6690000e-01,
		   1.6490000e-01,
		   1.6170000e-01,
		   1.5880000e-01,
		   1.5710000e-01,
		   1.5400000e-01,
		   1.5060000e-01,
		   1.4800000e-01,
		   1.4640000e-01,
		   1.4330000e-01,
		   1.4040000e-01,
		   1.3810000e-01,
		   1.3650000e-01,
		   1.3230000e-01,
		   1.2950000e-01,
		   1.2710000e-01,
		   1.2370000e-01,
		   1.2130000e-01,
		   1.2080000e-01,
		   1.1830000e-01,
		   1.1460000e-01,
		   1.1170000e-01,
		   1.1270000e-01,
		   1.1000000e-01,
		   1.0650000e-01,
		   1.0360000e-01,
		   1.0330000e-01,
		   1.0110000e-01,
		   9.8000000e-02,
		   9.5000000e-02,
		   9.4200000e-02,
		   9.2400000e-02,
		   9.1400000e-02,
		   8.8800000e-02,
		   8.5600000e-02,
		   8.2100000e-02,
		   8.1800000e-02,
		   7.7300000e-02,
		   7.2100000e-02,
		   7.3100000e-02,
		   7.7400000e-02,
		   7.6600000e-02,
		   7.6500000e-02,
		   7.2500000e-02,
		   7.2000000e-02,
		   7.2300000e-02,
		   7.0800000e-02,
		   6.9100000e-02,
		   6.9100000e-02,
		   6.8300000e-02,
		   6.7300000e-02,
		   6.5400000e-02,
		   6.6400000e-02,
		   6.5700000e-02,
		   6.4500000e-02,
		   6.3700000e-02,
		   6.4900000e-02,
		   6.4200000e-02,
		   6.2900000e-02,
		   6.3200000e-02,
		   6.2500000e-02,
		   6.0900000e-02,
		   6.1800000e-02,
		   6.2200000e-02,
		   6.0600000e-02,
		   6.0200000e-02,
		   6.1200000e-02,
		   6.1200000e-02,
		   6.0300000e-02,
		   6.1400000e-02,
		   6.1100000e-02,
		   5.9300000e-02,
		   5.9100000e-02,
		   5.9100000e-02,
		   5.8900000e-02,
		   5.8300000e-02,
		   5.9000000e-02,
		   6.0200000e-02,
		   5.8100000e-02,
		   5.6700000e-02,
		   5.8400000e-02,
		   5.9300000e-02,
		   5.7100000e-02,
		   5.6900000e-02,
		   5.8500000e-02,
		   5.8300000e-02,
		   5.7300000e-02,
		   5.7500000e-02,
		   5.7700000e-02,
		   5.7900000e-02,
		   5.7200000e-02,
		   5.8000000e-02,
		   5.7700000e-02,
		   5.7500000e-02,
		   5.8400000e-02,
		   5.8500000e-02,
		   5.6900000e-02,
		   5.7700000e-02,
		   5.9400000e-02,
		   5.6000000e-02,
		   5.4900000e-02,
		   5.6900000e-02,
		   5.8100000e-02,
		   5.6200000e-02,
		   5.7100000e-02,
		   5.7900000e-02,
		   5.8200000e-02,
		   5.6200000e-02,
		   5.7700000e-02,
		   5.7400000e-02,
		   5.6600000e-02,
		   5.6600000e-02,
		   5.8700000e-02,
		   5.7700000e-02,
		   5.6400000e-02,
		   5.6400000e-02,
		   5.8100000e-02,
		   5.8100000e-02,
		   5.7300000e-02,
		   5.7700000e-02,
		   5.7400000e-02,
		   5.8100000e-02,
		   5.9400000e-02,
		   5.9400000e-02,
		   5.7900000e-02,
		   5.8700000e-02,
		   5.9400000e-02,
		   5.9700000e-02,
		   5.8800000e-02,
		   5.8100000e-02,
		   5.9400000e-02,
		   5.8200000e-02,
		   5.7400000e-02,
		   6.0100000e-02,
		   6.0500000e-02,
		   5.9700000e-02,
		   6.1200000e-02,
		   6.1800000e-02,
		   6.0700000e-02,
		   6.0500000e-02,
		   6.1800000e-02,
		   6.2000000e-02,
		   6.1000000e-02,
		   5.9600000e-02,
		   5.9800000e-02,
		   5.8700000e-02,
		   5.8200000e-02,
		   6.0100000e-02,
		   6.1200000e-02,
		   6.0900000e-02,
		   6.2100000e-02,
		   6.3100000e-02,
		   6.1900000e-02,
		   6.1400000e-02,
		   6.2100000e-02,
		   6.1500000e-02,
		   6.1100000e-02,
		   6.2000000e-02,
		   6.3000000e-02,
		   6.3400000e-02,
		   6.2500000e-02,
		   6.2100000e-02,
		   6.2400000e-02,
		   6.0900000e-02,
		   6.0100000e-02,
		   5.8800000e-02,
		   5.6400000e-02,
		   5.7300000e-02,
		   6.4000000e-02,
		   6.3000000e-02,
		   6.2900000e-02,
		   6.3600000e-02,
		   6.1900000e-02,
		   6.2000000e-02,
		   6.3400000e-02,
		   6.2100000e-02,
		   6.1000000e-02,
		   6.2200000e-02,
		   6.2900000e-02
	#	   6.2100000e-02, 701
	#	   6.2600000e-02,
	#	   6.2600000e-02,
	#	   6.3300000e-02,
	#	   6.3100000e-02,
	#	   6.3300000e-02,
	#	   6.4000000e-02,
	#	   6.4500000e-02,
	#	   6.3200000e-02,
	#	   6.5800000e-02,
	#	   6.6600000e-02,
	#	   6.4800000e-02,
	#	   6.4700000e-02,
	#	   6.6600000e-02,
	#	   6.6900000e-02,
	#	   6.6800000e-02,
	#	   6.6000000e-02,
	#	   6.6500000e-02,
	#	   6.7800000e-02,
	#	   6.8000000e-02,
	#	   6.8500000e-02,
	#	   6.8900000e-02,
	#	   6.9100000e-02,
	#	   7.0500000e-02,
	#	   7.2300000e-02,
	#	   7.0800000e-02,
	#	   7.2200000e-02,
	#	   7.3700000e-02,
	#	   7.5000000e-02,
	#	   7.4500000e-02,
	#	   7.5500000e-02,
	#	   7.5700000e-02,
	#	   7.6100000e-02,
	#	   7.6200000e-02,
	#	   7.8000000e-02,
	#	   7.8400000e-02,
	#	   7.9700000e-02,
	#	   7.9900000e-02,
	#	   8.0800000e-02,
	#	   8.1100000e-02,
	#	   8.2000000e-02,
	#	   8.3500000e-02,
	#	   8.5300000e-02,
	#	   8.5300000e-02,
	#	   8.6900000e-02,
	#	   8.8500000e-02,
	#	   9.0600000e-02,
	#	   9.2200000e-02,
	#	   9.4700000e-02,
	#	   9.6100000e-02,
	#	   9.7200000e-02,
	#	   1.0050000e-01,
	#	   1.0150000e-01,
	#	   1.0240000e-01,
	#	   1.0380000e-01,
	#	   1.0650000e-01,
	#	   1.0870000e-01,
	#	   1.0930000e-01,
	#	   1.1110000e-01,
	#	   1.1290000e-01,
	#	   1.1280000e-01,
	#	   1.1290000e-01,
	#	   1.1580000e-01,
	#	   1.1680000e-01,
	#	   1.1660000e-01,
	#	   1.1840000e-01,
	#	   1.2150000e-01,
	#	   1.2210000e-01,
	#	   1.2320000e-01,
	#	   1.2330000e-01,
	#	   1.2360000e-01,
	#	   1.2230000e-01,
	#	   1.2260000e-01,
	#	   1.2400000e-01,
	#	   1.2350000e-01,
	#	   1.2320000e-01,
	#	   1.2510000e-01,
	#	   1.2570000e-01,
	#	   1.2600000e-01,
	#	   1.2710000e-01,
	#	   1.2640000e-01,
	#	   1.2360000e-01,
	#	   1.2700000e-01,
	#	   1.2900000e-01,
	#	   1.2920000e-01,
	#	   1.2720000e-01,
	#	   1.2650000e-01,
	#	   1.2770000e-01,
	#	   1.2920000e-01,
	#	   1.2890000e-01,
	#	   1.3080000e-01,
	#	   1.3250000e-01,
	#	   1.3220000e-01,
	#	   1.3320000e-01,
	#	   1.3620000e-01,
	#	   1.3710000e-01,
	#	   1.3690000e-01,
	#	   1.3630000e-01,
	#	   1.3890000e-01,
	#	   1.3960000e-01
		])

		m5b5_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   1.0170000e-01,
		   1.0280000e-01,
		   1.0690000e-01,
		   1.1290000e-01,
		   1.1250000e-01,
		   1.1730000e-01,
		   1.2380000e-01,
		   1.3000000e-01,
		   1.3250000e-01,
		   1.4030000e-01,
		   1.4550000e-01,
		   1.5330000e-01,
		   1.5820000e-01,
		   1.6690000e-01,
		   1.7220000e-01,
		   1.7960000e-01,
		   1.8430000e-01,
		   1.9320000e-01,
		   1.9610000e-01,
		   2.0000000e-01,
		   2.0620000e-01,
		   2.1160000e-01,
		   2.1100000e-01,
		   2.1580000e-01,
		   2.2020000e-01,
		   2.2340000e-01,
		   2.2430000e-01,
		   2.2920000e-01,
		   2.2980000e-01,
		   2.3120000e-01,
		   2.3340000e-01,
		   2.3530000e-01,
		   2.3520000e-01,
		   2.3480000e-01,
		   2.3620000e-01,
		   2.3720000e-01,
		   2.3480000e-01,
		   2.3610000e-01,
		   2.3870000e-01,
		   2.4230000e-01,
		   2.4230000e-01,
		   2.4390000e-01,
		   2.4470000e-01,
		   2.4550000e-01,
		   2.4620000e-01,
		   2.4580000e-01,
		   2.4510000e-01,
		   2.4720000e-01,
		   2.4840000e-01,
		   2.5040000e-01,
		   2.5090000e-01,
		   2.5290000e-01,
		   2.5560000e-01,
		   2.5630000e-01,
		   2.5610000e-01,
		   2.5950000e-01,
		   2.6210000e-01,
		   2.6100000e-01,
		   2.6050000e-01,
		   2.6480000e-01,
		   2.6590000e-01,
		   2.6520000e-01,
		   2.6680000e-01,
		   2.6950000e-01,
		   2.7020000e-01,
		   2.7060000e-01,
		   2.7250000e-01,
		   2.7410000e-01,
		   2.7430000e-01,
		   2.7390000e-01,
		   2.7520000e-01,
		   2.7770000e-01,
		   2.7700000e-01,
		   2.7920000e-01,
		   2.7880000e-01,
		   2.7850000e-01,
		   2.8130000e-01,
		   2.8340000e-01,
		   2.8200000e-01,
		   2.8450000e-01,
		   2.8780000e-01,
		   2.8960000e-01,
		   2.8920000e-01,
		   2.9060000e-01,
		   2.9310000e-01,
		   2.9410000e-01,
		   2.9210000e-01,
		   2.9470000e-01,
		   2.9540000e-01,
		   2.9440000e-01,
		   2.9590000e-01,
		   2.9730000e-01,
		   2.9650000e-01,
		   2.9640000e-01,
		   2.9790000e-01,
		   2.9860000e-01,
		   2.9780000e-01,
		   2.9860000e-01,
		   3.0010000e-01,
		   3.0030000e-01,
		   2.9910000e-01,
		   2.9890000e-01,
		   2.9970000e-01,
		   2.9720000e-01,
		   2.9690000e-01,
		   2.9790000e-01,
		   2.9810000e-01,
		   2.9670000e-01,
		   2.9760000e-01,
		   2.9750000e-01,
		   2.9560000e-01,
		   2.9260000e-01,
		   2.9470000e-01,
		   2.9400000e-01,
		   2.9290000e-01,
		   2.9360000e-01,
		   2.9410000e-01,
		   2.9290000e-01,
		   2.9170000e-01,
		   2.9100000e-01,
		   2.9160000e-01,
		   2.9040000e-01,
		   2.9000000e-01,
		   2.8950000e-01,
		   2.8720000e-01,
		   2.8470000e-01,
		   2.8560000e-01,
		   2.8400000e-01,
		   2.8100000e-01,
		   2.8140000e-01,
		   2.8030000e-01,
		   2.7680000e-01,
		   2.7550000e-01,
		   2.7460000e-01,
		   2.7370000e-01,
		   2.7170000e-01,
		   2.6970000e-01,
		   2.6780000e-01,
		   2.6560000e-01,
		   2.6200000e-01,
		   2.6190000e-01,
		   2.5890000e-01,
		   2.5460000e-01,
		   2.5230000e-01,
		   2.5240000e-01,
		   2.4840000e-01,
		   2.4550000e-01,
		   2.4160000e-01,
		   2.3910000e-01,
		   2.3750000e-01,
		   2.3550000e-01,
		   2.3190000e-01,
		   2.2880000e-01,
		   2.2470000e-01,
		   2.2370000e-01,
		   2.1980000e-01,
		   2.1640000e-01,
		   2.1160000e-01,
		   2.0970000e-01,
		   2.0300000e-01,
		   1.9830000e-01,
		   1.9500000e-01,
		   1.9290000e-01,
		   1.8730000e-01,
		   1.8570000e-01,
		   1.8280000e-01,
		   1.8030000e-01,
		   1.7420000e-01,
		   1.7390000e-01,
		   1.6950000e-01,
		   1.6570000e-01,
		   1.6160000e-01,
		   1.6060000e-01,
		   1.5750000e-01,
		   1.5420000e-01,
		   1.5020000e-01,
		   1.4810000e-01,
		   1.4420000e-01,
		   1.4110000e-01,
		   1.3600000e-01,
		   1.3170000e-01,
		   1.2830000e-01,
		   1.2990000e-01,
		   1.2500000e-01,
		   1.2490000e-01,
		   1.2310000e-01,
		   1.2000000e-01,
		   1.1750000e-01,
		   1.1670000e-01,
		   1.1410000e-01,
		   1.1160000e-01,
		   1.1040000e-01,
		   1.0740000e-01,
		   1.0540000e-01,
		   1.0630000e-01,
		   1.0580000e-01,
		   1.0320000e-01,
		   1.0250000e-01,
		   1.0370000e-01,
		   1.0240000e-01,
		   9.9300000e-02,
		   1.0000000e-01,
		   1.0020000e-01,
		   9.7600000e-02,
		   9.5700000e-02,
		   9.5100000e-02,
		   9.5500000e-02,
		   9.3500000e-02,
		   9.2900000e-02,
		   9.4500000e-02,
		   9.3300000e-02,
		   9.1800000e-02,
		   9.1500000e-02,
		   9.0100000e-02,
		   8.9300000e-02,
		   8.9600000e-02,
		   9.0400000e-02,
		   8.7600000e-02,
		   8.6300000e-02,
		   8.6500000e-02,
		   8.7100000e-02,
		   8.5000000e-02,
		   8.3600000e-02,
		   8.5200000e-02,
		   8.4400000e-02,
		   8.2200000e-02,
		   8.3500000e-02,
		   8.4100000e-02,
		   8.3000000e-02,
		   8.2800000e-02,
		   8.4300000e-02,
		   8.3300000e-02,
		   8.2700000e-02,
		   8.3300000e-02,
		   8.4400000e-02,
		   8.3500000e-02,
		   8.2400000e-02,
		   8.3100000e-02,
		   8.3300000e-02,
		   8.1800000e-02,
		   8.3000000e-02,
		   8.2300000e-02,
		   8.0700000e-02,
		   8.1300000e-02,
		   8.3300000e-02,
		   8.2800000e-02,
		   8.1000000e-02,
		   8.2500000e-02,
		   8.3300000e-02,
		   8.2600000e-02,
		   8.1800000e-02,
		   8.2400000e-02,
		   8.3100000e-02,
		   8.2700000e-02,
		   8.2200000e-02,
		   8.2100000e-02,
		   8.1100000e-02,
		   8.1300000e-02,
		   8.2400000e-02,
		   8.3100000e-02,
		   8.1500000e-02,
		   8.1200000e-02,
		   8.3200000e-02,
		   8.4000000e-02,
		   8.2200000e-02,
		   8.2500000e-02,
		   8.2800000e-02,
		   8.3200000e-02,
		   8.4300000e-02,
		   8.3900000e-02,
		   8.3500000e-02,
		   8.4200000e-02,
		   8.4700000e-02,
		   8.3300000e-02,
		   8.3000000e-02,
		   8.3500000e-02,
		   8.5200000e-02,
		   8.4100000e-02,
		   8.3500000e-02,
		   8.5700000e-02,
		   8.5700000e-02,
		   8.4500000e-02,
		   8.5500000e-02,
		   8.5800000e-02,
		   8.5400000e-02,
		   8.3400000e-02,
		   8.4100000e-02,
		   8.2900000e-02,
		   8.2300000e-02,
		   8.2800000e-02,
		   8.5100000e-02,
		   8.6500000e-02,
		   8.5200000e-02,
		   8.5600000e-02,
		   8.6900000e-02,
		   8.5500000e-02,
		   8.4500000e-02,
		   8.5500000e-02,
		   8.6100000e-02,
		   8.6700000e-02,
		   8.7500000e-02,
		   8.6500000e-02,
		   8.4700000e-02,
		   8.5700000e-02,
		   8.6800000e-02,
		   8.5600000e-02,
		   8.3000000e-02,
		   8.0200000e-02,
		   7.9700000e-02,
		   8.2400000e-02,
		   8.7700000e-02,
		   8.5500000e-02,
		   8.5200000e-02,
		   8.4600000e-02,
		   8.6100000e-02,
		   8.6100000e-02,
		   8.5500000e-02,
		   8.7200000e-02,
		   8.7100000e-02,
		   8.5700000e-02,
		   8.6300000e-02
	#	   8.6900000e-02, 701
	#	   8.5200000e-02,
	#	   8.5600000e-02,
	#	   8.6600000e-02,
	#	   8.6300000e-02,
	#	   8.5700000e-02,
	#	   8.5900000e-02,
	#	   8.6500000e-02,
	#	   8.7900000e-02,
	#	   8.6600000e-02,
	#	   8.7200000e-02,
	#	   8.6800000e-02,
	#	   8.7000000e-02,
	#	   8.7400000e-02,
	#	   8.7800000e-02,
	#	   8.8000000e-02,
	#	   8.9200000e-02,
	#	   9.0200000e-02,
	#	   9.0500000e-02,
	#	   8.9300000e-02,
	#	   9.1200000e-02,
	#	   9.2900000e-02,
	#	   9.3800000e-02,
	#	   9.3500000e-02,
	#	   9.5100000e-02,
	#	   9.6000000e-02,
	#	   9.7800000e-02,
	#	   9.9200000e-02,
	#	   1.0100000e-01,
	#	   1.0170000e-01,
	#	   1.0150000e-01,
	#	   1.0240000e-01,
	#	   1.0410000e-01,
	#	   1.0360000e-01,
	#	   1.0540000e-01,
	#	   1.0910000e-01,
	#	   1.1090000e-01,
	#	   1.1120000e-01,
	#	   1.1660000e-01,
	#	   1.1880000e-01,
	#	   1.2090000e-01,
	#	   1.2320000e-01,
	#	   1.2640000e-01,
	#	   1.2760000e-01,
	#	   1.3160000e-01,
	#	   1.3390000e-01,
	#	   1.3800000e-01,
	#	   1.4130000e-01,
	#	   1.4500000e-01,
	#	   1.4860000e-01,
	#	   1.5140000e-01,
	#	   1.5230000e-01,
	#	   1.5620000e-01,
	#	   1.5950000e-01,
	#	   1.6240000e-01,
	#	   1.6500000e-01,
	#	   1.6990000e-01,
	#	   1.7140000e-01,
	#	   1.7450000e-01,
	#	   1.7620000e-01,
	#	   1.7830000e-01,
	#	   1.7970000e-01,
	#	   1.8500000e-01,
	#	   1.8910000e-01,
	#	   1.9010000e-01,
	#	   1.8840000e-01,
	#	   1.9050000e-01,
	#	   1.9380000e-01,
	#	   1.9370000e-01,
	#	   1.9380000e-01,
	#	   1.9450000e-01,
	#	   1.9520000e-01,
	#	   1.9560000e-01,
	#	   1.9800000e-01,
	#	   1.9870000e-01,
	#	   1.9750000e-01,
	#	   1.9880000e-01,
	#	   2.0070000e-01,
	#	   1.9920000e-01,
	#	   2.0030000e-01,
	#	   2.0020000e-01,
	#	   2.0060000e-01,
	#	   2.0420000e-01,
	#	   2.0440000e-01,
	#	   2.0420000e-01,
	#	   2.0450000e-01,
	#	   2.0550000e-01,
	#	   2.0760000e-01,
	#	   2.0980000e-01,
	#	   2.0980000e-01,
	#	   2.1010000e-01,
	#	   2.1190000e-01,
	#	   2.1340000e-01,
	#	   2.1210000e-01,
	#	   2.1070000e-01,
	#	   2.1230000e-01,
	#	   2.1580000e-01,
	#	   2.1680000e-01,
	#	   2.1790000e-01,
	#	   2.1920000e-01
		])

		m5b6_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   1.2370000e-01,
		   1.2890000e-01,
		   1.3430000e-01,
		   1.4120000e-01,
		   1.4270000e-01,
		   1.4870000e-01,
		   1.5720000e-01,
		   1.6710000e-01,
		   1.7290000e-01,
		   1.8120000e-01,
		   1.8970000e-01,
		   2.0200000e-01,
		   2.0880000e-01,
		   2.2130000e-01,
		   2.3010000e-01,
		   2.4100000e-01,
		   2.4820000e-01,
		   2.6280000e-01,
		   2.6920000e-01,
		   2.7740000e-01,
		   2.8590000e-01,
		   2.9580000e-01,
		   2.9980000e-01,
		   3.0910000e-01,
		   3.1340000e-01,
		   3.1940000e-01,
		   3.2180000e-01,
		   3.2700000e-01,
		   3.3040000e-01,
		   3.3470000e-01,
		   3.3630000e-01,
		   3.3930000e-01,
		   3.3970000e-01,
		   3.3980000e-01,
		   3.4170000e-01,
		   3.4330000e-01,
		   3.4050000e-01,
		   3.4370000e-01,
		   3.4550000e-01,
		   3.4640000e-01,
		   3.4640000e-01,
		   3.4770000e-01,
		   3.4940000e-01,
		   3.5160000e-01,
		   3.5250000e-01,
		   3.5290000e-01,
		   3.5320000e-01,
		   3.5320000e-01,
		   3.5400000e-01,
		   3.5580000e-01,
		   3.5520000e-01,
		   3.5570000e-01,
		   3.5910000e-01,
		   3.6120000e-01,
		   3.6030000e-01,
		   3.6180000e-01,
		   3.6460000e-01,
		   3.6530000e-01,
		   3.6530000e-01,
		   3.6880000e-01,
		   3.6960000e-01,
		   3.7060000e-01,
		   3.7300000e-01,
		   3.7410000e-01,
		   3.7410000e-01,
		   3.7580000e-01,
		   3.7920000e-01,
		   3.8100000e-01,
		   3.7950000e-01,
		   3.8020000e-01,
		   3.8330000e-01,
		   3.8400000e-01,
		   3.8400000e-01,
		   3.8790000e-01,
		   3.9000000e-01,
		   3.8970000e-01,
		   3.9150000e-01,
		   3.9310000e-01,
		   3.9420000e-01,
		   3.9700000e-01,
		   3.9850000e-01,
		   4.0040000e-01,
		   4.0030000e-01,
		   4.0170000e-01,
		   4.0490000e-01,
		   4.0530000e-01,
		   4.0340000e-01,
		   4.0750000e-01,
		   4.0950000e-01,
		   4.0810000e-01,
		   4.0900000e-01,
		   4.0960000e-01,
		   4.0940000e-01,
		   4.1030000e-01,
		   4.1100000e-01,
		   4.1260000e-01,
		   4.1280000e-01,
		   4.1290000e-01,
		   4.1430000e-01,
		   4.1460000e-01,
		   4.1110000e-01,
		   4.0930000e-01,
		   4.1130000e-01,
		   4.1150000e-01,
		   4.1160000e-01,
		   4.1210000e-01,
		   4.1130000e-01,
		   4.0850000e-01,
		   4.0930000e-01,
		   4.0830000e-01,
		   4.0810000e-01,
		   4.0730000e-01,
		   4.0750000e-01,
		   4.0600000e-01,
		   4.0400000e-01,
		   4.0390000e-01,
		   4.0420000e-01,
		   4.0200000e-01,
		   4.0150000e-01,
		   4.0260000e-01,
		   4.0270000e-01,
		   3.9910000e-01,
		   3.9860000e-01,
		   3.9740000e-01,
		   3.9520000e-01,
		   3.9240000e-01,
		   3.9280000e-01,
		   3.9130000e-01,
		   3.8890000e-01,
		   3.8800000e-01,
		   3.8690000e-01,
		   3.8450000e-01,
		   3.8210000e-01,
		   3.8010000e-01,
		   3.8030000e-01,
		   3.7720000e-01,
		   3.7460000e-01,
		   3.7300000e-01,
		   3.7080000e-01,
		   3.6780000e-01,
		   3.6720000e-01,
		   3.6310000e-01,
		   3.5990000e-01,
		   3.5770000e-01,
		   3.5570000e-01,
		   3.5190000e-01,
		   3.4960000e-01,
		   3.4480000e-01,
		   3.4370000e-01,
		   3.4070000e-01,
		   3.3670000e-01,
		   3.3350000e-01,
		   3.3010000e-01,
		   3.2480000e-01,
		   3.2560000e-01,
		   3.2320000e-01,
		   3.1880000e-01,
		   3.1460000e-01,
		   3.1410000e-01,
		   3.0980000e-01,
		   3.0570000e-01,
		   3.0100000e-01,
		   2.9830000e-01,
		   2.9430000e-01,
		   2.9090000e-01,
		   2.8670000e-01,
		   2.8400000e-01,
		   2.7810000e-01,
		   2.7590000e-01,
		   2.7220000e-01,
		   2.6800000e-01,
		   2.6370000e-01,
		   2.6220000e-01,
		   2.5790000e-01,
		   2.5160000e-01,
		   2.4700000e-01,
		   2.4490000e-01,
		   2.3920000e-01,
		   2.3450000e-01,
		   2.2830000e-01,
		   2.2470000e-01,
		   2.2190000e-01,
		   2.2460000e-01,
		   2.1960000e-01,
		   2.2020000e-01,
		   2.1340000e-01,
		   2.0860000e-01,
		   2.0670000e-01,
		   2.0600000e-01,
		   2.0330000e-01,
		   2.0120000e-01,
		   1.9780000e-01,
		   1.9510000e-01,
		   1.9230000e-01,
		   1.9110000e-01,
		   1.8870000e-01,
		   1.8570000e-01,
		   1.8190000e-01,
		   1.8160000e-01,
		   1.8270000e-01,
		   1.7870000e-01,
		   1.7750000e-01,
		   1.7800000e-01,
		   1.7480000e-01,
		   1.7400000e-01,
		   1.7260000e-01,
		   1.7160000e-01,
		   1.7010000e-01,
		   1.7000000e-01,
		   1.6760000e-01,
		   1.6640000e-01,
		   1.6490000e-01,
		   1.6430000e-01,
		   1.6280000e-01,
		   1.6030000e-01,
		   1.5830000e-01,
		   1.5920000e-01,
		   1.5580000e-01,
		   1.5400000e-01,
		   1.5420000e-01,
		   1.5500000e-01,
		   1.5210000e-01,
		   1.5010000e-01,
		   1.5130000e-01,
		   1.5010000e-01,
		   1.4850000e-01,
		   1.4840000e-01,
		   1.4700000e-01,
		   1.4650000e-01,
		   1.4680000e-01,
		   1.4660000e-01,
		   1.4690000e-01,
		   1.4550000e-01,
		   1.4410000e-01,
		   1.4580000e-01,
		   1.4570000e-01,
		   1.4260000e-01,
		   1.4370000e-01,
		   1.4320000e-01,
		   1.4180000e-01,
		   1.4330000e-01,
		   1.4180000e-01,
		   1.4000000e-01,
		   1.4310000e-01,
		   1.4300000e-01,
		   1.4060000e-01,
		   1.3950000e-01,
		   1.3890000e-01,
		   1.3960000e-01,
		   1.3930000e-01,
		   1.3800000e-01,
		   1.3780000e-01,
		   1.3900000e-01,
		   1.3760000e-01,
		   1.3770000e-01,
		   1.3900000e-01,
		   1.3820000e-01,
		   1.3650000e-01,
		   1.3780000e-01,
		   1.3870000e-01,
		   1.3740000e-01,
		   1.3710000e-01,
		   1.3760000e-01,
		   1.3680000e-01,
		   1.3760000e-01,
		   1.3820000e-01,
		   1.3740000e-01,
		   1.3680000e-01,
		   1.3710000e-01,
		   1.3820000e-01,
		   1.3850000e-01,
		   1.3780000e-01,
		   1.3880000e-01,
		   1.3960000e-01,
		   1.3850000e-01,
		   1.3700000e-01,
		   1.3970000e-01,
		   1.3910000e-01,
		   1.3830000e-01,
		   1.4030000e-01,
		   1.4140000e-01,
		   1.4070000e-01,
		   1.4110000e-01,
		   1.4120000e-01,
		   1.4190000e-01,
		   1.4060000e-01,
		   1.4000000e-01,
		   1.3880000e-01,
		   1.3910000e-01,
		   1.3950000e-01,
		   1.4210000e-01,
		   1.4340000e-01,
		   1.4270000e-01,
		   1.4400000e-01,
		   1.4450000e-01,
		   1.4150000e-01,
		   1.4190000e-01,
		   1.4150000e-01,
		   1.4080000e-01,
		   1.4140000e-01,
		   1.4210000e-01,
		   1.4130000e-01,
		   1.4040000e-01,
		   1.4010000e-01,
		   1.3990000e-01,
		   1.4070000e-01,
		   1.3820000e-01,
		   1.3530000e-01,
		   1.3610000e-01,
		   1.3930000e-01,
		   1.4270000e-01,
		   1.4010000e-01,
		   1.3860000e-01,
		   1.3860000e-01,
		   1.3990000e-01,
		   1.3790000e-01,
		   1.3650000e-01,
		   1.3780000e-01,
		   1.3750000e-01,
		   1.3580000e-01,
		   1.3670000e-01
	#	   1.3620000e-01,
	#	   1.3500000e-01,
	#	   1.3460000e-01,
	#	   1.3450000e-01,
	#	   1.3380000e-01,
	#	   1.3340000e-01,
	#	   1.3300000e-01,
	#	   1.3410000e-01,
	#	   1.3460000e-01,
	#	   1.3280000e-01,
	#	   1.3520000e-01,
	#	   1.3520000e-01,
	#	   1.3490000e-01,
	#	   1.3620000e-01,
	#	   1.3650000e-01,
	#	   1.3480000e-01,
	#	   1.3600000e-01,
	#	   1.3780000e-01,
	#	   1.3820000e-01,
	#	   1.3800000e-01,
	#	   1.3960000e-01,
	#	   1.3990000e-01,
	#	   1.4110000e-01,
	#	   1.4100000e-01,
	#	   1.4340000e-01,
	#	   1.4520000e-01,
	#	   1.4550000e-01,
	#	   1.4670000e-01,
	#	   1.5050000e-01,
	#	   1.5100000e-01,
	#	   1.5100000e-01,
	#	   1.5370000e-01,
	#	   1.5800000e-01,
	#	   1.5840000e-01,
	#	   1.6070000e-01,
	#	   1.6360000e-01,
	#	   1.6600000e-01,
	#	   1.6790000e-01,
	#	   1.7260000e-01,
	#	   1.7670000e-01,
	#	   1.8040000e-01,
	#	   1.8320000e-01,
	#	   1.8760000e-01,
	#	   1.9150000e-01,
	#	   1.9530000e-01,
	#	   1.9780000e-01,
	#	   2.0330000e-01,
	#	   2.0510000e-01,
	#	   2.0890000e-01,
	#	   2.1480000e-01,
	#	   2.1980000e-01,
	#	   2.2160000e-01,
	#	   2.2560000e-01,
	#	   2.2820000e-01,
	#	   2.3280000e-01,
	#	   2.3630000e-01,
	#	   2.3940000e-01,
	#	   2.4270000e-01,
	#	   2.4620000e-01,
	#	   2.4680000e-01,
	#	   2.4850000e-01,
	#	   2.5070000e-01,
	#	   2.5300000e-01,
	#	   2.5660000e-01,
	#	   2.5760000e-01,
	#	   2.5740000e-01,
	#	   2.6170000e-01,
	#	   2.6400000e-01,
	#	   2.6650000e-01,
	#	   2.6670000e-01,
	#	   2.6770000e-01,
	#	   2.6640000e-01,
	#	   2.6680000e-01,
	#	   2.6790000e-01,
	#	   2.6980000e-01,
	#	   2.6940000e-01,
	#	   2.6880000e-01,
	#	   2.6920000e-01,
	#	   2.7000000e-01,
	#	   2.7000000e-01,
	#	   2.7090000e-01,
	#	   2.7220000e-01,
	#	   2.7380000e-01,
	#	   2.7480000e-01,
	#	   2.7470000e-01,
	#	   2.7580000e-01,
	#	   2.7920000e-01,
	#	   2.8110000e-01,
	#	   2.8000000e-01,
	#	   2.7780000e-01,
	#	   2.7990000e-01,
	#	   2.8610000e-01,
	#	   2.8300000e-01,
	#	   2.8020000e-01,
	#	   2.8210000e-01,
	#	   2.8390000e-01,
	#	   2.8630000e-01,
	#	   2.8850000e-01,
	#	   2.9240000e-01,
	#	   2.9260000e-01
		])

		m5b7_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   1.3740000e-01,
		   1.3940000e-01,
		   1.4420000e-01,
		   1.5350000e-01,
		   1.5790000e-01,
		   1.6730000e-01,
		   1.7420000e-01,
		   1.8450000e-01,
		   1.9270000e-01,
		   2.0660000e-01,
		   2.1520000e-01,
		   2.3140000e-01,
		   2.4230000e-01,
		   2.5830000e-01,
		   2.6820000e-01,
		   2.8490000e-01,
		   2.9900000e-01,
		   3.1830000e-01,
		   3.2830000e-01,
		   3.4780000e-01,
		   3.6040000e-01,
		   3.7630000e-01,
		   3.8700000e-01,
		   4.0250000e-01,
		   4.1090000e-01,
		   4.2260000e-01,
		   4.3090000e-01,
		   4.4000000e-01,
		   4.4430000e-01,
		   4.5100000e-01,
		   4.5730000e-01,
		   4.6440000e-01,
		   4.6600000e-01,
		   4.6830000e-01,
		   4.7000000e-01,
		   4.7150000e-01,
		   4.7090000e-01,
		   4.7480000e-01,
		   4.7630000e-01,
		   4.7760000e-01,
		   4.8120000e-01,
		   4.8430000e-01,
		   4.8270000e-01,
		   4.8390000e-01,
		   4.8680000e-01,
		   4.8700000e-01,
		   4.8510000e-01,
		   4.8910000e-01,
		   4.9090000e-01,
		   4.9150000e-01,
		   4.9230000e-01,
		   4.9630000e-01,
		   4.9930000e-01,
		   4.9990000e-01,
		   5.0030000e-01,
		   5.0380000e-01,
		   5.0440000e-01,
		   5.0470000e-01,
		   5.0840000e-01,
		   5.1110000e-01,
		   5.1070000e-01,
		   5.1370000e-01,
		   5.1660000e-01,
		   5.1800000e-01,
		   5.1910000e-01,
		   5.1900000e-01,
		   5.2040000e-01,
		   5.2150000e-01,
		   5.2230000e-01,
		   5.2570000e-01,
		   5.2620000e-01,
		   5.2470000e-01,
		   5.2760000e-01,
		   5.3110000e-01,
		   5.3040000e-01,
		   5.3190000e-01,
		   5.3570000e-01,
		   5.3760000e-01,
		   5.3730000e-01,
		   5.4120000e-01,
		   5.4350000e-01,
		   5.4620000e-01,
		   5.4740000e-01,
		   5.4780000e-01,
		   5.4900000e-01,
		   5.4990000e-01,
		   5.5150000e-01,
		   5.5310000e-01,
		   5.5350000e-01,
		   5.5530000e-01,
		   5.5750000e-01,
		   5.5690000e-01,
		   5.5690000e-01,
		   5.5780000e-01,
		   5.5950000e-01,
		   5.5850000e-01,
		   5.5850000e-01,
		   5.6000000e-01,
		   5.5980000e-01,
		   5.5950000e-01,
		   5.5950000e-01,
		   5.5980000e-01,
		   5.5900000e-01,
		   5.5740000e-01,
		   5.5870000e-01,
		   5.6080000e-01,
		   5.5780000e-01,
		   5.5560000e-01,
		   5.5760000e-01,
		   5.5580000e-01,
		   5.5450000e-01,
		   5.5580000e-01,
		   5.5500000e-01,
		   5.5230000e-01,
		   5.5300000e-01,
		   5.5110000e-01,
		   5.5010000e-01,
		   5.5000000e-01,
		   5.4960000e-01,
		   5.4680000e-01,
		   5.4610000e-01,
		   5.4490000e-01,
		   5.4620000e-01,
		   5.4440000e-01,
		   5.4080000e-01,
		   5.4070000e-01,
		   5.3990000e-01,
		   5.3520000e-01,
		   5.3440000e-01,
		   5.3460000e-01,
		   5.3220000e-01,
		   5.2910000e-01,
		   5.2930000e-01,
		   5.2450000e-01,
		   5.2160000e-01,
		   5.1770000e-01,
		   5.1760000e-01,
		   5.1570000e-01,
		   5.1290000e-01,
		   5.1080000e-01,
		   5.0980000e-01,
		   5.0470000e-01,
		   5.0120000e-01,
		   4.9760000e-01,
		   4.9520000e-01,
		   4.8950000e-01,
		   4.8770000e-01,
		   4.8440000e-01,
		   4.8030000e-01,
		   4.7630000e-01,
		   4.7570000e-01,
		   4.7070000e-01,
		   4.6580000e-01,
		   4.6090000e-01,
		   4.5980000e-01,
		   4.5560000e-01,
		   4.5150000e-01,
		   4.4560000e-01,
		   4.4250000e-01,
		   4.3790000e-01,
		   4.3410000e-01,
		   4.2890000e-01,
		   4.2460000e-01,
		   4.1940000e-01,
		   4.1840000e-01,
		   4.1100000e-01,
		   4.0500000e-01,
		   3.9980000e-01,
		   3.9720000e-01,
		   3.9190000e-01,
		   3.8830000e-01,
		   3.8120000e-01,
		   3.7720000e-01,
		   3.7100000e-01,
		   3.6720000e-01,
		   3.6130000e-01,
		   3.5620000e-01,
		   3.4820000e-01,
		   3.4560000e-01,
		   3.3740000e-01,
		   3.3040000e-01,
		   3.2840000e-01,
		   3.2890000e-01,
		   3.2250000e-01,
		   3.2080000e-01,
		   3.1500000e-01,
		   3.1170000e-01,
		   3.0830000e-01,
		   3.0440000e-01,
		   2.9900000e-01,
		   2.9770000e-01,
		   2.9480000e-01,
		   2.9100000e-01,
		   2.8560000e-01,
		   2.8440000e-01,
		   2.7990000e-01,
		   2.7790000e-01,
		   2.7740000e-01,
		   2.7510000e-01,
		   2.7190000e-01,
		   2.6880000e-01,
		   2.6510000e-01,
		   2.6360000e-01,
		   2.6040000e-01,
		   2.5960000e-01,
		   2.5870000e-01,
		   2.5670000e-01,
		   2.5390000e-01,
		   2.5530000e-01,
		   2.5520000e-01,
		   2.5090000e-01,
		   2.4980000e-01,
		   2.4890000e-01,
		   2.4500000e-01,
		   2.4310000e-01,
		   2.3880000e-01,
		   2.3730000e-01,
		   2.3550000e-01,
		   2.3560000e-01,
		   2.3430000e-01,
		   2.3320000e-01,
		   2.3170000e-01,
		   2.3020000e-01,
		   2.2740000e-01,
		   2.2480000e-01,
		   2.2420000e-01,
		   2.2450000e-01,
		   2.2270000e-01,
		   2.2220000e-01,
		   2.2090000e-01,
		   2.2010000e-01,
		   2.1900000e-01,
		   2.1860000e-01,
		   2.1890000e-01,
		   2.1760000e-01,
		   2.1610000e-01,
		   2.1660000e-01,
		   2.1730000e-01,
		   2.1620000e-01,
		   2.1560000e-01,
		   2.1650000e-01,
		   2.1550000e-01,
		   2.1420000e-01,
		   2.1590000e-01,
		   2.1630000e-01,
		   2.1540000e-01,
		   2.1500000e-01,
		   2.1320000e-01,
		   2.1240000e-01,
		   2.1320000e-01,
		   2.1460000e-01,
		   2.1330000e-01,
		   2.1120000e-01,
		   2.1230000e-01,
		   2.1290000e-01,
		   2.1300000e-01,
		   2.1110000e-01,
		   2.1270000e-01,
		   2.1430000e-01,
		   2.1300000e-01,
		   2.1130000e-01,
		   2.1180000e-01,
		   2.1170000e-01,
		   2.1210000e-01,
		   2.1420000e-01,
		   2.1570000e-01,
		   2.1520000e-01,
		   2.1510000e-01,
		   2.1530000e-01,
		   2.1540000e-01,
		   2.1540000e-01,
		   2.1420000e-01,
		   2.1630000e-01,
		   2.1690000e-01,
		   2.1670000e-01,
		   2.1810000e-01,
		   2.1870000e-01,
		   2.1800000e-01,
		   2.1880000e-01,
		   2.1980000e-01,
		   2.1960000e-01,
		   2.2090000e-01,
		   2.2210000e-01,
		   2.2070000e-01,
		   2.1960000e-01,
		   2.1850000e-01,
		   2.1900000e-01,
		   2.1760000e-01,
		   2.1720000e-01,
		   2.1860000e-01,
		   2.2080000e-01,
		   2.2130000e-01,
		   2.2140000e-01,
		   2.2170000e-01,
		   2.2150000e-01,
		   2.1970000e-01,
		   2.2190000e-01,
		   2.2120000e-01,
		   2.1940000e-01,
		   2.2000000e-01,
		   2.2190000e-01,
		   2.2050000e-01,
		   2.1890000e-01,
		   2.1820000e-01,
		   2.1870000e-01,
		   2.1730000e-01,
		   2.1620000e-01,
		   2.1410000e-01,
		   2.1250000e-01,
		   2.1220000e-01,
		   2.1690000e-01,
		   2.1670000e-01,
		   2.1610000e-01,
		   2.1540000e-01,
		   2.1430000e-01,
		   2.1360000e-01,
		   2.1340000e-01,
		   2.1280000e-01,
		   2.1140000e-01,
		   2.1130000e-01,
		   2.1200000e-01
	#	   2.0900000e-01,
	#	   2.0910000e-01,
	#	   2.0810000e-01,
	#	   2.0750000e-01,
	#	   2.0750000e-01,
	#	   2.0740000e-01,
	#	   2.0780000e-01,
	#	   2.0780000e-01,
	#	   2.0710000e-01,
	#	   2.0940000e-01,
	#	   2.0960000e-01,
	#	   2.0890000e-01,
	#	   2.0980000e-01,
	#	   2.1110000e-01,
	#	   2.1110000e-01,
	#	   2.1110000e-01,
	#	   2.1180000e-01,
	#	   2.1300000e-01,
	#	   2.1480000e-01,
	#	   2.1620000e-01,
	#	   2.1820000e-01,
	#	   2.1790000e-01,
	#	   2.2110000e-01,
	#	   2.2440000e-01,
	#	   2.2600000e-01,
	#	   2.2490000e-01,
	#	   2.2730000e-01,
	#	   2.2940000e-01,
	#	   2.3260000e-01,
	#	   2.3440000e-01,
	#	   2.3740000e-01,
	#	   2.3970000e-01,
	#	   2.4450000e-01,
	#	   2.4760000e-01,
	#	   2.5280000e-01,
	#	   2.5440000e-01,
	#	   2.6030000e-01,
	#	   2.6540000e-01,
	#	   2.7170000e-01,
	#	   2.7300000e-01,
	#	   2.7710000e-01,
	#	   2.8080000e-01,
	#	   2.8610000e-01,
	#	   2.8970000e-01,
	#	   2.9610000e-01,
	#	   3.0130000e-01,
	#	   3.1020000e-01,
	#	   3.1280000e-01,
	#	   3.2100000e-01,
	#	   3.2620000e-01,
	#	   3.3110000e-01,
	#	   3.3690000e-01,
	#	   3.4340000e-01,
	#	   3.4420000e-01,
	#	   3.4880000e-01,
	#	   3.5370000e-01,
	#	   3.5940000e-01,
	#	   3.6150000e-01,
	#	   3.6440000e-01,
	#	   3.6740000e-01,
	#	   3.7160000e-01,
	#	   3.7280000e-01,
	#	   3.7620000e-01,
	#	   3.7950000e-01,
	#	   3.8360000e-01,
	#	   3.8360000e-01,
	#	   3.8560000e-01,
	#	   3.8710000e-01,
	#	   3.8790000e-01,
	#	   3.8970000e-01,
	#	   3.9130000e-01,
	#	   3.8840000e-01,
	#	   3.8950000e-01,
	#	   3.9300000e-01,
	#	   3.9190000e-01,
	#	   3.9000000e-01,
	#	   3.9300000e-01,
	#	   3.9290000e-01,
	#	   3.9360000e-01,
	#	   3.9700000e-01,
	#	   3.9990000e-01,
	#	   3.9880000e-01,
	#	   3.9900000e-01,
	#	   3.9760000e-01,
	#	   3.9670000e-01,
	#	   3.9700000e-01,
	#	   3.9970000e-01,
	#	   4.0530000e-01,
	#	   4.0750000e-01,
	#	   4.0810000e-01,
	#	   4.1010000e-01,
	#	   4.1180000e-01,
	#	   4.0900000e-01,
	#	   4.1010000e-01,
	#	   4.1030000e-01,
	#	   4.1370000e-01,
	#	   4.1770000e-01,
	#	   4.1870000e-01,
	#	   4.2150000e-01,
	#	   4.2100000e-01
		])

		# 7.5B5-8/4
		m75b5_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   1.0600000e-01,
		   1.1030000e-01,
		   1.1210000e-01,
		   1.1460000e-01,
		   1.1770000e-01,
		   1.2290000e-01,
		   1.2820000e-01,
		   1.3580000e-01,
		   1.4250000e-01,
		   1.4990000e-01,
		   1.5300000e-01,
		   1.5970000e-01,
		   1.6710000e-01,
		   1.7480000e-01,
		   1.7860000e-01,
		   1.8680000e-01,
		   1.9340000e-01,
		   1.9930000e-01,
		   2.0340000e-01,
		   2.1080000e-01,
		   2.1410000e-01,
		   2.1790000e-01,
		   2.2220000e-01,
		   2.2590000e-01,
		   2.2710000e-01,
		   2.2880000e-01,
		   2.3090000e-01,
		   2.3210000e-01,
		   2.3290000e-01,
		   2.3440000e-01,
		   2.3710000e-01,
		   2.3770000e-01,
		   2.3710000e-01,
		   2.3820000e-01,
		   2.3810000e-01,
		   2.3570000e-01,
		   2.3650000e-01,
		   2.3800000e-01,
		   2.3650000e-01,
		   2.3780000e-01,
		   2.3960000e-01,
		   2.3870000e-01,
		   2.3820000e-01,
		   2.3950000e-01,
		   2.3970000e-01,
		   2.4020000e-01,
		   2.3860000e-01,
		   2.3900000e-01,
		   2.4030000e-01,
		   2.3920000e-01,
		   2.3960000e-01,
		   2.4140000e-01,
		   2.4020000e-01,
		   2.4020000e-01,
		   2.4210000e-01,
		   2.4300000e-01,
		   2.4080000e-01,
		   2.4020000e-01,
		   2.4180000e-01,
		   2.4240000e-01,
		   2.4130000e-01,
		   2.4190000e-01,
		   2.4210000e-01,
		   2.4230000e-01,
		   2.4340000e-01,
		   2.4270000e-01,
		   2.4220000e-01,
		   2.4220000e-01,
		   2.4260000e-01,
		   2.4360000e-01,
		   2.4240000e-01,
		   2.4060000e-01,
		   2.4250000e-01,
		   2.4140000e-01,
		   2.4050000e-01,
		   2.4170000e-01,
		   2.4190000e-01,
		   2.4060000e-01,
		   2.4180000e-01,
		   2.4160000e-01,
		   2.4000000e-01,
		   2.4080000e-01,
		   2.4260000e-01,
		   2.4080000e-01,
		   2.4010000e-01,
		   2.3920000e-01,
		   2.4060000e-01,
		   2.4060000e-01,
		   2.3820000e-01,
		   2.3840000e-01,
		   2.4030000e-01,
		   2.3900000e-01,
		   2.3810000e-01,
		   2.3860000e-01,
		   2.3780000e-01,
		   2.3650000e-01,
		   2.3750000e-01,
		   2.3760000e-01,
		   2.3530000e-01,
		   2.3600000e-01,
		   2.3620000e-01,
		   2.3700000e-01,
		   2.3540000e-01,
		   2.3540000e-01,
		   2.3410000e-01,
		   2.3360000e-01,
		   2.3220000e-01,
		   2.3220000e-01,
		   2.3130000e-01,
		   2.3160000e-01,
		   2.3130000e-01,
		   2.3230000e-01,
		   2.2960000e-01,
		   2.2930000e-01,
		   2.3240000e-01,
		   2.3130000e-01,
		   2.2900000e-01,
		   2.3070000e-01,
		   2.3110000e-01,
		   2.2940000e-01,
		   2.2920000e-01,
		   2.2760000e-01,
		   2.2720000e-01,
		   2.2790000e-01,
		   2.2770000e-01,
		   2.2720000e-01,
		   2.2620000e-01,
		   2.2370000e-01,
		   2.2550000e-01,
		   2.2600000e-01,
		   2.2360000e-01,
		   2.2210000e-01,
		   2.2380000e-01,
		   2.2220000e-01,
		   2.2040000e-01,
		   2.2130000e-01,
		   2.2130000e-01,
		   2.1940000e-01,
		   2.1940000e-01,
		   2.1940000e-01,
		   2.1810000e-01,
		   2.1670000e-01,
		   2.1640000e-01,
		   2.1490000e-01,
		   2.1240000e-01,
		   2.1190000e-01,
		   2.1150000e-01,
		   2.1060000e-01,
		   2.0980000e-01,
		   2.0840000e-01,
		   2.0830000e-01,
		   2.0700000e-01,
		   2.0490000e-01,
		   2.0540000e-01,
		   2.0380000e-01,
		   2.0190000e-01,
		   2.0190000e-01,
		   2.0050000e-01,
		   1.9730000e-01,
		   1.9760000e-01,
		   1.9690000e-01,
		   1.9320000e-01,
		   1.9110000e-01,
		   1.9070000e-01,
		   1.9000000e-01,
		   1.8770000e-01,
		   1.8510000e-01,
		   1.8320000e-01,
		   1.8350000e-01,
		   1.8120000e-01,
		   1.8050000e-01,
		   1.8020000e-01,
		   1.7740000e-01,
		   1.7310000e-01,
		   1.7250000e-01,
		   1.7050000e-01,
		   1.6750000e-01,
		   1.6410000e-01,
		   1.6170000e-01,
		   1.5790000e-01,
		   1.5500000e-01,
		   1.5700000e-01,
		   1.5940000e-01,
		   1.5810000e-01,
		   1.5580000e-01,
		   1.5210000e-01,
		   1.5240000e-01,
		   1.5100000e-01,
		   1.4840000e-01,
		   1.4710000e-01,
		   1.4740000e-01,
		   1.4420000e-01,
		   1.4300000e-01,
		   1.4130000e-01,
		   1.3980000e-01,
		   1.3740000e-01,
		   1.3610000e-01,
		   1.3550000e-01,
		   1.3440000e-01,
		   1.3410000e-01,
		   1.3460000e-01,
		   1.3350000e-01,
		   1.3220000e-01,
		   1.2990000e-01,
		   1.3130000e-01,
		   1.3010000e-01,
		   1.2890000e-01,
		   1.2880000e-01,
		   1.2870000e-01,
		   1.2550000e-01,
		   1.2410000e-01,
		   1.2450000e-01,
		   1.2370000e-01,
		   1.2240000e-01,
		   1.2340000e-01,
		   1.2170000e-01,
		   1.2020000e-01,
		   1.1950000e-01,
		   1.2090000e-01,
		   1.2060000e-01,
		   1.1940000e-01,
		   1.1800000e-01,
		   1.1950000e-01,
		   1.1920000e-01,
		   1.1780000e-01,
		   1.1710000e-01,
		   1.1760000e-01,
		   1.1570000e-01,
		   1.1550000e-01,
		   1.1590000e-01,
		   1.1490000e-01,
		   1.1490000e-01,
		   1.1500000e-01,
		   1.1450000e-01,
		   1.1390000e-01,
		   1.1290000e-01,
		   1.1370000e-01,
		   1.1420000e-01,
		   1.1340000e-01,
		   1.1350000e-01,
		   1.1430000e-01,
		   1.1460000e-01,
		   1.1480000e-01,
		   1.1520000e-01,
		   1.1360000e-01,
		   1.1240000e-01,
		   1.1280000e-01,
		   1.1290000e-01,
		   1.1190000e-01,
		   1.1170000e-01,
		   1.1320000e-01,
		   1.1240000e-01,
		   1.1210000e-01,
		   1.1300000e-01,
		   1.1320000e-01,
		   1.1390000e-01,
		   1.1450000e-01,
		   1.1390000e-01,
		   1.1410000e-01,
		   1.1380000e-01,
		   1.1410000e-01,
		   1.1450000e-01,
		   1.1290000e-01,
		   1.1290000e-01,
		   1.1510000e-01,
		   1.1620000e-01,
		   1.1400000e-01,
		   1.1510000e-01,
		   1.1570000e-01,
		   1.1510000e-01,
		   1.1560000e-01,
		   1.1660000e-01,
		   1.1630000e-01,
		   1.1660000e-01,
		   1.1730000e-01,
		   1.1830000e-01,
		   1.1760000e-01,
		   1.1700000e-01,
		   1.1800000e-01,
		   1.1960000e-01,
		   1.1890000e-01,
		   1.1800000e-01,
		   1.1860000e-01,
		   1.1720000e-01,
		   1.1690000e-01,
		   1.1670000e-01,
		   1.1540000e-01,
		   1.1440000e-01,
		   1.1640000e-01,
		   1.1660000e-01,
		   1.1590000e-01,
		   1.1760000e-01,
		   1.1950000e-01,
		   1.1810000e-01,
		   1.1600000e-01,
		   1.1590000e-01,
		   1.1670000e-01,
		   1.1650000e-01,
		   1.1530000e-01,
		   1.1500000e-01,
		   1.1490000e-01,
		   1.1360000e-01,
		   1.1190000e-01,
		   1.1200000e-01,
		   1.1070000e-01,
		   1.0900000e-01,
		   1.0980000e-01,
		   1.0780000e-01,
		   1.0370000e-01,
		   1.0680000e-01,
		   1.1060000e-01,
		   1.0880000e-01,
		   1.0990000e-01,
		   1.1050000e-01,
		   1.0850000e-01,
		   1.0820000e-01,
		   1.0840000e-01,
		   1.0780000e-01,
		   1.0690000e-01,
		   1.0780000e-01,
		   1.0730000e-01
	#	   1.0710000e-01,
	#	   1.0630000e-01,
	#	   1.0650000e-01,
	#	   1.0670000e-01,
	#	   1.0660000e-01,
	#	   1.0700000e-01,
	#	   1.0800000e-01,
	#	   1.0750000e-01,
	#	   1.0680000e-01,
	#	   1.0780000e-01,
	#	   1.0820000e-01,
	#	   1.0730000e-01,
	#	   1.0840000e-01,
	#	   1.0920000e-01,
	#	   1.0790000e-01,
	#	   1.1030000e-01,
	#	   1.1100000e-01,
	#	   1.1050000e-01,
	#	   1.1100000e-01,
	#	   1.1300000e-01,
	#	   1.1280000e-01,
	#	   1.1260000e-01,
	#	   1.1310000e-01,
	#	   1.1520000e-01,
	#	   1.1650000e-01,
	#	   1.1720000e-01,
	#	   1.1820000e-01,
	#	   1.2060000e-01,
	#	   1.2110000e-01,
	#	   1.2150000e-01,
	#	   1.2400000e-01,
	#	   1.2570000e-01,
	#	   1.2710000e-01,
	#	   1.2960000e-01,
	#	   1.3260000e-01,
	#	   1.3260000e-01,
	#	   1.3580000e-01,
	#	   1.3730000e-01,
	#	   1.3910000e-01,
	#	   1.4110000e-01,
	#	   1.4460000e-01,
	#	   1.4620000e-01,
	#	   1.4840000e-01,
	#	   1.5100000e-01,
	#	   1.5200000e-01,
	#	   1.5360000e-01,
	#	   1.5930000e-01,
	#	   1.6220000e-01,
	#	   1.6740000e-01,
	#	   1.6830000e-01,
	#	   1.7000000e-01,
	#	   1.7360000e-01,
	#	   1.7560000e-01,
	#	   1.7540000e-01,
	#	   1.7870000e-01,
	#	   1.8060000e-01,
	#	   1.8240000e-01,
	#	   1.8540000e-01,
	#	   1.8800000e-01,
	#	   1.8840000e-01,
	#	   1.9020000e-01,
	#	   1.9110000e-01,
	#	   1.9310000e-01,
	#	   1.9400000e-01,
	#	   1.9540000e-01,
	#	   1.9700000e-01,
	#	   1.9790000e-01,
	#	   1.9840000e-01,
	#	   1.9880000e-01,
	#	   1.9920000e-01,
	#	   1.9900000e-01,
	#	   1.9710000e-01,
	#	   1.9760000e-01,
	#	   1.9830000e-01,
	#	   1.9940000e-01,
	#	   1.9990000e-01,
	#	   2.0030000e-01,
	#	   2.0040000e-01,
	#	   1.9930000e-01,
	#	   2.0050000e-01,
	#	   2.0060000e-01,
	#	   2.0000000e-01,
	#	   2.0070000e-01,
	#	   2.0110000e-01,
	#	   2.0200000e-01,
	#	   2.0420000e-01,
	#	   2.0330000e-01,
	#	   2.0540000e-01,
	#	   2.0390000e-01,
	#	   2.0090000e-01,
	#	   2.0190000e-01,
	#	   2.0210000e-01,
	#	   2.0180000e-01,
	#	   2.0450000e-01,
	#	   2.0540000e-01,
	#	   2.0580000e-01,
	#	   2.0630000e-01,
	#	   2.0540000e-01,
	#	   2.0420000e-01,
	#	   2.0490000e-01
		])

		m75b6_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   1.1270000e-01,
		   1.1810000e-01,
		   1.2130000e-01,
		   1.2540000e-01,
		   1.3100000e-01,
		   1.3940000e-01,
		   1.4490000e-01,
		   1.5260000e-01,
		   1.6040000e-01,
		   1.7200000e-01,
		   1.7800000e-01,
		   1.8970000e-01,
		   1.9840000e-01,
		   2.1030000e-01,
		   2.1790000e-01,
		   2.2990000e-01,
		   2.3830000e-01,
		   2.5230000e-01,
		   2.6060000e-01,
		   2.7200000e-01,
		   2.7840000e-01,
		   2.8660000e-01,
		   2.9420000e-01,
		   3.0290000e-01,
		   3.0450000e-01,
		   3.0970000e-01,
		   3.1610000e-01,
		   3.1930000e-01,
		   3.1990000e-01,
		   3.2420000e-01,
		   3.2760000e-01,
		   3.2920000e-01,
		   3.2880000e-01,
		   3.2930000e-01,
		   3.2920000e-01,
		   3.2880000e-01,
		   3.2880000e-01,
		   3.3040000e-01,
		   3.3080000e-01,
		   3.3160000e-01,
		   3.3260000e-01,
		   3.3220000e-01,
		   3.3070000e-01,
		   3.3060000e-01,
		   3.3330000e-01,
		   3.3260000e-01,
		   3.3040000e-01,
		   3.3210000e-01,
		   3.3280000e-01,
		   3.3220000e-01,
		   3.3370000e-01,
		   3.3300000e-01,
		   3.3230000e-01,
		   3.3510000e-01,
		   3.3560000e-01,
		   3.3570000e-01,
		   3.3530000e-01,
		   3.3480000e-01,
		   3.3510000e-01,
		   3.3580000e-01,
		   3.3490000e-01,
		   3.3630000e-01,
		   3.3650000e-01,
		   3.3400000e-01,
		   3.3450000e-01,
		   3.3620000e-01,
		   3.3420000e-01,
		   3.3310000e-01,
		   3.3510000e-01,
		   3.3540000e-01,
		   3.3440000e-01,
		   3.3400000e-01,
		   3.3640000e-01,
		   3.3570000e-01,
		   3.3540000e-01,
		   3.3510000e-01,
		   3.3450000e-01,
		   3.3320000e-01,
		   3.3430000e-01,
		   3.3420000e-01,
		   3.3320000e-01,
		   3.3260000e-01,
		   3.3440000e-01,
		   3.3510000e-01,
		   3.3270000e-01,
		   3.3090000e-01,
		   3.3310000e-01,
		   3.3190000e-01,
		   3.2990000e-01,
		   3.3190000e-01,
		   3.3260000e-01,
		   3.3240000e-01,
		   3.3250000e-01,
		   3.3160000e-01,
		   3.3090000e-01,
		   3.2920000e-01,
		   3.2940000e-01,
		   3.2950000e-01,
		   3.2770000e-01,
		   3.2630000e-01,
		   3.2720000e-01,
		   3.2810000e-01,
		   3.2600000e-01,
		   3.2670000e-01,
		   3.2800000e-01,
		   3.2640000e-01,
		   3.2560000e-01,
		   3.2680000e-01,
		   3.2430000e-01,
		   3.2260000e-01,
		   3.2280000e-01,
		   3.2340000e-01,
		   3.2170000e-01,
		   3.2110000e-01,
		   3.2260000e-01,
		   3.2240000e-01,
		   3.2030000e-01,
		   3.1970000e-01,
		   3.2140000e-01,
		   3.2100000e-01,
		   3.1990000e-01,
		   3.1960000e-01,
		   3.1890000e-01,
		   3.1870000e-01,
		   3.1880000e-01,
		   3.1870000e-01,
		   3.1640000e-01,
		   3.1540000e-01,
		   3.1670000e-01,
		   3.1640000e-01,
		   3.1390000e-01,
		   3.1300000e-01,
		   3.1330000e-01,
		   3.1340000e-01,
		   3.1070000e-01,
		   3.0930000e-01,
		   3.0980000e-01,
		   3.0850000e-01,
		   3.0670000e-01,
		   3.0660000e-01,
		   3.0560000e-01,
		   3.0310000e-01,
		   3.0340000e-01,
		   3.0250000e-01,
		   2.9990000e-01,
		   2.9970000e-01,
		   2.9910000e-01,
		   2.9750000e-01,
		   2.9660000e-01,
		   2.9490000e-01,
		   2.9490000e-01,
		   2.9340000e-01,
		   2.9110000e-01,
		   2.8860000e-01,
		   2.8770000e-01,
		   2.8610000e-01,
		   2.8590000e-01,
		   2.8620000e-01,
		   2.8400000e-01,
		   2.8190000e-01,
		   2.8140000e-01,
		   2.7800000e-01,
		   2.7610000e-01,
		   2.7520000e-01,
		   2.7510000e-01,
		   2.7160000e-01,
		   2.6930000e-01,
		   2.6820000e-01,
		   2.6720000e-01,
		   2.6400000e-01,
		   2.6280000e-01,
		   2.6050000e-01,
		   2.5870000e-01,
		   2.5610000e-01,
		   2.5460000e-01,
		   2.5180000e-01,
		   2.4760000e-01,
		   2.4450000e-01,
		   2.4260000e-01,
		   2.3700000e-01,
		   2.3200000e-01,
		   2.3340000e-01,
		   2.3530000e-01,
		   2.3390000e-01,
		   2.3180000e-01,
		   2.2780000e-01,
		   2.2750000e-01,
		   2.2500000e-01,
		   2.2160000e-01,
		   2.1980000e-01,
		   2.2050000e-01,
		   2.1700000e-01,
		   2.1460000e-01,
		   2.1360000e-01,
		   2.1310000e-01,
		   2.1080000e-01,
		   2.0940000e-01,
		   2.0840000e-01,
		   2.0740000e-01,
		   2.0540000e-01,
		   2.0460000e-01,
		   2.0310000e-01,
		   2.0140000e-01,
		   1.9990000e-01,
		   1.9970000e-01,
		   1.9710000e-01,
		   1.9580000e-01,
		   1.9550000e-01,
		   1.9620000e-01,
		   1.9400000e-01,
		   1.9260000e-01,
		   1.9300000e-01,
		   1.9140000e-01,
		   1.8900000e-01,
		   1.9050000e-01,
		   1.8860000e-01,
		   1.8590000e-01,
		   1.8490000e-01,
		   1.8530000e-01,
		   1.8470000e-01,
		   1.8310000e-01,
		   1.8060000e-01,
		   1.8060000e-01,
		   1.8040000e-01,
		   1.7880000e-01,
		   1.7830000e-01,
		   1.7810000e-01,
		   1.7630000e-01,
		   1.7640000e-01,
		   1.7790000e-01,
		   1.7570000e-01,
		   1.7460000e-01,
		   1.7550000e-01,
		   1.7400000e-01,
		   1.7280000e-01,
		   1.7310000e-01,
		   1.7340000e-01,
		   1.7260000e-01,
		   1.7170000e-01,
		   1.7190000e-01,
		   1.7180000e-01,
		   1.7130000e-01,
		   1.7220000e-01,
		   1.7170000e-01,
		   1.7050000e-01,
		   1.6820000e-01,
		   1.7000000e-01,
		   1.6980000e-01,
		   1.6870000e-01,
		   1.6970000e-01,
		   1.7220000e-01,
		   1.7120000e-01,
		   1.6960000e-01,
		   1.7040000e-01,
		   1.7030000e-01,
		   1.7090000e-01,
		   1.7060000e-01,
		   1.7050000e-01,
		   1.7070000e-01,
		   1.7120000e-01,
		   1.7130000e-01,
		   1.7240000e-01,
		   1.7130000e-01,
		   1.7220000e-01,
		   1.7410000e-01,
		   1.7490000e-01,
		   1.7250000e-01,
		   1.7400000e-01,
		   1.7460000e-01,
		   1.7610000e-01,
		   1.7760000e-01,
		   1.7730000e-01,
		   1.7590000e-01,
		   1.7740000e-01,
		   1.7840000e-01,
		   1.7970000e-01,
		   1.7860000e-01,
		   1.7910000e-01,
		   1.8010000e-01,
		   1.8260000e-01,
		   1.8160000e-01,
		   1.7940000e-01,
		   1.8150000e-01,
		   1.8130000e-01,
		   1.8030000e-01,
		   1.7950000e-01,
		   1.7940000e-01,
		   1.7960000e-01,
		   1.8100000e-01,
		   1.7820000e-01,
		   1.7730000e-01,
		   1.7760000e-01,
		   1.7960000e-01,
		   1.7880000e-01,
		   1.7800000e-01,
		   1.7810000e-01,
		   1.7840000e-01,
		   1.7670000e-01,
		   1.7520000e-01,
		   1.7570000e-01,
		   1.7590000e-01,
		   1.7450000e-01,
		   1.7400000e-01,
		   1.7520000e-01,
		   1.7400000e-01,
		   1.7050000e-01,
		   1.6960000e-01,
		   1.6690000e-01,
		   1.6270000e-01,
		   1.6430000e-01,
		   1.6840000e-01,
		   1.6750000e-01,
		   1.6780000e-01,
		   1.6540000e-01,
		   1.6460000e-01,
		   1.6360000e-01,
		   1.6300000e-01,
		   1.6180000e-01,
		   1.6170000e-01,
		   1.6290000e-01,
		   1.6180000e-01
	#	   1.6010000e-01,
	#	   1.5990000e-01,
	#	   1.6030000e-01,
	#	   1.5960000e-01,
	#	   1.5840000e-01,
	#	   1.5890000e-01,
	#	   1.5840000e-01,
	#	   1.5780000e-01,
	#	   1.5780000e-01,
	#	   1.6020000e-01,
	#	   1.6010000e-01,
	#	   1.6000000e-01,
	#	   1.6120000e-01,
	#	   1.6170000e-01,
	#	   1.6240000e-01,
	#	   1.6490000e-01,
	#	   1.6450000e-01,
	#	   1.6550000e-01,
	#	   1.6890000e-01,
	#	   1.7140000e-01,
	#	   1.6950000e-01,
	#	   1.6880000e-01,
	#	   1.7150000e-01,
	#	   1.7290000e-01,
	#	   1.7590000e-01,
	#	   1.7600000e-01,
	#	   1.7960000e-01,
	#	   1.8250000e-01,
	#	   1.8320000e-01,
	#	   1.8540000e-01,
	#	   1.8700000e-01,
	#	   1.8740000e-01,
	#	   1.9040000e-01,
	#	   1.9420000e-01,
	#	   1.9840000e-01,
	#	   1.9900000e-01,
	#	   2.0310000e-01,
	#	   2.0570000e-01,
	#	   2.1120000e-01,
	#	   2.1350000e-01,
	#	   2.1660000e-01,
	#	   2.1970000e-01,
	#	   2.2550000e-01,
	#	   2.2790000e-01,
	#	   2.3090000e-01,
	#	   2.3080000e-01,
	#	   2.3560000e-01,
	#	   2.4110000e-01,
	#	   2.4730000e-01,
	#	   2.4720000e-01,
	#	   2.5130000e-01,
	#	   2.5500000e-01,
	#	   2.5520000e-01,
	#	   2.5700000e-01,
	#	   2.6150000e-01,
	#	   2.6330000e-01,
	#	   2.6500000e-01,
	#	   2.6690000e-01,
	#	   2.7130000e-01,
	#	   2.7060000e-01,
	#	   2.7070000e-01,
	#	   2.7190000e-01,
	#	   2.7550000e-01,
	#	   2.7610000e-01,
	#	   2.7620000e-01,
	#	   2.7800000e-01,
	#	   2.7860000e-01,
	#	   2.7810000e-01,
	#	   2.7770000e-01,
	#	   2.7830000e-01,
	#	   2.8100000e-01,
	#	   2.8010000e-01,
	#	   2.8200000e-01,
	#	   2.8290000e-01,
	#	   2.8070000e-01,
	#	   2.7940000e-01,
	#	   2.8390000e-01,
	#	   2.8450000e-01,
	#	   2.8340000e-01,
	#	   2.8230000e-01,
	#	   2.8600000e-01,
	#	   2.8430000e-01,
	#	   2.8420000e-01,
	#	   2.8500000e-01,
	#	   2.8340000e-01,
	#	   2.8250000e-01,
	#	   2.8340000e-01,
	#	   2.8560000e-01,
	#	   2.8800000e-01,
	#	   2.8770000e-01,
	#	   2.8620000e-01,
	#	   2.8400000e-01,
	#	   2.7920000e-01,
	#	   2.7990000e-01,
	#	   2.8370000e-01,
	#	   2.8610000e-01,
	#	   2.8640000e-01,
	#	   2.8620000e-01,
	#	   2.8480000e-01,
	#	   2.8420000e-01
		])

		m75b7_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   1.2100000e-01,
		   1.2790000e-01,
		   1.3060000e-01,
		   1.3610000e-01,
		   1.4210000e-01,
		   1.5010000e-01,
		   1.5540000e-01,
		   1.6770000e-01,
		   1.7620000e-01,
		   1.8870000e-01,
		   1.9810000e-01,
		   2.1200000e-01,
		   2.2130000e-01,
		   2.3810000e-01,
		   2.5080000e-01,
		   2.6620000e-01,
		   2.7620000e-01,
		   2.9560000e-01,
		   3.0980000e-01,
		   3.2980000e-01,
		   3.3950000e-01,
		   3.5620000e-01,
		   3.6860000e-01,
		   3.8400000e-01,
		   3.9560000e-01,
		   4.0950000e-01,
		   4.1580000e-01,
		   4.2480000e-01,
		   4.3020000e-01,
		   4.3790000e-01,
		   4.4210000e-01,
		   4.4800000e-01,
		   4.4930000e-01,
		   4.5040000e-01,
		   4.5210000e-01,
		   4.5190000e-01,
		   4.5110000e-01,
		   4.5350000e-01,
		   4.5450000e-01,
		   4.5690000e-01,
		   4.5730000e-01,
		   4.5820000e-01,
		   4.6010000e-01,
		   4.6030000e-01,
		   4.5840000e-01,
		   4.5850000e-01,
		   4.5960000e-01,
		   4.5890000e-01,
		   4.5780000e-01,
		   4.5850000e-01,
		   4.6010000e-01,
		   4.6020000e-01,
		   4.6000000e-01,
		   4.6110000e-01,
		   4.6130000e-01,
		   4.6230000e-01,
		   4.6370000e-01,
		   4.6430000e-01,
		   4.6360000e-01,
		   4.6330000e-01,
		   4.6380000e-01,
		   4.6520000e-01,
		   4.6240000e-01,
		   4.6400000e-01,
		   4.6490000e-01,
		   4.6220000e-01,
		   4.6090000e-01,
		   4.6240000e-01,
		   4.6170000e-01,
		   4.6160000e-01,
		   4.6310000e-01,
		   4.6280000e-01,
		   4.6180000e-01,
		   4.6200000e-01,
		   4.6130000e-01,
		   4.6270000e-01,
		   4.6180000e-01,
		   4.5960000e-01,
		   4.6070000e-01,
		   4.6120000e-01,
		   4.5940000e-01,
		   4.6110000e-01,
		   4.6220000e-01,
		   4.6010000e-01,
		   4.6110000e-01,
		   4.6240000e-01,
		   4.6000000e-01,
		   4.6080000e-01,
		   4.6200000e-01,
		   4.6080000e-01,
		   4.6100000e-01,
		   4.6150000e-01,
		   4.6220000e-01,
		   4.6050000e-01,
		   4.5950000e-01,
		   4.5920000e-01,
		   4.6050000e-01,
		   4.6090000e-01,
		   4.6000000e-01,
		   4.5880000e-01,
		   4.5830000e-01,
		   4.5630000e-01,
		   4.5650000e-01,
		   4.5800000e-01,
		   4.5690000e-01,
		   4.5570000e-01,
		   4.5560000e-01,
		   4.5510000e-01,
		   4.5380000e-01,
		   4.5380000e-01,
		   4.5290000e-01,
		   4.5350000e-01,
		   4.5290000e-01,
		   4.5230000e-01,
		   4.5210000e-01,
		   4.5160000e-01,
		   4.5130000e-01,
		   4.5210000e-01,
		   4.4940000e-01,
		   4.5030000e-01,
		   4.5190000e-01,
		   4.5000000e-01,
		   4.4640000e-01,
		   4.4890000e-01,
		   4.4930000e-01,
		   4.4560000e-01,
		   4.4520000e-01,
		   4.4590000e-01,
		   4.4450000e-01,
		   4.4380000e-01,
		   4.4260000e-01,
		   4.4170000e-01,
		   4.4090000e-01,
		   4.3780000e-01,
		   4.3770000e-01,
		   4.3780000e-01,
		   4.3540000e-01,
		   4.3250000e-01,
		   4.3280000e-01,
		   4.3150000e-01,
		   4.3120000e-01,
		   4.2990000e-01,
		   4.2910000e-01,
		   4.2680000e-01,
		   4.2590000e-01,
		   4.2390000e-01,
		   4.2140000e-01,
		   4.1920000e-01,
		   4.1980000e-01,
		   4.1760000e-01,
		   4.1580000e-01,
		   4.1510000e-01,
		   4.1480000e-01,
		   4.1290000e-01,
		   4.1080000e-01,
		   4.0980000e-01,
		   4.0940000e-01,
		   4.0490000e-01,
		   4.0370000e-01,
		   4.0270000e-01,
		   4.0100000e-01,
		   3.9720000e-01,
		   3.9740000e-01,
		   3.9360000e-01,
		   3.9060000e-01,
		   3.8940000e-01,
		   3.8870000e-01,
		   3.8350000e-01,
		   3.8180000e-01,
		   3.7850000e-01,
		   3.7630000e-01,
		   3.7250000e-01,
		   3.6960000e-01,
		   3.6680000e-01,
		   3.6430000e-01,
		   3.5970000e-01,
		   3.5760000e-01,
		   3.5450000e-01,
		   3.5040000e-01,
		   3.4370000e-01,
		   3.4170000e-01,
		   3.4040000e-01,
		   3.4000000e-01,
		   3.3540000e-01,
		   3.3230000e-01,
		   3.2850000e-01,
		   3.2710000e-01,
		   3.2410000e-01,
		   3.2200000e-01,
		   3.2090000e-01,
		   3.1840000e-01,
		   3.1450000e-01,
		   3.1410000e-01,
		   3.1140000e-01,
		   3.0730000e-01,
		   3.0620000e-01,
		   3.0520000e-01,
		   3.0260000e-01,
		   3.0020000e-01,
		   2.9760000e-01,
		   2.9670000e-01,
		   2.9580000e-01,
		   2.9540000e-01,
		   2.9380000e-01,
		   2.9240000e-01,
		   2.8960000e-01,
		   2.8910000e-01,
		   2.8700000e-01,
		   2.8580000e-01,
		   2.8390000e-01,
		   2.8380000e-01,
		   2.8140000e-01,
		   2.7970000e-01,
		   2.7800000e-01,
		   2.7680000e-01,
		   2.7470000e-01,
		   2.7360000e-01,
		   2.7240000e-01,
		   2.7180000e-01,
		   2.6990000e-01,
		   2.6970000e-01,
		   2.6650000e-01,
		   2.6420000e-01,
		   2.6370000e-01,
		   2.6480000e-01,
		   2.6280000e-01,
		   2.6040000e-01,
		   2.5970000e-01,
		   2.6040000e-01,
		   2.5930000e-01,
		   2.5740000e-01,
		   2.5760000e-01,
		   2.5680000e-01,
		   2.5600000e-01,
		   2.5720000e-01,
		   2.5530000e-01,
		   2.5250000e-01,
		   2.5340000e-01,
		   2.5390000e-01,
		   2.5520000e-01,
		   2.5390000e-01,
		   2.5360000e-01,
		   2.5390000e-01,
		   2.5500000e-01,
		   2.5330000e-01,
		   2.5260000e-01,
		   2.5270000e-01,
		   2.5170000e-01,
		   2.5190000e-01,
		   2.5290000e-01,
		   2.5270000e-01,
		   2.5290000e-01,
		   2.5450000e-01,
		   2.5370000e-01,
		   2.5200000e-01,
		   2.5320000e-01,
		   2.5480000e-01,
		   2.5470000e-01,
		   2.5500000e-01,
		   2.5500000e-01,
		   2.5490000e-01,
		   2.5380000e-01,
		   2.5450000e-01,
		   2.5640000e-01,
		   2.5640000e-01,
		   2.5710000e-01,
		   2.5890000e-01,
		   2.6030000e-01,
		   2.5880000e-01,
		   2.5770000e-01,
		   2.6040000e-01,
		   2.6100000e-01,
		   2.6010000e-01,
		   2.6250000e-01,
		   2.6360000e-01,
		   2.6280000e-01,
		   2.6380000e-01,
		   2.6490000e-01,
		   2.6520000e-01,
		   2.6640000e-01,
		   2.6720000e-01,
		   2.6720000e-01,
		   2.6580000e-01,
		   2.6500000e-01,
		   2.6620000e-01,
		   2.6640000e-01,
		   2.6490000e-01,
		   2.6390000e-01,
		   2.6380000e-01,
		   2.6360000e-01,
		   2.6390000e-01,
		   2.6510000e-01,
		   2.6430000e-01,
		   2.6250000e-01,
		   2.6340000e-01,
		   2.6370000e-01,
		   2.6240000e-01,
		   2.6200000e-01,
		   2.6170000e-01,
		   2.5960000e-01,
		   2.5790000e-01,
		   2.5760000e-01,
		   2.5790000e-01,
		   2.5620000e-01,
		   2.5440000e-01,
		   2.5420000e-01,
		   2.5220000e-01,
		   2.4880000e-01,
		   2.4730000e-01,
		   2.4870000e-01,
		   2.4980000e-01,
		   2.4860000e-01,
		   2.4880000e-01,
		   2.4690000e-01,
		   2.4640000e-01,
		   2.4600000e-01,
		   2.4390000e-01,
		   2.4250000e-01,
		   2.4270000e-01,
		   2.4230000e-01,
		   2.4060000e-01
	#	   2.3970000e-01,
	#	   2.3890000e-01,
	#	   2.3880000e-01,
	#	   2.3850000e-01,
	#	   2.3780000e-01,
	#	   2.3650000e-01,
	#	   2.3730000e-01,
	#	   2.3860000e-01,
	#	   2.3910000e-01,
	#	   2.3940000e-01,
	#	   2.3970000e-01,
	#	   2.3810000e-01,
	#	   2.3890000e-01,
	#	   2.4060000e-01,
	#	   2.4210000e-01,
	#	   2.4420000e-01,
	#	   2.4500000e-01,
	#	   2.4470000e-01,
	#	   2.4770000e-01,
	#	   2.4960000e-01,
	#	   2.5060000e-01,
	#	   2.5310000e-01,
	#	   2.5580000e-01,
	#	   2.5570000e-01,
	#	   2.5930000e-01,
	#	   2.6150000e-01,
	#	   2.6350000e-01,
	#	   2.6510000e-01,
	#	   2.6820000e-01,
	#	   2.7180000e-01,
	#	   2.7580000e-01,
	#	   2.7600000e-01,
	#	   2.7940000e-01,
	#	   2.8290000e-01,
	#	   2.8700000e-01,
	#	   2.9000000e-01,
	#	   2.9580000e-01,
	#	   2.9830000e-01,
	#	   3.0510000e-01,
	#	   3.1010000e-01,
	#	   3.1510000e-01,
	#	   3.1660000e-01,
	#	   3.2180000e-01,
	#	   3.2760000e-01,
	#	   3.3170000e-01,
	#	   3.3430000e-01,
	#	   3.4080000e-01,
	#	   3.4340000e-01,
	#	   3.4840000e-01,
	#	   3.5320000e-01,
	#	   3.6010000e-01,
	#	   3.6480000e-01,
	#	   3.6780000e-01,
	#	   3.7100000e-01,
	#	   3.7540000e-01,
	#	   3.7830000e-01,
	#	   3.7980000e-01,
	#	   3.8290000e-01,
	#	   3.8620000e-01,
	#	   3.8480000e-01,
	#	   3.8840000e-01,
	#	   3.9150000e-01,
	#	   3.9310000e-01,
	#	   3.9660000e-01,
	#	   3.9890000e-01,
	#	   3.9960000e-01,
	#	   4.0000000e-01,
	#	   4.0270000e-01,
	#	   4.0370000e-01,
	#	   4.0480000e-01,
	#	   4.0400000e-01,
	#	   4.0350000e-01,
	#	   4.0390000e-01,
	#	   4.0470000e-01,
	#	   4.0330000e-01,
	#	   4.0590000e-01,
	#	   4.0710000e-01,
	#	   4.0780000e-01,
	#	   4.0740000e-01,
	#	   4.0800000e-01,
	#	   4.0860000e-01,
	#	   4.0800000e-01,
	#	   4.1010000e-01,
	#	   4.1090000e-01,
	#	   4.0700000e-01,
	#	   4.0970000e-01,
	#	   4.0970000e-01,
	#	   4.0920000e-01,
	#	   4.1010000e-01,
	#	   4.0980000e-01,
	#	   4.0890000e-01,
	#	   4.0930000e-01,
	#	   4.0860000e-01,
	#	   4.0630000e-01,
	#	   4.0680000e-01,
	#	   4.0880000e-01,
	#	   4.1270000e-01,
	#	   4.1300000e-01,
	#	   4.1240000e-01,
	#	   4.1290000e-01
		])

		m75b8_1nm = np.array([
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   np.nan,
		   1.3290000e-01,
		   1.3800000e-01,
		   1.4380000e-01,
		   1.5450000e-01,
		   1.5850000e-01,
		   1.6880000e-01,
		   1.7640000e-01,
		   1.8640000e-01,
		   1.9440000e-01,
		   2.0970000e-01,
		   2.1980000e-01,
		   2.3720000e-01,
		   2.5050000e-01,
		   2.7020000e-01,
		   2.8200000e-01,
		   3.0510000e-01,
		   3.2180000e-01,
		   3.4540000e-01,
		   3.6080000e-01,
		   3.8680000e-01,
		   4.0560000e-01,
		   4.3230000e-01,
		   4.4890000e-01,
		   4.7400000e-01,
		   4.9280000e-01,
		   5.1640000e-01,
		   5.2970000e-01,
		   5.5040000e-01,
		   5.6070000e-01,
		   5.7590000e-01,
		   5.8710000e-01,
		   5.9940000e-01,
		   6.0420000e-01,
		   6.1250000e-01,
		   6.1620000e-01,
		   6.1970000e-01,
		   6.2230000e-01,
		   6.2650000e-01,
		   6.2700000e-01,
		   6.3120000e-01,
		   6.3430000e-01,
		   6.3700000e-01,
		   6.3870000e-01,
		   6.3770000e-01,
		   6.4150000e-01,
		   6.4390000e-01,
		   6.4110000e-01,
		   6.4250000e-01,
		   6.4570000e-01,
		   6.4470000e-01,
		   6.4420000e-01,
		   6.4720000e-01,
		   6.4680000e-01,
		   6.4690000e-01,
		   6.4830000e-01,
		   6.5140000e-01,
		   6.5210000e-01,
		   6.5510000e-01,
		   6.5680000e-01,
		   6.5670000e-01,
		   6.5630000e-01,
		   6.5610000e-01,
		   6.5740000e-01,
		   6.5860000e-01,
		   6.5730000e-01,
		   6.5640000e-01,
		   6.5950000e-01,
		   6.5660000e-01,
		   6.5690000e-01,
		   6.6050000e-01,
		   6.6130000e-01,
		   6.5810000e-01,
		   6.6040000e-01,
		   6.5980000e-01,
		   6.5810000e-01,
		   6.5900000e-01,
		   6.6130000e-01,
		   6.6090000e-01,
		   6.6110000e-01,
		   6.6090000e-01,
		   6.6010000e-01,
		   6.6070000e-01,
		   6.5960000e-01,
		   6.5940000e-01,
		   6.6070000e-01,
		   6.6020000e-01,
		   6.6050000e-01,
		   6.6080000e-01,
		   6.5970000e-01,
		   6.5920000e-01,
		   6.6160000e-01,
		   6.6080000e-01,
		   6.5830000e-01,
		   6.5910000e-01,
		   6.6090000e-01,
		   6.5850000e-01,
		   6.5790000e-01,
		   6.5830000e-01,
		   6.5650000e-01,
		   6.5730000e-01,
		   6.5750000e-01,
		   6.5500000e-01,
		   6.5410000e-01,
		   6.5380000e-01,
		   6.5280000e-01,
		   6.5510000e-01,
		   6.5300000e-01,
		   6.5180000e-01,
		   6.5410000e-01,
		   6.5270000e-01,
		   6.4990000e-01,
		   6.5160000e-01,
		   6.5140000e-01,
		   6.4840000e-01,
		   6.4920000e-01,
		   6.4890000e-01,
		   6.4680000e-01,
		   6.4820000e-01,
		   6.4930000e-01,
		   6.4630000e-01,
		   6.4550000e-01,
		   6.4540000e-01,
		   6.4570000e-01,
		   6.4530000e-01,
		   6.4410000e-01,
		   6.4230000e-01,
		   6.4300000e-01,
		   6.3990000e-01,
		   6.3870000e-01,
		   6.3880000e-01,
		   6.3660000e-01,
		   6.3320000e-01,
		   6.3420000e-01,
		   6.3430000e-01,
		   6.3180000e-01,
		   6.2930000e-01,
		   6.3050000e-01,
		   6.2620000e-01,
		   6.2410000e-01,
		   6.2460000e-01,
		   6.2350000e-01,
		   6.1940000e-01,
		   6.1620000e-01,
		   6.1150000e-01,
		   6.1000000e-01,
		   6.0800000e-01,
		   6.0530000e-01,
		   6.0390000e-01,
		   6.0180000e-01,
		   5.9750000e-01,
		   5.9800000e-01,
		   5.9470000e-01,
		   5.8950000e-01,
		   5.8920000e-01,
		   5.8890000e-01,
		   5.8390000e-01,
		   5.8110000e-01,
		   5.7690000e-01,
		   5.7430000e-01,
		   5.7130000e-01,
		   5.6870000e-01,
		   5.6430000e-01,
		   5.6260000e-01,
		   5.5940000e-01,
		   5.5750000e-01,
		   5.5470000e-01,
		   5.5150000e-01,
		   5.4480000e-01,
		   5.4450000e-01,
		   5.4090000e-01,
		   5.3640000e-01,
		   5.3010000e-01,
		   5.2840000e-01,
		   5.2610000e-01,
		   5.2250000e-01,
		   5.1530000e-01,
		   5.1010000e-01,
		   5.0600000e-01,
		   5.0380000e-01,
		   4.9790000e-01,
		   4.9440000e-01,
		   4.9280000e-01,
		   4.9090000e-01,
		   4.8690000e-01,
		   4.8730000e-01,
		   4.8160000e-01,
		   4.7750000e-01,
		   4.7550000e-01,
		   4.7230000e-01,
		   4.6800000e-01,
		   4.6580000e-01,
		   4.6410000e-01,
		   4.6090000e-01,
		   4.5780000e-01,
		   4.5580000e-01,
		   4.5060000e-01,
		   4.4760000e-01,
		   4.4580000e-01,
		   4.4370000e-01,
		   4.4190000e-01,
		   4.3910000e-01,
		   4.3760000e-01,
		   4.3680000e-01,
		   4.3360000e-01,
		   4.3140000e-01,
		   4.3070000e-01,
		   4.2760000e-01,
		   4.2390000e-01,
		   4.2560000e-01,
		   4.2420000e-01,
		   4.1960000e-01,
		   4.1890000e-01,
		   4.1800000e-01,
		   4.1520000e-01,
		   4.1450000e-01,
		   4.1120000e-01,
		   4.0980000e-01,
		   4.0760000e-01,
		   4.0650000e-01,
		   4.0500000e-01,
		   4.0400000e-01,
		   4.0260000e-01,
		   4.0080000e-01,
		   4.0010000e-01,
		   3.9830000e-01,
		   3.9720000e-01,
		   3.9680000e-01,
		   3.9300000e-01,
		   3.9120000e-01,
		   3.9070000e-01,
		   3.9140000e-01,
		   3.8830000e-01,
		   3.8840000e-01,
		   3.8910000e-01,
		   3.8750000e-01,
		   3.8640000e-01,
		   3.8670000e-01,
		   3.8660000e-01,
		   3.8550000e-01,
		   3.8550000e-01,
		   3.8590000e-01,
		   3.8530000e-01,
		   3.8190000e-01,
		   3.8170000e-01,
		   3.8320000e-01,
		   3.8190000e-01,
		   3.8010000e-01,
		   3.8110000e-01,
		   3.8100000e-01,
		   3.8100000e-01,
		   3.8210000e-01,
		   3.8050000e-01,
		   3.7850000e-01,
		   3.8100000e-01,
		   3.8080000e-01,
		   3.8010000e-01,
		   3.7990000e-01,
		   3.8030000e-01,
		   3.8170000e-01,
		   3.8180000e-01,
		   3.8040000e-01,
		   3.8050000e-01,
		   3.8200000e-01,
		   3.8120000e-01,
		   3.8300000e-01,
		   3.8550000e-01,
		   3.8520000e-01,
		   3.8480000e-01,
		   3.8660000e-01,
		   3.8690000e-01,
		   3.8750000e-01,
		   3.8780000e-01,
		   3.8880000e-01,
		   3.8940000e-01,
		   3.9030000e-01,
		   3.9210000e-01,
		   3.9360000e-01,
		   3.9410000e-01,
		   3.9430000e-01,
		   3.9480000e-01,
		   3.9400000e-01,
		   3.9460000e-01,
		   3.9540000e-01,
		   3.9560000e-01,
		   3.9460000e-01,
		   3.9380000e-01,
		   3.9550000e-01,
		   3.9390000e-01,
		   3.9260000e-01,
		   3.9450000e-01,
		   3.9420000e-01,
		   3.9310000e-01,
		   3.9340000e-01,
		   3.9260000e-01,
		   3.9310000e-01,
		   3.9150000e-01,
		   3.9160000e-01,
		   3.9010000e-01,
		   3.8890000e-01,
		   3.8680000e-01,
		   3.8860000e-01,
		   3.8810000e-01,
		   3.8570000e-01,
		   3.8440000e-01,
		   3.8350000e-01,
		   3.8160000e-01,
		   3.8070000e-01,
		   3.7830000e-01,
		   3.7590000e-01,
		   3.7550000e-01,
		   3.7890000e-01,
		   3.7780000e-01,
		   3.7560000e-01,
		   3.7420000e-01,
		   3.7220000e-01,
		   3.6960000e-01,
		   3.6910000e-01,
		   3.6610000e-01,
		   3.6410000e-01,
		   3.6340000e-01,
		   3.6350000e-01
	#	   3.6150000e-01,
	#	   3.6070000e-01,
	#	   3.5750000e-01,
	#	   3.5640000e-01,
	#	   3.5530000e-01,
	#	   3.5560000e-01,
	#	   3.5410000e-01,
	#	   3.5430000e-01,
	#	   3.5290000e-01,
	#	   3.5390000e-01,
	#	   3.5510000e-01,
	#	   3.5380000e-01,
	#	   3.5340000e-01,
	#	   3.5510000e-01,
	#	   3.5610000e-01,
	#	   3.5680000e-01,
	#	   3.5760000e-01,
	#	   3.5850000e-01,
	#	   3.5980000e-01,
	#	   3.6060000e-01,
	#	   3.6410000e-01,
	#	   3.6510000e-01,
	#	   3.6810000e-01,
	#	   3.6900000e-01,
	#	   3.7320000e-01,
	#	   3.7430000e-01,
	#	   3.7820000e-01,
	#	   3.8070000e-01,
	#	   3.8360000e-01,
	#	   3.8630000e-01,
	#	   3.8920000e-01,
	#	   3.9260000e-01,
	#	   3.9710000e-01,
	#	   3.9970000e-01,
	#	   4.0660000e-01,
	#	   4.1020000e-01,
	#	   4.1730000e-01,
	#	   4.2250000e-01,
	#	   4.3160000e-01,
	#	   4.3570000e-01,
	#	   4.4130000e-01,
	#	   4.4650000e-01,
	#	   4.5570000e-01,
	#	   4.6070000e-01,
	#	   4.6830000e-01,
	#	   4.7400000e-01,
	#	   4.8380000e-01,
	#	   4.8960000e-01,
	#	   5.0300000e-01,
	#	   5.0890000e-01,
	#	   5.1750000e-01,
	#	   5.2580000e-01,
	#	   5.3690000e-01,
	#	   5.4210000e-01,
	#	   5.5050000e-01,
	#	   5.5510000e-01,
	#	   5.6450000e-01,
	#	   5.6890000e-01,
	#	   5.7320000e-01,
	#	   5.7700000e-01,
	#	   5.8320000e-01,
	#	   5.8580000e-01,
	#	   5.8990000e-01,
	#	   5.9260000e-01,
	#	   5.9760000e-01,
	#	   6.0020000e-01,
	#	   6.0630000e-01,
	#	   6.0880000e-01,
	#	   6.1070000e-01,
	#	   6.1120000e-01,
	#	   6.1460000e-01,
	#	   6.1320000e-01,
	#	   6.1700000e-01,
	#	   6.2220000e-01,
	#	   6.1910000e-01,
	#	   6.1670000e-01,
	#	   6.2160000e-01,
	#	   6.2360000e-01,
	#	   6.1950000e-01,
	#	   6.2130000e-01,
	#	   6.2230000e-01,
	#	   6.2340000e-01,
	#	   6.2460000e-01,
	#	   6.2640000e-01,
	#	   6.2290000e-01,
	#	   6.2500000e-01,
	#	   6.2670000e-01,
	#	   6.2950000e-01,
	#	   6.3130000e-01,
	#	   6.2930000e-01,
	#	   6.3040000e-01,
	#	   6.3330000e-01,
	#	   6.3140000e-01,
	#	   6.2860000e-01,
	#	   6.2920000e-01,
	#	   6.3180000e-01,
	#	   6.3300000e-01,
	#	   6.3020000e-01,
	#	   6.2490000e-01,
	#	   6.2770000e-01
		])
	
	xvalues = np.empty(41)
	for i in range(0, 41):
		xvalues[i] = i*10 + 300
	xvalues_1nm = np.empty(401)
	for i in range(0, 401):
		xvalues_1nm[i] = i + 300
	
	if (args.mv2):
		m5b4 = m5b4_1nm
		m5b5 = m5b5_1nm
		m5b6 = m5b6_1nm
		m5b7 = m5b7_1nm
		m75b5 = m75b5_1nm
		m75b6 = m75b6_1nm
		m75b7 = m75b7_1nm
		m75b8 = m75b8_1nm
	if (args.mv3):
		m5b4 = np.interp(xvalues_1nm, xvalues, m5b4)
		m5b5 = np.interp(xvalues_1nm, xvalues, m5b5)
		m5b6 = np.interp(xvalues_1nm, xvalues, m5b6)
		m5b7 = np.interp(xvalues_1nm, xvalues, m5b7)
		m75b5 = np.interp(xvalues_1nm, xvalues, m75b5)
		m75b6 = np.interp(xvalues_1nm, xvalues, m75b6)
		m75b7 = np.interp(xvalues_1nm, xvalues, m75b7)
		m75b8 = np.interp(xvalues_1nm, xvalues, m75b8)
		
		for i in range(120): # 300-419
			m5b4[i] = m5b4_1nm[i]
			m5b5[i] = m5b5_1nm[i]
			m5b6[i] = m5b6_1nm[i]
			m5b7[i] = m5b7_1nm[i]
			m75b5[i] = m75b5_1nm[i]
			m75b6[i] = m75b6_1nm[i]
			m75b7[i] = m75b7_1nm[i]
			m75b8[i] = m75b8_1nm[i]
	
	# another alternate version: just remove the reflectance at 380 and 390 nm
	# 16/03/2025 -- fixed an error that was causing 380 to be left in: range was
	# set to (9,10), but 380 has index 8.
	if (args.mv4):
		for i in range(8,10): # 380-390
			m5b4[i] = np.nan
			m5b5[i] = np.nan
			m5b6[i] = np.nan
			m5b7[i] = np.nan
			m75b5[i] = np.nan
			m75b6[i] = np.nan
			m75b7[i] = np.nan
			m75b8[i] = np.nan
			m5gy5[i] = np.nan
			m5gy6[i] = np.nan
			m5gy7[i] = np.nan
			m5gy8[i] = np.nan
			m10yr5[i] = np.nan
			m10yr6[i] = np.nan
			m10yr7[i] = np.nan
			m10yr8[i] = np.nan
		
	# plot spectra
	if (args.mv2 or args.mv3):
		xvalues_1nm = np.empty(401)
		for i in range(0, 401):
			xvalues_1nm[i] = i + 300
		xvalues_alt = xvalues_1nm
	else: xvalues_alt = xvalues
	
	plt.subplot(2, 2, 1)
	plt.plot(xvalues_alt, m5b4, color='deepskyblue', label="5B4/6")
	plt.plot(xvalues_alt, m5b5, color='deepskyblue', label="5B5/6")
	plt.plot(xvalues_alt, m5b6, color='deepskyblue', label="5B6/6")
	plt.plot(xvalues_alt, m5b7, color='deepskyblue', label="5B7/6")
	#plt.legend()
	#plt.title("5B")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 2)
	plt.plot(xvalues_alt, m75b5, color='dodgerblue', label="7.5B5/4")
	plt.plot(xvalues_alt, m75b6, color='dodgerblue', label="7.5B6/4")
	plt.plot(xvalues_alt, m75b7, color='dodgerblue', label="7.5B7/4")
	plt.plot(xvalues_alt, m75b8, color='dodgerblue', label="7.5B8/4")
	#plt.legend()
	#plt.title("7.5B")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 3)
	plt.plot(xvalues, m10yr5, color='orange', label="10YR5/10")
	plt.plot(xvalues, m10yr6, color='orange', label="10YR6/10")
	plt.plot(xvalues, m10yr7, color='orange', label="10YR7/10")
	plt.plot(xvalues, m10yr8, color='orange', label="10YR8/10")
	#plt.legend()
	#plt.title("10YR")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 4)
	plt.plot(xvalues, m5gy5, color='chartreuse', label="5GY5/10")
	plt.plot(xvalues, m5gy6, color='chartreuse', label="5GY6/10")
	plt.plot(xvalues, m5gy7, color='chartreuse', label="5GY7/10")
	plt.plot(xvalues, m5gy8, color='chartreuse', label="5GY8/10")
	#plt.legend()
	#plt.title("5GY")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	plt.show()
	
	# plot linear interpolations (for checking)
	#plt.subplot(2, 2, 1)
	#plt.plot(x1nm, m5b4_1nm, marker='^', color='deepskyblue', mec='k', label="5B4/6")
	#plt.plot(x1nm, m5b5_1nm, marker='^', color='deepskyblue', mec='k', label="5B5/6")
	#plt.plot(x1nm, m5b6_1nm, marker='^', color='deepskyblue', mec='k', label="5B6/6")
	#plt.plot(x1nm, m5b7_1nm, marker='^', color='deepskyblue', mec='k', label="5B7/6")
	#plt.legend()
	#plt.title("5B")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	#plt.subplot(2, 2, 2)
	#plt.plot(x1nm, m75b5_1nm, marker='o', color='dodgerblue', mec='k', label="7.5B5/4")
	#plt.plot(x1nm, m75b6_1nm, marker='o', color='dodgerblue', mec='k', label="7.5B6/4")
#	plt.plot(x1nm, m75b7_1nm, marker='o', color='dodgerblue', mec='k', label="7.5B7/4")
#	plt.plot(x1nm, m75b8_1nm, marker='o', color='dodgerblue', mec='k', label="7.5B8/4")
#	#plt.legend()
#	plt.title("7.5B")
#	#plt.xlabel("Wavelength (nm)")
#	#plt.ylabel("Reflectance")
#	#plt.show()
#	plt.subplot(2, 2, 3)
#	plt.plot(x1nm, m10yr5_1nm, marker='s', color='orange', mec='k', label="10YR5/10")
#	plt.plot(x1nm, m10yr6_1nm, marker='s', color='orange', mec='k', label="10YR6/10")
#	plt.plot(x1nm, m10yr7_1nm, marker='s', color='orange', mec='k', label="10YR7/10")
#	plt.plot(x1nm, m10yr8_1nm, marker='s', color='orange', mec='k', label="10YR8/10")
#	#plt.legend()
#	plt.title("10YR")
#	#plt.xlabel("Wavelength (nm)")
#	#plt.ylabel("Reflectance")
#	#plt.show()
#	plt.subplot(2, 2, 4)
#	plt.plot(x1nm, m5gy5_1nm, marker='D', color='chartreuse', mec='k', label="5GY5/10")
#	plt.plot(x1nm, m5gy6_1nm, marker='D', color='chartreuse', mec='k', label="5GY6/10")
#	plt.plot(x1nm, m5gy7_1nm, marker='D', color='chartreuse', mec='k', label="5GY7/10")
#	plt.plot(x1nm, m5gy8_1nm, marker='D', color='chartreuse', mec='k', label="5GY8/10")
#	#plt.legend()
#	plt.title("5GY")
#	#plt.xlabel("Wavelength (nm)")
#	#plt.ylabel("Reflectance")
#	plt.show()
	
	# rendering
	m5bcs = np.empty((3, 4))
	print("5B4/6")
	m5bdata = spectral_rendering(m5b4)
	m5bcs[0][0] = m5bdata[0]
	m5bcs[1][0] = m5bdata[1]
	m5bcs[2][0] = m5bdata[2]
	print("5B5/6")
	m5bdata = spectral_rendering(m5b5)
	m5bcs[0][1] = m5bdata[0]
	m5bcs[1][1] = m5bdata[1]
	m5bcs[2][1] = m5bdata[2]
	print("5B6/6")
	m5bdata = spectral_rendering(m5b6)
	m5bcs[0][2] = m5bdata[0]
	m5bcs[1][2] = m5bdata[1]
	m5bcs[2][2] = m5bdata[2]
	print("5B7/6")
	m5bdata = spectral_rendering(m5b7)
	m5bcs[0][3] = m5bdata[0]
	m5bcs[1][3] = m5bdata[1]
	m5bcs[2][3] = m5bdata[2]
	
	m75bcs = np.empty((3, 4))
	print("7.5B5/4")
	m75bdata = spectral_rendering(m75b5)
	m75bcs[0][0] = m75bdata[0]
	m75bcs[1][0] = m75bdata[1]
	m75bcs[2][0] = m75bdata[2]
	print("7.5B6/4")
	m75bdata = spectral_rendering(m75b6)
	m75bcs[0][1] = m75bdata[0]
	m75bcs[1][1] = m75bdata[1]
	m75bcs[2][1] = m75bdata[2]
	print("7.5B7/4")
	m75bdata = spectral_rendering(m75b7)
	m75bcs[0][2] = m75bdata[0]
	m75bcs[1][2] = m75bdata[1]
	m75bcs[2][2] = m75bdata[2]
	print("7.5B8/4")
	m75bdata = spectral_rendering(m75b8)
	m75bcs[0][3] = m75bdata[0]
	m75bcs[1][3] = m75bdata[1]
	m75bcs[2][3] = m75bdata[2]
	
	m10yrcs = np.empty((3, 4))
	print("10YR5/10")
	m10yrdata = spectral_rendering(m10yr5)
	m10yrcs[0][0] = m10yrdata[0]
	m10yrcs[1][0] = m10yrdata[1]
	m10yrcs[2][0] = m10yrdata[2]
	print("10YR6/10")
	m10yrdata = spectral_rendering(m10yr6)
	m10yrcs[0][1] = m10yrdata[0]
	m10yrcs[1][1] = m10yrdata[1]
	m10yrcs[2][1] = m10yrdata[2]
	print("10YR7/10")
	m10yrdata = spectral_rendering(m10yr7)
	m10yrcs[0][2] = m10yrdata[0]
	m10yrcs[1][2] = m10yrdata[1]
	m10yrcs[2][2] = m10yrdata[2]
	print("10YR8/10")
	m10yrdata = spectral_rendering(m10yr8)
	m10yrcs[0][3] = m10yrdata[0]
	m10yrcs[1][3] = m10yrdata[1]
	m10yrcs[2][3] = m10yrdata[2]
	
	m5gycs = np.empty((3, 4))
	print("5GY5/10")
	m5gydata = spectral_rendering(m5gy5)
	m5gycs[0][0] = m5gydata[0]
	m5gycs[1][0] = m5gydata[1]
	m5gycs[2][0] = m5gydata[2]
	print("5GY6/10")
	m5gydata = spectral_rendering(m5gy6)
	m5gycs[0][1] = m5gydata[0]
	m5gycs[1][1] = m5gydata[1]
	m5gycs[2][1] = m5gydata[2]
	print("5GY7/10")
	m5gydata = spectral_rendering(m5gy7)
	m5gycs[0][2] = m5gydata[0]
	m5gycs[1][2] = m5gydata[1]
	m5gycs[2][2] = m5gydata[2]
	print("5GY8/10")
	m5gydata = spectral_rendering(m5gy8)
	m5gycs[0][3] = m5gydata[0]
	m5gycs[1][3] = m5gydata[1]
	m5gycs[2][3] = m5gydata[2]
	
	# color space plots
	if (l1 == m1):
		plt.plot(m5bcs[1], m5bcs[0], marker='^', linestyle='', color='deepskyblue', mec='k', label="5B")
		plt.plot(m75bcs[1], m75bcs[0], marker='o', linestyle='', color='dodgerblue', mec='k', label="7.5B")
		plt.plot(m10yrcs[1], m10yrcs[0], marker='s', linestyle='', color='orange', mec='k', label="10YR")
		plt.plot(m5gycs[1], m5gycs[0], marker='D', linestyle='', color='chartreuse', mec='k', label="5GY")
		plt.xlabel("Chromaticity ((L-S)/(L+S))")
		plt.ylabel("Brightness (L)")
		plt.legend()
		plt.show()
	else:
		plt.plot(m5bcs[1], m5bcs[2], marker='^', linestyle='', color='deepskyblue', mec='k', label="5B")
		plt.plot(m75bcs[1], m75bcs[2], marker='o', linestyle='', color='dodgerblue', mec='k', label="7.5B")
		plt.plot(m10yrcs[1], m10yrcs[2], marker='s', linestyle='', color='orange', mec='k', label="10YR")
		plt.plot(m5gycs[1], m5gycs[2], marker='D', linestyle='', color='chartreuse', mec='k', label="5GY")
		plt.xlabel("(L-S)/(L+M+S)")
		plt.ylabel("(M-0.5(L+S))/(L+M+S)")
		xborder = np.array([-math.sqrt(1/2), 0, math.sqrt(1/2), -math.sqrt(1/2)])
		yborder = np.array([-math.sqrt(2/3)/2, math.sqrt(2/3), -math.sqrt(2/3)/2, -math.sqrt(2/3)/2])
		plt.plot(xborder, yborder, '-k')
		plt.text(-math.sqrt(1/2) - 0.05, -math.sqrt(2/3)/2 - 0.025, 'S')
		plt.text(0 - 0.025, math.sqrt(2/3) + 0.0125, 'M')
		plt.text(math.sqrt(1/2) + 0.0125, -math.sqrt(2/3)/2 - 0.025, 'L')
		plt.legend()
		plt.show()
	
	# color contrast
	# Here we compare the expected performance with the study's definition of chance level,
	# 62.5%.
	# 11/02/2025 -- still not sure what definition of significance we should be using.
	# Note that using alternative=greater here will never produce a value less than 0.05
	# whether we compare against 62.5% (10 of 16) or 50% (8 of 16). If using 62.5%, the
	# minimum value (for no distinguishable pairs) is 0.2272491455078125, whereas if
	# using 50%, the minimum is 0.5981903076171875.
	# This has been fixed.
	box = np.empty((16, 4))
	labels = ["5B", "7.5B", "5GY", "10YR"]
	
	# 5B vs. 10YR
	print("5B vs. 10YR")
	color_disc(m5b4, m5b5, m5b6, m5b7, m10yr5, m10yr6, m10yr7, m10yr8, correct=41, trials=64, order=0, alternative='greater')
	
	# 7.5B vs. 10YR
	print("7.5B vs. 10YR")
	color_disc(m75b5, m75b6, m75b7, m75b8, m10yr5, m10yr6, m10yr7, m10yr8, correct=41, trials=64, order=1, alternative='greater')
	
	# 5GY vs. 10YR
	print("5GY vs. 10YR")
	color_disc(m5gy5, m5gy6, m5gy7, m5gy8, m10yr5, m10yr6, m10yr7, m10yr8, correct=41, trials=64, order=2, alternative='greater')
	
	# 10YR vs. 10YR, for comparison
	print("10YR vs. 10YR")
	color_disc(m10yr5, m10yr6, m10yr7, m10yr8, m10yr5, m10yr6, m10yr7, m10yr8, correct=41, trials=64, order=3, alternative='greater')
	
	# plot contrast values
	plt.boxplot(x=box, tick_labels=labels)
	plt.ylabel("ΔS (JND)")
	plt.show()
	
	# find optimal values for pigments
	# This heats up Myotis pretty badly while it's running, but it doesn't take long and is
	# vastly preferable to doing it by hand for both time and computation costs.
	if (args.mopt1):
		# optimize dichromacy for blue vs. orange with L held constant at specified value and
		# S varied between specified value and L-1
		print("optimizing 5B")
		median_5b = np.empty(194)
		for i in range(194):
			color_disc(m5b4, m5b5, m5b6, m5b7, m10yr5, m10yr6, m10yr7, m10yr8, order=0, sw=(args.sw + i), output=False)
			median_5b[i] = color_disc(m5b4, m5b5, m5b6, m5b7, m10yr5, m10yr6, m10yr7, m10yr8, order=0, sw=(args.sw + i), output=False)[0]
			print(str(args.sw + i) + " nm: " + str(median_5b[i]))
		
		print("")
		print("optimizing 7.5B")
		median_75b = np.empty(194)
		for i in range(194):
			median_75b[i] = color_disc(m75b5, m75b6, m75b7, m75b8, m10yr5, m10yr6, m10yr7, m10yr8, order=0, sw=args.sw + i, output=False)[0]
			print(str(args.sw + i) + " nm: " + str(median_75b[i]))
		
		xvalues = np.empty(194)
		for i in range(194):
			xvalues[i] = i + args.sw
		
		plt.plot(xvalues, median_5b, 'k', label="5B")
		plt.plot(xvalues, median_75b, '--k', label="7.5B")
		plt.xlabel("λmax of S cone (nm)")
		plt.ylabel("Median contrast (JND)")
		plt.legend()
		plt.show()
		
	if (args.mopt2):
		# optimize trichromacy for green and orange vs. orange
		print("optimizing 5GY")
		median_5gy = np.empty(194)
		min_5gy = np.empty(194)
		max_5gy = np.empty(194)
		for i in range(194):
			values = color_disc(m5gy5, m5gy6, m5gy7, m5gy8, m10yr5, m10yr6, m10yr7, m10yr8, order=0, mw=(args.sw + i), output=False)
			median_5gy[i] = values[0]
			min_5gy[i] = values[1]
			max_5gy[i] = values[2]
			print(str(args.sw + i) + " nm: " + str(values))
		
		print("")
		print("optimizing 10YR")
		median_10yr = np.empty(194)
		#min_10yr = np.empty(194) we don't need this, it's always 0 by definition
		max_10yr = np.empty(194)
		for i in range(194):
			values = color_disc(m10yr5, m10yr6, m10yr7, m10yr8, m10yr5, m10yr6, m10yr7, m10yr8, order=0, mw=(args.sw + i), output=False)
			median_10yr[i] = values[0]
			#min_10yr[i] = values[1]
			max_10yr[i] = values[2]
			print(str(args.sw + i) + " nm: " + str(values))
		
		xvalues = np.empty(194)
		for i in range(194):
			xvalues[i] = i + args.sw
		
		plt.plot(xvalues, median_5gy, 'k', label="5GY")
		plt.plot(xvalues, min_5gy, color='gray')
		plt.plot(xvalues, max_5gy, color='gray')
		plt.plot(xvalues, median_10yr, '--k', label="10YR")
		#plt.plot(xvalues, min_10yr, '--', color='gray')
		plt.plot(xvalues, max_10yr, '--', color='gray')
		plt.xlabel("λmax of M cone (nm)")
		plt.ylabel("Median contrast (JND)")
		plt.legend()
		plt.show()
		
	if (args.mopt3):
		# compare positions in color space on the LM axis
		# This seems to be less computationally intensive than finding the contrast
		# values.
		print("optimizing 5GY")
		lm_5gy5 = np.empty(193)
		lm_5gy6 = np.empty(193)
		lm_5gy7 = np.empty(193)
		lm_5gy8 = np.empty(193)
		for i in range(193):
			lm_5gy5[i] = spectral_rendering(m5gy5, output=False, mw=args.sw + i)
			lm_5gy6[i] = spectral_rendering(m5gy6, output=False, mw=args.sw + i)
			lm_5gy7[i] = spectral_rendering(m5gy7, output=False, mw=args.sw + i)
			lm_5gy8[i] = spectral_rendering(m5gy8, output=False, mw=args.sw + i)
			print(str(args.sw + i) + " nm: " + str(lm_5gy5[i]) + ", " + str(lm_5gy6[i]) + ", " + str(lm_5gy7[i]) + ", " + str(lm_5gy8[i]))
		
		print("")
		print("optimizing 10YR")
		lm_10yr5 = np.empty(193)
		lm_10yr6 = np.empty(193)
		lm_10yr7 = np.empty(193)
		lm_10yr8 = np.empty(193)
		for i in range(193):
			lm_10yr5[i] = spectral_rendering(m10yr5, output=False, mw=args.sw + i)
			lm_10yr6[i] = spectral_rendering(m10yr6, output=False, mw=args.sw + i)
			lm_10yr7[i] = spectral_rendering(m10yr7, output=False, mw=args.sw + i)
			lm_10yr8[i] = spectral_rendering(m10yr8, output=False, mw=args.sw + i)
			print(str(args.sw + i) + " nm: " + str(lm_10yr5[i]) + ", " + str(lm_10yr6[i]) + ", " + str(lm_10yr7[i]) + ", " + str(lm_10yr8[i]))
		
		xvalues = np.empty(193)
		for i in range(193):
			xvalues[i] = i + args.sw
		
		plt.plot(xvalues, lm_5gy5, 'k', label="5GY5/10")
		plt.plot(xvalues, lm_5gy6, '--k', label="5GY6/10")
		plt.plot(xvalues, lm_5gy7, ':k', label="5GY7/10")
		plt.plot(xvalues, lm_5gy8, '-.k', label="5GY8/10")
		plt.plot(xvalues, lm_10yr5, color='gray', label="10YR5/10")
		plt.plot(xvalues, lm_10yr6, '--', color='gray', label="10YR6/10")
		plt.plot(xvalues, lm_10yr7, ':', color='gray', label="10YR7/10")
		plt.plot(xvalues, lm_10yr8, '-.', color='gray', label="10YR8/10")
		plt.xlabel("λmax of M cone (nm)")
		plt.ylabel("L-M chromaticity")
		plt.legend()
		plt.show()

# simplified ERG example
if (args.erg):
	x1nm = np.empty(301)
	for i in range(301): x1nm[i] = i + 370
	x10nm = np.empty(31)
	for i in range(31): x10nm[i] = i*10 + 370

	y_example = np.empty(301)
	for i in range(301): y_example[i] = vpt(i + 370, 560)
	
	y_shifted = np.empty(31)
	newpeakx = 0
	newpeaky = 0
	for i in range(31):
		y_shifted[i] = vpt(i*10 + 370, 560) / mouse_filter_data[i + 6]
		#print(y_shifted[i])
		if (y_shifted[i] > newpeaky):
			newpeaky = y_shifted[i]
			newpeakx = i*10 + 370
	y_shifted = y_shifted / newpeaky # set peak to 1
	print("New peak (mouse): " + str(newpeakx))
	
	y_shifted1 = np.empty(31)
	newpeakx = 0
	newpeaky = 0
	for i in range(31):
		y_shifted1[i] = vpt(i*10 + 370, 560) / mouse_filter_data[i + 6]**2
		#print(y_shifted1[i])
		if (y_shifted1[i] > newpeaky):
			newpeaky = y_shifted1[i]
			newpeakx = i*10 + 370
	y_shifted1 = y_shifted1 / newpeaky # set peak to 1
	print("New peak (mouse^2): " + str(newpeakx))

	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.yscale("log")
	plt.plot(x1nm, y_example, "k")
	plt.plot(x10nm, y_shifted, '--k')
	plt.plot(x10nm, y_shifted1, ':k')
	plt.show()
	
	# function for scipy.optimize.curve_fit
	def vpt_fit2(xdata, t1, t2, scalet1, scalet2):
		scalet1 = abs(scalet1)
		scalet2 = abs(scalet2)
		if (scalet2 < 0):
			scalet2 = 0
		ydata = np.empty(xdata.shape[0])
		for i in range(xdata.shape[0]):
			value = scalet1*vpt(xdata[i], t1) + scalet2*vpt(xdata[i], t2)
			if (value >= 0):
				ydata[i] = value
			else:
				ydata[i] = 0
		return(ydata)
	
	# example best fits
	best_fit = scipy.optimize.curve_fit(vpt_fit2, x10nm, y_shifted, p0=[560, 360, 1, 1])
	total = best_fit[0][2] + best_fit[0][3]
	print("LWS: " + str(best_fit[0][0]) + " nm (" + str(best_fit[0][2]/total * 100) + "%)")
	print("SWS: " + str(best_fit[0][1]) + " nm (" + str(best_fit[0][3]/total * 100) + "%)")
	uvs = np.empty(301)
	for i in range(301): uvs[i] = vpt(i + 370, best_fit[0][1]) * best_fit[0][3]
	lws = np.empty(301)
	for i in range(301): lws[i] = vpt(i + 370, best_fit[0][0]) * best_fit[0][2]
	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.yscale("log")
	plt.ylim(0.01, 1.1)
	plt.plot(x10nm, y_shifted, 'k')
	plt.plot(x1nm, uvs, ':k')
	plt.plot(x1nm, lws, ':k')
	plt.show()
	
	best_fit = scipy.optimize.curve_fit(vpt_fit2, x10nm, y_shifted1, p0=[560, 360, 1, 1])
	total = best_fit[0][2] + best_fit[0][3]
	print("LWS: " + str(best_fit[0][0]) + " nm (" + str(best_fit[0][2]/total * 100) + "%)")
	print("SWS: " + str(best_fit[0][1]) + " nm (" + str(best_fit[0][3]/total * 100) + "%)")
	uvs = np.empty(301)
	for i in range(301): uvs[i] = vpt(i + 370, best_fit[0][1]) * best_fit[0][3]
	lws = np.empty(301)
	for i in range(301): lws[i] = vpt(i + 370, best_fit[0][0]) * best_fit[0][2]
	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.yscale("log")
	plt.ylim(0.01, 1.1)
	plt.plot(x10nm, y_shifted, 'k')
	plt.plot(x1nm, uvs, ':k')
	plt.plot(x1nm, lws, ':k')
	plt.show()

# print execution time
print("%s seconds" % (time.time() - start_time))
