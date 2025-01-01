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
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--lw", type=int, default=0, help="longwave cone sensitivity")
parser.add_argument("-m", "--mw", type=int, default=0, help="mediumwave cone sensitivity")
parser.add_argument("-s", "--sw", type=int, default=0, help="shortwave cone sensitivity")
parser.add_argument("--rod", type=int, default=0, help="rod sensitivity")
parser.add_argument("--weber", type=float, default=0.018, help="Weber fraction for L cones")
parser.add_argument("--weberm", type=float, default=0, help="override Weber fraction for M cones")
parser.add_argument("--webers", type=float, default=0, help="override Weber fraction for S cones")
parser.add_argument("--weberb", type=float, default=0.14, help="Weber fraction for brightness")
parser.add_argument("--lp", type=float, default=9, help="number of L cones per receptive field")
parser.add_argument("--mp", type=float, default=4.5, help="number of M cones per receptive field")
parser.add_argument("--sp", type=float, default=1, help="number of S cones per receptive field")
parser.add_argument("--filter", type=str, default="none", help="type of lens filtering")
parser.add_argument("--qn", help="use quantum noise in color differences", action="store_true")
parser.add_argument("--lb", type=float, default=0.68990272, help="contribution of L cones to perception of brightness")
parser.add_argument("--mb", type=float, default=0.34832189, help="contribution of M cones to perception of brightness")
parser.add_argument("--sb", type=float, default=0.0, help="contribution of S cones to perception of brightness")
parser.add_argument("--rb", type=float, default=0.0, help="contribution of rods to perception of brightness")
parser.add_argument("--novk", help="disable von Kries transform", action="store_true")
parser.add_argument("--peak", help="show wavelength with maximum sensitivity", action="store_true")
parser.add_argument("--primaries", help="show primary and secondary wavelengths", action="store_true")
parser.add_argument("--white", default="e", help="reference white")
parser.add_argument("--lighting", help="standard light sources", action="store_true")
parser.add_argument("--leaves", help="several types of leaves", action="store_true")
parser.add_argument("--flowers", help="several types of flowers", action="store_true")
parser.add_argument("--sky", help="daylight and a clear blue sky", action="store_true")
parser.add_argument("--kodak", help="Kodak Wratten camera filters", action="store_true")
parser.add_argument("--kcheck", help="check camera filter brightness matches", action="store_true")
parser.add_argument("--check", help="show colormath spectral rendering for comparison", action="store_true")
parser.add_argument("--qcheck", help="show absolute quantum catches", action="store_true")
parser.add_argument("--vpt", help="plot visual pigment templates", action="store_true")
parser.add_argument("--triangle", help="plot color triangle", action="store_true")
parser.add_argument("--triangle1", help="plot logarithmic noise-scaled color triangle", action="store_true")
parser.add_argument("--hexagon", help="plot color hexagon", action="store_true")
parser.add_argument("--wd", help="plot wavelength discrimination function", action="store_true")
parser.add_argument("--blackbody", type=int, default=0, help="plot black body curve for a given temperature in Kelvin")
parser.add_argument("--munsell", help="Munsell color cards", action="store_true")
parser.add_argument("--convergence", type=float, default=1, help="cone-to-ganglion cell convergence ratio")
parser.add_argument("--ff", type=float, default=22.4, help="flicker fusion frequency (Hz)")
parser.add_argument("--osd", type=float, default=1.2, help="outer segment diameter (um)")
parser.add_argument("--fl", type=float, default=6.81, help="focal length (mm)")
parser.add_argument("--pupil", type=float, default=6, help="pupil diameter (mm)")
parser.add_argument("--od", type=float, default=0.015, help="optical density of pigment (um)")
parser.add_argument("--osl", type=float, default=30, help="outer segment length (um)")
parser.add_argument("--dist", type=float, default=5, help="distance from light source (cm)")
parser.add_argument("--width", type=float, default=3.8, help="width of light/object (cm)")
parser.add_argument("--height", type=float, default=0, help="height of light/object (cm)")
parser.add_argument("--radius", type=float, default=1, help="radius of light source (cm)")
parser.add_argument("--watts", type=int, default=12, help="number of watts produced by light source")
parser.add_argument("--lux", type=int, default=115, help="illuminance in lux")
#parser.add_argument("--start", type=int, default=300, help="start of wavelength interval")
#parser.add_argument("--end", type=int, default=900, help="end of wavelength interval")
parser.add_argument("--ct", type=int, default=2856, help="color temperature of incandescent/blackbody illuminant")
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
# and 68-73%, 20-25% and 7% in the fat-tailed dunnart ("Diversity of Color Vision [sic]: Not All
# Australian Marsupials are Trichromatic"). In Norway rats, the percentage of S cones is 11-12%
# ("Cone-based vision of rats for ultraviolet and visible lights") and in the house mouse the
# percentage of exclusive non-coexpressing S cones is 26% ("Number and Distribution of Mouse
# Retinal Cone Photoreceptors").

# function approximating human lens filtering, for real this time (Lamb 1995)
def template_filter(w, a, b, c, d, e, f):
	try:
		density = a*math.exp((b - w) / c) + d*math.exp((e - w) / f)
	except OverflowError:
		print("Warning (template_filter): math overflow, returning 0")
		return 0
	if (density < 0):
		return 0
	return math.exp(-density)
def human_filter(w):
	return template_filter(w, 1.1, 400, 15, 0.11, 500, 80)

# array version for curve fitting
def filter_fit(xdata, a, b, c, d, e, f):
	ydata = np.empty(xdata.shape[0])
	for i in range(xdata.shape[0]):
		ydata[i] = template_filter(xdata[i], a, b, c, d, e, f)
	return(ydata)

# Data for Thylamys elegans from Palacios et al. (2010). They say the measurements were
# every 20 nm, but the graph shows every 4 nm. Numbered because I can't count.
opossum_filter_data = np.array([
	4.198057132, # 300
	4.940696766, # 304
	5.078876161, # 308
	5.501061962, # 312
	6.584859955, # 316
	7.209278612, # 320
	7.901294748, # 324
	8.620954224, # 328
	9.398772899, # 332
	10.347077491, # 336
	10.347077491, # 340
	11.277550076, # 344
	11.7416299, # 348
	12.36134807, # 352
	12.706348891, # 356
	13.514198455, # 360
	13.671217103, # 364
	13.883186681, # 368
	14.331523604, # 372
	14.679061196, # 376
	14.956837592, # 380
	15.445837458, # 384
	15.658142785, # 388
	16.148261815, # 392
	16.148261815, # 396
	16.469760201, # 400
	17.049486926, # 404
	17.049486926, # 408
	17.049486926, # 412
	17.754783802, # 416
	18.383417975, # 420
	18.665186055, # 424
	18.465415359, # 428
	18.800194486, # 432
	18.800194486, # 436
	19.319971349, # 440
	19.747790274, # 444
	19.747790274, # 448
	19.747790274, # 452
	19.554921087, # 456
	19.96248315, # 460
	19.627330969, # 464
	19.76767408, # 468
	19.860751183, # 472
	19.921297931, # 476
	20.081897902, # 480
	21.329242998, # 484
	21.266830976, # 488
	20.999089747, # 492
	21.050795104, # 496
	20.877623195, # 500
	20.972341738, # 504
	21.270598827, # 508
	21.485291704, # 512
	21.705468482, # 516
	21.848684112, # 520
	21.723113961, # 526
	22.91506046, # 530
	22.780275861, # 534
	22.620944276, # 536
	23.162768661, # 540
	23.162768661, # 544
	22.572782937, # 548
	22.581960078, # 552
	22.620571221, # 556
	22.597516452, # 560
	23.17746701, # 564
	23.074242824, # 568
	23.545671821, # 572
	23.609016479, # 576
	23.112518217, # 580
	23.144973961, # 584
	23.489564422, # 588
	23.638674313, # 592
	24.366951349, # 596
	23.92103928, # 600
	24.628835622, # 604
	24.580935422, # 608
	24.381015504, # 612
	23.953420413, # 616
	24.122488721, # 620
	24.43510841, # 624
	24.616338296, # 628
	24.61712171, # 632
	24.756233741, # 636
	24.676549295, # 640
	24.652897639, # 644
	24.798276986, # 648
	25.008567819, # 652
	25.276943241, # 656
	25.865996329, # 660
	25.24340564, # 664
	25.300445676, # 668
	25.404042916, # 672
	25.44354939, # 676
	25.638134628, # 680
	25.358679486, # 684
	25.557368324, # 688
	25.510661898, # 692
	26.081547231, # 696
	26.750956263 # 700
])

# normalize to 1 at 700
opossum_filter_data = opossum_filter_data / opossum_filter_data[100]

xvalues = np.empty(101)
for i in range(0, 101):
	xvalues[i] = i*4 + 300

opossum_fit = scipy.optimize.curve_fit(filter_fit, xvalues, opossum_filter_data, p0=[1.1, 400, 15, 0.11, 500, 80])

# mouse (Jacobs & Williams 2007)
mouse_filter_data = np.array([
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

xvalues = np.empty(40)
for i in range(0, 40):
	xvalues[i] = i*10 + 310

mouse_fit = scipy.optimize.curve_fit(filter_fit, xvalues, mouse_filter_data/100, p0=[1.1, 400, 15, 0.11, 500, 80])

mousesquared_fit = scipy.optimize.curve_fit(filter_fit, xvalues, (mouse_filter_data/100)**2, p0=[1.1, 400, 15, 0.11, 500, 80])

# functions
def opossum_filter(w):
	return template_filter(w, *opossum_fit[0])

def mouse_filter(w):
	return template_filter(w, *mouse_fit[0])

def mousesquared_filter(w):
	return template_filter(w, *mousesquared_fit[0])

def lens_filter(w):
	if (args.filter == "human"):
		return human_filter(w)
	elif (args.filter == "opossum"):
		return opossum_filter(w)
	elif (args.filter == "mouse"):
		return mouse_filter(w)
	elif (args.filter == "mousesquared"):
		return mousesquared_filter(w)
	return 1

def sensitivity(w, l=l1, m=m1, s=s1):
	if (args.filter == "human"):
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
		print("Warning (blackbody): math overflow, returning 1")
		return 1
	return value

# normal distribution
def normal(mu, std_dev, x):
	return (1 / std_dev*math.sqrt(2*math.pi)) * math.exp((-1/2)*((x - mu)/std_dev)**2)

# light sources -- standardize to 300-900
#d65 = np.zeros(61)
#for i in range(4, 54):
#	d65[i] = spectral_constants.REF_ILLUM_TABLE["d65"][i-4]

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

d65 = np.zeros(51)
for i in range(0, 51):
	d65[i] = d65_1nm[i*10]
	#print(d65[i])

e = np.empty(51)
for i in range(0, 51):
	e[i] = 100

a = np.zeros(51)
for i in range(4, 51):
	a[i] = spectral_constants.REF_ILLUM_TABLE["a"][i-4]

# incandescent lighting with specified color temperature
incandescent = np.empty(51)
for i in range(0, 51):
	w = i*10 + 300
	# normalize to 100% at 560 nm
	incandescent[i] = 100*blackbody(w, args.ct) / blackbody(560, args.ct)
# 1-nm resolution
incandescent_1nm = np.empty(501)
for i in range(501):
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

# find sensitivity of a given photoreceptor type to a wavelength using visual pigment templates
# from Govardovskii et al. (2000)
def vpt(w, lmax):
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
		print("Warning (vpt): math overflow, clipping to 2.2250738585072014e-308")
		return 2.2250738585072014e-308
	
	return alpha + beta

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

# estimate hue, saturation and lightness for a spectral power distribution
# Lens filtering is now used for colors properly with a von Kries transform.
def spectral_rendering(table, light_source=wp):
	table_l = 0
	table_m = 0
	table_s = 0
	brightness = 0
	
	for i in range(0, table.shape[0]):
		w = i*10 + 300 # wavelength
		if (table[i] > 0): # don't add zeroes or nans
			table_l += vpt(w, l1) * table[i] * light_source[i] * lens_filter(w)
			# remove either M or S for dichromacy
			if (m1 != l1):
				table_m += vpt(w, m1) * table[i] * light_source[i] * lens_filter(w)
			if (s1 != m1):
				table_s += vpt(w, s1) * table[i] * light_source[i] * lens_filter(w)
			
			# brightness
			brightness += sensitivity(w) * table[i] * light_source[i]
	
	# normalize according to provided light source
	n = 0
	wpl = 0
	wpm = 0
	wps = 0
	for i in range(0, table.shape[0]):
		w = i*10 + 300
		n += sensitivity(w) * light_source[i]
		wpl += vpt(w, l1) * light_source[i] * lens_filter(w)
		wpm += vpt(w, m1) * light_source[i] * lens_filter(w)
		wps += vpt(w, s1) * light_source[i] * lens_filter(w)
	
	# von Kries transform
	if (args.novk == False):
		table_l = table_l / wpl
		table_m = table_m / wpm
		table_s = table_s / wps
	
	# express brightness of reflected/emitted light proportional to the light source
	brightness = brightness / n
	
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
	if (l1 != m1 != s1):
		print("Cone response: l=" + str(table_l) + ", m=" + str(table_m) + ", s=" + str(table_s))
		
		# This should be the nearest wavelength on a line from the white point.
		match = 300
		x0 = math.sqrt(1/2)*(table_l - table_s) / (table_l + table_m + table_s)
		y0 = math.sqrt(2/3)*(table_m - 0.5*(table_l + table_s)) / (table_l + table_m + table_s)
		wx = math.sqrt(1/2)*(wpl - wps) / (wpl + wpm + wps)
		wy = math.sqrt(2/3)*(wpm - 0.5*(wpl + wps)) / (wpl + wpm + wps)
		
		# adjust to chosen visible range
		while (lens_filter(match) == 0):
			match += 1
		
		for i in range(300, 800):
			matchl1 = vpt(match, l1)*lens_filter(match) / wpl
			matchm1 = vpt(match, m1)*lens_filter(match) / wpm
			matchs1 = vpt(match, s1)*lens_filter(match) / wps
			matchl2 = vpt(i, l1)*lens_filter(i) / wpl
			matchm2 = vpt(i, m1)*lens_filter(i) / wpm
			matchs2 = vpt(i, s1)*lens_filter(i) / wps
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
		
		lms = np.array([
			[table_l],
			[table_m],
			[table_s]
		])
		
		if (args.check):
			matrix = np.matmul(lms_to_xyz, lms)
			xyz = XYZColor(*matrix, illuminant="e")
			print("Color coordinates:")
			print("CIE XYZ: " + str(xyz.get_value_tuple()))
			print("sRGB: " + str(convert_color(xyz, sRGBColor).get_value_tuple()))
			print("CIE LAB: " + str(convert_color(xyz, LabColor).get_value_tuple()))
	elif (l1 == m1 != s1):
		# typical dichromacy
		print("Cone response: l/m=" + str(table_l) + ", s=" + str(table_s))
		
		match = 300
		for i in range(300, 800):
			diff0 = abs(table_l / table_s - vpt(match, l1) / vpt(match, s1))
			diff1 = abs(table_l / table_s - vpt(i, l1) / vpt(i, s1))
			if (diff1 < diff0):
				match = i
		print("Matching wavelength: " + str(match))
		
		pos = (table_l - table_s) / (table_l + table_s)
		dist = pos - (wpl - wps) / (wpl + wps)
		print("Position on color line: " + str(pos))
	elif (l1 != m1 == s1):
		# tritanopia
		print("Cone response: l=" + str(table_l) + ", m=" + str(table_m))
	elif (l1 == m1 == s1):
		# monochromacy
		print("Cone response: " + str(table_l))
	
	# estimate brightness
	print("Brightness relative to light source: " + str(brightness))
	
	# "SpectralColor" conversion for comparison to check if our estimates are reasonable
	if (args.check):
		spectral = SpectralColor(spec_340nm=table[4],
		spec_350nm=table[4+1],
		spec_360nm=table[4+2],
		spec_370nm=table[4+3],
		spec_380nm=table[4+4],
		spec_390nm=table[4+5],
		spec_400nm=table[4+6],
		spec_410nm=table[4+7],
		spec_420nm=table[4+8],
		spec_430nm=table[4+9],
		spec_440nm=table[4+10],
		spec_450nm=table[4+11],
		spec_460nm=table[4+12],
		spec_470nm=table[4+13],
		spec_480nm=table[4+14],
		spec_490nm=table[4+15],
		spec_500nm=table[4+16],
		spec_510nm=table[4+17],
		spec_520nm=table[4+18],
		spec_530nm=table[4+19],
		spec_540nm=table[4+20],
		spec_550nm=table[4+21],
		spec_560nm=table[4+22],
		spec_570nm=table[4+23],
		spec_580nm=table[4+24],
		spec_590nm=table[4+25],
		spec_600nm=table[4+26],
		spec_610nm=table[4+27],
		spec_620nm=table[4+28],
		spec_630nm=table[4+29],
		spec_640nm=table[4+30],
		spec_650nm=table[4+31],
		spec_660nm=table[4+32],
		spec_670nm=table[4+33],
		spec_680nm=table[4+34],
		spec_690nm=table[4+35],
		spec_700nm=table[4+36],
		spec_710nm=table[4+37],
		spec_720nm=table[4+38],
		spec_730nm=table[4+39],
		spec_740nm=table[4+40],
		spec_750nm=table[4+41],
		spec_760nm=table[4+42],
		spec_770nm=table[4+43],
		spec_780nm=table[4+44],
		spec_790nm=table[4+45],
		spec_800nm=table[4+46],
		spec_810nm=table[4+47],
		spec_820nm=table[4+48],
		spec_830nm=table[4+49], illuminant='e')
		
		print("colormath conversion:")
		print(convert_color(spectral, XYZColor))
		print(convert_color(spectral, sRGBColor))
		print(convert_color(spectral, LabColor))
	
	# break up text
	print("")
	
	# return brightness and color for comparisons
	if (l1 != m1):
		return(brightness, posx, posy)
	elif (m1 != s1): # add zeroes for compatibility
		return(brightness, pos, 0)
	else:
		return(brightness, 0, 0)

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
# with linear interpolation. This is probably better than hacking together a vague
# estimate of the "true size" of the integral.
def color_contrast(table1, table2, quantum_noise=args.qn):
	# interpolate 1-nm intervals
	x10nm = np.empty(51)
	for i in range(51):
		x10nm[i] = i*10 + 300
	
	x1nm = np.empty(501)
	for i in range(501):
		x1nm[i] = i + 300
	
	table1 = np.interp(x1nm, x10nm, table1)
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
			ql1 += vpt(w, l1) * table1[i] * wp_1nm[i] * lens_filter(w)
			qm1 += vpt(w, m1) * table1[i] * wp_1nm[i] * lens_filter(w)
			qs1 += vpt(w, s1) * table1[i] * wp_1nm[i] * lens_filter(w)
		wpl += vpt(w, l1) * wp_1nm[i] * lens_filter(w)
		wpm += vpt(w, m1) * wp_1nm[i] * lens_filter(w)
		wps += vpt(w, s1) * wp_1nm[i] * lens_filter(w)
	for i in range(0, table2.shape[0]):
		w = i + 300
		if (table2[i] > 0):
			ql2 += vpt(w, l1) * table2[i] * wp_1nm[i] * lens_filter(w)
			qm2 += vpt(w, m1) * table2[i] * wp_1nm[i] * lens_filter(w)
			qs2 += vpt(w, s1) * table2[i] * wp_1nm[i] * lens_filter(w)
	
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
		T = 1 / lens_filter(800) # set to 1 at 800 nm
		
		# select illuminant
		if (args.white == "i"):
			photons = ia
		elif (args.white == "d65"):
			photons = d65a
	
		aql1 = 0
		aqm1 = 0
		aqs1 = 0
		aql2 = 0
		aqm2 = 0
		aqs2 = 0
		for i in range(table1.shape[0]):
			w = i + 300
			if (table1[i] > 0):
				aql1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, l1)*l)) * vpt(w, l1) * table1[i] * photons[i] * lens_filter(w)
				aqm1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, m1)*l)) * vpt(w, m1) * table1[i] * photons[i] * lens_filter(w) 
				aqs1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, s1)*l)) * vpt(w, s1) * table1[i] * photons[i] * lens_filter(w)
		for i in range(table2.shape[0]):
			w = i + 300
			if (table2[i] > 0):
				aql2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, l1)*l)) * vpt(w, l1) * table2[i] * photons[i] * lens_filter(w)
				aqm2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, m1)*l)) * vpt(w, m1) * table2[i] * photons[i] * lens_filter(w)
				aqs2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, s1)*l)) * vpt(w, s1) * table2[i] * photons[i] * lens_filter(w)
		
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
		w = i*10 + 300
		wpt += sensitivity(w, l1) * wp[i]
		#ql1 += vpt(w, l1) * table1[i] * wp[i] * lens_filter(w)
		#qm1 += vpt(w, m1) * table1[i] * wp[i] * lens_filter(w)
		#qr1 += vpt(w, rod) * table1[i] * wp[i] * lens_filter(w)
		q1 += sensitivity(w) * table1[i] * wp[i]
	for i in range(table2.shape[0]):
		w = i*10 + 300
		#ql2 += vpt(w, l1) * table2[i] * wp[i] * lens_filter(w)
		#qm2 += vpt(w, m1) * table2[i] * wp[i] * lens_filter(w)
		#qr2 += vpt(w, rod) * table2[i] * wp[i] * lens_filter(w)
		q2 += sensitivity(w) * table2[i] * wp[i]
	
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

	plt.plot(xvalues, lvalues/lmax, 'gray', label="L (" + str(args.lw) + ")")
	# don't plot redundant curves
	if (l1 != m1):
		plt.plot(xvalues, mvalues/mmax, 'silver', label="M (" + str(args.mw) + ")")
	if (m1 != s1):
		plt.plot(xvalues, svalues/smax, 'k', label="S (" + str(args.sw) + ")")
	if (args.rod != 0):
		plt.plot(xvalues, rvalues/rmax, ':k', label="rod (" + str(args.rod) + ")")
	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.legend()
	plt.show()

# show luminous efficiency function and lens filtering
if (args.peak):
	xvalues = np.empty(500)
	yvalues = np.empty(500)
	yvalues1 = np.empty(500)
	ms = 300
	for i in range(300, 800):
		if (sensitivity(i) > (sensitivity(ms))):
			ms = i
		xvalues[i-300] = i
		yvalues[i-300] = sensitivity(i)
		yvalues1[i-300] = lens_filter(i)
	print("Maximum sensitivity: " + str(ms))
	print("")
	
	plt.plot(xvalues, yvalues/sensitivity(ms), 'k')
	if (args.filter == "none"):
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
# Now a real Maxwell triangle with Cartesian coordinates.
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
	xvalues = np.empty(37)
	yvalues = np.empty(37)
	labels = np.empty(37)
	
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
	
	for i in range(0, 37):
		w = i*10 + 340
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
	for i in range(0, 37):
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
	
	xvalues = np.empty(350)
	yvalues = np.empty(350)
	for i in range(350, 700):
		xvalues[i-350] = i
		
		# derivatives -- pretty sure the lens filtering cancels
		#scale1 = (vpt(i+1, l1) + vpt(i+1, m1) + vpt(i+1, s1))
		scale1 = 1
		dl1 = math.log(1 + kl*vpt(i+1, l1)) / scale1
		dm1 = math.log(1 + km*vpt(i+1, m1)) / scale1
		ds1 = math.log(1 + ks*vpt(i+1, s1)) / scale1
		#scale2 = (vpt(i-1, l1) + vpt(i-1, m1) + vpt(i-1, s1))
		scale2 = 1
		dl2 = math.log(1 + kl*vpt(i-1, l1)) / scale2
		dm2 = math.log(1 + km*vpt(i-1, m1)) / scale2
		ds2 = math.log(1 + ks*vpt(i-1, s1)) / scale2
		# Not totally sure if this is supposed to be logarithmic or linear. The
		# contrast functions are logarithmic.
		dll = (dl1 - dl2) / 2
		dlm = (dm1 - dm2) / 2
		dls = (ds1 - ds2) / 2
		dfl = (kl / (1 + kl*vpt(i, l1))) * dll
		dfm = (km / (1 + km*vpt(i, m1))) * dlm
		dfs = (ks / (1 + ks*vpt(i, s1))) * dls
		
		v = math.sqrt(((wl*wm)**2 + (wl*ws)**2 + (wm*ws)**2) / (wl**2*(dfm - dfs)**2 + wm**2*(dfl - dfs)**2 + ws**2*(dfm - dfl)**2))
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
		
		if (v > 100):
			yvalues[i-350] = float('inf')
		else:
			yvalues[i-350] = v
	
	plt.plot(xvalues, yvalues, 'k')
	plt.xlabel("Wavelength (nm)")
	plt.ylabel('Δλ (nm)')
	plt.show()
	
# white
# The "light source" is specified as E because these are not reflective objects
# and multiplying them by some other SPD does not make sense. Same with the sky
# and daylight spectra below.
	
if (args.lighting):
	print("White (E, equal-energy)")
	spectral_rendering(e/100, light_source=e)
	
	print("White (D65)")
	spectral_rendering(d65/100, light_source=e)
	
	print("Incandescent lighting (A)")
	spectral_rendering(a/100, light_source=e)
	
	print("Incandescent lighting (approximated with 2856 K blackbody spectrum)")
	spectral_rendering(incandescent/100, light_source=e)
	
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
		0.0, # 300
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.01, # 400
		0.01,
		0.02,
		0.02,
		0.02,
		0.03,
		0.03,
		0.03,
		0.04,
		0.04,
		0.05, # 500
		0.08,
		0.19,
		0.23,
		0.25,
		0.27,
		0.25,
		0.23,
		0.19,
		0.17,
		0.17, # 600
		0.15,
		0.14,
		0.14,
		0.10,
		0.07,
		0.05,
		0.04,
		0.10,
		0.30,
		0.35, # 700
		0.45,
		0.55,
		0.60,
		0.63,
		0.65,
		0.65,
		0.65,
		0.65,
		0.65,
		0.65, # 800
		0.65,
		0.65,
		0.65,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0 # 900
	])
	print("Green maple leaf")
	spectral_rendering(leaf_table)
	
	red_leaf_table = np.array([
		0.0, # 300
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0, # 400
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0, # 500
		0.01,
		0.01,
		0.02,
		0.02,
		0.03,
		0.05,
		0.06,
		0.08,
		0.10,
		0.10, # 600
		0.09,
		0.08,
		0.05,
		0.10,
		0.30,
		0.35,
		0.45,
		0.55,
		0.60,
		0.63, # 700
		0.65,
		0.65,
		0.65,
		0.65,
		0.65,
		0.65,
		0.65,
		0.65,
		0.65,
		0.0, # 800
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0 # 900
	])
	print("Red maple leaf")
	spectral_rendering(red_leaf_table)
	
	# color/brightness contrast test
	print("Color contrast between green and red maple leaves: " + str(color_contrast(leaf_table, red_leaf_table)))
	print("Brightness contrast between green and red maple leaves: " + str(brightness_contrast(leaf_table, red_leaf_table)))
	print("")
	
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
		0.80,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0
	])
	print("Corn leaf")
	spectral_rendering(corn_table)
		
	# hydrangea
	# estimated from https://spectralevolution.com/wp-content/uploads/2024/04/RT_hydrang_ref.jpg
	hydrangea_table = np.array([
		0.2,
		0.2,
		0.2,
		0.2,
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
		0.525,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0
	])
	print("Hydrangea leaf")
	spectral_rendering(hydrangea_table)
	
# flowers
# estimated from https://www.researchgate.net/profile/Robert-Gegear/publication/222035355/figure/fig1/AS:632345413578753@1527774303378/The-spectral-reflectance-curves-of-coloured-flowers-and-the-green-array-background-used.png
	
if (args.flowers):
	flower0_table = np.array([
		0.05,
		0.08,
		0.08,
		0.08,
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
		0.08,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0
	])
	print("Blue flower")
	spectral_rendering(flower0_table)
		
	flower1_table = np.array([
		0.08,
		0.08,
		0.08,
		0.08,
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
		0.1,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0
	])
	print("Yellow flower")
	spectral_rendering(flower1_table)
	
	# green flower
	flower2_table = np.array([
		0.0,
		0.0,
		0.0,
		0.0,
		0.120045133,
		0.12964042,
		0.135836894,
		0.139659595,
		0.144615233,
		0.146226008,
		0.1469813,
		0.14724334,
		0.148160479,
		0.149686477,
		0.151150818,
		0.154565044,
		0.167196912,
		0.183543582,
		0.205685959,
		0.238541146,
		0.293554125,
		0.352829107,
		0.380420374,
		0.387256534,
		0.38048203,
		0.366046711,
		0.341592216,
		0.317623267,
		0.293184186,
		0.269269186,
		0.24610177,
		0.226148197,
		0.20576303,
		0.197354629,
		0.199111839,
		0.221816831,
		0.264459982,
		0.330794628,
		0.401653318,
		0.45626553,
		0.482384749,
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
		0.0
	])
	print("Green flower")
	spectral_rendering(flower2_table)
	
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
		0.9375,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0
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
		0.75,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0
	])
	print("Mature Banksia attenuata flower")
	spectral_rendering(banksia1_table)
	
	print("Color contrast between blue and yellow flower: " + str(color_contrast(flower0_table, flower1_table)))
	print("Brightness contrast between blue and yellow flower: " + str(brightness_contrast(flower0_table, flower1_table)))
	print("Color contrast between Banksia flowers: " + str(color_contrast(banksia0_table, banksia1_table)))
	print("Brightness contrast between Banksia flowers: " + str(brightness_contrast(banksia0_table, banksia1_table)))

# sky and daylight
# These are both normalized to "1" at 550 nm. If I try to make "sky" relative to "daylight",
# the values are way too low and probably not meaningful.
if (args.sky):
	daylight_table = np.empty(61)
	for i in range(0, 61):
		w = i*10 + 300
		# normalize to 1 at 550
		daylight_table[i] = blackbody(w / 1000000000, 5800) / blackbody(550 / 1000000000, 5800)
	
	print("Black body spectrum approximating daylight")
	spectral_rendering(daylight_table, light_source=e)
	
	sky_table = np.empty(61)
	for i in range(0, 61):
		w = i*10 + 300
		sky_table[i] = blackbody(w / 1000000000, 5800) * w**-4 / (blackbody(550 / 1000000000, 5800) * 550**-4)
		
	print("Black body spectrum with Rayleigh scattering approximating a blue sky")
	spectral_rendering(sky_table, light_source=e)
	
	sky1_table = np.empty(61)
	for i in range(0, 61):
		w = i*10 + 300
		sky1_table[i] = d65[i] * w**-4 * 5e8
	
	print("D65 with Rayleigh scattering")
	spectral_rendering(sky1_table, light_source=e)

# Kodak Wratten camera filters
# Note these are 300-900 nm rather than 340-830 nm. I removed 800-900 nm because it's probably
# not relevant here.

# optical density (-log10 transmission, not -ln(transmission))
# The names of these end in "d" for density so I don't confuse them with the transmission
# arrays (ending in "t").
# yellow 15
yellow15d = np.array([
	np.nan, # 300
	2.07925,
	1.87288,
	2.39146,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan, # 400
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan, # 500
	1.61358,
	0.58700,
	0.21129,
	0.08958,
	0.05254,
	0.04196,
	0.04196,
	0.04196,
	0.04196,
	0.03667, # 600
	0.03667,
	0.03667,
	0.03667,
	0.04196,
	0.03667,
	0.03667,
	0.03667,
	0.04196,
	0.03667,
	0.03667, # 700
	0.04196,
	0.03137,
	0.04196,
	0.03667,
	0.03667,
	0.03667,
	0.03667,
	0.03667,
	0.03137,
	0.03667 # 800
])

# red 25
red25d = np.array([
	2.41792, # 300
	2.29092,
	2.08983,
	1.81996,
	1.77763,
	1.93108,
	2.19037,
	2.62429,
	np.nan,
	np.nan,
	np.nan, # 400
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan, # 500
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	1.94167,
	0.51292,
	0.13721, # 600
	0.06313,
	0.04725,
	0.04725,
	0.04196,
	0.04196,
	0.04196,
	0.03667,
	0.04196,
	0.03667,
	0.03667, # 700
	0.03667,
	0.03667,
	0.03667,
	0.03667,
	0.03667,
	0.03138,
	0.04196,
	0.03138,
	0.03667,
	0.03667 # 800
])

# blue 47
blue47d = np.array([
	2.80421, # 300
	2.94708,
	2.75658,
	2.53962,
	2.60312,
	2.51317,
	2.20625,
	1.84642,
	1.54479,
	1.24317,
	0.90979, # 400
	0.57112,
	0.37533,
	0.31183,
	0.30654,
	0.33829,
	0.39650,
	0.49175,
	0.61875,
	0.81983,
	1.10029, # 500
	1.51304,
	2.11100,
	2.92062,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan, # 600
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	2.96296,
	2.42850, # 700
	2.01575,
	1.64004,
	1.28021,
	0.91508,
	0.61875,
	0.39121,
	0.24833,
	0.15308,
	0.11075,
	0.08429 # 800
])

# green 58
green58d = np.array([
	np.nan, # 300
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan, # 400
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	2.15333,
	1.47600,
	0.94154,
	0.50763, # 500
	0.30654,
	0.25892,
	0.28008,
	0.34358,
	0.44942,
	0.59758,
	0.81454,
	1.13733,
	1.57125,
	2.07925, # 600
	2.60842,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	3.00000, # 700
	2.24858,
	1.48129,
	0.93625,
	0.58171,
	0.36475,
	0.23246,
	0.16367,
	0.12133,
	0.09488,
	0.07900 # 800
])

# neutral density filters
# 0.1
nd01d = np.array([
	0.33728, # 300
	0.30024,
	0.25790,
	0.22615,
	0.21028,
	0.18911,
	0.17324,
	0.15736,
	0.15207,
	0.14678,
	0.13620, # 400
	0.13090,
	0.12561,
	0.12032,
	0.12032,
	0.11503,
	0.11503,
	0.11503,
	0.10974,
	0.10974,
	0.10974, # 500
	0.10974,
	0.10974,
	0.10974,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445, # 600
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445,
	0.10445, # 700
	0.09915,
	0.09386,
	0.08857,
	0.08328,
	0.08328,
	0.08328,
	0.07799,
	0.07799,
	0.07799,
	0.07799 # 800
])

# 0.2
nd02d = np.array([
	0.52360, # 300
	0.47068,
	0.39660,
	0.35426,
	0.33839,
	0.30664,
	0.28547,
	0.26960,
	0.26431,
	0.25372,
	0.24843, # 400
	0.23785,
	0.22726,
	0.22197,
	0.21668,
	0.21139,
	0.21139,
	0.20610,
	0.20610,
	0.20610,
	0.20081, # 500
	0.20081,
	0.20081,
	0.20081,
	0.20081,
	0.20081,
	0.20081,
	0.20081,
	0.20081,
	0.19551,
	0.19551, # 600
	0.19551,
	0.20081,
	0.20081,
	0.20610,
	0.20610,
	0.20610,
	0.20610,
	0.20610,
	0.20081,
	0.19551, # 700
	0.19022,
	0.16906,
	0.16376,
	0.14789,
	0.14260,
	0.13731,
	0.13731,
	0.13201,
	0.12672,
	0.12672 # 800
])

# 0.3
nd03d = np.array([
	0.84947, # 300
	0.76480,
	0.63780,
	0.56901,
	0.54255,
	0.48964,
	0.44730,
	0.42614,
	0.41555,
	0.39968,
	0.38909, # 400
	0.37322,
	0.36264,
	0.34676,
	0.34147,
	0.33089,
	0.32559,
	0.32030,
	0.32030,
	0.31501,
	0.31501, # 500
	0.31501,
	0.31501,
	0.31501,
	0.31501,
	0.31501,
	0.31501,
	0.32030,
	0.31501,
	0.30443,
	0.30443, # 600
	0.30972,
	0.31501,
	0.32030,
	0.32559,
	0.32559,
	0.32030,
	0.32559,
	0.32559,
	0.32559,
	0.31501, # 700
	0.29914,
	0.27797,
	0.25151,
	0.23034,
	0.21447,
	0.20918,
	0.19859,
	0.19859,
	0.19330,
	0.18801 # 800
])

# 0.4
nd04d = np.array([
	1.03346, # 300
	0.94350,
	0.80592,
	0.73184,
	0.69480,
	0.63659,
	0.58367,
	0.55721,
	0.53605,
	0.52017,
	0.50430, # 400
	0.48313,
	0.46725,
	0.45138,
	0.44080,
	0.43021,
	0.42492,
	0.41963,
	0.41434,
	0.41434,
	0.40905, # 500
	0.40905,
	0.40905,
	0.40905,
	0.40375,
	0.40375,
	0.40905,
	0.40905,
	0.40375,
	0.39317,
	0.39317, # 600
	0.39317,
	0.40375,
	0.40905,
	0.41963,
	0.41963,
	0.41434,
	0.41434,
	0.41963,
	0.41434,
	0.40375, # 700
	0.38259,
	0.35613,
	0.32438,
	0.29792,
	0.28205,
	0.26617,
	0.26088,
	0.25030,
	0.24500,
	0.23971 # 800
])

# 0.5
nd05d = np.array([
	1.13191, # 300
	1.05782,
	0.96257,
	0.90436,
	0.86203,
	0.79324,
	0.72445,
	0.68741,
	0.66624,
	0.64507,
	0.62391, # 400
	0.59745,
	0.58157,
	0.56041,
	0.54982,
	0.53395,
	0.52336,
	0.51807,
	0.51278,
	0.50749,
	0.50220, # 500
	0.50220,
	0.49691,
	0.49691,
	0.49691,
	0.49161,
	0.49161,
	0.49161,
	0.48632,
	0.48103,
	0.47574, # 600
	0.48103,
	0.49161,
	0.49691,
	0.50220,
	0.50220,
	0.49691,
	0.49691,
	0.50220,
	0.49691,
	0.48632, # 700
	0.45986,
	0.42811,
	0.39636,
	0.36991,
	0.34874,
	0.33816,
	0.32757,
	0.31699,
	0.31170,
	0.30641 # 800
])

# 0.6
nd06d = np.array([
	1.33333, # 300
	1.25396,
	1.13225,
	1.05816,
	1.02112,
	0.93116,
	0.85179,
	0.80945,
	0.78300,
	0.75654,
	0.73537, # 400
	0.70362,
	0.68245,
	0.66129,
	0.64541,
	0.62954,
	0.61895,
	0.60837,
	0.60308,
	0.60308,
	0.59779, # 500
	0.59779,
	0.59779,
	0.59779,
	0.59250,
	0.59250,
	0.59779,
	0.60308,
	0.59250,
	0.57662,
	0.57133, # 600
	0.57662,
	0.59250,
	0.60308,
	0.61366,
	0.61366,
	0.61366,
	0.60837,
	0.61366,
	0.60837,
	0.60308, # 700
	0.57133,
	0.52900,
	0.48137,
	0.43904,
	0.41258,
	0.38612,
	0.37554,
	0.35966,
	0.35437,
	0.34908 # 800
])

# 0.7
nd07d = np.array([
	1.54847, # 300
	1.45852,
	1.32622,
	1.24156,
	1.19922,
	1.09868,
	1.00872,
	0.96110,
	0.92935,
	0.89760,
	0.86585, # 400
	0.83410,
	0.80235,
	0.78118,
	0.76002,
	0.74414,
	0.72827,
	0.71768,
	0.70710,
	0.70181,
	0.69652, # 500
	0.69652,
	0.69652,
	0.69122,
	0.68593,
	0.68593,
	0.69122,
	0.69122,
	0.68064,
	0.66477,
	0.65947, # 600
	0.66477,
	0.68064,
	0.69122,
	0.70181,
	0.70181,
	0.70181,
	0.69652,
	0.70181,
	0.69652,
	0.68593, # 700
	0.64889,
	0.60656,
	0.54835,
	0.51131,
	0.47427,
	0.45310,
	0.43722,
	0.42135,
	0.41606,
	0.40547 # 800
])

# 0.8
nd08d = np.array([
	1.73990, # 300
	1.63406,
	1.49648,
	1.40976,
	1.36214,
	1.25101,
	1.14518,
	1.09226,
	1.05522,
	1.01881,
	0.98706, # 400
	0.95002,
	0.91827,
	0.88652,
	0.86535,
	0.84419,
	0.82831,
	0.81773,
	0.80715,
	0.79656,
	0.79127, # 500
	0.79127,
	0.79127,
	0.78598,
	0.78069,
	0.77540,
	0.78069,
	0.78069,
	0.77010,
	0.75148,
	0.74618, # 600
	0.75677,
	0.76735,
	0.78323,
	0.78852,
	0.78852,
	0.78852,
	0.78323,
	0.79381,
	0.78852,
	0.77264, # 700
	0.73560,
	0.68798,
	0.62977,
	0.57685,
	0.53981,
	0.50806,
	0.49218,
	0.47631,
	0.46573,
	0.46043 # 800
])

# 0.9
nd09d = np.array([
	1.88770, # 300
	1.76868,
	1.62620,
	1.53316,
	1.46437,
	1.34735,
	1.23623,
	1.17273,
	1.13189,
	1.09335,
	1.06160, # 400
	1.01927,
	0.98752,
	0.96009,
	0.93363,
	0.91247,
	0.89659,
	0.88601,
	0.87543,
	0.86484,
	0.85955, # 500
	0.85426,
	0.85426,
	0.85128,
	0.84368,
	0.83967,
	0.83838,
	0.84368,
	0.83309,
	0.81873,
	0.81193, # 600
	0.81921,
	0.83174,
	0.84762,
	0.85820,
	0.85820,
	0.85291,
	0.85052,
	0.85291,
	0.85291,
	0.83785, # 700
	0.80193,
	0.74707,
	0.68887,
	0.63595,
	0.59566,
	0.56671,
	0.54804,
	0.53216,
	0.52525,
	0.51629 # 800
])

# 1.0
nd10d = np.array([
	2.16924, # 300
	2.03191,
	1.87489,
	1.76947,
	1.69268,
	1.55559,
	1.43246,
	1.35751,
	1.31309,
	1.27284,
	1.23051, # 400
	1.18817,
	1.14892,
	1.11741,
	1.08966,
	1.06646,
	1.04862,
	1.03686,
	1.02574,
	1.01687,
	1.01158, # 500
	1.00628,
	1.00628,
	1.00628,
	0.99998,
	0.99570,
	0.99570,
	1.00099,
	0.98982,
	0.96924,
	0.96395, # 600
	0.96924,
	0.98512,
	1.00628,
	1.01687,
	1.01687,
	1.01158,
	1.01158,
	1.01158,
	1.01158,
	0.99570, # 700
	0.94808,
	0.88458,
	0.81578,
	0.75758,
	0.71289,
	0.68007,
	0.65886,
	0.64116,
	0.63058,
	0.62143 # 800
])

# 2.0
nd20d = np.array([
	np.nan, # 300
	4.00000,
	3.81479,
	3.73012,
	3.55021,
	3.24858,
	2.98400,
	2.87287,
	2.75117,
	2.66650,
	2.58183, # 400
	2.47071,
	2.39133,
	2.31196,
	2.24846,
	2.20083,
	2.15850,
	2.11617,
	2.08971,
	2.07383,
	2.05267, # 500
	2.03679,
	2.04208,
	2.03679,
	2.01033,
	1.99975,
	2.00504,
	2.01562,
	1.98387,
	1.92037,
	1.90450, # 600
	1.91508,
	1.94154,
	1.97858,
	2.00504,
	2.00504,
	1.99446,
	1.98917,
	1.98387,
	1.98917,
	1.96271, # 700
	1.89392,
	1.76162,
	1.61875,
	1.48117,
	1.38062,
	1.30654,
	1.25362,
	1.21658,
	1.19012,
	1.17425 # 800
])

# 3.0
nd30d = np.array([
	np.nan, # 300
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	3.96296,
	3.86242, # 400
	3.68779,
	3.56608,
	3.47083,
	3.37029,
	3.30150,
	3.23271,
	3.17979,
	3.14804,
	3.10042,
	3.08454, # 500
	3.06337,
	3.05279,
	3.03692,
	3.01575,
	2.98929,
	2.98929,
	2.98929,
	2.96283,
	2.91521,
	2.88346, # 600
	2.88875,
	2.93108,
	2.98400,
	3.02633,
	3.03162,
	3.01046,
	2.99987,
	2.99987,
	3.01575,
	2.97871, # 700
	2.87287,
	2.67708,
	2.46012,
	2.24846,
	2.08971,
	1.97858,
	1.89921,
	1.84100,
	1.80396,
	1.76692 # 800
])

# 4.0
nd40d = np.array([
	np.nan, # 300
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan,
	np.nan, # 400
	4.91108,
	4.76816,
	4.62525,
	4.49186,
	4.40612,
	4.34895,
	4.24415,
	4.20604,
	4.16475,
	4.13617, # 500
	4.11711,
	4.10758,
	4.10123,
	4.06312,
	4.03454,
	4.03454,
	4.05360,
	3.99325,
	3.90115,
	3.86940, # 600
	3.87575,
	3.93291,
	4.00913,
	4.05677,
	4.06630,
	4.04089,
	4.02184,
	4.03136,
	4.03772,
	3.98055, # 700
	3.83446,
	3.58039,
	3.27234,
	3.01827,
	2.80866,
	2.65940,
	2.56412,
	2.48155,
	2.43391,
	2.39262 # 800
])

# transmission
yellow15t = np.empty(51)
for i in range(51):
	if (np.isnan(yellow15d[i])): # nan check
		yellow15t[i] = 0
	else:
		yellow15t[i] = 10**(-yellow15d[i])
red25t = np.empty(51)
for i in range(51):
	if (np.isnan(red25d[i])):
		red25t[i] = 0
	else:
		red25t[i] = 10**(-red25d[i])
blue47t = np.empty(51)
for i in range(51):
	if (np.isnan(blue47d[i])):
		blue47t[i] = 0
	else:
		blue47t[i] = 10**(-blue47d[i])
green58t = np.empty(51)
for i in range(51):
	if (np.isnan(green58d[i])):
		green58t[i] = 0
	else:
		green58t[i] = 10**(-green58d[i])
# neutral density filters -- these do not contain nans except 2.0, 3.0 and 4.0
nd01t = np.empty(51)
for i in range(51):
	nd01t[i] = 10**(-nd01d[i])
nd02t = np.empty(51)
for i in range(51):
	nd02t[i] = 10**(-nd02d[i])
nd03t = np.empty(51)
for i in range(51):
	nd03t[i] = 10**(-nd03d[i])
nd04t = np.empty(51)
for i in range(51):
	nd04t[i] = 10**(-nd04d[i])
nd05t = np.empty(51)
for i in range(51):
	nd05t[i] = 10**(-nd05d[i])
nd06t = np.empty(51)
for i in range(51):
	nd06t[i] = 10**(-nd06d[i])
nd07t = np.empty(51)
for i in range(51):
	nd07t[i] = 10**(-nd07d[i])
nd08t = np.empty(51)
for i in range(51):
	nd08t[i] = 10**(-nd08d[i])
nd09t = np.empty(51)
for i in range(51):
	nd09t[i] = 10**(-nd09d[i])
nd10t = np.empty(51)
for i in range(51):
	nd10t[i] = 10**(-nd10d[i])
nd20t = np.empty(51)
for i in range(51):
	if (np.isnan(nd20d[i])):
		nd20t[i] = 0
	else:
		nd20t[i] = 10**(-nd20d[i])
nd30t = np.empty(51)
for i in range(51):
	if (np.isnan(nd30d[i])):
		nd30t[i] = 0
	else:
		nd30t[i] = 10**(-nd30d[i])
nd40t = np.empty(51)
for i in range(51):
	if (np.isnan(nd40d[i])):
		nd40t[i] = 0
	else:
		nd40t[i] = 10**(-nd40d[i])

if (args.kcheck):
	# plot test
	#xvalues = np.empty(51)
	#for i in range(51):
#		xvalues[i] = i*10 + 300
#	plt.plot(xvalues, red25, 's-r', mec='k', label="red 25")
#	plt.plot(xvalues, yellow15, 'D-y', mec='k', label="yellow 15")
#	plt.plot(xvalues, green58, '^-g', mec='k', label="green 58")
#	plt.plot(xvalues, blue47, 'o-b', mec='k', label="blue 47")
#	plt.xlabel("Wavelength (nm)")
#	plt.ylabel("Optical density")
#	plt.legend()
#	plt.show()
	
	# plot test 2
#	plt.subplot(2, 2, 1)
#	plt.plot(xvalues, red25, 'r')
#	plt.title("Red 25")
#	plt.subplot(2, 2, 2)
#	plt.plot(xvalues, yellow15, 'y')
#	plt.title("Yellow 15")
#	plt.subplot(2, 2, 3)
#	plt.plot(xvalues, green58, 'g')
#	plt.title("Green 58")
#	plt.subplot(2, 2, 4)
#	plt.plot(xvalues, blue47, 'b')
#	plt.title("Blue 47")
	#plt.show()
	
#	plt.subplot(2, 2, 1)
#	plt.plot(xvalues, nd01, color='0.5')
#	plt.plot(xvalues, nd02, color='0.25')
#	plt.plot(xvalues, nd03, 'k')
#	plt.title("0.1-0.3")
#	plt.subplot(2, 2, 2)
#	plt.plot(xvalues, nd04, color='0.5')
#	plt.plot(xvalues, nd05, color='0.25')
#	plt.plot(xvalues, nd06, 'k')
#	plt.title("0.4-0.6")
#	plt.subplot(2, 2, 3)
#	plt.plot(xvalues, nd07, color='0.5')
#	plt.plot(xvalues, nd08, color='0.25')
#	plt.plot(xvalues, nd09, 'k')
#	plt.title("0.7-0.9")
#	plt.subplot(2, 2, 4)
#	plt.plot(xvalues, nd10, color='0.6')
#	plt.plot(xvalues, nd20, color='0.4')
#	plt.plot(xvalues, nd30, color='0.2')
#	plt.plot(xvalues, nd40, 'k')
#	plt.title("1.0-4.0")
#	#plt.show()
	
	# transmission
#	red25t = np.empty(51)
#	for i in range(51):
#		red25t[i] = math.exp(-red25[i])
#	yellow15t = np.empty(51)
#	for i in range(51):
#		yellow15t[i] = math.exp(-yellow15[i])
#	green58t = np.empty(51)
#	for i in range(51):
#		green58t[i] = math.exp(-green58[i])
#	blue47t = np.empty(51)
#	for i in range(51):
#		blue47t[i] = math.exp(-blue47[i])
#	plt.plot(xvalues, red25t*100, 's-r', mec='k', label="red 25")
#	plt.plot(xvalues, yellow15t*100, 'D-y', mec='k', label="yellow 15")
#	plt.plot(xvalues, green58t*100, '^-g', mec='k', label="green 58")
#	plt.plot(xvalues, blue47t*100, 'o-b', mec='k', label="blue 47")
#	plt.xlabel("Wavelength (nm)")
#	plt.ylabel("Transmission (%)")
#	plt.legend()
#	#plt.show()
	
	# "dummy" array for gray as we assume no color filter
	placeholder = np.empty(51)
	for i in range(51):
		placeholder[i] = 1
	
	# brightness tests -- wrap this up in a function so we don't have to copy it a million
	# times
	def brightness_tests(f):
		f_01 = np.empty(51)
		for i in range(51):
			f_01[i] = f[i] * nd01t[i]
		f_02 = np.empty(51)
		for i in range(51):
			f_02[i] = f[i] * nd02t[i]
		f_03 = np.empty(51)
		for i in range(51):
			f_03[i] = f[i] * nd03t[i]
		f_04 = np.empty(51)
		for i in range(51):
			f_04[i] = f[i] * nd04t[i]
		f_05 = np.empty(51)
		for i in range(51):
			f_05[i] = f[i] * nd05t[i]
		f_06 = np.empty(51)
		for i in range(51):
			f_06[i] = f[i] * nd06t[i]
		f_07 = np.empty(51)
		for i in range(51):
			f_07[i] = f[i] * nd07t[i]
		f_08 = np.empty(51)
		for i in range(51):
			f_08[i] = f[i] * nd08t[i]
		f_09 = np.empty(51)
		for i in range(51):
			f_09[i] = f[i] * nd09t[i]
		# 1.0 + others
		f_10 = np.empty(51)
		for i in range(51):
			f_10[i] = f[i] * nd10t[i]
		f_1001 = np.empty(51)
		for i in range(51):
			f_1001[i] = f[i] * nd10t[i] * nd01t[i]
		f_1002 = np.empty(51)
		for i in range(51):
			f_1002[i] = f[i] * nd10t[i] * nd02t[i]
		f_1003 = np.empty(51)
		for i in range(51):
			f_1003[i] = f[i] * nd10t[i] * nd03t[i]
		f_1004 = np.empty(51)
		for i in range(51):
			f_1004[i] = f[i] * nd10t[i] * nd04t[i]
		f_1005 = np.empty(51)
		for i in range(51):
			f_1005[i] = f[i] * nd10t[i] * nd05t[i]
		f_1006 = np.empty(51)
		for i in range(51):
			f_1006[i] = f[i] * nd10t[i] * nd06t[i]
		f_1007 = np.empty(51)
		for i in range(51):
			f_1007[i] = f[i] * nd10t[i] * nd07t[i]
		f_1008 = np.empty(51)
		for i in range(51):
			f_1008[i] = f[i] * nd10t[i] * nd08t[i]
		f_1009 = np.empty(51)
		for i in range(51):
			f_1009[i] = f[i] * nd10t[i] * nd09t[i]
		# 2.0
		f_20 = np.empty(51)
		for i in range(51):
			f_20[i] = f[i] * nd20t[i]
		f_2001 = np.empty(51)
		for i in range(51):
			f_2001[i] = f[i] * nd20t[i] * nd01t[i]
		f_2002 = np.empty(51)
		for i in range(51):
			f_2002[i] = f[i] * nd20t[i] * nd02t[i]
		f_2003 = np.empty(51)
		for i in range(51):
			f_2003[i] = f[i] * nd20t[i] * nd03t[i]
		f_2004 = np.empty(51)
		for i in range(51):
			f_2004[i] = f[i] * nd20t[i] * nd04t[i]
		f_2005 = np.empty(51)
		for i in range(51):
			f_2005[i] = f[i] * nd20t[i] * nd05t[i]
		f_2006 = np.empty(51)
		for i in range(51):
			f_2006[i] = f[i] * nd20t[i] * nd06t[i]
		f_2007 = np.empty(51)
		for i in range(51):
			f_2007[i] = f[i] * nd20t[i] * nd07t[i]
		f_2008 = np.empty(51)
		for i in range(51):
			f_2008[i] = f[i] * nd20t[i] * nd08t[i]
		f_2009 = np.empty(51)
		for i in range(51):
			f_2009[i] = f[i] * nd20t[i] * nd09t[i]
		# 3.0
		f_30 = np.empty(51)
		for i in range(51):
			f_30[i] = f[i] * nd30t[i]
		f_3001 = np.empty(51)
		for i in range(51):
			f_3001[i] = f[i] * nd30t[i] * nd01t[i]
		f_3002 = np.empty(51)
		for i in range(51):
			f_3002[i] = f[i] * nd30t[i] * nd02t[i]
		f_3003 = np.empty(51)
		for i in range(51):
			f_3003[i] = f[i] * nd30t[i] * nd03t[i]
		f_3004 = np.empty(51)
		for i in range(51):
			f_3004[i] = f[i] * nd30t[i] * nd04t[i]
		f_3005 = np.empty(51)
		for i in range(51):
			f_3005[i] = f[i] * nd30t[i] * nd05t[i]
		f_3006 = np.empty(51)
		for i in range(51):
			f_3006[i] = f[i] * nd30t[i] * nd06t[i]
		f_3007 = np.empty(51)
		for i in range(51):
			f_3007[i] = f[i] * nd30t[i] * nd07t[i]
		f_3008 = np.empty(51)
		for i in range(51):
			f_3008[i] = f[i] * nd30t[i] * nd08t[i]
		f_3009 = np.empty(51)
		for i in range(51):
			f_3009[i] = f[i] * nd30t[i] * nd09t[i]
		# 4.0
		f_40 = np.empty(51)
		for i in range(51):
			f_40[i] = f[i] * nd40t[i]
		f_4001 = np.empty(51)
		for i in range(51):
			f_4001[i] = f[i] * nd40t[i] * nd01t[i]
		f_4002 = np.empty(51)
		for i in range(51):
			f_4002[i] = f[i] * nd40t[i] * nd02t[i]
		f_4003 = np.empty(51)
		for i in range(51):
			f_4003[i] = f[i] * nd40t[i] * nd03t[i]
		f_4004 = np.empty(51)
		for i in range(51):
			f_4004[i] = f[i] * nd40t[i] * nd04t[i]
		f_4005 = np.empty(51)
		for i in range(51):
			f_4005[i] = f[i] * nd40t[i] * nd05t[i]
		f_4006 = np.empty(51)
		for i in range(51):
			f_4006[i] = f[i] * nd40t[i] * nd06t[i]
		f_4007 = np.empty(51)
		for i in range(51):
			f_4007[i] = f[i] * nd40t[i] * nd07t[i]
		f_4008 = np.empty(51)
		for i in range(51):
			f_4008[i] = f[i] * nd40t[i] * nd08t[i]
		f_4009 = np.empty(51)
		for i in range(51):
			f_4009[i] = f[i] * nd40t[i] * nd09t[i]
		
		print("Brightness contrast with blue:")
		print("0.1: " + str(brightness_contrast(f_01, blue47t)))
		print("0.2: " + str(brightness_contrast(f_02, blue47t)))
		print("0.3: " + str(brightness_contrast(f_03, blue47t)))
		print("0.4: " + str(brightness_contrast(f_04, blue47t)))
		print("0.5: " + str(brightness_contrast(f_05, blue47t)))
		print("0.6: " + str(brightness_contrast(f_06, blue47t)))
		print("0.7: " + str(brightness_contrast(f_07, blue47t)))
		print("0.8: " + str(brightness_contrast(f_08, blue47t)))
		print("0.9: " + str(brightness_contrast(f_09, blue47t)))
		print("1.0: " + str(brightness_contrast(f_10, blue47t)))
		print("1.0 + 0.1: " + str(brightness_contrast(f_1001, blue47t)))
		print("1.0 + 0.2: " + str(brightness_contrast(f_1002, blue47t)))
		print("1.0 + 0.3: " + str(brightness_contrast(f_1003, blue47t)))
		print("1.0 + 0.4: " + str(brightness_contrast(f_1004, blue47t)))
		print("1.0 + 0.5: " + str(brightness_contrast(f_1005, blue47t)))
		print("1.0 + 0.6: " + str(brightness_contrast(f_1006, blue47t)))
		print("1.0 + 0.7: " + str(brightness_contrast(f_1007, blue47t)))
		print("1.0 + 0.8: " + str(brightness_contrast(f_1008, blue47t)))
		print("1.0 + 0.9: " + str(brightness_contrast(f_1009, blue47t)))
		print("2.0: " + str(brightness_contrast(f_20, blue47t)))
		print("2.0 + 0.1: " + str(brightness_contrast(f_2001, blue47t)))
		print("2.0 + 0.2: " + str(brightness_contrast(f_2002, blue47t)))
		print("2.0 + 0.3: " + str(brightness_contrast(f_2003, blue47t)))
		print("2.0 + 0.4: " + str(brightness_contrast(f_2004, blue47t)))
		print("2.0 + 0.5: " + str(brightness_contrast(f_2005, blue47t)))
		print("2.0 + 0.6: " + str(brightness_contrast(f_2006, blue47t)))
		print("2.0 + 0.7: " + str(brightness_contrast(f_2007, blue47t)))
		print("2.0 + 0.8: " + str(brightness_contrast(f_2008, blue47t)))
		print("2.0 + 0.9: " + str(brightness_contrast(f_2009, blue47t)))
		print("3.0: " + str(brightness_contrast(f_30, blue47t)))
		print("3.0 + 0.1: " + str(brightness_contrast(f_3001, blue47t)))
		print("3.0 + 0.2: " + str(brightness_contrast(f_3002, blue47t)))
		print("3.0 + 0.3: " + str(brightness_contrast(f_3003, blue47t)))
		print("3.0 + 0.4: " + str(brightness_contrast(f_3004, blue47t)))
		print("3.0 + 0.5: " + str(brightness_contrast(f_3005, blue47t)))
		print("3.0 + 0.6: " + str(brightness_contrast(f_3006, blue47t)))
		print("3.0 + 0.7: " + str(brightness_contrast(f_3007, blue47t)))
		print("3.0 + 0.8: " + str(brightness_contrast(f_3008, blue47t)))
		print("3.0 + 0.9: " + str(brightness_contrast(f_3009, blue47t)))
		print("4.0: " + str(brightness_contrast(f_40, blue47t)))
		print("4.0 + 0.1: " + str(brightness_contrast(f_4001, blue47t)))
		print("4.0 + 0.2: " + str(brightness_contrast(f_4002, blue47t)))
		print("4.0 + 0.3: " + str(brightness_contrast(f_4003, blue47t)))
		print("4.0 + 0.4: " + str(brightness_contrast(f_4004, blue47t)))
		print("4.0 + 0.5: " + str(brightness_contrast(f_4005, blue47t)))
		print("4.0 + 0.6: " + str(brightness_contrast(f_4006, blue47t)))
		print("4.0 + 0.7: " + str(brightness_contrast(f_4007, blue47t)))
		print("4.0 + 0.8: " + str(brightness_contrast(f_4008, blue47t)))
		print("4.0 + 0.9: " + str(brightness_contrast(f_4009, blue47t)))
	
	# red
	print("Red 25")
	brightness_tests(red25t)
	
	# yellow
	print("")
	print("Yellow 15")
	brightness_tests(yellow15t)
	
	# green
	print("")
	print("Green 58")
	brightness_tests(green58t)
	
	# gray
	print("")
	print("Gray")
	brightness_tests(placeholder)

# These are used by both --blackbody and --kodak, so we keep them out here.
# red 25
red25_0 = np.empty(51)
for i in range(51):
	red25_0[i] = red25t[i] * nd10t[i] * nd05t[i]
red25_03 = np.empty(51)
for i in range(51):
	red25_03[i] = red25_0[i] * nd03t[i]
red25_07 = np.empty(51)
for i in range(51):
	red25_07[i] = red25_0[i] * nd07t[i]
red25_10 = np.empty(51)
for i in range(51):
	red25_10[i] = red25_0[i] * nd10t[i]

# yellow 15
yellow15_0 = np.empty(51)
for i in range(51):
	yellow15_0[i] = yellow15t[i] * nd20t[i]
yellow15_03 = np.empty(51)
for i in range(51):
	yellow15_03[i] = yellow15_0[i] * nd03t[i]
yellow15_07 = np.empty(51)
for i in range(51):
	yellow15_07[i] = yellow15_0[i] * nd07t[i]
yellow15_10 = np.empty(51)
for i in range(51):
	yellow15_10[i] = yellow15_0[i] * nd10t[i]

# green 58
green58_0 = np.empty(51)
for i in range(51):
	green58_0[i] = green58t[i] * nd10t[i] * nd03t[i]
green58_03 = np.empty(51)
for i in range(51):
	green58_03[i] = green58_0[i] * nd03t[i]
green58_07 = np.empty(51)
for i in range(51):
	green58_07[i] = green58_0[i] * nd07t[i]
green58_10 = np.empty(51)
for i in range(51):
	green58_10[i] = green58_0[i] * nd10t[i]

# blue 47
blue47_0 = blue47t # keep this name for convenience
blue47_03 = np.empty(51)
for i in range(51):
	blue47_03[i] = blue47_0[i] * nd03t[i]
blue47_07 = np.empty(51)
for i in range(51):
	blue47_07[i] = blue47_0[i] * nd07t[i]
blue47_10 = np.empty(51)
for i in range(51):
	blue47_10[i] = blue47_0[i] * nd10t[i]

# gray
gray_0 = np.empty(51)
for i in range(51):
	gray_0[i] = nd20t[i] * nd01t[i]
gray_03 = np.empty(51)
for i in range(51):
	gray_03[i] = gray_0[i] * nd03t[i]
gray_07 = np.empty(51)
for i in range(51):
	gray_07[i] = gray_0[i] * nd07t[i]
gray_10 = np.empty(51)
for i in range(51):
	gray_10[i] = gray_0[i] * nd10t[i]

# used by both --kodak and --munsell
def color_disc(f0, f1, f2, f3, s0, s1, s2, s3, correct=35, trials=40, order=0, alternative='greater'):
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
			contrast = color_contrast(flevel, slevel)
			box[counter][order] = contrast
			#print(box[counter][order])
			print(str(i) + " vs. " + str(j) + ": " + str(contrast))
			if (contrast < 1):
				s += 1
			else:
				d += 1
			counter += 1
	print("Different: " + str(d))
	print("Same: " + str(s))
	if (l1 == m1):
		binomial = binomtest(correct, trials, p=(d + s/2)/16, alternative=alternative)
	else:
		binomial = binomtest(correct, trials, p=(d + s/2)/16, alternative=alternative)
	print("P-value: " + str(binomial.pvalue))
	
	# median contrast
	#print(box)
	#print([box[0][order], box[1][order], box[2][order], box[3][order], box[4][order], box[5][order], box[6][order], box[7][order], box[8][order], box[9][order], box[10][order], box[11][order], box[12][order], box[13][order], box[14][order], box[15][order]])
	mc = statistics.median([box[0][order], box[1][order], box[2][order], box[3][order], box[4][order], box[5][order], box[6][order], box[7][order], box[8][order], box[9][order], box[10][order], box[11][order], box[12][order], box[13][order], box[14][order], box[15][order]])
	print("Median contrast: " + str(mc))
	print("")

if (args.kodak):
	
	# plot
	xvalues = np.empty(51)
	for i in range(51):
		xvalues[i] = i*10 + 300
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
		print("P-value: " + str(binomial.pvalue))
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
	plt.ylabel("Perceived brightness")
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
	
	# integral approximations
	#integral = 0
	#for i in range(30, 90):
	#	integral += blackbody(10*i, args.blackbody)
	#print("Integral using 10-nm bins from 300 to 900: " + str(integral) + " (" + str(100*integral / (energy)) + "% of actual value)")
	#integral = 0
	#for i in range(38, 73):
	#	integral += blackbody(10*i, args.blackbody)
	#print("Integral using 10-nm bins from 380 to 730: " + str(integral) + " (" + str(100*integral / (energy)) + "% of actual value)")
	#for i in range(1, 100000):
	#	integral += blackbody(10*i, args.blackbody)
	#print("Integral using 10-nm bins from 10 to 1000000: " + str(integral) + " (" + str(100*integral / (energy)) + "% of actual value)")
	#integral = 0
	#for i in range(1, 100000):
	#	integral += blackbody(i, args.blackbody)
	#print("Integral using 1-nm bins from 1 to 100000: " + str(integral) + " (" + str(100*integral / (energy)) + "% of actual value)")
	
	# integral approximations with SciPy
	#integral = quad(blackbody, 0, 100000, args=args.blackbody)
	#print("SciPy integral from 0 to 100000: " + str(integral[0]) + " (" + str(100*integral[0] / (energy)) + "% of actual value)")
	#integral = quad(blackbody, 300, 900, args=args.blackbody)
	#print("SciPy integral from 300 to 900: " + str(integral[0]) + " (" + str(100*integral[0] / (energy)) + "% of actual value)")
	#integral = quad(blackbody, 340, 830, args=args.blackbody)
	#print("SciPy integral from 340 to 830: " + str(integral[0]) + " (" + str(100*integral[0] / (energy)) + "% of actual value)")
	#integral = quad(blackbody, 380, 730, args=args.blackbody)
	#print("SciPy integral from 380 to 730: " + str(integral[0]) + " (" + str(100*integral[0] / (energy)) + "% of actual value)")
	
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
	
	# calculate the "underestimation factor"
	# No, forget this. We can do better.
	radiance_full = quad(blackbody, 300, 700, args=args.blackbody)
#	radiance_10nm = 0
#	for i in range(30, 81):
#		radiance_10nm += blackbody(i*10, args.blackbody)
#	luminance_full = 0
#	for i in range(300, 801):
#		luminance_full += blackbody(i, args.blackbody) * (0.68990272*vpt(i, 565) + 0.34832189*vpt(i, 535))*human_filter(i) * 683.002
#	luminance_10nm = 0
#	for i in range(30, 81):
#		luminance_10nm += blackbody(i*10, args.blackbody) * (0.68990272*vpt(i*10, 565) + 0.34832189*vpt(i*10, 535))*human_filter(i*10) * 683.002
	#visible_total = 0
	#visible = np.empty(41)
	#for i in range(41):
	#	visible[i] = blackbody(i+300, args.blackbody) * watts_m2 / (energy)
	#	visible_total += visible[i]
	visible = radiance_full[0]
	x10nm = np.empty(51)
	for i in range(51):
		x10nm[i] = i*10 + 300
	x1nm = np.empty(501)
	for i in range(501):
		x1nm[i] = i + 300
	# interpolate
	red25_0i = np.interp(x1nm, x10nm, red25_0)
	red25_03i = np.interp(x1nm, x10nm, red25_03)
	red25_07i = np.interp(x1nm, x10nm, red25_07)
	red25_10i = np.interp(x1nm, x10nm, red25_10)
	yellow15_0i = np.interp(x1nm, x10nm, yellow15_0)
	yellow15_03i = np.interp(x1nm, x10nm, yellow15_03)
	yellow15_07i = np.interp(x1nm, x10nm, yellow15_07)
	yellow15_10i = np.interp(x1nm, x10nm, yellow15_10)
	green58_0i = np.interp(x1nm, x10nm, green58_0)
	green58_03i = np.interp(x1nm, x10nm, green58_03)
	green58_07i = np.interp(x1nm, x10nm, green58_07)
	green58_10i = np.interp(x1nm, x10nm, green58_10)
	blue47_0i = np.interp(x1nm, x10nm, blue47_0)
	blue47_03i = np.interp(x1nm, x10nm, blue47_03)
	blue47_07i = np.interp(x1nm, x10nm, blue47_07)
	blue47_10i = np.interp(x1nm, x10nm, blue47_10)
	gray_0i = np.interp(x1nm, x10nm, gray_0)
	gray_03i = np.interp(x1nm, x10nm, gray_03)
	gray_07i = np.interp(x1nm, x10nm, gray_07)
	gray_10i = np.interp(x1nm, x10nm, gray_10)
	visible_r0 = 0
	for i in range(300, 701):
		visible_r0 += blackbody(i, args.blackbody) * red25_0i[i-300]
	visible_r03 = 0
	for i in range(300, 701):
		visible_r03 += blackbody(i, args.blackbody) * red25_03i[i-300]
	visible_r07 = 0
	for i in range(300, 701):
		visible_r07 += blackbody(i, args.blackbody) * red25_07i[i-300]
	visible_r10 = 0
	for i in range(300, 701):
		visible_r10 += blackbody(i, args.blackbody) * red25_10i[i-300]
	visible_y0 = 0
	for i in range(300, 701):
		visible_y0 += blackbody(i, args.blackbody) * yellow15_0i[i-300]
	visible_y03 = 0
	for i in range(300, 701):
		visible_y03 += blackbody(i, args.blackbody) * yellow15_03i[i-300]
	visible_y07 = 0
	for i in range(300, 701):
		visible_y07 += blackbody(i, args.blackbody) * yellow15_07i[i-300]
	visible_y10 = 0
	for i in range(300, 701):
		visible_y10 += blackbody(i, args.blackbody) * yellow15_10i[i-300]
	visible_g0 = 0
	for i in range(300, 701):
		visible_g0 += blackbody(i, args.blackbody) * green58_0i[i-300]
	visible_g03 = 0
	for i in range(300, 701):
		visible_g03 += blackbody(i, args.blackbody) * green58_03i[i-300]
	visible_g07 = 0
	for i in range(300, 701):
		visible_g07 += blackbody(i, args.blackbody) * green58_07i[i-300]
	visible_g10 = 0
	for i in range(300, 701):
		visible_g10 += blackbody(i, args.blackbody) * green58_10i[i-300]
	visible_b0 = 0
	for i in range(300, 701):
		visible_b0 += blackbody(i, args.blackbody) * blue47_0i[i-300]
	visible_b03 = 0
	for i in range(300, 701):
		visible_b03 += blackbody(i, args.blackbody) * blue47_03i[i-300]
	visible_b07 = 0
	for i in range(300, 701):
		visible_b07 += blackbody(i, args.blackbody) * blue47_07i[i-300]
	visible_b10 = 0
	for i in range(300, 701):
		visible_b10 += blackbody(i, args.blackbody) * blue47_10i[i-300]
	visible_gray0 = 0
	for i in range(300, 701):
		visible_gray0 += blackbody(i, args.blackbody) * gray_0i[i-300]
	visible_gray03 = 0
	for i in range(300, 701):
		visible_gray03 += blackbody(i, args.blackbody) * gray_03i[i-300]
	visible_gray07 = 0
	for i in range(300, 701):
		visible_gray07 += blackbody(i, args.blackbody) * gray_07i[i-300]
	visible_gray10 = 0
	for i in range(300, 701):
		visible_gray10 += blackbody(i, args.blackbody) * gray_10i[i-300]
	# scale
	#visible_filtered *= luminance_full / luminance_10nm
	#print(luminance_full / luminance_10nm)
	
	# As above, this is radiance, so we multiply by the sr value we found earlier to get W/m^2
	#visible_total = visible_total * sr
	#visible_filtered_total = visible_filtered_total * sr
	
	# lx/cd
	visible_bcd0 = 0
	for i in range(300, 801):
		visible_bcd0 += blackbody(i, args.blackbody) * blue47_0i[i-300] *  luminosity[i-300] * 683.002
	visible_bcd03 = 0
	for i in range(300, 801):
		visible_bcd03 += blackbody(i, args.blackbody) * blue47_03i[i-300] *  luminosity[i-300] * 683.002
	visible_bcd07 = 0
	for i in range(300, 801):
		visible_bcd07 += blackbody(i, args.blackbody) * blue47_07i[i-300] *  luminosity[i-300] * 683.002
	visible_bcd10 = 0
	for i in range(300, 801):
		visible_bcd10 += blackbody(i, args.blackbody) * blue47_10i[i-300] *  luminosity[i-300] * 683.002
	#visible_filtered_cd *= luminance_full / luminance_10nm
	visible_lx0 = visible_bcd0 * math.pi
	visible_lx03 = visible_bcd03 * math.pi
	visible_lx07 = visible_bcd07 * math.pi
	visible_lx10 = visible_bcd10 * math.pi
	# other colors
	# red
	visible_rcd0 = 0
	for i in range(300, 801):
		visible_rcd0 += blackbody(i, args.blackbody) * red25_0i[i-300] *  luminosity[i-300] * 683.002
	visible_rcd03 = 0
	for i in range(300, 801):
		visible_rcd03 += blackbody(i, args.blackbody) * red25_03i[i-300] *  luminosity[i-300] * 683.002
	visible_rcd07 = 0
	for i in range(300, 801):
		visible_rcd07 += blackbody(i, args.blackbody) * red25_07i[i-300] *  luminosity[i-300] * 683.002
	visible_rcd10 = 0
	for i in range(300, 801):
		visible_rcd10 += blackbody(i, args.blackbody) * red25_10i[i-300] *  luminosity[i-300] * 683.002
	# yellow
	visible_ycd0 = 0
	for i in range(300, 801):
		visible_ycd0 += blackbody(i, args.blackbody) * yellow15_0i[i-300] *  luminosity[i-300] * 683.002
	visible_ycd03 = 0
	for i in range(300, 801):
		visible_ycd03 += blackbody(i, args.blackbody) * yellow15_03i[i-300] *  luminosity[i-300] * 683.002
	visible_ycd07 = 0
	for i in range(300, 801):
		visible_ycd07 += blackbody(i, args.blackbody) * yellow15_07i[i-300] *  luminosity[i-300] * 683.002
	visible_ycd10 = 0
	for i in range(300, 801):
		visible_ycd10 += blackbody(i, args.blackbody) * yellow15_10i[i-300] *  luminosity[i-300] * 683.002
	# green
	visible_gcd0 = 0
	for i in range(300, 801):
		visible_gcd0 += blackbody(i, args.blackbody) * green58_0i[i-300] *  luminosity[i-300] * 683.002
	visible_gcd03 = 0
	for i in range(300, 801):
		visible_gcd03 += blackbody(i, args.blackbody) * green58_03i[i-300] *  luminosity[i-300] * 683.002
	visible_gcd07 = 0
	for i in range(300, 801):
		visible_gcd07 += blackbody(i, args.blackbody) * green58_07i[i-300] *  luminosity[i-300] * 683.002
	visible_gcd10 = 0
	for i in range(300, 801):
		visible_gcd10 += blackbody(i, args.blackbody) * green58_10i[i-300] *  luminosity[i-300] * 683.002
	# gray
	visible_graycd0 = 0
	for i in range(300, 801):
		visible_graycd0 += blackbody(i, args.blackbody) * gray_0i[i-300] *  luminosity[i-300] * 683.002
	visible_graycd03 = 0
	for i in range(300, 801):
		visible_graycd03 += blackbody(i, args.blackbody) * gray_03i[i-300] *  luminosity[i-300] * 683.002
	visible_graycd07 = 0
	for i in range(300, 801):
		visible_graycd07 += blackbody(i, args.blackbody) * gray_07i[i-300] *  luminosity[i-300] * 683.002
	visible_graycd10 = 0
	for i in range(300, 801):
		visible_graycd10 += blackbody(i, args.blackbody) * gray_10i[i-300] *  luminosity[i-300] * 683.002
	
	# candela equivalents weighted by a custom luminosity function in place of CIE
	# to compare brightness for another species
	
	# blue
	visible_bcde0 = 0
	for i in range(300, 801):
		visible_bcde0 += blackbody(i, args.blackbody) * blue47_0i[i-300] *  sensitivity(i) * 683.002
	visible_bcde03 = 0
	for i in range(300, 801):
		visible_bcde03 += blackbody(i, args.blackbody) * blue47_03i[i-300] * sensitivity(i) * 683.002
	visible_bcde07 = 0
	for i in range(300, 801):
		visible_bcde07 += blackbody(i, args.blackbody) * blue47_07i[i-300] * sensitivity(i) * 683.002
	visible_bcde10 = 0
	for i in range(300, 801):
		visible_bcde10 += blackbody(i, args.blackbody) * blue47_10i[i-300] * sensitivity(i) * 683.002
	# red
	visible_rcde0 = 0
	for i in range(300, 801):
		visible_rcde0 += blackbody(i, args.blackbody) * red25_0i[i-300] * sensitivity(i) * 683.002
	visible_rcde03 = 0
	for i in range(300, 801):
		visible_rcde03 += blackbody(i, args.blackbody) * red25_03i[i-300] * sensitivity(i) * 683.002
	visible_rcde07 = 0
	for i in range(300, 801):
		visible_rcde07 += blackbody(i, args.blackbody) * red25_07i[i-300] * sensitivity(i) * 683.002
	visible_rcde10 = 0
	for i in range(300, 801):
		visible_rcde10 += blackbody(i, args.blackbody) * red25_10i[i-300] * sensitivity(i) * 683.002
	# yellow
	visible_ycde0 = 0
	for i in range(300, 801):
		visible_ycde0 += blackbody(i, args.blackbody) * yellow15_0i[i-300] * sensitivity(i) * 683.002
	visible_ycde03 = 0
	for i in range(300, 801):
		visible_ycde03 += blackbody(i, args.blackbody) * yellow15_03i[i-300] * sensitivity(i) * 683.002
	visible_ycde07 = 0
	for i in range(300, 801):
		visible_ycde07 += blackbody(i, args.blackbody) * yellow15_07i[i-300] * sensitivity(i) * 683.002
	visible_ycde10 = 0
	for i in range(300, 801):
		visible_ycde10 += blackbody(i, args.blackbody) * yellow15_10i[i-300] * sensitivity(i) * 683.002
	# green
	visible_gcde0 = 0
	for i in range(300, 801):
		visible_gcde0 += blackbody(i, args.blackbody) * green58_0i[i-300] * sensitivity(i) * 683.002
	visible_gcde03 = 0
	for i in range(300, 801):
		visible_gcde03 += blackbody(i, args.blackbody) * green58_03i[i-300] * sensitivity(i) * 683.002
	visible_gcde07 = 0
	for i in range(300, 801):
		visible_gcde07 += blackbody(i, args.blackbody) * green58_07i[i-300] * sensitivity(i) * 683.002
	visible_gcde10 = 0
	for i in range(300, 801):
		visible_gcde10 += blackbody(i, args.blackbody) * green58_10i[i-300] * sensitivity(i) * 683.002
	# gray
	visible_graycde0 = 0
	for i in range(300, 801):
		visible_graycde0 += blackbody(i, args.blackbody) * gray_0i[i-300] * sensitivity(i) * 683.002
	visible_graycde03 = 0
	for i in range(300, 801):
		visible_graycde03 += blackbody(i, args.blackbody) * gray_03i[i-300] * sensitivity(i) * 683.002
	visible_graycde07 = 0
	for i in range(300, 801):
		visible_graycde07 += blackbody(i, args.blackbody) * gray_07i[i-300] * sensitivity(i) * 683.002
	visible_graycde10 = 0
	for i in range(300, 801):
		visible_graycde10 += blackbody(i, args.blackbody) * gray_10i[i-300] * sensitivity(i) * 683.002
	
	print("Radiance from 300-800 nm (W/m^2/sr): " + str(visible))
	print("Irradiance from 300-800 nm (W/m^2): " + str(visible*sr))
	print("Irradiance from 300-800 nm (uW/cm^2): " + str(visible*sr*1000000/100**2))
#	print("Irradiance from 300-800 nm filtered through blue 47 (W/m^2): " + str(visible_filtered*sr))
#	print("Irradiance from 300-800 nm filtered through blue 47 (uW/cm^2): " + str(visible_filtered*sr*1000000/100**2))
	#print("Luminance from 300-800 nm (cd/m^2): " + str(visible_filtered_cd))
	#print("Illuminance from 300-800 nm (lx): " + str(visible_filtered_lx))
	
	# scaled to specified watt number
	print("")
	print("Scaled radiance/irradiance:")
	print("Radiance from 300-800 nm (W/m^2/sr): " + str(scale*visible))
	radiance = 0
	for i in range(501):
		radiance += ia[i]
	print("Radiance from 300-800 nm (photons/sec/sr): " + str(radiance))
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
	#print("Irradiance from 300-800 nm filtered through blue 47 (uW/cm^2): " + str(scale*visible_filtered*sr*1000000/100**2))
	print("Luminance from 300-800 nm for blue 47 (cd/m^2):")
	print("0: " + str(scale*visible_bcd0))
	print("0.3: " + str(scale*visible_bcd03))
	print("0.7: " + str(scale*visible_bcd07))
	print("1.0: " + str(scale*visible_bcd10))
	print("Luminance from 300-800 nm for red 25 (cd/m^2):")
	print("0: " + str(scale*visible_rcd0))
	print("0.3: " + str(scale*visible_rcd03))
	print("0.7: " + str(scale*visible_rcd07))
	print("1.0: " + str(scale*visible_rcd10))
	print("Luminance from 300-800 nm for yellow 15 (cd/m^2):")
	print("0: " + str(scale*visible_ycd0))
	print("0.3: " + str(scale*visible_ycd03))
	print("0.7: " + str(scale*visible_ycd07))
	print("1.0: " + str(scale*visible_ycd10))
	print("Luminance from 300-800 nm for green 58 (cd/m^2):")
	print("0: " + str(scale*visible_gcd0))
	print("0.3: " + str(scale*visible_gcd03))
	print("0.7: " + str(scale*visible_gcd07))
	print("1.0: " + str(scale*visible_gcd10))
	print("Luminance from 300-800 nm for gray (cd/m^2):")
	print("0: " + str(scale*visible_graycd0))
	print("0.3: " + str(scale*visible_graycd03))
	print("0.7: " + str(scale*visible_graycd07))
	print("1.0: " + str(scale*visible_graycd10))
	print("Illuminance (lx):")
	print("0: " + str(scale*visible_lx0))
	print("0.3: " + str(scale*visible_lx03))
	print("0.7: " + str(scale*visible_lx07))
	print("1.0: " + str(scale*visible_lx10))
	# alternatively it may depend on the size of the circle depending on what a
	# foot lambert is
	#print("Luminance/illuminance scaled to size of circle:")
	#print("Luminance (lx/sr):")
	#print("0: " + str(scale*visible_filtered_lx0/sr))
	#print("0.3: " + str(scale*visible_filtered_lx03/sr))
	#print("0.7: " + str(scale*visible_filtered_lx07/sr))
	#print("1.0: " + str(scale*visible_filtered_lx10/sr))
	#print("Illuminance (cd*sr/m^2):")
	#print("0: " + str(scale*visible_filtered_cd0*sr))
	#print("0.3: " + str(scale*visible_filtered_cd03*sr))
	#print("0.7: " + str(scale*visible_filtered_cd07*sr))
	#print("1.0: " + str(scale*visible_filtered_cd10*sr))
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
	plt.subplot(1, 2, 2)
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_rcde0, scale*visible_rcde03, scale*visible_rcde07, scale*visible_rcde10], 'sr', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_ycde0, scale*visible_ycde03, scale*visible_ycde07, scale*visible_ycde10], 'Dy', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_gcde0, scale*visible_gcde03, scale*visible_gcde07, scale*visible_gcde10], '^g', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_bcde0, scale*visible_bcde03, scale*visible_bcde07, scale*visible_bcde10], 'ob', mec='k')
	plt.plot([0, 0.3, 0.7, 1.0], [scale*visible_graycde0, scale*visible_graycde03, scale*visible_graycde07, scale*visible_graycde10], marker='v', linestyle='', color='gray', mec='k')
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
		0.056283,
		0.058632,
		0.062006,
		0.066158, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.115894,
		0.119029,
		0.127256,
		0.144525, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.196282,
		0.203308,
		0.214796,
		0.231492, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.23233,
		0.223291,
		0.227917,
		0.250112, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.148741,
		0.149692,
		0.155841,
		0.172085, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.225293,
		0.229032,
		0.236621,
		0.246813, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.350633,
		0.353927,
		0.362951,
		0.379432, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.462779,
		0.469476,
		0.477307,
		0.489983, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.309138,
		0.305939,
		0.302531, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.462856,
		0.461125,
		0.458641,
		0.455757, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.671756,
		0.673282,
		0.672417,
		0.670845, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.884059,
		0.885425,
		0.882944,
		0.880923, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.09505,
		0.100989,
		0.100258,
		0.098844, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.145016,
		0.144087,
		0.140965, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.286497,
		0.29664,
		0.295497,
		0.293521, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
		0.475862,
		0.491195,
		0.491365,
		0.487952, # 730
		np.nan, # 740
		np.nan, # 750
		np.nan, # 760
		np.nan, # 770
		np.nan, # 780
		np.nan, # 790
		np.nan # 800
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
	
	# plot spectra
	xvalues = np.empty(51)
	for i in range(0, 51):
		xvalues[i] = i*10 + 300
	
	plt.subplot(2, 2, 1)
	plt.plot(xvalues, m5b4, color='deepskyblue', label="5B4/6")
	plt.plot(xvalues, m5b5, color='deepskyblue', label="5B5/6")
	plt.plot(xvalues, m5b6, color='deepskyblue', label="5B6/6")
	plt.plot(xvalues, m5b7, color='deepskyblue', label="5B7/6")
	#plt.legend()
	plt.title("5B")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 2)
	plt.plot(xvalues, m75b5, color='dodgerblue', label="7.5B5/4")
	plt.plot(xvalues, m75b6, color='dodgerblue', label="7.5B6/4")
	plt.plot(xvalues, m75b7, color='dodgerblue', label="7.5B7/4")
	plt.plot(xvalues, m75b8, color='dodgerblue', label="7.5B8/4")
	#plt.legend()
	plt.title("7.5B")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 3)
	plt.plot(xvalues, m10yr5, color='orange', label="10YR5/10")
	plt.plot(xvalues, m10yr6, color='orange', label="10YR6/10")
	plt.plot(xvalues, m10yr7, color='orange', label="10YR7/10")
	plt.plot(xvalues, m10yr8, color='orange', label="10YR8/10")
	#plt.legend()
	plt.title("10YR")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 4)
	plt.plot(xvalues, m5gy5, color='chartreuse', label="5GY5/10")
	plt.plot(xvalues, m5gy6, color='chartreuse', label="5GY6/10")
	plt.plot(xvalues, m5gy7, color='chartreuse', label="5GY7/10")
	plt.plot(xvalues, m5gy8, color='chartreuse', label="5GY8/10")
	#plt.legend()
	plt.title("5GY")
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
	box = np.empty((16, 3))
	
	# 5B vs. 10YR
	print("5B vs. 10YR")
	color_disc(m5b4, m5b5, m5b6, m5b7, m10yr5, m10yr6, m10yr7, m10yr8, correct=8, trials=16, order=0, alternative='less')
	
	# 7.5B vs. 10YR
	print("7.5B vs. 10YR")
	color_disc(m75b5, m75b6, m75b7, m75b8, m10yr5, m10yr6, m10yr7, m10yr8, correct=8, trials=16, order=0, alternative='less')
	
	# 5GY vs. 10YR
	print("5GY vs. 10YR")
	color_disc(m5gy5, m5gy6, m5gy7, m5gy8, m10yr5, m10yr6, m10yr7, m10yr8, correct=8, trials=16, order=0, alternative='less')

# print execution time
print("%s seconds" % (time.time() - start_time))
