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
parser.add_argument("--lp", type=float, default=0.62, help="proportion of L cones")
parser.add_argument("--mp", type=float, default=0.31, help="proportion of M cones")
parser.add_argument("--sp", type=float, default=0.07, help="proportion of S cones")
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
parser.add_argument("--convergence", type=float, default=1, help="cone-to-ganglion cell convergence ratio")
parser.add_argument("--ff", type=float, default=22.4, help="flicker fusion frequency (Hz)")
parser.add_argument("--osd", type=float, default=1.2, help="outer segment diameter (um)")
parser.add_argument("--fl", type=float, default=22.3, help="focal length (mm)")
parser.add_argument("--pupil", type=float, default=11, help="pupil diameter (mm)")
parser.add_argument("--od", type=float, default=0.015, help="optical density of pigment (um)")
parser.add_argument("--osl", type=float, default=13.4, help="outer segment length (um)")
parser.add_argument("--dist", type=float, default=5, help="distance from light source (cm)")
parser.add_argument("--width", type=float, default=3.8, help="width of object (cm)")
parser.add_argument("--watts", type=int, default=12, help="number of watts produced by light source")
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
	density = a*math.exp((b - w) / c) + d*math.exp((e - w) / f)
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

# function
def opossum_filter(w):
	return template_filter(w, *opossum_fit[0])

def lens_filter(w):
	if (args.filter == "human"):
		return human_filter(w)
	elif (args.filter == "opossum"):
		return opossum_filter(w)
	return 1

def sensitivity(w, l=l1, m=m1, s=s1):
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
	return 1e-9*2*h*c**2*l**-5 / (math.exp((h*c) / (k*l*t)) - 1)

# normal distribution
def normal(mu, std_dev, x):
	return (1 / std_dev*math.sqrt(2*math.pi)) * math.exp((-1/2)*((x - mu)/std_dev)**2)

# light sources -- standardize to 300-900
d65 = np.zeros(61)
for i in range(4, 54):
	d65[i] = spectral_constants.REF_ILLUM_TABLE["d65"][i-4]

e = np.empty(61)
for i in range(0, 61):
	e[i] = 100

a = np.zeros(61)
for i in range(4, 54):
	a[i] = spectral_constants.REF_ILLUM_TABLE["a"][i-4]

# incandescent lighting based on CIE A illuminant
incandescent = np.empty(61)
for i in range(0, 61):
	w = i*10 + 300
	# normalize to 100% at 560 nm
	incandescent[i] = 100*blackbody(w, 2856) / blackbody(560, 2856)

# absolute number of photons
# Since the output of blackbody() is joules/m^3/sec/sr, this should be photons/m^3/sec/sr.
# I think what we really want is square meters. (No, the meters cancel, see the quantum
# noise calculations)
ia = np.empty(61)
# scale it down to what amount of energy would be produced by the specified number of
# watts

# find total energy produced using Stefan-Boltzmann constant
sigma = 5.670374419e-8
energy = sigma * 2856**4

# then find the amount of energy we "should" have in terms of the sphere corresponding to
# the distance we chose
sphere_area = 4*math.pi*args.dist**2 # sphere surface area
watts_cm2 = args.watts / sphere_area
watts_m2 = watts_cm2 * 100**2
scale = watts_m2 / energy # scaling factor -- units should cancel (W/m^2)

for i in range(0, 61):
	w = i*10 + 300 # nanometers
	ia[i] = scale*1e-9*w*blackbody(w, 2856) / (h*c) # meters * joules/m^2/nm/sec/sr / joule*sec * m/s
	#print(ia[i])

if (args.blackbody != 0):
	#xvalues = yvalues = np.empty(61) don't do this, it won't do what you think it does
	xvalues = np.empty(61)
	yvalues = np.empty(61)
	for i in range(61):
		w = i*10 + 300
		xvalues[i] = w
		yvalues[i] = blackbody(w, args.blackbody)
	plt.plot(xvalues, yvalues)
	#plt.plot(xvalues, ia)
	#plt.plot(xvalues, a)
	#plt.plot(xvalues, incandescent)
	#plt.plot(xvalues, d65)
	print(sphere_area)
	print(watts_cm2)
	print(watts_m2)
	print(energy)
	print(sigma * args.blackbody**4)
	plt.show()
	
	# how much energy is within the (human-)visible spectrum? Incandescent bulbs usually
	# have 10%. https://physlab.org/wp-content/uploads/2016/03/Planck_ref8.pdf
	#energy = 0
	#for i in range(380, 780): # 380-780 nm
		#energy += blackbody(i, args.blackbody)
		#print(energy)
	#print(str(energy) + " J/m^3/s/sr = W/m^3/sr")
	#print(0.6 / energy)
	#print(1e-9*300 / (h*c))
	#print((0.6 / energy) * 1e-9*300 / (h*c))
	#print(blackbody(300, args.blackbody))
	#print(blackbody(300, args.blackbody) * scale * 1e-9*300 / (h*c))
	#print(blackbody(600, args.blackbody))
	#print(blackbody(600, args.blackbody) * scale * 1e-9*600 / (h*c))
	#print(h*c)

# white point
if (args.white == "d65"):
	wp = d65
if (args.white == "a"):
	wp = a
elif (args.white == "i"):
	wp = incandescent
else:
	wp = e

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
		print("Warning: math overflow, clipping to 2.2250738585072014e-308")
		return 2.2250738585072014e-308
	
	return alpha + beta

# estimate hue, saturation and lightness for a spectral power distribution
# Lens filtering is now used for colors properly by normalizing each signal according to
# the light source.
def spectral_rendering(table, light_source=wp/100):
	table_l = 0
	table_m = 0
	table_s = 0
	brightness = 0
	
	for i in range(0, light_source.shape[0]):
		w = i*10 + 300 # wavelength
		table_l += vpt(w, l1) * table[i] * light_source[i] * lens_filter(w)
		# remove either M or S for dichromacy
		if (m1 != l1):
			table_m += vpt(w, m1) * table[i] * light_source[i] * lens_filter(w)
		if (s1 != m1):
			table_s += vpt(w, s1) * table[i] * light_source[i] * lens_filter(w)
		
		# brightness
		brightness += sensitivity(w) * table[i] * light_source[i] * lens_filter(w)
	
	# normalize according to provided light source
	n = 0
	wpl = 0
	wpm = 0
	wps = 0
	for i in range(0, light_source.shape[0]):
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
	# reported agree reasonably well with what the hue should be.
	if (l1 != m1 != s1):
		print("Cone response: l=" + str(table_l) + ", m=" + str(table_m) + ", s=" + str(table_s))
		
		# This should be the nearest wavelength on a line from the white point.
		match = 300
		x0 = (table_l - table_s) / (table_l + table_m + table_s)
		y0 = (table_m - table_l - table_s) / (table_l + table_m + table_s)
		white_l = 0
		white_m = 0
		white_s = 0
		for i in range(wp.shape[0]):
			w = i*10 + 300
			white_l += vpt(w, l1) * wp[i]
			white_m += vpt(w, m1) * wp[i]
			white_s += vpt(w, s1) * wp[i]
		wx = (white_l - white_s) / (white_l + white_m + white_s)
		wy = (white_m - white_l - white_s) / (white_l + white_m + white_s)
		
		for i in range(300, 800):
			x1 = (vpt(match, l1) - vpt(match, s1)) / (vpt(match, l1) + vpt(match, m1) + vpt(match, s1))
			y1 = (vpt(match, m1) - vpt(match, l1) - vpt(match, s1)) / (vpt(match, l1) + vpt(match, m1) + vpt(match, s1))
			x2 = (vpt(i, l1) - vpt(i, s1)) / (vpt(i, l1) + vpt(i, m1) + vpt(i, s1))
			y2 = (vpt(i, m1) - vpt(i, l1) - vpt(i, s1)) / (vpt(i, l1) + vpt(i, m1) + vpt(i, s1))
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
		white_l = 0
		white_s = 0
		for i in range(wp.shape[0]):
			w = i*10 + 300
			if (w >= 300):
				white_l += vpt(w, l1) * wp[i]
				white_s += vpt(w, s1) * wp[i]
		dist = pos - (white_l - white_s) / (white_l + white_s)
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
# Neither of these models do both of these.
def color_contrast(table1, table2, quantum_noise=args.qn):
	
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
	for i in range(0, wp.shape[0]):
		w = i*10 + 300
		ql1 += vpt(w, l1) * table1[i] * wp[i] * lens_filter(w)
		qm1 += vpt(w, m1) * table1[i] * wp[i] * lens_filter(w)
		qs1 += vpt(w, s1) * table1[i] * wp[i] * lens_filter(w)
		wpl += vpt(w, l1) * wp[i] * lens_filter(w)
		wpm += vpt(w, m1) * wp[i] * lens_filter(w)
		wps += vpt(w, s1) * wp[i] * lens_filter(w)
	for i in range(0, wp.shape[0]):
		w = i*10 + 300
		ql2 += vpt(w, l1) * table2[i] * wp[i] * lens_filter(w)
		qm2 += vpt(w, m1) * table2[i] * wp[i] * lens_filter(w)
		qs2 += vpt(w, s1) * table2[i] * wp[i] * lens_filter(w)
	
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
		# I am not using the terms for transduction efficiency, oil droplets or the various
		# refractive indices that appear in the two chicken studies because introducing
		# more probably wrong estimates than I need to is a bad idea and I had a hard
		# enough time figuring out the rest of this. This is like 1000% more math than I
		# usually bother with. Accounting for oil droplets would be even harder than it
		# sounds because in opossums some cones have them and others don't.
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
		
		# Estimates used:
		# Where necessary, I've converted centimeters, millimeters and micrometers to meters.
		# v = 3000 * 450: cone density for Didelphis aurita (Ahnelt et al. 2009), retinal
		# area for cats (Vaney 1985). This number is probably an overestimate because
		# 3000/mm^2 is the maximum density, not the average density. (Not needed, see
		# above)
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
		# f = 22.3 / 1000: 22.3 mm, cats (axial length, Thomasy et al. 2016 -- note axial
		# length and focal length are not necessarily the same and I don't really know what
		# focal length is)
		# D = 11 / 1000: 11 mm, maximum size of rabbit pupil (Peiffer et al. 1994)
		# O: Thylamys elegans (see the lens filter functions)
		# k = 0.015: mice (Yin et al. 2013)
		# l = 13.4: mice
		# sphere = 4*math.pi*5**2: the filters were located 5 cm away from the light bulbs
		# sr = math.pi*3.8**2 / sphere: 3.8 cm diameter circle
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
		# transmittance of ocular media: 80% (estimate for chickens)
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
		sphere = 4*math.pi*args.dist**2
		sr = math.pi*args.width**2 / sphere
		
		# lens filtering
		T = 0.27 / lens_filter(700)
	
		aql1 = 0
		aqm1 = 0
		aqs1 = 0
		aql2 = 0
		aqm2 = 0
		aqs2 = 0
		for i in range(wp.shape[0]):
			w = i*10 + 300
			aql1 += v*(math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, l1)*l)) * vpt(w, l1) * table1[i] * ia[i]*sr * lens_filter(w)
			aqm1 += v*(math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, m1)*l)) * vpt(w, m1) * table1[i] * ia[i]*sr * lens_filter(w) 
			aqs1 += v*(math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, s1)*l)) * vpt(w, s1) * table1[i] * ia[i]*sr * lens_filter(w)
		for i in range(wp.shape[0]):
			w = i*10 + 300
			aql2 += v*(math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, l1)*l)) * vpt(w, l1) * table2[i] * ia[i]*sr * lens_filter(w)
			aqm2 += v*(math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, m1)*l)) * vpt(w, m1) * table2[i] * ia[i]*sr * lens_filter(w) 
			aqs2 += v*(math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, s1)*l)) * vpt(w, s1) * table2[i] * ia[i]*sr * lens_filter(w)
		
		# wavelength interval (see The Optics of Life)
		aql1 = aql1 * 600
		aqm1 = aqm1 * 600
		aqs1 = aqs1 * 600
		aql2 = aql2 * 600
		aqm2 = aqm2 * 600
		aqs2 = aqs2 * 600
		
		el = math.sqrt((1 / aql1 + 1 / aql2) + 2*wl**2)
		em = math.sqrt((1 / aqm1 + 1 / aqm2) + 2*wm**2)
		es = math.sqrt((1 / aqs1 + 1 / aqs2) + 2*ws**2)
		
		if (args.qcheck):
			print("aq1:")
			print(aql1)
			print(aqm1)
			print(aqs1)
			print("aq2:")
			print(aql2)
			print(aqm2)
			print(aqs2)
			print("noise:")
			print(el)
			print(em)
			print(es)
			print("Weber fractions:")
			print(wl)
			print(wm)
			print(ws)
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
	wpt = 0
	
	ql1 = 0 # L cones
	qm1 = 0 # M cones
	qr1 = 0 # rods
	ql2 = 0
	qm2 = 0
	qr2 = 0
	for i in range(wp.shape[0]):
		w = i*10 + 300
		wpt += sensitivity(w, l1) * wp[i]
		ql1 += vpt(w, l1) * table1[i] * wp[i] * lens_filter(w)
		qm1 += vpt(w, m1) * table1[i] * wp[i] * lens_filter(w)
		qr1 += vpt(w, rod) * table1[i] * wp[i] * lens_filter(w)
	for i in range(wp.shape[0]):
		w = i*10 + 300
		ql2 += vpt(w, l1) * table2[i] * wp[i] * lens_filter(w)
		qm2 += vpt(w, m1) * table2[i] * wp[i] * lens_filter(w)
		qr2 += vpt(w, rod) * table2[i] * wp[i] * lens_filter(w)
	
	if (l1 != m1):
		df = math.log((args.lb*ql1 + args.mb*qm1 + args.rb*qr1) / (args.lb*ql2 + args.mb*qm2 + args.rb*qr2))
	else:
		df = math.log((args.lb*ql1 + args.rb*qr1) / (args.lb*ql2 + args.rb*qr2))
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
		if (args.filter == "human"):
			lvalues[i] = vpt(i+300, l1) * human_filter(i+300)
			mvalues[i] = vpt(i+300, m1) * human_filter(i+300)
			svalues[i] = vpt(i+300, s1) * human_filter(i+300)
			rvalues[i] = vpt(i+300, rod) * human_filter(i+300)
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
	
	plt.plot(xvalues, yvalues/ms, 'k')
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
	xvalues = np.empty(37)
	yvalues = np.empty(37)
	labels = np.empty(37)
	
	for i in range(0, 37):
		w = i*10 + 340
		labels[i] = w
		el += vpt(w, l1) * e[i]
		em += vpt(w, m1) * e[i]
		es += vpt(w, s1) * e[i]
		d65l += vpt(w, l1) * d65[i]
		d65m += vpt(w, m1) * d65[i]
		d65s += vpt(w, s1) * d65[i]
		l = vpt(w, l1) / (vpt(w, l1) + vpt(w, m1) + vpt(w, s1))
		m = vpt(w, m1) / (vpt(w, l1) + vpt(w, m1) + vpt(w, s1))
		s = vpt(w, s1) / (vpt(w, l1) + vpt(w, m1) + vpt(w, s1))
		xvalues[i] = math.sqrt(1/2)*(l - s)
		yvalues[i] = math.sqrt(2/3)*(m - (l + s)/2)
	plt.plot(xvalues, yvalues, '-k')
	for i in range(0, 17):
		plt.plot(xvalues[i*2], yvalues[i*2], 'ok')
		plt.text(xvalues[i*2], yvalues[i*2], labels[i*2])
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
	
# wavelength discrimination take three
# This is mostly right but for honeybees and butterflies gives an extra minimum
# where a maximum should be (really an extra peak in the first trough) and way too
# high values on either side of it. I cut the values off at 500 because otherwise
# the graph is unreadable.
# It looks much better if I include chromatic adaptation and use the "i" light. Also
# better if the white values are divided by 100?
if (args.wd):
	# background light
	wpl = 0
	wpm = 0
	wps = 0
	if (args.white == 'i'):
		start = 300
	else:
		start = 340
	
	for i in range(wp.shape[0]):
		w = i*10 + start
		wpl += vpt(w, l1) * wp[i]
		wpm += vpt(w, m1) * wp[i]
		wps += vpt(w, s1) * wp[i]
	kl = 1
	km = kl / wpm
	ks = kl / wps
	
	xvalues = np.empty(350)
	yvalues = np.empty(350)
	for i in range(350, 700):
		xvalues[i-350] = i
		
		# derivatives
		dl1 = vpt(i+0.01, l1) / (vpt(i+0.01, l1) + vpt(i+0.01, m1) + vpt(i+0.01, s1))
		dm1 = vpt(i+0.01, l1) / (vpt(i+0.01, l1) + vpt(i+0.01, m1) + vpt(i+0.01, s1))
		ds1 = vpt(i+0.01, s1) / (vpt(i+0.01, l1) + vpt(i+0.01, m1) + vpt(i+0.01, s1))
		dl2 = vpt(i-0.01, l1) / (vpt(i-0.01, l1) + vpt(i-0.01, m1) + vpt(i-0.01, s1))
		dm2 = vpt(i-0.01, l1) / (vpt(i-0.01, l1) + vpt(i-0.01, m1) + vpt(i-0.01, s1))
		ds2 = vpt(i-0.01, s1) / (vpt(i-0.01, l1) + vpt(i-0.01, m1) + vpt(i-0.01, s1))
		dll = (dl1 - dl2) / 0.02
		dlm = (dm1 - dm2) / 0.02
		dls = (ds1 - ds2) / 0.02
		dfl = (kl / (1 + kl*vpt(i, l1))) * dll
		dfm = (km / (1 + km*vpt(i, m1))) * dlm
		dfs = (ks / (1 + ks*vpt(i, s1))) * dls
		
		v = math.sqrt(((wl*wm)**2 + (wl*ws)**2 + (wm*ws)**2) / (wl**2*(dfm - dfs)**2 + wm**2*(dfl - dfs)**2 + ws**2*(dfm - dfl)**2))
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
# Note these are 300-900 nm rather than 340-830 nm.
if (args.kodak or args.kcheck):

	# optical density (-log transmission)
	red25 = np.array([
		2.423441558, # 300
		2.300623377,
		2.095532468,
		1.826415584,
		1.785675325,
		1.919655844,
		2.181538961,
		2.589188312,
		3.0,
		3.0,
		3.0, # 400
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0, # 500
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		1.813064935,
		0.525675325,
		0.12974026, # 600
		0.061480519,
		0.048883117,
		0.048883117,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385, # 700
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385, # 800
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385,
		0.0385 # 900
	])

	yellow15 = np.array([
		3.0, # 300
		2.012928025,
		1.90243502,
		2.381175769,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0, # 400
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0, # 500
		1.442115433,
		0.538040586,
		0.197090581,
		0.087322195,
		0.051651715,
		0.044886722,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986, # 600
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986, # 700
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986, # 800
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986,
		0.039887986 # 900
	])
	
	green58 = np.array([
		3.0, # 300
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0, # 400
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		2.099148329,
		1.44207606,
		0.903494563,
		0.505296555, # 500
		0.304288077,
		0.265710889,
		0.284540692,
		0.353630693,
		0.450726311,
		0.603871935,
		0.824304706,
		1.116218359,
		1.557238985,
		2.058303034, # 600
		2.586623128,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0, # 700
		2.194344166,
		1.485783951,
		0.889530526,
		0.575549687,
		0.367349181,
		0.239760136,
		0.16994641,
		0.121469681,
		0.10261403,
		0.085315682, # 800
		0.07716084,
		0.071661813,
		0.069529406,
		0.067338842,
		0.065083659,
		0.062027209,
		0.062027209,
		0.062027209,
		0.062027209,
		0.062027209 # 900
	])
	
	blue47 = np.array([
		2.806412992, #300
		2.941547894,
		2.72582778,
		2.540505092,
		2.601416432,
		2.502024554,
		2.174584034,
		1.826738344,
		1.532238615,
		1.211871628,
		0.899497151, # 400
		0.577505776,
		0.382462642,
		0.312497438,
		0.310284128,
		0.33718038,
		0.402880755,
		0.493788277,
		0.61670467,
		0.81348221,
		1.084716293, # 500
		1.532238615,
		2.067976245,
		2.782500178,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0, # 600
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		3.0,
		2.939224565,
		2.454930225, # 700
		2.015037999,
		1.644936243,
		1.272291122,
		0.941084091,
		0.625933268,
		0.393477421,
		0.256006247,
		0.15697031,
		0.118942788,
		0.088875418, # 800
		0.065046736,
		0.065046736,
		0.050944971,
		0.050944971,
		0.050944971,
		0.043864966,
		0.043864966,
		0.043864966,
		0.043864966,
		0.043864966 # 900
	])
	
	# neutral density filters
	nd01 = np.array([
		0.338821199, # 300
		0.296074593,
		0.250350004,
		0.229619894,
		0.214813592,
		0.194578741,
		0.171481681,
		0.159318902,
		0.149683871,
		0.144904947,
		0.139180529, # 400
		0.126323102,
		0.12428418,
		0.120341407,
		0.120341407,
		0.113311308,
		0.113311308,
		0.113311308,
		0.113311308,
		0.110474826,
		0.110474826, # 500
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759, # 600
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.107882759,
		0.103804915,
		0.103804915, # 700
		0.099592001,
		0.096472515,
		0.086831052,
		0.084386919,
		0.082309406,
		0.081621189,
		0.081621189,
		0.077298933,
		0.075253579,
		0.075253579, # 800
		0.075253579,
		0.075253579,
		0.075253579,
		0.075253579,
		0.073652028,
		0.073652028,
		0.073652028,
		0.073652028,
		0.073652028,
		0.073652028 # 900
	])
	
	nd02 = np.array([
		0.512531049, # 300
		0.468346895,
		0.39609636,
		0.360025696,
		0.33735546,
		0.30890364,
		0.288263383,
		0.265792291,
		0.262085653,
		0.254447537,
		0.247426124, # 400
		0.235850107,
		0.231963597,
		0.22562955,
		0.218222698,
		0.212762313,
		0.209428266,
		0.200126338,
		0.200126338,
		0.200126338,
		0.200126338, # 500
		0.197325482,
		0.197325482,
		0.197325482,
		0.197325482,
		0.197325482,
		0.197325482,
		0.197325482,
		0.197325482,
		0.197325482,
		0.197325482, # 600
		0.197325482,
		0.197325482,
		0.202066381,
		0.202066381,
		0.202066381,
		0.202066381,
		0.202066381,
		0.202066381,
		0.202066381,
		0.196194861, # 700
		0.186533191,
		0.172336188,
		0.163124197,
		0.150077088,
		0.14179015,
		0.13417773,
		0.13417773,
		0.13417773,
		0.130040685,
		0.130040685, # 800
		0.126777302,
		0.125852248,
		0.125852248,
		0.125852248,
		0.121316916,
		0.121316916,
		0.121316916,
		0.121316916,
		0.121316916,
		0.121316916 # 900
	])
	
	nd03 = np.array([
		0.854971744, # 300
		0.760632686,
		0.631502736,
		0.57084613,
		0.540373666,
		0.49574778,
		0.457375285,
		0.425115222,
		0.41742022,
		0.396359877,
		0.392637316, # 400
		0.377010247,
		0.363689757,
		0.355629546,
		0.346345568,
		0.337343505,
		0.332083225,
		0.32595157,
		0.32595157,
		0.317667109,
		0.317667109, # 500
		0.317667109,
		0.317667109,
		0.317667109,
		0.317667109,
		0.317667109,
		0.319204828,
		0.320345303,
		0.312272279,
		0.307524571,
		0.306268768, # 600
		0.309100733,
		0.312797666,
		0.323478406,
		0.32595157,
		0.32595157,
		0.32595157,
		0.32595157,
		0.32595157,
		0.32595157,
		0.315347717, # 700
		0.301418546,
		0.277737674,
		0.252845848,
		0.235014715,
		0.218529086,
		0.210834084,
		0.20336974,
		0.197148386,
		0.197148386,
		0.193701332, # 800
		0.192157206,
		0.192157206,
		0.184346875,
		0.184346875,
		0.184346875,
		0.179906712,
		0.179906712,
		0.178779051,
		0.178779051,
		0.178779051 # 900
	])
	
	nd04 = np.array([
		1.029570513, # 300
		0.926320513,
		0.800365385,
		0.727153846,
		0.692576923,
		0.633108974,
		0.584423077,
		0.550814103,
		0.535929487,
		0.518224359,
		0.502814103, # 400
		0.481301282,
		0.467788462,
		0.453423077,
		0.442448718,
		0.433044872,
		0.423653846,
		0.414519231,
		0.414519231,
		0.414519231,
		0.40699359, # 500
		0.40699359,
		0.40699359,
		0.40699359,
		0.40699359,
		0.40699359,
		0.40699359,
		0.40699359,
		0.403378205,
		0.394410256,
		0.39400641, # 600
		0.397480769,
		0.403339744,
		0.40699359,
		0.414371795,
		0.414371795,
		0.414371795,
		0.414371795,
		0.414371795,
		0.414371795,
		0.40699359, # 700
		0.384628205,
		0.360371795,
		0.328782051,
		0.304942308,
		0.304942308,
		0.269974359,
		0.2575,
		0.252064103,
		0.247923077,
		0.241929487, # 800
		0.239083333,
		0.239083333,
		0.235544872,
		0.230679487,
		0.230679487,
		0.226410256,
		0.223673077,
		0.223673077,
		0.221519231,
		0.221519231 # 900
	])
	
	nd05 = np.array([
		1.135265525, # 300
		1.045072805,
		0.963321199,
		0.900211991,
		0.855880086,
		0.791775161,
		0.724541756,
		0.688085653,
		0.658124197,
		0.645398287,
		0.621770878, # 400
		0.601650964,
		0.578762313,
		0.563171306,
		0.543481799,
		0.532817987,
		0.524762313,
		0.520985011,
		0.50933833,
		0.506126338,
		0.501995717, # 500
		0.501211991,
		0.499618844,
		0.497640257,
		0.493773019,
		0.493773019,
		0.493773019,
		0.493773019,
		0.488987152,
		0.477520343,
		0.475252677, # 600
		0.48,
		0.488139186,
		0.497640257,
		0.499618844,
		0.499618844,
		0.499618844,
		0.499618844,
		0.499618844,
		0.495436831,
		0.485678801, # 700
		0.461762313,
		0.429359743,
		0.399199143,
		0.371948608,
		0.352310493,
		0.334856531,
		0.324880086,
		0.317897216,
		0.310837259,
		0.305248394, # 800
		0.303745182,
		0.29832334,
		0.29832334,
		0.293087794,
		0.288250535,
		0.284479657,
		0.283188437,
		0.277477516,
		0.277477516,
		0.277477516 # 900
	])
	
	nd06 = np.array([
		1.323256959, # 300
		1.24437045,
		1.112158458,
		1.052364026,
		1.01943469,
		0.930738758,
		0.849501071,
		0.811066381,
		0.785858672,
		0.760342612,
		0.735012848, # 400
		0.707254818,
		0.680884368,
		0.663353319,
		0.643246253,
		0.624706638,
		0.620100642,
		0.607291221,
		0.60003212,
		0.60003212,
		0.60003212, # 500
		0.60003212,
		0.60003212,
		0.60003212,
		0.596119914,
		0.596036403,
		0.597674518,
		0.60003212,
		0.592265525,
		0.577554604,
		0.570970021, # 600
		0.576796574,
		0.591115632,
		0.60003212,
		0.609366167,
		0.609366167,
		0.609366167,
		0.609366167,
		0.609366167,
		0.609366167,
		0.60003212, # 700
		0.573738758,
		0.532130621,
		0.486558887,
		0.444276231,
		0.412862955,
		0.392421842,
		0.375372591,
		0.363835118,
		0.352856531,
		0.352856531, # 800
		0.347852248,
		0.340734475,
		0.340734475,
		0.330089936,
		0.330089936,
		0.325117773,
		0.325117773,
		0.31866167,
		0.31866167,
		0.31866167 # 900
	])
	
	nd07 = np.array([
		1.543503212, #300
		1.451524625,
		1.323366167,
		1.241633833,
		1.193717345,
		1.090805139,
		1.010550321,
		0.957025696,
		0.925586724,
		0.899152034,
		0.87072591, # 400
		0.839049251,
		0.805477516,
		0.780366167,
		0.759982869,
		0.747122056,
		0.730226981,
		0.719203426,
		0.705571734,
		0.701961456,
		0.695447537, # 500
		0.695447537,
		0.695447537,
		0.695447537,
		0.687057816,
		0.688278373,
		0.689184154,
		0.695447537,
		0.683293362,
		0.666122056,
		0.660648822, # 600
		0.669186296,
		0.679316916,
		0.690070664,
		0.695447537,
		0.695447537,
		0.695447537,
		0.695447537,
		0.695447537,
		0.695447537,
		0.688593148, # 700
		0.655336188,
		0.614961456,
		0.556265525,
		0.517053533,
		0.473152034,
		0.454246253,
		0.438289079,
		0.423411135,
		0.417147752,
		0.409561028, # 800
		0.404595289,
		0.401164882,
		0.39382227,
		0.391252677,
		0.388271949,
		0.3808394,
		0.37882227,
		0.376593148,
		0.371094218,
		0.367901499 # 900
	])
	
	nd08 = np.array([
		1.733531049, # 300
		1.63137045,
		1.490531049,
		1.411541756,
		1.355505353,
		1.235119914,
		1.145216274,
		1.091762313,
		1.048374732,
		1.02243469,
		0.98766167, # 400
		0.953184154,
		0.918231263,
		0.890357602,
		0.864443255,
		0.845788009,
		0.832438972,
		0.819957173,
		0.8071606,
		0.802959315,
		0.794511777, # 500
		0.793284797,
		0.788980728,
		0.788980728,
		0.777841542,
		0.774044968,
		0.777841542,
		0.783680942,
		0.774237687,
		0.754149893,
		0.748156317, # 600
		0.756526767,
		0.767029979,
		0.780327623,
		0.791922912,
		0.791922912,
		0.791922912,
		0.785498929,
		0.791922912,
		0.787149893,
		0.777025696, # 700
		0.739586724,
		0.690353319,
		0.631766595,
		0.583124197,
		0.543809422,
		0.515890792,
		0.4934197,
		0.478111349,
		0.468025696,
		0.462610278, # 800
		0.454940043,
		0.449062099,
		0.443620985,
		0.438077088,
		0.435648822,
		0.429147752,
		0.42369379,
		0.418432548,
		0.418432548,
		0.418432548 # 900
	])
	
	nd09 = np.array([
		1.876477516, # 300
		1.76406424,
		1.621329764,
		1.535800857,
		1.460492505,
		1.34320985,
		1.242269807,
		1.171824411,
		1.129361884,
		1.094267666,
		1.060477516, # 400
		1.019139186,
		0.988773019,
		0.959029979,
		0.932331906,
		0.914678801,
		0.896865096,
		0.887453961,
		0.876655246,
		0.864635974,
		0.857293362, # 500
		0.852321199,
		0.852321199,
		0.850079229,
		0.845023555,
		0.839473233,
		0.838670236,
		0.838021413,
		0.834199143,
		0.815942184,
		0.814085653, # 600
		0.819404711,
		0.829824411,
		0.844143469,
		0.852321199,
		0.852321199,
		0.852321199,
		0.848511777,
		0.852321199,
		0.852321199,
		0.835059957, # 700
		0.801295503,
		0.752672377,
		0.687449679,
		0.632319058,
		0.593569593,
		0.568162741,
		0.545788009,
		0.529477516,
		0.519655246,
		0.51079015, # 800
		0.504089936,
		0.496740899,
		0.492256959,
		0.489873662,
		0.484888651,
		0.479794433,
		0.472554604,
		0.47008137,
		0.465494647,
		0.465494647 # 900
	])
	
	nd10 = np.array([
		2.158094421, # 300
		2.02001073,
		1.860347639,
		1.766027897,
		1.68323176,
		1.554251073,
		1.424053648,
		1.355980687,
		1.311186695,
		1.269270386,
		1.232182403, # 400
		1.186332618,
		1.150751073,
		1.118201717,
		1.092251073,
		1.069648069,
		1.048821888,
		1.035334764,
		1.024345494,
		1.013658798,
		1.00543133, # 500
		1.00543133,
		1.00543133,
		1.00543133,
		0.999199571,
		0.999199571,
		0.999199571,
		0.999199571,
		0.985667382,
		0.969347639,
		0.961075107, # 600
		0.969412017,
		0.98622103,
		0.999199571,
		1.008012876,
		1.016349785,
		1.008012876,
		1.008012876,
		1.008012876,
		1.008012876,
		0.986774678, # 700
		0.948154506,
		0.884697425,
		0.816296137,
		0.757435622,
		0.712487124,
		0.6777103,
		0.655281116,
		0.640944206,
		0.629040773,
		0.620652361, # 800
		0.611948498,
		0.603263948,
		0.600251073,
		0.595545064,
		0.584961373,
		0.579225322,
		0.574680258,
		0.567914163,
		0.564830472,
		0.560774678 # 900
	])
	
	nd20 = np.array([
		3.0, # 300
		3.0,
		2.862978678,
		2.791944563,
		2.647950959,
		2.418345416,
		2.243904051,
		2.156552239,
		2.064997868,
		2.000859275,
		1.938486141, # 400
		1.858381663,
		1.784085288,
		1.739296375,
		1.687420043,
		1.654081023,
		1.623345416,
		1.593639659,
		1.57201919,
		1.554255864,
		1.541648188, # 500
		1.529712154,
		1.529712154,
		1.526603412,
		1.526603412,
		1.500850746,
		1.510228145,
		1.510228145,
		1.488102345,
		1.448168443,
		1.42961194, # 600
		1.436123667,
		1.455249467,
		1.481014925,
		1.504816631,
		1.504816631,
		1.493034115,
		1.490008529,
		1.490008529,
		1.490008529,
		1.473486141, # 700
		1.414791045,
		1.323473348,
		1.210765458,
		1.12065032,
		1.04263113,
		0.980597015,
		0.941353945,
		0.91608742,
		0.895823028,
		0.880324094, # 800
		0.865970149,
		0.855428571,
		0.84817484,
		0.837108742,
		0.82473774,
		0.816882729,
		0.808157783,
		0.799995736,
		0.789978678,
		0.783383795 # 900
	])
	
	if (args.kcheck):
		# plot test
		xvalues = np.empty(61)
		for i in range(0, 61):
			xvalues[i] = i*10 + 300
		plt.plot(xvalues, red25, 'o-r', mec='k', label="red 25")
		plt.plot(xvalues, yellow15, 's-y', mec='k', label="yellow 15")
		plt.plot(xvalues, green58, '^-g', mec='k', label="green 58")
		plt.plot(xvalues, blue47, 'D-b', mec='k', label="blue 47")
		plt.xlabel("Wavelength (nm)")
		plt.ylabel("Optical density")
		plt.legend()
		plt.show()
		
		# transmission
		red25t = np.empty(61)
		for i in range(0, 61):
			red25t[i] = math.exp(-red25[i]) * 100
		yellow15t = np.empty(61)
		for i in range(0, 61):
			yellow15t[i] = math.exp(-yellow15[i]) * 100
		green58t = np.empty(61)
		for i in range(0, 61):
			green58t[i] = math.exp(-green58[i]) * 100
		blue47t = np.empty(61)
		for i in range(0, 61):
			blue47t[i] = math.exp(-blue47[i]) * 100
		plt.plot(xvalues, red25t, 'o-r', mec='k', label="red 25")
		plt.plot(xvalues, yellow15t, 's-y', mec='k', label="yellow 15")
		plt.plot(xvalues, green58t, '^-g', mec='k', label="green 58")
		plt.plot(xvalues, blue47t, 'D-b', mec='k', label="blue 47")
		plt.xlabel("Wavelength (nm)")
		plt.ylabel("Transmission (%)")
		plt.legend()
		plt.show()
		
		# red brightness tests
		red25_01 = np.empty(61)
		for i in range(0, 61):
			red25_01[i] = red25[i] + nd01[i]
		red25_01t = np.empty(61)
		for i in range(0, 61):
			red25_01t[i] = math.exp(-red25_01[i])
		red25_02 = np.empty(61)
		for i in range(0, 61):
			red25_02[i] = red25[i] + nd02[i]
		red25_02t = np.empty(61)
		for i in range(0, 61):
			red25_02t[i] = math.exp(-red25_02[i])
		red25_03 = np.empty(61)
		for i in range(0, 61):
			red25_03[i] = red25[i] + nd03[i]
		red25_03t = np.empty(61)
		for i in range(0, 61):
			red25_03t[i] = math.exp(-red25_03[i])
		red25_04 = np.empty(61)
		for i in range(0, 61):
			red25_04[i] = red25[i] + nd04[i]
		red25_04t = np.empty(61)
		for i in range(0, 61):
			red25_04t[i] = math.exp(-red25_04[i])
		red25_05 = np.empty(61)
		for i in range(0, 61):
			red25_05[i] = red25[i] + nd05[i]
		red25_05t = np.empty(61)
		for i in range(0, 61):
			red25_05t[i] = math.exp(-red25_05[i])
		red25_06 = np.empty(61)
		for i in range(0, 61):
			red25_06[i] = red25[i] + nd06[i]
		red25_06t = np.empty(61)
		for i in range(0, 61):
			red25_06t[i] = math.exp(-red25_06[i])
		red25_07 = np.empty(61)
		for i in range(0, 61):
			red25_07[i] = red25[i] + nd07[i]
		red25_07t = np.empty(61)
		for i in range(0, 61):
			red25_07t[i] = math.exp(-red25_07[i])
		red25_08 = np.empty(61)
		for i in range(0, 61):
			red25_08[i] = red25[i] + nd08[i]
		red25_08t = np.empty(61)
		for i in range(0, 61):
			red25_08t[i] = math.exp(-red25_08[i])
		red25_09 = np.empty(61)
		for i in range(0, 61):
			red25_09[i] = red25[i] + nd09[i]
		red25_09t = np.empty(61)
		for i in range(0, 61):
			red25_09t[i] = math.exp(-red25_09[i])
		red25_10 = np.empty(61)
		for i in range(0, 61):
			red25_10[i] = red25[i] + nd10[i]
		red25_10t = np.empty(61)
		for i in range(0, 61):
			red25_10t[i] = math.exp(-red25_10[i])
		red25_1001 = np.empty(61)
		for i in range(0, 61):
			red25_1001[i] = red25[i] + nd10[i] + nd01[i]
		red25_1001t = np.empty(61)
		for i in range(0, 61):
			red25_1001t[i] = math.exp(-red25_1001[i])
		red25_1002 = np.empty(61)
		for i in range(0, 61):
			red25_1002[i] = red25[i] + nd10[i] + nd02[i]
		red25_1002t = np.empty(61)
		for i in range(0, 61):
			red25_1002t[i] = math.exp(-red25_1002[i])
		red25_1003 = np.empty(61)
		for i in range(0, 61):
			red25_1003[i] = red25[i] + nd10[i] + nd03[i]
		red25_1003t = np.empty(61)
		for i in range(0, 61):
			red25_1003t[i] = math.exp(-red25_1003[i])
		red25_1004 = np.empty(61)
		for i in range(0, 61):
			red25_1004[i] = red25[i] + nd10[i] + nd04[i]
		red25_1004t = np.empty(61)
		for i in range(0, 61):
			red25_1004t[i] = math.exp(-red25_1004[i])
		red25_20 = np.empty(61)
		for i in range(0, 61):
			red25_20[i] = red25[i] + nd20[i]
		red25_20t = np.empty(61)
		for i in range(0, 61):
			red25_20t[i] = math.exp(-red25_20[i])
		
		# yellow brightness tests
		yellow15_01 = np.empty(61)
		for i in range(0, 61):
			yellow15_01[i] = yellow15[i] + nd01[i]
		yellow15_01t = np.empty(61)
		for i in range(0, 61):
			yellow15_01t[i] = math.exp(-yellow15_01[i])
		yellow15_02 = np.empty(61)
		for i in range(0, 61):
			yellow15_02[i] = yellow15[i] + nd02[i]
		yellow15_02t = np.empty(61)
		for i in range(0, 61):
			yellow15_02t[i] = math.exp(-yellow15_02[i])
		yellow15_03 = np.empty(61)
		for i in range(0, 61):
			yellow15_03[i] = yellow15[i] + nd03[i]
		yellow15_03t = np.empty(61)
		for i in range(0, 61):
			yellow15_03t[i] = math.exp(-yellow15_03[i])
		yellow15_04 = np.empty(61)
		for i in range(0, 61):
			yellow15_04[i] = yellow15[i] + nd04[i]
		yellow15_04t = np.empty(61)
		for i in range(0, 61):
			yellow15_04t[i] = math.exp(-yellow15_04[i])
		yellow15_05 = np.empty(61)
		for i in range(0, 61):
			yellow15_05[i] = yellow15[i] + nd05[i]
		yellow15_05t = np.empty(61)
		for i in range(0, 61):
			yellow15_05t[i] = math.exp(-yellow15_05[i])
		yellow15_06 = np.empty(61)
		for i in range(0, 61):
			yellow15_06[i] = yellow15[i] + nd06[i]
		yellow15_06t = np.empty(61)
		for i in range(0, 61):
			yellow15_06t[i] = math.exp(-yellow15_06[i])
		yellow15_07 = np.empty(61)
		for i in range(0, 61):
			yellow15_07[i] = yellow15[i] + nd07[i]
		yellow15_07t = np.empty(61)
		for i in range(0, 61):
			yellow15_07t[i] = math.exp(-yellow15_07[i])
		yellow15_08 = np.empty(61)
		for i in range(0, 61):
			yellow15_08[i] = yellow15[i] + nd08[i]
		yellow15_08t = np.empty(61)
		for i in range(0, 61):
			yellow15_08t[i] = math.exp(-yellow15_08[i])
		yellow15_09 = np.empty(61)
		for i in range(0, 61):
			yellow15_09[i] = yellow15[i] + nd09[i]
		yellow15_09t = np.empty(61)
		for i in range(0, 61):
			yellow15_09t[i] = math.exp(-yellow15_09[i])
		yellow15_10 = np.empty(61)
		for i in range(0, 61):
			yellow15_10[i] = yellow15[i] + nd10[i]
		yellow15_10t = np.empty(61)
		for i in range(0, 61):
			yellow15_10t[i] = math.exp(-yellow15_10[i])
		yellow15_1001 = np.empty(61)
		for i in range(0, 61):
			yellow15_1001[i] = yellow15[i] + nd10[i] + nd01[i]
		yellow15_1001t = np.empty(61)
		for i in range(0, 61):
			yellow15_1001t[i] = math.exp(-yellow15_1001[i])
		yellow15_1002 = np.empty(61)
		for i in range(0, 61):
			yellow15_1002[i] = yellow15[i] + nd10[i] + nd02[i]
		yellow15_1002t = np.empty(61)
		for i in range(0, 61):
			yellow15_1002t[i] = math.exp(-yellow15_1002[i])
		yellow15_1010 = np.empty(61)
		for i in range(0, 61):
			yellow15_1010[i] = yellow15[i] + nd10[i] + nd10[i]
		yellow15_1010t = np.empty(61)
		for i in range(0, 61):
			yellow15_1010t[i] = math.exp(-yellow15_1010[i])
		yellow15_20 = np.empty(61)
		for i in range(0, 61):
			yellow15_20[i] = yellow15[i] + nd20[i]
		yellow15_20t = np.empty(61)
		for i in range(0, 61):
			yellow15_20t[i] = math.exp(-yellow15_20[i])
		yellow15_2001 = np.empty(61)
		for i in range(0, 61):
			yellow15_2001[i] = yellow15[i] + nd20[i] + nd01[i]
		yellow15_2001t = np.empty(61)
		for i in range(0, 61):
			yellow15_2001t[i] = math.exp(-yellow15_2001[i])
		yellow15_2002 = np.empty(61)
		for i in range(0, 61):
			yellow15_2002[i] = yellow15[i] + nd20[i] + nd02[i]
		yellow15_2002t = np.empty(61)
		for i in range(0, 61):
			yellow15_2002t[i] = math.exp(-yellow15_2002[i])
		yellow15_2003 = np.empty(61)
		for i in range(0, 61):
			yellow15_2003[i] = yellow15[i] + nd20[i] + nd03[i]
		yellow15_2003t = np.empty(61)
		for i in range(0, 61):
			yellow15_2003t[i] = math.exp(-yellow15_2003[i])
		yellow15_2004 = np.empty(61)
		for i in range(0, 61):
			yellow15_2004[i] = yellow15[i] + nd20[i] + nd04[i]
		yellow15_2004t = np.empty(61)
		for i in range(0, 61):
			yellow15_2004t[i] = math.exp(-yellow15_2004[i])
		yellow15_2005 = np.empty(61)
		for i in range(0, 61):
			yellow15_2005[i] = yellow15[i] + nd20[i] + nd05[i]
		yellow15_2005t = np.empty(61)
		for i in range(0, 61):
			yellow15_2005t[i] = math.exp(-yellow15_2005[i])
		yellow15_2006 = np.empty(61)
		for i in range(0, 61):
			yellow15_2006[i] = yellow15[i] + nd20[i] + nd06[i]
		yellow15_2006t = np.empty(61)
		for i in range(0, 61):
			yellow15_2006t[i] = math.exp(-yellow15_2006[i])
		yellow15_2007 = np.empty(61)
		for i in range(0, 61):
			yellow15_2007[i] = yellow15[i] + nd20[i] + nd07[i]
		yellow15_2007t = np.empty(61)
		for i in range(0, 61):
			yellow15_2007t[i] = math.exp(-yellow15_2007[i])
		
		# green brightness tests
		green58_01 = np.empty(61)
		for i in range(0, 61):
			green58_01[i] = green58[i] + nd01[i]
		green58_01t = np.empty(61)
		for i in range(0, 61):
			green58_01t[i] = math.exp(-green58_01[i])
		green58_02 = np.empty(61)
		for i in range(0, 61):
			green58_02[i] = green58[i] + nd02[i]
		green58_02t = np.empty(61)
		for i in range(0, 61):
			green58_02t[i] = math.exp(-green58_02[i])
		green58_03 = np.empty(61)
		for i in range(0, 61):
			green58_03[i] = green58[i] + nd03[i]
		green58_03t = np.empty(61)
		for i in range(0, 61):
			green58_03t[i] = math.exp(-green58_03[i])
		green58_04 = np.empty(61)
		for i in range(0, 61):
			green58_04[i] = green58[i] + nd04[i]
		green58_04t = np.empty(61)
		for i in range(0, 61):
			green58_04t[i] = math.exp(-green58_04[i])
		green58_05 = np.empty(61)
		for i in range(0, 61):
			green58_05[i] = green58[i] + nd05[i]
		green58_05t = np.empty(61)
		for i in range(0, 61):
			green58_05t[i] = math.exp(-green58_05[i])
		green58_06 = np.empty(61)
		for i in range(0, 61):
			green58_06[i] = green58[i] + nd06[i]
		green58_06t = np.empty(61)
		for i in range(0, 61):
			green58_06t[i] = math.exp(-green58_06[i])
		green58_07 = np.empty(61)
		for i in range(0, 61):
			green58_07[i] = green58[i] + nd07[i]
		green58_07t = np.empty(61)
		for i in range(0, 61):
			green58_07t[i] = math.exp(-green58_07[i])
		green58_08 = np.empty(61)
		for i in range(0, 61):
			green58_08[i] = green58[i] + nd08[i]
		green58_08t = np.empty(61)
		for i in range(0, 61):
			green58_08t[i] = math.exp(-green58_08[i])
		green58_09 = np.empty(61)
		for i in range(0, 61):
			green58_09[i] = green58[i] + nd09[i]
		green58_09t = np.empty(61)
		for i in range(0, 61):
			green58_09t[i] = math.exp(-green58_09[i])
		green58_10 = np.empty(61)
		for i in range(0, 61):
			green58_10[i] = green58[i] + nd10[i]
		green58_10t = np.empty(61)
		for i in range(0, 61):
			green58_10t[i] = math.exp(-green58_10[i])
		green58_1001 = np.empty(61)
		for i in range(0, 61):
			green58_1001[i] = green58[i] + nd10[i] + nd01[i]
		green58_1001t = np.empty(61)
		for i in range(0, 61):
			green58_1001t[i] = math.exp(-green58_1001[i])
		green58_1002 = np.empty(61)
		for i in range(0, 61):
			green58_1002[i] = green58[i] + nd10[i] + nd02[i]
		green58_1002t = np.empty(61)
		for i in range(0, 61):
			green58_1002t[i] = math.exp(-green58_1002[i])
		green58_1003 = np.empty(61)
		for i in range(0, 61):
			green58_1003[i] = green58[i] + nd10[i] + nd03[i]
		green58_1003t = np.empty(61)
		for i in range(0, 61):
			green58_1003t[i] = math.exp(-green58_1003[i])
		green58_1004 = np.empty(61)
		for i in range(0, 61):
			green58_1004[i] = green58[i] + nd10[i] + nd04[i]
		green58_1004t = np.empty(61)
		for i in range(0, 61):
			green58_1004t[i] = math.exp(-green58_1004[i])
		green58_20 = np.empty(61)
		for i in range(0, 61):
			green58_20[i] = green58[i] + nd20[i]
		green58_20t = np.empty(61)
		for i in range(0, 61):
			green58_20t[i] = math.exp(-green58_20[i])
		
		# blue
		blue47t = np.empty(61)
		for i in range(0, 61):
			blue47t[i] = math.exp(-blue47[i])
		
		print("Brightness difference between blue and red:")
		print("0.1: " + str(brightness_contrast(red25_01t, blue47t)))
		print("0.2: " + str(brightness_contrast(red25_02t, blue47t)))
		print("0.3: " + str(brightness_contrast(red25_03t, blue47t)))
		print("0.4: " + str(brightness_contrast(red25_04t, blue47t)))
		print("0.5: " + str(brightness_contrast(red25_05t, blue47t)))
		print("0.6: " + str(brightness_contrast(red25_06t, blue47t)))
		print("0.7: " + str(brightness_contrast(red25_07t, blue47t)))
		print("0.8: " + str(brightness_contrast(red25_08t, blue47t)))
		print("0.9: " + str(brightness_contrast(red25_09t, blue47t)))
		print("1.0: " + str(brightness_contrast(red25_10t, blue47t)))
		print("1.0 + 0.1: " + str(brightness_contrast(red25_1001t, blue47t)))
		print("1.0 + 0.2: " + str(brightness_contrast(red25_1002t, blue47t)))
		print("1.0 + 0.3: " + str(brightness_contrast(red25_1003t, blue47t)))
		print("1.0 + 0.4: " + str(brightness_contrast(red25_1004t, blue47t)))
		print("2.0: " + str(brightness_contrast(red25_20t, blue47t)))
		print("Brightness difference between blue and yellow:")
		print("0.1: " + str(brightness_contrast(yellow15_01t, blue47t)))
		print("0.2: " + str(brightness_contrast(yellow15_02t, blue47t)))
		print("0.3: " + str(brightness_contrast(yellow15_03t, blue47t)))
		print("0.4: " + str(brightness_contrast(yellow15_04t, blue47t)))
		print("0.5: " + str(brightness_contrast(yellow15_05t, blue47t)))
		print("0.6: " + str(brightness_contrast(yellow15_06t, blue47t)))
		print("0.7: " + str(brightness_contrast(yellow15_07t, blue47t)))
		print("0.8: " + str(brightness_contrast(yellow15_08t, blue47t)))
		print("0.9: " + str(brightness_contrast(yellow15_09t, blue47t)))
		print("1.0: " + str(brightness_contrast(yellow15_10t, blue47t)))
		print("1.0 + 0.1: " + str(brightness_contrast(yellow15_1001t, blue47t)))
		print("1.0 + 0.2: " + str(brightness_contrast(yellow15_1002t, blue47t)))
		print("1.0 + 1.0: " + str(brightness_contrast(yellow15_1010t, blue47t)))
		print("2.0: " + str(brightness_contrast(yellow15_20t, blue47t)))
		print("2.0 + 0.1: " + str(brightness_contrast(yellow15_2001t, blue47t)))
		print("2.0 + 0.2: " + str(brightness_contrast(yellow15_2002t, blue47t)))
		print("2.0 + 0.3: " + str(brightness_contrast(yellow15_2003t, blue47t)))
		print("2.0 + 0.4: " + str(brightness_contrast(yellow15_2004t, blue47t)))
		print("2.0 + 0.5: " + str(brightness_contrast(yellow15_2005t, blue47t)))
		print("2.0 + 0.6: " + str(brightness_contrast(yellow15_2006t, blue47t)))
		print("2.0 + 0.7: " + str(brightness_contrast(yellow15_2007t, blue47t)))
		print("Brightness difference between blue and green:")
		print("0.1: " + str(brightness_contrast(green58_01t, blue47t)))
		print("0.2: " + str(brightness_contrast(green58_02t, blue47t)))
		print("0.3: " + str(brightness_contrast(green58_03t, blue47t)))
		print("0.4: " + str(brightness_contrast(green58_04t, blue47t)))
		print("0.5: " + str(brightness_contrast(green58_05t, blue47t)))
		print("0.6: " + str(brightness_contrast(green58_06t, blue47t)))
		print("0.7: " + str(brightness_contrast(green58_07t, blue47t)))
		print("0.8: " + str(brightness_contrast(green58_08t, blue47t)))
		print("0.9: " + str(brightness_contrast(green58_09t, blue47t)))
		print("1.0: " + str(brightness_contrast(green58_10t, blue47t)))
		print("1.0 + 0.1: " + str(brightness_contrast(green58_1001t, blue47t)))
		print("1.0 + 0.2: " + str(brightness_contrast(green58_1002t, blue47t)))
		print("1.0 + 0.3: " + str(brightness_contrast(green58_1003t, blue47t)))
		print("1.0 + 0.4: " + str(brightness_contrast(green58_1004t, blue47t)))
		print("2.0: " + str(brightness_contrast(green58_20t, blue47t)))
	
	if (args.kodak):
		# red 25
		red25e = np.empty(61)
		for i in range(0, 61):
			red25e[i] = red25[i] + nd10[i] + nd01[i]
		red25e_0 = np.empty(61)
		for i in range(0, 61):
			red25e_0[i] = math.exp(-red25e[i])
		red25e_03 = np.empty(61)
		for i in range(0, 61):
			red25e_03[i] = math.exp(-(red25e[i] + nd03[i]))
		red25e_07 = np.empty(61)
		for i in range(0, 61):
			red25e_07[i] = math.exp(-(red25e[i] + nd07[i]))
		red25e_10 = np.empty(61)
		for i in range(0, 61):
			red25e_10[i] = math.exp(-(red25e[i] + nd10[i]))
		
		# yellow 15
		yellow15e = np.empty(61)
		for i in range(0, 61):
			yellow15e[i] = yellow15[i] + nd20[i] + nd05[i]
		yellow15e_0 = np.empty(61)
		for i in range(0, 61):
			yellow15e_0[i] = math.exp(-yellow15e[i])
		yellow15e_03 = np.empty(61)
		for i in range(0, 61):
			yellow15e_03[i] = math.exp(-(yellow15e[i] + nd03[i]))
		yellow15e_07 = np.empty(61)
		for i in range(0, 61):
			yellow15e_07[i] = math.exp(-(yellow15e[i] + nd07[i]))
		yellow15e_10 = np.empty(61)
		for i in range(0, 61):
			yellow15e_10[i] = math.exp(-(yellow15e[i] + nd10[i]))
		
		# green 58
		green58e = np.empty(61)
		for i in range(0, 61):
			green58e[i] = green58[i] + nd10[i] + nd03[i]
		green58e_0 = np.empty(61)
		for i in range(0, 61):
			green58e_0[i] = math.exp(-green58e[i])
		green58e_03 = np.empty(61)
		for i in range(0, 61):
			green58e_03[i] = math.exp(-(green58e[i] + nd03[i]))
		green58e_07 = np.empty(61)
		for i in range(0, 61):
			green58e_07[i] = math.exp(-(green58e[i] + nd07[i]))
		green58e_10 = np.empty(61)
		for i in range(0, 61):
			green58e_10[i] = math.exp(-(green58e[i] + nd10[i]))
		
		# blue 47
		blue47_0 = np.empty(61)
		for i in range(0, 61):
			blue47_0[i] = math.exp(-blue47[i])
		blue47_03 = np.empty(61)
		for i in range(0, 61):
			blue47_03[i] = math.exp(-(blue47[i] + nd03[i]))
		blue47_07 = np.empty(61)
		for i in range(0, 61):
			blue47_07[i] = math.exp(-(blue47[i] + nd07[i]))
		blue47_10 = np.empty(61)
		for i in range(0, 61):
			blue47_10[i] = math.exp(-(blue47[i] + nd10[i]))
		
		# brightness levels and color space coordinates
		red25l = np.empty(4)
		red25cs = np.empty(4)
		red25cs1 = np.empty(4)
		print("Red 25 (+ 1.0 + 0.2)")
		red25data = spectral_rendering(red25e_0)
		red25l[0] = red25data[0]
		red25cs[0] = red25data[1]
		red25cs1[0] = red25data[2]
		print("+ 0.3")
		red25data = spectral_rendering(red25e_03)
		red25l[1] = red25data[0]
		red25cs[1] = red25data[1]
		red25cs1[1] = red25data[2]
		print("+ 0.7")
		red25data = spectral_rendering(red25e_07)
		red25l[2] = red25data[0]
		red25cs[2] = red25data[1]
		red25cs1[2] = red25data[2]
		print("+ 1.0")
		red25data = spectral_rendering(red25e_10)
		red25l[3] = red25data[0]
		red25cs[3] = red25data[1]
		red25cs1[3] = red25data[2]
		
		yellow15l = np.empty(4)
		yellow15cs = np.empty(4)
		yellow15cs1 = np.empty(4)
		print("Yellow 15 (+ 2.0 + 0.6)")
		yellow15data = spectral_rendering(yellow15e_0)
		yellow15l[0] = yellow15data[0]
		yellow15cs[0] = yellow15data[1]
		yellow15cs1[0] = yellow15data[2]
		print("+ 0.3")
		yellow15data = spectral_rendering(yellow15e_03)
		yellow15l[1] = yellow15data[0]
		yellow15cs[1] = yellow15data[1]
		yellow15cs1[1] = yellow15data[2]
		print("+ 0.7")
		yellow15data = spectral_rendering(yellow15e_07)
		yellow15l[2] = yellow15data[0]
		yellow15cs[2] = yellow15data[1]
		yellow15cs1[2] = yellow15data[2]
		print("+ 1.0")
		yellow15data = spectral_rendering(yellow15e_10)
		yellow15l[3] = yellow15data[0]
		yellow15cs[3] = yellow15data[1]
		yellow15cs1[3] = yellow15data[2]
		
		green58l = np.empty(4)
		green58cs = np.empty(4)
		green58cs1 = np.empty(4)
		print("Green 58 (+ 1.0 + 0.3)")
		green58data = spectral_rendering(green58e_0)
		green58l[0] = green58data[0]
		green58cs[0] = green58data[1]
		green58cs1[0] = green58data[2]
		print("+ 0.3")
		green58data = spectral_rendering(green58e_03)
		green58l[1] = green58data[0]
		green58cs[1] = green58data[1]
		green58cs1[1] = green58data[2]
		print("+ 0.7")
		green58data = spectral_rendering(green58e_07)
		green58l[2] = green58data[0]
		green58cs[2] = green58data[1]
		green58cs1[2] = green58data[2]
		print("+ 1.0")
		green58data = spectral_rendering(green58e_10)
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
		
		# brightness differences between colors
		# We want to know both the direction and the strength of the difference
		# to judge whether it could be reliably used as a cue.
		
		# red-yellow
		ry = 0 # red brighter
		yr = 0 # yellow brighter
		equal = 0
		for i in range(0, 4):
			if (i == 0):
				rlevel = red25e
			elif (i == 1):
				rlevel = red25e_03
			elif (i == 2):
				rlevel = red25e_07
			elif (i == 3):
				rlevel = red25e_10
			for j in range(0, 4):
				if (j == 0):
					ylevel = yellow15e
				elif (j == 1):
					ylevel = yellow15e_03
				elif (j == 2):
					ylevel = yellow15e_07
				elif (j == 3):
					ylevel = yellow15e_10
				diff = red25l[i] - yellow15l[j]
				contrast = brightness_contrast(rlevel, ylevel)
				if (contrast >= 1):
					if (diff > 0):
						ry += 1
					elif (diff < 0):
						yr += 1
				else:
					equal += 1
		print("R-Y brightness differences:")
		print("Red brighter: " + str(ry))
		print("Yellow brighter: " + str(yr))
		print("Equal: " + str(equal))
		
		# red-green
		rg = 0 # red brighter
		gr = 0 # green brighter
		equal = 0
		for i in range(0, 4):
			if (i == 0):
				rlevel = red25e
			elif (i == 1):
				rlevel = red25e_03
			elif (i == 2):
				rlevel = red25e_07
			elif (i == 3):
				rlevel = red25e_10
			for j in range(0, 4):
				if (j == 0):
					glevel = green58e
				elif (j == 1):
					glevel = green58e_03
				elif (j == 2):
					glevel = green58e_07
				elif (j == 3):
					glevel = green58e_10
				diff = red25l[i] - green58l[j]
				contrast = brightness_contrast(rlevel, glevel)
				if (contrast >= 1):
					if (diff > 0):
						rg += 1
					elif (diff < 0):
						gr += 1
				else:
					equal += 1
		print("R-G brightness differences:")
		print("Red brighter: " + str(rg))
		print("Green brighter: " + str(gr))
		print("Equal: " + str(equal))
		
		# red-blue
		rb = 0 # red brighter
		br = 0 # blue brighter
		equal = 0
		for i in range(0, 4):
			if (i == 0):
				rlevel = red25e
			elif (i == 1):
				rlevel = red25e_03
			elif (i == 2):
				rlevel = red25e_07
			elif (i == 3):
				rlevel = red25e_10
			for j in range(0, 4):
				if (j == 0):
					blevel = blue47
				elif (j == 1):
					blevel = blue47_03
				elif (j == 2):
					blevel = blue47_07
				elif (j == 3):
					blevel = blue47_10
				diff = red25l[i] - blue47l[j]
				contrast = brightness_contrast(rlevel, blevel)
				if (contrast >= 1):
					if (diff > 0):
						rb += 1
					elif (diff < 0):
						br += 1
				else:
					equal += 1
		print("R-B brightness differences:")
		print("Red brighter: " + str(rb))
		print("Blue brighter: " + str(br))
		print("Equal: " + str(equal))
		
		# yellow-green
		yg = 0 # yellow brighter
		gy = 0 # green brighter
		equal = 0
		for i in range(0, 4):
			if (i == 0):
				ylevel = yellow15e
			elif (i == 1):
				ylevel = yellow15e_03
			elif (i == 2):
				ylevel = yellow15e_07
			elif (i == 3):
				ylevel = yellow15e_10
			for j in range(0, 4):
				if (j == 0):
					glevel = green58e
				elif (j == 1):
					glevel = green58e_03
				elif (j == 2):
					glevel = green58e_07
				elif (j == 3):
					glevel = green58e_10
				diff = yellow15l[i] - green58l[j]
				contrast = brightness_contrast(ylevel, glevel)
				if (contrast >= 1):
					if (diff > 0):
						yg += 1
					elif (diff < 0):
						gy += 1
				else:
					equal += 1
		print("Y-G brightness differences:")
		print("Yellow brighter: " + str(yg))
		print("Green brighter: " + str(gy))
		print("Equal: " + str(equal))
		
		# yellow-blue
		yb = 0 # yellow brighter
		by = 0 # blue brighter
		equal = 0
		for i in range(0, 4):
			if (i == 0):
				ylevel = yellow15e
			elif (i == 1):
				ylevel = yellow15e_03
			elif (i == 2):
				ylevel = yellow15e_07
			elif (i == 3):
				ylevel = yellow15e_10
			for j in range(0, 4):
				if (j == 0):
					blevel = blue47
				elif (j == 1):
					blevel = blue47_03
				elif (j == 2):
					blevel = blue47_07
				elif (j == 3):
					blevel = blue47_10
				diff = yellow15l[i] - blue47l[j]
				contrast = brightness_contrast(ylevel, blevel)
				if (contrast >= 1):
					if (diff > 0):
						yb += 1
					elif (diff < 0):
						by += 1
				else:
					equal += 1
		print("Y-B brightness differences:")
		print("Yellow brighter: " + str(yb))
		print("Blue brighter: " + str(by))
		print("Equal: " + str(equal))
		
		# green-blue
		gb = 0 # green brighter
		bg = 0 # blue brighter
		equal = 0
		for i in range(0, 4):
			if (i == 0):
				glevel = green58e
			elif (i == 1):
				glevel = green58e_03
			elif (i == 2):
				glevel = green58e_07
			elif (i == 3):
				glevel = green58e_10
			for j in range(0, 4):
				if (j == 0):
					blevel = blue47
				elif (j == 1):
					blevel = blue47_03
				elif (j == 2):
					blevel = blue47_07
				elif (j == 3):
					blevel = blue47_10
				diff = green58l[i] - blue47l[j]
				contrast = brightness_contrast(glevel, blevel)
				if (contrast >= 1):
					if (diff > 0):
						gb += 1
					elif (diff < 0):
						bg += 1
				else:
					equal += 1
		print("G-B brightness differences:")
		print("Green brighter: " + str(gb))
		print("Blue brighter: " + str(bg))
		print("Equal: " + str(equal))
		
		# significance level of brightness choices assuming a normal distribution
		integral = 0
		for i in range(0, 16):
			integral += normal(8, 8/3, i)
		tail = 0
		for i in range(max(ry, yr), 16):
			tail += normal(8, 8/3, i)
		
		print("R-Y significance: " + str(tail / integral))
		
		integral = 0
		for i in range(0, 16):
			integral += normal(8, 8/3, i)
		tail = 0
		for i in range(max(rg, gr), 16):
			tail += normal(8, 8/3, i)
		
		print("R-G significance: " + str(tail / integral))
		
		integral = 0
		for i in range(0, 16):
			integral += normal(8, 8/3, i)
		tail = 0
		for i in range(max(rb, br), 16):
			tail += normal(8, 8/3, i)
		
		print("R-B significance: " + str(tail / integral))
		
		integral = 0
		for i in range(0, 16):
			integral += normal(8, 8/3, i)
		tail = 0
		for i in range(max(yg, gy), 16):
			tail += normal(8, 8/3, i)
		
		print("Y-G significance: " + str(tail / integral))
		
		integral = 0
		for i in range(0, 16):
			integral += normal(8, 8/3, i)
		tail = 0
		for i in range(max(yb, by), 16):
			tail += normal(8, 8/3, i)
		
		print("Y-B significance: " + str(tail / integral))
		
		integral = 0
		for i in range(0, 16):
			integral += normal(8, 8/3, i)
		tail = 0
		for i in range(max(gb, bg), 16):
			tail += normal(8, 8/3, i)
		
		print("G-B significance: " + str(tail / integral))
		
		# plot
		x = np.array([
			"B+1.0",
			"G+1.0",
			"Y+1.0",
			"R+1.0",
			"B+0.7",
			"G+0.7",
			"Y+0.7",
			"R+0.7",
			"B+0.3",
			"G+0.3",
			"Y+0.3",
			"R+0.3",
			"B",
			"G",
			"Y",
			"R"
		])
		y = np.array([
			blue47l[3],
			green58l[3],
			yellow15l[3],
			red25l[3],
			blue47l[2],
			green58l[2],
			yellow15l[2],
			red25l[2],
			blue47l[1],
			green58l[1],
			yellow15l[1],
			red25l[1],
			blue47l[0],
			green58l[0],
			yellow15l[0],
			red25l[0]
		])
		x = np.array([0.0, 0.3, 0.7, 1.0])
		#plt.barh(x, y, color="black")
		#plt.xlabel("Perceived brightness")
		plt.plot(x, red25l, 'o-r', mec='k', label="red 25")
		plt.plot(x, yellow15l, 's-y', mec='k', label="yellow 15")
		plt.plot(x, green58l, '^-g', mec='k', label="green 58")
		plt.plot(x, blue47l, 'D-b', mec='k', label="blue 47")
		plt.xlabel("Filter optical density")
		plt.ylabel("Perceived brightness")
		plt.yscale("log")
		plt.legend()
		plt.show()
		
		x = np.array(["R-Y", "R-G", "R-B", "Y-G", "Y-B", "G-B"])
		y = np.empty(6)
		y[0] = (max(ry, yr) / 16) * 100
		y[1] = (max(rg, gr) / 16) * 100
		y[2] = (max(rb, br) / 16) * 100
		y[3] = (max(yg, gy) / 16) * 100
		y[4] = (max(yb, by) / 16) * 100
		y[5] = (max(gb, bg) / 16) * 100
		plt.bar(x, y, color="black")
		plt.ylabel("% correct")
		plt.show()
		
		# color differences -- this is the hard part...
		
		# first we plot them in a two-dimensional color space with (L-S)/(L+S) on
		# the x-axis and brightness on the y-axis. Not sure how much this tells us
		# though.
		if (l1 == m1):
			plt.plot(red25cs, red25l, 'or', mec='k', label="red 25")
			plt.plot(yellow15cs, yellow15l, 'sy', mec='k', label="yellow 15")
			plt.plot(green58cs, green58l, '^g', mec='k', label="green 58")
			plt.plot(blue47cs, blue47l, 'Db', mec='k', label="blue 47")
			plt.xlabel("Chromaticity ((L-S)/(L+S))")
			plt.ylabel("Brightness (L)")
			plt.legend()
			plt.show()
		else:
			plt.plot(red25cs, red25cs1, 'or', mec='k', label="red 25")
			plt.plot(yellow15cs, yellow15cs1, 'sy', mec='k', label="yellow 15")
			plt.plot(green58cs, green58cs1, '^g', mec='k', label="green 58")
			plt.plot(blue47cs, blue47cs1, 'Db', mec='k', label="blue 47")
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
		
		# Next we try to assess whether they're distinguishable. There appears
		# to be no chromatic overlap, so we won't bother with the direction.
		
		y = np.empty(6)
		
		# R-Y
		d = 0 # different
		n = 0 # near threshold
		s = 0 # same
		
		print("R-Y")
		for i in range(0, 4):
			if (i == 0):
				rlevel = red25e
			elif (i == 1):
				rlevel = red25e_03
			elif (i == 2):
				rlevel = red25e_07
			elif (i == 3):
				rlevel = red25e_10
			for j in range(0, 4):
				if (j == 0):
					ylevel = yellow15e
				elif (j == 1):
					ylevel = yellow15e_03
				elif (j == 2):
					ylevel = yellow15e_07
				elif (j == 3):
					ylevel = yellow15e_10
				contrast = color_contrast(rlevel, ylevel)
				print("R" + str(i) + " vs. Y" + str(j) + ": " + str(contrast))
				if (contrast >= 4):
					d += 1
				elif (contrast >= 1):
					n += 1
				else:
					s += 1
		print("Different: " + str(d))
		print("Close: " + str(n))
		print("Same: " + str(s))
		print("")
		y[0] = d+n
		
		# R-G
		d = 0
		n = 0
		s = 0
		
		print("R-G")
		for i in range(0, 4):
			if (i == 0):
				rlevel = red25e
			elif (i == 1):
				rlevel = red25e_03
			elif (i == 2):
				rlevel = red25e_07
			elif (i == 3):
				rlevel = red25e_10
			for j in range(0, 4):
				if (j == 0):
					glevel = green58e
				elif (j == 1):
					glevel = green58e_03
				elif (j == 2):
					glevel = green58e_07
				elif (j == 3):
					glevel = green58e_10
				contrast = color_contrast(rlevel, glevel)
				print("R" + str(i) + " vs. G" + str(j) + ": " + str(contrast))
				if (contrast >= 4):
					d += 1
				elif (contrast >= 1):
					n += 1
				else:
					s += 1
		print("Different: " + str(d))
		print("Close: " + str(n))
		print("Same: " + str(s))
		print("")
		y[1] = d+n
		
		# R-B
		d = 0
		n = 0
		s = 0
		
		print("R-B")
		for i in range(0, 4):
			if (i == 0):
				rlevel = red25e
			elif (i == 1):
				rlevel = red25e_03
			elif (i == 2):
				rlevel = red25e_07
			elif (i == 3):
				rlevel = red25e_10
			for j in range(0, 4):
				if (j == 0):
					blevel = blue47
				elif (j == 1):
					blevel = blue47_03
				elif (j == 2):
					blevel = blue47_07
				elif (j == 3):
					blevel = blue47_10
				contrast = color_contrast(rlevel, blevel)
				print("R" + str(i) + " vs. B" + str(j) + ": " + str(contrast))
				if (contrast >= 4):
					d += 1
				elif (contrast >= 1):
					n += 1
				else:
					s += 1
		print("Different: " + str(d))
		print("Close: " + str(n))
		print("Same: " + str(s))
		print("")
		y[2] = d+n
		
		# Y-G
		d = 0
		n = 0
		s = 0
		
		print("Y-G")
		for i in range(0, 4):
			if (i == 0):
				ylevel = yellow15e
			elif (i == 1):
				ylevel = yellow15e_03
			elif (i == 2):
				ylevel = yellow15e_07
			elif (i == 3):
				ylevel = yellow15e_10
			for j in range(0, 4):
				if (j == 0):
					glevel = green58e
				elif (j == 1):
					glevel = green58e_03
				elif (j == 2):
					glevel = green58e_07
				elif (j == 3):
					glevel = green58e_10
				contrast = color_contrast(ylevel, glevel)
				print("Y" + str(i) + " vs. G" + str(j) + ": " + str(contrast))
				if (contrast >= 4):
					d += 1
				elif (contrast >= 1):
					n += 1
				else:
					s += 1
		print("Different: " + str(d))
		print("Close: " + str(n))
		print("Same: " + str(s))
		print("")
		y[3] = d+n
		
		# Y-B
		d = 0
		n = 0
		s = 0
		
		print("Y-B")
		for i in range(0, 4):
			if (i == 0):
				ylevel = yellow15e
			elif (i == 1):
				ylevel = yellow15e_03
			elif (i == 2):
				ylevel = yellow15e_07
			elif (i == 3):
				ylevel = yellow15e_10
			for j in range(0, 4):
				if (j == 0):
					blevel = blue47
				elif (j == 1):
					blevel = blue47_03
				elif (j == 2):
					blevel = blue47_07
				elif (j == 3):
					blevel = blue47_10
				contrast = color_contrast(rlevel, blevel)
				print("Y" + str(i) + " vs. B" + str(j) + ": " + str(contrast))
				if (contrast >= 4):
					d += 1
				elif (contrast >= 1):
					n += 1
				else:
					s += 1
		print("Different: " + str(d))
		print("Close: " + str(n))
		print("Same: " + str(s))
		print("")
		y[4] = d+n
		
		# G-B
		d = 0
		n = 0
		s = 0
		
		print("G-B")
		for i in range(0, 4):
			if (i == 0):
				glevel = green58e
			elif (i == 1):
				glevel = green58e_03
			elif (i == 2):
				glevel = green58e_07
			elif (i == 3):
				glevel = green58e_10
			for j in range(0, 4):
				if (j == 0):
					blevel = blue47
				elif (j == 1):
					blevel = blue47_03
				elif (j == 2):
					blevel = blue47_07
				elif (j == 3):
					blevel = blue47_10
				contrast = color_contrast(glevel, blevel)
				print("G" + str(i) + " vs. B" + str(j) + ": " + str(contrast))
				if (contrast >= 4):
					d += 1
				elif (contrast >= 1):
					n += 1
				else:
					s += 1
		print("Different: " + str(d))
		print("Close: " + str(n))
		print("Same: " + str(s))
		print("")
		y[5] = d+n
		
		# percentages
		for i in range(6):
			y[i] = y[i]*100 / 16
		
		# plot
		plt.bar(x, y, color="black")
		plt.ylabel("% correct")
		plt.show()

# print execution time
print("%s seconds" % (time.time() - start_time))
