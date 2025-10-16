"""
This file contains the following:
* list of arguments
* visual pigments: standard visual pigment template (vpt()) and human cone fundamentals
* ocular media transmittance for various species
* blackbody: blackbody()
* color/brightness analysis: color_contrast(), brightness_contrast() and spec2rgb()
* illuminants: D65, E, incandescent (identical to CIE A by default but can be customized)
* color matching functions
* Maxwell triangle color space plots: triangle()
* color and brightness discrimination: color_disc() and brightness_disc()

Other stuff has been spun off to the following files, which import this one. I'd like to add additional arguments to those but
don't see how to do that.
* erg.py: analysis of ERG data from Jacobs & Williams (2010)
* illuminants.py: test illuminants
* optimal.py: optimal colors
* peak.py: Gaussian peak spectra
* ramp.py: ramp function spectra
* sky.py: sky/daylight spectra
* wratten.py: Kodak Wratten filters
* munsell.py: Munsell color cards
* image_processing.py: tool for creating false-color images
** This one can't be part of the main file because OpenCV and matplotlib don't like each other. That's
why the main file has the part for generating and plotting the source->target curves.

wratten.py and munsell.py just produce plots by default. To see the color rendering and trial modeling, you have to use
--wratten/--munsell. This is because they also contain other features that are accessed with other arguments (see the list).

The image processing functions that used to be here have also been moved into another file because they
didn't get along with matplotlib.
"""

# this list doesn't recursively import
import math
import numpy as np
import time
import argparse
import colormath
from colormath.color_objects import sRGBColor, SpectralColor
from colormath.color_conversions import convert_color
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binomtest
import statistics
from scipy.integrate import quad
import time
import csv

# execution time
start_time = time.time()

# arguments
# The absolute quantum catch arguments mostly have default values specific to my analysis of
# opossums. I may change this.
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--receptors", "--target", nargs="+", type=int, default=[562], help="receptor peak sensitivities")
parser.add_argument("--rcsv", nargs="+", default="none", help="custom receptors from CSV file(s)")
parser.add_argument("--rod", type=int, default=0, help="rod peak sensitivity")
parser.add_argument("--weber", type=float, default=0.05, help="base Weber fraction")
parser.add_argument("--weberr2", type=float, default=0, help="override Weber fraction for receptor 2")
parser.add_argument("--weberr3", type=float, default=0, help="override Weber fraction for receptor 3")
parser.add_argument("--weberr4", type=float, default=0, help="override Weber fraction for receptor 4")
parser.add_argument("--weberb", type=float, default=0.11, help="achromatic Weber fraction")
parser.add_argument("--r1p", type=float, default=1, help="number of receptor 1 per receptive field")
parser.add_argument("--r2p", type=float, default=1, help="number of receptor 2 per receptive field")
parser.add_argument("--r3p", type=float, default=1, help="number of receptor 3 per receptive field")
parser.add_argument("--r4p", type=float, default=1, help="number of receptor 4 per receptive field")
parser.add_argument("--media", "--filter", type=str, default="none", help="ocular media transmittance")
parser.add_argument("--qn", help="use quantum noise in color differences", action="store_true")
parser.add_argument("--r1s", type=float, default=1, help="contribution of receptor 1 to spectral sensitivity")
parser.add_argument("--rs", type=float, default=0, help="contribution of rods to spectral sensitivity")
parser.add_argument("--vonkries", help="enable von Kries transform", action="store_true")
parser.add_argument("--luminosity", help="show spectral sensitivity function", action="store_true")
parser.add_argument("--white", default="d65", help="set illuminant")
parser.add_argument("--wratten", help="Kodak Wratten camera filters", action="store_true")
parser.add_argument("--wv2", help="use tabulated Kodak Wratten data from Kodak Photographic Filters Handbook (1990)", action="store_true")
parser.add_argument("--wcheck", help="check camera filter brightness matches", action="store_true")
parser.add_argument("--wopt1", help="test varying spectral sensitivity (probability)", action="store_true")
parser.add_argument("--wopt2", help="test varying spectral sensitivity (order)", action="store_true")
parser.add_argument("--wopt3", help="test varying spectral sensitivity (contrast)", action="store_true")
parser.add_argument("--qcheck", help="show absolute quantum catches", action="store_true")
parser.add_argument("--vpt", help="plot visual pigment templates", action="store_true")
parser.add_argument("--csplot", "--triangle", nargs="*", help="create a Maxwell-style color space plot")
parser.add_argument("--r1name", default='', help="name of receptor 1 (used in triangle/line plots)")
parser.add_argument("--r2name", default='', help="name of receptor 2 (used in triangle/line plots)")
parser.add_argument("--r3name", default='', help="name of receptor 3 (used in triangle/line plots)")
parser.add_argument("--fcmode", help="switch default false color mode", default="none")
parser.add_argument("--s2tmode", help="switch default color transform mode", default="c")
parser.add_argument("--blackbody", type=int, default=0, help="test Wratten filters with a specified blackbody temperature")
parser.add_argument("--munsell", help="Munsell color cards", action="store_true")
parser.add_argument("--mv2", help="use alternative Munsell spectra", action="store_true")
parser.add_argument("--mv3", help="use alternative Munsell spectra for <420nm only", action="store_true")
parser.add_argument("--mv4", help="remove wavelengths <400nm from Munsell spectra", action="store_true")
parser.add_argument("--mopt1", help="find optimal visual pigments for Munsell cards (dichromacy)", action="store_true")
parser.add_argument("--mopt2", help="find optimal visual pigments for Munsell cards (trichromacy)", action="store_true")
parser.add_argument("--mopt3", help="find L-M chromatic spread/overlap of Munsell cards", action="store_true")
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
parser.add_argument("--primaries", nargs="+", type=int, default=[700, 546.1, 435.8, 0, 0], help="primary wavelengths for color-matching functions")
parser.add_argument("-p", "--preset", help="use specified preset vision system", default="none")
parser.add_argument("-v", "--verbose", help="", action="store_true")
parser.add_argument("-i", "--image", default="none", help="path to input image")
parser.add_argument("-u", "--uvimage", default="none", help="path to input UV image")
parser.add_argument("--source", nargs="+", help="source phenotype")
parser.add_argument("--s2tplot", help="source->target transform plot", action="store_true")
parser.add_argument("--bound", default=-np.inf, type=float, help="lower bound for source->target coefficients")
parser.add_argument("--ingamma", default=-1, type=float, help="input gamma")
parser.add_argument("--outgamma", default=-1, type=float, help="output gamma")
parser.add_argument("--interactions", default=1, type=int, help="interactions between terms for curve fitting")
parser.add_argument("--luminance", help="type of luminance adjustment", default="default")
args = parser.parse_args()

# X values representing wavelengths 300-700 nm, used by multiple functions
x_10nm = np.empty(41)
for i in range(41): x_10nm[i] = i*10 + 300
x_1nm = np.empty(401)
for i in range(401): x_1nm[i] = i + 300

# receptor peak sensitivities
receptors = args.receptors
rod = args.rod
dimension = len(args.receptors)

# Weber fractions
wr1 = args.weber
wr2 = wr1 * math.sqrt(args.r1p) / math.sqrt(args.r2p)
wr3 = wr1 * math.sqrt(args.r1p) / math.sqrt(args.r3p)
wr4 = wr1 * math.sqrt(args.r1p) / math.sqrt(args.r4p)
if (args.weberr2 > 0): wr2 = args.weberr2
if (args.weberr3 > 0): ws = args.weberr3
if (args.weberr4 > 0): ws = args.weberr4

# read CSV files
# This function allows for using arbitrary CSV files as visual pigments, ocular media, illuminants,
# or reflectance spectra. At the moment it can only handle files that contain a single spectrum.
# More would be nice but not sure how yet.
def csv2spec(filename):
	if (args.verbose): print("Reading CSV file: " + filename)
	with open(filename, newline='') as csvfile:
		reader = csv.DictReader(csvfile, fieldnames=['wl', 'value'])
		wl = []
		values = []
		# This bit uses an exception (can be either ValueError or TypeError) to
		# ignore headers and other extra text. This kind of thing is considered bad
		# coding practice but seems to be the only way to get the behavior I want
		# when reading PlotDigitizer files.
		for row in reader:
			try:
				wl_float = float(row['wl'])
				value_float = float(row['value'])
				wl.append(wl_float)
				values.append(value_float)
			except:
				if (args.warnings):
					print("The following row of " + filename + " was ignored: "
					+ str(row))
	
	# interpolate to 300-700
	values_1nm = np.interp(x_1nm, wl, values)

	# remove extrapolation
	wl_min = min(*wl)
	wl_max = max(*wl)
	
	if (wl_min > 300):
		for i in range(math.floor(wl_min) - 300): values_1nm[i] = 0
	if (wl_max < 700):
		for i in range(math.ceil(wl_max) - 300, 401): values_1nm[i] = 0
	
	return values_1nm

# types of ocular media

# mouse (Jacobs & Williams 2007)
# Transmission below 310 nm is not provided and has been set to 0. This means we're effectively
# only considering wavelengths >=310, but it's easier to keep track of the numbers if they begin
# with a multiple of 100.
mouse_lens = np.array([
	np.nan, # 300 -- prevent interpolation
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

# divide by 100
mouse_lens = mouse_lens / 100

# 1-nm intervals
mouse_1nm = np.interp(x_1nm, x_10nm, mouse_lens)
for i in range(mouse_1nm.shape[0]):
	if (np.isnan(mouse_1nm[i])): mouse_1nm[i] = 0 # remove nans

# 10-nm intervals
mouse_10nm = mouse_lens
mouse_10nm[0] = 0 # remove nan

# choose media type
media_1nm = np.ones(401)
media_10nm = np.ones(41)
if (args.media == "mouse"):
	media_1nm = mouse_1nm
	media_10nm = mouse_10nm
elif (args.media == "human"):
	for i in range(401):
		w = i + 300
		try:
			density = 1.1*math.exp((400 - w) / 15) + 0.11*math.exp((500 - w) / 80)
			media_1nm[i] = 10**(-density)
		except OverflowError:
			media_1nm[i] = 0
	for i in range(41):
		w = i*10 + 300
		try:
			density = 1.1*math.exp((400 - w) / 15) + 0.11*math.exp((500 - w) / 80)
			media_10nm[i] = 10**(-density)
		except OverflowError:
			media_10nm[i] = 0
elif (args.media != "none"):
	media_1nm = csv2spec(args.media)
	media_10nm = np.interp(x_10nm, x_1nm, values)
	
	# set to 1 at 700 nm
	media_10nm = media_10nm / media_10nm[40]
	media_1nm = media_1nm / media_1nm[400]

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

d65_10nm = np.empty(41)
for i in range(41): d65_10nm[i] = d65_1nm[i*10]

e_10nm = np.full(41, 100)
e_1nm = np.full(401, 100)

# incandescent lighting with specified color temperature
incandescent_10nm = np.empty(41)
for i in range(41):
	w = i*10 + 300
	incandescent_10nm[i] = blackbody(w, args.ct)
incandescent_10nm = 100*incandescent_10nm / max(*incandescent_10nm) # 0-100
# 1-nm resolution
incandescent_1nm = np.empty(401)
for i in range(401):
	w = i + 300
	incandescent_1nm[i] = blackbody(w, args.ct)
incandescent_1nm = 100*incandescent_1nm / max(*incandescent_1nm) # 0-100

# absolute number of photons
# Since the output of blackbody() is joules/m^3/sec/sr, this should be photons/m^3/sec/sr.
# I think what we really want is square meters. (No, the meters cancel, see the quantum
# noise calculations)
# Changing this to 1-nm intervals because 10-nm isn't very useful.
incandescent_photons = np.empty(501)
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
if (args.height > 0): area = args.width * args.height
else: area = math.pi*(args.width/2)**2
# steradians -- not divided by the whole sphere surface area, just the square of the radius
sr = area / args.dist**2
watts_cm2 = args.watts / (sphere_area)
watts_m2 = watts_cm2 * 100**2
scale = watts_m2 / energy # scaling factor -- units should cancel (W/m^2)
#scale = args.watts / power
#print("W/m^2 scaling factor: " + str(scale))

for i in range(501):
	w = i + 300 # nanometers
	incandescent_photons[i] = scale*1e-9*w*blackbody(w, args.ct) / (h*c) # meters * joules/m^2/nm/sec/sr / joule*sec * m/s

# absolute number of photons, part 2: scaling D65 to a specified value in lux
# For values not included in Python's D65, we approximate this with a black body with
# temperature 6504 K. (Never mind, I found the real thing)

# human photopic luminosity function (2019 CIE standard)
# No values are provided below 360 nm, so we assume 0.
cie_luminosity = np.array([
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

# integral of D65 x human photopic luminosity function (we hardcode this instead of using
# sensitivity() because that varies with settings for the specified species)
# The units we want at the end are W/m^2, so we multiply this by the lux-to-watts scaling
# factor for 555 nm (683.002 lm/W). This will end up on the bottom.
integral = sum(d65_1nm * cie_luminosity) * 683.002

# lux scaling factor (units = lx / lx/nm = nm)
lxscale = args.lux / integral
#print(lxscale)

# scaled version (W/m^2 * m / joule*sec * m/s * nm = joule/sec/m^2 * m / joule/sec * m/s * nm * sr)
# Per nm is implied. The values we just found do not include steradians (probably), so we divide
# by the sr value we found earlier. (No we don't)
# Also changing this to 1 nm.
d65_photons = np.empty(531)
for i in range(d65_photons.shape[0]):
	w = i + 300
	d65_photons[i] = lxscale*1e-9*w*d65_1nm[i] / (h*c*math.pi)

e_photons = np.empty(531)
for i in range(e_photons.shape[0]):
	w = i + 300
	e_photons[i] = lxscale*1e-9*w*100 / (h*c*math.pi)

# white point
# Selecting the absolute number of photons is now done here to avoid going through an if statement
# every time color_contrast() is run.
if (args.white == "d65"):
	wp_10nm = d65_10nm
	wp_1nm = np.empty(401)
	for i in range(401): wp_1nm[i] = d65_1nm[i]
	wp_photons = d65_photons
elif (args.white == "e"):
	wp_10nm = e_10nm
	wp_1nm = e_1nm
	wp_photons = e_photons
elif (args.white == "a" or args.white == "i"):
	wp_10nm = incandescent_10nm
	wp_1nm = incandescent_1nm
	wp_photons = incandescent_photons
# custom illuminant
else:
	wp_1nm = csv2spec(args.white)
	wp_10nm = np.interp(x_10nm, x_1nm, wp_1nm)

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

# find sensitivity of a given photoreceptor type to a wavelength using visual pigment templates
# from Govardovskii et al. (2000). l is short for lambda. See also https://pmc.ncbi.nlm.nih.gov/articles/PMC2962788/
# Palacios et al. 2010 use a different equation for the beta-band peak (123 + 0.429 * alpha
# peak).
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
		if (args.warnings): print("Warning (vpt): OverflowError, clipping to 2.2250738585072014e-308")
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

# presets

# CIE 2-deg fundamentals
cie2_l = np.zeros(401)
cie2_m = np.zeros(401)
cie2_s = np.zeros(401)
for i in range(390, 701):
	cie2_l[i-300] = 10**(hcf2deg[int(round(i))-390][1])
	cie2_m[i-300] = 10**(hcf2deg[int(round(i))-390][2])
	cie2_s[i-300] = 10**(hcf2deg[int(round(i))-390][3])
cie2 = [cie2_l, cie2_m, cie2_s]

# CIE 10-deg fundamentals
cie10_l = np.zeros(401)
cie10_m = np.zeros(401)
cie10_s = np.zeros(401)
cie10_luminosity = np.zeros(401)
for i in range(390, 701):
	cie10_l[i-300] = 10**(hcf10deg[int(round(i))-390][1])
	cie10_m[i-300] = 10**(hcf10deg[int(round(i))-390][2])
	cie10_s[i-300] = 10**(hcf10deg[int(round(i))-390][3])
	cie10_luminosity[i-300] = cie_luminosity[i-300] # just so it's the right length
cie10 = [cie10_l, cie10_m, cie10_s]

# typical UVS and VS birds, digitized with WebPlotDigitizer
# source: https://www.nature.com/articles/s41467-018-08142-5
# licensed under Creative Commons 4.0 http://creativecommons.org/licenses/by/4.0/
bird_l = csv2spec('bird-l.csv')
bird_m = csv2spec('bird-m.csv')
bird_su = csv2spec('bird-su.csv')
bird_sv = csv2spec('bird-sv.csv')
bird_uv = csv2spec('bird-uv.csv')
bird_v = csv2spec('bird-v.csv')

uvbird = [bird_l,bird_m,bird_su,bird_uv]
vbird = [bird_l,bird_m,bird_sv,bird_v]

# set receptors and luminosity
rlist = []
r1_1nm = np.zeros(401)
r2_1nm = np.zeros(401)
r3_1nm = np.zeros(401)
r4_1nm = np.zeros(401)
r5_1nm = np.zeros(401)
luminosity = np.empty(401)

if (args.preset == "cie2"):
	rlist = cie2
	luminosity = cie10_luminosity
# CIE 10-deg fundamentals
elif (args.preset == "cie10"):
	rlist = cie10
	luminosity = cie10_luminosity
# UVS bird
elif (args.preset == "uvbird"): rlist = uvbird
# VS bird
elif (args.preset == "vbird"): rlist = vbird
elif (args.rcsv != "none"):
	dimension = len(args.rcsv)
	for i in range(len(args.rcsv)):
		rlist.append(csv2spec(args.rcsv[i]))
else:
	rod = args.rod
	for i in range(401):
		w = i+300
		if (args.rod > 0):
			luminosity[i] = (args.r1s*vpt(w, receptors[0]) + args.rs*vpt(w, rod)) * media_1nm[w-300]
		else:
			luminosity[i] = (args.r1s*vpt(w, receptors[0])) * media_1nm[w-300]

	for i in range(dimension):
		# whole list
		yvalues = np.empty(401)
		peak = args.receptors[i] # peak sensitivity
		for j in range(401):
			yvalues[j] = vpt(j+300, peak)
		rlist.append(yvalues)

dimension = len(rlist)
if (dimension > 4): r5_1nm = rlist[4]
if (dimension > 3): r4_1nm = rlist[3]
if (dimension > 2): r3_1nm = rlist[2]
if (dimension > 1): r2_1nm = rlist[1]
r1_1nm = rlist[0]

"""
CMFs (color-matching functions)

The names "red", "green" and "blue" for the first 3 primaries have been removed, but the default
values are still the CIE RGB primaries (700, 546.1 and 435.8). Since the second two are not
integers, they have to be rounded to 546 and 436 so we can look up array indexes. (It's faster
that way and probably doesn't matter for precision.)

Other recommended values for trichromatic primaries:
* marsupials: 620, 450, 360 (Arrese et al. 2006 doi.org/10.1016/j.cub.2006.02.036)
* humans/primates:
** 630, 532, 467 (Rec. 2020 https://en.wikipedia.org/wiki/Rec._2020)
** 610, 550, 465 (approximate dominant wavelengths of sRGB primaries: see https://en.wikipedia.org/wiki/File:SRGB_chromaticity_CIE1931.svg)
** other color spaces: see https://en.wikipedia.org/wiki/RGB_color_spaces

Dichromatic, trichromatic, tetrachromatic and pentachromatic CMFs are supported.
"""
# primaries
p1 = round(args.primaries[0])
p2 = round(args.primaries[1])
if (len(args.primaries) > 2): p3 = round(args.primaries[2])
if (len(args.primaries) > 3): p4 = round(args.primaries[3])
if (len(args.primaries) > 4): p5 = round(args.primaries[4])

# color coding

def cmf():
	if (dimension > 4):
		# sensitivity to primaries
		matrix_a = np.array([
			[r1_1nm[p1-300]*media_1nm[p1-300], r1_1nm[p2-300]*media_1nm[p2-300], r1_1nm[p3-300]*media_1nm[p3-300], r1_1nm[p4-300]*media_1nm[p4-300], r1_1nm[p5-300]*media_1nm[p5-300]],
			[r2_1nm[p1-300]*media_1nm[p1-300], r2_1nm[p2-300]*media_1nm[p2-300], r2_1nm[p3-300]*media_1nm[p3-300], r2_1nm[p4-300]*media_1nm[p4-300], r2_1nm[p5-300]*media_1nm[p5-300]],
			[r3_1nm[p1-300]*media_1nm[p1-300], r3_1nm[p2-300]*media_1nm[p2-300], r3_1nm[p3-300]*media_1nm[p3-300], r3_1nm[p4-300]*media_1nm[p4-300], r3_1nm[p5-300]*media_1nm[p5-300]],
			[r4_1nm[p1-300]*media_1nm[p1-300], r4_1nm[p2-300]*media_1nm[p2-300], r4_1nm[p3-300]*media_1nm[p3-300], r4_1nm[p4-300]*media_1nm[p4-300], r4_1nm[p5-300]*media_1nm[p5-300]],
			[r5_1nm[p1-300]*media_1nm[p1-300], r5_1nm[p2-300]*media_1nm[p2-300], r5_1nm[p3-300]*media_1nm[p3-300], r5_1nm[p4-300]*media_1nm[p4-300], r5_1nm[p5-300]*media_1nm[p5-300]],
		])

		# sensitivity to each wavelength from 300 to 700 nm in 1 nm increments
		matrix_c = np.empty((5, 401)) # 400x5 matrix -- this is height x width, not width x height

		for i in range(401):
		# p1 row
			matrix_c[0][i] = r1_1nm[i] * media_1nm[i]
		# p2 row
			matrix_c[1][i] = r2_1nm[i] * media_1nm[i]
		# p3 row
			matrix_c[2][i] = r3_1nm[i] * media_1nm[i]
		# p4 row
			matrix_c[3][i] = r4_1nm[i] * media_1nm[i]
		# p4 row
			matrix_c[4][i] = r5_1nm[i] * media_1nm[i]

		# initial color match matrix
		matrix_m = np.matmul(np.linalg.inv(matrix_a), matrix_c) # A x M = C, so M = A^-1 x C

		# adjust to reference white
		#watts = 0
		#for i in range(401): watts += wp_1nm[i] # scale to 1 watt
		p1w = sum(matrix_m[0] * wp_1nm)
		p2w = sum(matrix_m[1] * wp_1nm)
		p3w = sum(matrix_m[2] * wp_1nm)
		p4w = sum(matrix_m[3] * wp_1nm)
		p5w = sum(matrix_m[4] * wp_1nm)
		matrix_w = [
			[1/p1w, 0, 0, 0, 0],
			[0, 1/p2w, 0, 0, 0],
			[0, 0, 1/p3w, 0, 0],
			[0, 0, 0, 1/p4w, 0],
			[0, 0, 0, 0, 1/p5w],
		]

		# final color match matrix
		matrix_cmf = np.matmul(matrix_w, matrix_m)
	elif (dimension == 4):
		# sensitivity to primaries
		matrix_a = np.array([
			[r1_1nm[p1-300]*media_1nm[p1-300], r1_1nm[p2-300]*media_1nm[p2-300], r1_1nm[p3-300]*media_1nm[p3-300], r1_1nm[p4-300]*media_1nm[p4-300]],
			[r2_1nm[p1-300]*media_1nm[p1-300], r2_1nm[p2-300]*media_1nm[p2-300], r2_1nm[p3-300]*media_1nm[p3-300], r2_1nm[p4-300]*media_1nm[p4-300]],
			[r3_1nm[p1-300]*media_1nm[p1-300], r3_1nm[p2-300]*media_1nm[p2-300], r3_1nm[p3-300]*media_1nm[p3-300], r3_1nm[p4-300]*media_1nm[p4-300]],
			[r4_1nm[p1-300]*media_1nm[p1-300], r4_1nm[p2-300]*media_1nm[p2-300], r4_1nm[p3-300]*media_1nm[p3-300], r4_1nm[p4-300]*media_1nm[p4-300]],
		])

		# sensitivity to each wavelength from 300 to 700 nm in 1 nm increments
		matrix_c = np.empty((4, 401)) # 400x4 matrix -- this is height x width, not width x height

		for i in range(401):
		# p1 row
			matrix_c[0][i] = r1_1nm[i] * media_1nm[i]
		# p2 row
			matrix_c[1][i] = r2_1nm[i] * media_1nm[i]
		# p3 row
			matrix_c[2][i] = r3_1nm[i] * media_1nm[i]
		# p4 row
			matrix_c[3][i] = r4_1nm[i] * media_1nm[i]

		# initial color match matrix
		matrix_m = np.matmul(np.linalg.inv(matrix_a), matrix_c) # A x M = C, so M = A^-1 x C

		# adjust to reference white
		#watts = 0
		#for i in range(401): watts += wp_1nm[i] # scale to 1 watt
		p1w = sum(matrix_m[0] * wp_1nm)
		p2w = sum(matrix_m[1] * wp_1nm)
		p3w = sum(matrix_m[2] * wp_1nm)
		p4w = sum(matrix_m[3] * wp_1nm)
		matrix_w = [
			[1/p1w, 0, 0, 0],
			[0, 1/p2w, 0, 0],
			[0, 0, 1/p3w, 0],
			[0, 0, 0, 1/p4w],
		]

		# final color match matrix
		matrix_cmf = np.matmul(matrix_w, matrix_m)
	elif (dimension == 3):
		# sensitivity to primaries
		matrix_a = np.array([
			[r1_1nm[p1-300]*media_1nm[p1-300], r1_1nm[p2-300]*media_1nm[p2-300], r1_1nm[p3-300]*media_1nm[p3-300]],
			[r2_1nm[p1-300]*media_1nm[p1-300], r2_1nm[p2-300]*media_1nm[p2-300], r2_1nm[p3-300]*media_1nm[p3-300]],
			[r3_1nm[p1-300]*media_1nm[p1-300], r3_1nm[p2-300]*media_1nm[p2-300], r3_1nm[p3-300]*media_1nm[p3-300]],
		])

		# sensitivity to each wavelength from 300 to 700 nm in 1 nm increments
		matrix_c = np.empty((3, 401)) # 400x3 matrix -- this is height x width, not width x height

		for i in range(401):
		# red row
			matrix_c[0][i] = r1_1nm[i] * media_1nm[i]
		# green row
			matrix_c[1][i] = r2_1nm[i] * media_1nm[i]
		# blue row
			matrix_c[2][i] = r3_1nm[i] * media_1nm[i]

		# initial color match matrix
		matrix_m = np.matmul(np.linalg.inv(matrix_a), matrix_c) # A x M = C, so M = A^-1 x C

		# adjust to reference white
		#watts = 0
		#for i in range(401): watts += wp_1nm[i] # scale to 1 watt
		rw = sum(matrix_m[0] * wp_1nm)
		gw = sum(matrix_m[1] * wp_1nm)
		bw = sum(matrix_m[2] * wp_1nm)
		matrix_w = [
			[1/rw, 0, 0],
			[0, 1/gw, 0],
			[0, 0, 1/bw]
		]

		# final color match matrix
		matrix_cmf = np.matmul(matrix_w, matrix_m)
	elif (dimension == 2):

		# sensitivity to primaries
		matrix_a = np.array([
			[r1_1nm[p1-300] * media_1nm[p1-300], r1_1nm[p2-300] * media_1nm[p2-300]],
			[r2_1nm[p1-300] * media_1nm[p1-300], r2_1nm[p2-300] * media_1nm[p2-300]],
		])

		# sensitivity to each wavelength from 300 to 700 nm in 1 nm increments
		matrix_c = np.empty((2, 401)) # 400x2 matrix

		for i in range(401):
		# red row
			matrix_c[0][i] = r1_1nm[i] * media_1nm[i]
		# blue row
			matrix_c[1][i] = r2_1nm[i] * media_1nm[i]

		# initial color match matrix
		matrix_m = np.matmul(np.linalg.inv(matrix_a), matrix_c) # A x M = C, so M = A^-1 x C

		# adjust to reference white
		#watts = 0
		#for i in range(401): watts += wp_1nm[i] # scale to 1 watt
		rw = sum(matrix_m[0] * wp_1nm)
		gw = 0
		bw = sum(matrix_m[1] * wp_1nm)
		matrix_w = [
			[1/rw, 0],
			[0, 1/bw]
		]

		# final color match matrix
		matrix_cmf = np.matmul(matrix_w, matrix_m)
		
	return matrix_cmf

if (args.cmf):
	if (dimension < 2):
		print("CMFs are not meaningful for monochromacy.")
	else:
		matrix_cmf = cmf()[0]
		# plot curves
		plt.plot(x_1nm, matrix_cmf[0])
		if (dimension > 1): plt.plot(x_1nm, matrix_cmf[1], color='k')
		if (dimension > 2): plt.plot(x_1nm, matrix_cmf[2], color='k')
		if (dimension > 3): plt.plot(x_1nm, matrix_cmf[3], color='k')
		if (dimension > 4): plt.plot(x_1nm, matrix_cmf[4], color='k')
		plt.show()

# set source receptors (for image processing)
if (args.image != "none" or args.s2tplot):
	srlist = []
	if (args.source == "cie2"): srlist = cie2
	# CIE 10-deg fundamentals
	elif (args.source == "cie10"): srlist = cie10
	# UVS bird
	elif (args.source == "uvbird"): srlist = uvbird
	# VS bird
	elif (args.source == "vbird"): srlist = vbird
	else:
		for i in range(len(args.source)):
			srlist.append(csv2spec(args.source[i]))

	sr1 = srlist[0]
	sr2 = srlist[1]
	sr3 = srlist[2]
	sr4 = np.zeros(401)
	if ((len(srlist)) == 4): sr4 = srlist[3]

	# constants for luminance adjustment
	# sRGB -> XYZ Y: http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
	lum_r = 0.2126729
	lum_g = 0.7151522
	lum_b = 0.0721750
	
	# contribution of each channel to the target's luminance vision
	lum_r1 = sum(sr1 * luminosity)
	lum_g1 = sum(sr2 * luminosity)
	lum_b1 = sum(sr3 * luminosity)
	lum_uv = sum(sr4 * luminosity)
	lum_max = lum_r1 + lum_g1 + lum_b1 + lum_uv
	
	# adapt to illuminant
	sr1 *= wp_1nm
	sr2 *= wp_1nm
	sr3 *= wp_1nm
	sr4 *= wp_1nm
	wpsr1 = sum(sr1)
	wpsr2 = sum(sr2)
	wpsr3 = sum(sr3)
	wpsr4 = sum(sr4)
	sr1 /= wpsr1
	sr2 /= wpsr2
	sr3 /= wpsr3
	if (wpsr4 > 0): sr4 /= wpsr4
	
	# declare variables
	tr1 = np.zeros(401)
	tr2 = np.zeros(401)
	tr3 = np.zeros(401)
	wptr1 = 1
	wptr2 = 1
	wptr3 = 1

	# use CMFs instead of raw receptors
	# No von Kries transform because they're already adapted.
	if (args.fcmode == 'cmf'):
		matrix_cmf = cmf()
		tr1 = matrix_cmf[0]
		if (dimension > 1): tr2 = matrix_cmf[1]
		if (dimension > 2): tr3 = matrix_cmf[2]
		
		# cut off negative numbers -- they mess with chromaticity and represent something
		# that can't be displayed anyway
		for i in range(401):
			if (tr1[i] < 0): tr1[i] = 0
			if (tr2[i] < 0): tr2[i] = 0
			if (tr3[i] < 0): tr3[i] = 0
	else:
		tr1 = r1_1nm * media_1nm * wp_1nm
		tr2 = r2_1nm * media_1nm * wp_1nm
		tr3 = r3_1nm * media_1nm * wp_1nm

		# von Kries transform
		wptr1 = sum(tr1)
		wptr2 = sum(tr2)
		wptr3 = sum(tr3)
		tr1 /= wptr1
		if (wptr2 > 0): tr2 /= wptr2
		if (wptr3 > 0): tr3 /= wptr3
	if (args.s2tmode == "c"):
		# convert to chromaticity
		sr1_c = np.zeros(401)
		sr2_c = np.zeros(401)
		sr3_c = np.zeros(401)
		sr4_c = np.zeros(401)
		tr1_c = np.zeros(401)
		tr2_c = np.zeros(401)
		tr3_c = np.zeros(401)
		for i in range(401):
			stotal = sr1[i] + sr2[i] + sr3[i] + sr4[i]
			ttotal = tr1[i] + tr2[i] + tr3[i]
			if (stotal > 0):
				sr1_c[i] = sr1[i] / stotal
				sr2_c[i] = sr2[i] / stotal
				sr3_c[i] = sr3[i] / stotal
				sr4_c[i] = sr4[i] / stotal
			if (ttotal > 0):
				tr1_c[i] = tr1[i] / ttotal
				tr2_c[i] = tr2[i] / ttotal
				tr3_c[i] = tr3[i] / ttotal
	elif (args.s2tmode == "s"):
		# 0 - 65535
		smax = max(*sr1,*sr2,*sr3,*sr4)
		tmax = max(*tr1,*tr2,*tr3)
		sr1_c = sr1 * 65535 / smax
		sr2_c = sr2 * 65535 / smax
		sr3_c = sr3 * 65535 / smax
		sr4_c = sr4 * 65535 / smax
		tr1_c = tr1 * 65535 / smax
		tr2_c = tr2 * 65535 / smax
		tr3_c = tr3 * 65535 / smax

	# find best fit to target
	def source2target(xdata, r, g, b, u, rg=0, rb=0, ru=0, gb=0, gu=0, bu=0, rgb=0, rgu=0, rbu=0, gbu=0, r2=0, g2=0, b2=0, u2=0):
		first_order = r*sr1_c + g*sr2_c + b*sr3_c + u*sr4_c
		second_order = rg*(sr1_c*sr2_c) + rb*(sr1_c*sr3_c) + ru*(sr1_c*sr4_c) + gb*(sr2_c*sr3_c) + gu*(sr2_c*sr4_c) + bu*(sr3_c*sr4_c)
		third_order = rgb*(sr1_c*sr2_c*sr3_c) + rgu*(sr1_c*sr2_c*sr4_c) + rbu*(sr1_c*sr3_c*sr4_c) + gbu*(sr2_c*sr3_c*sr4_c) + r2*(sr1_c**2) + g2*(sr2_c**2) + b2*(sr3_c**2) + u2*(sr4_c**2)
		if (args.interactions == 3): ydata = first_order + second_order + third_order
		elif (args.interactions == 2): ydata = first_order + second_order
		else: ydata = first_order
		return(ydata)

	if (args.interactions == 1): p0 = [1,1,1,1]
	elif (args.interactions == 2): p0 = [1,1,1,1,0,0,0,0,0,0]
	elif (args.interactions == 3): p0 = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

	popt, pcov = scipy.optimize.curve_fit(source2target, x_1nm, tr1_c, p0=p0, bounds=[args.bound, np.inf])
	if (args.verbose): print(popt)
	s2t_r1 = popt
	popt, pcov = scipy.optimize.curve_fit(source2target, x_1nm, tr2_c, p0=p0, bounds=[args.bound, np.inf])
	if (args.verbose): print(popt)
	s2t_r2 = popt
	popt, pcov = scipy.optimize.curve_fit(source2target, x_1nm, tr3_c, p0=p0, bounds=[args.bound, np.inf])
	if (args.verbose): print(popt)
	s2t_r3 = popt
	
	if (args.s2tplot):
		plt.plot(x_1nm, tr1_c, color='k')
		plt.plot(x_1nm, source2target(x_1nm, *s2t_r1), '--', color='k')
		if (dimension > 1):
			plt.plot(x_1nm, tr2_c, color='k')
			plt.plot(x_1nm, source2target(x_1nm, *s2t_r2), '--', color='k')
			#plt.plot(x_1nm, sr1_c, ':r')
			#plt.plot(x_1nm, sr2_c, ':g')
			#plt.plot(x_1nm, sr1_c*sr2_c, ':y')
		if (dimension > 2):
			plt.plot(x_1nm, tr3_c, color='k')
			plt.plot(x_1nm, source2target(x_1nm, *s2t_r3), '--', color='k')
			#plt.plot(x_1nm, sr3_c, ':b')
			#plt.plot(x_1nm, sr1_c*sr3_c, ':m')
			#plt.plot(x_1nm, sr2_c*sr3_c, ':c')
		plt.show()

# gamma functions -- see https://en.wikipedia.org/wiki/Gamma_correction
# The default value is -1 and produces a standard sRGB gamma. Set to 1 to disable.
def gamma_decode(v): # decoding/linearization
	if (args.ingamma == -1):
		if v <= 0.04045: return np.float64(v / 12.92)
		else: return np.float64(math.pow((v + 0.055) / 1.055, 2.4))
	else:
		if (v <= 0): return 0
		return np.float64(math.pow(v, args.ingamma))

def gamma_encode(v): # encoding
	if (args.outgamma == -1):
		if v <= 0.0031308: return v * 12.92
		else: return np.float64(1.055 * math.pow(v, 1 / 2.4) - 0.055)
	else:
		if (v <= 0): return 0
		return np.float64(math.pow(v, 1 / args.outgamma))

"""
generate an RGB color from a spectral power distribution
This can produce either a "true color" using colormath or a "false color"
using the specified visual system. Out-of-gamut colors are scaled down if a
value is more than 1 and clipped if a value is less than 0.
"""
def spec2rgb(table, reflect=True, mode=args.fcmode):
	rgb_r = 0
	rgb_g = 0
	rgb_b = 0
	
	# interpolate
	if(len(table) < 401): table = np.interp(x_1nm, x_10nm, table)
	
	# combine SPD with illuminant
	if (reflect == True): table = table * wp_1nm
	
	# remove nans
	for i in range(table.shape[0]):
		if (np.isnan(table[i])): table[i] = 0
	
	if (mode == "lms"):
		# LMS
		for i in range(table.shape[0]):
			rgb_r += r1_1nm[i] * table[i] * media_1nm[i]
			rgb_g += r2_1nm[i] * table[i] * media_1nm[i]
			rgb_b += r3_1nm[i] * table[i] * media_1nm[i]
	
		# von Kries transform
		if (args.vonkries):
			wpr1 = 0
			wpr2 = 0
			wpr3 = 0
			for i in range(401):
				wpr1 += r1_1nm[i] * wp_1nm[i] * media_1nm[i]
				wpr2 += r2_1nm[i] * wp_1nm[i] * media_1nm[i]
				wpr3 += r3_1nm[i] * wp_1nm[i] * media_1nm[i]
		
			rgb_r = rgb_r / wpr1
			rgb_g = rgb_g / wpr2
			rgb_b = rgb_b / wpr3

		# gamma correction
		rgb_r = gamma_encode(rgb_r)
		rgb_g = gamma_encode(rgb_g)
		rgb_b = gamma_encode(rgb_b)

		# render
		if (dimension == 2): rgb_g = rgb_r
		rgb_tuple = (rgb_r, rgb_g, rgb_b)

	elif (mode == "cmf"):
		# CMF
		matrix_cmf = cmf()
		for i in range(401):
			if (dimension > 2):
				rgb_r += matrix_cmf[0][i] * table[i]
				rgb_g += matrix_cmf[1][i] * table[i]
				rgb_b += matrix_cmf[2][i] * table[i]
			elif (dimension == 2):
				rgb_r += matrix_cmf[0][i] * table[i]
				rgb_b += matrix_cmf[1][i] * table[i]
				rgb_g = rgb_r
		
		# gamma correction
		rgb_r = gamma_encode(rgb_r)
		rgb_g = gamma_encode(rgb_g)
		rgb_b = gamma_encode(rgb_b)

		# render
		rgb_tuple = (rgb_r, rgb_g, rgb_b)
	
	else:
		# colormath true color
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
		spec_830nm=0,
		illuminant='d65')
		
		# If an illuminant is set, we'll also set it here. Note that colormath
		# only adapts the CMFs and not the SPD itself, so we've already adjusted
		# it.
		if (args.white == 'd65'): spectral.set_illuminant('d65')
		elif (args.white == 'a' or args.white == 'i'): spectral.set_illuminant('a')
		elif (args.white == 'e'): spectral.set_illuminant('e')
		
		srgb = convert_color(spectral, sRGBColor)
		rgb_tuple = srgb.get_value_tuple()
		rgb_r = rgb_tuple[0]
		rgb_g = rgb_tuple[1]
		rgb_b = rgb_tuple[2]
	
	if (args.verbose):
		print("Mode: " + mode)
		print("RGB: " + str(rgb_tuple))

	rgb_max = max(*rgb_tuple,1)
	rgb_min = min(*rgb_tuple)
	if ((rgb_max > 1 or rgb_min < 0) and args.warnings):
		print("Warning (spec2rgb): color is out of gamut. For details, use --verbose/-v.")
	rgb_tuple = (max(rgb_r/rgb_max, 0), max(rgb_g/rgb_max, 0), max(rgb_b/rgb_max, 0))

	return(rgb_tuple)

"""
contrast sensitivity (Vorobyev and Osorio (1998); Vorobyev and Osorio (2001))
Without quantum noise, this tends to make unrealistic predictions when the value
for one or more cone signals is very small.
"""
def color_contrast(table1, table2, quantum_noise=args.qn, r1 = r1_1nm, r2 = r2_1nm, r3 = r3_1nm, r4 = r4_1nm):
	if (dimension > 4):
		print("Warning: color_contrast() does not support dimensions higher than 4. Any "
		+ "receptors after r4 will be ignored.")
	# interpolate 1-nm intervals if we're provided with 10-nm
	if (len(table1) < 401): table1 = np.interp(x_1nm, x_10nm, table1)
	if (len(table2) < 401): table2 = np.interp(x_1nm, x_10nm, table2)
	
	# background light
	wpr1 = 0
	wpr2 = 0
	wpr3 = 0
	wpr4 = 0
	
	qr11 = 0
	qr21 = 0
	qr31 = 0
	qr41 = 0
	qr12 = 0
	qr22 = 0
	qr32 = 0
	qr42 = 0
	for i in range(table1.shape[0]):
		if (table1[i] > 0): # zero/nan check
			qr11 += r1[i] * table1[i] * wp_1nm[i] * media_1nm[i]
			qr21 += r2[i] * table1[i] * wp_1nm[i] * media_1nm[i]
			qr31 += r3[i] * table1[i] * wp_1nm[i] * media_1nm[i]
			qr41 += r2[i] * table1[i] * wp_1nm[i] * media_1nm[i]
		wpr1 += r1[i] * wp_1nm[i] * media_1nm[i]
		wpr2 += r2[i] * wp_1nm[i] * media_1nm[i]
		wpr3 += r3[i] * wp_1nm[i] * media_1nm[i]
		wpr4 += r2[i] * wp_1nm[i] * media_1nm[i]
	for i in range(table2.shape[0]):
		w = i + 300
		if (table2[i] > 0):
			qr12 += r1[i] * table2[i] * wp_1nm[i] * media_1nm[i]
			qr22 += r2[i] * table2[i] * wp_1nm[i] * media_1nm[i]
			qr32 += r3[i] * table2[i] * wp_1nm[i] * media_1nm[i]
			qr42 += r2[i] * table2[i] * wp_1nm[i] * media_1nm[i]
	
	# normalize
	if (args.vonkries):
		qr11 = qr11 / wpr1
		qr21 = qr21 / wpr2
		qr31 = qr31 / wpr3
		qr41 = qr41 / wpr4
		qr12 = qr12 / wpr1
		qr22 = qr22 / wpr2
		qr32 = qr32 / wpr3
		qr42 = qr42 / wpr4
	
	# differences
	dfr1 = math.log(qr11 / qr12)
	dfr2 = math.log(qr21 / qr22)
	dfr3 = math.log(qr31 / qr32)
	dfr4 = math.log(qr41 / qr42)
	
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
		# area of retina: 1094 mm^2 (human), 80 mm^2 (rat; Mayhew & Astle 1997)
		# total number of cones: ~6 million (human)
		# outer segment length: 30 um (chicken), 6.2 um (cat; Fisher, Pfeffer & Anderson
		# 1983), 7 um (gray squirrel; "), 20 um (honey possum), 25.1-40 um (human),
		# 15.2 um (goldfish red single cones, Harosi & Flamarique 2012), 5 um (tammar
		# wallaby and fat-tailed dunnart; Ebeling, Natoli & Hemmi 2010)
		# outer segment width: 1.73 um (chicken), 6.1 um (goldfish red single cones)
		# optical density: 0.015 (coral reef triggerfish, birds)
		# focal length: 8300 um = 8.3 mm (chickens), 2.6 mm (mice; Geng et al. 2011),
		# 22.3 mm (human; "), 2.5 mm (triggerfish)
		# pupil size: 4900-3500 um = 4.9-3.5 mm (chicken), 2 mm (mouse), 6 mm (human)
		# transmittance of ocular media: 80% (general estimate from The Optics of Light)
		# integration time/flicker fusion: 50-12 ms = 20-83 Hz (chicken),
		# 70-80 Hz = 14-12 ms (cat and dog), 18 Hz = 56 ms (rat; Gil, Valente & Shemesh 2024),
		# 14 Hz = 71 ms (mice; "), 60 Hz = 17 ms (human), 25 Hz = 40 ms (some fish)
		
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
		#T = 0.8 / ocular_media(700)
		#T = 1 / ocular_media(700) # set to 1 at 700 nm
		T = 1
	
		aqr11 = 0
		aqr21 = 0
		aqr31 = 0
		aqr41 = 0
		aqr12 = 0
		aqr22 = 0
		aqr32 = 0
		aqr42 = 0
		for i in range(table1.shape[0]):
			w = i + 300
			if (table1[i] > 0):
				aqr11 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r1[i]*l)) * r1[i] * table1[i] * wp_photons[i] * media_1nm[i]
				aqr21 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r2[i]*l)) * r2[i] * table1[i] * wp_photons[i] * media_1nm[i]
				aqr31 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r3[i]*l)) * r3[i] * table1[i] * wp_photons[i] * media_1nm[i]
				aqr41 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r4[i]*l)) * r4[i] * table1[i] * wp_photons[i] * media_1nm[i]
		for i in range(table2.shape[0]):
			w = i + 300
			if (table2[i] > 0):
				aqr12 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r1[i]*l)) * r1[i] * table2[i] * wp_photons[i] * media_1nm[i]
				aqr22 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r2[i]*l)) * r2[i] * table2[i] * wp_photons[i] * media_1nm[i]
				aqr32 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r3[i]*l)) * r3[i] * table2[i] * wp_photons[i] * media_1nm[i]
				aqr42 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*r4[i]*l)) * r4[i] * table2[i] * wp_photons[i] * media_1nm[i]
		
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
			aqr11 = aqr11*v*args.r1p
			aqr12 = aqr12*v*args.r1p
			aqr21 = aqr21*v*args.r2p
			aqr22 = aqr22*v*args.r2p
			aqr31 = aqr31*v*args.r3p
			aqr32 = aqr32*v*args.r3p
			aqr41 = aqr41*v*args.r4p
			aqr42 = aqr42*v*args.r4p
		
		er1 = math.sqrt((1 / aqr11 + 1 / aqr12) + 2*wr1**2)
		er2 = math.sqrt((1 / aqr21 + 1 / aqr22) + 2*wr2**2)
		er3 = math.sqrt((1 / aqr31 + 1 / aqr32) + 2*wr3**2)
		er4 = math.sqrt((1 / aqr41 + 1 / aqr42) + 2*wr4**2)
		
		if (args.qcheck):
			print("r11: " + str(aqr11) + ", r21: " + str(aqr21)
			+ ", r31: " + str(aqr31) + ", r41: " + str(aqr41))
			print("r12: " + str(aqr12) + ", r22: " + str(aqr22)
			+ ", r32: " + str(aqr32) + ", r42: " + str(aqr42))
			print("er1: " + str(er1) + ", er2: " + str(er2)
			+ ", er3: " + str(er3) + ", er4: " + str(er4))
			print("wr1: " + str(wr1) + ", wr2: " + str(wr2)
			+ ", wr3: " + str(wr3) + ", wr4: " + str(wr4))
	else:
		er1 = wr1
		er2 = wr2
		er3 = wr3
		er4 = wr4
	
	if (dimension > 3):
		delta_s = ((er3*er4)**2*(dfr1 - dfr2)**2
		+ (er2*er4)**2*(dfr1 - dfr3)**2
		+ (er2*er3)**2*(dfr1 - dfr4)**2
		+ (er1*er4)**2*(dfr2 - dfr3)**2
		+ (er1*er3)**2*(dfr2 - dfr4)**2
		+ (er1*er2)**2*(dfr3 - dfr4)**2) / ((er1*er2*er3)**2 + (er1*er2*er4)**2 + (er1*er3*er4)**2 + (er2*er3*er4)**2)
	elif (dimension == 3):
		delta_s = (er3**2*(dfr1 - dfr2)**2
		+ er2**2*(dfr1 - dfr3)**2
		+ er1**2*(dfr3 - dfr2)**2) / ((er1*er2)**2 + (er1*er3)**2 + (er2*er3)**2)
	else:
		delta_s = (dfr1 - dfr2)**2 / (er1**2 + er2**2)
	return math.sqrt(delta_s)

# brightness contrast based on https://journals.biologists.com/jeb/article/207/14/2471/14748/Interspecific-and-intraspecific-views-of-color
# I've changed it to a signed value because we care about the direction.
def brightness_contrast(table1, table2, r1=luminosity):
	# interpolate 1-nm intervals if we're provided with 10-nm
	if (len(table1) < 401): table1 = np.interp(x_1nm, x_10nm, table1)
	if (len(table2) < 401): table2 = np.interp(x_1nm, x_10nm, table2)
	
	# We're redoing this to use sensitivity() so we can use the photopic luminosity
	# function for humans.
	q1 = 0
	q2 = 0
	for i in range(401):
		q1 += r1[i] * table1[i] * wp_1nm[i]
		q2 += r1[i] * table2[i] * wp_1nm[i]
	
	df = math.log(q1 / q2)
	delta_s = df / args.weberb # not an absolute value
	return delta_s

"""
plot visual pigment templates
I changed the line colors back to red, green and blue. I've been trying to avoid this kind
of color coding, but in this case it's generally understood which visual pigment template
is which based on its position on the graph, and the shades-of-gray "color" scheme was too
pale for my liking.

09/10/2025 -- you can now plot an arbitrary number of templates and the color for each is
chosen based on their peak sensitivity
"""
if (args.vpt):
	# rod
	if (args.rod != 0):
		yvalues = np.empty(401)
		ymax = 0
		for i in range(401):
			x_1nm[i] = i + 300
			yvalues[i] = vpt(i+300, rod) * media_1nm[i]
			ymax = max(ymax, yvalues[i])
		plt.plot(x_1nm, yvalues/ymax, ':k')
	
	# non-rod receptors
	for i in range(dimension):
		yvalues = rlist[i] * media_1nm
		ymax = max(*yvalues)
		
		plt.plot(x_1nm, yvalues/ymax, color='k')

	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.show()

# show luminous efficiency function and lens filtering
if (args.luminosity):
	ms = 300
	msy = 0
	for i in range(401):
		if (luminosity[i] > msy):
			msy = luminosity[i]
			ms = i+300
	print("Maximum sensitivity: " + str(ms))
	print("")
	
	plt.plot(x_1nm, luminosity, 'k')
	if (args.media != "none"):
		plt.plot(x_1nm, media_1nm, ':k')
	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.show()

"""
plot color triangle
See this paper for information on Maxwell triangles: https://www.researchgate.net/publication/8064482_Animal_colour_vision_-_Behavioural_tests_and_physiological_concepts
In the literature, Maxwell triangles may have either M on the left and S on top, as in the above
paper, or the reverse (S left, M top). triangle() defaults to the former because it simplifies
having multiple dimensions of plots. If you want a different orientation, you can specify another
order for --receptors and change the labels with --r1name, --r2name and --r3name. Note the default
label for the left side is "M" for both lines and triangles.

The name "triangle" is now slightly misleading because it produces a line or tetrahedron when
provided with 2 or 4+ receptors.

Numbers are off by default because they overlap with everything and are almost unreadable.

Axes are off by default because these plots don't usually include them and the Y-axis isn't
meaningful for dichromacy. This doesn't hide the color bar for the 4th dimension.
"""
def triangle(spectra=[], reflect=True, markers=[], colors=[], text=[], legend=False, numbers=False, mec='k', gamut=False, gamutcolor='k', gamutedge='-', achro=True, axes=False):
	# default list elements
	# Checking if they're blank doesn't always work if you call the function twice in a row.
	if (len(markers) < len(spectra)):
		for i in range(len(spectra)): markers.append('o')
	if (len(colors) < len(spectra)):
		for i in range(len(spectra)): colors.append('k')
	if (len(text) < len(spectra)):
		for i in range(len(spectra)): text.append('')
	
	# plot triangle and visible gamut
	wr1 = 0
	wr2 = 0
	wr3 = 0
	wr4 = 0
	wr5 = 0
	xvalues = np.empty(41)
	yvalues = np.zeros(41)
	zvalues = np.zeros(41)
	wvalues = np.zeros(41)
	labels = np.empty(41)
	
	# triangle
	if (dimension > 3):
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		
		# vertices
		r1v = [1/math.sqrt(2), -math.sqrt(2)/(2*math.sqrt(3)), -1/(2*math.sqrt(3))]
		r2v = [-1/math.sqrt(2), -math.sqrt(2)/(2*math.sqrt(3)), -1/(2*math.sqrt(3))]
		r3v = [0, math.sqrt(2/3), -1/(2*math.sqrt(3))]
		r4v = [0, 0, math.sqrt(3)/2]
		
		# LMS triangle
		ax.plot([r1v[0], r2v[0], r3v[0], r1v[0]], [r1v[1], r2v[1], r3v[1], r1v[1]],
		[r1v[2], r2v[2], r3v[2], r1v[2]], '-k')
		# LMU triangle
		ax.plot([r1v[0], r2v[0], r4v[0], r1v[0]], [r1v[1], r2v[1], r4v[1], r1v[1]],
		[r1v[2], r2v[2], r4v[2], r1v[2]], '-k')
		# LSU triangle
		ax.plot([r1v[0], r3v[0], r4v[0], r1v[0]], [r1v[1], r3v[1], r4v[1], r1v[1]],
		[r1v[2], r3v[2], r4v[2], r1v[2]], '-k')
		# MSU triangle
		ax.plot([r2v[0], r3v[0], r4v[0], r2v[0]], [r2v[1], r3v[1], r4v[1], r2v[1]],
		[r2v[2], r3v[2], r4v[2], r2v[2]], '-k')
		
	elif (dimension == 3):
		xborder = np.array([-math.sqrt(1/2), 0, math.sqrt(1/2), -math.sqrt(1/2)])
		yborder = np.array([-math.sqrt(2/3)/2, math.sqrt(2/3), -math.sqrt(2/3)/2, -math.sqrt(2/3)/2])
		plt.plot(xborder, yborder, '-k')
		plt.text(-math.sqrt(1/2) - 0.05, -math.sqrt(2/3)/2 - 0.025, args.r2name)
		plt.text(0 - 0.025, math.sqrt(2/3) + 0.0125, args.r3name)
		plt.text(math.sqrt(1/2) + 0.0125, -math.sqrt(2/3)/2 - 0.025, args.r1name)
	else:
		plt.plot([-1, 1], [0, 0], '-k')
		plt.text(-1.1, -.025, args.r2name)
		plt.text(1.1, -.025, args.r1name)
		plt.ylim(-1, 1)

	# gamut
	for i in range(41):
		w = i*10 + 300
		wr1 += r1_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		wr2 += r2_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		wr3 += r3_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		wr4 += r4_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		wr5 += r5_1nm[i*10] * wp_10nm[i] * media_10nm[i]
	for i in range(41):
		w = i*10 + 300
		labels[i] = w
		r1 = r1_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		r2 = r2_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		r3 = r3_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		r4 = r4_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		r5 = r5_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		if (args.vonkries):
			r1 = r1 / wr1
			if (wr2 > 0): r2 = r2 / wr2
			if (wr3 > 0): r3 = r3 / wr3
			if (wr4 > 0): r4 = r4 / wr4
			if (wr5 > 0): r5 = r5 / wr5
		total = r1 + r2 + r3 + r4 + r5
		if (total > 0): # avoid "invalid value encountered in scalar divide" warning
			r1 = r1 / total
			r2 = r2 / total
			r3 = r3 / total
			r4 = r4 / total
			r5 = r5 / total
		if (dimension > 2):
			xvalues[i] = math.sqrt(1/2)*(r1 - r2)
			yvalues[i] = math.sqrt(2/3)*(r3 - (r1 + r2)/2)
			if (dimension > 3): zvalues[i] = math.sqrt(3)/2*(r4 - (r1 + r2 + r3)/3)
			# we're not combining this with the others, so no coefficient
			if (dimension > 4): wvalues[i] = (r5 - (r1 + r2 + r3 + r4)/4)
		else:
			if (total > 0): xvalues[i] = (r1 - r2)
			else: xvalues[i] = 0
		if (total > 0): # Hide zeroes. Line segments connecting them are still plotted.
			if (dimension > 3):
				if (gamut): ax.plot(xvalues[i], yvalues[i], zvalues[i], 'o', color=gamutcolor)
				# use int() to cut off the ".0" at the end
				if (numbers): ax.text(xvalues[i], yvalues[i], zvalues[i], int(labels[i]))
			else:
				if (gamut): plt.plot(xvalues[i], yvalues[i], 'o', color=gamutcolor)
				if (numbers): plt.text(xvalues[i], yvalues[i], int(labels[i]))
	
	if (gamut):
		if (dimension > 3): plt.plot(xvalues, yvalues, zvalues, linestyle=gamutedge, color=gamutcolor)
		else: plt.plot(xvalues, yvalues, linestyle=gamutedge, color=gamutcolor)
	
	# plot point at origin
	if (achro):
		if (dimension > 3): plt.plot(0, 0, 0, linestyle='', color='0.8', marker='s')
		else: plt.plot(0, 0, linestyle='', color='0.8', marker='s')
	
	
	# plot spectra
	if (len(spectra) > 0):
		wpr1 = 0
		wpr2 = 0
		wpr3 = 0
		wpr4 = 0
		wpr5 = 0
		posxlist = []
		posylist = []
		poszlist = []
		poswlist = []
		for i in range(401):
			wpr1 += r1_1nm[i] * wp_1nm[i] * media_1nm[i]
			wpr2 += r2_1nm[i] * wp_1nm[i] * media_1nm[i]
			wpr3 += r3_1nm[i] * wp_1nm[i] * media_1nm[i]
			wpr4 += r4_1nm[i] * wp_1nm[i] * media_1nm[i]
			wpr5 += r5_1nm[i] * wp_1nm[i] * media_1nm[i]
		
		for i in range(len(spectra)):
			table = spectra[i]
			table_r1 = 0
			table_r2 = 0
			table_r3 = 0
			table_r4 = 0
			table_r5 = 0
			
			# interpolate
			if(len(table) < 401): table = np.interp(x_1nm, x_10nm, table)
			
			# combine SPD with illuminant
			if (reflect == True):
				for j in range(401): table[j] *= wp_1nm[j]
			
			# remove nans
			for j in range(401):
				if (np.isnan(table[j])): table[j] = 0
			
			for j in range(401):
				table_r1 += r1_1nm[j] * table[j] * media_1nm[j]
				table_r2 += r2_1nm[j] * table[j] * media_1nm[j]
				table_r3 += r3_1nm[j] * table[j] * media_1nm[j]
				table_r4 += r4_1nm[j] * table[j] * media_1nm[j]
				table_r5 += r5_1nm[j] * table[j] * media_1nm[j]
			
			# von Kries transform
			if (args.vonkries):
				table_r1 = table_r1 / wpr1
				if (wpr2 > 0): table_r2 = table_r2 / wpr2
				if (wpr3 > 0): table_r3 = table_r3 / wpr3
				if (wpr4 > 0): table_r4 = table_r4 / wpr4
				if (wpr5 > 0): table_r5 = table_r5 / wpr5

			# declare variables
			posx = 0
			posy = 0
			posz = 0
			posw = 0

			if (dimension > 3):
				total = table_r1 + table_r2 + table_r3 + table_r4 + table_r5
				if (total > 0):
					posx = math.sqrt(1/2)*(table_r1 - table_r2) / total
					posy = math.sqrt(2/3)*(table_r3 - (table_r1 + table_r2)/2) / total
					posz = math.sqrt(3)/2*(table_r4 - (table_r1 + table_r2 + table_r3)/3) / total
					posw = (table_r5 - (table_r1 + table_r2 + table_r3 + table_r4)/4) / total
				else:
					posx = 0
					posy = 0
					posz = 0
					posw = 0
				if (args.verbose):
					print("Relative quantum catches: r1=" + str(table_r1) + ", r2=" + str(table_r2) + ", r3=" + str(table_r3) + ", r4=" + str(table_r4) + ", r5=" + str(table_r5))
					print("Position in color tetrahedron: " + str((posx, posy, posz, posw)))
			elif (dimension == 3):
				total = table_r1 + table_r2 + table_r3
				if (total > 0):
					posx = math.sqrt(1/2)*(table_r1 - table_r2) / total
					posy = math.sqrt(2/3)*(table_r3 - (table_r1 + table_r2)/2) / total
				else:
					posx = 0
					posy = 0
				if (args.verbose):
					print("Relative quantum catches: r1=" + str(table_r1) + ", r2=" + str(table_r2) + ", r3=" + str(table_r3))
					print("Position in color triangle: " + str((posx, posy)))
			elif (dimension == 2):
				if (table_r1 + table_r2 > 0): posx = (table_r1 - table_r2) / (table_r1 + table_r2)
				else: posx = 0
				if (args.verbose):
					print("Relative quantum catches: r1=" + str(table_r1) + ", r2=" + str(table_r2))
					print("Position on color line: " + str(posx))
			elif (dimension == 1 and args.verbose):
				# monochromacy
				print("Relative quantum catch: " + str(table_r1))

			if (dimension == 4):
				if (text[i] == ''):
					ax.plot(posx, posy, posz, marker=markers[i], mec=mec, color=colors[i], linestyle='')
				else:
					ax.plot(posx, posy, posz, marker=markers[i], mec=mec, color=colors[i], label=text[i], linestyle='')
			elif (dimension < 4):
				if (text[i] == ''):
					plt.plot(posx, posy, marker=markers[i], mec=mec, color=colors[i], linestyle='')
				else:
					plt.plot(posx, posy, marker=markers[i], mec=mec, color=colors[i], label=text[i], linestyle='')

			posxlist.append(posx)
			posylist.append(posy)
			poszlist.append(posz)
			poswlist.append(posw)
	
		# cmap needs the whole list at once
		if (dimension > 4):
			img = ax.scatter(posxlist, posylist, poszlist, c=poswlist, cmap=plt.viridis())
			fig.colorbar(img)
	
	if (dimension > 3):
		# positioning text is awkward, I have a better idea...
		ax.plot(*r1v, 'o', color='r')
		ax.plot(*r2v, 'o', color='g')
		ax.plot(*r3v, 'o', color='b')
		ax.plot(*r4v, 'o', color='purple')
	
	if (axes == False): plt.axis('off')
	if (legend and dimension < 5): plt.legend()
	plt.show()

if (type(args.csplot) == list):
	if (args.csplot == []):
		triangle(numbers=True, gamut=True)
	else:
		spectra = []
		colors = []
		for i in range(len(args.csplot)):
			spectrum = csv2spec(args.csplot[i])
			spectra.append(spectrum)
			colors.append(spec2rgb(spectrum))
		triangle(spectra, colors=colors)

"""
function for modeling color discrimination trials

This outputs a probability for both success (defined as at least the criterion) and failure (less
than the criterion) because the behavioral results include both. We use different thresholds
for "success" for the different sets of discriminations. In Friedman (1967) the criterion was
35 out of 40 for the color pairs and 68 of 80 for the color/gray pairs. In Gutierrez et al. (2011)
significant performance is defined as more than 62.5% and the total number of trials for each
pair was 64, so the criterion would be 41 of 64. This produces considerably different
significance thresholds when comparing "success" or "failure" against expected values:

* 35/40 and 68/80: success 0-12 (out of 16), failure 15-16
* 41/64: success 0-8, failure 12-16

For contrast where we don't consider overlap and expect chance (50%) performance for deltaS < 1,
this means that under the first criterion "success" is defined as at least 10 of 16 distinguishable
pairs and "failure" is defined as at most 13 of 16, whereas under the second criterion "success"
is at least 1 and "failure" is at most 7. There are also several intermediate values where
the probability of both success and failure exceeds 0.05. With the first two criteria, success
becomes more likely than failure at 14 (12 distinguishable), and with the second this occurs at 10.5
(5 distinguishable).

This is now also used directly for optimization graphs, so all the parts relevant only to
testing one set of pigments have been made conditional to reduce noise and computation.
"""
def color_disc(first, second, correct=35, trials=40, r1 = 0, r2 = 0, r3 = 0, output=True):
	# custom receptors
	custom_r1 = r1_1nm
	custom_r2 = r2_1nm
	custom_r3 = r3_1nm
	if (r1 > 0):
		for i in range(401): custom_r1[i] = vpt(i+300, r1)
	if (r2 > 0):
		for i in range(401): custom_r2[i] = vpt(i+300, r1)
	if (r3 > 0):
		for i in range(401): custom_r3[i] = vpt(i+300, r3)
	plt.show()
	d = 0 # different
	s = 0 # same
	counter = 0
	box = np.empty(16)
	
	for i in range(4):
		for j in range(4):
			contrast = color_contrast(first[i], second[j], r1 = custom_r1, r2 = custom_r2, r3 = custom_r3)
			box[counter] = contrast
			if (output):
				print(str(i) + " vs. " + str(j) + ": " + str(contrast))
				if (contrast < 1): s += 1
				else: d += 1
			counter += 1
	if (output):
		print("Different: " + str(d))
		print("Same: " + str(s))
		binomial = binomtest(correct, trials, p=(d + s/2)/16, alternative='greater')
		binomial2 = binomtest(correct-1, trials, p=(d + s/2)/16, alternative='less')
		print("P-value (success): " + str(binomial.pvalue))
		print("P-value (failure): " + str(binomial2.pvalue))
	
	# median contrast
	mc = statistics.median(box)
	minc = min(*box)
	maxc = max(*box)
	if (output):
		print("Median contrast: " + str(mc))
		print("")
	return [mc, minc, maxc, box]


"""
like color_disc except for brightness
We want to know both the direction and the strength of the difference
to judge whether it could be reliably used as a cue.
"""
def brightness_disc(first, second, correct=35, trials=40, r1=0, output=True, customize=False):
	# custom L cone
	custom_r1 = luminosity
	if (r1 > 0):
		for i in range(401): custom_r1[i] = vpt(i+300, r1) * media_1nm[i]
	fb = 0 # first brighter
	sb = 0 # second brighter
	# first set
	for i in range(4):
		# second set
		for j in range(4):
			contrast = brightness_contrast(first[i], second[j], r1=custom_r1)
			if (contrast >= 1): fb += 1
			elif (contrast <= -1): sb += 1
	equal = 16 - fb - sb
	
	binomial = binomtest(correct, trials, (max(fb, sb) + equal/2) / 16, alternative='greater')
	binomial2 = binomtest(correct-1, trials, (max(fb, sb) + equal/2) / 16, alternative='less')
	if (output):
		print("First brighter: " + str(fb))
		print("Second brighter: " + str(sb))
		print("Equal: " + str(equal))
		print("Expected % correct: " + str(100*(max(fb, sb) + equal/2) / 16) + "% (" + str(100*max(fb, sb) / 16) + "-" + str(100*(max(fb, sb) + equal) / 16) + "%)")
		print("P-value (success): " + str(binomial.pvalue))
		print("P-value (failure): " + str(binomial2.pvalue))
		print("")

	return [binomial.pvalue, binomial2.pvalue]

# print execution time
#print("%s seconds" % (time.time() - start_time))
