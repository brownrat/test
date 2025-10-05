"""
This file contains the following:
* list of arguments
* visual pigments: standard visual pigment template (vpt()) and human cone fundamentals
* ocular media transmittance for various species
* blackbody: blackbody()
* color/brightness analysis: color_contrast(), brightness_contrast() and spectral_rendering()
* illuminants: D65, E, A, custom incandescent similar to A
* color matching functions
* Maxwell triangle color space plots
* color discrimination: color_disc()

Other stuff has been spun off to the following files, which import this one. I'd like to add additional arguments to those but
don't see how to do that.
* erg.py: analysis of ERG data from Jacobs & Williams (2010)
* illuminants.py: test illuminants
* step.py: step function spectra
* peak.py: Gaussian peak spectra
* ramp.py: ramp function spectra
* sky.py: sky/daylight spectra
* wratten.py: Kodak Wratten filters
* munsell.py: Munsell color cards

wratten.py and munsell.py just produce plots by default. To see the color rendering and trial modeling, you have to use
--wratten/--munsell. This is because they also contain other features that are accessed with other arguments (see the list).
"""

# this list doesn't recursively import
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

# arguments
# 07/03/2025 -- replaced some default values:
# * weberb: 0.14 -> 0.11. 0.11 is the value determined for humans in a study of manatees (https://pubmed.ncbi.nlm.nih.gov/9202447/).
# 0.14 was found in a previous study (Cornsweet & Pinsker), but the number seems to only appear in citations
# and not in the actual paper. All I can find in the paper is graphs.
# * lp/mp: 9/4.5 -> 10/5. The ratio of L to S cones in Didelphis aurita and other marsupials is closer to 10:1 than
# 9:1. The L:M ratio is kept at 2:1.
# * ls/ms: 0.68990272/0.34832189 -> 1/0. The former values are for humans. Since we're now modeling human luminance vision with the CIE
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
parser.add_argument("--media", type=str, default="none", help="ocular media transmittance")
parser.add_argument("--qn", help="use quantum noise in color differences", action="store_true")
parser.add_argument("--ls", type=float, default=1, help="contribution of L cones to spectral sensitivity") # formerly 0.68990272
parser.add_argument("--ms", type=float, default=0, help="contribution of M cones to spectral sensitivity")
parser.add_argument("--ss", type=float, default=0, help="contribution of S cones to spectral sensitivity")
parser.add_argument("--rb", type=float, default=0, help="contribution of rods to spectral sensitivity")
parser.add_argument("--novk", help="disable von Kries transform", action="store_true")
parser.add_argument("--luminosity", help="show wavelength with maximum sensitivity", action="store_true")
parser.add_argument("--white", default="e", help="reference white")
parser.add_argument("--wratten", help="Kodak Wratten camera filters", action="store_true")
parser.add_argument("--wv2", help="use tabulated Kodak Wratten data from Kodak Photographic Filters Handbook (1990)", action="store_true")
parser.add_argument("--wcheck", help="check camera filter brightness matches", action="store_true")
parser.add_argument("--wopt1", help="test varying rod spectral sensitivity", action="store_true")
parser.add_argument("--render", help="show colormath spectral rendering for comparison", action="store_true")
parser.add_argument("--qcheck", help="show absolute quantum catches", action="store_true")
parser.add_argument("--vpt", help="plot visual pigment templates", action="store_true")
parser.add_argument("--triangle", help="plot color triangle", action="store_true")
parser.add_argument("--blackbody", type=int, default=0, help="plot black body curve for a given temperature in Kelvin")
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
parser.add_argument("-r", "--red", type=float, default=700, help="'red' primary")
parser.add_argument("-g", "--green", type=float, default=546.1, help="'green' primary")
parser.add_argument("-b", "--blue", type=float, default=435.8, help="'blue' primary")
parser.add_argument("--q10deg", help="use 10-degree instead of 2-degree cone fundamentals", action="store_true")
parser.add_argument("--e2deg", help="use 2-degree energy cone fundamentals", action="store_true")
parser.add_argument("-v", "--verbose", help="", action="store_true")
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

# types of ocular media

# human lens approximation (Lamb 1995)
# 21/03/2025 -- optical density uses 10 not e, again.
# We don't need this for modeling a standard human because we have the CIE luminosity function
# and human cone fundamentals. I'm only leaving it in for customization.
def human_filter(w):
	try:
		density = 1.1*math.exp((400 - w) / 15) + 0.11*math.exp((500 - w) / 80)
	except OverflowError:
		if (args.warnings): print("Warning (template_filter): math overflow, returning 0")
		return 0
	if (density < 0):
		return 0
	return 10**(-density)

# mouse (Jacobs & Williams 2007)
# Transmission below 310 nm is not provided and has been set to 0. This means we're effectively
# only considering wavelengths >=310, but it's easier to keep track of the numbers if they begin
# with a multiple of 100.
mouse_lens = np.array([
	0, # 300
	0, #309 (this is here so we don't interpolate 301-309)
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

# Thylamys elegans lens (Palacios et al. 2010, fig. 7)
# As with the ERG data this is a dot plot, apparently with intervals of 4 nm (but the
# text says measurements were made every 20 nm). Relative to the value at 700 nm,
# the 50% transmission wavelength is between 372 and 376.
thylamys_lens = np.array([
	[300,1.99085],
	[303.883,2.84897],
	[308.091,3.02060],
	[311.974,3.36384],
	[316.181,4.56522],
	[320.065,5.28604],
	[324.272,6.00687],
	[328.155,6.76201],
	[332.362,7.58581],
	[336.246,8.64989],
	[340.129,8.68421],
	[344.337,9.61098],
	[348.220,10.0915],
	[352.104,10.7094],
	[356.311,11.1899],
	[360.194,12.1167],
	[364.078,12.2540],
	[367.961,12.4600],
	[372.006,12.9405],
	[376.133,13.3009],
	[380.016,13.6098],
	[383.900,14.1505],
	[388.026,14.4079],
	[392.152,14.8970],
	[396.036,14.9485],
	[399.919,15.2574],
	[404.045,15.8238],
	[407.929,15.8753],
	[412.055,15.9010],
	[416.181,16.6476],
	[420.065,17.2912],
	[423.948,17.7031],
	[428.074,17.4199],
	[432.201,17.8061],
	[436.084,17.8318],
	[439.968,18.4239],
	[444.094,18.7843],
	[447.977,18.8616],
	[452.104,18.8358],
	[455.987,18.6556],
	[460.113,19.0675],
	[463.997,18.7071],
	[467.880,18.8873],
	[472.006,18.9645],
	[476.133,19.1705],
	[480.016,19.2477],
	[484.142,20.6121],
	[487.783,20.5606],
	[491.909,20.2002],
	[496.036,20.2260],
	[500.162,19.9943],
	[504.045,20.1745],
	[507.929,20.5092],
	[512.055,20.6894],
	[516.181,20.9468],
	[520.065,21.1785],
	[523.948,21.0498],
	[528.074,22.2855],
	[531.958,22.1567],
	[536.084,21.8993],
	[540.210,22.5429],
	[544.094,22.5944],
	[547.977,22.0023],
	[552.104,21.9251],
	[555.987,22.0023],
	[559.871,21.9508],
	[564.240,22.5686],
	[568.123,22.4399],
	[572.006,22.9548],
	[575.890,22.9805],
	[580.016,22.5429],
	[583.900,22.5429],
	[588.026,22.9548],
	[591.909,23.0578],
	[596.278,23.7786],
	[599.919,23.3410],
	[604.045,24.1133],
	[607.929,24.0875],
	[612.055,23.8816],
	[616.181,23.4182],
	[620.065,23.5984],
	[624.191,23.9331],
	[628.074,24.1648],
	[632.201,24.1133],
	[636.084,24.3192],
	[639.968,24.2162],
	[644.094,24.1905],
	[647.977,24.3192],
	[651.618,24.6281],
	[655.744,24.9113],
	[659.871,25.4519],
	[663.754,24.7569],
	[667.880,24.8341],
	[672.006,24.9113],
	[675.890,24.9886],
	[679.773,25.2975],
	[683.900,24.8084],
	[688.026,25.0915],
	[691.909,25.1430],
	[696.036,25.7094],
	[700.162,26.4302],
])

# brushtail possum whole eye (Vlahos 2012, fig. II.7)
# This is slightly different from the transmission curve for just the lens. I didn't
# use that because they overlap in the graph and the whole eye curve (blue) is
# on top. Since it's a continuous curve rather than a dot plot, I again made the
# arbitrary choice to sample it at 5-nm intervals. The wavelength of 50% transmission
# is at about 370, as reported.
brushtail_eye = np.array([
	[300.000,0.330716],
	[305.186,0.800026],
	[310.084,2.35957],
	[314.881,6.87857],
	[319.678,11.0341],
	[325.246,13.7362],
	[330.431,15.8151],
	[334.847,17.7898],
	[340.223,20.3879],
	[345.405,24.3359],
	[349.243,27.3489],
	[354.615,33.1661],
	[360.178,39.3988],
	[365.358,45.4237],
	[369.960,51.8637],
	[375.330,59.1347],
	[380.317,65.1595],
	[385.304,70.8728],
	[389.717,75.3397],
	[394.708,78.5608],
	[399.892,81.5742],
	[405.460,83.8609],
	[410.453,85.4205],
	[415.062,86.7722],
	[420.440,88.2281],
	[425.241,89.3723],
	[429.659,90.5163],
	[435.037,91.3491],
	[439.839,91.7663],
	[445.410,92.1839],
	[450.022,91.9780],
	[455.211,90.8377],
	[459.631,89.6972],
	[464.628,88.4530],
	[470.394,86.7938],
	[475.390,86.0688],
	[480.001,85.5514],
	[484.805,85.1379],
	[490.186,84.3093],
	[495.181,83.9997],
	[499.984,84.5208],
	[504.978,85.0419],
	[510.357,85.5632],
	[515.160,85.8766],
	[519.770,86.2938],
	[524.572,87.0226],
	[529.950,87.8554],
	[535.136,88.9997],
	[540.322,90.0401],
	[544.548,90.8725],
	[549.925,92.2246],
	[554.727,92.9533],
	[559.914,93.6822],
	[565.292,94.5151],
	[570.286,95.2439],
	[574.704,95.9725],
	[579.699,96.3898],
	[585.078,96.9111],
	[590.457,96.7056],
	[594.876,96.9150],
	[599.678,97.3322],
	[604.865,97.6457],
	[610.244,97.8555],
	[615.240,97.9613],
	[620.234,98.1709],
	[625.229,98.2767],
	[630.225,98.1748],
	[635.220,98.3845],
	[640.407,98.3865],
	[644.633,98.6997],
	[650.013,98.7018],
	[655.007,99.1191],
	[660.002,99.4325],
	[665.382,99.1231],
	[669.801,99.1248],
	[675.180,99.3346],
	[680.175,99.3365],
	[685.170,99.5462],
	[690.165,99.7558],
	[695.353,99.5501],
	[700.156,99.6559],
])

# separate X and Y values
thylamys_x = np.empty(thylamys_lens.shape[0])
thylamys_y = np.empty(thylamys_lens.shape[0])
for i in range(thylamys_lens.shape[0]):
	thylamys_x[i] = round(thylamys_lens[i][0]) # again, dot plot
	thylamys_y[i] = thylamys_lens[i][1]
brushtail_x = np.empty(brushtail_eye.shape[0])
brushtail_y = np.empty(brushtail_eye.shape[0])
for i in range(brushtail_eye.shape[0]):
	brushtail_x[i] = brushtail_eye[i][0] # not a dot plot
	brushtail_y[i] = brushtail_eye[i][1]
# special case for mouse
mouse_x = np.empty(42)
mouse_x[0] = 300
mouse_x[1] = 309
for i in range(2, 42): mouse_x[i] = (i-1)*10 + 300

# interpolate
x_10nm = np.empty(41)
for i in range(41): x_10nm[i] = i*10 + 300
x_1nm = np.empty(401)
for i in range(401): x_1nm[i] = i + 300

# divide by 100
mouse_lens = mouse_lens / 100
brushtail_y = brushtail_y / 100

# 10-nm intervals
mouse_10nm = np.empty(41)
mouse_10nm[0] == 0
for i in range(1, 41): mouse_10nm[i] = mouse_lens[i+1] # remove the extra 0 at 309
thylamys_10nm = np.interp(x_10nm, thylamys_x, thylamys_y)
thylamys_10nm = thylamys_10nm / thylamys_10nm[40] # value at 700 nm
brushtail_10nm = np.interp(x_10nm, brushtail_x, brushtail_y)

# 1-nm intervals
mouse_1nm = np.interp(x_1nm, mouse_x, mouse_lens)
thylamys_1nm = np.interp(x_1nm, thylamys_x, thylamys_y)
thylamys_1nm = thylamys_1nm / thylamys_1nm[400] # value at 700 nm
brushtail_1nm = np.interp(x_1nm, brushtail_x, brushtail_y)

def ocular_media(w):
	if (args.media == "human"): return human_filter(w)
	elif (args.media == "mouse"):
		w = round(w)
		return mouse_1nm[w-300]
	elif (args.media == "thylamys" or args.media == "opossum"):
		w = round(w)
		return thylamys_1nm[w-300]
	elif (args.media == "brushtail"):
		w = round(w)
		return brushtail_1nm[w-300]
	return 1

# relative wavelength sensitivity
def sensitivity(w, l=l1, m=m1, s=s1):
	if (args.media == "human2"): return luminosity[w-300]
	return (args.ls*vpt(w, l) + args.ms*vpt(w, m) + args.ss*vpt(w, s) + args.rb*vpt(w, rod)) * ocular_media(w)

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
		[vpt(r, l1)*ocular_media(r), vpt(g, l1)*ocular_media(g), vpt(b, l1)*ocular_media(b)],
		[vpt(r, m1)*ocular_media(r), vpt(g, m1)*ocular_media(g), vpt(b, m1)*ocular_media(b)],
		[vpt(r, s1)*ocular_media(r), vpt(g, s1)*ocular_media(g), vpt(b, s1)*ocular_media(b)]
	])

	# sensitivity of cones to each wavelength from 300 to 800 nm in 1 nm increments
	matrix_c = np.empty((3, 401)) # 500x3 matrix -- this is height x width, not width x height

	# red row
	for i in range(401):
		matrix_c[0][i] = vpt(i + 300, l1)*ocular_media(i+300)
	# green row
		matrix_c[1][i] = vpt(i + 300, m1)*ocular_media(i+300)
	# blue row
		matrix_c[2][i] = vpt(i + 300, s1)*ocular_media(i+300)

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
	matrix_w = [
		[1/rw, 0, 0],
		[0, 1/gw, 0],
		[0, 0, 1/bw]
	]

	# final color match matrix
	matrix_cmf = np.matmul(matrix_w, matrix_m)

if (args.cmf):
	print("LMS->RGB matrix: " + str(np.matmul(matrix_w, np.linalg.inv(matrix_a))))
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

# determine hue, saturation and lightness for a spectral power distribution
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
			if (reflect == True): table_l += vpt(w, lw) * table[i] * light_source[i] * ocular_media(w)
			else: table_l += vpt(w, lw) * table[i] * ocular_media(w)
			# remove either M or S for dichromacy
			if (mw != lw):
				if (reflect == True): table_m += vpt(w, mw) * table[i] * light_source[i] * ocular_media(w)
				else: table_m += vpt(w, mw) * table[i] * ocular_media(w)
			if (sw != mw):
				if (reflect == True): table_s += vpt(w, sw) * table[i] * light_source[i] * ocular_media(w)
				else: table_s += vpt(w, sw) * table[i] * ocular_media(w)
			
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
		wpl += vpt(w, lw) * light_source[i] * ocular_media(w)
		wpm += vpt(w, mw) * light_source[i] * ocular_media(w)
		wps += vpt(w, sw) * light_source[i] * ocular_media(w)
	
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

	# declare posx and posy
	posx = -1
	posy = -1
	
	if (output):
		# render color
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
			table_l_gamma = gamma(table_l)
			table_m_gamma = gamma(table_m)
			table_s_gamma = gamma(table_s)
			print("RGB=LMS false color: (" + str(round(table_l_gamma * 255)) + ", " + str(round(table_m_gamma * 255)) + ", " + str(round(table_s_gamma * 255)) + ")")
			try:
				hsl = convert_color(sRGBColor(table_l, table_m, table_s), HSLColor)
				print("LMS->HSL: " + str(hsl.get_value_tuple()))
			except ZeroDivisionError: print("Warning: HSL conversion failed (ZeroDivisionError)")
			
			# CMF false color rendering. Out-of-gamut colors are scaled down if
			# too large and clipped toward the white point if negative.
			rgb_max = max(table_r,table_g,table_b)
			rgb_min = min(table_r,table_g,table_b,0)
			rgb_max2 = max(table_r-rgb_min,table_g-rgb_min,table_b-rgb_min,rgb_max,1)
			print("CMF false color: (" + str(round(table_r * 255)) + ", " + str(round(table_g * 255)) + ", " + str(round(table_b * 255)) + ")")
			print("CMF false color (within gamut): (" + str((table_r-rgb_min) * 255 / rgb_max2) + ", " + str((table_g-rgb_min) * 255 / rgb_max2) + ", " + str((table_b-rgb_min) * 255 / rgb_max2) + ")")
			try:
				hsl = convert_color(sRGBColor(table_r-rgb_min, table_g-rgb_min, table_b-rgb_min), HSLColor)
				print("CMF->HSL: " + str(hsl.get_value_tuple()))
			except ZeroDivisionError: print("Warning: HSL conversion failed (ZeroDivisionError)")
				
			
			# This should be the nearest wavelength on a line from the white point.
			match = 300
			x0 = math.sqrt(1/2)*(table_l - table_s) / (table_l + table_m + table_s)
			y0 = math.sqrt(2/3)*(table_m - 0.5*(table_l + table_s)) / (table_l + table_m + table_s)
			wx = math.sqrt(1/2)*(wpl - wps) / (wpl + wpm + wps)
			wy = math.sqrt(2/3)*(wpm - 0.5*(wpl + wps)) / (wpl + wpm + wps)
			
			# adjust to chosen visible range
			while (ocular_media(match) == 0):
				match += 1
			
			for i in range(300, 701):
				matchl1 = vpt(match, lw)*ocular_media(match) / wpl
				matchm1 = vpt(match, mw)*ocular_media(match) / wpm
				matchs1 = vpt(match, sw)*ocular_media(match) / wps
				matchl2 = vpt(i, lw)*ocular_media(i) / wpl
				matchm2 = vpt(i, mw)*ocular_media(i) / wpm
				matchs2 = vpt(i, sw)*ocular_media(i) / wps
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
		
		# "SpectralColor" conversion for comparison
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
		# break up text
		print("")
	
	# return various items:
	# 0: brightness
	# 1: color triangle x
	# 2: color triangle y
	# 3: L-M position
	return([brightness, posx, posy, lm])

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
			ql1 += vpt(w, lw) * table1[i] * wp_1nm[i] * ocular_media(w)
			qm1 += vpt(w, mw) * table1[i] * wp_1nm[i] * ocular_media(w)
			qs1 += vpt(w, sw) * table1[i] * wp_1nm[i] * ocular_media(w)
		wpl += vpt(w, lw) * wp_1nm[i] * ocular_media(w)
		wpm += vpt(w, mw) * wp_1nm[i] * ocular_media(w)
		wps += vpt(w, sw) * wp_1nm[i] * ocular_media(w)
	for i in range(0, table2.shape[0]):
		w = i + 300
		if (table2[i] > 0):
			ql2 += vpt(w, lw) * table2[i] * wp_1nm[i] * ocular_media(w)
			qm2 += vpt(w, mw) * table2[i] * wp_1nm[i] * ocular_media(w)
			qs2 += vpt(w, sw) * table2[i] * wp_1nm[i] * ocular_media(w)
	
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
		T = 1 / ocular_media(700) # set to 1 at 700 nm
		
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
				aql1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, lw)*l)) * vpt(w, lw) * table1[i] * photons[i] * ocular_media(w)
				aqm1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, mw)*l)) * vpt(w, mw) * table1[i] * photons[i] * ocular_media(w) 
				aqs1 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, sw)*l)) * vpt(w, sw) * table1[i] * photons[i] * ocular_media(w)
		for i in range(table2.shape[0]):
			w = i + 300
			if (table2[i] > 0):
				aql2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, lw)*l)) * vpt(w, lw) * table2[i] * photons[i] * ocular_media(w)
				aqm2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, mw)*l)) * vpt(w, mw) * table2[i] * photons[i] * ocular_media(w)
				aqs2 += (math.pi/4)**2*(d / f)**2*D**2*K*T*dt*(1 - math.exp(-k*vpt(w, sw)*l)) * vpt(w, sw) * table2[i] * photons[i] * ocular_media(w)
		
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
# I've changed it to a signed value because we care about the direction.
def brightness_contrast(table1, table2, lw=args.lw):
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
	
	# We're redoing this to use sensitivity() so we can use the photopic luminosity
	# function for humans.
	q1 = 0
	q2 = 0
	for i in range(table1.shape[0]):
		w = i + 300
		wpt += sensitivity(w, lw) * wp_1nm[i]
		q1 += sensitivity(w, lw) * table1[i] * wp_1nm[i]
	for i in range(table2.shape[0]):
		w = i + 300
		q2 += sensitivity(w, lw) * table2[i] * wp_1nm[i]
	
	df = math.log(q1 / q2)
	delta_s = df / args.weberb # not an absolute value
	return delta_s

# print specified information
print("L: " + str(l1) + ", M: " + str(m1) + ", S: " + str(s1) + ", rod: " + str(rod))
print("")

"""
plot visual pigment templates
I changed the line colors back to red, green and blue. I've been trying to avoid this kind
of color coding, but in this case it's generally understood which visual pigment template
is which based on its position on the graph, and the shades-of-gray "color" scheme was too
pale for my liking.
"""
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
		if (args.media != "none"):
			lvalues[i] = vpt(i+300, l1) * ocular_media(i+300)
			mvalues[i] = vpt(i+300, m1) * ocular_media(i+300)
			svalues[i] = vpt(i+300, s1) * ocular_media(i+300)
			rvalues[i] = vpt(i+300, rod) * ocular_media(i+300)
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
		yvalues1[i-300] = ocular_media(i)
	print("Maximum sensitivity: " + str(ms))
	print("")
	
	plt.plot(xvalues, yvalues/sensitivity(ms), 'k')
	if (args.media != "none"):
		plt.plot(xvalues, yvalues1, ':k')
	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.show()

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
		wl += vpt(w, l1) * wp[i] * ocular_media(w)
		wm += vpt(w, m1) * wp[i] * ocular_media(w)
		ws += vpt(w, s1) * wp[i] * ocular_media(w)
		el += vpt(w, l1) * e[i] * ocular_media(w)
		em += vpt(w, m1) * e[i] * ocular_media(w)
		es += vpt(w, s1) * e[i] * ocular_media(w)
		d65l += vpt(w, l1) * d65[i] * ocular_media(w)
		d65m += vpt(w, m1) * d65[i] * ocular_media(w)
		d65s += vpt(w, s1) * d65[i] * ocular_media(w)
	
	for i in range(0, 41):
		w = i*10 + 300
		labels[i] = w
		l = vpt(w, l1) * ocular_media(w)
		m = vpt(w, m1) * ocular_media(w)
		s = vpt(w, s1) * ocular_media(w)
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
def color_disc(f0, f1, f2, f3, s0, s1, s2, s3, correct=35, trials=40, alternative='greater', lw = args.lw, mw = args.mw, sw = args.sw, output=True):
	d = 0 # different
	s = 0 # same
	counter = 0
	box = np.empty(16)
	
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
			box[counter] = contrast
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
	mc = statistics.median(box)
	minc = min(*box)
	maxc = max(*box)
	if (output):
		print("Median contrast: " + str(mc))
		print("")
		#return
	return [mc, minc, maxc, box]


"""
like color_disc except for brightness
We want to know both the direction and the strength of the difference
to judge whether it could be reliably used as a cue.
"""
def brightness_disc(f0, f1, f2, f3, s0, s1, s2, s3, correct=35, trials=40, lw=args.lw, output=True): 
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
			contrast = brightness_contrast(flevel, slevel, lw=lw)
			if (contrast >= 1): fb += 1
			elif (contrast <= -1): sb += 1
			else: equal += 1
	
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
