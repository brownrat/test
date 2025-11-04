"""
This file contains the following:
* list of arguments
* visual pigments: standard visual pigment template (vpt()) and human cone
fundamentals
* ocular media transmittance for various species
* blackbody: blackbody()
* color/brightness analysis: color_contrast(), brightness_contrast() and
spec2rgb()
* illuminants: D65, E, blackbody (identical to CIE A by default but can be
customized)
* color matching functions
* Maxwell triangle color space plots: triangle()
* color and brightness discrimination: color_disc() and brightness_disc()

Other stuff has been spun off to the following files, which import this one.
I'd like to add additional arguments to those but don't see how to do that.

* erg.py: analysis of ERG data from Jacobs & Williams (2010)
* illuminants.py: test illuminants
* optimal.py: optimal colors
* peak.py: Gaussian peak spectra
* ramp.py: ramp function spectra
* sky.py: sky/daylight spectra
* wratten.py: Kodak Wratten filters
* munsell.py: Munsell color cards
* image_processing.py: false-color images

wratten.py and munsell.py have different parts accessed with different
arguments (see list), so they do nothing by default.
"""

# this list doesn't recursively import
import math
import numpy as np
import time
import argparse
import colormath
from colormath.color_objects import sRGBColor, SpectralColor, LabColor
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
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--receptors", "--target", nargs="+", type=int, default=[562], help="receptor peak sensitivities")
parser.add_argument("--rcsv", nargs="+", default="none", help="custom receptors from CSV file(s)")
parser.add_argument("--weber", type=float, default=0.05, help="base Weber fraction")
parser.add_argument("--weberr2", type=float, default=0, help="override Weber fraction for receptor 2")
parser.add_argument("--weberr3", type=float, default=0, help="override Weber fraction for receptor 3")
parser.add_argument("--weberr4", type=float, default=0, help="override Weber fraction for receptor 4")
parser.add_argument("--weberb", type=float, default=0.11, help="achromatic Weber fraction")
parser.add_argument("--rf", nargs="+", type=float, default=[1,1,1,1,1], help="number of receptors of each type per receptive field")
parser.add_argument("--media", "--filter", type=str, default="none", help="ocular media transmittance")
parser.add_argument("--qn", help="use quantum noise in color differences", action="store_true")
parser.add_argument("--r1s", type=float, default=1, help="contribution of receptor 1 to spectral sensitivity")
parser.add_argument("--rs", type=float, default=0, help="contribution of rods to spectral sensitivity")
parser.add_argument("--vonkries", help="enable von Kries transform", action="store_true")
parser.add_argument("--sensitivity", help="show spectral sensitivity function", action="store_true")
parser.add_argument("--white", default="d65", help="set illuminant")
parser.add_argument("--wratten", help="Kodak Wratten camera filters", action="store_true")
parser.add_argument("--wv2", help="use tabulated Kodak Wratten data from Kodak Photographic Filters Handbook (1990)", action="store_true")
parser.add_argument("--ndswap", help="use flat transmittance instead of neutral density curves", action="store_true")
parser.add_argument("--grayswap", help="set gray to yellow/green average", action="store_true")
parser.add_argument("--wcheck", help="find Wratten filter brightness matches", action="store_true")
parser.add_argument("--wcheck2", help="plot brightness matches and show radiance/luminance", action="store_true")
parser.add_argument("--wopt1", help="test varying spectral sensitivity (probability)", action="store_true")
parser.add_argument("--wopt2", help="test varying spectral sensitivity (order)", action="store_true")
parser.add_argument("--wopt3", help="test varying spectral sensitivity (contrast)", action="store_true")
parser.add_argument("--wopt4", help="test varying spectral sensitivity (contrast)", action="store_true")
parser.add_argument("--vpt", help="plot visual pigment templates", action="store_true")
parser.add_argument("--plot", nargs="+", help="plot spectra from CSV files")
parser.add_argument("--ylabel", default="Reflectance", help="Y-axis label")
parser.add_argument("--colorplot", help="use generated colors for plots", action="store_true")
parser.add_argument("--csplot", "--triangle", nargs="*", help="create a Maxwell-style color space plot")
parser.add_argument("--labels", nargs="+", default=['','',''], help="labels to use in triangle and line plots")
parser.add_argument("--fcmode", help="switch default false color mode", default="none")
parser.add_argument("--s2tmode", help="switch default color transform mode", default="c")
parser.add_argument("--munsell", help="Munsell color cards", action="store_true")
parser.add_argument("--mcheck", help="plots and brightness", action="store_true")
parser.add_argument("--mv2", help="use alternative Munsell spectra", action="store_true")
parser.add_argument("--mv3", help="use alternative Munsell spectra for <420nm only", action="store_true")
parser.add_argument("--mv4", help="remove wavelengths <400nm from Munsell spectra", action="store_true")
parser.add_argument("--mopt1", help="find optimal visual pigments for Munsell cards (dichromacy)", action="store_true")
parser.add_argument("--mopt2", help="find optimal visual pigments for Munsell cards (trichromacy)", action="store_true")
parser.add_argument("--mopt3", help="find L-M chromatic spread/overlap of Munsell cards", action="store_true")
parser.add_argument("--summation", type=float, default=1, help="spatial summation")
parser.add_argument("--dt", type=float, default=1/22.4, help="integration time (s)")
parser.add_argument("-d", "--osd", type=float, default=1, help="outer segment diameter (um)")
parser.add_argument("-f", "--fl", type=float, default=6.81, help="focal length (mm)")
parser.add_argument("--pupil", type=float, default=6, help="pupil diameter (mm)")
parser.add_argument("-k", "--od", type=float, default=0.009, help="optical density of pigment (um^-1)")
parser.add_argument("-l", "--osl", type=float, default=30, help="outer segment length (um)")
parser.add_argument("--dist", type=float, default=5, help="distance from light source (cm)")
parser.add_argument("--sdist", type=float, default=0, help="distance from stimulus (cm)")
parser.add_argument("--width", type=float, default=0, help="width of stimulus (cm)")
parser.add_argument("--height", type=float, default=0, help="height of stimulus (cm)")
parser.add_argument("--square", help="make stimulus a square/rectangle instead of circle/ellipse", action="store_true")
parser.add_argument("--watts", type=float, default=12, help="illuminant power in watts")
parser.add_argument("--w2p", help="compare watt and photon illuminant spectra", action="store_true")
parser.add_argument("--lux", "--lx", type=float, default=-1, help="illuminance in lux")
parser.add_argument("--edi", type=float, default=-1, help="illuminance in EDI (equivalent daylight illuminance)")
parser.add_argument("--ct", type=int, default=2856, help="color temperature of blackbody illuminant")
parser.add_argument("--selfscreen", help="add pigment self-screening", action="store_true")
parser.add_argument("--subjects", type=int, default=1, help="number of subjects")
parser.add_argument("--warnings", help="show warnings", action="store_true")
parser.add_argument("--cmf", help="color-matching functions", action="store_true")
parser.add_argument("--primaries", nargs="+", type=int, default=[-1], help="primary wavelengths for color-matching functions")
parser.add_argument("-p", "--preset", help="use specified preset vision system", default="none")
parser.add_argument("-v", "--verbose", help="", action="store_true")
parser.add_argument("-i", "--image", default="none", help="path to input image")
parser.add_argument("-u", "--uvimage", default="none", help="path to input UV image")
parser.add_argument("--outimage", default="result.png", help="path to save output image")
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
dimension = len(args.receptors)

# Weber fractions
wr1 = args.weber
# set number per receptive field to 1 when not specified
rf = [1,1,1,1,1]
for i in range(len(args.rf)): rf[i] = args.rf[i]
wr2 = wr1 * math.sqrt(rf[0]) / math.sqrt(rf[1])
wr3 = wr1 * math.sqrt(rf[0]) / math.sqrt(rf[2])
wr4 = wr1 * math.sqrt(rf[0]) / math.sqrt(rf[3])
if (args.weberr2 > 0): wr2 = args.weberr2
if (args.weberr3 > 0): wr3 = args.weberr3
if (args.weberr4 > 0): wr4 = args.weberr4

"""
read CSV files

This function allows for using arbitrary CSV files as visual pigments, ocular
media, illuminants, or reflectance spectra. By default, only the first two
columns are read and the first column is interpreted as wavelength, and it
returns the second column interpolated to 1-nm increments from 300-700 nm.
Additional columns can be read by setting numcols to a value higher than 2.
It would be nice to automatically detect how many columns there are, but that
works best if the first row is a header and the rest are values, and the files
I'm working with aren't like that. Interpolation can be disabled by setting
interp=False, in which case all the columns will be returned as-is.

If the range of data doesn't cover 300-700, values not included are filled with
0. If the data is the log of what you really want (e.g. absorption vs.
transmission), you should set log=True so negative infinity is used instead
(10**-inf = 0, but 10**0 = 1).
"""
def csv2spec(filename, numcols=2, log=False, interp=True):
	if (args.verbose): print("Reading CSV file: " + filename)
	# use negative infinity as "zero" for logarithmic data
	zero = 0
	if (log): zero = -np.inf
	with open(filename, newline='') as csvfile:
		fieldnames = []
		cols = []
		for i in range(numcols):
			name = 'col' + str(i)
			fieldnames.append(name)
			cols.append([])
		reader = csv.DictReader(csvfile, fieldnames=fieldnames)
		# This bit uses an exception (can be either ValueError or
		# TypeError) to ignore headers and other extra text. This kind
		# of thing is considered bad coding practice but seems to be the
		# only way to get the behavior I want when reading PlotDigitizer
		# files.
		for row in reader:
			try:
				if (len(row) > 1 and row[fieldnames[numcols-1]] != None):
					for i in range(min(len(row), numcols)):
						if (row[fieldnames[i]] == ''):
							cols[i].append(zero)
						else:
							f = float(row[fieldnames[i]])
							cols[i].append(f)
			except:
				if (args.warnings):
					print("The following row of " + filename + " was ignored: "
					+ str(row))

	# interpolate to 300-700
	if (interp):
		wl = cols[0]
		values = []
		for i in range(1, numcols):
			interp = np.interp(x_1nm, wl, cols[i])
			values.append(interp)

		# remove extrapolation
		wl_min = min(*wl)
		wl_max = max(*wl)
		
		if (wl_min > 300):
			for i in range(math.floor(wl_min) - 300):
				for j in range(numcols-1): values[j][i] = zero
		if (wl_max < 700):
			for i in range(math.ceil(wl_max) - 300, 401):
				for j in range(numcols-1): values[j][i] = zero

		if (numcols == 2): return values[0]
		return values

	return cols

# types of ocular media
# The presets don't include transmission/absorption below a certain wavelength
# (310 nm for mouse, 390 nm for human), so it's set to 0 by csv2spec().

# mouse (Jacobs & Williams 2007)
mouse_1nm = csv2spec('mouse.csv')

# divide by 100
mouse_1nm /= 100

# 10-nm intervals
mouse_10nm = mouse_1nm[::10]

# human macular pigment and ocular media
# http://www.cvrl.org/maclens.htm (CVRL functions, 1 nm)
mac = csv2spec('macss_1.csv', log=True)
lens = csv2spec('lensss_1.csv', log=True)
for i in range(401):
	mac[i] = 10**-abs(mac[i])
	lens[i] = 10**-abs(lens[i])

# choose media type
media_1nm = np.ones(401)
media_10nm = np.ones(41)
if (args.media == "mouse"):
	media_1nm = mouse_1nm
	media_10nm = mouse_10nm
elif (args.media == "cvrl"): media_1nm = lens
elif (args.media == "cvrlmac"): media_1nm = mac*lens
elif (args.media != "none"):
	media_1nm = csv2spec(args.media)
	media_1nm /= media_1nm[400]
	media_10nm = media_1nm[::10]

# quantal human cone fundamentals -- http://www.cvrl.org/cones.htm
cie2 = csv2spec('ss2_10q_1.csv', numcols=4, log=True)
cie2_l = np.zeros(401)
cie2_m = np.zeros(401)
cie2_s = np.zeros(401)
for i in range(401):
	cie2_l[i] = 10**cie2[0][i]
	cie2_m[i] = 10**cie2[1][i]
	cie2_s[i] = 10**cie2[2][i]

cie10 = csv2spec('ss10q_1.csv', numcols=4, log=True)
cie10_l = np.zeros(401)
cie10_m = np.zeros(401)
cie10_s = np.zeros(401)
for i in range(401):
	cie10_l[i] = 10**cie10[0][i]
	cie10_m[i] = 10**cie10[1][i]
	cie10_s[i] = 10**cie10[2][i]

"""
find sensitivity of a given photoreceptor type to a wavelength using visual
pigment templates from Govardovskii et al. (2000). l is short for lambda.
See also https://pmc.ncbi.nlm.nih.gov/articles/PMC2962788/

Palacios et al. 2010 use a different equation for the beta-band peak
(123 + 0.429 * alpha peak).
"""
def vpt(w, lmax, od=-1):
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
	"""
	self-screening

	For what self-screening is, see: https://pmc.ncbi.nlm.nih.gov/articles/PMC3269789/

	This is the same as equation 4.12 from The Optics of Life (p. 110).
	A similar equation is given in Lamb (1995): log10[1 - S(1 - 10^-D)]/(-D),
	where S is the pigment absorption and D is the optical density. Here he's
	dividing the log of sensitivity by the density, so the division isn't
	really part of the sensitivity. I previously had a mistake where the
	whole thing was multiplied by the non-screened value, which meant the
	resulting curve was narrower rather than broader. Look again at the
	official equations and you'll see we don't need to do this.
	"""
	value = alpha + beta
	if (od == -1):
		if (args.selfscreen): return (1 - math.exp(-args.od*value*args.osl))
		return value
	else: return (1 - math.exp(-od*value))

# presets

# CIE 2-deg fundamentals
cie2 = [cie2_l, cie2_m, cie2_s]

# CIE 10-deg fundamentals
cie10 = [cie10_l, cie10_m, cie10_s]

# typical UVS and VS birds, digitized with PlotDigitizer
# source: https://www.nature.com/articles/s41467-018-08142-5
# licensed under Creative Commons 4.0 http://creativecommons.org/licenses/by/4.0/
bird_l = csv2spec('bird-l.csv')
bird_m = csv2spec('bird-m.csv')
bird_su = csv2spec('bird-su.csv')
bird_sv = csv2spec('bird-sv.csv')
bird_u = csv2spec('bird-u.csv')
bird_v = csv2spec('bird-v.csv')

# set max values to 1
bird_l /= max(*bird_l)
bird_m /= max(*bird_m)
bird_su /= max(*bird_su)
bird_sv /= max(*bird_sv)
bird_u /= max(*bird_u)
bird_v /= max(*bird_v)

uvbird = [bird_l,bird_m,bird_su,bird_u]
vbird = [bird_l,bird_m,bird_sv,bird_v]

# set receptors and luminosity
rlist = []
r1_1nm = np.zeros(401)
r2_1nm = np.zeros(401)
r3_1nm = np.zeros(401)
r4_1nm = np.zeros(401)
luminosity = np.ones(401)

# CIE 2-deg fundamentals
if (args.preset == "cie2"):
	rlist = cie2
# CIE 10-deg fundamentals
elif (args.preset == "cie10"):
	rlist = cie10
# UVS bird
elif (args.preset == "uvbird"): rlist = uvbird
# VS bird
elif (args.preset == "vbird"): rlist = vbird
elif (args.rcsv != "none"):
	dimension = len(args.rcsv)
	for i in range(len(args.rcsv)):
		rlist.append(csv2spec(args.rcsv[i]))
else:
	for i in range(401):
		w = i+300
		luminosity[i] = (args.r1s*vpt(w, receptors[0])) * media_1nm[w-300]

	for i in range(dimension):
		# whole list
		yvalues = np.empty(401)
		peak = args.receptors[i] # peak sensitivity
		for j in range(401):
			yvalues[j] = vpt(j+300, peak)
		rlist.append(yvalues)

dimension = len(rlist)
if (dimension > 3): r4_1nm = rlist[3]
if (dimension > 2): r3_1nm = rlist[2]
if (dimension > 1): r2_1nm = rlist[1]
r1_1nm = rlist[0]

# black body -- SI units, returns energy in joules/m^2/nm/sec/steradian = watts/m^2/nm/steradian
# The "1e-9" converts from meters to nanometers, as this is not cubic meters but
# "per wavelength".
# https://yceo.yale.edu/sites/default/files/files/ComputingThePlanckFunction.pdf
h = 6.626068e-34 # Planck's constant (joule sec)
k = 1.38066e-23 # Boltzmann's constant (joule/deg)
c = 2.997925e+8 # speed of light (m/s)

def blackbody(nm, t):
	wl = nm / 1000000000 # wavelength (m)
	try:
		value = 1e-9*2*h*c**2*wl**-5 / (math.exp((h*c) / (k*wl*t)) - 1)
	except OverflowError:
		if (args.warnings): print("Warning (blackbody): math overflow, returning 1")
		return 1
	return value

# watts to photons/second (The Optics of Life, p. 12)
watts2photons = np.empty(401)
for i in range(401): watts2photons[i] = 1e-9*(i+300) / (h*c)

# human photopic luminosity function
# https://cie.co.at/datatable/cie-spectral-luminous-efficiency-photopic-vision
cie_photopic_full = csv2spec('CIE_sle_photopic.csv', interp=False)
# add zeroes
for i in range(60):
	cie_photopic_full[0].insert(0, 359-i)
	cie_photopic_full[1].insert(0, 0)
cie_photopic = np.array(cie_photopic_full[1][:401])

# The luminous efficiency function relates to watts rather than photons, so
# we convert it for compatibility.
if (args.preset == ('cie2' or 'cie10')):
	luminosity = cie_photopic / watts2photons

luminosity /= max(*luminosity)

# illuminants

# E for easy
e_1nm = np.ones(401)

# D65 300-830
# http://cie.co.at/datatable/cie-standard-illuminant-d65
d65_full = csv2spec("CIE_std_illum_D65.csv", interp=False)
d65_1nm = np.array(d65_full[1][:401])

# CIE A
a_1nm = csv2spec("CIE_std_illum_A_1nm.csv")

# black body with specified color temperature
blackbody_1nm = np.empty(401)
for i in range(401):
	blackbody_1nm[i] = blackbody(i+300, args.ct)

# white point
if (args.white == "d65"):
	wp_1nm = d65_1nm
elif (args.white == "e"):
	wp_1nm = e_1nm
elif (args.white == "a"):
	wp_1nm = a_1nm
elif (args.white == "i" or args.white == "blackbody"):
	wp_1nm = blackbody_1nm
# custom illuminant
else:
	wp_1nm = csv2spec(args.white)

# absolute numbers of photons

"""
Here we calculate the absolute number of photons per wavelength produced
by the blackbody illuminant. Note the difference between cubic meters
and square meter-wavelength (see notes on quantum noise). The premise is
that an incandescent bulb is a "gray body" with an arbitrary color
temperature and an emissivity determined by the energy going into the
filament. That's not really how they work, but color temperature doesn't
exactly depend on wattage and this is the most convenient way to model
that.
"""
sphere_area = 4*math.pi*args.dist**2 # sphere surface area
watts_cm2 = args.watts / (sphere_area)
watts_m2 = watts_cm2 * 100**2

if (args.white == ('i' or 'blackbody')):
	# scale it down to what amount of energy would be produced by the specified
	# number of watts

	# find total energy produced using Stefan-Boltzmann constant
	# We divide by pi steradians to get the radiance and make the units equivalent to
	# blackbody() (not 4pi, this gives the wrong value -- see p. 61 here:
	# https://eodg.atm.ox.ac.uk/user/grainger/research/book/protected/Chapter3.pdf)
	# Except not, because we're actually looking at irradiance. Irradiance goes
	# down with the square of the distance, but radiance doesn't. We then get the
	# radiance of the thing lit by that irradiance. Confusing, right? See also:
	# The Optics of Life, p. 20-24
	sigma = 5.670374419e-8
	energy = sigma * args.ct**4

# We can also scale other illuminants to specified watts, but here we can't
# use information outside the range we're given. This is why I don't use the
# "real" CIE A.
else:
	energy = sum(wp_1nm)

scale = watts_m2 / energy
if (args.verbose): print("W/m^2 scaling factor: " + str(scale))
wp_watts = scale * wp_1nm

# integral of D65 x human photopic luminosity function
# The units we want at the end are W/m^2, so we multiply this by the
# lux-to-watts scaling factor for 555 nm (683.002 lm/W). This will end
# up on the bottom.
integral = sum(np.array(d65_full[1]) * np.array(cie_photopic_full[1])) * 683.002

# By default, the illuminant is scaled in watts. (We need a default.)
# This can be changed to lux by setting --lux or EDI (equivalent daylight
# illuminance) by setting --edi. For an explanation of EDI, see:
# https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-02038-1
if (args.lux != -1 or args.edi != -1):
	# lux scaling factor
	# I previously wrote that this has units of nanometers, but the integral we
	# just found doesn't include them because it's implicitly multiplied by the
	# wavelength bin. The scaling factor is unitless.
	lxscale = args.lux / integral
	if (args.verbose): print("Lux scaling factor: " + str(lxscale))

	# scaled version
	# units: (W/m^2 * m) / (joule*sec * m/s * sr) = (s * m^2 * sr)^-1
	d65_watts = lxscale * d65_1nm / math.pi

	# We can also scale other illuminants this way.
	if (args.white == 'd65'): wp_watts = d65_watts
	else:
		if (args.lux != -1):
			lxscale = args.lux / (sum(wp_1nm * cie_photopic) * 683.002)
			wp_watts = lxscale * wp_1nm / math.pi
		else:
			d65_lx = args.edi / integral
			d65_opic = d65_lx * sum(d65_1nm*watts2photons*luminosity)
			wp_opic = sum(wp_1nm*watts2photons*luminosity)
			ediscale = d65_opic / wp_opic
			wp_watts = ediscale * wp_1nm / math.pi

# We just found watts above because the original spectrum is assumed to be watts.
# Watts->photons can be done separately.
# Even though we don't need the absolute number for relative quantum catches, we
# still need to have photons and not watts, so we replace the original illuminant
# array with this.
wp_1nm = wp_watts * watts2photons
wp_10nm = wp_1nm[::10]

# show comparison of illuminant spectrum in watts and photons (for testing)
if (args.w2p):
	# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/two_scales.html
	fig, ax1 = plt.subplots()
	color='tab:red'
	ax1.set_xlabel('wavelength (nm)')
	ax1.set_ylabel('watts', color=color)
	ax1.plot(x_1nm, wp_watts, color=color)
	ax1.tick_params(axis='y', labelcolor=color)

	ax2 = ax1.twinx()

	color = 'tab:blue'
	ax2.set_ylabel('photons', color=color)
	ax2.plot(x_1nm, wp_1nm, ':', color=color)
	ax2.tick_params(axis='y', labelcolor=color)

	fig.tight_layout()
	plt.show()

"""
CMFs (color-matching functions)

For how this works, see: https://horizon-lab.org/colorvis/cone2cmf.html

In humans, CMFs are related to perception of "red", "green" and "blue"
unique hues that derive from opponent processing, as mentioned here:
https://pmc.ncbi.nlm.nih.gov/articles/PMC8084791/ The original human
CMFs can be found at http://www.cvrl.org/.

By default, the primary wavelengths will be automatically chosen as
whichever are the closest to stimulating only one receptor type. This
means they won't match what's best for a digital display because the
longest and shortest ones tend to be placed way at the ends of the
spectrum (limited to the range 300-700). All I care about here is the
shape of the curves. They also don't directly relate to unique hues
because the "corners" of human color space are closer to teal and violet
than blue and green.

Some recommended values for trichromatic primaries:
* marsupials: 620, 450, 360 (Arrese et al. 2006 doi.org/10.1016/j.cub.2006.02.036)
* humans/primates:
** 700, 546.1, 435.8 (CIE 1931 https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_RGB_color_space)
** 630, 532, 467 (Rec. 2020 https://en.wikipedia.org/wiki/Rec._2020)
** 612, 549, 465 (dominant wavelengths of sRGB primaries: see:
*** https://en.wikipedia.org/wiki/File:SRGB_chromaticity_CIE1931.svg
*** https://clarkvision.com/articles/color-spaces/
** other color spaces: see https://en.wikipedia.org/wiki/RGB_color_spaces

Dichromatic, trichromatic and tetrachromatic CMFs are supported.
"""
if (args.cmf or (args.fcmode == 'cmf')):
	# primaries
	p1 = 0
	p2 = 0
	p3 = 0
	p4 = 0

	# find "optimal" primaries
	if (args.primaries[0] == -1):
		r1_cur = 0
		r2_cur = 0
		r3_cur = 0
		r1_max = 0
		r2_max = 0
		r3_max = 0
		for i in range(401):
			total = r1_1nm[i] + r2_1nm[i] + r3_1nm[i]
			if (total > 0):
				r1_cur = r1_1nm[i]/total
				r2_cur = r2_1nm[i]/total
				r3_cur = r3_1nm[i]/total
			if (r1_cur > r1_max):
				r1_max = r1_cur
				p1 = i+300
			if (r2_cur > r2_max):
				r2_max = r2_cur
				p2 = i+300
			if (r3_cur > r3_max):
				r3_max = r3_cur
				p3 = i+300
		if (args.verbose):
			print("Optimized primary wavelengths: "
				+ str(p1) + ", " + str(p2)
				+ ", " + str(p3))
	else:
		p1 = round(args.primaries[0])
		p2 = round(args.primaries[1])
		if (len(args.primaries) > 2): p3 = round(args.primaries[2])
		if (len(args.primaries) > 3): p4 = round(args.primaries[3])

def cmf():
	if (dimension > 3):
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
		matrix_cmf = cmf()
		# plot curves
		plt.plot(x_1nm, matrix_cmf[0], color='k')
		if (dimension > 1): plt.plot(x_1nm, matrix_cmf[1], color='k')
		if (dimension > 2): plt.plot(x_1nm, matrix_cmf[2], color='k')
		if (dimension > 3): plt.plot(x_1nm, matrix_cmf[3], color='k')
		if (dimension > 4): plt.plot(x_1nm, matrix_cmf[4], color='k')
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
generate an RGB color from a spectrum

This can produce either a "true color" using colormath or a "false color"
using the specified visual system. Out-of-gamut colors are scaled down if a
value is more than 1 and clipped if a value is less than 0.
"""
def spec2rgb(spec, reflect=True, mode=args.fcmode, scale=1, average=False):
	rgb_r = 0
	rgb_g = 0
	rgb_b = 0
	
	# interpolate
	if(len(spec) < 401): spec = np.interp(x_1nm, x_10nm, spec)
	
	# combine SPD with illuminant
	if (reflect and mode != 'none'): spec = spec * wp_1nm
	
	# remove nans
	for i in range(401):
		if (np.isnan(spec[i])): spec[i] = 0
	
	if (mode == "lms"):
		# LMS
		for i in range(spec.shape[0]):
			rgb_r += r1_1nm[i] * spec[i] * media_1nm[i]
			rgb_g += r2_1nm[i] * spec[i] * media_1nm[i]
			rgb_b += r3_1nm[i] * spec[i] * media_1nm[i]
	
		# von Kries transform
		if (args.vonkries):
			wpr1 = 0
			wpr2 = 0
			wpr3 = 0
			for i in range(401):
				wpr1 += r1_1nm[i] * wp_1nm[i] * media_1nm[i]
				wpr2 += r2_1nm[i] * wp_1nm[i] * media_1nm[i]
				wpr3 += r3_1nm[i] * wp_1nm[i] * media_1nm[i]
		
			rgb_r /= wpr1
			if (wpr2 > 0): rgb_g /= wpr2
			if (wpr3 > 0): rgb_b /= wpr3

		# LMS->RGB
		if (dimension < 3):
			rgb_b = rgb_g
			rgb_g = rgb_r
		if (dimension < 2):
			rgb_b = rgb_r

		# gamma correction
		rgb_r = gamma_encode(rgb_r)
		rgb_g = gamma_encode(rgb_g)
		rgb_b = gamma_encode(rgb_b)

	elif (mode == "cmf"):
		# CMF
		matrix_cmf = cmf()

		rgb_r = sum(matrix_cmf[0] * spec)
		if (dimension > 2):
			rgb_g = sum(matrix_cmf[1] * spec)
			rgb_b = sum(matrix_cmf[2] * spec)
		elif (dimension == 2):
			rgb_g = rgb_r
			rgb_b = sum(matrix_cmf[1] * spec)
				
		# gamma correction
		rgb_r = gamma_encode(rgb_r)
		rgb_g = gamma_encode(rgb_g)
		rgb_b = gamma_encode(rgb_b)
	
	else:
		# colormath true color
		# SpectralColor objects are limited to 340-830 nm in 10-nm
		# intervals. We have to choose whether to use every tenth
		# value or the average of the 10 values surrounding it.
		downsample = np.zeros(37)
		if (average):
			for i in range(37):
				v = 40+i*10
				mean = np.mean(spec[v-5:v+5])
				if (mean > 0): downsample[i] = mean
		else: downsample = spec[40::10]
		spectral = SpectralColor(spec_340nm=downsample[0],
		spec_350nm=downsample[1],
		spec_360nm=downsample[2],
		spec_370nm=downsample[3],
		spec_380nm=downsample[4],
		spec_390nm=downsample[5],
		spec_400nm=downsample[6],
		spec_410nm=downsample[7],
		spec_420nm=downsample[8],
		spec_430nm=downsample[9],
		spec_440nm=downsample[10],
		spec_450nm=downsample[11],
		spec_460nm=downsample[12],
		spec_470nm=downsample[13],
		spec_480nm=downsample[14],
		spec_490nm=downsample[15],
		spec_500nm=downsample[16],
		spec_510nm=downsample[17],
		spec_520nm=downsample[18],
		spec_530nm=downsample[19],
		spec_540nm=downsample[20],
		spec_550nm=downsample[21],
		spec_560nm=downsample[22],
		spec_570nm=downsample[23],
		spec_580nm=downsample[24],
		spec_590nm=downsample[25],
		spec_600nm=downsample[26],
		spec_610nm=downsample[27],
		spec_620nm=downsample[28],
		spec_630nm=downsample[29],
		spec_640nm=downsample[30],
		spec_650nm=downsample[31],
		spec_660nm=downsample[32],
		spec_670nm=downsample[33],
		spec_680nm=downsample[34],
		spec_690nm=downsample[35],
		spec_700nm=downsample[36],
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
		
		srgb = convert_color(spectral, sRGBColor)
		rgb_tuple = srgb.get_value_tuple()
		rgb_r = rgb_tuple[0]
		rgb_g = rgb_tuple[1]
		rgb_b = rgb_tuple[2]

	# scale
	rgb_r *= scale
	rgb_g *= scale
	rgb_b *= scale
	rgb_tuple = (rgb_r, rgb_g, rgb_b)
	
	if (args.verbose):
		print("Mode: " + mode)
		print("RGB: " + str(rgb_tuple))

	# prevent out-of-gamut colors -- matplotlib doesn't accept them
	rgb_max = max(*rgb_tuple,1)
	rgb_min = min(*rgb_tuple)
	if ((rgb_max > 1 or rgb_min < 0) and args.warnings):
		print("Warning (spec2rgb): color is out of gamut. For details, use --verbose/-v.")
	rgb_tuple = (max(rgb_r/rgb_max, 0), max(rgb_g/rgb_max, 0), max(rgb_b/rgb_max, 0))

	return(rgb_tuple)

# find radiance/luminance of a spectrum
def spec2radiance(spec, reflect=True):
	# interpolate
	if(len(spec) < 401): spec = np.interp(x_1nm, x_10nm, spec)
	for i in range(401):
		if (not spec[i] > 0): spec[i] = 0

	radiance = spec
	if (reflect): radiance *= wp_watts
	print("Radiance (W/m^2/sr): " + str(sum(radiance)))

	# cd/m^2
	cdm2 = sum(cie_photopic*radiance)*683.002
	print("Luminance (cd/m^2): " + str(cdm2))

	# EDI cd/m^2 equivalent
	d65_lx = 1 / integral
	# 1 unit = 1 lx D65
	d65_opic = d65_lx * sum(d65_1nm*watts2photons*luminosity)
	spec_opic = sum(radiance*watts2photons*luminosity)
	edi = spec_opic / d65_opic
	print("Luminance (EDI/sr): " + str(edi))

	return([sum(radiance), cdm2, edi])

"""
quantum noise

Based on these equations:
* https://journals.biologists.com/jeb/article/218/2/184/14274/Bird-colour-vision-behavioural-thresholds-reveal
* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5043323/
* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1691374/pdf/12816638.pdf
* The Optics of Life, p. 107-112

Meaning of variables:
* v: number of cones per receptive field (spatial summation)
* d: receptor diameter (um)
* f: focal length (um, cancels d)
* D: pupil diameter (um)
* O/Ai/T: transmittance of ocular media ("Ai" is for oil droplets).
This is built in by multiplying by the ocular media array we
created earlier.
* T/dt: integration time (ms), inverse of flicker fusion frequency
(Hz)
* ui/k: optical density of pigment (um^-1)
* l: length of outer segment (um, cancels ui)

This set of variables is (more or less) standard. Terms used in
some papers that I don't include are oil droplets, refractive
indices, dark noise, retinal area covered by cones, and degrees
of the visual field covered by the stimulus. We don't know
enough about the first three here, and the term for visual
angle is usually set as the angle of a single receptive field
determined by the focal length and receptor width. (The (pi/4)^2
comes from this and the pupil size.) For flexibility, I've also
added the option to do it the other way by setting --width and
--dist (and optionally --height). In this case --summation
should be set to the retinal area covered by cones. (The
absolute retinal area isn't relevant.)

The units are confusing: the radiance looks like it has units of
photons/m^3/s/sr, but the "cubic meters" are really per square
meter per wavelength (nm). See The Optics of Life p. 14 and 112
and this StackExchange question:
https://physics.stackexchange.com/questions/690117/what-does-spectral-flux-density-mean-per-wavelength
The Optics of Life advises to multiply the integral by the
wavelength interval, i.e. the bin size, if you find it as a sum
like I do here. This step isn't necessary because we interpolate
to 1 nm so the bin size is 1.

These are the estimates I'm currently using. Where necessary,
I've converted centimeters, millimeters and micrometers to meters.

* dt = 1/22.4: 22.4 Hz = 45 ms, brushtail possum (Signal, Temple & Foster 2011)
* d = 1 / 1000000: 1 um, Didelphis virginiana (Walls 1939)
** This is the width given for rod outer segments (p. 79). Walls doesn't
give the exact width of cones but says "The outer segment is
filamentous, even more slender than that of a rod" (p. 82). He only
describes the width of oil-droplet single cones, so I'm inferring
other cones are similar. In Marmosa mexicana, though, droplet-free
single cones are supposed to be wider.
** Rod width is also confirmed by Kolb & Wang (1985). We don't have the
width of cone outer segments because they cut through at the inner
segment.
* f = 6.81 / 1000: Didelphis aurita (Oswaldo-Cruz, Hokoc & Sousa 1979)
* D = 6 / 1000: Didelphis aurita (")
* k = 0.009: honey possum (Arrese et al. 2002)
* l = 30: Didelphis aurita (Ahnelt et al. 1995)
** Rods and cones are said to have similar dimensions. Note that
the lengths given for rods by Walls (1939) range from 19.7-25.4 (subtracting
the inner segment from the total length).
* K = 0.5 (see The Optics of Life)

Values provided for size:
* 3.8 cm circle (Friedman 1967)
** No distance, they could get as close as they wanted. The choice chamber
was 60 cm long, but there's no mention of a divider.
* 12.5 x 7.5 cm rectangle (Gutierrez et al. 2011)
** As above. We're also told the windows were 19.4 cm apart, which would have
been the limit on comparing them.

Some values for other species for comparison:

* area of retina: 1094 mm^2 (human), 80 (rat; Mayhew & Astle 1997), 450 (cat;
Vaney 1985)
* total number of cones: ~6 million (human)
* outer segment length: 30 um (chicken), 6.2 um (cat; Fisher, Pfeffer & Anderson
1983), 7 um (gray squirrel; "), 25.1-40 um (human), 15.2 um (goldfish red single
cones, Harosi & Flamarique 2012), 5 um (tammar wallaby and fat-tailed dunnart;
Ebeling, Natoli & Hemmi 2010)
** I previously listed 20 for the honey possum (Sumner, Arrese & Partridge),
but after checking their reference, it looks like they used the whole
length instead of the outer segment. Oops.
* outer segment width: 1.73 um (chicken), 6.1 um (goldfish red single
cones), 1.2 um (mouse; Nikonov et al. 2006, Fu & Yau 2007)
* optical density: 0.015 (coral reef triggerfish, birds)
* focal length: 8300 um = 8.3 mm (chickens), 2.6 mm (mouse; Geng et al. 2011),
22.3 mm (human; "), 2.5 mm (triggerfish)
* pupil size: 4900-3500 um = 4.9-3.5 mm (chicken), 2 mm (mouse), 6 mm (human)
* transmittance of ocular media: 80% (general estimate from The Optics of Light)
* integration time/flicker fusion: 50-12 ms = 20-83 Hz (chicken),
70-80 Hz = 14-12 ms (cat and dog), 18 Hz = 56 ms (rat; Gil, Valente & Shemesh 2024),
14 Hz = 71 ms (mouse; "), 60 Hz = 17 ms (human), 25 Hz = 40 ms (some fish)

There was a mistake in here. The formula already includes the absorption
curve (p. 110, 4.12), but then I added it a second time. Look at the final
equation again (p. 111, 4.13). This means we only need the sensitivity array
and k and l can be folded in to it.
"""
abs_r1 = np.zeros(401)
abs_r2 = np.zeros(401)
abs_r3 = np.zeros(401)
abs_r4 = np.zeros(401)

def abs_catch(w, r, rf):
	v = args.summation
	d = args.osd / 1000000
	f = args.fl / 1000
	D = args.pupil / 1000
	K = 0.5
	dt = args.dt
##	k = args.od
##	l = args.osl

	# solid angle
	if (args.width == 0): # receptor
		R = ((d/f)**2)*(math.pi/4)
	else: # stimulus
		if (args.height == 0):
			R = (args.width**2)/(args.sdist**2)
		else:
			R = args.width*args.height/(args.sdist**2)
		if (not args.square): R *= math.pi

	i = w-300
	aq = (math.pi/4)*(R)*(D**2)*K*dt*r[i] * wp_1nm[i] * media_1nm[i]

	"""
	Set the number of receptors per receptive field and degree of
	convergence/spatial summation. By default, the number is
	taken from the list given to the rf (receptive field)
	argument. Setting --summation to a nonzero number
	multiplies the catches by that, and setting it to 0
	sets the receptor numbers to 1 (for the catches but not the
	Weber fractions).
	"""
	if (v != 0): aq *= v*rf
	
	return aq

# fill arrays
for i in range(401):
	w = i+300
	abs_r1[i] = abs_catch(w, r1_1nm, rf[0])
	if (dimension > 1): abs_r2[i] = abs_catch(w, r2_1nm, rf[1])
	if (dimension > 2): abs_r3[i] = abs_catch(w, r3_1nm, rf[2])
	if (dimension > 3): abs_r4[i] = abs_catch(w, r4_1nm, rf[3])

"""
contrast sensitivity (Vorobyev and Osorio (1998); Vorobyev and Osorio (2001))

Without quantum noise, this tends to make unrealistic predictions when the value
for one or more cone signals is very small.
"""
def color_contrast(spec1, spec2, quantum_noise=args.qn, r1=r1_1nm, r2=r2_1nm, abs1_r1=abs_r1, abs1_r2=abs_r2):
	# interpolate 1-nm intervals if we're provided with 10-nm
	if (len(spec1) < 401): spec1 = np.interp(x_1nm, x_10nm, spec1)
	if (len(spec2) < 401): spec2 = np.interp(x_1nm, x_10nm, spec2)
	
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
	for i in range(401):
		if (spec1[i] > 0): # zero/nan check
			qr11 += r1[i] * spec1[i] * wp_1nm[i] * media_1nm[i]
			qr21 += r2[i] * spec1[i] * wp_1nm[i] * media_1nm[i]
			qr31 += r3_1nm[i] * spec1[i] * wp_1nm[i] * media_1nm[i]
			qr41 += r4_1nm[i] * spec1[i] * wp_1nm[i] * media_1nm[i]
		else: spec1[i] = 0 # replace nans and negative numbers with 0
		wpr1 += r1[i] * wp_1nm[i] * media_1nm[i]
		wpr2 += r2[i] * wp_1nm[i] * media_1nm[i]
		wpr3 += r3_1nm[i] * wp_1nm[i] * media_1nm[i]
		wpr4 += r4_1nm[i] * wp_1nm[i] * media_1nm[i]
	for i in range(401):
		w = i + 300
		if (spec2[i] > 0):
			qr12 += r1[i] * spec2[i] * wp_1nm[i] * media_1nm[i]
			qr22 += r2[i] * spec2[i] * wp_1nm[i] * media_1nm[i]
			qr32 += r3_1nm[i] * spec2[i] * wp_1nm[i] * media_1nm[i]
			qr42 += r4_1nm[i] * spec2[i] * wp_1nm[i] * media_1nm[i]
		else: spec2[i] = 0
	
	# normalize
	if (args.vonkries):
		qr11 = qr11 / wpr1
		if (dimension > 1): qr21 = qr21 / wpr2
		if (dimension > 2): qr31 = qr31 / wpr3
		if (dimension > 3): qr41 = qr41 / wpr4
		qr12 = qr12 / wpr1
		if (dimension > 1): qr22 = qr22 / wpr2
		if (dimension > 2): qr32 = qr32 / wpr3
		if (dimension > 3): qr42 = qr42 / wpr4
	
	# differences
	dfr1 = math.log(qr11 / qr12)
	dfr2 = math.log(qr21 / qr22)
	if (dimension > 2): dfr3 = math.log(qr31 / qr32)
	if (dimension > 3): dfr4 = math.log(qr41 / qr42)
	
	# quantum noise
	if (quantum_noise):
		aqr11 = sum(abs1_r1 * spec1)
		aqr21 = sum(abs1_r2 * spec1)
		aqr31 = sum(abs_r3 * spec1)
		aqr41 = sum(abs_r4 * spec1)
		aqr12 = sum(abs1_r1 * spec2)
		aqr22 = sum(abs1_r2 * spec2)
		aqr32 = sum(abs_r3 * spec2)
		aqr42 = sum(abs_r4 * spec2)
		
		er1 = math.sqrt((1 / aqr11 + 1 / aqr12) + 2*wr1**2)
		er2 = math.sqrt((1 / aqr21 + 1 / aqr22) + 2*wr2**2)
		if (dimension > 2): er3 = math.sqrt((1 / aqr31 + 1 / aqr32) + 2*wr3**2)
		if (dimension > 3): er4 = math.sqrt((1 / aqr41 + 1 / aqr42) + 2*wr4**2)
		
		if (args.verbose):
			catch1 = "Absolute catch (1): r1=" + str(aqr11) + ", r2=" + str(aqr21)
			catch2 = "Absolute catch (2): r1=" + str(aqr12) + ", r2=" + str(aqr22)
			noise = "Noise terms: er1=" + str(er1) + ", er2=" + str(er2)
			weber = "Weber fractions: wr1=" + str(wr1) + ", wr2=" + str(wr2)

			if (dimension > 2):
				catch1 += ", r3=" + str(aqr31)
				catch2 += ", r3=" + str(aqr32)
				noise += ", er3=" + str(er3)
				weber += ", wr3=" + str(wr3)
			if (dimension > 3):
				catch1 += ", r4=" + str(aqr41)
				catch2 += ", r4=" + str(aqr42)
				noise += ", er4=" + str(er4)
				weber += ", wr4=" + str(wr4)

			print(catch1)
			print(catch2)
			print(noise)
			print(weber)
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
		+ (er1*er2)**2*(dfr3 - dfr4)**2) / ((er1*er2*er3)**2 + (er1*er2*er4)**2
		+ (er1*er3*er4)**2 + (er2*er3*er4)**2)
	elif (dimension == 3):
		delta_s = (er3**2*(dfr1 - dfr2)**2
		+ er2**2*(dfr1 - dfr3)**2
		+ er1**2*(dfr3 - dfr2)**2) / ((er1*er2)**2 + (er1*er3)**2 + (er2*er3)**2)
	else:
		delta_s = (dfr1 - dfr2)**2 / (er1**2 + er2**2)
	return math.sqrt(delta_s)

# brightness contrast based on https://journals.biologists.com/jeb/article/207/14/2471/14748/Interspecific-and-intraspecific-views-of-color
# I've changed it to a signed value because we care about the direction.
def brightness_contrast(spec1, spec2, r1=luminosity, compare=-1):
	# interpolate 1-nm intervals if we're provided with 10-nm
	if (len(spec1) < 401): spec1 = np.interp(x_1nm, x_10nm, spec1)
	if (len(spec2) < 401): spec2 = np.interp(x_1nm, x_10nm, spec2)

	q1 = sum(r1 * spec1 * wp_1nm)
	if (compare != -1): q2 = compare
	else: q2 = sum(r1 * spec2 * wp_1nm)
	
	df = math.log(q1 / q2)
	delta_s = df / args.weberb # not an absolute value
	return delta_s

"""
plot any spectra

Like --vpt, --plot can plot anything you want, except the original Y values are
left as is.
"""
if (args.plot):
	for i in range(len(args.plot)):
		yvalues = csv2spec(args.plot[i])
		if (args.colorplot): color=spec2rgb(yvalues)
		else: color='k'
		plt.plot(x_1nm, yvalues, color=color)

	plt.xlabel("Wavelength (nm)")
	plt.ylabel(args.ylabel)
	plt.show()

"""
plot visual pigment sensitivity curves

You can plot an arbitrary number of either standard visual pigment templates
(with -r/--receptors), preset curves (--preset) or CSV files (--rcsv). All of
them will be scaled to 1 at their maximum value. They're all plotted as solid
black lines because the default set of Python colors can be counterintuitive
and it's too awkward to assign colors to them that "make sense".
"""
if (args.vpt):
	for i in range(dimension):
		yvalues = rlist[i] * media_1nm
		ymax = max(*yvalues)

		if (args.colorplot):
			spec = np.zeros(401)
			peak = yvalues.argmax(axis=0)
			spec[peak] = 100
			color = spec2rgb(spec, average=True)
		else: color = 'k'
		
		plt.plot(x_1nm, yvalues/ymax, color=color)

	plt.xlabel("Wavelength (nm)")
	plt.ylabel("Relative sensitivity")
	plt.show()

# show luminous efficiency function and lens filtering
# 16/10/2025 changed the argument name from luminosity to sensitivity
# because we have another argument named luminance.
if (args.sensitivity):
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
See this paper for information on Maxwell triangles:
https://www.researchgate.net/publication/8064482_Animal_colour_vision_-_Behavioural_tests_and_physiological_concepts
Maxwell triangles usually are labeled with L on the right and either M on
the left and S on top or the reverse. triangle() puts r1 on the right, r2 on
the left and r3 on top and does not have labels by default. You can add any
labels you want with the argument --labels.

The name "triangle" is slightly misleading because it produces a line or
tetrahedron when provided with 2 or 4+ receptors. (Everything after 4 is
ignored.)

Numbers are off by default because they overlap with everything and are almost
unreadable.

Axes are off by default because these plots don't usually include them and the
Y-axis isn't meaningful for dichromacy.
"""
def triangle(spectra=[], reflect=True, markers=[], colors=[], text=[], legend=False, numbers=False, mec='k', gamut=False, gamutcolor='k', gamutedge='-', achro=True, axes=False):
	# default list elements
	# Checking if they're blank doesn't always work if you call the
	# function twice in a row.
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
		plt.text(-math.sqrt(1/2) - 0.05, -math.sqrt(2/3)/2 - 0.025, args.labels[1])
		plt.text(0 - 0.025, math.sqrt(2/3) + 0.0125, args.labels[2])
		plt.text(math.sqrt(1/2) + 0.0125, -math.sqrt(2/3)/2 - 0.025, args.labels[0])
	else:
		plt.plot([-1, 1], [0, 0], '-k')
		plt.text(-1.1, -.025, args.labels[1])
		plt.text(1.1, -.025, args.labels[0])
		plt.ylim(-1, 1)

	# gamut
	for i in range(41):
		w = i*10 + 300
		wr1 += r1_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		wr2 += r2_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		wr3 += r3_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		wr4 += r4_1nm[i*10] * wp_10nm[i] * media_10nm[i]
	for i in range(41):
		w = i*10 + 300
		labels[i] = w
		r1 = r1_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		r2 = r2_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		r3 = r3_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		r4 = r4_1nm[i*10] * wp_10nm[i] * media_10nm[i]
		if (args.vonkries):
			r1 = r1 / wr1
			if (wr2 > 0): r2 = r2 / wr2
			if (wr3 > 0): r3 = r3 / wr3
			if (wr4 > 0): r4 = r4 / wr4
		total = r1 + r2 + r3 + r4
		if (total > 0): # avoid "invalid value encountered in scalar divide" warning
			r1 = r1 / total
			r2 = r2 / total
			r3 = r3 / total
			r4 = r4 / total
		if (dimension > 2):
			xvalues[i] = math.sqrt(1/2)*(r1 - r2)
			yvalues[i] = math.sqrt(2/3)*(r3 - (r1 + r2)/2)
			if (dimension > 3): zvalues[i] = math.sqrt(3)/2*(r4 - (r1 + r2 + r3)/3)
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
		posxlist = []
		posylist = []
		poszlist = []
		for i in range(401):
			wpr1 += r1_1nm[i] * wp_1nm[i] * media_1nm[i]
			wpr2 += r2_1nm[i] * wp_1nm[i] * media_1nm[i]
			wpr3 += r3_1nm[i] * wp_1nm[i] * media_1nm[i]
			wpr4 += r4_1nm[i] * wp_1nm[i] * media_1nm[i]
		
		for i in range(len(spectra)):
			# dummy multiplier to make sure we don't modify the original array
			spec = spectra[i]*1
			spec_r1 = 0
			spec_r2 = 0
			spec_r3 = 0
			spec_r4 = 0
			
			# interpolate
			if(len(spec) < 401): spec = np.interp(x_1nm, x_10nm, spec)
			
			# combine SPD with illuminant
			if (reflect):
				for j in range(401): spec[j] *= wp_1nm[j]
			
			# remove nans
			for j in range(401):
				if (np.isnan(spec[j])): spec[j] = 0
			
			for j in range(401):
				spec_r1 += r1_1nm[j] * spec[j] * media_1nm[j]
				spec_r2 += r2_1nm[j] * spec[j] * media_1nm[j]
				spec_r3 += r3_1nm[j] * spec[j] * media_1nm[j]
				spec_r4 += r4_1nm[j] * spec[j] * media_1nm[j]
			
			# von Kries transform
			if (args.vonkries):
				spec_r1 = spec_r1 / wpr1
				if (wpr2 > 0): spec_r2 = spec_r2 / wpr2
				if (wpr3 > 0): spec_r3 = spec_r3 / wpr3
				if (wpr4 > 0): spec_r4 = spec_r4 / wpr4

			# declare variables
			posx = 0
			posy = 0
			posz = 0

			if (dimension > 3):
				total = spec_r1 + spec_r2 + spec_r3 + spec_r4
				if (total > 0):
					posx = math.sqrt(1/2)*(spec_r1 - spec_r2) / total
					posy = math.sqrt(2/3)*(spec_r3 - (spec_r1 + spec_r2)/2) / total
					posz = math.sqrt(3)/2*(spec_r4 - (spec_r1 + spec_r2 + spec_r3)/3) / total
				else:
					posx = 0
					posy = 0
					posz = 0
				if (args.verbose):
					print("Relative catches: r1=" + str(spec_r1) + ", r2=" + str(spec_r2) + ", r3=" + str(spec_r3) + ", r4=" + str(spec_r4))
					print("Position in color tetrahedron: " + str((posx, posy, posz)))
			elif (dimension == 3):
				total = spec_r1 + spec_r2 + spec_r3
				if (total > 0):
					posx = math.sqrt(1/2)*(spec_r1 - spec_r2) / total
					posy = math.sqrt(2/3)*(spec_r3 - (spec_r1 + spec_r2)/2) / total
				else:
					posx = 0
					posy = 0
				if (args.verbose):
					print("Relative catches: r1=" + str(spec_r1) + ", r2=" + str(spec_r2) + ", r3=" + str(spec_r3))
					print("Position in color triangle: " + str((posx, posy)))
			elif (dimension == 2):
				if (spec_r1 + spec_r2 > 0): posx = (spec_r1 - spec_r2) / (spec_r1 + spec_r2)
				else: posx = 0
				if (args.verbose):
					print("Relative catches: r1=" + str(spec_r1) + ", r2=" + str(spec_r2))
					print("Position on color line: " + str(posx))
			elif (dimension == 1 and args.verbose):
				# monochromacy
				print("Relative catch: " + str(spec_r1))

			if (dimension > 3):
				if (text[i] == ''):
					ax.plot(posx, posy, posz, marker=markers[i], mec=mec, color=colors[i], linestyle='')
				else:
					ax.plot(posx, posy, posz, marker=markers[i], mec=mec, color=colors[i], label=text[i], linestyle='')
			else:
				if (text[i] == ''):
					plt.plot(posx, posy, marker=markers[i], mec=mec, color=colors[i], linestyle='')
				else:
					plt.plot(posx, posy, marker=markers[i], mec=mec, color=colors[i], label=text[i], linestyle='')

			posxlist.append(posx)
			posylist.append(posy)
			poszlist.append(posz)
	
	if (dimension > 3):
		# positioning text is awkward, I have a better idea...
		ax.plot(*r1v, 'o', color='r')
		ax.plot(*r2v, 'o', color='g')
		ax.plot(*r3v, 'o', color='b')
		ax.plot(*r4v, 'o', color='purple')
	
	if (not axes): plt.axis('off')
	if (legend): plt.legend()
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
			if (args.colorplot): colors.append(spec2rgb(spectrum))
			else: colors.append('k')
		triangle(spectra, colors=colors)

"""
function for modeling color discrimination trials

This outputs a probability for both success (defined as at least the criterion) and failure (less
than the criterion) because the behavioral results include both. We use
different thresholds for "success" for the different sets of discriminations.
In Friedman (1967) the criterion was 35 out of 40 for the color pairs and 68
of 80 for the color/gray pairs. In Gutierrez et al. (2011) significant
performance is defined as more than 62.5% and the total number of trials
for each pair was 64, so the criterion would be 41 of 64. This produces
considerably different significance thresholds when comparing "success" or
"failure" against expected values:

* 35/40 and 68/80: success 0-12 (out of 16), failure 15-16
* 41/64: success 0-8, failure 12-16

For contrast where we don't consider overlap and expect chance (50%)
performance for deltaS < 1, this means that under the first criterion
"success" is defined as at least 10 of 16 distinguishable pairs and "failure"
is defined as at most 13 of 16, whereas under the second criterion "success"
is at least 1 and "failure" is at most 7. There are also several intermediate
values where the probability of both success and failure exceeds 0.05. With
the first two criteria, success becomes more likely than failure at 14
(12 distinguishable), and with the second this occurs at 10.5 (5
distinguishable).

This is now also used directly for optimization graphs, so all the parts
relevant only to testing one set of pigments have been made conditional to
reduce noise and computation.

The running time of color_disc() and brightness_disc() is O(n^2) because we
look at every possible combination. If we want to compare two sets of 4 colors,
as in the above cases, we run through the nested for loop 16 times. There's
no way to change this, so running time can only be cut down by changing what it
does on each pass.

The number of subjects can be factored in with --subjects to get the total
number of trials completed by all of them. (I might could use a
distributive numeral.) This doesn't make much difference for the above
studies because they only had 2.
"""
def color_disc(first, second, correct=35, trials=40, r1=r1_1nm, r2=r2_1nm, abs1_r1=abs_r1, abs1_r2=abs_r2, output=True):
	d = 0 # different
	s = 0 # same
	counter = 0
	combine = len(first)*len(second)
	box = np.empty(combine)
	
	for i in range(len(first)):
		for j in range(len(second)):
			contrast = color_contrast(first[i], second[j], r1=r1, r2=r2, abs1_r1=abs1_r1, abs1_r2=abs1_r2)
			box[counter] = contrast
			if (output):
				print(str(i) + " vs. " + str(j) + ": " + str(contrast))
				if (contrast < 1): s += 1
				else: d += 1
			counter += 1
	if (output):
		print("Different: " + str(d))
		print("Same: " + str(s))
		binomial = binomtest(correct*args.subjects, trials*args.subjects, p=(d + s/2)/combine, alternative='greater')
		binomial2 = binomtest((correct-1)*args.subjects, trials*args.subjects, p=(d + s/2)/combine, alternative='less')
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
def brightness_disc(first, second, correct=35*args.subjects, trials=40*args.subjects, r1=luminosity, output=True, compare=[-1,-1,-1,-1]):
	fb = 0 # first brighter
	sb = 0 # second brighter
	combine = len(first)*len(second)
	# first set
	for i in range(len(first)):
		# second set
		for j in range(len(second)):
			contrast = brightness_contrast(first[i], second[j], r1=r1, compare=compare[j])
			if (contrast >= 1): fb += 1
			elif (contrast <= -1): sb += 1
	equal = combine - fb - sb
	
	binomial = binomtest(correct*args.subjects, trials*args.subjects, (max(fb, sb) + equal/2) / combine, alternative='greater')
	binomial2 = binomtest((correct-1)*args.subjects, trials*args.subjects, (max(fb, sb) + equal/2) / combine, alternative='less')
	if (output):
		print("First brighter: " + str(fb))
		print("Second brighter: " + str(sb))
		print("Equal: " + str(equal))
		print("Expected % correct: " + str(100*(max(fb, sb) + equal/2) / combine)
		      + "% (" + str(100*max(fb, sb) / combine)
		      + "-" + str(100*(max(fb, sb) + equal) / combine) + "%)")
		print("P-value (success): " + str(binomial.pvalue))
		print("P-value (failure): " + str(binomial2.pvalue))
		print("")

	return [binomial.pvalue, binomial2.pvalue]

# Optimize visual pigments for a certain task. Now a function so I don't
# have to keep rewriting it.
def color_opt(first, second, r2=300):
	if (len(args.receptors) > 1): r2 = args.receptors[1]
	w_range = args.receptors[0] - r2 + 1 # inclusive
	mediany = np.empty(w_range)
	miny = np.empty(w_range)
	maxy = np.empty(w_range)
	best_median = 0
	best_mediany = 0
	best_min = 0
	best_miny = 0
	best_max = 0
	best_maxy = 0
	
	for i in range(w_range):
		w = r2 + i
		r2_cur = np.empty(401)
		for j in range(401): r2_cur[j] = vpt(j+300, w)
		abs_r2 = np.empty(401)
		for j in range(401):
			abs_r2[j] = abs_catch(j+300, r2_cur, rf[1])
		if (len(first) > 1 and len(second) > 1):
			values = color_disc(first, second, r2=r2_cur*media_1nm, abs1_r2=abs_r2, output=False)
			if (args.verbose): print(str(w) + " nm: " + str(median[i]))
			mediany[i] = values[0]
			miny[i] = values[1]
			maxy[i] = values[2]
			if (mediany[i] > best_mediany):
				best_mediany = mediany[i]
				best_median = w
			if (miny[i] > best_miny):
				best_miny = miny[i]
				best_min = w
			if (maxy[i] > best_maxy):
				best_maxy = maxy[i]
				best_max = w
		else:
			mediany[i] = color_contrast(first[0], second[0], r2=r2_cur*media_1nm, abs1_r2=abs_r2)
			if (mediany[i] > best_mediany):
				best_mediany = mediany[i]
				best_median = w

	print("Best median: " + str(best_median))
	print("Best min: " + str(best_min))
	print("Best max: " + str(best_max))
	return(mediany, miny, maxy, best_median, best_min, best_max)

# the other big ones
def color_vary(spectrum, r2=300):
	if (len(args.receptors) > 1): r2 = args.receptors[1]
	w_range = args.receptors[0] - r2 + 1 # inclusive

	color = np.empty(w_range)
	for i in range(w_range):
		w = r2 + i
		r2_cur = np.empty(401)
		for j in range(401): r2_cur[j] = vpt(j+300, w)
		if (args.verbose): print(w)
		wpr1 = 0
		wpr2 = 0
		r1_color = 0
		r2_color = 0

		for j in range(401):
			wpr1 += r1_1nm[j] * media_1nm[j] * wp_1nm[j]
			wpr2 += r2_cur[j] * media_1nm[j] * wp_1nm[j]
			if (spectrum[j] > 0):
				r1_color += r1_1nm[j] * media_1nm[j] * spectrum[j] * wp_1nm[j]
				r2_color += r2_cur[j] * media_1nm[j] * spectrum[j] * wp_1nm[j]
		r1_color /= wpr1
		r2_color /= wpr2

		if ((r1_color + r2_color) > 0):
			color[i] = (r1_color - r2_color) / (r1_color + r2_color)
		else: color[i] = 0

	return color

def brightness_vary(spectrum, r2=300):
	if (len(args.receptors) > 1): r2 = args.receptors[1]
	w_range = args.receptors[0] - r2 + 1 # inclusive
	
	brightness = np.empty(w_range)
	for i in range(w_range):
		w = r2 + i
		if (args.verbose): print(w)
		r2_cur = np.empty(401)
		for j in range(401): r2_cur[j] = vpt(j+300, w)
		
		r2_brightness = 0
		for j in range(spectrum.shape[0]):
			r2_brightness += r2_cur[j] * media_1nm[j] * spectrum[j] * wp_1nm[j]
		
		brightness[i] = r2_brightness

	return brightness

# print execution time
#print("%s seconds" % (time.time() - start_time))
