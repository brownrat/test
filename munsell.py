"""
Munsell color cards: 5B, 7.5B, 10YR, 5GY

Reflectance spectra from here: https://www.munsellcolourscienceforpainters.com/MunsellResources/SpectralReflectancesOf2007MunsellBookOfColorGlossy.txt
I can't find any Munsell reflectance spectra that go below 380 nm, so we'l
just have to make do with this.

Note Python variables can't begin with a digit, so the arrays begin with "m"
for Munsell.

Lightness/chroma values used by Gutierrez et al. 2011: 10YR5-8/10, 5B4-7/6,
7.5B5-8/4, 5GY5-8/10

Lighting specifications for absolute quantum catches:
* card dimensions: 12.5 x 7.5 cm
* illuminant: D65 light bulb 1.5 m above the test box, filtered with baking
paper to provide an illuminance of 115 lx
** I can't find any measurements of transmission for baking paper, so we'll
have to assume it transmits all wavelengths equally. Bleached paper probably
doesn't absorb much UV because the lignin is removed.
** How to scale the illuminant to a value in lux is described in central.py.
"""
import central as c
import numpy as np
import matplotlib.pyplot as plt
import math
import time
args = c.args
r1 = 562
r2 = 362
if (len(args.receptors) > 1):
	r1 = args.receptors[0]
	r2 = args.receptors[1]

# 2007 Munsell Book of Color measured with ColorMunki spectrophotometer
# https://www.munsellcolorscienceforpainters.com/MunsellResources/SpectralReflectancesOf2007MunsellBookOfColorGlossy.txt
# These are arranged in rows rather than columns. I pasted every row I wanted
# into a new WordPad document and saved it as a plain text file.
glossy = c.csv2spec('munsell-glossy-2012.csv', numcols=37, h=True)
m5b4 = glossy[0]
m5b5 = glossy[1]
m5b6 = glossy[2]
m5b7 = glossy[3]
m75b5 = glossy[4]
m75b6 = glossy[5]
m75b7 = glossy[6]
m75b8 = glossy[7]
m10yr5 = glossy[8]
m10yr6 = glossy[9]
m10yr7 = glossy[10]
m10yr8 = glossy[11]
m5gy5 = glossy[12]
m5gy6 = glossy[13]
m5gy7 = glossy[14]
m5gy8 = glossy[15]

# cs.joensuu.fi no longer exists. The Perkin-Elmer measurements of matte
# colors can be found here:
# https://sites.uef.fi/spectral/databases-software/munsell-colors-matt-spectrofotometer-measured/
# These are used in this paper:
# https://www.researchgate.net/publication/6764034_Transforming_reflectance_spectra_into_Munsell_color_space_by_using_prime_colors#pf3
# 10YR5/10 and 5GY5-7/10 are missing from this set. Maybe they only come in
# glossy versions?

if (args.matte):
	matte = c.csv2spec("joensuu-matte.csv", numcols=9)
	m5b7 = matte[0]
	m5b6 = matte[1]
	m5b5 = matte[2]
	m5b4 = matte[3]
	m75b8 = matte[4]
	m75b7 = matte[5]
	m75b6 = matte[6]
	m75b5 = matte[7]

# Joensuu measurements of glossy colors, including all 16:
# https://sites.uef.fi/spectral/databases-software/munsell-colors-glossy-all-spectrofotometer-measured/
# This one is tricky because it has 1600 columns. I replaced spaces with commas
# in Notepad, imported it into Google Sheets (renamed to .csv so Google
# accepts it), and pasted the columns I wanted into a new spreadsheet. This
# tool was very useful for finding the right columns:
# https://www.mathsisfun.com/numbers/convert-base.html
if (args.glossy2):
    glossy2 = c.csv2spec('joensuu-glossy.csv', numcols=17)
    m5b4 = glossy2[0]
    m5b5 = glossy2[1]
    m5b6 = glossy2[2]
    m5b7 = glossy2[3]
    m75b5 = glossy2[4]
    m75b6 = glossy2[5]
    m75b7 = glossy2[6]
    m75b8 = glossy2[7]
    m10yr5 = glossy2[8]
    m10yr6 = glossy2[9]
    m10yr7 = glossy2[10]
    m10yr8 = glossy2[11]
    m5gy5 = glossy2[12]
    m5gy6 = glossy2[13]
    m5gy7 = glossy2[14]
    m5gy8 = glossy2[15]

# another alternate version: just remove the reflectance at 380 and 390 nm
# 16/03/2025 -- fixed an error that was causing 380 to be left in: range was
# set to (9,10), but 380 has index 8.
if (args.nouv):
	for i in range(80,91): # 380-390
		m5b4[i] = 0
		m5b5[i] = 0
		m5b6[i] = 0
		m5b7[i] = 0
		m75b5[i] = 0
		m75b6[i] = 0
		m75b7[i] = 0
		m75b8[i] = 0
		m5gy5[i] = 0
		m5gy6[i] = 0
		m5gy7[i] = 0
		m5gy8[i] = 0
		m10yr5[i] = 0
		m10yr6[i] = 0
		m10yr7[i] = 0
		m10yr8[i] = 0

# colors
colors = []
col5b4 = c.spec2rgb(m5b4)
colors.append(col5b4)
col5b5 = c.spec2rgb(m5b5)
colors.append(col5b5)
col5b6 = c.spec2rgb(m5b6)
colors.append(col5b6)
col5b7 = c.spec2rgb(m5b7)
colors.append(col5b7)
col75b5 = c.spec2rgb(m75b5)
colors.append(col75b5)
col75b6 = c.spec2rgb(m75b6)
colors.append(col75b6)
col75b7 = c.spec2rgb(m75b7)
colors.append(col75b7)
col75b8 = c.spec2rgb(m75b8)
colors.append(col75b8)
col10yr5 = c.spec2rgb(m10yr5)
colors.append(col10yr5)
col10yr6 = c.spec2rgb(m10yr6)
colors.append(col10yr6)
col10yr7 = c.spec2rgb(m10yr7)
colors.append(col10yr7)
col10yr8 = c.spec2rgb(m10yr8)
colors.append(col10yr8)
col5gy5 = c.spec2rgb(m5gy5)
colors.append(col5gy5)
col5gy6 = c.spec2rgb(m5gy6)
colors.append(col5gy6)
col5gy7 = c.spec2rgb(m5gy7)
colors.append(col5gy7)
col5gy8 = c.spec2rgb(m5gy8)
colors.append(col5gy8)

if (args.mcheck):
	# plot spectra
	plt.subplot(2, 2, 1)
	plt.plot(c.x_1nm, m5b4, color=col5b4)
	plt.plot(c.x_1nm, m5b5, color=col5b5)
	plt.plot(c.x_1nm, m5b6, color=col5b6)
	plt.plot(c.x_1nm, m5b7, color=col5b7)
	plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("5B")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 2)
	plt.plot(c.x_1nm, m75b5, color=col75b5)
	plt.plot(c.x_1nm, m75b6, color=col75b6)
	plt.plot(c.x_1nm, m75b7, color=col75b7)
	plt.plot(c.x_1nm, m75b8, color=col75b8)
	plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("7.5B")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 3)
	plt.plot(c.x_1nm, m10yr5, color=col10yr5)
	plt.plot(c.x_1nm, m10yr6, color=col10yr6)
	plt.plot(c.x_1nm, m10yr7, color=col10yr7)
	plt.plot(c.x_1nm, m10yr8, color=col10yr8)
	#plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("10YR")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	#plt.show()
	plt.subplot(2, 2, 4)
	plt.plot(c.x_1nm, m5gy5, color=col5gy5)
	plt.plot(c.x_1nm, m5gy6, color=col5gy6)
	plt.plot(c.x_1nm, m5gy7, color=col5gy7)
	plt.plot(c.x_1nm, m5gy8, color=col5gy8)
	#plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("5GY")
	#plt.xlabel("Wavelength (nm)")
	#plt.ylabel("Reflectance")
	plt.show()

	# radiance
	print("5B4/6")
	rad5b4 = c.spec2radiance(m5b4)
	print("5B5/6")
	rad5b5 = c.spec2radiance(m5b5)
	print("5B6/6")
	rad5b6 = c.spec2radiance(m5b6)
	print("5B7/6")
	rad5b7 = c.spec2radiance(m5b7)
	print("7.5B5/4")
	rad75b5 = c.spec2radiance(m75b5)
	print("7.5B6/4")
	rad75b6 = c.spec2radiance(m75b6)
	print("7.5B7/4")
	rad75b7 = c.spec2radiance(m75b7)
	print("7.5B8/4")
	rad75b8 = c.spec2radiance(m75b8)
	print("10YR5/10")
	rad10yr5 = c.spec2radiance(m10yr5)
	print("10YR6/10")
	rad10yr6 = c.spec2radiance(m10yr6)
	print("10YR7/10")
	rad10yr7 = c.spec2radiance(m10yr7)
	print("10YR8/10")
	rad10yr8 = c.spec2radiance(m10yr8)
	print("5GY5/10")
	rad5gy5 = c.spec2radiance(m5gy5)
	print("5GY6/10")
	rad5gy6 = c.spec2radiance(m5gy6)
	print("5GY7/10")
	rad5gy7 = c.spec2radiance(m5gy7)
	print("5GY8/10")
	rad5gy8 = c.spec2radiance(m5gy8)
	
	x = [0,1,2,3]
	plt.subplot(2,1,1)
	plt.plot(x[0], rad5b4[1], '^', color=col5b4, mec='k')
	plt.plot(x[1], rad5b5[1], '^', color=col5b5, mec='k')
	plt.plot(x[2], rad5b6[1], '^', color=col5b6, mec='k')
	plt.plot(x[3], rad5b7[1], '^', color=col5b7, mec='k', label="5B")
	plt.plot(x[0], rad75b5[1], 'o', color=col75b5, mec='k')
	plt.plot(x[1], rad75b6[1], 'o', color=col75b6, mec='k')
	plt.plot(x[2], rad75b7[1], 'o', color=col75b7, mec='k')
	plt.plot(x[3], rad75b8[1], 'o', color=col75b8, mec='k', label="7.5B")
	plt.plot(x[0], rad10yr5[1], 's', color=col10yr5, mec='k')
	plt.plot(x[1], rad10yr6[1], 's', color=col10yr6, mec='k')
	plt.plot(x[2], rad10yr7[1], 's', color=col10yr7, mec='k')
	plt.plot(x[3], rad10yr8[1], 's', color=col10yr8, mec='k', label="10YR")
	plt.plot(x[0], rad5gy5[1], 'D', color=col5gy5, mec='k')
	plt.plot(x[1], rad5gy6[1], 'D', color=col5gy6, mec='k')
	plt.plot(x[2], rad5gy7[1], 'D', color=col5gy7, mec='k')
	plt.plot(x[3], rad5gy8[1], 'D', color=col5gy8, mec='k', label="5GY")
	#plt.xlabel("Filter optical density")
	plt.ylabel("Luminance (cd/m^2)")
	plt.yscale('log')
	plt.legend()
	plt.subplot(2,1,2)
	plt.plot(x[0], rad5b4[2], '^', color=col5b4, mec='k')
	plt.plot(x[1], rad5b5[2], '^', color=col5b5, mec='k')
	plt.plot(x[2], rad5b6[2], '^', color=col5b6, mec='k')
	plt.plot(x[3], rad5b7[2], '^', color=col5b7, mec='k', label="5B")
	plt.plot(x[0], rad75b5[2], 'o', color=col75b5, mec='k')
	plt.plot(x[1], rad75b6[2], 'o', color=col75b6, mec='k')
	plt.plot(x[2], rad75b7[2], 'o', color=col75b7, mec='k')
	plt.plot(x[3], rad75b8[2], 'o', color=col75b8, mec='k', label="7.5B")
	plt.plot(x[0], rad10yr5[2], 's', color=col10yr5, mec='k')
	plt.plot(x[1], rad10yr6[2], 's', color=col10yr6, mec='k')
	plt.plot(x[2], rad10yr7[2], 's', color=col10yr7, mec='k')
	plt.plot(x[3], rad10yr8[2], 's', color=col10yr8, mec='k', label="10YR")
	plt.plot(x[0], rad5gy5[2], 'D', color=col5gy5, mec='k')
	plt.plot(x[1], rad5gy6[2], 'D', color=col5gy6, mec='k')
	plt.plot(x[2], rad5gy7[2], 'D', color=col5gy7, mec='k')
	plt.plot(x[3], rad5gy8[2], 'D', color=col5gy8, mec='k', label="5GY")
	plt.xlabel("Filter optical density")
	plt.ylabel("Luminance (EDI/sr)")
	plt.yscale('log')
	#plt.legend()
	plt.show()

# consolidate brightness levels
m5b_all = [m5b4, m5b5, m5b6, m5b7]
m75b_all = [m75b5, m75b6, m75b7, m75b8]
m5gy_all = [m5gy5, m5gy6, m5gy7, m5gy8]
m10yr_all = [m10yr5, m10yr6, m10yr7, m10yr8]

# trials
if (args.munsell):
	# color space plots
	
	spectra = [m5b4, m5b5, m5b6, m5b7,
		m75b5, m75b6, m75b7, m75b8,
		m10yr5, m10yr6, m10yr7, m10yr8,
		m5gy5, m5gy6, m5gy7, m5gy8]
	markers = []
	for i in range(4): markers.append('^')
	for i in range(4): markers.append('o')
	for i in range(4): markers.append('s')
	for i in range(4): markers.append('D')
	labels = []
	for i in range(len(spectra)): labels.append('')
	labels[3] = '5B'
	labels[7] = '7.5B'
	labels[11] = '10YR'
	labels[15] = '5GY'
	c.triangle(
		spectra=spectra,
		markers=markers,
		colors=colors,
		text=labels,
		legend=True,
		gamut=True,
		gamutcolor="0.7",
		gamutedge='',
		)

	# color contrast
	# Here we compare the expected performance with the study's definition of chance level,
	# 62.5%.
	# 11/02/2025 -- still not sure what definition of significance we should be using.
	# Note that using alternative=greater here will never produce a value less than 0.05
	# whether we compare against 62.5% (10 of 16) or 50% (8 of 16). If using 62.5%, the
	# minimum value (for no distinguishable pairs) is 0.2272491455078125, whereas if
	# using 50%, the minimum is 0.5981903076171875.
	# This has been fixed.

	# 5B vs. 10YR
	print("5B vs. 10YR")
	trial_5b = c.color_disc(m5b_all, m10yr_all, correct=41, trials=64)

	# 7.5B vs. 10YR
	print("7.5B vs. 10YR")
	trial_75b = c.color_disc(m75b_all, m10yr_all, correct=41, trials=64)

	# 5GY vs. 10YR
	print("5GY vs. 10YR")
	trial_5gy = c.color_disc(m5gy_all, m10yr_all, correct=41, trials=64)

	# 10YR vs. 10YR, for comparison
	print("10YR vs. 10YR")
	trial_10yr = c.color_disc(m10yr_all, m10yr_all, correct=41, trials=64)

	# plot contrast values
	boxes = [trial_5b[3], trial_75b[3], trial_5gy[3], trial_10yr[3]]
	labels = ["5B", "7.5B", "5GY", "10YR"]
	bplot = plt.boxplot(x=boxes, tick_labels=labels, patch_artist=True)
	colors = [col5b7, col75b8, col5gy8, col10yr8]
	for patch, color in zip(bplot['boxes'], colors):
		patch.set_facecolor(color)
	plt.ylabel("ΔS")
	plt.show()

# find optimal values for pigments
# This heats up Myotis pretty badly while it's running, but it doesn't take long and is
# vastly preferable to doing it by hand for both time and computation costs.
w_range = r1 - r2 + 1
if (args.mopt1):
	# optimize dichromacy for blue vs. orange with L held constant at specified value and
	# S varied between specified value and L-1
	print("optimizing 5B")

	median_5b = c.color_opt(m5b_all, m10yr_all, r2=r2)[0]
	
	print("")
	print("optimizing 7.5B")

	median_75b = c.color_opt(m75b_all, m10yr_all, r2=r2)[0]
	xvalues = np.empty(w_range)
	for i in range(w_range): xvalues[i] = i + r2
	
	plt.plot(xvalues, median_5b, 'k', label="5B")
	plt.plot(xvalues, median_75b, '--k', label="7.5B")
	plt.xlabel("λmax (nm)")
	plt.ylabel("Median contrast")
	plt.legend()
	plt.show()
	
if (args.mopt2):
	# optimize trichromacy for green and orange vs. orange
	print("optimizing 5GY")
	values_5gy = c.color_opt(m5gy_all, m10yr_all, r2=r2)
	print("optimizing 10YR")
	values_10yr = c.color_opt(m10yr_all, m10yr_all, r2=r2)
	
	xvalues = np.empty(w_range)
	for i in range(w_range): xvalues[i] = i + r2

	# compare
	min_contrast = 0
	for i in range(w_range):
		w = r2 + i
		if (values_5gy[1][i] > values_10yr[2][i]
			and values_5gy[1][i-1] < values_10yr[2][i-1]):
			min_contrast = w
	print("5GY > 10YR: " + str(min_contrast))
	
	plt.plot(xvalues, values_5gy[0], 'k', label="5GY")
	plt.plot(xvalues, values_5gy[1], color='gray')
	plt.plot(xvalues, values_5gy[2], color='gray')
	plt.plot(xvalues, values_10yr[0], '--k', label="10YR")
	plt.plot(xvalues, values_10yr[2], '--', color='gray')
	plt.xlabel("λmax (nm)")
	plt.ylabel("Median contrast")
	plt.legend()
	plt.show()
	
if (args.mopt3):
	# compare positions in color space on the LM axis
	# This seems to be less computationally intensive than finding the contrast
	# values.
	# 06/10/2025 -- rewritten to not use spectral_rendering(), saves 25 seconds:
	# 26.77785038948059−2.019402503967285

	# execution time
	#start_time = time.time()
	xvalues = np.empty(r1 - r2 + 1)
	for i in range(r1 - r2 + 1): xvalues[i] = i+r2

	gy5 = c.color_vary(m5gy5_1nm, r2=r2)
	gy6 = c.color_vary(m5gy6_1nm, r2=r2)
	gy7 = c.color_vary(m5gy7_1nm, r2=r2)
	gy8 = c.color_vary(m5gy8_1nm, r2=r2)
	yr5 = c.color_vary(m10yr5_1nm, r2=r2)
	yr6 = c.color_vary(m10yr6_1nm, r2=r2)
	yr7 = c.color_vary(m10yr7_1nm, r2=r2)
	yr8 = c.color_vary(m10yr8_1nm, r2=r2)

	# overlap
	cross = 0
	for i in range(w_range):
		if (gy5[i] < yr8[i] and gy5[i-1] > yr8[i-1]):
			cross = i + r2
			break

	print("Overlap begins/ends: " + str(cross) + " nm")
	
	plt.plot(xvalues, gy5, 'k', label="5GY5/10")
	plt.plot(xvalues, gy6, '--k', label="5GY6/10")
	plt.plot(xvalues, gy7, ':k', label="5GY7/10")
	plt.plot(xvalues, gy8, '-.k', label="5GY8/10")
	plt.plot(xvalues, yr5, color='gray', label="10YR5/10")
	plt.plot(xvalues, yr6, '--', color='gray', label="10YR6/10")
	plt.plot(xvalues, yr7, ':', color='gray', label="10YR7/10")
	plt.plot(xvalues, yr8, '-.', color='gray', label="10YR8/10")
	plt.xlabel("λmax (nm)")
	plt.ylabel("Chromaticity")
	plt.legend()
	plt.show()
