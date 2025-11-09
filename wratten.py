"""
Kodak Wratten camera filters

These were obtained from the graphs at https://www.kodak.com/en/motion/page/wratten-2-filters/
using Plot Digitizer (https://plotdigitizer.sourceforge.net/). The wavelength
intervals are roughly 5 nm (getting exact integer X values usually isn't
possible).
"""

import central as c
import numpy as np
import math
import matplotlib.pyplot as plt
import statistics
import time
args = c.args
r1 = args.receptors[0]
r2 = 362
if (len(args.receptors) > 1): r2 = args.receptors[1]

# execution time
start_time = time.time()

"""
Alternate data for 400-700 nm from Kodak Photographic Filters Handbook (1990).
There are also graphs going down to 300 that we could scan, but there's no
information for neutral density filters other than 1.0. The biggest difference
below 400 nm evident in these and earlier published graphs is that the Wratten
2 glass version of red 25 transmits a small amount of light in this range
whereas all other sources show no transmission at all. The other four are
consistently shown as either transmitting UV (yellow, blue, ND 1.0) or
not (green).

Yellow, red, blue, green and ND 96 can be found on pages 104, 111, 121, 124
and 132 respectively. When redoing these with Google Sheets, I realized
yellow goes straight from 630 to 650 without a value for 640, so I separated
it.
"""
if (args.w1990):
	yellow15_1nm = c.csv2spec('handbook-filters-yellow.csv')
	others = c.csv2spec('handbook-filters.csv', numcols=5)
	red25_1nm = others[0]
	blue47_1nm = others[1]
	green58_1nm = others[2]
	nd10_1nm = others[3]

	# divide by 100
	yellow15_1nm /= 100
	red25_1nm /= 100
	green58_1nm /= 100
	blue47_1nm /= 100
	nd10_1nm /= 100

	# fill in other ND filters
	nd01_1nm = nd10_1nm * (10**(-.1 + 1))
	nd02_1nm = nd10_1nm * (10**(-.2 + 1))
	nd03_1nm = nd10_1nm * (10**(-.3 + 1))
	nd04_1nm = nd10_1nm * (10**(-.4 + 1))
	nd05_1nm = nd10_1nm * (10**(-.5 + 1))
	nd06_1nm = nd10_1nm * (10**(-.6 + 1))
	nd07_1nm = nd10_1nm * (10**(-.7 + 1))
	nd08_1nm = nd10_1nm * (10**(-.8 + 1))
	nd09_1nm = nd10_1nm * (10**(-.9 + 1))
	nd20_1nm = nd10_1nm * (10**(-2 + 1))
	nd30_1nm = nd10_1nm * (10**(-3 + 1))
	nd40_1nm = nd10_1nm * (10**(-4 + 1))

# Otherwise use the Wratten 2 digitized data. The colors are in two parts
# because the transmission goes to 0 in the middle.
else:
	# yellow 15
	yellow15_part1 = c.csv2spec('yellow15-1.csv')
	yellow15_part2 = c.csv2spec('yellow15-2.csv')
	yellow15_1nm = yellow15_part1 + yellow15_part2

	# red 25
	red25_part1 = c.csv2spec('red25-1.csv')
	red25_part2 = c.csv2spec('red25-2.csv')
	red25_1nm = red25_part1 + red25_part2

	# blue 47
	blue47_part1 = c.csv2spec('blue47-1.csv')
	blue47_part2 = c.csv2spec('blue47-2.csv')
	blue47_1nm = blue47_part1 + blue47_part2

	# green 58
	green58_part1 = c.csv2spec('green58-1.csv')
	green58_part2 = c.csv2spec('green58-2.csv')
	green58_1nm = green58_part1 + green58_part2

	# neutral density 0.1
	nd01_1nm = c.csv2spec('nd01.csv')
	# 0.2
	nd02_1nm = c.csv2spec('nd02-fixed.csv')
	# 0.3
	nd03_1nm = c.csv2spec('nd03.csv')
	# 0.4
	nd04_1nm = c.csv2spec('nd04.csv')
	# 0.5
	nd05_1nm = c.csv2spec('nd05.csv')
	# 0.6
	nd06_1nm = c.csv2spec('nd06.csv')
	# 0.7
	nd07_1nm = c.csv2spec('nd07.csv')
	# 0.8
	nd08_1nm = c.csv2spec('nd08.csv')
	# 0.9
	nd09_1nm = c.csv2spec('nd09.csv')
	# 1.0
	nd10_1nm = c.csv2spec('nd10.csv')
	# 2.0
	nd20_1nm = c.csv2spec('nd20.csv')
	# 3.0
	nd30_1nm = c.csv2spec('nd30.csv')
	# 4.0
	nd40_1nm = c.csv2spec('nd40.csv')

# other replacement test
if (args.ndswap):
	nd01_1nm = np.full(c.range1, 10**(-.1))
	nd02_1nm = np.full(c.range1, 10**(-.2))
	nd03_1nm = np.full(c.range1, 10**(-.3))
	nd04_1nm = np.full(c.range1, 10**(-.4))
	nd05_1nm = np.full(c.range1, 10**(-.5))
	nd06_1nm = np.full(c.range1, 10**(-.6))
	nd07_1nm = np.full(c.range1, 10**(-.7))
	nd08_1nm = np.full(c.range1, 10**(-.8))
	nd09_1nm = np.full(c.range1, 10**(-.9))
	nd10_1nm = np.full(c.range1, 10**(-1))
	nd20_1nm = np.full(c.range1, 10**(-2))
	nd30_1nm = np.full(c.range1, 10**(-3))
	nd40_1nm = np.full(c.range1, 10**(-4))

"""
Find brightness matches for human vision. Since we begin with colors that are much
brighter than blue 47, I intended to choose the "brightest" ND filter or pair of filters that produces
an achromatic contrast less than 1, but in practice there is be only one match for each
color. This depends on the specified vision system and illuminant like everything else, so
we have to set --media human2 --white i (L, M and S don't matter). The current matches are:
* red: 1.0 + 0.5
* yellow: 2.0
* green: 1.0 + 0.3
* gray: 2.0 + 0.1
When using the 1990 data, the matches are slightly brighter:
* red: 1.0 + 0.3
* yellow: 1.0 + 0.8
* green: 1.0 + 0.2
* gray: 1.0 + 0.9
"""
if (args.wcheck):	
	# "dummy" array for gray as we assume no color filter
	placeholder = np.ones(c.range1)
	
	# brightness tests -- wrap this up in a function so we don't have to copy it a million
	# times
	def brightness_tests(f):
		f_01 = np.empty(c.range1)
		for i in range(c.range1):
			f_01[i] = f[i] * nd01_1nm[i]
		f_02 = np.empty(c.range1)
		for i in range(c.range1):
			f_02[i] = f[i] * nd02_1nm[i]
		f_03 = np.empty(c.range1)
		for i in range(c.range1):
			f_03[i] = f[i] * nd03_1nm[i]
		f_04 = np.empty(c.range1)
		for i in range(c.range1):
			f_04[i] = f[i] * nd04_1nm[i]
		f_05 = np.empty(c.range1)
		for i in range(c.range1):
			f_05[i] = f[i] * nd05_1nm[i]
		f_06 = np.empty(c.range1)
		for i in range(c.range1):
			f_06[i] = f[i] * nd06_1nm[i]
		f_07 = np.empty(c.range1)
		for i in range(c.range1):
			f_07[i] = f[i] * nd07_1nm[i]
		f_08 = np.empty(c.range1)
		for i in range(c.range1):
			f_08[i] = f[i] * nd08_1nm[i]
		f_09 = np.empty(c.range1)
		for i in range(c.range1):
			f_09[i] = f[i] * nd09_1nm[i]
		# 1.0 + others
		f_10 = np.empty(c.range1)
		for i in range(c.range1):
			f_10[i] = f[i] * nd10_1nm[i]
		f_1001 = np.empty(c.range1)
		for i in range(c.range1):
			f_1001[i] = f[i] * nd10_1nm[i] * nd01_1nm[i]
		f_1002 = np.empty(c.range1)
		for i in range(c.range1):
			f_1002[i] = f[i] * nd10_1nm[i] * nd02_1nm[i]
		f_1003 = np.empty(c.range1)
		for i in range(c.range1):
			f_1003[i] = f[i] * nd10_1nm[i] * nd03_1nm[i]
		f_1004 = np.empty(c.range1)
		for i in range(c.range1):
			f_1004[i] = f[i] * nd10_1nm[i] * nd04_1nm[i]
		f_1005 = np.empty(c.range1)
		for i in range(c.range1):
			f_1005[i] = f[i] * nd10_1nm[i] * nd05_1nm[i]
		f_1006 = np.empty(c.range1)
		for i in range(c.range1):
			f_1006[i] = f[i] * nd10_1nm[i] * nd06_1nm[i]
		f_1007 = np.empty(c.range1)
		for i in range(c.range1):
			f_1007[i] = f[i] * nd10_1nm[i] * nd07_1nm[i]
		f_1008 = np.empty(c.range1)
		for i in range(c.range1):
			f_1008[i] = f[i] * nd10_1nm[i] * nd08_1nm[i]
		f_1009 = np.empty(c.range1)
		for i in range(c.range1):
			f_1009[i] = f[i] * nd10_1nm[i] * nd09_1nm[i]
		# 2.0
		f_20 = np.empty(c.range1)
		for i in range(c.range1):
			f_20[i] = f[i] * nd20_1nm[i]
		f_2001 = np.empty(c.range1)
		for i in range(c.range1):
			f_2001[i] = f[i] * nd20_1nm[i] * nd01_1nm[i]
		f_2002 = np.empty(c.range1)
		for i in range(c.range1):
			f_2002[i] = f[i] * nd20_1nm[i] * nd02_1nm[i]
		f_2003 = np.empty(c.range1)
		for i in range(c.range1):
			f_2003[i] = f[i] * nd20_1nm[i] * nd03_1nm[i]
		f_2004 = np.empty(c.range1)
		for i in range(c.range1):
			f_2004[i] = f[i] * nd20_1nm[i] * nd04_1nm[i]
		f_2005 = np.empty(c.range1)
		for i in range(c.range1):
			f_2005[i] = f[i] * nd20_1nm[i] * nd05_1nm[i]
		f_2006 = np.empty(c.range1)
		for i in range(c.range1):
			f_2006[i] = f[i] * nd20_1nm[i] * nd06_1nm[i]
		f_2007 = np.empty(c.range1)
		for i in range(c.range1):
			f_2007[i] = f[i] * nd20_1nm[i] * nd07_1nm[i]
		f_2008 = np.empty(c.range1)
		for i in range(c.range1):
			f_2008[i] = f[i] * nd20_1nm[i] * nd08_1nm[i]
		f_2009 = np.empty(c.range1)
		for i in range(c.range1):
			f_2009[i] = f[i] * nd20_1nm[i] * nd09_1nm[i]
		# 3.0
		f_30 = np.empty(c.range1)
		for i in range(c.range1):
			f_30[i] = f[i] * nd30_1nm[i]
		f_3001 = np.empty(c.range1)
		for i in range(c.range1):
			f_3001[i] = f[i] * nd30_1nm[i] * nd01_1nm[i]
		f_3002 = np.empty(c.range1)
		for i in range(c.range1):
			f_3002[i] = f[i] * nd30_1nm[i] * nd02_1nm[i]
		f_3003 = np.empty(c.range1)
		for i in range(c.range1):
			f_3003[i] = f[i] * nd30_1nm[i] * nd03_1nm[i]
		f_3004 = np.empty(c.range1)
		for i in range(c.range1):
			f_3004[i] = f[i] * nd30_1nm[i] * nd04_1nm[i]
		f_3005 = np.empty(c.range1)
		for i in range(c.range1):
			f_3005[i] = f[i] * nd30_1nm[i] * nd05_1nm[i]
		f_3006 = np.empty(c.range1)
		for i in range(c.range1):
			f_3006[i] = f[i] * nd30_1nm[i] * nd06_1nm[i]
		f_3007 = np.empty(c.range1)
		for i in range(c.range1):
			f_3007[i] = f[i] * nd30_1nm[i] * nd07_1nm[i]
		f_3008 = np.empty(c.range1)
		for i in range(c.range1):
			f_3008[i] = f[i] * nd30_1nm[i] * nd08_1nm[i]
		f_3009 = np.empty(c.range1)
		for i in range(c.range1):
			f_3009[i] = f[i] * nd30_1nm[i] * nd09_1nm[i]
		# 4.0
		f_40 = np.empty(c.range1)
		for i in range(c.range1):
			f_40[i] = f[i] * nd40_1nm[i]
		f_4001 = np.empty(c.range1)
		for i in range(c.range1):
			f_4001[i] = f[i] * nd40_1nm[i] * nd01_1nm[i]
		f_4002 = np.empty(c.range1)
		for i in range(c.range1):
			f_4002[i] = f[i] * nd40_1nm[i] * nd02_1nm[i]
		f_4003 = np.empty(c.range1)
		for i in range(c.range1):
			f_4003[i] = f[i] * nd40_1nm[i] * nd03_1nm[i]
		f_4004 = np.empty(c.range1)
		for i in range(c.range1):
			f_4004[i] = f[i] * nd40_1nm[i] * nd04_1nm[i]
		f_4005 = np.empty(c.range1)
		for i in range(c.range1):
			f_4005[i] = f[i] * nd40_1nm[i] * nd05_1nm[i]
		f_4006 = np.empty(c.range1)
		for i in range(c.range1):
			f_4006[i] = f[i] * nd40_1nm[i] * nd06_1nm[i]
		f_4007 = np.empty(c.range1)
		for i in range(c.range1):
			f_4007[i] = f[i] * nd40_1nm[i] * nd07_1nm[i]
		f_4008 = np.empty(c.range1)
		for i in range(c.range1):
			f_4008[i] = f[i] * nd40_1nm[i] * nd08_1nm[i]
		f_4009 = np.empty(c.range1)
		for i in range(c.range1):
			f_4009[i] = f[i] * nd40_1nm[i] * nd09_1nm[i]
		
		print("Brightness contrast with blue:")
		print("0.1: " + str(c.brightness_contrast(f_01, blue47_1nm)))
		print("0.2: " + str(c.brightness_contrast(f_02, blue47_1nm)))
		print("0.3: " + str(c.brightness_contrast(f_03, blue47_1nm)))
		print("0.4: " + str(c.brightness_contrast(f_04, blue47_1nm)))
		print("0.5: " + str(c.brightness_contrast(f_05, blue47_1nm)))
		print("0.6: " + str(c.brightness_contrast(f_06, blue47_1nm)))
		print("0.7: " + str(c.brightness_contrast(f_07, blue47_1nm)))
		print("0.8: " + str(c.brightness_contrast(f_08, blue47_1nm)))
		print("0.9: " + str(c.brightness_contrast(f_09, blue47_1nm)))
		print("1.0: " + str(c.brightness_contrast(f_10, blue47_1nm)))
		print("1.0 + 0.1: " + str(c.brightness_contrast(f_1001, blue47_1nm)))
		print("1.0 + 0.2: " + str(c.brightness_contrast(f_1002, blue47_1nm)))
		print("1.0 + 0.3: " + str(c.brightness_contrast(f_1003, blue47_1nm)))
		print("1.0 + 0.4: " + str(c.brightness_contrast(f_1004, blue47_1nm)))
		print("1.0 + 0.5: " + str(c.brightness_contrast(f_1005, blue47_1nm)))
		print("1.0 + 0.6: " + str(c.brightness_contrast(f_1006, blue47_1nm)))
		print("1.0 + 0.7: " + str(c.brightness_contrast(f_1007, blue47_1nm)))
		print("1.0 + 0.8: " + str(c.brightness_contrast(f_1008, blue47_1nm)))
		print("1.0 + 0.9: " + str(c.brightness_contrast(f_1009, blue47_1nm)))
		print("2.0: " + str(c.brightness_contrast(f_20, blue47_1nm)))
		print("2.0 + 0.1: " + str(c.brightness_contrast(f_2001, blue47_1nm)))
		print("2.0 + 0.2: " + str(c.brightness_contrast(f_2002, blue47_1nm)))
		print("2.0 + 0.3: " + str(c.brightness_contrast(f_2003, blue47_1nm)))
		print("2.0 + 0.4: " + str(c.brightness_contrast(f_2004, blue47_1nm)))
		print("2.0 + 0.5: " + str(c.brightness_contrast(f_2005, blue47_1nm)))
		print("2.0 + 0.6: " + str(c.brightness_contrast(f_2006, blue47_1nm)))
		print("2.0 + 0.7: " + str(c.brightness_contrast(f_2007, blue47_1nm)))
		print("2.0 + 0.8: " + str(c.brightness_contrast(f_2008, blue47_1nm)))
		print("2.0 + 0.9: " + str(c.brightness_contrast(f_2009, blue47_1nm)))
		print("3.0: " + str(c.brightness_contrast(f_30, blue47_1nm)))
		print("3.0 + 0.1: " + str(c.brightness_contrast(f_3001, blue47_1nm)))
		print("3.0 + 0.2: " + str(c.brightness_contrast(f_3002, blue47_1nm)))
		print("3.0 + 0.3: " + str(c.brightness_contrast(f_3003, blue47_1nm)))
		print("3.0 + 0.4: " + str(c.brightness_contrast(f_3004, blue47_1nm)))
		print("3.0 + 0.5: " + str(c.brightness_contrast(f_3005, blue47_1nm)))
		print("3.0 + 0.6: " + str(c.brightness_contrast(f_3006, blue47_1nm)))
		print("3.0 + 0.7: " + str(c.brightness_contrast(f_3007, blue47_1nm)))
		print("3.0 + 0.8: " + str(c.brightness_contrast(f_3008, blue47_1nm)))
		print("3.0 + 0.9: " + str(c.brightness_contrast(f_3009, blue47_1nm)))
		print("4.0: " + str(c.brightness_contrast(f_40, blue47_1nm)))
		print("4.0 + 0.1: " + str(c.brightness_contrast(f_4001, blue47_1nm)))
		print("4.0 + 0.2: " + str(c.brightness_contrast(f_4002, blue47_1nm)))
		print("4.0 + 0.3: " + str(c.brightness_contrast(f_4003, blue47_1nm)))
		print("4.0 + 0.4: " + str(c.brightness_contrast(f_4004, blue47_1nm)))
		print("4.0 + 0.5: " + str(c.brightness_contrast(f_4005, blue47_1nm)))
		print("4.0 + 0.6: " + str(c.brightness_contrast(f_4006, blue47_1nm)))
		print("4.0 + 0.7: " + str(c.brightness_contrast(f_4007, blue47_1nm)))
		print("4.0 + 0.8: " + str(c.brightness_contrast(f_4008, blue47_1nm)))
		print("4.0 + 0.9: " + str(c.brightness_contrast(f_4009, blue47_1nm)))
	
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

# These are used by several different parts of this file, so we keep them
# out here.

# red 25
red25_0 = np.empty(c.range1)
for i in range(c.range1):
	if (args.w1990):
		red25_0[i] = red25_1nm[i] * nd10_1nm[i] * nd03_1nm[i]
	else:
		red25_0[i] = red25_1nm[i] * nd10_1nm[i] * nd05_1nm[i]
red25_03 = np.empty(c.range1)
for i in range(c.range1):
	red25_03[i] = red25_0[i] * nd03_1nm[i]
red25_07 = np.empty(c.range1)
for i in range(c.range1):
	red25_07[i] = red25_0[i] * nd07_1nm[i]
red25_10 = np.empty(c.range1)
for i in range(c.range1):
	red25_10[i] = red25_0[i] * nd10_1nm[i]

# yellow 15
yellow15_0 = np.empty(c.range1)
for i in range(c.range1):
	if (args.w1990 and not args.ndswap):
		yellow15_0[i] = yellow15_1nm[i] * nd10_1nm[i] * nd08_1nm[i]
	elif (args.ndswap):
		yellow15_0[i] = yellow15_1nm[i] * nd10_1nm[i] * nd09_1nm[i]
	else: yellow15_0[i] = yellow15_1nm[i] * nd20_1nm[i]
yellow15_03 = np.empty(c.range1)
for i in range(c.range1):
	yellow15_03[i] = yellow15_0[i] * nd03_1nm[i]
yellow15_07 = np.empty(c.range1)
for i in range(c.range1):
	yellow15_07[i] = yellow15_0[i] * nd07_1nm[i]
yellow15_10 = np.empty(c.range1)
for i in range(c.range1):
	yellow15_10[i] = yellow15_0[i] * nd10_1nm[i]

# green 58
green58_0 = np.empty(c.range1)
for i in range(c.range1):
	if (args.w1990 and not args.ndswap):
		green58_0[i] = green58_1nm[i] * nd10_1nm[i] * nd02_1nm[i]
	else: green58_0[i] = green58_1nm[i] * nd10_1nm[i] * nd03_1nm[i]
green58_03 = np.empty(c.range1)
for i in range(c.range1):
	green58_03[i] = green58_0[i] * nd03_1nm[i]
green58_07 = np.empty(c.range1)
for i in range(c.range1):
	green58_07[i] = green58_0[i] * nd07_1nm[i]
green58_10 = np.empty(c.range1)
for i in range(c.range1):
	green58_10[i] = green58_0[i] * nd10_1nm[i]

# blue 47
blue47_0 = blue47_1nm # keep this name for convenience
blue47_03 = np.empty(c.range1)
for i in range(c.range1):
	blue47_03[i] = blue47_0[i] * nd03_1nm[i]
blue47_07 = np.empty(c.range1)
for i in range(c.range1):
	blue47_07[i] = blue47_0[i] * nd07_1nm[i]
blue47_10 = np.empty(c.range1)
for i in range(c.range1):
	blue47_10[i] = blue47_0[i] * nd10_1nm[i]

# gray
gray_0 = np.empty(c.range1)
for i in range(c.range1):
          if (args.w1990 and not args.ndswap):
              gray_0[i] = nd10_1nm[i] * nd09_1nm[i]
          elif (args.ndswap): gray_0[i] = nd20_1nm[i]
          else: gray_0[i] = nd20_1nm[i] * nd01_1nm[i]

# hypothetical "gray" with a luminance exactly halfway between yellow and green
compare = [-1,-1,-1,-1]
if (args.grayswap): compare = [
        math.sqrt(sum(green58_0*c.luminosity*c.wp_1nm) * sum(yellow15_0*c.luminosity*c.wp_1nm)),
        math.sqrt(sum(green58_03*c.luminosity*c.wp_1nm) * sum(yellow15_03*c.luminosity*c.wp_1nm)),
        math.sqrt(sum(green58_07*c.luminosity*c.wp_1nm) * sum(yellow15_07*c.luminosity*c.wp_1nm)),
        math.sqrt(sum(green58_10*c.luminosity*c.wp_1nm) * sum(yellow15_10*c.luminosity*c.wp_1nm))
        ]
gray_03 = np.empty(c.range1)
for i in range(c.range1):
	gray_03[i] = gray_0[i] * nd03_1nm[i]
gray_07 = np.empty(c.range1)
for i in range(c.range1):
	gray_07[i] = gray_0[i] * nd07_1nm[i]
gray_10 = np.empty(c.range1)
for i in range(c.range1):
	gray_10[i] = gray_0[i] * nd10_1nm[i]

# colors
scale = 6
cr0 = c.spec2rgb(red25_0, scale=scale)
cr03 = c.spec2rgb(red25_03, scale=scale)
cr07 = c.spec2rgb(red25_07, scale=scale)
cr10 = c.spec2rgb(red25_10, scale=scale)
cy0 = c.spec2rgb(yellow15_0, scale=scale)
cy03 = c.spec2rgb(yellow15_03, scale=scale)
cy07 = c.spec2rgb(yellow15_07, scale=scale)
cy10 = c.spec2rgb(yellow15_10, scale=scale)
cg0 = c.spec2rgb(green58_0, scale=scale)
cg03 = c.spec2rgb(green58_03, scale=scale)
cg07 = c.spec2rgb(green58_07, scale=scale)
cg10 = c.spec2rgb(green58_10, scale=scale)
cb0 = c.spec2rgb(blue47_0, scale=scale)
cb03 = c.spec2rgb(blue47_03, scale=scale)
cb07 = c.spec2rgb(blue47_07, scale=scale)
cb10 = c.spec2rgb(blue47_10, scale=scale)
cgray0 = c.spec2rgb(gray_0, scale=scale)
cgray03 = c.spec2rgb(gray_03, scale=scale)
cgray07 = c.spec2rgb(gray_07, scale=scale)
cgray10 = c.spec2rgb(gray_10, scale=scale)

# other checks
if (args.wcheck2):
	# plot
	plt.subplot(3, 2, 1)
	plt.plot(c.x_1nm, red25_0*100, color=cr0)
	plt.plot(c.x_1nm, red25_03*100, color=cr03)
	plt.plot(c.x_1nm, red25_07*100, color=cr07)
	plt.plot(c.x_1nm, red25_10*100, color=cr10)
	# problem with overlapping text
	plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("Red 25")
	plt.subplot(3, 2, 2)
	plt.plot(c.x_1nm, yellow15_0*100, color=cy0)
	plt.plot(c.x_1nm, yellow15_03*100, color=cy03)
	plt.plot(c.x_1nm, yellow15_07*100, color=cy07)
	plt.plot(c.x_1nm, yellow15_10*100, color=cy10)
	plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("Yellow 15")
	plt.subplot(3, 2, 3)
	plt.plot(c.x_1nm, green58_0*100, color=cg0)
	plt.plot(c.x_1nm, green58_03*100, color=cg03)
	plt.plot(c.x_1nm, green58_07*100, color=cg07)
	plt.plot(c.x_1nm, green58_10*100, color=cg10)
	plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("Green 58")
	plt.subplot(3, 2, 4)
	plt.plot(c.x_1nm, blue47_0*100, color=cb0)
	plt.plot(c.x_1nm, blue47_03*100, color=cb03)
	plt.plot(c.x_1nm, blue47_07*100, color=cb07)
	plt.plot(c.x_1nm, blue47_10*100, color=cb10)
	plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("Blue 47")
	plt.subplot(3, 2, 5)
	plt.plot(c.x_1nm, gray_0*100, color=cgray0)
	plt.plot(c.x_1nm, gray_03*100, color=cgray03)
	plt.plot(c.x_1nm, gray_07*100, color=cgray07)
	plt.plot(c.x_1nm, gray_10*100, color=cgray10)
	plt.gca().axes.get_xaxis().set_ticklabels([])
	plt.title("Gray")
	plt.show()

	# radiance/luminance
	print("Red 25")
	radr0 = c.spec2radiance(red25_0)
	print("+ 0.3")
	radr03 = c.spec2radiance(red25_03)
	print("+ 0.7")
	radr07 = c.spec2radiance(red25_07)
	print("+ 1.0")
	radr10 = c.spec2radiance(red25_10)
	print("")
	print("Yellow 15")
	rady0 = c.spec2radiance(yellow15_0)
	print("+ 0.3")
	rady03 = c.spec2radiance(yellow15_03)
	print("+ 0.7")
	rady07 = c.spec2radiance(yellow15_07)
	print("+ 1.0")
	rady10 = c.spec2radiance(yellow15_10)
	print("")
	print("Green 58")
	radg0 = c.spec2radiance(green58_0)
	print("+ 0.3")
	radg03 = c.spec2radiance(green58_03)
	print("+ 0.7")
	radg07 = c.spec2radiance(green58_07)
	print("+ 1.0")
	radg10 = c.spec2radiance(green58_10)
	print("")
	print("Blue 47")
	radb0 = c.spec2radiance(blue47_0)
	print("+ 0.3")
	radb03 = c.spec2radiance(blue47_03)
	print("+ 0.7")
	radb07 = c.spec2radiance(blue47_07)
	print("+ 1.0")
	radb10 = c.spec2radiance(blue47_10)
	print("")
	print("Gray")
	radgray0 = c.spec2radiance(gray_0)
	print("+ 0.3")
	radgray03 = c.spec2radiance(gray_03)
	print("+ 0.7")
	radgray07 = c.spec2radiance(gray_07)
	print("+ 1.0")
	radgray10 = c.spec2radiance(gray_10)

	# relative brightness
	x = [0.0, 0.3, 0.7, 1.0]
	plt.subplot(2,1,1)
	plt.plot(x[0], radr0[1], 's', color=cr0, mec='k', label="Red")
	plt.plot(x[1], radr03[1], 's', color=cr03, mec='k')
	plt.plot(x[2], radr07[1], 's', color=cr07, mec='k')
	plt.plot(x[3], radr10[1], 's', color=cr10, mec='k')
	plt.plot(x[0], rady0[1], 'D', color=cy0, mec='k', label="Yellow")
	plt.plot(x[1], rady03[1], 'D', color=cy03, mec='k')
	plt.plot(x[2], rady07[1], 'D', color=cy07, mec='k')
	plt.plot(x[3], rady10[1], 'D', color=cy10, mec='k')
	plt.plot(x[0], radg0[1], '^', color=cg0, mec='k', label="Green")
	plt.plot(x[1], radg03[1], '^', color=cg03, mec='k')
	plt.plot(x[2], radg07[1], '^', color=cg07, mec='k')
	plt.plot(x[3], radg10[1], '^', color=cg10, mec='k')
	plt.plot(x[0], radb0[1], 'o', color=cb0, mec='k', label="Blue")
	plt.plot(x[1], radb03[1], 'o', color=cb03, mec='k')
	plt.plot(x[2], radb07[1], 'o', color=cb07, mec='k')
	plt.plot(x[3], radb10[1], 'o', color=cb10, mec='k')
	plt.plot(x[0], radgray0[1], 'v', color=cgray0, mec='k', label="Gray")
	plt.plot(x[1], radgray03[1], 'v', color=cgray03, mec='k')
	plt.plot(x[2], radgray07[1], 'v', color=cgray07, mec='k')
	plt.plot(x[3], radgray10[1], 'v', color=cgray10, mec='k')
	#plt.xlabel("Filter optical density")
	plt.ylabel("Luminance (cd/m^2)")
	plt.yscale('log')
	plt.legend()
	plt.subplot(2,1,2)
	plt.plot(x[0], radr0[2], 's', color=cr0, mec='k', label="Red")
	plt.plot(x[1], radr03[2], 's', color=cr03, mec='k')
	plt.plot(x[2], radr07[2], 's', color=cr07, mec='k')
	plt.plot(x[3], radr10[2], 's', color=cr10, mec='k')
	plt.plot(x[0], rady0[2], 'D', color=cy0, mec='k', label="Yellow")
	plt.plot(x[1], rady03[2], 'D', color=cy03, mec='k')
	plt.plot(x[2], rady07[2], 'D', color=cy07, mec='k')
	plt.plot(x[3], rady10[2], 'D', color=cy10, mec='k')
	plt.plot(x[0], radg0[2], '^', color=cg0, mec='k', label="Green")
	plt.plot(x[1], radg03[2], '^', color=cg03, mec='k')
	plt.plot(x[2], radg07[2], '^', color=cg07, mec='k')
	plt.plot(x[3], radg10[2], '^', color=cg10, mec='k')
	plt.plot(x[0], radb0[2], 'o', color=cb0, mec='k', label="Blue")
	plt.plot(x[1], radb03[2], 'o', color=cb03, mec='k')
	plt.plot(x[2], radb07[2], 'o', color=cb07, mec='k')
	plt.plot(x[3], radb10[2], 'o', color=cb10, mec='k')
	plt.plot(x[0], radgray0[1], 'v', color=cgray0, mec='k', label="Gray")
	plt.plot(x[1], radgray03[1], 'v', color=cgray03, mec='k')
	plt.plot(x[2], radgray07[1], 'v', color=cgray07, mec='k')
	plt.plot(x[3], radgray10[1], 'v', color=cgray10, mec='k')
	plt.xlabel("Filter optical density")
	plt.ylabel("Luminance (EDI/sr)")
	plt.yscale('log')
	#plt.legend()
	plt.show()

# arrays containing all 4 variants of the colors
red25_all = [red25_0, red25_03, red25_07, red25_10]
yellow15_all = [yellow15_0, yellow15_03, yellow15_07, yellow15_10]
green58_all = [green58_0, green58_03, green58_07, green58_10]
blue47_all = [blue47_0, blue47_03, blue47_07, blue47_10]
gray_all = [gray_0, gray_03, gray_07, gray_10]

if (args.wratten):
	# red-yellow
	print("R-Y")
	c.brightness_disc(red25_all, yellow15_all)
	
	# red-green
	print("R-G")
	c.brightness_disc(red25_all, green58_all)
	
	# red-blue
	print("R-B")
	c.brightness_disc(red25_all, blue47_all)
	
	# yellow-green
	print("Y-G")
	c.brightness_disc(yellow15_all, green58_all)
	
	# yellow-blue
	print("Y-B")
	c.brightness_disc(yellow15_all, blue47_all)
	
	# green-blue
	print("G-B")
	c.brightness_disc(green58_all, blue47_all)
	
	# colors vs. gray
	
	# red
	print("Red vs. gray")
	c.brightness_disc(red25_all, gray_all, correct=68, trials=80, compare=compare)
	
	# yellow
	print("Yellow vs. gray")
	c.brightness_disc(yellow15_all, gray_all, correct=68, trials=80, compare=compare)
	
	# green
	print("Green vs. gray")
	c.brightness_disc(green58_all, gray_all, correct=68, trials=80, compare=compare)
	
	# blue
	print("Blue vs. gray")
	c.brightness_disc(blue47_all, gray_all, correct=68, trials=80, compare=compare)

	# print execution time
	print("%s seconds" % (time.time() - start_time))
	
	# color differences -- this is the hard part...
	if (c.dimension > 1):
		# color space plot
		spectra = [red25_10, red25_07, red25_03, red25_0,
			yellow15_10, yellow15_07, yellow15_03, yellow15_0,
			green58_10, green58_07, green58_03, green58_0,
			blue47_10, blue47_07, blue47_03, blue47_0,
			gray_10, gray_07, gray_03, gray_0]
		markers = []
		for i in range(4): markers.append('s')
		for i in range(4): markers.append('D')
		for i in range(4): markers.append('^')
		for i in range(4): markers.append('o')
		for i in range(4): markers.append('v')
		colors = []
		for i in range(len(spectra)): colors.append(c.spec2rgb(spectra[i],scale=scale))
		labels = []
		for i in range(len(spectra)): labels.append('')
		labels[3] = 'Red'
		labels[7] = 'Yellow'
		labels[11] = 'Green'
		labels[15] = 'Blue'
		labels[19] = 'Gray'
		c.triangle(
			spectra=spectra,
			markers=markers,
			colors=colors,
			text=labels,
			legend=True,
			gamut=True,
			gamutcolor="0.7",
			gamutedge=''
			)
		
		# Next we try to assess whether they're distinguishable. We don't
		# check the direction here.

		# R-Y
		print("R-Y")
		ry = c.color_disc(red25_all, yellow15_all)
		
		# R-G
		print("R-G")
		rg = c.color_disc(red25_all, green58_all)
		
		# R-B
		print("R-B")
		rb = c.color_disc(red25_all, blue47_all)
		
		# Y-G
		print("Y-G")
		yg = c.color_disc(yellow15_all, green58_all)
		
		# Y-B
		print("Y-B")
		yb = c.color_disc(yellow15_all, blue47_all)
		
		# G-B
		print("G-B")
		gb = c.color_disc(green58_all, blue47_all)

		# collect box plot data	
		labels = ["R-Y", "R-G", "R-B", "Y-G", "Y-B", "G-B"]
		boxes = [ry[3], rg[3], rb[3], yg[3], yb[3], gb[3]]
		
		# medians
		print("R-Y median contrast: " + str(ry[0]))
		print("R-G median contrast: " + str(rg[0]))
		print("R-B median contrast: " + str(rb[0]))
		print("Y-G median contrast: " + str(yg[0]))
		print("Y-B median contrast: " + str(yb[0]))
		print("G-B median contrast: " + str(gb[0]))
		print("Lowest median contrast: " + str(min(ry[0], rg[0], rb[0], yg[0], yb[0], gb[0])))
		print("Highest median contrast: " + str(max(ry[0], rg[0], rb[0], yg[0], yb[0], gb[0])))
		print("")
		
		# plot contrast values
		plt.boxplot(x=boxes, tick_labels=labels)
		plt.ylabel("ΔS (JND)")
		plt.show()
		
		# and again for colors vs gray. Is gray distinguishable from yellow
		# and green? If not, we have a major problem.
		
		# red
		print("Red vs. gray")
		r_gray = c.color_disc(red25_all, gray_all, correct=68, trials=80)
		
		# yellow
		print("Yellow vs. gray")
		y_gray = c.color_disc(yellow15_all, gray_all, correct=68, trials=80)
		
		# green
		print("Green vs. gray")
		g_gray = c.color_disc(green58_all, gray_all, correct=68, trials=80)
		
		# blue
		print("Blue vs. gray")
		b_gray = c.color_disc(blue47_all, gray_all, correct=68, trials=80)
		
		
		labels = ["red", "yellow", "green", "blue"]
		boxes = [r_gray[3], y_gray[3], g_gray[3], b_gray[3]]
		
		plt.boxplot(x=boxes, tick_labels=labels)
		plt.ylabel("ΔS (JND)")
		plt.show()
		# Yes, it is! This came as a bit of a surprise considering the rods and cones apparently
		# have very similar responses to these.

# variables used in tests
w_range = r1 - r2 + 1 # inclusive
xvalues = np.empty(w_range)
for i in range(w_range): xvalues[i] = i + r2

"""
Test 1: probability of success
I previously tested varying the relative contribution of rods and cones, but
other papers just measure rod and cone contrast separately. This is simpler
and addresses the issue of whether 493 nm is accurate for living rods or
they're red-shifted closer to other marsupials.

This is really slow. brightness_disc and color_disc potentially have O(n^2)
running time because they have a nested for loop, but in practice it's always
exactly 16 passes because there are 8 items to compare. (If there were 2, it
would be 1 pass; 4, 4; 6, 9; 8, 16; 10, 25; etc. n^2/2) Meanwhile, the code
below is O(n). This means each round of the for loop is doing 16 comparisons *
10 color pairs = 160. Is there any way we can not do this? I can't think of
one. Instead I recommend saving computing power by setting -s to 480 nm
because values below this are implausible for both L cone and rod pigments (in
terrestrial vertebrates).

Considerably better now that the custom luminosity is set once as an array in
brightness_disc() rather than doing all the calculations every time we run
brightness_contrast(). The latest improvement took it from ~67 s to ~24 s.
"""
if (args.wopt1):
	rod_ry = np.empty(w_range)
	rod_rg = np.empty(w_range)
	rod_rb = np.empty(w_range)
	rod_yg = np.empty(w_range)
	rod_yb = np.empty(w_range)
	rod_gb = np.empty(w_range)
	rod_rgray = np.empty(w_range)
	rod_ygray = np.empty(w_range)
	rod_ggray = np.empty(w_range)
	rod_bgray = np.empty(w_range)
	plus_ry = 0
	plus_rg = 0
	plus_rb = 0
	plus_yg = 0
	plus_yb = 0
	plus_gb = 0
	plus_rgray = 0
	plus_ygray = 0
	plus_ggray = 0
	plus_bgray = 0
	minus_ry = 0
	minus_rg = 0
	minus_rb = 0
	minus_yg = 0
	minus_yb = 0
	minus_gb = 0
	minus_rgray = 0
	minus_ygray = 0
	minus_ggray = 0
	minus_bgray = 0

	for i in range(w_range):
		w = r2 + i
		if (args.verbose): print(w)
		luminosity = np.empty(c.range1)
		for j in range(c.range1): luminosity[j] = c.vpt(j+300, w) * c.media_1nm[j]
		compare = [-1,-1,-1,-1]
		if (args.grayswap): compare = [
			math.sqrt(sum(green58_0*luminosity*c.wp_1nm) * sum(yellow15_0*luminosity*c.wp_1nm)),
			math.sqrt(sum(green58_03*luminosity*c.wp_1nm) * sum(yellow15_03*luminosity*c.wp_1nm)),
			math.sqrt(sum(green58_07*luminosity*c.wp_1nm) * sum(yellow15_07*luminosity*c.wp_1nm)),
			math.sqrt(sum(green58_10*luminosity*c.wp_1nm) * sum(yellow15_10*luminosity*c.wp_1nm))
			]
		rod_ry[i] = c.brightness_disc(red25_all, yellow15_all, r1=luminosity, output=False)[0]
		rod_rg[i] = c.brightness_disc(red25_all, green58_all, r1=luminosity, output=False)[0]
		rod_rb[i] = c.brightness_disc(red25_all, blue47_all, r1=luminosity, output=False)[0]
		rod_yg[i] = c.brightness_disc(yellow15_all, green58_all, r1=luminosity, output=False)[0]
		rod_yb[i] = c.brightness_disc(yellow15_all, blue47_all, r1=luminosity, output=False)[0]
		rod_gb[i] = c.brightness_disc(green58_all, blue47_all, r1=luminosity, output=False)[0]
		rod_rgray[i] = c.brightness_disc(red25_all, gray_all, correct=68, trials=80, r1=luminosity, output=False, compare=compare)[0]
		rod_ggray[i] = c.brightness_disc(green58_all, gray_all, correct=68, trials=80, r1=luminosity, output=False, compare=compare)[0]
		rod_ygray[i] = c.brightness_disc(yellow15_all, gray_all, correct=68, trials=80, r1=luminosity, output=False, compare=compare)[0]
		rod_bgray[i] = c.brightness_disc(blue47_all, gray_all, correct=68, trials=80, r1=luminosity, output=False, compare=compare)[0]
		if (args.verbose):
			print("RY: " + str(rod_ry[i]))
			print("RG: " + str(rod_rg[i]))
			print("RB: " + str(rod_rb[i]))
			print("YG: " + str(rod_yg[i]))
			print("YB: " + str(rod_yb[i]))
			print("GB: " + str(rod_gb[i]))
			print("R-gray: " + str(rod_rgray[i]))
			print("Y-gray: " + str(rod_ygray[i]))
			print("G-gray: " + str(rod_ggray[i]))
			print("B-gray: " + str(rod_bgray[i]))
		
		# threshold crossings
		if ((rod_ry[i] < 0.05) and (rod_ry[i-1] > 0.05)): minus_ry = w
		if ((rod_rg[i] < 0.05) and (rod_rg[i-1] > 0.05)): minus_rg = w
		if ((rod_rb[i] < 0.05) and (rod_rb[i-1] > 0.05)): minus_rb = w
		if ((rod_yg[i] < 0.05) and (rod_yg[i-1] > 0.05)): minus_yg = w
		if ((rod_yb[i] < 0.05) and (rod_yb[i-1] > 0.05)): minus_yb = w
		if ((rod_gb[i] < 0.05) and (rod_gb[i-1] > 0.05)): minus_gb = w
		if ((rod_rgray[i] < 0.05) and (rod_rgray[i-1] > 0.05)): minus_rgray = w
		if ((rod_ygray[i] < 0.05) and (rod_ygray[i-1] > 0.05)): minus_ygray = w
		if ((rod_ggray[i] < 0.05) and (rod_ggray[i-1] > 0.05)): minus_ggray = w
		if ((rod_bgray[i] < 0.05) and (rod_bgray[i-1] > 0.05)): minus_bgray = w
		if ((rod_ry[i] > 0.05) and (rod_ry[i-1] < 0.05)): plus_ry = w
		if ((rod_rg[i] > 0.05) and (rod_rg[i-1] < 0.05)): plus_rg = w
		if ((rod_rb[i] > 0.05) and (rod_rb[i-1] < 0.05)): plus_rb = w
		if ((rod_yg[i] > 0.05) and (rod_yg[i-1] < 0.05)): plus_yg = w
		if ((rod_yb[i] > 0.05) and (rod_yb[i-1] < 0.05)): plus_yb = w
		if ((rod_gb[i] > 0.05) and (rod_gb[i-1] < 0.05)): plus_gb = w
		if ((rod_rgray[i] > 0.05) and (rod_rgray[i-1] < 0.05)): plus_rgray = w
		if ((rod_ygray[i] > 0.05) and (rod_ygray[i-1] < 0.05)): plus_ygray = w
		if ((rod_ggray[i] > 0.05) and (rod_ggray[i-1] < 0.05)): plus_ggray = w
		if ((rod_bgray[i] > 0.05) and (rod_bgray[i-1] < 0.05)): plus_bgray = w
		
	print("Threshold crossings (- => +):")
	print("R-Y: " + str(plus_ry))
	print("R-G: " + str(plus_rg))
	print("R-B: " + str(plus_rb))
	print("Y-G: " + str(plus_yg))
	print("Y-B: " + str(plus_yb))
	print("G-B: " + str(plus_gb))
	print("R-gray: " + str(plus_rgray))
	print("Y-gray: " + str(plus_ygray))
	print("G-gray: " + str(plus_ggray))
	print("B-gray: " + str(plus_bgray))
	print("Threshold crossings (+ => -):")
	print("R-Y: " + str(minus_ry))
	print("R-G: " + str(minus_rg))
	print("R-B: " + str(minus_rb))
	print("Y-G: " + str(minus_yg))
	print("Y-B: " + str(minus_yb))
	print("G-B: " + str(minus_gb))
	print("R-gray: " + str(minus_rgray))
	print("Y-gray: " + str(minus_ygray))
	print("G-gray: " + str(minus_ggray))
	print("B-gray: " + str(minus_bgray))

	# print execution time
	print("%s seconds" % (time.time() - start_time))
	
	plt.plot(xvalues, rod_ry, 'k', label="RY")
	plt.plot(xvalues, rod_rg, '--k', label="RG")
	plt.plot(xvalues, rod_rb, '-.k', label="RB")
	plt.plot(xvalues, rod_yg, 'k', dashes=[1,5], label="YG")
	plt.plot(xvalues, rod_yb, 'k', dashes=[8,1], label="YB")
	plt.plot(xvalues, rod_gb, 'k', dashes=[8,1,1,1,1,1], label="GB")
	plt.plot([r2, r1], [0.05, 0.05], ':', color="gray")
	plt.xlabel("λmax (nm)")
	plt.ylabel("Probability of success")
	plt.legend()
	plt.show()
	
	plt.plot(xvalues, rod_rgray, color=cr0, label="R-gray")
	plt.plot(xvalues, rod_ygray, '--', color=cy0,label="Y-gray")
	plt.plot(xvalues, rod_ggray, '-.', color=cg0, label="G-gray")
	plt.plot(xvalues, rod_bgray, color=cb0, dashes=[8,1], label="B-gray")
	plt.plot([r2, r1], [0.05, 0.05], ':', color="gray")
	plt.xlabel("λmax (nm)")
	plt.ylabel("Probability of success")
	plt.legend()
	plt.show()

"""
Test 2: order

The opossums are described as ordering the colors as red-yellow-gray-green-blue. This doesn't tell us
whether they used color, but it does tell us (a) they used visual cues, (b) there were probably no more
than two pigments involved (dichromatic species prefer relative discrimination), (c) the contrast must
be such that they appear in this order.

This should be quicker than 1. (It is.)
"""
if (args.wopt2):
	# 2a: brightness contrast with 1 pigment
	print("brightness order")

	r = c.brightness_vary(red25_0)
	y = c.brightness_vary(yellow15_0)
	g = c.brightness_vary(green58_0)
	b = c.brightness_vary(blue47_0)
	if (args.grayswap):
		gray = np.empty(w_range)
		for i in range(w_range): gray[i] = math.sqrt(y[i] * g[i])
	else: gray = c.brightness_vary(gray_0)

	# order changes
	ry = 0
	yg = 0
	ygray = 0
	grayg = 0
	gb = 0
	r2 = args.receptors[1]
	for i in range(w_range):
		if (r[i] < y[i] and r[i-1] > y[i-1]): ry = i + r2
		if (y[i] < g[i] and y[i-1] > g[i-1]): yg = i + r2
		if (y[i] < gray[i] and y[i-1] > gray[i-1]): ygray = i + r2
		if (gray[i] < g[i] and gray[i-1] > g[i-1]): grayg = i + r2
		if (g[i] < b[i] and g[i-1] > b[i-1]): gb = i + r2

	print("Order changes (+ => -):")
	print("R-Y: " + str(ry))
	print("Y-G: " + str(yg))
	print("G-B: " + str(gb))
	print("Y-gray: " + str(ygray))
	print("Gray-G: " + str(grayg))
	
	plt.plot(xvalues, r, color=cr0, label="Red")
	plt.plot(xvalues, y, '--', color=cy0, label="Yellow")
	plt.plot(xvalues, g, '-.', color=cg0, label="Green")
	plt.plot(xvalues, b, color=cb0, dashes=[8,1], label="Blue")
	plt.plot(xvalues, gray, color=cgray0, dashes=[8,1,1,1,1,1], label="Gray")
	plt.xlabel("λmax (nm)")
	plt.ylabel("Brightness (arbitrary units)")
	plt.yscale("log")
	plt.legend()
	plt.show()

	# 2b: color contrast with 2 pigments
	# This is similar to the above but slightly different because it takes both into account.
	print("color order")

	r = c.color_vary(red25_0)
	y = c.color_vary(yellow15_0)
	g = c.color_vary(green58_0)
	b = c.color_vary(blue47_0)
	gray = c.color_vary(gray_0)

	# order changes
	ry = 0
	yg = 0
	ygray = 0
	grayg = 0
	gb = 0
	r2 = args.receptors[1]
	for i in range(w_range):
		if (r[i] > y[i] and r[i-1] < y[i-1]): ry = i + r2
		if (y[i] > g[i] and y[i-1] < g[i-1]): yg = i + r2
		if (y[i] > gray[i] and y[i-1] < gray[i-1]): ygray = i + r2
		if (gray[i] > g[i] and gray[i-1] < g[i-1]): grayg = i + r2
		if (g[i] > b[i] and g[i-1] < b[i-1]): gb = i + r2

	print("Order changes (- => +):")
	print("R-Y: " + str(ry))
	print("Y-G: " + str(yg))
	print("G-B: " + str(gb))
	print("Y-gray: " + str(ygray))
	print("Gray-G: " + str(grayg))
	    
	plt.plot(xvalues, r, color=cr0, label="Red")
	plt.plot(xvalues, y, '--', color=cy0, label="Yellow")
	plt.plot(xvalues, g, '-.', color=cg0, label="Green")
	plt.plot(xvalues, b, color=cb0, dashes=[8,1], label="Blue")
	plt.plot(xvalues, gray, color=cgray0, dashes=[8,1,1,1,1,1], label="Gray")
	plt.xlabel("λmax (nm)")
	plt.ylabel("Chromaticity")
	plt.legend()
	plt.show()
	
	# 2c
	print("brightness order")

	r = c.brightness_vary(red25_0)
	r3 = c.brightness_vary(red25_03)
	r7 = c.brightness_vary(red25_07)
	r10 = c.brightness_vary(red25_10)
	y = c.brightness_vary(yellow15_0)
	y3 = c.brightness_vary(yellow15_03)
	y7 = c.brightness_vary(yellow15_07)
	y10 = c.brightness_vary(yellow15_10)
	g = c.brightness_vary(green58_0)
	g3 = c.brightness_vary(green58_03)
	g7 = c.brightness_vary(green58_07)
	g10 = c.brightness_vary(green58_10)
	b = c.brightness_vary(blue47_0)
	b3 = c.brightness_vary(blue47_03)
	b7 = c.brightness_vary(blue47_07)
	b10 = c.brightness_vary(blue47_10)
	if (args.grayswap):
		gray = np.empty(w_range)
		gray3 = np.empty(w_range)
		gray7 = np.empty(w_range)
		gray10 = np.empty(w_range)
		for i in range(w_range):
			gray[i] = math.sqrt(y[i] * g[i])
			gray3[i] = math.sqrt(y3[i] * g3[i])
			gray7[i] = math.sqrt(y7[i] * g7[i])
			gray10[i] = math.sqrt(y10[i] * g10[i])
	else:
		gray = c.brightness_vary(gray_0)
		gray3 = c.brightness_vary(gray_03)
		gray7 = c.brightness_vary(gray_07)
		gray10 = c.brightness_vary(gray_10)
	
	plt.plot(xvalues, r, color=cr0, label="Red")
	plt.plot(xvalues, r3, color=cr03)
	plt.plot(xvalues, r7, color=cr07)
	plt.plot(xvalues, r10, color=cr10)
	plt.plot(xvalues, y, '--', color=cy0, label="Yellow")
	plt.plot(xvalues, y3, '--', color=cy03)
	plt.plot(xvalues, y7, '--', color=cy07)
	plt.plot(xvalues, y10, '--', color=cy10)
	plt.plot(xvalues, g, '-.', color=cg0, label="Green")
	plt.plot(xvalues, g3, '-.', color=cg03)
	plt.plot(xvalues, g7, '-.', color=cg07)
	plt.plot(xvalues, g10, '-.', color=cg10)
	plt.plot(xvalues, b, color=cb0, dashes=[8,1], label="Blue")
	plt.plot(xvalues, b3, color=cb03, dashes=[8,1])
	plt.plot(xvalues, b7, color=cb07, dashes=[8,1])
	plt.plot(xvalues, b10, color=cb10, dashes=[8,1])
	plt.plot(xvalues, gray, color=cgray0, dashes=[8,1,1,1,1,1], label="Gray")
	plt.plot(xvalues, gray3, color=cgray03, dashes=[8,1,1,1,1,1])
	plt.plot(xvalues, gray7, color=cgray07, dashes=[8,1,1,1,1,1])
	plt.plot(xvalues, gray10, color=cgray10, dashes=[8,1,1,1,1,1])
	plt.xlabel("λmax (nm)")
	plt.ylabel("Brightness (arbitrary units)")
	plt.yscale("log")
	plt.legend()
	plt.show()

	# 2d
	print("color order")

	r = c.color_vary(red25_0)
	r3 = c.color_vary(red25_03)
	r7 = c.color_vary(red25_07)
	r10 = c.color_vary(red25_10)
	y = c.color_vary(yellow15_0)
	y3 = c.color_vary(yellow15_03)
	y7 = c.color_vary(yellow15_07)
	y10 = c.color_vary(yellow15_10)
	g = c.color_vary(green58_0)
	g3 = c.color_vary(green58_03)
	g7 = c.color_vary(green58_07)
	g10 = c.color_vary(green58_10)
	b = c.color_vary(blue47_0)
	b3 = c.color_vary(blue47_03)
	b7 = c.color_vary(blue47_07)
	b10 = c.color_vary(blue47_10)
	gray = c.color_vary(gray_0)
	gray3 = c.color_vary(gray_03)
	gray7 = c.color_vary(gray_07)
	gray10 = c.color_vary(gray_10)
    
	for i in range(w_range):
		if (min(r[i],r3[i],r7[i],r10[i]) > max(y[i],y3[i],y7[i],y10[i])
			and min(r[i-1],r3[i-1],r7[i-1],r10[i-1]) < max(y[i-1],y3[i-1],y7[i-1],y10[i-1])):
			ry = i + r2
		if (min(y[i],y3[i],y7[i],y10[i]) > max(g[i],g3[i],g7[i],g10[i])
			and min(y[i-1],y3[i-1],y7[i-1],y10[i-1]) < max(g[i-1],g3[i-1],g7[i-1],g10[i-1])):
			yg = i + r2
		if (min(y[i],y3[i],y7[i],y10[i]) > max(gray[i],gray3[i],gray7[i],gray10[i])
			and min(y[i-1],y3[i-1],y7[i-1],y10[i-1]) < max(gray[i-1],gray3[i-1],gray7[i-1],gray10[i-1])):
			ygray = i + r2
		if (min(gray[i],gray3[i],gray7[i],gray10[i]) > max(g[i],g3[i],g7[i],g10[i])
			and min(gray[i-1],gray3[i-1],gray7[i-1],gray10[i-1]) < max(g[i-1],g3[i-1],g7[i-1],g10[i-1])):
			grayg = i + r2
		if (min(g[i],g3[i],g7[i],g10[i]) > max(b[i],b3[i],b7[i],b10[i])
			and min(g[i-1],g3[i-1],g7[i-1],g10[i-1]) < max(b[i-1],b3[i-1],b7[i-1],b10[i-1])):
			gb = i + r2

	print("Order changes (- => +):")
	print("R-Y: " + str(ry))
	print("Y-G: " + str(yg))
	print("G-B: " + str(gb))
	print("Y-gray: " + str(ygray))
	print("Gray-G: " + str(grayg))

	plt.plot(xvalues, r10, color=cr10)
	plt.plot(xvalues, r7, color=cr07)
	plt.plot(xvalues, r3, color=cr03)
	plt.plot(xvalues, r, color=cr0, label="Red")
	plt.plot(xvalues, y10, '--', color=cy10)
	plt.plot(xvalues, y7, '--', color=cy07)
	plt.plot(xvalues, y3, '--', color=cy03)
	plt.plot(xvalues, y, '--', color=cy0, label="Yellow")
	plt.plot(xvalues, g10, '-.', color=cg10)
	plt.plot(xvalues, g7, '-.', color=cg07)
	plt.plot(xvalues, g3, '-.', color=cg03)
	plt.plot(xvalues, g, '-.', color=cg0, label="Green")
	plt.plot(xvalues, b10, color=cb10, dashes=[8,1])
	plt.plot(xvalues, b7, color=cb07, dashes=[8,1])
	plt.plot(xvalues, b3, color=cb03, dashes=[8,1])
	plt.plot(xvalues, b, color=cb0, dashes=[8,1], label="Blue")
	plt.plot(xvalues, gray10, color=cgray10, dashes=[8,1,1,1,1,1])
	plt.plot(xvalues, gray7, color=cgray07, dashes=[8,1,1,1,1,1])
	plt.plot(xvalues, gray3, color=cgray03, dashes=[8,1,1,1,1,1])
	plt.plot(xvalues, gray, color=cgray0, dashes=[8,1,1,1,1,1], label="Gray")
	plt.xlabel("λmax (nm)")
	plt.ylabel("Chromaticity")
	plt.legend()
	plt.show()

"""
Test 3: discrimination

Here we test different values of a second pigment. This is set up so there are
only two because we assume S/UV cones are only negligibly involved.

As expected, this one runs slower than a dog in January.
(07/10: 776.850195646286 s) Not much we can do about that unless I use only
the brightest variant instead of all four.

I shaved it down to 165.02264022827148. It still isn't fast, but it's more
reasonable.

With the quantum noise calculations also delegated to arrays, the running time
improves further: 131.77082228660583 seconds. Python uses 30-40% of the
CPU while running it.
"""
if (args.wopt3):
	print("RY")
	ry = c.color_opt(red25_all, yellow15_all)
	print("RG")
	rg = c.color_opt(red25_all, green58_all)
	print("RB")
	rb = c.color_opt(red25_all, blue47_all)
	print("YG")
	yg = c.color_opt(yellow15_all, green58_all)
	print("YB")
	yb = c.color_opt(yellow15_all, blue47_all)
	print("GB")
	gb = c.color_opt(green58_all, blue47_all)
	print("R-gray")
	rgray = c.color_opt(red25_all, gray_all)
	print("Y-gray")
	ygray = c.color_opt(yellow15_all, gray_all)
	print("G-gray")
	ggray = c.color_opt(green58_all, gray_all)
	print("B-gray")
	bgray = c.color_opt(blue47_all, gray_all)

	# print execution time
	print("%s seconds" % (time.time() - start_time))

	lowest_col = 0
	lowest_coly = 0
	lowest_gray = 0
	lowest_grayy = 0
	lowest = 0
	lowest_y = 0
	for i in range(w_range):
		cur_x = i + args.receptors[1]
		cur_col = min(ry[0][i],rg[0][i],rb[0][i],yg[0][i],yb[0][i],gb[0][i])
		if (cur_col > lowest_coly):
			lowest_coly = cur_col
			lowest_col = cur_x
		cur_gray = min(rgray[0][i],ggray[0][i],ygray[0][i],bgray[0][i])
		if (cur_gray > lowest_grayy):
			lowest_grayy = cur_gray
			lowest_gray = cur_x
		cur = min(cur_col, cur_gray)
		if (cur > lowest_y):
			lowest_y = cur
			lowest = cur_x
	print("Best min (color pairs): " + str(lowest_col))
	print("Best min (gray): " + str(lowest_gray))
	print("Best min (all): " + str(lowest))

	plt.subplot(2,1,1)
	plt.plot(xvalues, ry[0], color=cr0, label="R-Y")
	plt.plot(xvalues, rg[0], '--', color=cr0, label="R-G")
	plt.plot(xvalues, rb[0], '-.', color=cr0, label="R-B")
	plt.plot(xvalues, yg[0], ':', color=cy0, label="Y-G")
	plt.plot(xvalues, yb[0], color=cy0, dashes=[8,1], label="Y-B")
	plt.plot(xvalues, gb[0], color=cg0, dashes=[8,1,1,1,1,1], label="G-B")
	plt.axvline(color='k', x=lowest_col)
	plt.axvline(color='k', linestyle=":", x=lowest)
	#plt.xlabel("λmax (nm)")
	plt.ylabel("Contrast")
	plt.legend()
	plt.subplot(2,1,2)
	plt.plot(xvalues, rgray[0], color=cr0, label="R-gray")
	plt.plot(xvalues, ygray[0], '--', color=cy0, label="Y-gray")
	plt.plot(xvalues, ggray[0], '-.', color=cg0, label="G-gray")
	plt.plot(xvalues, bgray[0], color=cb0, dashes=[8,1], label="B-gray")
	plt.axvline(color='k', x=lowest_gray)
	plt.axvline(color='k', linestyle=":", x=lowest)
	plt.xlabel("λmax (nm)")
	plt.ylabel("Contrast")
	plt.legend()
	plt.show()

# darkest versions
if (args.wopt4):
	print("RY")
	ry = c.color_opt([red25_10], [yellow15_10])
	print("RG")
	rg = c.color_opt([red25_10], [green58_10])
	print("RB")
	rb = c.color_opt([red25_10], [blue47_10])
	print("YG")
	yg = c.color_opt([yellow15_10], [green58_10])
	print("YB")
	yb = c.color_opt([yellow15_10], [blue47_10])
	print("GB")
	gb = c.color_opt([green58_10], [blue47_10])
	print("R-gray")
	rgray = c.color_opt([red25_10], [gray_10])
	print("Y-gray")
	ygray = c.color_opt([yellow15_10], [gray_10])
	print("G-gray")
	ggray = c.color_opt([green58_10], [gray_10])
	print("B-gray")
	bgray = c.color_opt([blue47_10], [gray_10])

	# print execution time
	print("%s seconds" % (time.time() - start_time))

	lowest_col = 0
	lowest_coly = 0
	lowest_gray = 0
	lowest_grayy = 0
	lowest = 0
	lowest_y = 0
	for i in range(w_range):
		cur_x = i + args.receptors[1]
		cur_col = min(ry[0][i],rg[0][i],rb[0][i],yg[0][i],yb[0][i],gb[0][i])
		if (cur_col > lowest_coly):
			lowest_coly = cur_col
			lowest_col = cur_x
		cur_gray = min(rgray[0][i],ggray[0][i],ygray[0][i],bgray[0][i])
		if (cur_gray > lowest_grayy):
			lowest_grayy = cur_gray
			lowest_gray = cur_x
		cur = min(cur_col, cur_gray)
		if (cur > lowest_y):
			lowest_y = cur
			lowest = cur_x
	print("Best min (color pairs): " + str(lowest_col))
	print("Best min (gray): " + str(lowest_gray))
	print("Best min (all): " + str(lowest))

	plt.subplot(2,1,1)
	plt.plot(xvalues, ry[0], color=cr0, label="R-Y")
	plt.plot(xvalues, rg[0], '--', color=cr0, label="R-G")
	plt.plot(xvalues, rb[0], '-.', color=cr0, label="R-B")
	plt.plot(xvalues, yg[0], ':', color=cy0, label="Y-G")
	plt.plot(xvalues, yb[0], color=cy0, dashes=[8,1], label="Y-B")
	plt.plot(xvalues, gb[0], color=cg0, dashes=[8,1,1,1,1,1], label="G-B")
	plt.axvline(color='k', x=lowest_col)
	plt.axvline(color='k', linestyle=":", x=lowest)
	#plt.xlabel("λmax (nm)") runs into the bottom plot
	plt.ylabel("Contrast")
	plt.legend()
	plt.subplot(2,1,2)
	plt.plot(xvalues, rgray[0], color=cr0, label="R-gray")
	plt.plot(xvalues, ygray[0], '--', color=cy0, label="Y-gray")
	plt.plot(xvalues, ggray[0], '-.', color=cg0, label="G-gray")
	plt.plot(xvalues, bgray[0], color=cb0, dashes=[8,1], label="B-gray")
	plt.axvline(color='k', x=lowest_gray)
	plt.axvline(color='k', linestyle=":", x=lowest)
	plt.xlabel("λmax (nm)")
	plt.ylabel("Contrast")
	plt.legend()
	plt.show()

