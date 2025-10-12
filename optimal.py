"""
Three types of optimal colors: see https://en.wikipedia.org/wiki/Gamut#Optimal_colors

* long-pass step functions (0 up to a certain wavelength and then 1)
* short-pass step functions (1 then 0)
* monochromatic single wavelengths (set to 100 because otherwise they come out really dark)
"""

import central as c
import numpy as np
import matplotlib.pyplot as plt

# "long-pass" step functions: 0 for wavelengths less than the step and 1 for longer wavelengths.
# Many white, yellow, orange and red flowers have reflectance spectra like these, and they're
# also found in paints and camera filters.
print("Long-pass step functions")
longpass = []
# 350-nm step
step350 = np.ones(401)
for i in range(49): step350[i] = 0
longpass.append(step350)

# 400-nm step
step400 = np.ones(401)
for i in range(99): step400[i] = 0
longpass.append(step400)

# 450-nm step
step450 = np.ones(401)
for i in range(149): step450[i] = 0
longpass.append(step450)

# 500-nm step
step500 = np.ones(401)
for i in range(199): step500[i] = 0
longpass.append(step500)

# 550-nm step
step550 = np.ones(401)
for i in range(249): step550[i] = 0
longpass.append(step550)

# 600-nm step
step600 = np.ones(401)
for i in range(299): step600[i] = 0
longpass.append(step600)

# 650-nm step
step650 = np.ones(401)
for i in range(349): step650[i] = 0
longpass.append(step650)

# plots
colors = []
text = []
for i in range(len(longpass)):
	colors.append(c.spec2rgb(longpass[i]))
	text.append(i*50 + 350)

for i in range(len(longpass)):
	plt.plot(c.x_1nm, longpass[i], color=colors[i])
plt.show()

c.triangle(
	spectra=longpass,
	colors=colors,
	gamut=True,
	gamutcolor="0.7",
	gamutedge='',
	text=text,
	legend=True
	)

# "short-pass" step functions, the reverse of the above. These form the other edge of the
# optimal color solid (see Wikipedia) but don't seem to occur in nature or have any
# applications. Blue/green surfaces usually are the "band-pass" type with one or two
# peaks.
print("Short-pass step functions")
shortpass = []
# 350-nm step
step350 = np.zeros(401)
for i in range(49): step350[i] = 1
shortpass.append(step350)

# 400-nm step
step400 = np.zeros(401)
for i in range(99): step400[i] = 1
shortpass.append(step400)

# 450-nm step
step450 = np.zeros(401)
for i in range(149): step450[i] = 1
shortpass.append(step450)

# 500-nm step
step500 = np.zeros(401)
for i in range(199): step500[i] = 1
shortpass.append(step500)

# 550-nm step
step550 = np.zeros(401)
for i in range(249): step550[i] = 1
shortpass.append(step550)

# 600-nm step
step600 = np.zeros(401)
for i in range(299): step600[i] = 1
shortpass.append(step600)

# 650-nm step
step650 = np.zeros(401)
for i in range(349): step650[i] = 1
shortpass.append(step650)

# plots
colors = []
text = []
for i in range(len(shortpass)):
	colors.append(c.spec2rgb(shortpass[i]))
	text.append(i*50 + 350)

for i in range(len(shortpass)):
	plt.plot(c.x_1nm, shortpass[i], color=colors[i])
plt.show()

c.triangle(
	spectra=shortpass,
	colors=colors,
	gamut=True,
	gamutcolor="0.7",
	gamutedge='',
	text=text,
	legend=True
	)

# single wavelengths
print("Single wavelengths")

mono = []

# 300
wl300 = np.zeros(401)
wl300[0] = 100
mono.append(wl300)

mono = []
# 310
wl310 = np.zeros(401)
wl310[10] = 100
mono.append(wl310)

# 320
wl320 = np.zeros(401)
wl320[20] = 100
mono.append(wl320)

# 330
wl330 = np.zeros(401)
wl330[30] = 100
mono.append(wl330)

# 340
wl340 = np.zeros(401)
wl340[40] = 100
mono.append(wl340)

# 350
wl350 = np.zeros(401)
wl350[50] = 100
mono.append(wl350)

# 360
wl360 = np.zeros(401)
wl360[60] = 100
mono.append(wl360)

# 370
wl370 = np.zeros(401)
wl370[70] = 100
mono.append(wl370)

# 380
wl380 = np.zeros(401)
wl380[80] = 100
mono.append(wl380)

# 390
wl390 = np.zeros(401)
wl390[90] = 100
mono.append(wl390)

# 400
wl400 = np.zeros(401)
wl400[100] = 100
mono.append(wl400)

# 410
wl410 = np.zeros(401)
wl410[110] = 100
mono.append(wl410)

# 420
wl420 = np.zeros(401)
wl420[120] = 100
mono.append(wl420)

# 430
wl430 = np.zeros(401)
wl430[130] = 100
mono.append(wl430)

# 440
wl440 = np.zeros(401)
wl440[140] = 100
mono.append(wl440)

# 450
wl450 = np.zeros(401)
wl450[150] = 100
mono.append(wl450)

# 460
wl460 = np.zeros(401)
wl460[160] = 100
mono.append(wl460)

# 470
wl470 = np.zeros(401)
wl470[170] = 100
mono.append(wl470)

# 480
wl480 = np.zeros(401)
wl480[180] = 100
mono.append(wl480)

# 490
wl490 = np.zeros(401)
wl490[190] = 100
mono.append(wl490)

# 500
wl500 = np.zeros(401)
wl500[200] = 100
mono.append(wl500)

# 510
wl510 = np.zeros(401)
wl510[210] = 100
mono.append(wl510)

# 520
wl520 = np.zeros(401)
wl520[220] = 100
mono.append(wl520)

# 530
wl530 = np.zeros(401)
wl530[230] = 100
mono.append(wl530)

# 540
wl540 = np.zeros(401)
wl540[240] = 100
mono.append(wl540)

# 550
wl550 = np.zeros(401)
wl550[250] = 100
mono.append(wl550)

# 560
wl560 = np.zeros(401)
wl560[260] = 100
mono.append(wl560)

# 570
wl570 = np.zeros(401)
wl570[270] = 100
mono.append(wl570)

# 580
wl580 = np.zeros(401)
wl580[280] = 100
mono.append(wl580)

# 590
wl590 = np.zeros(401)
wl590[290] = 100
mono.append(wl590)

# 600
wl600 = np.zeros(401)
wl600[300] = 100
mono.append(wl600)

# 610
wl610 = np.zeros(401)
wl610[310] = 100
mono.append(wl610)

# 620
wl620 = np.zeros(401)
wl620[320] = 100
mono.append(wl620)

# 630
wl630 = np.zeros(401)
wl630[330] = 100
mono.append(wl630)

# 640
wl640 = np.zeros(401)
wl640[340] = 100
mono.append(wl640)

# 650
wl650 = np.zeros(401)
wl650[350] = 100
mono.append(wl650)

# 660
wl660 = np.zeros(401)
wl660[360] = 100
mono.append(wl660)

# 670
wl670 = np.zeros(401)
wl670[370] = 100
mono.append(wl670)

# 680
wl680 = np.zeros(401)
wl680[380] = 100
mono.append(wl680)

# 690
wl690 = np.zeros(401)
wl690[390] = 100
mono.append(wl690)

# 700
wl700 = np.zeros(401)
wl700[400] = 100
mono.append(wl700)

# plots
colors = []
text = []
for i in range(len(mono)):
	colors.append(c.spec2rgb(mono[i]))
	text.append(i*10 + 300)

for i in range(len(mono)):
	plt.plot(c.x_1nm, mono[i], color=colors[i])
plt.show()

c.triangle(
	spectra=mono,
	colors=colors,
	#gamut=True,
	#gamutcolor="0.7",
	#gamutedge='',
	#text=text,
	#legend=True
	)
