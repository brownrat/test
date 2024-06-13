# extra data related to or derived from colorvision.py

import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import scipy

parser = argparse.ArgumentParser()
parser.add_argument("--performance", help="performance graphs", action="store_true")
parser.add_argument("--erg", help="ERG graphs", action="store_true")
args = parser.parse_args()

# copied from colorvision.py
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

# curve fitting

# fit 1 visual pigment template
def vpt_fit1(xdata, scale, shift, t1):
	ydata = np.empty(xdata.shape[0])
	for i in range(xdata.shape[0]):
		ydata[i] = scale*math.log(vpt(xdata[i], t1)) + shift
	return(ydata)

# fit 2 (SciPy doesn't seem to want to figure out what number of templates provides the
# best fit, so I just have to pick one)
def vpt_fit2(xdata, scale, shift, t1, t2, scalet1, scalet2):
	# bounds
	#if (t1 < 530):
	#	t1 = 530
	#if (t1 > 570):
	#	t1 = 570
	#if (t2 < 350):
	#	t2 = 350
	#if (t2 > 515):
	#	t2 = 515
	scalet1 = abs(scalet1)
	scalet2 = abs(scalet2)
	if (scalet2 < 0):
		scalet2 = 0
	ydata = np.empty(xdata.shape[0])
	for i in range(xdata.shape[0]):
		value = scalet1*vpt(xdata[i], t1) + scalet2*vpt(xdata[i], t2)
		if (value >= 0):
			ydata[i] = scale*math.log(value) + shift
		else:
			ydata[i] = 0
	return(ydata)

# fit 1-3
def vpt_fit3(xdata, scale, shift, t1, t2, t3, scalet1, scalet2, scalet3):
	# bounds
	if (t1 < 555):
		t1 = 555
	if (t1 > 562):
		t1 = 562
	if (t2 < 500):
		t2 = 500
	if (t2 > 510):
		t2 = 510
	if (t3 < 350):
		t3 = 350
	if (t3 > 365):
		t3 = 365
	scalet1 = abs(scalet1)
	scalet2 = abs(scalet2)
	scalet3 = abs(scalet3)
	ydata = np.empty(xdata.shape[0])
	for i in range(xdata.shape[0]):
		value = scalet1*vpt(xdata[i], t1) + scalet2*vpt(xdata[i], t2) + scalet3*vpt(xdata[i], t3)
		if (value >= 0):
			ydata[i] = scale*math.log(value) + shift
		else:
			ydata[i] = 0
	return(ydata)

def template_filter(w, a, b, c, d, e, f):
	density = a*math.exp((b - w) / c) + d*math.exp((e - w) / f)
	if (density < 0):
		return 0
	return math.exp(-density)

def filter_fit(xdata, a, b, c, d, e, f):
	ydata = np.empty(xdata.shape[0])
	for i in range(xdata.shape[0]):
		ydata[i] = template_filter(xdata[i], a, b, c, d, e, f)
	return(ydata)

# Data for Thylamys elegans from Palacios et al. (2010)
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

if (args.performance):
	# upper limit of performance with varying levels of dark adaptation
	x = np.array([
		0.0,
		0.05,
		0.1,
		0.15,
		0.2,
		0.25,
		0.3,
		0.35,
		0.4,
		0.45,
		0.5,
		0.55,
		0.6,
		0.65,
		0.7,
		0.75,
		0.8,
		0.85,
		0.9,
		0.95,
		1.0
	])

	y = np.array([
		62.5, # R-Y, Y-G, G-B
		62.5, # R-Y, Y-G, Y-B, G-B
		62.5, # all except R-G
		62.5, # all equal
		62.5, # all equal
		62.5, # all equal
		62.5, # all except G-B
		62.5, # R-Y, R-G, Y-G
		62.5, # R-Y, R-G, Y-G
		62.5, # R-Y, R-G, Y-G
		62.5, # R-Y, R-G, R-B, Y-G
		62.5, # R-Y
		62.5, # R-Y
		62.5, # R-Y
		62.5, # R-Y
		68.75, # R-Y
		75, # R-Y
		81.25, # R-Y, Y-G
		81.25, # R-Y
		87.5, # R-Y
		93.75 # R-Y, Y-G, G-B
	])

	y1 = np.array([
		37.5, # G-B
		37.5, # G-B
		37.5, # Y-G, Y-B, G-B
		37.5, # R-B, Y-G, Y-B, G-B
		37.5, # R-G, Y-G, Y-B, G-B
		37.5, # R-G, Y-G, Y-B
		37.5, # Y-G
		37.5, # Y-G
		37.5, # Y-G
		37.5, # Y-G
		37.5, # Y-G
		43.75, # R-Y
		43.75, # R-Y
		37.5, # R-Y
		56.25, # R-Y
		56.25, # R-Y
		56.25, # R-Y
		62.5, # R-Y
		68.75, # R-Y
		75.0, # R-Y
		87.5 # R-Y
	])

	y2 = np.empty(21)
	for i in range(21):
		y2[i] = 87.5

	# p-values -- lowest
	pvalues = np.array([
		6.91306013322901e-07, # G-B
		6.91306013322901e-07, # G-B
		6.91306013322901e-07, # Y-G, Y-B, G-B
		6.91306013322901e-07, # R-B, Y-G, Y-B, G-B
		6.91306013322901e-07, # R-G, Y-G, Y-B, G-B
		6.91306013322901e-07, # R-G, Y-G, Y-B
		6.91306013322901e-07, # Y-G
		6.91306013322901e-07, # Y-G
		6.91306013322901e-07, # Y-G
		6.91306013322901e-07, # Y-G
		6.91306013322901e-07, # Y-G
		4.263225520892315e-06, # R-Y
		4.263225520892315e-06, # R-Y
		6.91306013322901e-07, # R-Y
		0.00010750156966200634, # R-Y
		0.00044783204074632933, # R-Y
		0.001661061926555766, # R-Y
		0.016325647807321812, # R-Y
		0.04327398257201562, # R-Y
		0.21328292249639952, # R-Y
		0.8321331599980754 # R-Y
	])

	# p-values -- highest
	pvalues1 = np.array([
		0.001661061926555766, # R-G, R-B
		0.00044783204074632933, # R-G, R-B
		0.00044783204074632933, # R-G
		0.00010750156966200634, # R-G
		4.263225520892315e-06, # R-Y, R-B
		0.00010750156966200634, # G-B
		0.00044783204074632933, # G-B
		0.005505225640202571, # R-B
		0.04327398257201562, # R-B
		0.04327398257201562, # R-B
		0.21328292249639952, # R-B
		0.3899284990417625, # R-B, Y-B
		0.8321331599980754, # R-B
		0.9633061393856144, # R-B
		1.0, # R-B
		1.0, # R-B
		1.0, # R-B, Y-B
		1.0, # R-B, Y-B
		1.0, # R-G, R-B, Y-B
		1.0, # R-G, R-B, Y-B
		1.0 # R-G, R-B, Y-B
	])

	y3 = np.empty(21)
	for i in range(21):
		y3[i] = 0.05

	plt.plot(x*100, y, 'o-k', label="upper limit of expected performance")
	#plt.plot(x*100, y1, 's-k', mfc='w', label="lower limit of expected performance")
	plt.plot(x*100, y2, ':k', label="lower limit of actual performance")
	plt.ylabel("% correct")
	plt.xlabel("Rod contribution to spectral sensitivity (%)")
	plt.legend()
	plt.show()

	plt.plot(x*100, pvalues, 'o-k', label="lowest P-value")
	#plt.plot(x*100, pvalues1, 's-k', mfc='w', label="highest P-value")
	plt.plot(x*100, y3, ':k', label="0.05 threshold")
	plt.ylabel("% correct")
	plt.xlabel("Rod contribution to spectral sensitivity (%)")
	plt.legend()
	plt.show()

	# performance limit with various types of S-cones
	x = np.array([
		360,
		365,
		370,
		375,
		380,
		385,
		390,
		395,
		400,
		405,
		410,
		415,
		420,
		425,
		430,
		435,
		440,
		445,
		450,
		455,
		460,
		465,
		470,
		475,
		480,
		485,
		490,
		495,
		500,
		505,
		510,
		515,
		520,
		525,
		530,
		535,
		540,
		545,
		550,
		555,
		560
		#565
	])

	# lowest % distinguishable
	y = np.array([
		0.375, #360, R-G
		0.3125, #370, R-G
		0.25, #380, R-G
		0.25, #390, R-G
		0.1875, #400, R-G
		0.375, #410, R-G
		0.375, #420, R-G
		0.8125, #430, R-Y, G-B
		0.5, #440, R-Y
		0.4375, #450, R-Y
		0.4375, #460, R-Y
		1, #470, all
		0.8125, #480, R-Y
		0.8125, #490, R-Y
		0.8125, #500, R-G, Y-G, Y-B
		0.8125, #510, R-G, Y-G, Y-B
		0.625, #520, Y-G
		0.375, #530, G-B
		0.375, #540, G-B
		0.0625, #550, G-B
		0, #560, all
	])

	# lowest median contrast (w=0.05, lp=0.9, sp=0.1)
	y1 = np.array([
		0.3266314163271572, #360
		0.3123140968136736,
		0.294692481663287, #370
		0.2749841053971982,
		0.25358900702795295, #380
		0.2297576011301315,
		0.20101546212040633, #390
		0.16217117434685305,
		0.1385978582741992, #400
		0.22413398793797112,
		0.360140825448527, #410
		0.5290701835481114,
		0.8246340141211987, #420
		1.2044545505984863,
		1.2437681847173452, #430
		1.1628075198472771,
		1.0110039459276385, #440
		0.7811668198932582,
		0.47352629369445715, #450
		0.10597334789894955,
		0.45124712598197514, #460
		0.9370662009444986,
		1.4408740678430836, #470
		1.957620198861704,
		2.4692178598487162, #480
		2.9206392178065466,
		3.2840400911191026, #490
		3.404534843345984,
		2.9894431838138678, #500
		2.5680310349127566,
		2.1519652754687524, #510
		1.7670201880157301,
		1.4155645749064338, #520
		1.0997111775551862,
		0.8213182784333408, #530
		0.5815076585767733,
		0.38204817592999074, #540
		0.22508031347466514,
		0.10895973688862494, # 550
		0.03463295053117855,
		0.0017301516821745318 # 560
		#0.012404390197030556
	])

	# highest median contrast
	y2 = np.array([
		2.1621469641940374, #360
		2.4865347637322404,
		2.8332275461018823, #370
		3.2041900331533317,
		3.5998111446102925, #380
		4.02695746566881,
		4.498676795965365, #390
		4.941442562133805,
		5.377072697866122, #400
		5.84698102314721,
		6.315881356235241, #410
		6.777888967503436,
		7.227211583413264, #420
		7.658275829179262,
		8.065854580697074, #430
		8.358370551868125,
		8.704125234996223, #440
		9.013430252007048,
		9.282905654763674, #450
		9.509452518332456,
		9.6899987429787, #460
		9.821099808930427,
		9.898382365798803, #470
		9.915897607223247,
		9.865624483556992, #480
		9.737613107339428,
		9.52141349397467, #490
		9.209095476225826,
		8.79900787867561, #500
		8.297984561938536,
		7.7196473719410275, #510
		7.078991912572704,
		6.387063699478136, #520
		5.650353488917943,
		4.875368545505786, #530
		4.073547854820973,
		3.26185713958644, #540
		2.459406468283915,
		1.68246206224998, # 550
		0.945721970920892,
		0.265284992872467 # 560
	])

	#plt.plot(x, y*100, 'o-k')
	#plt.ylabel("Minimum distinguishable color/brightness pairs (%)")
	#plt.xlabel("λmax of S cone (nm)")
	#plt.show()

	# lowest median contrast for M cones given S=362 and L=562 (w=0.06, lp=0.62, mp=0.31, sp=0.07)
	y3 = np.array([
		0.4577810630009087, #360
		0.44776081789168903,
		0.4357331354341544, #370
		0.4227265819906353,
		0.40919044233894086, #380
		0.39489345912858065,
		0.3788622999185072, #390
		0.3596327561653114,
		0.3372274074867613, #400
		0.32125224453533163,
		0.39655260606543763, #410
		0.5996616557723204,
		0.9083700248880022, #420
		1.263912820556233,
		1.451097969967503, #430
		1.3824769148088656,
		1.2523591523899231, #440
		1.0632256495333963,
		0.8438788737748868, #450
		0.6893491353807499,
		0.8220135600947149, #460
		1.2432181350217457,
		1.8058928030480033, #470
		2.389418289871686,
		2.920852805745188, #480
		3.484594743888741,
		3.960957270273307, #490
		4.317496733216242,
		4.534780539468157, #500
		4.176826303773539,
		3.8055488298720848, #510
		3.5127580000497227,
		3.2020868450037656, #520
		2.967152190750146,
		2.7944747449506826, #530
		2.429652156671667,
		2.055459685638413, #540
		1.6386336705431992,
		1.249677500317743, # 550
		0.8262518446018268,
		0.508541829016269 # 560
	])

	plt.plot(x, y1, 'o-k')
	#plt.plot(x, y2, 's-k', mfc='w', label="highest median contrast")
	plt.ylabel("Lowest median contrast (JND)")
	plt.xlabel("λmax of S cone (nm)")
	plt.plot([362, 362], [0, 4], ':k')
	plt.plot([493, 493], [0, 4], ':k')
	plt.text(363, 3.5, 'D. aurita S cone')
	plt.text(494, 3.5, 'D. virginiana rod')
	plt.show()

	plt.plot(x, y3, 'o-k')
	plt.ylabel("Lowest median contrast (JND)")
	plt.xlabel("λmax of M cone (nm)")
	plt.plot([362, 362], [0, 5], ':k')
	plt.plot([493, 493], [0, 5], ':k')
	plt.text(363, 4.75, 'D. aurita S cone')
	plt.text(494, 4.75, 'D. virginiana rod')
	plt.show()

	x1 = np.array([
		490,
		491,
		492,
		493,
		494,
		495,
		496,
		497,
		498,
		499,
		500,
		501,
		502,
		503,
		504,
		505,
		506,
		507,
		508,
		509,
		510
	])

	# S cones
	y4 = np.array([
		3.2840400911191026, #490
		3.3441619681478674,
		3.3997664789917073,
		3.450742959490288,
		3.496998101636657,
		3.404534843345984, # 495
		3.3089980109980814,
		3.2223337300604356,
		3.1533525970922387,
		3.0760519243187847,
		2.9894431838138678, #500
		2.9215088715949395,
		2.831513041100504,
		2.7425852550251006,
		2.6547500073427583,
		2.5680310349127566, # 505
		2.482451379458902,
		2.398033441796644,
		2.3147990289765064,
		2.232769394968497,
		2.1519652754687524 #510
	])

	# M cones
	y5 = np.array([
		3.960957270273307, #490
		4.042684666768339,
		4.119372624390497,
		4.190838916757469,
		4.256924920553885,
		4.317496733216242, # 495
		4.372446095712421,
		4.421691067385867,
		4.465176403960991,
		4.502873596627392,
		4.534780539468157, #500
		4.532630556593004,
		4.427249951609493,
		4.32396550077768,
		4.2498452529733335,
		4.176826303773539, # 505
		4.10310829512685,
		4.0075637129664665,
		3.9144899584563957,
		3.849336456485055,
		3.8055488298720848 #510
	])

	plt.plot(x1, y4, 'o-k')
	plt.plot(x1, y5, 's-k', mfc='w')
	plt.ylabel("Lowest median contrast (JND)")
	plt.xlabel("λmax of M/S cone (nm)")
	plt.plot([493, 493], [2, 5], ':k')
	plt.text(493.5, 4.75, 'D. virginiana rod')
	plt.show()

# ERG graphs (Jacobs & Williams 2010)

if (args.erg):
	# average for several animals, 460-650 (figure 1)
	erg0 = np.array([
		-0.764865291, # 460
		-0.65598016, # 470
		-0.538437606, # 480
		-0.403219479, # 490
		-0.311450795, # 500
		-0.216994702, # 510
		-0.153506933, # 520
		-0.102355991, # 530
		-0.070846579, # 540
		-0.036036524, # 550
		0, # 560
		-0.011074287, # 570
		-0.051024687, # 580
		-0.090145418, # 590
		-0.128905422, # 600
		-0.219339421, # 610
		-0.344727765, # 620
		-0.553966858, # 630
		-0.74387104, # 640
		-0.984349002 # 650
	])

	x0 = np.array([
		460,
		470,
		480,
		490,
		500,
		510,
		520,
		530,
		540,
		550,
		560,
		570,
		580,
		590,
		600,
		610,
		620,
		630,
		640,
		650
	])

	# This is not a real best-fit curve, just an approximation of what the original image shows
	x0_detail = np.empty(190)
	#best_fit0 = np.empty(190)
	for i in range(190):
		x0_detail[i] = i+460
	#	best_fit0[i] = 0.44*math.log(vpt(i+460, 562.4)) - 0.01

	#plt.plot(x0, erg0, 'ok')
	#plt.plot(x0_detail, best_fit0, 'k')
	#plt.show()
	
	# actual best fit curve
	best_fit = scipy.optimize.curve_fit(vpt_fit1, x0, erg0, p0=[1, 0, 560])
	print(best_fit[0])
	best_fit0 = vpt_fit1(x0_detail, *best_fit[0])
	plt.plot(x0, erg0, 'ok')
	plt.plot(x0_detail, best_fit0, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# one animal, 370-670 (figure 2)
	erg1 = np.array([
		-0.676985896, # 370
		-0.687185684, # 380
		-0.69775213, # 390
		-0.776750484, # 400
		-0.804949897, # 410
		-0.850848941, # 420
		-0.86864857, # 430
		-0.926680694, # 440
		-0.918614196, # 450
		-0.777917127, # 460
		-0.724351576, # 470
		-0.577321306, # 480
		-0.372758901, # 490
		-0.309926877, # 500
		-0.211162267, # 510
		-0.11886419, # 520
		-0.09409804, # 530
		-0.06303202, # 540
		-0.034332618, # 550
		0, # 560
		-0.009233141, # 570
		-0.078398367, # 580
		-0.152363492, # 590
		-0.204395742, # 600
		-0.287694006, # 610
		-0.364825733, # 620
		-0.5392221, # 630
		-0.729518135, # 640
		-1.026578613, # 650
		-1.261140393, # 660
		-1.5664007 # 670
	])
	
	x1 = np.empty(31)
	for i in range(31):
		x1[i] = i*10 + 370
	
	x1_detail = np.empty(300)
	for i in range(300):
		x1_detail[i] = i+370
	best_fit = scipy.optimize.curve_fit(vpt_fit1, x1, erg1, p0=[1, 0, 560])
	print(best_fit[0])
	best_fit1 = vpt_fit1(x1_detail, *best_fit[0])
	plt.plot(x1, erg1, 'ok')
	plt.plot(x1_detail, best_fit1, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# one animal, 390-670 nm, sensitivity differences (figure 3)
	erg2 = np.array([
		0.009655977, # 390
		-0.019895044, #400
		-0.03419242, #410
		-0.020058309, #420
		-0.059731778, #430
		-0.020571429, #440
		0.111673469, #450
		0.141271137, #460
		0.009166181, #470
		0.080769679, #480
		0.13187172, #490
		0.156921283, #500
		0.06071137, #510
		0.049772595, #520
		0.071626822, #530
		0.081002915, #540
		0.008419825, #550
		0.029690962, #560
		-0.009632653, #570
		-0.001399417, #580
		-0.019988338, #590
		0.029854227, #600
		-0.059335277, #610
		0.020011662, #620
		0.025889213, #630
		0.029527697, #640
		0.111953353, #650
		-0.010542274, #660
		-0.009212828 #670
	])
	
	x2 = np.empty(29)
	for i in range(29):
		x2[i] = i*10 + 390
	
	plt.plot(x2, erg2, 'ok')
	plt.show()
	
	# correction for lens filtering
	#print(opossum_fit[0])
	filter_x = np.empty(101)
	for i in range(101):
		filter_x[i] = i*4 + 300
	filter_x1 = np.empty(400)
	for i in range(400):
		filter_x1[i] = i + 300
	plt.plot(filter_x, opossum_filter_data, 'ok')
	filter_y = filter_fit(filter_x1, *opossum_fit[0])
	plt.plot(filter_x1, filter_y, 'k')
	#plt.show()
	
	#print(mouse_fit[0])
	filter_x2 = np.empty(40)
	for i in range(40):
		filter_x2[i] = i*10 + 310
	filter_x3 = np.empty(400)
	for i in range(400):
		filter_x3[i] = i + 300
	plt.plot(filter_x2, mouse_filter_data/100, 'or')
	filter_y1 = filter_fit(filter_x3, *mouse_fit[0])
	plt.plot(filter_x3, filter_y1, 'r')
	plt.show()
	
	erg0_adjusted = np.empty(20)
	peak = 0
	for i in range(20):
		#erg0_adjusted[i] = math.exp(erg0[i]) / template_filter(i*10 + 450, *opossum_fit[0])
		erg0_adjusted[i] = math.exp(erg0[i]) / mouse_filter_data[i+14]**2
		peak = max(peak, erg0_adjusted[i])
	# adjust downward and convert back to log format
	for i in range(20):
		erg0_adjusted[i] = math.log(erg0_adjusted[i] / peak)
	plt.plot(x0, erg0_adjusted, 'ok')
	plt.plot(x0_detail, best_fit0, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	erg1_adjusted = np.empty(31)
	peak = 0
	for i in range(31):
		#erg1_adjusted[i] = math.exp(erg1[i]) / template_filter(i*10 + 370, *opossum_fit[0])
		erg1_adjusted[i] = math.exp(erg1[i]) / mouse_filter_data[i+6]**2
		peak = max(peak, erg1_adjusted[i])
	# adjust downward and convert back to log format
	for i in range(31):
		erg1_adjusted[i] = math.log(erg1_adjusted[i] / peak)
	plt.plot(x1, erg1_adjusted, 'ok')
	plt.plot(x1_detail, best_fit1, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# more curve fitting
	best_fit = scipy.optimize.curve_fit(vpt_fit1, x0, erg0_adjusted, p0=[1, 0, 560])
	best_fit2 = vpt_fit1(x0_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x0, erg0_adjusted, 'ok')
	plt.plot(x0_detail, best_fit2, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	best_fit = scipy.optimize.curve_fit(vpt_fit1, x1, erg1_adjusted, p0=[1, 0, 560])
	best_fit3 = vpt_fit1(x1_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x1, erg1_adjusted, 'ok')
	plt.plot(x1_detail, best_fit3, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# 2 templates
	best_fit = scipy.optimize.curve_fit(vpt_fit2, x0, erg0_adjusted, p0=[1, 0, 560, 500, 1, 1])
	best_fit2 = vpt_fit2(x0_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x0, erg0_adjusted, 'ok')
	plt.plot(x0_detail, best_fit2, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	#best_fit = scipy.optimize.curve_fit(vpt_fit2, x1, erg1_adjusted, p0=[1, 0, 560, 500, 1, 1])
	#best_fit3 = vpt_fit2(x1_detail, *best_fit[0])
	#print(best_fit[0])
	#plt.plot(x1, erg1_adjusted, 'ok')
	#plt.plot(x1_detail, best_fit3, 'k')
	#plt.show()
	
	#best_fit = scipy.optimize.curve_fit(vpt_fit2, x0, erg0_adjusted, p0=[1, 0, 560, 360, 1, 1])
	#best_fit2 = vpt_fit2(x0_detail, *best_fit[0])
	#print(best_fit[0])
	#plt.plot(x0, erg0_adjusted, 'ok')
	#plt.plot(x0_detail, best_fit2, 'k')
	#plt.show()
	
	best_fit = scipy.optimize.curve_fit(vpt_fit2, x1, erg1_adjusted, p0=[1, 0, 560, 360, 1, 1])
	best_fit3 = vpt_fit2(x1_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x1, erg1_adjusted, 'ok')
	plt.plot(x1_detail, best_fit3, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# 3 templates
	#best_fit = scipy.optimize.curve_fit(vpt_fit3, x0, erg0_adjusted, p0=[1, 0, 560, 500, 360, 1, 1, 1])
	#best_fit2 = vpt_fit3(x0_detail, *best_fit[0])
	#print(best_fit[0])
	#plt.plot(x0, erg0_adjusted, 'ok')
	#plt.plot(x0_detail, best_fit2, 'k')
	#plt.show()
	
	best_fit = scipy.optimize.curve_fit(vpt_fit3, x1, erg1_adjusted, p0=[1, 0, 560, 500, 360, 1, 1, 1])
	best_fit3 = vpt_fit3(x1_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x1, erg1_adjusted, 'ok')
	plt.plot(x1_detail, best_fit3, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# chromatic adaptation
	erg0_adapted = np.empty(20) # 460-650
	peak = 0
	for i in range(20):
		erg0_adapted[i] = erg0[i] + erg2[i+7]
		peak = max(peak, erg0_adapted[i])
	# adjust downward
	for i in range(20):
		erg0_adapted[i] = erg0_adapted[i] - peak
	plt.plot(x0, erg0_adapted, 'ok')
	plt.plot(x0_detail, best_fit0, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	erg1_adapted = np.empty(29) # 390-670
	peak = 0
	for i in range(29):
		erg1_adapted[i] = erg1[i+2] + erg2[i]
		peak = max(peak, erg1_adapted[i])
	# adjust downward
	for i in range(29):
		erg1_adapted[i] = erg1_adapted[i] - peak
	plt.plot(x2, erg1_adapted, 'ok')
	plt.plot(x1_detail, best_fit1, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	erg0_adapted_adjusted = np.empty(20)
	peak = 0
	for i in range(20):
		erg0_adapted_adjusted[i] = erg0_adjusted[i] + erg2[i+7]
		peak = max(peak, erg0_adapted_adjusted[i])
	# adjust downward
	for i in range(20):
		erg0_adapted_adjusted[i] = erg0_adapted_adjusted[i] - peak
	plt.plot(x0, erg0_adapted_adjusted, 'ok')
	plt.plot(x0_detail, best_fit0, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	erg1_adapted_adjusted = np.empty(29)
	peak = 0
	for i in range(29):
		erg1_adapted_adjusted[i] = erg1_adjusted[i+2] + erg2[i]
		peak = max(peak, erg1_adapted_adjusted[i])
	# adjust downward
	for i in range(29):
		erg1_adapted_adjusted[i] = erg1_adapted_adjusted[i] - peak
	plt.plot(x2, erg1_adapted_adjusted, 'ok')
	plt.plot(x1_detail, best_fit1, 'k')
	plt.show()
	
	# further curve fitting
	best_fit = scipy.optimize.curve_fit(vpt_fit1, x0, erg0_adapted_adjusted, p0=[1, 0, 560])
	best_fit2 = vpt_fit1(x0_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x0, erg0_adjusted, 'ok')
	plt.plot(x0_detail, best_fit2, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	best_fit = scipy.optimize.curve_fit(vpt_fit1, x2, erg1_adapted_adjusted, p0=[1, 0, 560])
	best_fit3 = vpt_fit1(x1_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x2, erg1_adapted_adjusted, 'ok')
	plt.plot(x1_detail, best_fit3, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# 2 templates
	best_fit = scipy.optimize.curve_fit(vpt_fit2, x0, erg0_adapted_adjusted, p0=[1, 0, 560, 500, 1, 1])
	best_fit2 = vpt_fit2(x0_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x0, erg0_adjusted, 'ok')
	plt.plot(x0_detail, best_fit2, 'k')
	plt.show()
	
	#best_fit = scipy.optimize.curve_fit(vpt_fit2, x2, erg1_adapted_adjusted, p0=[1, 0, 560, 500, 1, 1])
	#best_fit3 = vpt_fit2(x1_detail, *best_fit[0])
	#print(best_fit[0])
	#plt.plot(x2, erg1_adapted_adjusted, 'ok')
	#plt.plot(x1_detail, best_fit3, 'k')
	#plt.show()
	
	#best_fit = scipy.optimize.curve_fit(vpt_fit2, x0, erg0_adapted_adjusted, p0=[1, 0, 560, 360, 1, 1])
	#best_fit2 = vpt_fit2(x0_detail, *best_fit[0])
	#print(best_fit[0])
	#plt.plot(x0, erg0_adjusted, 'ok')
	#plt.plot(x0_detail, best_fit2, 'k')
	#plt.show()
	
	best_fit = scipy.optimize.curve_fit(vpt_fit2, x2, erg1_adapted_adjusted, p0=[1, 0, 560, 360, 1, 1])
	best_fit3 = vpt_fit2(x1_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x2, erg1_adapted_adjusted, 'ok')
	plt.plot(x1_detail, best_fit3, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
	
	# 3 templates
	#best_fit = scipy.optimize.curve_fit(vpt_fit3, x0, erg0_adapted_adjusted, p0=[1, 0, 560, 500, 360, 1, 1, 1])
	#best_fit2 = vpt_fit3(x0_detail, *best_fit[0])
	#print(best_fit[0])
	#plt.plot(x0, erg0_adjusted, 'ok')
	#plt.plot(x0_detail, best_fit2, 'k')
	#plt.show()
	
	best_fit = scipy.optimize.curve_fit(vpt_fit3, x2, erg1_adapted_adjusted, p0=[1, 0, 560, 500, 360, 1, 0.5, 0.1])
	best_fit3 = vpt_fit3(x1_detail, *best_fit[0])
	print(best_fit[0])
	plt.plot(x2, erg1_adapted_adjusted, 'ok')
	plt.plot(x1_detail, best_fit3, 'k')
	plt.ylabel("Log relative sensitivity")
	plt.xlabel("Wavelength (nm)")
	plt.show()
