"""
This file analyzes the ERG data from Jacobs & Williams (2009) to test their conclusion that they recorded signals
from only one visual pigment with peak sensitivity of 562 nm. First we show how the graphs could differ if
ocular media transmission were known, and then we reconstruct the effects of chromatic adaptation by
combining their fig. 3 with fig. 1 or fig. 2. Which species' ocular media we use is specified by the --media
argument.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import scipy
import central as c
args = c.args

# simplify printing and interpreting curve_fit output
def opt_nm(popt, pcov, i=0):
    return str(round(popt[i], 1)) + " nm ± " + str(round(np.sqrt(np.diag(pcov))[i], 1))
def opt_percent(popt, pcov, i=0):
    try: return str(round(popt[i]/total * 100)) + "% ± " + str(round(np.sqrt(np.diag(pcov))[i]/total * 100))
    except OverflowError: return "Could not calculate percentage: one or more parameters is infinite"
def opt_raw():
        if (args.verbose):
                print("Raw output:")
                print("popt: " + str(popt))
                print("pcov: " + str(pcov))
                print("condition number: " + str(np.linalg.cond(pcov)))
                print("diagonal elements: " + str(np.diag(pcov)))
                print("standard deviations: " + str(np.sqrt(np.diag(pcov))))

# figure 1
fig1 = np.array([
        [460.123,-0.766813],
        [469.975,-0.657009],
        [480.248,-0.540837],
        [490.108,-0.403942],
        [500.164,-0.311671],
        [510.012,-0.217803],
        [520.270,-0.154221],
        [530.315,-0.103385],
        [540.355,-0.0716718],
        [550.187,-0.0367684],
        [560.438,-0.000278095],
        [570.257,-0.0115901],
        [580.277,-0.0515907],
        [590.297,-0.0931850],
        [600.318,-0.131592],
        [610.324,-0.222589],
        [620.322,-0.347052],
        [630.086,-0.557568],
        [640.066,-0.747370],
        [650.031,-0.989762],
])

for i in range(20): plt.plot(round(fig1[i][0]), fig1[i][1], 'o')
plt.show()

# figure 2
fig2 = np.array([
        [370.193,-0.677121],
        [380.275,-0.687269],
        [390.078,-0.699926],
        [400.152,-0.779149],
        [410.232,-0.806880],
        [420.171,-0.853447],
        [430.112,-0.871130],
        [440.189,-0.930258],
        [450.134,-0.921566],
        [460.236,-0.779752],
        [470.186,-0.728361],
        [480.009,-0.579009],
        [490.118,-0.378168],
        [500.350,-0.315476],
        [510.446,-0.216362],
        [520.262,-0.124781],
        [530.209,-0.0997632],
        [540.296,-0.0709791],
        [550.104,-0.0421922],
        [560.333,-0.00587418],
        [570.135,-0.0135080],
        [580.351,-0.0826846],
        [590.145,-0.158136],
        [600.223,-0.209729],
        [610.156,-0.295229],
        [620.371,-0.371941],
        [630.153,-0.547864],
        [640.073,-0.736347],
        [650.120,-1.03284],
        [660.035,-1.26528],
        [670.081,-1.56805],
])

for i in range(31): plt.plot(round(fig2[i][0]), fig2[i][1], 'o')
plt.show()

# figure 3
fig3 = np.array([
        [390.000,0.00892070],
        [400.123,-0.0216754],
        [410.253,-0.0354776],
        [420.018,-0.0217969],
        [430.139,-0.0615532],
        [439.913,-0.0219184],
        [450.096,0.107790],
        [460.243,0.138262],
        [469.954,0.00843243],
        [480.115,0.0785991],
        [490.269,0.128919],
        [500.038,0.153286],
        [510.138,0.0585685],
        [520.269,0.0478197],
        [530.036,0.0676073],
        [540.175,0.0782324],
        [549.909,0.00794416],
        [560.052,0.0292562],
        [569.797,-0.0120246],
        [580.311,-0.00140178],
        [590.063,-0.0228353],
        [600.217,0.0274842],
        [609.944,-0.0611245],
        [620.108,0.0166757],
        [630.245,0.0227207],
        [640.006,0.0257145],
        [650.170,0.105041],
        [659.886,-0.0141014],
        [670.022,-0.0126366],
])

for i in range(29): plt.plot(round(fig3[i][0]), fig3[i][1], 'o')
plt.show()

# plots
for i in range(20): plt.plot(round(fig1[i][0]), 10**fig1[i][1], 'ok')
for i in range(31): plt.plot(round(fig2[i][0]), 10**fig2[i][1], 'or')
for i in range(29): plt.plot(round(fig3[i][0]), 10**fig3[i][1], 'ob')
plt.show()
# fig. 3 combined with the others for context
for i in range(29): plt.plot(round(fig3[i][0]), 10**(fig3[i][1] + fig2[i+2][1]), 'ok')
for i in range(20): plt.plot(round(fig3[i+7][0]), 10**(fig3[i+7][1] + fig1[i][1]), 'or')
plt.show()

x = np.empty(401)
for i in range(401): x[i] = i+300

# conversion to linear sensitivity
# figure 1
x_fig1 = np.empty(20)
# we round x because it's a dot plot that only really has integer x values
for i in range(20): x_fig1[i] = round(fig1[i][0])
y_fig1 = np.empty(20)
for i in range(20): y_fig1[i] = 10**fig1[i][1]

# figure 1 + 3
y_fig13 = np.empty(20)
for i in range(20): y_fig13[i] = 10**(fig1[i][1] + fig3[i+7][1])

# figure 2
x_fig2 = np.empty(31)
for i in range(31): x_fig2[i] = round(fig2[i][0])
y_fig2 = np.empty(31)
for i in range(31): y_fig2[i] = 10**fig2[i][1]

# figure 2 + 3
x_fig23 = np.empty(29)
for i in range(29): x_fig23[i] = round(fig2[i+2][0])
y_fig23 = np.empty(29)
for i in range(29): y_fig23[i] = 10**(fig2[i+2][1] + fig3[i][1])

# choose ocular media
if (args.filter == "mouse"): media_10nm = c.mouse_10nm
elif (args.filter == "thylamys"): media_10nm = c.thylamys_10nm
elif (args.filter == "brushtail"): media_10nm = c.brushtail_10nm

# remove unnecessary parameters from fitting functions -- we "needed" these because of the
# log->linear issue. You should (I think) never need to shift the whole thing up and
# down the y-axis. That just doesn't make sense.
def vpt_fit1(xdata, t1, scalet1):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1)
                if (value >= 0): ydata[i] = value
                else: ydata[i] = 0
        return(ydata)

def vpt_fit2(xdata, t1, t2, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1) + scalet2*c.vpt(xdata[i], t2)
                if (value >= 0): ydata[i] = value
                else: ydata[i] = 0
        return(ydata)

def vpt_fit3(xdata, t1, t2, t3, scalet1, scalet2, scalet3):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1) + scalet2*c.vpt(xdata[i], t2) + scalet3*c.vpt(xdata[i], t3)
                if (value >= 0):
                        ydata[i] = value
                else:
                        ydata[i] = 0
        return(ydata)

# also test the difference between unfixed and fixed templates
def vpt_fit2_fixs(xdata, t1, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1) + scalet2*c.vpt(xdata[i], args.sw)
                if (value >= 0): ydata[i] = value
                else: ydata[i] = 0
        return(ydata)

# In these two, the first parameter is named t2 not t1 so the relationship with the
# coefficients is intuitive.
def vpt_fit2_fixl(xdata, t2, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], args.lw) + scalet2*c.vpt(xdata[i], t2)
                if (value >= 0): ydata[i] = value
                else: ydata[i] = 0
        return(ydata)

def vpt_fit3_fixls(xdata, t2, scalet1, scalet2, scalet3):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], args.lw) + scalet2*c.vpt(xdata[i], t2) + scalet3*c.vpt(xdata[i], args.sw)
                if (value >= 0): ydata[i] = value
                else: ydata[i] = 0
        return(ydata)

def vpt_fit3_fixall(xdata, scalet1, scalet2, scalet3):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], args.lw) + scalet2*c.vpt(xdata[i], args.mw) + scalet3*c.vpt(xdata[i], args.sw)
                if (value >= 0): ydata[i] = value
                else: ydata[i] = 0
        return(ydata)

# fixed L and M
def vpt_fit2_fixall(xdata, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], args.lw) + scalet2*c.vpt(xdata[i], args.mw)
                if (value >= 0): ydata[i] = value
                else: ydata[i] = 0
        return(ydata)

# 1: As before, we try fitting one template to the original data to make sure our
# results aren't too far off. SciPy comes up with 563.6 nm (rounding off to one
# decimal place) for fig. 1 and 560.5 nm for fig. 2. Jacobs & Williams' fits are
# 562.4 and 561.6, so these may not be quite right but are pretty close.
print("Original data from fig. 1")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig1, y_fig1, p0=[560, 1])
opt_raw()
print("LWS: " + opt_nm(popt, pcov))
x_1nm = np.empty(401)
lws = np.empty(401)
for i in range(401):
        x_1nm[i] = i + 300
        lws[i] = c.vpt(i + 300, popt[0]) * popt[1]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig1, y_fig1, 'ok')
plt.plot(x_1nm, lws, 'k')
plt.show()

print("Original data from fig. 2")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig2, y_fig2, p0=[560, 1])
opt_raw()
print("LWS: " + opt_nm(popt, pcov))
x_1nm = np.empty(401)
lws = np.empty(401)
for i in range(401):
        x_1nm[i] = i + 300
        lws[i] = c.vpt(i + 300, popt[0]) * popt[1]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig2, y_fig2, 'ok')
plt.plot(x_1nm, lws, 'k')
plt.show()

# 2: custom ocular media
# fig. 2
print("Custom ocular media (fig. 2)")
media_fig2 = np.empty(31)
for i in range(31):
        media_fig2[i] = y_fig2[i] / media_10nm[i+7]
popt, pcov, infodict, mesg, ier = scipy.optimize.curve_fit(vpt_fit2, x_fig2, media_fig2, p0=[args.lw, 360, 1, 0.1], full_output=True)
opt_raw()
total = popt[2] + popt[3]
print("LWS: " + opt_nm(popt, pcov) + " (" + opt_percent(popt, pcov, 2) + ")")
print("SWS: " + opt_nm(popt, pcov, 1) + " (" + opt_percent(popt, pcov, 3) + ")")
curve = vpt_fit2(x_1nm, *popt)
lws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = c.vpt(i + 300, popt[0]) * popt[2]
        sws[i] = c.vpt(i + 300, popt[1]) * popt[3]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig2, media_fig2, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, sws, ':k')
plt.show()

# We're "expecting" to find we can fit a UV pigment in, but let's also compare just
# one template.
# result: about the same value (559.2 nm). The curve deviates in the UV, but this isn't
# really visible in the residual plot where the errors in this region are smaller than
# at longer wavelengths.
# Further reading: https://statisticsbyjim.com/regression/check-residual-plots-regression-analysis/
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig2, media_fig2, p0=[args.lw, 1])
opt_raw()
print("LWS: " + opt_nm(popt, pcov))
curve = vpt_fit1(x_1nm, *popt)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig2, media_fig2, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()

# two types of residual plots: vs. independent variable and vs. fitted values
# The second plot shows what may be some heteroskedasticity (cone shape),
# but probably not enough to worry about. That word is hard to spell.
# https://www.statology.org/heteroscedasticity-regression/
fitted = vpt_fit1(x_fig2, *popt)
residual = fitted - media_fig2
plt.plot(x_fig2, residual, 'o')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Residuals")
plt.show()
plt.plot(fitted, residual, 'o')
plt.xlabel("Fitted values")
plt.ylabel("Residuals")
plt.show()

# fig. 2 + 3
print("fig. 2 + 3")
media_fig23 = np.empty(29)
for i in range(29):
        media_fig23[i] = y_fig23[i] / media_10nm[i+9]
popt, pcov = scipy.optimize.curve_fit(vpt_fit2, x_fig23, media_fig23, p0=[args.lw, 360, 1, 0.1])
opt_raw()
total = popt[2] + popt[3]
print("LWS: " + opt_nm(popt, pcov) + " (" + opt_percent(popt, pcov, 2) + ")")
print("SWS: " + opt_nm(popt, pcov, 1) + " (" + opt_percent(popt, pcov, 3) + ")")
curve = vpt_fit2(x_1nm, *popt)
lws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = c.vpt(i + 300, popt[0]) * popt[2]
        sws[i] = c.vpt(i + 300, popt[1]) * popt[3]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, sws, ':k')
plt.show()

fitted = vpt_fit2(x_fig23, *popt)
residual = fitted - media_fig23
plt.plot(x_fig23, residual, 'o')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Residuals")
plt.show()

# fixed S
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixs, x_fig23, media_fig23, p0=[args.lw, 1, 0.1])
opt_raw()
total = popt[1] + popt[2]
print("LWS: " + opt_nm(popt, pcov) + " (" + opt_percent(popt, pcov, 1) + ")")
print("SWS: " + str(args.sw) + " nm (fixed) (" + opt_percent(popt, pcov, 2) + ")")
curve = vpt_fit2_fixs(x_1nm, *popt)
lws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = c.vpt(i + 300, popt[0]) * popt[1]
        sws[i] = c.vpt(i + 300, args.sw) * popt[2]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, sws, ':k')
plt.show()

# fixed L, second variable (S doesn't work well so we're just doing L and M)
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixl, x_fig23, media_fig23, p0=[args.mw, 1, 0.1])
opt_raw()
total = popt[1] + popt[2]
print("LWS: " + str(args.lw) + " nm (fixed) (" + opt_percent(popt, pcov, 1) + ")")
print("MWS: " + opt_nm(popt, pcov) + " (" + opt_percent(popt, pcov, 2) + ")")

curve = vpt_fit2_fixl(x_1nm, *popt)
lws = np.empty(401)
mws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = c.vpt(i + 300, args.lw) * popt[1]
        mws[i] = c.vpt(i + 300, popt[0]) * popt[2]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, mws, '-.k')
plt.show()

# fixed L and M
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixall, x_fig23, media_fig23, p0=[1, 0.1])
opt_raw()
total = popt[0] + popt[1]
print("LWS: " + str(args.lw) + " nm (fixed) (" + opt_percent(popt, pcov) + ")")
print("MWS: " + str(args.mw) + " nm (fixed) (" + opt_percent(popt, pcov, 1) + ")")
curve = vpt_fit2_fixall(x_1nm, *popt)
lws = np.empty(401)
mws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = c.vpt(i + 300, args.lw) * popt[0]
        mws[i] = c.vpt(i + 300, args.mw) * popt[1]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, mws, '-.k')
plt.show()

# fig. 1
print("fig. 1")
media_fig1 = np.empty(20)
for i in range(20):
        media_fig1[i] = y_fig1[i] / media_10nm[i+16]
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig1, media_fig1, p0=[args.lw, 1])
opt_raw()
print("LWS: " + opt_nm(popt, pcov))
curve = vpt_fit1(x_1nm, *popt)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig1, media_fig1, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()

# fig. 1+3
print("fig. 1 + 3")
media_fig13 = np.empty(20)
for i in range(20):
        media_fig13[i] = y_fig13[i] / media_10nm[i+16]
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig1, media_fig13, p0=[args.lw, 1])
opt_raw()
print("LWS: " + opt_nm(popt, pcov))
curve = vpt_fit1(x_1nm, *popt)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig1, media_fig13, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()

# 2 templates, L fixed
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixl, x_fig1, media_fig13, p0=[args.mw, 1, 0.1])
opt_raw()
total = popt[1] + popt[2]
print("LWS: " + str(args.lw) + " nm (fixed) (" + opt_percent(popt, pcov, 1) + ")")
print("MWS: " + opt_nm(popt, pcov) + " (" + opt_percent(popt, pcov, 2) + ")")
curve = vpt_fit2_fixl(x_1nm, *popt)
lws = np.empty(401)
mws = np.empty(401)
for i in range(401):
        lws[i] = c.vpt(i + 300, args.lw) * popt[1]
        mws[i] = c.vpt(i + 300, popt[0]) * popt[2]
plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative sensitivity")
plt.yscale("log")
plt.ylim(0.01)
plt.plot(x_fig1, media_fig13, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, mws, '-.k')
plt.show()
