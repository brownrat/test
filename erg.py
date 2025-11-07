"""
This file analyzes the ERG data from Jacobs & Williams (2010) to test their
conclusion that they recorded signals from only one visual pigment with peak
sensitivity of 562 nm. First we show how the graphs could differ if ocular
media transmission were known, and then we reconstruct the effects of chromatic
adaptation by combining their fig. 3 with fig. 1 or fig. 2. Which species'
ocular media we use is specified by the --media argument.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import scipy
import central as c
args = c.args
lw = args.receptors[0]
if (len(args.receptors) > 1): mw = args.receptors[1]
if (len(args.receptors) > 2): sw = args.receptors[2]

# simplify printing and interpreting curve_fit output
# I'm still not quite sure what the square root of the diagonal covariance
# elements represents. This person on Stack Overflow isn't either:
# https://stackoverflow.com/questions/44723194/what-exactly-is-the-variance-on-the-parameters-of-scipy-curve-fit-python
def opt(popt, pcov, i=0, r=1):
    return str(round(popt[i], r)) + " ± " + str(round(np.sqrt(np.diag(pcov))[i], r))
def opt_percent(start=1):
    total = sum(popt[start:])
    string = "Percentages: "
    for i in range(start, len(popt)):
        try:
            string += str(round(popt[i]/total * 100)) + "% "
        except OverflowError:
            print("Could not calculate percentage: "
                  + "one or more parameters is infinite")
    print(string)
def opt_raw():
        if (args.verbose):
                print("Raw output:")
                print("popt: " + str(popt))
                print("pcov: " + str(pcov))
                print("condition number: " + str(np.linalg.cond(pcov)))
                print("diagonal elements: " + str(np.diag(pcov)))
                print("standard deviations: " + str(np.sqrt(np.diag(pcov))))

# root mean square error
def rmse(data, fit, p, start=0, interval=1):
    s = 0
    for i in range(data.shape[0]): s += (data[i] - fit[start+i*interval])**2
    s /= (data.shape[0] - p)
    s = math.sqrt(s)
    print("Root mean square error: " + str(s))

# figure 1
fig1 = c.csv2spec('fig1-jw2010.csv', interp=False)

# figure 2
fig2 = c.csv2spec('fig2-jw2010.csv', interp=False)

# figure 3
fig3 = c.csv2spec('fig3-jw2010.csv', interp=False)

x = np.empty(401)
for i in range(401): x[i] = i+300

# linearize or not
zero = 0
if (args.log): zero = -np.inf

ylabel = "Relative sensitivity"
if (args.log): ylabel = "Log relative sensitivity"

def ylim():
    if (args.log): plt.ylim(-2)

def exp(value):
    if (args.log): return value
    return 10**value

def log(value):
    if (args.log):
        if (value < 0): return -np.inf
        return math.log(value, 10)
    return value

def divide(x, y):
    if (args.log): return x - y
    return x / y

# figure 1
x_fig1 = np.empty(20)
# we round x because it's a dot plot that only really has integer x values
for i in range(20): x_fig1[i] = round(fig1[0][i])
y_fig1 = np.empty(20)
for i in range(20): y_fig1[i] = exp(fig1[1][i])

# figure 2
x_fig2 = np.empty(31)
for i in range(31): x_fig2[i] = round(fig2[0][i])
y_fig2 = np.empty(31)
for i in range(31): y_fig2[i] = exp(fig2[1][i])

# figure 3
x_fig3 = np.empty(29)
for i in range(29): x_fig3[i] = round(fig3[0][i])
y_fig3 = np.empty(29)
for i in range(29): y_fig3[i] = exp(fig3[1][i])

# figure 1 + 3
y_fig13 = np.empty(20)
for i in range(20): y_fig13[i] = exp(fig1[1][i] + fig3[1][i+7])

# figure 2 + 3
x_fig23 = np.empty(29)
for i in range(29): x_fig23[i] = round(fig2[0][i+2])
y_fig23 = np.empty(29)
for i in range(29): y_fig23[i] = exp(fig2[1][i+2] + fig3[1][i])

# choose ocular media
media_10nm = c.media_10nm

# fitting functions
def vpt_fit1(xdata, t1, scalet1):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

def vpt_fit2(xdata, t1, t2, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1) + scalet2*c.vpt(xdata[i], t2)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

def vpt_fit3(xdata, t1, t2, t3, scalet1, scalet2, scalet3):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1) + scalet2*c.vpt(xdata[i], t2) + scalet3*c.vpt(xdata[i], t3)
                if (value >= 0):
                        ydata[i] = log(value)
                else:
                        ydata[i] = zero
        return(ydata)

# also test the difference between unfixed and fixed templates
def vpt_fit2_fixs(xdata, t1, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1) + scalet2*c.vpt(xdata[i], sw)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

# In these two, the first parameter is named t2 not t1 so the relationship with the
# coefficients is intuitive.
def vpt_fit2_fixl(xdata, t2, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], lw) + scalet2*c.vpt(xdata[i], t2)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

def vpt_fit3_fixls(xdata, t2, scalet1, scalet2, scalet3):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], lw) + scalet2*c.vpt(xdata[i], t2) + scalet3*c.vpt(xdata[i], sw)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

def vpt_fit3_fixall(xdata, scalet1, scalet2, scalet3):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], lw) + scalet2*c.vpt(xdata[i], mw) + scalet3*c.vpt(xdata[i], sw)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

# fixed L and M
def vpt_fit2_fixall(xdata, scalet1, scalet2):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], lw) + scalet2*c.vpt(xdata[i], mw)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

# self-screening
def vpt_fit1_od(xdata, t1, scalet1, od):
        ydata = np.empty(xdata.shape[0])
        for i in range(xdata.shape[0]):
                value = scalet1*c.vpt(xdata[i], t1, od=od)
                if (value >= 0): ydata[i] = log(value)
                else: ydata[i] = zero
        return(ydata)

"""
1: As before, we try fitting one template to the original data to make sure
our results aren't too far off. SciPy comes up with 563.6 nm (rounding off to
one decimal place) for fig. 1 and 560.5 nm for fig. 2. Jacobs & Williams' fits
are 562.4 and 561.6, so these may not be quite right but are pretty close.
"""
print("Original data from fig. 1")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig1, y_fig1, p0=[560, 1])
opt_raw()
print("LWS: " + opt(popt, pcov))
x_1nm = np.empty(401)
lws = np.empty(401)
for i in range(401):
        x_1nm[i] = i + 300
        lws[i] = log(c.vpt(i + 300, popt[0]) * popt[1])
rmse(y_fig1, lws, 2, start=160, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
##plt.yscale("log")
ylim()
plt.plot(x_fig1, y_fig1, 'ok')
plt.plot(x_1nm, lws, 'k')
plt.show()
print("")

print("Original data from fig. 2")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig2, y_fig2, p0=[560, 1])
opt_raw()
print("LWS: " + opt(popt, pcov))
x_1nm = np.empty(401)
lws = np.empty(401)
for i in range(401):
        x_1nm[i] = i + 300
        lws[i] = log(c.vpt(i + 300, popt[0]) * popt[1])
rmse(y_fig2, lws, 2, start=70, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
##plt.yscale("log")
ylim()
plt.plot(x_fig2, y_fig2, 'ok')
plt.plot(x_1nm, lws, 'k')
plt.show()
print("")

# linear regression on fig. 3
print("Original data from fig. 3")
line = scipy.stats.linregress(x_fig3, y_fig3)
print(f"R-squared: {line.rvalue**2:.4f}")
print("Slope: " + str(line.slope))
print(f"P-value: {line.pvalue:.4f}")
fit = line.intercept + line.slope*x_fig3
plt.plot(x_fig3, y_fig3, 'ok')
plt.plot(x_fig3, fit, 'k')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Log sensitivity difference")
plt.show()
# residual
plt.plot(x_fig3, fit - y_fig3, 'ok')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Residuals")
plt.show()
print("")

print("fig. 3 limited to >460 nm")
x_fig3b = x_fig3[7:]
y_fig3b = y_fig3[7:]
line = scipy.stats.linregress(x_fig3b, y_fig3b)
print(f"R-squared: {line.rvalue**2:.4f}")
print("Slope: " + str(line.slope))
print(f"P-value: {line.pvalue:.4f}")
fit = line.intercept + line.slope*x_fig3b
plt.plot(x_fig3b, y_fig3b, 'ok')
plt.plot(x_fig3b, fit, 'k')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Log sensitivity difference")
plt.show()
# residual
plt.plot(x_fig3b, fit - y_fig3b, 'ok')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Residuals")
plt.show()

# fig. 1
print("Modified fig. 1")
media_fig1 = np.empty(20)
for i in range(20):
        media_fig1[i] = divide(y_fig1[i], log(media_10nm[i+16]))
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig1, media_fig1, p0=[lw, 1])
opt_raw()
print("LWS: " + opt(popt, pcov))
curve = vpt_fit1(x_1nm, *popt)
rmse(media_fig1, curve, 2, start=160, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig1, media_fig1, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")

# fig. 1+3
print("fig. 1 + 3")
media_fig13 = np.empty(20)
for i in range(20):
        media_fig13[i] = divide(y_fig13[i], log(media_10nm[i+16]))
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig1, media_fig13, p0=[lw, 1])
opt_raw()
print("LWS: " + opt(popt, pcov))
curve = vpt_fit1(x_1nm, *popt)
rmse(media_fig13, curve, 2, start=160, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig1, media_fig13, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")

# 2 templates, L fixed
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixl, x_fig1, media_fig13, p0=[mw, 1, 0.1])
opt_raw()
print("LWS: " + str(lw) + " nm (fixed) (" + opt(popt, pcov, 1, 2) + ")")
print("MWS: " + opt(popt, pcov) + " (" + opt(popt, pcov, 2, 2) + ")")
opt_percent()
curve = vpt_fit2_fixl(x_1nm, *popt)
rmse(media_fig13, curve, 3, start=160, interval=10)
lws = np.empty(401)
mws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, lw) * popt[1])
        mws[i] = log(c.vpt(i + 300, popt[0]) * popt[2])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig1, media_fig13, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, mws, '-.k')
plt.show()
print("")

# fix both
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixall, x_fig1, media_fig13, p0=[1, 0.1])
opt_raw()
print("LWS: " + str(lw) + " nm (fixed) (" + opt(popt, pcov, r=2) + ")")
print("MWS: " + str(mw) + " nm (fixed) (" + opt(popt, pcov, 1, 2) + ")")
opt_percent(0)
curve = vpt_fit2_fixall(x_1nm, *popt)
rmse(media_fig13, curve, 2, start=160, interval=10)
lws = np.empty(401)
mws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, lw) * popt[0])
        mws[i] = log(c.vpt(i + 300, mw) * popt[1])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig1, media_fig13, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, mws, '-.k')
plt.show()
print("")

# fig. 2
print("Modified fig. 2")

# We're "expecting" to find we can fit a UV pigment in, but let's also compare
# just one template.
media_fig2 = np.empty(31)
for i in range(31):
        media_fig2[i] = divide(y_fig2[i], log(media_10nm[i+7]))
popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig2, media_fig2, p0=[lw, 1])
opt_raw()
print("LWS: " + opt(popt, pcov))
curve = vpt_fit1(x_1nm, *popt)
rmse(media_fig2, curve, 2, start=70, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig2, media_fig2, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")

# Residual plots: https://statisticsbyjim.com/regression/check-residual-plots-regression-analysis/
# two types of residual plots: vs. independent variable and vs. fitted values
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

popt, pcov, infodict, mesg, ier = scipy.optimize.curve_fit(vpt_fit2, x_fig2, media_fig2, p0=[lw, 360, 1, 0.1], full_output=True)
opt_raw()
print("LWS: " + opt(popt, pcov) + " (" + opt(popt, pcov, 2, 2) + ")")
print("SWS: " + opt(popt, pcov, 1) + " (" + opt(popt, pcov, 3, 2) + ")")
opt_percent(2)
curve = vpt_fit2(x_1nm, *popt)
rmse(media_fig2, curve, 4, start=70, interval=10)
lws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, popt[0]) * popt[2])
        sws[i] = log(c.vpt(i + 300, popt[1]) * popt[3])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig2, media_fig2, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, sws, ':k')
plt.show()
print("")

# fixed S
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixs, x_fig2, media_fig2, p0=[lw, 1, 0.1])
opt_raw()
print("LWS: " + opt(popt, pcov) + " (" + opt(popt, pcov, 1, 2) + ")")
print("SWS: " + str(sw) + " nm (fixed) (" + opt(popt, pcov, 2, 2) + ")")
opt_percent()
curve = vpt_fit2_fixs(x_1nm, *popt)
rmse(media_fig2, curve, 3, start=70, interval=10)
lws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, popt[0]) * popt[1])
        sws[i] = log(c.vpt(i + 300, sw) * popt[2])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig2, media_fig2, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, sws, ':k')
plt.show()
print("")

# fig. 2 + 3
print("fig. 2 + 3")
# just 1
media_fig23 = np.empty(29)
for i in range(29):
        media_fig23[i] = divide(y_fig23[i], log(media_10nm[i+9]))

popt, pcov = scipy.optimize.curve_fit(vpt_fit1, x_fig23, media_fig23, p0=[lw, 1])
opt_raw()
print("LWS: " + opt(popt, pcov))
curve = vpt_fit1(x_1nm, *popt)
rmse(media_fig23, curve, 1, start=90, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")

# L and S
popt, pcov = scipy.optimize.curve_fit(vpt_fit2, x_fig23, media_fig23, p0=[lw, 360, 1, 0.1])
opt_raw()
print("LWS: " + opt(popt, pcov) + " (" + opt(popt, pcov, 2, 2) + ")")
print("SWS: " + opt(popt, pcov, 1) + " (" + opt(popt, pcov, 3, 2) + ")")
opt_percent(2)
curve = vpt_fit2(x_1nm, *popt)
rmse(media_fig23, curve, 4, start=90, interval=10)
lws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, popt[0]) * popt[2])
        sws[i] = log(c.vpt(i + 300, popt[1]) * popt[3])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, sws, ':k')
plt.show()

##fitted = vpt_fit2(x_fig23, *popt)
##residual = fitted - media_fig23
##plt.plot(x_fig23, residual, 'o')
##plt.xlabel("Wavelength (nm)")
##plt.ylabel("Residuals")
##plt.show()
print("")

# fixed S
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixs, x_fig23, media_fig23, p0=[lw, 1, 0.1])
opt_raw()
total = popt[1] + popt[2]
print("LWS: " + opt(popt, pcov) + " (" + opt(popt, pcov, 1, 2) + ")")
print("SWS: " + str(sw) + " nm (fixed) (" + opt(popt, pcov, 2, 2) + ")")
opt_percent()
curve = vpt_fit2_fixs(x_1nm, *popt)
rmse(media_fig23, curve, 4, start=90, interval=10)
lws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, popt[0]) * popt[1])
        sws[i] = log(c.vpt(i + 300, sw) * popt[2])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, sws, ':k')
plt.show()
print("")

# fixed L, second variable (S doesn't work well so we're just doing L and M)
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixl, x_fig23, media_fig23, p0=[mw, 1, 0.1])
opt_raw()
print("LWS: " + str(lw) + " nm (fixed) (" + opt(popt, pcov, 1, 2) + ")")
print("MWS: " + opt(popt, pcov) + " (" + opt(popt, pcov, 2, 2) + ")")
opt_percent()
curve = vpt_fit2_fixl(x_1nm, *popt)
rmse(media_fig23, curve, 3, start=90, interval=10)
lws = np.empty(401)
mws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, lw) * popt[1])
        mws[i] = log(c.vpt(i + 300, popt[0]) * popt[2])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, mws, '-.k')
plt.show()
print("")

# fixed L and M
popt, pcov = scipy.optimize.curve_fit(vpt_fit2_fixall, x_fig23, media_fig23, p0=[1, 0.1])
opt_raw()
print("LWS: " + str(lw) + " nm (fixed) (" + opt(popt, pcov, r=2) + ")")
print("MWS: " + str(mw) + " nm (fixed) (" + opt(popt, pcov, 1, 2) + ")")
opt_percent(0)
curve = vpt_fit2_fixall(x_1nm, *popt)
rmse(media_fig23, curve, 2, start=90, interval=10)
lws = np.empty(401)
mws = np.empty(401)
sws = np.empty(401)
for i in range(401):
        lws[i] = log(c.vpt(i + 300, lw) * popt[0])
        mws[i] = log(c.vpt(i + 300, mw) * popt[1])
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.plot(x_1nm, lws, '--k')
plt.plot(x_1nm, mws, '-.k')
plt.show()
print("")

# self-screening
p0 = [lw, 1, 1]
# fig. 1
print("Modified fig. 1")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1_od, x_fig1, media_fig1, p0=p0)
opt_raw()
print("LWS: " + opt(popt, pcov))
print("OS optical density: " + opt(popt, pcov, 2, 3))
print("Predicted OS length (od="+str(args.od)+"): " + str(round(popt[2]/args.od, 3)))
print("Predicted absorption coefficient (osl="+str(args.osl)+"): " + str(round(popt[2]/args.osl, 3)))
curve = vpt_fit1_od(x_1nm, *popt)
rmse(media_fig1, curve, 3, start=160, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig1, media_fig1, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")

# fig. 1+3
print("fig. 1 + 3")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1_od, x_fig1, media_fig13, p0=p0)
opt_raw()
print("LWS: " + opt(popt, pcov))
print("OS optical density: " + opt(popt, pcov, 2, 3))
print("Predicted OS length (od="+str(args.od)+"): " + str(round(popt[2]/args.od, 3)))
print("Predicted absorption coefficient (osl="+str(args.osl)+"): " + str(round(popt[2]/args.osl, 3)))
curve = vpt_fit1_od(x_1nm, *popt)
rmse(media_fig13, curve, 3, start=160, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig1, media_fig13, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")

# 2
print("fig. 2")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1_od, x_fig2, media_fig2, p0=p0)
opt_raw()
print("LWS: " + opt(popt, pcov))
print("OS optical density: " + opt(popt, pcov, 2, 3))
print("Predicted OS length (od="+str(args.od)+"): " + str(round(popt[2]/args.od, 3)))
print("Predicted absorption coefficient (osl="+str(args.osl)+"): " + str(round(popt[2]/args.osl, 3)))
curve = vpt_fit1_od(x_1nm, *popt)
rmse(media_fig2, curve, 2, start=70, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig2, media_fig2, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")

# 2+3
print("fig. 2 + 3")
popt, pcov = scipy.optimize.curve_fit(vpt_fit1_od, x_fig23, media_fig23, p0=p0)
opt_raw()
print("LWS: " + opt(popt, pcov))
print("OS optical density: " + opt(popt, pcov, 2, 3))
print("Predicted OS length (od="+str(args.od)+"): " + str(round(popt[2]/args.od, 3)))
print("Predicted absorption coefficient (osl="+str(args.osl)+"): " + str(round(popt[2]/args.osl, 3)))
curve = vpt_fit1_od(x_1nm, *popt)
rmse(media_fig23, curve, 1, start=90, interval=10)
plt.xlabel("Wavelength (nm)")
plt.ylabel(ylabel)
#plt.yscale("log")
ylim()
plt.plot(x_fig23, media_fig23, 'ok')
plt.plot(x_1nm, curve, 'k')
plt.show()
print("")
