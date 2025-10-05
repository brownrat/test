"""
ramp functions
Many natural surfaces have this type of reflectance, especially those that appear
"brown"/"beige". Examples:
* soil and vegetation: https://www.researchgate.net/profile/Yingying-Yang-5/publication/340281560/figure/fig5/AS:878006243430404@1586344409333/a-A-comparison-of-spectral-reflectance-curves-between-soil-and-withered-vegetation.png
* different soil types: https://www.researchgate.net/publication/47426815/figure/fig1/AS:11431281212212810@1702556618056/Typical-spectral-reflectance-of-vegetation-soils-water-and-snow-The-original-data.tif
* fur: https://www.mdpi.com/2072-4292/8/4/273
* red vs. green leaves: http://2.bp.blogspot.com/-DfLsYNxVoMo/UbwHSrq_KxI/AAAAAAAAAYM/OZyOnBPRxf4/s1600/LeafSpectra.jpg
* human skin:
** https://www.researchgate.net/figure/Skin-reflectance-curves-for-individuals-of-European-East-Asian-South-Asian-and-African_fig8_5800680
** https://www.researchgate.net/figure/Skin-reflectance-in-the-wavelength-regions-of-UVB-280-320-nm-UVA-320-400-nm-and_fig1_228955239
* human hair: https://media.springernature.com/lw685/springer-static/image/art%3A10.1007%2Fs10439-019-02315-z/MediaObjects/10439_2019_2315_Fig5_HTML.png

The color of red/yellow leaves is mainly due to carotenoids, and carotenoid spectra
seem to be somewhere between a step and a ramp.

Some mammals with white fur (e.g. polar bear) absorb UV whereas others (e.g. Arctic
hare) do not. An example of non-absorbing white fur is the Honduran white bat (see
carotenoid example).
"""
import central as c
import numpy as np

# 300-nm ramp
ramp300 = np.empty(401)
for i in range(401): ramp300[i] = i/401
print("300-nm ramp function")
c.spectral_rendering(ramp300)

# 350-nm ramp
ramp350 = np.empty(401)
for i in range(49): ramp350[i] = 0
for i in range(50, 401): ramp350[i] = (i-50)/351
print("350-nm ramp function")
c.spectral_rendering(ramp350)

# 400-nm ramp
ramp400 = np.empty(401)
for i in range(99): ramp400[i] = 0
for i in range(100, 401): ramp400[i] = (i-100)/301
print("400-nm ramp function")
c.spectral_rendering(ramp400)

# 450-nm ramp
ramp450 = np.empty(401)
for i in range(149): ramp450[i] = 0
for i in range(150, 401): ramp450[i] = (i-150)/251
print("450-nm ramp function")
c.spectral_rendering(ramp450)

# 500-nm ramp
ramp500 = np.empty(401)
for i in range(199): ramp500[i] = 0
for i in range(200, 401): ramp500[i] = (i-200)/201
print("500-nm ramp function")
c.spectral_rendering(ramp500)

# 550-nm ramp
ramp550 = np.empty(401)
for i in range(249): ramp550[i] = 0
for i in range(250, 401): ramp550[i] = (i-250)/151
print("550-nm ramp function")
c.spectral_rendering(ramp550)

# 600-nm ramp
ramp600 = np.empty(401)
for i in range(299): ramp600[i] = 0
for i in range(300, 401): ramp600[i] = (i-300)/101
print("600-nm ramp function")
c.spectral_rendering(ramp600)
