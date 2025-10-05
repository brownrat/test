# peaked or "band-pass" spectra. Blue, green and purple objects usually have this type
# of reflectance, including green leaves and some flowers. Modeled with a normal distribution.
# See https://mathworld.wolfram.com/NormalDistribution.html
# The peak region is usually set to about 100 nm wide by choosing 50 as the standard
# deviation. The width at half-maximum is slightly narrower.

import central as c
import numpy as np
import math

# 350-nm peak
peak350 = np.empty(401)
for i in range(401):
        peak350[i] = math.exp(-(i+300-350)**2/(2*25)**2)
print("350-nm peak")
c.spectral_rendering(peak350)

# 400-nm peak
peak400 = np.empty(401)
for i in range(401):
        peak400[i] = math.exp(-(i+300-400)**2/(2*25)**2)
print("400-nm peak")
c.spectral_rendering(peak400)

# 450-nm peak
peak450 = np.empty(401)
for i in range(401):
        peak450[i] = math.exp(-(i+300-450)**2/(2*25)**2)
print("450-nm peak")
c.spectral_rendering(peak450)

# 500-nm peak
peak500 = np.empty(401)
for i in range(401):
        peak500[i] = math.exp(-(i+300-500)**2/(2*25)**2)
print("500-nm peak")
c.spectral_rendering(peak500)

# 550-nm peak
peak550 = np.empty(401)
for i in range(401):
        peak550[i] = math.exp(-(i+300-550)**2/(2*25)**2)
print("550-nm peak")
c.spectral_rendering(peak550)

# 600-nm peak
peak600 = np.empty(401)
for i in range(401):
        peak600[i] = math.exp(-(i+300-600)**2/(2*25)**2)
print("600-nm peak")
c.spectral_rendering(peak600)

# 650-nm peak
peak650 = np.empty(401)
for i in range(401):
        peak650[i] = math.exp(-(i+300-650)**2/(2*25)**2)
print("650-nm peak")
c.spectral_rendering(peak650)
