# step functions

import central as c
import numpy as np

# "long-pass" step functions: 0 for wavelengths less than the step and 1 for longer wavelengths.
# Many white, yellow, orange and red flowers have reflectance spectra like these, and they're
# also found in paints and camera filters.
print("Long-pass step functions")
# 350-nm step
step350 = np.empty(401)
for i in range(49): step350[i] = 0
for i in range(50, 401): step350[i] = 1
print("350-nm step function")
c.spectral_rendering(step350)

# 400-nm step
step400 = np.empty(401)
for i in range(99): step400[i] = 0
for i in range(100, 401): step400[i] = 1
print("400-nm step function")
c.spectral_rendering(step400)

# 450-nm step
step450 = np.empty(401)
for i in range(149): step450[i] = 0
for i in range(150, 401): step450[i] = 1
print("450-nm step function")
c.spectral_rendering(step450)

# 500-nm step
step500 = np.empty(401)
for i in range(199): step500[i] = 0
for i in range(200, 401): step500[i] = 1
print("500-nm step function")
c.spectral_rendering(step500)

# 550-nm step
step550 = np.empty(401)
for i in range(249): step550[i] = 0
for i in range(250, 401): step550[i] = 1
print("550-nm step function")
c.spectral_rendering(step550)

# 600-nm step
step600 = np.empty(401)
for i in range(299): step600[i] = 0
for i in range(299, 401): step600[i] = 1
print("600-nm step function")
c.spectral_rendering(step600)

# 650-nm step
step650 = np.empty(401)
for i in range(349): step650[i] = 0
for i in range(350, 401): step650[i] = 1
print("650-nm step function")
c.spectral_rendering(step650)

# "short-pass" step functions, the reverse of the above. These form the other edge of the
# optimal color solid (see Wikipedia) but don't seem to occur in nature or have any
# applications. Blue/green surfaces usually are the "band-pass" type with one or two
# peaks.
print("Short-pass step functions")
# 350-nm step
step350 = np.empty(401)
for i in range(49): step350[i] = 1
for i in range(50, 401): step350[i] = 0
print("350-nm step function")
c.spectral_rendering(step350)

# 400-nm step
step400 = np.empty(401)
for i in range(99): step400[i] = 1
for i in range(100, 401): step400[i] = 0
print("400-nm step function")
c.spectral_rendering(step400)

# 450-nm step
step450 = np.empty(401)
for i in range(149): step450[i] = 1
for i in range(150, 401): step450[i] = 0
print("450-nm step function")
c.spectral_rendering(step450)

# 500-nm step
step500 = np.empty(401)
for i in range(199): step500[i] = 1
for i in range(200, 401): step500[i] = 0
print("500-nm step function")
c.spectral_rendering(step500)

# 550-nm step
step550 = np.empty(401)
for i in range(249): step550[i] = 1
for i in range(250, 401): step550[i] = 0
print("550-nm step function")
c.spectral_rendering(step550)

# 600-nm step
step600 = np.empty(401)
for i in range(299): step600[i] = 1
for i in range(300, 401): step600[i] = 0
print("600-nm step function")
c.spectral_rendering(step600)

# 650-nm step
step650 = np.empty(401)
for i in range(349): step650[i] = 1
for i in range(350, 401): step650[i] = 0
print("650-nm step function")
c.spectral_rendering(step650)
