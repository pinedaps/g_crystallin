import numpy as np
# -----------------------------
# Data (mg/mL)
# -----------------------------
c_mgml = np.array([0.26893, 0.79647, 1.06413])
n = np.array([1.33141, 1.33143, 1.33175])

# Convert to g/mL
c = c_mgml / 1000.0

# -----------------------------
# Linear fit: n = a + b*c
# -----------------------------
c_bar = np.mean(c)
n_bar = np.mean(n)

b = np.sum((c - c_bar)*(n - n_bar)) / np.sum((c - c_bar)**2)
a = n_bar - b*c_bar

print("Intercept (n0) =", a)
print("dn/dc =", b, "mL/g")

# -----------------------------
# Residuals (keep your structure)
# -----------------------------
n_fit = a + b*c
residuals = n - n_fit

print("\nResiduals:", residuals)

RSS = np.sum(residuals**2)
RMS = np.sqrt(np.mean(residuals**2))

print("RMS error:", RMS)
