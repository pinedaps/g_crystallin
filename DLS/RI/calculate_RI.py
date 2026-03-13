import numpy as np

# Wavelengths in nm
lam = np.array([435.8, 546.1, 579.0])

# Corresponding refractive indices
# Buffer
#n = np.array([(1.34084+1.34013+1.34071)/3,(1.33495+1.33417+1.33473)/3,(1.33370+1.33347+1.33374)/3])
# 33 %
#n = np.array([(1.34063+1.34028+1.34104)/3,(1.33469+1.33482+1.33468)/3,(1.33367+1.33342+1.33348)/3])
# 67 %
#n = np.array([(1.34078+1.34130+1.34066)/3,(1.33477+1.33487+1.33491)/3,(1.33337+1.33366+1.33374)/3])
# 100 %
n = np.array([(1.34027+1.34125+1.33987)/3,(1.33486+1.33487+1.33495)/3,(1.33372+1.33356+1.33393)/3])

x = 1.0 / lam**2

x_bar = np.mean(x)
n_bar = np.mean(n)

B = np.sum((x - x_bar) * (n - n_bar)) / np.sum((x - x_bar)**2)
A = n_bar - B * x_bar

print("Fitted coefficients:")
print("A =", A)
print("B =", B)

n_fit = A + B / lam**2
residuals = n - n_fit

print("\nResiduals:", residuals)
print("RMS error:", np.sqrt(np.mean(residuals**2)))

lambda_new = 660 

n_new = A + B / lambda_new**2

print(f"\nPredicted refractive index at λ = {lambda_new}:")
print("n =", n_new)
