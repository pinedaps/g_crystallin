import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Heat capacities [J/K/mol] for aminoacids
 
df = pd.read_excel("Cp_amino_acids.xlsx")
T  = np.array(df.columns.values.tolist()[1:7])+273

coefficients = []
for i, row in df.iterrows():
	Cp        = (row[1:7].values).tolist()
	plt.plot(T, Cp, label=row["Residue"])	
	coef      = np.polyfit(T,Cp,3)
	poly1d_fn = np.poly1d(coef)
	plt.plot(T, poly1d_fn(T), '--k')
	coefficients.append((row["Residue"],float(coef[0]),float(coef[1]),float(coef[2]),float(coef[3])))
cf_df = pd.DataFrame(coefficients, columns=['Residue','a','b','c','d'])
cf_df.to_csv('Cp_coefficients.csv', index=False)
plt.xlabel("Temperature (T)")
plt.ylabel("Heat Capacity (Cp)")
plt.title("Fitting of Amino Acid Heat Capacities vs Temperature")
plt.xlim((270,430))
plt.legend(labelspacing=0.23)
plt.savefig('cp_aa.png')
plt.close()

df_test = pd.read_csv("Cp_coefficients.csv")
for i, row in df_test.iterrows():
        coef      = (row[1:].values).tolist()
        poly1d_fn = np.poly1d(coef)
        plt.plot(T, poly1d_fn(T), '--', label=row["Residue"])
plt.xlabel("Temperature (T)")
plt.ylabel("Heat Capacity (Cp)")
plt.title("Fitted Amino Acid Heat Capacities vs Temperature")
plt.xlim((270,430))
plt.legend(labelspacing=0.23)
plt.savefig('cp_aa_fit.png')

