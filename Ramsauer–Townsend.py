#%%
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import math
#%%
## Voltage values
V = np.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,])
Vp = np.array([0.2e-3,0.3e-3,0.6e-3,0.9e-3,1.2e-3,1.7e-3,2.2e-3,2.9e-3,3.5e-3,4.2e-3,5.0e-3,5.8e-3,6.6e-3,7.3e-3,8.0e-3,8.7e-3,9.4e-3,9.9e-3,10.4e-3,10.9e-3,11.4e-3,13.3e-3,12.3e-3,10.8e-3,9.9e-3,9.6e-3,9.5e-3,9.7e-3,9.9e-3,10.3e-3,10.9e-3,11.8e-3,12.7e-3,13.9e-3,15.3e-3,17.1e-3,19.4e-3,22.1e-3,26.1e-3,30.9e-3,37.2e-3,47.8e-3,67.2e-3])
Vp_star = np.array([1.5e-3,2.0e-3,2.6e-3,3.2e-3,3.9e-3,4.6e-3,5.4e-3,6.2e-3,6.9e-3,7.7e-3,8.6e-3,9.4e-3,10.3e-3,11.2e-3,12.1e-3,13.2e-3,14.2e-3,15.2e-3,16.2e-3,17.3e-3,18.4e-3,30.3e-3,42.6e-3,55.4e-3,69.3e-3,84.2e-3,100.2e-3,117.3e-3,135.2e-3,153.6e-3,172.2e-3,190.2e-3,0.207,0.225,0.245,0.268,0.292,0.314,0.333,0.349,0.363,0.382,0.404,])
Vs = np.array([1.1e-3,1.7e-3,2.4e-3,3.3e-3,4.5e-3,5.8e-3,7.1e-3,8.6e-3,10.1e-3,11.8e-3,13.6e-3,15.6e-3,17.5e-3,19.5e-3,21.6e-3,23.7e-3,25.8e-3,27.9e-3,30.2e-3,32.5e-3,34.9e-3,60.3e-3,88.8e-3,118.6e-3,148.9e-3,179.4e-3,0.209,0.241,0.274,0.308,0.343,0.379,0.418,0.457,0.498,0.541,0.586,0.631,0.677,0.727,0.779,0.842,0.933,])
Vs_star = np.array([2.7e-3,3.6e-3,4.7e-3,6.1e-3,7.5e-3,9.1e-3,10.9e-3,12.6e-3,14.4e-3,16.2e-3,18.1e-3,20.0e-3,22.1e-3,24.2e-3,26.4e-3,28.9e-3,31.1e-3,33.2e-3,35.6e-3,38.0e-3,40.6e-3,66.8e-3,96.4e-3,128.1e-3,161.0e-3,195.1e-3,0.229,0.267,0.305,0.346,0.390,0.437,0.485,0.536,0.588,0.642,0.698,0.754,0.810,0.867,0.924,0.979,1.031,])

## Plate current values
Ip = Vp / 10000
Istar_p = Vp_star / 10000

Is = Vs / 100
Istar_s = Vs_star / 100

## Uncertainties 
sigma_Vp = 0.1E-3
sigma_Vs = 0.1E-3
sigma_Vp_star = 0.1E-3
sigma_Vs_star = 0.1E-3
sigma_Ip = sigma_Vp / 10000
sigma_Is = sigma_Vs / 100
sigma_Istar_p = sigma_Vp_star / 10000
sigma_Istar_s = sigma_Vs_star / 100

## Plotting plate current vs voltage
plt.plot(V, Ip, 'bo-', label = '$I_p$')
plt.errorbar(V, Ip, xerr= 0.1, yerr = sigma_Ip, ls = 'none', capsize=5, label = 'Error on $I_p$')
plt.plot(V, Istar_p, 'ro-', label='$I_{p}^{*}$')
plt.errorbar(V, Istar_p, xerr= 0.1, yerr = sigma_Istar_p, ls = 'none', capsize=5, label = 'Error on $I_{p}^{*}$')

plt.title('Plate Current vs Voltage')
plt.xlabel('Voltage (V)')
plt.ylabel('Plate Current (A)')
plt.legend()
plt.grid()

plt.show()
#%%
## Probability of scattering
P = 1 - (Ip * Istar_s)/(Istar_p * Is)

## Probability error
def P_error(P, Ip, Istar_s, Istar_p, Is, sigma_Ip, sigma_Istar_s, sigma_Istar_p, sigma_Is):
    sigma_P = math.sqrt(((sigma_Ip * Istar_s / (Istar_p * Is)) ** 2) + ((Ip * sigma_Istar_s / (Istar_p * Is)) ** 2) + ((Ip * Istar_s * sigma_Istar_p / (Istar_p ** 2 * Is)) ** 2) + ((Ip * Istar_s * sigma_Is / (Istar_p * Is ** 2)) ** 2))
    return sigma_P

sigma_P_array = []
for i in range(len(P)):
    sigma_P = P_error(P[i], Ip[i], Istar_s[i], Istar_p[i], Is[i], sigma_Ip, sigma_Istar_s, sigma_Istar_p, sigma_Is)
    sigma_P_array.append(sigma_P)

## Plotting probability of scattering vs electron momentum 
plt.errorbar((V-Vs)**1/2, P, yerr = sigma_P,ls='none', capsize=7, label = 'Error on Probability')
plt.plot((V-Vs)**1/2, P, 'ro-', label='P')

plt.title('Probability of Scattering vs Electron Momentum')
plt.xlabel('Electron Momentum $(\sqrt{\mathrm{Volts}})$')
plt.ylabel('Probability of Scattering $P_s$') 
plt.legend()
plt.grid()

plt.show()

## Electron mmtm. corresponding to min and max probabilities
x_a = (V-Vs)**1/2
min_y_a = np.min(P)
min_y_a_index = np.argmin(P)
min_x_a = x_a[min_y_a_index]

max_y_a = np.max(P)
max_y_a_index = np.argmax(P)
max_x_a = x_a[max_y_a_index]

print(f"The minimum probability of scattering is: {min_y_a:.2f}")
print(f"Corresponding energy value: {min_x_a:.2f}")

print(f"The maximum probability of scattering is: {max_y_a:.2f}")
print(f"Corresponding energy value: {max_x_a:.2f}")

# %%
## Reciprocal mean free path vs electron momentum
MeanFreePathRecip = -(np.log(1-P))*(1/0.7)

plt.plot((V-Vs)**1/2, MeanFreePathRecip, 'ro-')

plt.title('Reciprocal Mean Free Path vs Electron Momentum')
plt.xlabel('Electron Momentum $(\sqrt{\mathrm{Volts}}$)')
plt.ylabel('Reciprocal Mean Free Path $(\mathrm{cm}^{-1})$') 
plt.grid()

plt.show()

# %%
## Extension with reversed polarity
Vs_star_reverse = np.array([2.3e-3,2.2e-3,2.0e-3,1.9e-3,1.7e-3,1.6e-3,1.5e-3,1.4e-3,1.4e-3,1.3e-3,1.2e-3,1.1e-3,1.0e-3,0.9e-3,0.8e-3,0.7e-3,0.7e-3,0.6e-3,0.6e-3,0.5e-3,0.4e-3,0.4e-3,0.3e-3,0.3e-3,0.3e-3,0.2e-3,0.2e-3,0.2e-3,0.2e-3,0.1e-3,0.1e-3,0.1e-3,0.1e-3,0.1e-3,0.1e-3])
Vrange = np.array([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34])

Is_star_reverse = Vs_star_reverse / 100

## Fit a line to the first slope
x1 = -Vrange[:20]
y1 = Is_star_reverse[:20]
m1, c1, r1, p1, std_err1 = linregress(x1, np.log10(y1))
line1 = 10 ** (m1 * x1 + c1)
m1_error1 = std_err1

## Fit a line to the second slope
x2 = -Vrange[10:]
y2 = Is_star_reverse[10:]
m2, c2, r2, p2, std_err2 = linregress(x2, np.log10(y2))
line2 = 10 ** (m2 * x2 + c2)
m2_error2 = std_err2

## Plot the graph with the lines and error bars
plt.plot(-Vrange, Is_star_reverse, 'r.', label = 'Experimental Data')
plt.plot(x1, line1, 'b--', label = f'Right slope = {m1:.2f}±{m1_error1:.2f}')
plt.plot(x2, line2, 'g--', label = f'Left slope = {m2:.2f}±{m2_error2:.2f}')

plt.title('$\log(I_{s}^{*})$ vs Potential Between Shield and Cathode')
plt.yscale('log')
plt.xlabel('Potential between shield and cathode')
plt.ylabel('$\log(I_{s}^{*})$')
plt.grid()
plt.legend()

Voltage_intersection = (c2 - c1) / (m1 - m2)
print(f"The intersection point is at {Voltage_intersection:.2f} V.")

# %%
Vc = 0.14
Vmean = 3 / (2 * m2)

plt.errorbar((V-Vs+Vc+Vmean)**1/2, P, yerr = sigma_P,ls='none', capsize=7, label = 'Error on Probability')
plt.plot((V-Vs+Vc+Vmean)**1/2, P, 'ro-', label='P')

plt.title('Probability of Scattering vs Electron Momentum')
plt.xlabel('Electron Momentum $\sqrt{V - V_{s} + V_{c} + \overline{V}} \, (\sqrt{\mathrm{Volts}})$')
plt.ylabel('Probability of Scattering $P_s$') 
plt.legend()
plt.grid()

plt.show()

## Electron mmtm. corresponding to to min and max probabilities (adjusted values)
x_b = (V-Vs+Vc+Vmean)**1/2
min_y_b = np.min(P)
min_y_b_index = np.argmin(P)
min_x_b = x_b[min_y_b_index]

max_y_b = np.max(P)
max_y_b_index = np.argmax(P)
max_x_b = x_b[max_y_b_index]

print(f"The minimum probability of scattering is: {min_y_b:.2f}")
print(f"Corresponding energy value: {min_x_b:.2f}")

print(f"The maximum probability of scattering is: {max_y_b:.2f}")
print(f"Corresponding energy value: {max_x_b:.2f}")

print(f"Difference in electron momenta for minima: {min_x_b-min_x_a:.2f}")
print(f"Difference in electron momenta for maxima: {max_x_b-max_x_a:.2f}")
# %%
print(c1) ## Intercept of line 1 
print(c2) ## Intercept of line 2 
print(Vmean)