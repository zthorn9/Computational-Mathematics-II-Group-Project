import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import label, sum

def cluster_sizes(L, p, nsamp):
    sizes = []
    structure = np.array([1, 1, 1])
    for i in range(nsamp):
        z = np.random.rand(L)
        m = z < p
        lw, num = label(m)
        if num > 0:
#             cluster_sizes = sum(m, lw, index=np.arange(1, num + 1))
            cluster_sizes = sum(m, lw, index=np.arange(lw.max() + 1))
            sizes.extend(cluster_sizes)
    return np.array(sizes)

L = 1000        
nsamp = 1000
# p values to investigate average cluster size against occupation probability p
p_vals = np.linspace(0.5, 0.99, 50)
S_measured = np.zeros(len(p_vals))
S_theoretical = np.zeros(len(p_vals))
s_xi_measured = np.zeros(len(p_vals))
# begin simulation
for i,p in enumerate(p_vals):
    s_xi_measured[i] = -1.0/np.log(p)
    sizes = cluster_sizes(L, p, nsamp)
    if sizes.size > 0:
        summed_s_sq = np.sum(sizes**2)
        summed_s = np.sum(sizes)
        S_measured[i] = summed_s_sq / summed_s
    # Calculate theoretical value
    S_theoretical[i] = (1 + p) / (1 - p)

fit_points = 15 

p_data = p_vals[-fit_points:]
S_data = S_measured[-fit_points:]
Xi_data = s_xi_measured[-fit_points:]
x_log = np.log(1 - p_data)
coefficients_gamma = np.polyfit(x_log, np.log(S_data), 1) 
gamma_sim = -coefficients_gamma[0] 
coefficients_nu = np.polyfit(x_log, np.log(Xi_data), 1) 
nu_sim = -coefficients_nu[0] 

print(f"Theoretical Gamma (γ) = 1.0")
print(f"Simulated Gamma (γ) ≈ {gamma_sim:.3f}")
print(f"Theoretical Nu (ν) = 1.0")
print(f"Simulated Nu (ν) ≈ {nu_sim:.3f}")

p_test = 0.9
sizes_n = cluster_sizes(L, p_test, nsamp)
if sizes_n.size > 0:
    s_max = sizes_n.max()
    n, sbins = np.histogram(sizes_n, bins = int(s_max) + 1, range = (0, s_max + 1))
    s_values = 0.5*(sbins[1:] + sbins[:-1])
    n_s_p_measured = n/(L*nsamp)
    s_xi = -1.0/np.log(p_test)
    n_s_p_theoretical = np.exp(-s_values/s_xi)*(1-p_test)**2
else:
    s_values = np.array([0])
    n_s_p_measured = np.array([0])
    n_s_p_theoretical = np.array([0])

print(f"Characteristic cluster size = {s_xi:.2f}")    

# Figure 1: log(n(s,p) against s
plt.semilogy(s_values[1:],n_s_p_measured[1:],'bo', label = 'Numerical Result n(s,p)')
plt.semilogy(s_values[1:], n_s_p_theoretical[1:], 'r-', label=f'Theory: $\\propto e^{{-s/s_\\xi}}$ (where $s_\\xi={s_xi:.2f}$)')
plt.xlabel('Cluster size s',fontsize = 15)
plt.ylabel('log(n(s,p))',fontsize = 15)
plt.title(f'Cluster Number density n(s,p) at p={p_test}',fontsize = 15)
plt.legend(fontsize = 15)
plt.grid(True)
plt.xlim(0, 50) 
plt.ylim(bottom=1e-10)
plt.show()

# Figure 2
# Plot 1: Log-Log plot to find exponent gamma
plt.subplot(1, 2, 1)
plt.loglog(1 - p_vals, S_measured, 'bo', label='Measured S')
plt.loglog(1 - p_vals, S_theoretical, 'r-', label='Theoretical S')
plt.xlabel('Log(1 - p)',fontsize = 15)
plt.ylabel('Log(S)',fontsize = 15)
plt.legend(fontsize = 15)
plt.grid(True)
plt.title('Log-Log plot of Average Cluster Size S against (1-p)',fontsize = 13)

# Plot 2: log plot of (1-p) vs S_xi
plt.subplot(1, 2, 2)
plt.loglog(1 - p_vals, s_xi_measured, 'bo', label='Numerical $s_\\xi$')
plt.loglog(1 - p_vals, 1/(1 - p_vals), 'r-', label='Theoretical $s_\\xi \\propto (1-p)^{-1}$')
plt.xlabel('Log(1 - p)',fontsize = 15)
plt.ylabel('Log($s_\\xi$)',fontsize = 15)
plt.title('Log-Log Plot of Characteristic Cluster $s_\\xi$ against (1-p)',fontsize = 13)
plt.legend(fontsize = 15)
plt.grid(True)

plt.show()