import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import label, sum

def cluster_sizes_2D(L, p, nsamp):
    n_s_tot = {}
    sites_tot = L * L
    structure = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]
    for _ in range(nsamp):
        z = np.random.rand(L,L)
        m = z < p
        lw, num = label(m, structure = structure)
        if num > 0:
            cluster_sizes = sum(m, lw, index=np.arange(1, num + 1))
            sizes_uni, counts = np.unique(cluster_sizes, return_counts = True)
            for s, N_s in zip(sizes_uni, counts):
                density = N_s/sites_tot
                n_s_tot.setdefault(s, []).append(density)
    avg_sizes = np.array(list(n_s_tot.keys()))
    avg_densities = np.array([np.mean(n_list) for n_list in n_s_tot.values()])
    sort_ind = np.argsort(avg_sizes)
    return avg_sizes[sort_ind], avg_densities[sort_ind]

def cluster_sizes_3D(L, p, nsamp):
    n_s_tot = {}
    sites_tot = L * L * L
    structure = [[[0,0,0],[0,1,0],[0,0,0]],
        [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
        [[0,0,0],[0,1,0],[0,0,0]]]
    for _ in range(nsamp):
        z = np.random.rand(L,L,L)
        m = z < p
        lw, num = label(m, structure = structure)
        if num > 0:
            cluster_sizes = sum(m, lw, index=np.arange(1, num + 1))
            sizes_uni, counts = np.unique(cluster_sizes, return_counts = True)
            for s, N_s in zip(sizes_uni, counts):
                density = N_s/sites_tot
                n_s_tot.setdefault(s, []).append(density)
    avg_sizes = np.array(list(n_s_tot.keys()))
    avg_densities = np.array([np.mean(n_list) for n_list in n_s_tot.values()])
    sort_ind = np.argsort(avg_sizes)
    return avg_sizes[sort_ind], avg_densities[sort_ind]
# define functions to find tau, s_xi and sigma for plotting to investigate n(s,p)
def tau_measure(sizes, densities):
    min = 3
    max = 70
    mask = (sizes >= min) & (sizes <= max) & (densities > 0)
    if np.sum(mask) < 5:
        return 2.05, None # the theoretical tau for 2D
    log_s = np.log(sizes[mask])
    log_n_s = np.log(densities[mask])
    coeffs = np.polyfit(log_s, log_n_s, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    tau_measured = -slope
    return tau_measured, (log_s, log_n_s, slope, intercept)

def s_xi_measure(sizes, densities, tau_measured):
    mask = densities > 0
    if np.sum(mask) == 0:
        return np.nan
    s = sizes[mask]
    n_s = densities[mask]
    y_scaled = np.log(n_s) + tau_measured * np.log(s)
    max_ind = np.argmax(y_scaled)
    s_xi_measured = s[max_ind]
    return s_xi_measured

def sigma_measure(p_vals, s_xi_vals, p_c):
    mask = (s_xi_vals > 0) & (np.abs(p_vals - p_c) > 1e-6)
    if np.sum(mask) < 3:
        return 0.45, None
    delta_p = np.abs(p_vals[mask] - p_c)
    s_xi = s_xi_vals[mask]
    log_delta_p = np.log(delta_p)
    log_s_xi = np.log(s_xi)
    coeffs = np.polyfit(log_delta_p, log_s_xi, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    sigma_measured =  -1.0 / slope
    return sigma_measured, (log_delta_p, log_s_xi, slope, intercept)
# # execute investigation for 2D site percolation
L_2d = 1000   
nsamp_2d = 1000
p_c_2d = 0.5927
tau_theoretical_2d = 2.05
sigma_theoretical_2d = 0.4
# p values close to p_c = 0.5927.... to investigate exponents
p_vals_2d = np.array([0.54, 0.57, p_c_2d, 0.61, 0.63])
s_xi_measured_2d = np.zeros(len(p_vals_2d))
all_dist_2d = {}
# execute investigation for 3D site percolation
L_3d = 200  
nsamp_3d = 1000
p_c_3d = 0.3116
tau_theoretical_3d = 2.2
sigma_theoretical_3d = 0.45
# p values near p_c = 0.3
p_vals_3d = np.array([0.28,0.285, 0.29,0.3, p_c_3d])
s_xi_measured_3d = np.zeros(len(p_vals_3d))
all_dist_3d = {}
beta_3d = 0.4
# begin simulation for tau (2d)
sizes_p_c_2d, densities_p_c_2d = cluster_sizes_2D(L_2d, p_c_2d, nsamp_2d)
all_dist_2d[p_c_2d] = (sizes_p_c_2d, densities_p_c_2d)
tau_measured_2d, tau_data_2d = tau_measure(sizes_p_c_2d, densities_p_c_2d)
if tau_data_2d:
    log_s_2d, log_n_s_2d, slope_2d, intercept_2d = tau_data_2d

# begin simulation for tau (3d)
sizes_p_c_3d, densities_p_c_3d = cluster_sizes_3D(L_3d, p_c_3d, nsamp_3d)
all_dist_3d[p_c_3d] = (sizes_p_c_3d, densities_p_c_3d)
tau_measured_3d, tau_data_3d = tau_measure(sizes_p_c_3d, densities_p_c_3d)
if tau_data_3d:
    log_s_3d, log_n_s_3d, slope_3d, intercept_3d = tau_data_3d

# begin simulation for sigma 2d
tau_s_xi_2d = tau_measured_2d if tau_data_2d else tau_theoretical_2d
for i,p in enumerate(p_vals_2d):
    if p == p_c_2d:
        s_xi_measured_2d[i] = np.nan
        continue
    sizes_2d, densities_2d = cluster_sizes_2D(L_2d, p, nsamp_2d)
    all_dist_2d[p] = (sizes_2d, densities_2d)
    s_xi_2d = s_xi_measure(sizes_2d, densities_2d, tau_s_xi_2d)
    s_xi_measured_2d[i] = s_xi_2d if not np.isnan(s_xi_2d) and s_xi_2d > 0 else 1.0    
p_data_2d = p_vals_2d
s_xi_data_2d = s_xi_measured_2d
sigma_measured_2d, sigma_data_2d = sigma_measure(p_data_2d, s_xi_data_2d, p_c_2d)
if sigma_data_2d:
    log_delta_p_2d, log_s_xi_2d, slope_2d, intercept_2d = sigma_data_2d

# begin simulation for sigma 3d
tau_s_xi_3d = tau_measured_3d if tau_data_3d else tau_theoretical_3d
for i,p in enumerate(p_vals_3d):
    if p == p_c_3d:
        s_xi_measured_3d[i] = np.nan
        continue
    sizes_3d, densities_3d = cluster_sizes_3D(L_3d, p, nsamp_3d)
    all_dist_3d[p] = (sizes_3d, densities_3d)
    s_xi_3d = s_xi_measure(sizes_3d, densities_3d, tau_s_xi_3d)
    s_xi_measured_3d[i] = s_xi_3d if not np.isnan(s_xi_3d) and s_xi_3d > 0 else 1.0
    
p_data_3d = p_vals_3d
s_xi_data_3d = s_xi_measured_3d
sigma_measured_3d, sigma_data_3d = sigma_measure(p_data_3d, s_xi_data_3d, p_c_3d)
if sigma_data_3d:
    log_delta_p_3d, log_s_xi_3d, slope_3d, intercept_3d = sigma_data_3d
    
# log log plot to measure sigma 2d
    plt.subplot(1,2,1) #2D
#     plt.figure(figsize = (10,6))
    plt.loglog(np.exp(log_delta_p_2d), np.exp(slope_2d * log_delta_p_2d + intercept_2d), 'r--', label=f'Fit: $\\propto |p-p_c|^{{-1/\\sigma}}$ (Sim. $\\sigma={sigma_measured_2d:.3f}$)')
    plt.loglog(np.abs(p_data_2d - p_c_2d), s_xi_data_2d, 'bo', label='$s_\\xi$')
    plt.title(f'Log-log plot to calculate $\\sigma$ in 2D', fontsize = 19)
    plt.xlabel(r'$|p-p_c|$ (log scale)', fontsize = 19)
    plt.ylabel(r'Characteristic Cluster Size (log scale)', fontsize = 19)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 14)
    plt.legend(loc = 'best',fontsize = 19)
    plt.grid(True, which="both", ls="--")
    print(sigma_measured_2d)
    
#  log log plot to measure sigma 3d
    plt.subplot(1,2,2)
#     plt.figure(figsize = (10,6))
    plt.loglog(np.exp(log_delta_p_3d), np.exp(slope_3d * log_delta_p_3d + intercept_3d), 'r--', label=f'Fit: $\\propto |p-p_c|^{{-1/\\sigma}}$ (Sim. $\\sigma={sigma_measured_3d:.3f}$)')
    plt.loglog(np.abs(p_data_3d - p_c_3d), s_xi_data_3d, 'bo', label='$s_\\xi$')
    plt.title(f'Log-log plot to calculate $\\sigma$ in 3D', fontsize = 19)
    plt.xlabel(r'$|p-p_c|$ (log scale)', fontsize = 19)
    plt.ylabel(r'Characteristic Cluster Size (log scale)', fontsize = 19)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 14)
    plt.legend(loc = 'upper center', fontsize = 19)
    plt.grid(True, which="both", ls="--")
    print(sigma_measured_3d)
    plt.show()
    
# 2d
p_col_2d = {0.54: 'mediumblue', 0.57: 'skyblue', p_c_2d: 'magenta', 0.61: 'darkgreen', 0.63: 'gold'}
for p_val_2d in sorted(all_dist_2d.keys()):
    sizes_2d, densities_2d = all_dist_2d[p_val_2d]
    mask_2d = densities_2d > 0
    col_2d = p_col_2d.get(p_val_2d, 'gray')
    marker_2d = 's' if p_val_2d == p_c_2d else 'o'
    plt.loglog(sizes_2d[mask_2d], densities_2d[mask_2d], marker_2d, ms=4, color=col_2d, alpha=0.7, label=f'$p = {p_val_2d:.3f}$')
    sizes_p_val_2d, densities_p_val_2d = cluster_sizes_2D(L_2d, p_val_2d, nsamp_2d)
    all_dist_2d[p_val_2d] = (sizes_p_val_2d, densities_p_val_2d)
    tau_measured_2d, tau_data_2d = tau_measure(sizes_p_val_2d, densities_p_val_2d)
    if tau_data_2d:
        log_s_2d, log_n_s_2d, slope_2d, intercept_2d = tau_data_2d
        plt.loglog(np.exp(log_s_2d), np.exp(slope_2d * log_s_2d + intercept_2d), 'r--', label=f'Fit: $\\propto s^{{-\\tau}}$ (Sim. $\\tau={tau_measured_2d:.3f}$)')
plt.title(f'Cluster Number density $n(s,p)$ for a range of $p$ ({L}x{L} Lattice)',fontsize = 19)
plt.xlabel(r'Cluster size $s$ (log scale)', fontsize = 19)
plt.ylabel(r'Cluster Number density $n(s,p)$ (log scale)', fontsize = 19)
plt.xlim(0,1e3)
plt.legend(loc = 'upper right', title='Occupation probability $p$', fontsize = 15, title_fontsize = 15)
plt.grid(True, which="both", ls="--")
plt.show()

# 3d
p_col_3d = {0.29: 'mediumblue', 0.3: 'skyblue', p_c_3d: 'magenta', 0.32: 'darkgreen', 0.325: 'gold'}
for p_val_3d in sorted(all_dist_3d.keys()):
    sizes_3d, densities_3d = all_dist_3d[p_val_3d]
    mask_3d = densities_3d > 0
    col_3d = p_col_3d.get(p_val_3d, 'gray')
    marker_3d = 's' if p_val_3d == p_c_3d else 'o'
    plt.loglog(sizes_3d[mask_3d], densities_3d[mask_3d], marker_3d, ms=4, color=col_3d, alpha=0.7, label=f'$p = {p_val_3d:.3f}$')
    sizes_p_val_3d, densities_p_val_3d = cluster_sizes_3D(L_3d, p_val_3d, nsamp_3d)
    all_dist_3d[p_val_3d] = (sizes_p_val_3d, densities_p_val_3d)
    tau_measured_3d, tau_data_3d = tau_measure(sizes_p_val_3d, densities_p_val_3d)
    if tau_data_3d:
        log_s_3d, log_n_s_3d, slope_3d, intercept_3d = tau_data_3d
        plt.loglog(np.exp(log_s_3d), np.exp(slope_3d * log_s_3d + intercept_3d), 'r--', label=f'Fit: $\\propto s^{{-\\tau}}$ (Sim. $\\tau={tau_measured_3d:.3f}$)')
plt.title(f'Cluster Number density $n(s,p)$ for a range of $p$ ({L}x{L}x{L} Lattice)',fontsize = 19)
plt.xlabel(r'Cluster size $s$ (log scale)', fontsize = 19)
plt.ylabel(r'Cluster Number density $n(s,p)$ (log scale)', fontsize = 19)
plt.xlim(0,1e3)
plt.legend(loc = 'upper right', title='Occupation probability $p$', fontsize = 15, title_fontsize = 15)
plt.grid(True, which="both", ls="--")
plt.show()
