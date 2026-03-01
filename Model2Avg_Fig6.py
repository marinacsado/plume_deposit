"""
@author: marina

This is a python script for Figure 6 in Casado-Anarte et al. 2025.

Calculates the rosion time of plume’s deposit as a function of mass flux and time of eruption
25 km away from the source for an average case.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Constants
A = 10**8  # Area of the plume source (m^2)
v_avg = 460  # m/s
s_r = 2.5e14   # Sputtering rate (ptcls/sm²)
radrate_O2 = 6.3e13  # Radiolysis rate (O2/sm²)
radrate_H2 = 8e13  # Radiolysis rate (H2/m²s )
density_ppcc = 3.3e10/7000    # Ptc/cm³ for a plume of m = 1 kg/s
H2O = 2.987e-26  # Kg per particle of H2O
O2 = 5.3137E-26  # Kg per particle of O2
H2 = 3.3543E-27  # Kg per particle of H2
distance_conversion = 10**6  # cubic meters per 1 cm³


# Converting the constants
radiolysis_O2 = radrate_O2 * O2  # Radiolysis conversion kg/sm²
radiolysis_H2 = radrate_H2 * H2  # Radiolysis conversion kg /m²
sputtering_rate = s_r * H2O  # Sputtering conversion (kg/sm²)

erosion_factors = radiolysis_H2 + radiolysis_O2 + sputtering_rate

rho = density_ppcc * H2O * distance_conversion


seconds_per_month = 2.628e6  # Conversion factor: seconds per month
seconds_per_year = 3.154e+7  # Conversion factor: seconds per year

# Define the ranges for eruption time and mass flux
eruption_times = np.logspace(np.log10(262.800288),
                             np.log10(3.154e+7), 1000)  # Eruption times from less than 1-12 months
mass_flux_rates = np.logspace(0, 4, 1000)  # Mass flux rates from 1-10,000 kg/s


# CALCULATIONS:
# Function to calculate the density based on mass flux rates
def calculate_density(rho, mass_flux_rates):
    rho_step = rho * mass_flux_rates
    return rho_step


# Calculate the mass flux above the surface of the falling particles
def calculate_mass_flux_falling(rho_step, A, v_avg):
    massflux_f = rho_step * A * v_avg
    return massflux_f


# Calculate total particles erupted
def calculate_total_particles_eruption(massflux_f, t_eruption):
    N = massflux_f * t_eruption  # Convert time to seconds
    return N


# Calculate maximum disappearance time
def calculate_max_disappearance_time(total_particles, erosion_factors):
    time_in_seconds = total_particles / (erosion_factors * A)
    Tmax = time_in_seconds / seconds_per_year  # Convert time to years
    return Tmax


# Create a meshgrid for plotting
mass_flux_grid, eruption_time_grid = np.meshgrid(mass_flux_rates,
                                                 eruption_times)

# Initialize the array for maximum disappearance times
Tmax_grid = np.zeros_like(mass_flux_grid, dtype=float)

# Calculate Tmax for each combination of mass flux rate and eruption time
for i in range(len(eruption_times)):
    for j in range(len(mass_flux_rates)):
        rho_step = calculate_density(rho, mass_flux_rates[j])
        massflux_f = calculate_mass_flux_falling(rho_step, A, v_avg)
        total_particles = calculate_total_particles_eruption(massflux_f,
                                                             eruption_times[i])
        Tmax_grid[i, j] = calculate_max_disappearance_time(total_particles,
                                                           erosion_factors)



# Labels added manually
plt.figure(figsize=(8,  6))
contour = plt.contour(mass_flux_grid, eruption_time_grid, Tmax_grid,
                      levels=[10**i for i in range(-5, 6)],
                      colors='black', linewidths=0.5)
plt.text(3, 2000, '1.0e-02 years', fontsize=8,
         color='black',  ha='left',  va='center')
plt.text(35, 17000, '1.0e-00 years', fontsize=8,
         color='black',  ha='left',  va='center')
plt.text(550, 130000, '1.0e+02years', fontsize=8,
         color='black',  ha='left',  va='center')
plt.text(1300, 5e6, '1.0e+04 years', fontsize=8,
         color='black',  ha='left',  va='center')
plt.pcolormesh(mass_flux_grid, eruption_time_grid,
               Tmax_grid, shading='auto', norm=LogNorm(), cmap='viridis')
plt.colorbar(label='Erosion Time (years)')

# Highlight the specific point as this is the plume that was observed
example_mass_flux_rate = 7000  # Example mass flux rate (in kg/s)
example_eruption_time = 25200  # Eruption of 7h (converted to months)
plt.plot(example_mass_flux_rate, example_eruption_time, 'o', color='#D41159')
plt.text(example_mass_flux_rate * 0.4, example_eruption_time * 0.8,
         '(7000 kg/s, ~7h)', color='black', ha='center', va='top')


# Reference lines in the y-axis text (commented out for this example)

plt.axhline(y=3600, color='grey',
            linestyle='-', label='1 h')  # 1 h reference line
plt.text(0.8, 3600, '1 h', color='grey', ha='right', va='bottom')

plt.axhline(y=86400, color='grey',
            linestyle='-', label='1 day')  # 1 day reference line
plt.text(0.8, 60000, '1 day', color='grey', ha='right', va='bottom')

plt.axhline(y=2.628e+6, color='grey',
            linestyle='-', label='1 month')  # 1 month reference line
plt.text(0.8, 2.628e+6, '1 month', color='grey', ha='right', va='bottom')

plt.axhline(y=1.577e+7, color='grey',
            linestyle='-', label='6 months')  # 6 months reference line
plt.text(0.8, 1.577e+7, '6 months', color='grey', ha='right', va='bottom')

# Add the contour line for Tmax = 28 years
c28 = plt.contour(mass_flux_grid, eruption_time_grid,
                  Tmax_grid, levels=[28], colors='white', linewidths=1.5)
plt.text(35, 5e5, '28 years', color='white', ha='left', va='center')

# Plotting
plt.xlabel('Mass Flux Rate (kg/s)')
plt.ylabel('Eruption Time (seconds)', labelpad=15)
plt.yscale('log') 
plt.xscale('log') 
plt.title('Erosion Time of Deposit 25 km Away from Plume Source ')
plt.show()

# Example calculations for print
rho_step_example = calculate_density(rho, example_mass_flux_rate)
massflux_f_example = calculate_mass_flux_falling(rho_step_example, A, v_avg)
total_particles = calculate_total_particles_eruption(massflux_f_example,
                                                     example_eruption_time)
max_disappearance_time = calculate_max_disappearance_time(total_particles,
                                                          erosion_factors)
print(f"Density for given mass flux:  {rho_step_example:.2E} kg/m³")
print(f"Mass flux falling to surface : {massflux_f_example:.2E} kg/s")
print(f"Total particles erupted: {total_particles:.2E} kg")
print(f"Max disappearance time: {max_disappearance_time:.2E} years")
