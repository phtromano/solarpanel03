import pandas as pd
import pvlib
from pvlib.location import Location
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import integrate
import datetime


timezone = "Europe/Copenhagen"
start_date = "2024-04-01"
end_date = "2024-04-30"
delta_time = "Min"  # "Min", "H", 

# Definition of Location object. Coordinates and elevation of Amager, Copenhagen (Denmark)
site = Location(
    latitude=55.78505061217175, longitude=12.51974478835592,
                    tz='Europe/Copenhagen', altitude = 155, name='Dtu 303'
)  # latitude, longitude, time_zone, altitude, name

# Definition of a time range of simulation
times = pd.date_range(
    start_date + " 00:00:00", end_date + " 23:59:00", inclusive="left", freq=delta_time, tz=timezone
)
#closed="left", - missing from date range
sunpos = site.get_solarposition(times)

chosen_date = "2024-04-20"
#print(sunpos.loc[chosen_date].zenith)


theta_sun = np.array(sunpos.loc[chosen_date].zenith)
phi_sun = np.array(sunpos.loc[chosen_date].azimuth)
array_size=np.shape(theta_sun)
theta_panel = 0
phi_panel = 180


def solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel):
    # convert angles from degrees to radians
    theta_sun = np.deg2rad(theta_sun)
    phi_sun = np.deg2rad(phi_sun)
    theta_panel = np.deg2rad(theta_panel)
    phi_panel = np.deg2rad(phi_panel)

    # calculate the unit vectors
    u_s = np.array([np.sin(theta_sun) * np.cos(phi_sun), np.sin(theta_sun) * np.sin(phi_sun), np.cos(theta_sun)])
    u_p = np.array([np.sin(theta_panel) * np.cos(phi_panel), np.sin(theta_panel) * np.sin(phi_panel), np.cos(theta_panel)])

    # calculate the dot product
    dot_products=u_p[0]*np.array(u_s[0])+u_p[1]*np.array(u_s[1])+u_p[2]*np.array(u_s[2])
    return dot_products

#print(solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel))

#CALCULATING POWER BY SOLAR FLUX
def solar_flux(theta_sun, phi_sun, theta_panel, phi_panel):
    dot_products=solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel)
    #flux_array=np.empty_like(dot_products)
    flux_array = np.array([])
    for element in dot_products:
        if element>0:
            flux_array = np.append(flux_array, 0.5 * element) #per squared meter
            #flux_array = 0.5 * element #from flux calculations
        else:
            flux_array = np.append(flux_array, 0)
    return flux_array

np.set_printoptions(threshold=np.inf) #printing the whole array
#print(solar_flux(theta_sun, phi_sun, theta_panel, phi_panel))

time = np.arange(0, 24*60, 1) #time in minutes

# Plots for solar zenith and solar azimuth angles
fig, (ax1) = plt.subplots(1, 1, figsize=(30, 10))
fig.suptitle("Solar Energy production Estimation in " + site.name + " " + chosen_date)

# plot for solar zenith angle
ax1.plot(solar_flux(theta_sun, phi_sun, theta_panel, phi_panel))
ax1.set_ylabel("Energy")
ax1.set_xlabel("Time UTC+2 (hour)")
# Get the current x-axis tick positions and labels
ticks = ax1.get_xticks()
tick_labels = ax1.get_xticklabels()
# Convert the tick positions to datetime objects
tick_datetimes = [mdates.num2date(tick) for tick in ticks]
# Add 2 hours to each datetime object
new_tick_datetimes = [tick_datetime + datetime.timedelta(hours=2) for tick_datetime in tick_datetimes]
# Convert the updated datetime objects back to numerical values
new_ticks = mdates.date2num(new_tick_datetimes)
# Set the new tick positions and labels
ax1.set_xticks(ticks)
ax1.set_xticklabels([tick_label.strftime('%H') for tick_label in new_tick_datetimes])

# PLOT FOR SOLAR FLUX OVER TIME
fig, ax2 = plt.subplots(figsize=(10, 5))
ax2.plot(time, solar_flux(theta_sun, phi_sun, theta_panel, phi_panel))
ax2.set_ylabel("Solar Flux (W/m^2)")
ax2.set_xlabel("UTC Time (minutes)")
ax2.set_title("Solar Flux over Time on 20th of April 2024 DTU 303")
#plt.show()

#CALCULATING ENERGY PRODUCTION
integral_value = integrate.simps(solar_flux(theta_sun, phi_sun, theta_panel, phi_panel),time, dx=60)

#print(integral_value)
int_values=np.array([])
#dot_products=np.array([])
for i in range(0, 90):
    theta_panel = i
    # dot_product=solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel)
    # dot_integrals_value = integrate.simps(dot_products,time, dx=60)
    # dot_products = np.append(dot_products, dot_integrals_value)
    integral_value = integrate.simps(solar_flux(theta_sun, phi_sun, theta_panel, phi_panel),time, dx=60)
    int_values = np.append(int_values, integral_value)


#PLOT ENERGY PRODUCTION OVER ANGLES
#print(int_values)
#print(dot_products)
angles = np.arange(0, 90, 1) #angles
# Plot for integrals over degrees
fig, ax3 = plt.subplots(figsize=(10, 5))
ax3.plot(angles, int_values)
ax3.set_ylabel("Energy production (Wh/m^2)")
ax3.set_xlabel("Theta_p angle (degrees)")
ax3.set_title("Daily energy production over theta_p angle")

plt.show()