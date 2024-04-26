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
start_date = "2024-01-01"
end_date = "2024-12-31"
delta_time = "h"  # "Min", "H", 

# Definition of Location object. Coordinates and elevation of Amager, Copenhagen (Denmark)
site = Location(
    latitude=55.78505061217175, longitude=12.51974478835592,
                    tz='Europe/Copenhagen', altitude = 155, name='Dtu 303'
)  # latitude, longitude, time_zone, altitude, name

# Definition of a time range of simulation
times = pd.date_range(
    start_date + " 00:00", end_date + " 00:00", inclusive="left", freq=delta_time, tz=timezone
)
sunpos = site.get_solarposition(times)

theta_sun = np.array(sunpos.zenith)
phi_sun = np.array(sunpos.azimuth)
theta_panel = 0
phi_panel = 90
#print(sunpos.zenith)

def solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel):
    # convert angles from degrees to radians
    theta_sun = np.deg2rad(theta_sun)
    phi_sun = np.deg2rad(phi_sun)
    theta_panel = np.deg2rad(theta_panel)
    phi_panel = np.deg2rad(phi_panel)

    # calculate the unit vectors
    u_s = np.array([np.sin(theta_sun) * np.cos(phi_sun), np.sin(theta_sun) * np.sin(phi_sun), np.cos(theta_sun)])
    u_p = np.array([np.sin(theta_panel) * np.cos(phi_panel), np.sin(theta_panel) * np.sin(phi_panel), np.cos(theta_panel)])

    dot_products=u_p[0]*np.array(u_s[0])+u_p[1]*np.array(u_s[1])+u_p[2]*np.array(u_s[2])
    return dot_products
#print (solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel))

def solar_flux(theta_sun, phi_sun, theta_panel, phi_panel):
    dot_products=solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel)
    #flux_array=np.empty_like(dot_products)
    flux_array = np.array([])
    for element in dot_products:
        if element>0:
            flux_array = np.append(flux_array, 0.5 * element)
        else:
            flux_array = np.append(flux_array, 0)
    return flux_array

np.set_printoptions(threshold=np.inf) #printing the whole array
#print(solar_flux(theta_sun, phi_sun, theta_panel, phi_panel))

time = np.arange(0, 365*24, 1) #time in minutes

#Plot for solar flux over time
fig, ax2 = plt.subplots(figsize=(10, 5))
ax2.plot(time, solar_flux(theta_sun, phi_sun, theta_panel, phi_panel))
ax2.set_ylabel("Solar Flux")
ax2.set_xlabel("Time (hours)")
ax2.set_title("Solar Flux over Time")
#plt.show()


#print(integral_value)

integral_value = integrate.simps(solar_flux(theta_sun, phi_sun, theta_panel, phi_panel),time, dx=3600)

#print(integral_value)
int_values=np.array([])
for i in range(0, 90):
    theta_panel = i
    integral_value = integrate.simps(solar_flux(theta_sun, phi_sun, theta_panel, phi_panel),time, dx=3600)
    int_values = np.append(int_values, integral_value)
#print(int_values)

angles = np.arange(0, 90, 1) #angles
# Plot for integrals over degrees
fig, ax3 = plt.subplots(figsize=(10, 5))
ax3.plot(angles, int_values)
ax3.set_ylabel("Energy production (Wh/m^2)")
ax3.set_xlabel("Theta_p angle (degrees)")
ax3.set_title("Yearly energy prodcution over theta_p angle")

plt.show()

max_prod=np.max(int_values)
optimal_angle=np.argmax(int_values)
print("The optimal angle is: ", optimal_angle, "degrees")
print(max_prod)