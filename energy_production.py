import pandas as pd
from pvlib.location import Location
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import integrate

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
# function to calculate the solar panel projection
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

def solar_flux(theta_sun, phi_sun, theta_panel, phi_panel):
    dot_products=solar_panel_projection(theta_sun, phi_sun, theta_panel, phi_panel)
    flux_array = np.array([])
    for element in dot_products:
        if element>0:
            flux_array = np.append(flux_array, 0.5 * element)
        else:
            flux_array = np.append(flux_array, 0)
    return flux_array

number_of_panels=20
panels_area=2.16
total_area=number_of_panels*panels_area
theta_panel = 55
phi_panel = 180
efficiency = 0.215
delta_time = "H" 
time = pd.date_range(start_date + " 00:00:00", end_date + " 00:00:00", inclusive="left", freq=delta_time)
sunpos = site.get_solarposition(time)
factor = efficiency * total_area
#calculate the energy production for each day
energy_production = []
for day in times.date.tolist():
    chosen_date = str(day)
    daily_energy = factor * integrate.simps(solar_flux(sunpos.loc[chosen_date].zenith, sunpos.loc[chosen_date].azimuth, theta_panel, phi_panel), np.arange(0, 24, 1), dx=3600)
    print(daily_energy)
    energy_production.append(daily_energy)
#plot the energy production
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(times.date.tolist(), energy_production)
ax.set_ylabel("Energy Production (kWh)")
ax.set_xlabel("Date")
ax.set_title("Energy Production per Day")
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
plt.xticks(rotation=45)
plt.show()