import pandas as pd
import pvlib
from pvlib.location import Location
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from scipy.signal import argrelextrema


timezone = "Europe/Copenhagen"
start_date = "2024-04-01"
end_date = "2024-06-30"
delta_time = "Min"  # "Min", "H", 


# latitude, longitude, time_zone, altitude, name
site = Location(latitude=55.78505061217175, longitude=12.51974478835592,
                    tz='Europe/Copenhagen', altitude = 155, name='Dtu 303') #Building 303

# Definition of a time range of simulation
times = pd.date_range(
    start_date + " 00:00:00", end_date + " 23:59:00", inclusive="left", freq=delta_time, tz=timezone
)

sunpos = site.get_solarposition(times)

chosen_date = "2024-04-22"

# Plots for solar zenith and solar azimuth angles
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(30, 10))
fig.suptitle("Solar Position Estimation in " + site.name + " " + chosen_date)

# plot for solar zenith angle
ax1.plot(sunpos.loc[chosen_date].zenith)
ax1.set_ylabel("Solar zenith angle (degree)")
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

# plot for solar azimuth angle
ax2.plot(sunpos.loc[chosen_date].azimuth)
ax2.set_ylabel("Solar azimuth angle (degree)")
ax2.set_xlabel("Time UTC+2 (hour)")
# Get the current x-axis tick positions and labels
ticks = ax2.get_xticks()
tick_labels = ax2.get_xticklabels()
# Convert the tick positions to datetime objects
tick_datetimes = [mdates.num2date(tick) for tick in ticks]
# Add 2 hours to each datetime object
new_tick_datetimes = [tick_datetime + datetime.timedelta(hours=2) for tick_datetime in tick_datetimes]
# Convert the updated datetime objects back to numerical values
new_ticks = mdates.date2num(new_tick_datetimes)
# Set the new tick positions and labels
ax2.set_xticks(ticks)
ax2.set_xticklabels([tick_label.strftime('%H') for tick_label in new_tick_datetimes])

# plot for solar elevation angle
ax3.plot(sunpos.loc[chosen_date].elevation)
ax3.set_ylabel("Solar elevation angle in orange (degree)")
ax3.set_xlabel("Time UTC+2 (hour)")
# Get the current x-axis tick positions and labels
ticks = ax3.get_xticks()
tick_labels = ax3.get_xticklabels()
# Convert the tick positions to datetime objects
tick_datetimes = [mdates.num2date(tick) for tick in ticks]
# Add 2 hours to each datetime object
new_tick_datetimes = [tick_datetime + datetime.timedelta(hours=2) for tick_datetime in tick_datetimes]
# Convert the updated datetime objects back to numerical values
new_ticks = mdates.date2num(new_tick_datetimes)
# Set the new tick positions and labels
ax3.set_xticks(ticks)
ax3.set_xticklabels([tick_label.strftime('%H') for tick_label in new_tick_datetimes])



# plt.show()

# def solar_elevation_angle(theta):
#     return 90 - theta

# chosen_date = "2024-04-20"
# print(sunpos.loc[chosen_date].zenith)
# print(sunpos.loc[chosen_date].elevation)
# print(sunpos.loc[chosen_date].azimuth)

# FINDING THE HIGHEST POINT OF SUN

chosen_date = "2024-06-20"
elevation_vector=np.array(sunpos.loc[chosen_date].elevation)
highest_angle=elevation_vector.max()

y_value = highest_angle  # The y value you want to find the corresponding x value for
x_values = sunpos.loc[chosen_date].index  # The x values from the plot
# Find the index of the x value where the y value is closest to the desired value
index = np.where(np.isclose(sunpos.loc[chosen_date].elevation, y_value))[0][0]
# Get the corresponding x value
x_value = x_values[index]
print("x value:", x_value)
print("y value:", y_value)


#FINDING SUNSET AND SUNRISE TIMES

chosen_date = "2024-04-20"
series = np.abs(sunpos.loc[chosen_date].apparent_elevation)

# Convert pandas Series to NumPy array
arr = series.to_numpy()
# Find indices of local minima
minima_indices = argrelextrema(arr, np.less)[0]
# Get values at local minima
minima_values = arr[minima_indices]

# print("Indices of local minima:", minima_indices)
# print("Values of local minima:", minima_values)
x_values = sunpos.loc[chosen_date].index  # The x values from the plot
for value in minima_indices:
    if sunpos.loc[chosen_date].apparent_elevation[value]<sunpos.loc[chosen_date].apparent_elevation[value+1]:
        sunrise = sunpos.loc[chosen_date].apparent_elevation[value]
        sunrise_i=value
        sunrise_value = x_values[sunrise_i]
        print("Sunrise time: ", sunrise_value)
    else:
        sunset = sunpos.loc[chosen_date].apparent_elevation[value]
        sunset_i=value
        sunset_value = x_values[sunset_i]
        print("Sunset time: ", sunset_value)

#FUNCTION FOR FINDING THE HIGHEST POINT FOR THE SUN

# import pytz
# from tzwhere import tzwhere
# #import elevation

# def get_timezone(lat, lon):
#     tz = tzwhere.tzwhere()
#     timezone_str = tz.tzNameAt(lat, lon)
#     return pytz.timezone(timezone_str)

def highest_point_of_the_sun(date, latitude, longitude):
    delta_time = "Min"  # "Min", "H", 
    #timezone = get_timezone(latitude, longitude)
    #altitude = elevation.ggm.load(latitude, longitude)
    site = Location(latitude, longitude)
    times = pd.date_range(
        date + " 00:00:00", date + " 23:59:00", inclusive="left", freq=delta_time
    )
    sunpos = site.get_solarposition(times)
    elevation_vector = np.array(sunpos.loc[date].elevation)
    highest_angle = elevation_vector.max()
    y_value = highest_angle  # The y value you want to find the corresponding x value for
    x_values = sunpos.loc[date].index  # The x values from the plot
    # Find the index of the x value where the y value is closest to the desired value
    index = np.where(np.isclose(sunpos.loc[date].elevation, y_value))[0][0]
    # Get the corresponding x value
    x_value = x_values[index]
    #return (str(x_value) + ' UTC')
    return (highest_angle)

#print(highest_point_of_the_sun("2024-04-20", 55.78505061217175, 12.51974478835592))

#FUNCTION CHANGING FROM ANGLES TO XYZ

# def angles_to_xyz(date, latitude, longitude):
#     chosen_date=date
#     site = Location(latitude, longitude)
#     sunpos = site.get_solarposition(times)
#     theta_p = sunpos.loc[chosen_date].zenith
#     phi_p = sunpos.loc[chosen_date].azimuth

def angles_to_xyz(zenith, azimuth):
    #r=pvlib.solarposition.nrel_earthsun_distance(times) * 149597870700
    # The unit normal vector to the solar panel expressed in terms of zenith and azimuth angles
    zenith=np.deg2rad(zenith)
    azimuth=np.deg2rad(azimuth)
    x = np.sin(zenith) * np.cos(azimuth)
    y = np.sin(zenith) * np.sin(azimuth)
    z = np.cos(zenith)
    return (x, y, z)

#print(angles_to_xyz(sunpos.loc[chosen_date].zenith[1], sunpos.loc[chosen_date].azimuth[1]))

#FUNCTION CHANGING FROM XYZ TO ANGLES

def xyz_to_angles(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    zenith = np.arccos(z / r)
    zenith = np.rad2deg(zenith)
    azimuth = np.arctan2(y, x)
    azimuth = np.rad2deg(azimuth)
    return (zenith, azimuth)

#print(xyz_to_angles(1, -4, -5))

