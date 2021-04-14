###########################################################################
# GAUSSIAN PLUME MODEL FOR TEACHING PURPOSES                              #
# PAUL CONNOLLY (UNIVERSITY OF MANCHESTER, 2017)                          #
# MODIFIED BY AIDEN KERR (UNIVERSITY OF BRISTOL, 2021)                    #
# USED TO ESTIMATE THE IMPACTS OF METHANE POINT SOURCES ON                #
# ATMOSPHERIC POLLUTANT LEVELS USING REAL-LIFE METEOROLOGICAL DATA        #
###########################################################################

import numpy as np
import pandas as pd
import sys
from scipy.special import erfcinv as erfcinv
import tqdm as tqdm

from gauss_func import gauss_func

import matplotlib.pyplot as plt
from matplotlib import rc


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


###########################################################################
# Do not change these variables                                           #
###########################################################################


# SECTION 0: Definitions (normally don't modify this section)
# view
PLAN_VIEW = 1
HEIGHT_SLICE = 2
SURFACE_TIME = 3
NO_PLOT = 4
BOTH = 5

# wind field
CONSTANT_WIND = 1
FLUCTUATING_WIND = 2
PREVAILING_WIND = 3
MET_WIND = 4

# number of stacks
ONE_STACK = 1
TWO_STACKS = 2
THREE_STACKS = 3

# stability of the atmosphere
CONSTANT_STABILITY = 1
ANNUAL_CYCLE = 2
stability_str = ['Very unstable', 'Moderately unstable', 'Slightly unstable',
                 'Neutral', 'Moderately stable', 'Very stable']
# Aerosol properties
HUMIDIFY = 2
DRY_AEROSOL = 1

SODIUM_CHLORIDE = 1
SULPHURIC_ACID = 2
ORGANIC_ACID = 3
AMMONIUM_NITRATE = 4
METHANE = 5
nu = [2., 2.5, 1., 2., 0]  # hydration number? - does methane hydrate in atmosphere? i think no
rho_s = [2160., 1840., 1500., 1725., 657.]  # density - for methane, 657 g/m3
Ms = [58.44e-3, 98e-3, 200e-3, 80e-3, 16e-3]  # mass of solvent (kg / mol)
Mw = 18e-3  # mass of water

dxy = 50  # resolution of the model in both x and y directions (m) - 100 standard!
dz = 10
radius = 10000
x = np.mgrid[-radius:radius + dxy:dxy]  # solve on a x by x km domain
y = x  # x-grid is same as y-grid
###########################################################################


# SECTION 1: Configuration
# Variables can be changed by the user+++++++++++++++++++++++++++++++++++++
RH = 0.90 # relative humididty!
aerosol_type = METHANE

dry_size = 60e-9
humidify = DRY_AEROSOL

stab1 = 5  # set from 1-6
stability_used = CONSTANT_STABILITY

output = BOTH  # for some reason, HEIGHT_SLICE isn't working - due to index out of range error
x_slice = 200  # position to take the slice in the x-direction
y_slice = 200  # position to plot concentrations vs time

wind = MET_WIND
wind_data = 'hourly_met_data.csv'
wind_spd_column = 'Wind Speed  Avg(m/s)'

wind_speed_input = 1
#stacks = 'NAEI SOURCES'
NAEI_data = 'Close_Sources.csv'
stacks = THREE_STACKS
stack_x = [-4894., 2536., -6294.]  # x pos. of each of the stacks
stack_y = [9306., 7806., 7306.]

map_img = plt.imread('Bristol10kmradiusmap.png')

Q = [0.747, 0.381, 0.0964]  # mass emitted per unit time (g / s) (23.57 tonnes per year = 0.747 g/s)
H = [10., 10., 10.]  # stack height, m
days = 5  # run the model for n days
start_date = '2019-01-01'  # for use of wind data timeseries
# --------------------------------------------------------------------------
times = np.mgrid[1:days * 24 + 1:1] / 24.

Dy = 10.
Dz = 10.

# SECTION 2: Act on the configuration information
if stacks == 'NAEI SOURCES':
    centre_loc = [358394, 173194]
    sources = pd.read_csv(NAEI_data)
    stack_x = sources['Easting'] - centre_loc[0]
    stack_y = sources['Northing'] - centre_loc[1]
    Q = sources['Emission'] * 1000000 / (365 * 24 * 60 * 60)  # converts emissions data to g/s
    print(str(Q))
    stacks = len(stack_x)
    print('Stacks: ' + str(stacks))
    H = [10.0] * stacks
    print(str(H))

# Decide which stability profile to use
if stability_used == CONSTANT_STABILITY:

    stability = stab1 * np.ones((days * 24, 1))
    stability_str = stability_str[stab1 - 1]
elif stability_used == ANNUAL_CYCLE:

    stability = np.round(2.5 * np.cos(times * 2. * np.pi / (365.)) + 3.5)
    stability_str = 'Annual cycle'
else:
    sys.exit()

# decide what kind of run to do, plan view or y-z slice, or time series
if output == PLAN_VIEW or output == SURFACE_TIME or output == NO_PLOT or output == BOTH:

    C1 = np.zeros((len(x), len(y), days * 24))  # array to store data, initialised to be zero

    [x, y] = np.meshgrid(x, y)  # x and y defined at all positions on the grid
    z = np.zeros(np.shape(x))  # z is defined to be at ground level.
elif output == HEIGHT_SLICE:
    z = np.mgrid[0:500 + dz:dz]  # z-grid

    C1 = np.zeros((len(y), len(z), days * 24))  # array to store data, initialised to be zero

    [y, z] = np.meshgrid(y, z)  # y and z defined at all positions on the grid
    x = x[x_slice] * np.ones(np.shape(y))  # x is defined to be x at x_slice
else:
    sys.exit()

# Set the wind based on input flags++++++++++++++++++++++++++++++++++++++++
wind_speed = wind_speed_input * np.ones((days * 24, 1))  # m/s - wind speed array (hourly)
if wind == CONSTANT_WIND:
    wind_dir = 0. * np.ones((days * 24, 1))
    wind_dir_str = 'Constant wind'
elif wind == FLUCTUATING_WIND:
    wind_dir = 360. * np.random.rand(days * 24, 1)
    wind_dir_str = 'Random wind'
elif wind == PREVAILING_WIND:
    wind_dir = -np.sqrt(2.) * erfcinv(2. * np.random.rand(24 * days, 1)) * 40.  # norminv(rand(days.*24,1),0,40)
    # note at this point you can add on the prevailing wind direction, i.e.
    # wind_dir=wind_dir+200
    wind_dir[np.where(wind_dir >= 360.)] = \
        np.mod(wind_dir[np.where(wind_dir >= 360)], 360)
    wind_dir_str = 'Prevailing wind'
elif wind == MET_WIND:
    print('Using MET data')
    end_date = pd.Timedelta(days - 1, unit='d') + pd.to_datetime(start_date)
    end_date = end_date.strftime('%Y-%m-%d')
    print('Start: ' + start_date + ' End: ' + end_date)
    met_data = pd.read_csv(wind_data, index_col='datetime', parse_dates=True)
    met_data = met_data.loc[start_date:end_date]
    print(met_data.to_string())
    wind_speed = met_data[wind_spd_column]
    wind_speed[wind_speed < 0.5] = None  # Replace all values lower than 2 as None
    wind_speed = wind_speed.to_list()
    print('Wind speed:' + str(wind_speed))
    wind_dir = met_data['Wind Direction (degrees)'].to_list()
    print('Wind dir:' + str(wind_dir))
    wind_dir_str = 'MET data wind'
else:
    sys.exit()
# --------------------------------------------------------------------------


# SECTION 3: Main loop
# For all times...
C1 = np.zeros((len(x), len(y), len(wind_dir)))  # Fill colormap with 1s
for i in tqdm.tqdm(range(0, len(wind_dir))):  # For each iteration (loading bar)
    for j in range(0, stacks):  # for each stack
        C = np.ones((len(x), len(y)))
        C = gauss_func(Q[j], wind_speed[i], wind_dir[i], x, y, z,
                       stack_x[j], stack_y[j], H[j], Dy, Dz, stability[i])
        C1[:, :, i] = C1[:, :, i] + C

# SECTION 4: Post process / output

# decide whether to humidify the aerosol and hence increase the mass
if humidify == DRY_AEROSOL:
    print('do not humidify')
elif humidify == HUMIDIFY:
    mass = np.pi / 6. * rho_s[aerosol_type] * dry_size ** 3.
    moles = mass / Ms[aerosol_type]

    nw = RH * nu[aerosol_type] * moles / (1. - RH)
    mass2 = nw * Mw + moles * Ms[aerosol_type - 1]
    C1 = C1 * mass2 / mass
else:
    sys.exit()

# convert C1 (color map from quantities) from g /cm3 to ppb
temp = 15 + 273  # Temp in K (temp in C + 273)
mol_volume = 22.41 * (temp / 273)  # from http://www.apis.ac.uk/unit-conversion
mol_mass_gcm3 = Ms[aerosol_type - 1] * 1000
ppb = C1 * (mol_volume / mol_mass_gcm3)


# output the plots
if wind == MET_WIND:
    title = 'Gaussian Modelling starting from ' + start_date + ' for ' + str(days) + ' days'
else:
    title = stability_str + '\n' + wind_dir_str + ', ' + str(wind_speed_input) + 'm/s' + ' for ' + str(days) + 'days'

if output == PLAN_VIEW or output == BOTH:

    BBox = (-radius, radius, -radius, radius)
    plt.figure()
    plt.ion()
    plt.imshow(map_img, zorder=0, extent=BBox, aspect='equal')
    plt.pcolor(x, y, np.mean(ppb, axis=2) * 1e6, cmap='jet', shading='auto', alpha=0.15)
    # plt.scatter(x_slice, y_slice, c='black', s=20)  # delete!!
    plt.clim((0, 10))
    plt.title(title)
    plt.xlabel('x (metres)')
    plt.ylabel('y (metres)')
    cb1 = plt.colorbar()
    cb1.set_label('ppb')
    plt.savefig('Plan_view_plot')
    plt.show()
    print('Plotted Plan View')

if output == HEIGHT_SLICE:
    plt.figure()
    plt.ion()

    plt.pcolor(y, z, np.mean(C1, axis=2) * 1e6, cmap='jet')
    plt.clim((0, 1e2))
    plt.xlabel('y (metres)')
    plt.ylabel('z (metres)')
    plt.title(title)
    cb1 = plt.colorbar()
    cb1.set_label('$\mu$ g m$^{-3}$')
    plt.savefig('Height_slice_plot')
    plt.show()

if output == SURFACE_TIME or output == BOTH:
    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
    ax1.plot(times, 1e6 * np.squeeze(ppb[y_slice, x_slice, :]))
    try:
        ax1.plot(times, smooth(1e6 * np.squeeze(ppb[y_slice, x_slice, :]), 24), 'r')
        ax1.legend(('Hourly mean', 'Daily mean'))
    except:
        print('Problem, cannot plot')
        sys.exit()

    ax1.set_xlabel('time (days)')
    plt.xticks()
    ax1.set_ylabel('ppb')
    ax1.set_title(title)

    ax2.plot(times, stability)
    ax2.set_xlabel('time (days)')
    ax2.set_ylabel('Stability parameter')
    f.savefig('Surface_time_plot')
    f.show()
    print('Surface Time plotted')

if output == NO_PLOT:
    print('don''t plot')




