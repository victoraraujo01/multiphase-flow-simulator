import numpy as np
import matplotlib.pyplot as plt

from oil import Oil
from water import Water
from gas import Gas
from mixture import Mixture
from beggs_and_brill import BeggsAndBrill
from tubing import Tubing
from ipr import IPR
import helpers

well_head_pressure = 150.0
temperature = 170.0

prod_glr = 600.0
water_cut = 0.0

x = []
y = []
x_ipr = []
y_ipr = []

tubing = Tubing(10000, 1.995, 90.0, 0.00015)
gas = Gas(0.7)
oil = Oil(25.0, gas)
water = Water(1.07, gas)
mixture = Mixture(gas, oil, water, water_cut, prod_glr, 0.)
bubble_point = helpers.bubble_point(temperature, mixture)

for flow_rate in np.arange(1, 6500, 10):
    tubing = Tubing(10000, 1.995, 90.0, 0.00015)
    gas = Gas(0.7)
    oil = Oil(25.0, gas)
    water = Water(1.07, gas)
    mixture = Mixture(gas, oil, water, water_cut, prod_glr, flow_rate)
    correlation = BeggsAndBrill()

    last_pressure = 0.0
    pressure = well_head_pressure
    last_depth = 0.0
    depth = 0.0

    delta_p = 0.5
    last_limit = 5.0
    limit = 10.0
    while depth < tubing.length:
        gas.update_conditions(pressure, temperature)
        oil.update_conditions(pressure, bubble_point, temperature, mixture.water_cut)
        water.update_conditions(pressure, bubble_point, temperature, mixture.water_cut)
        mixture.update_conditions(pressure, bubble_point, tubing.diameter)
        correlation.update_conditions(mixture, tubing)

        dp_dz = correlation.grav_pressure_gradient + correlation.fric_pressure_gradient
        last_depth = depth
        depth = depth + delta_p / abs(dp_dz)

        last_pressure = pressure
        pressure = pressure + delta_p

        if pressure >= limit:
            delta_p = limit / 10.0
            tmp_limit = limit
            limit = last_limit * 10.0
            last_limit = tmp_limit

    factor = (tubing.length - last_depth) / (depth - last_depth)
    pressure = last_pressure + factor * (pressure - last_pressure)

    x.append(flow_rate)
    y.append(pressure)

ipr = IPR.from_tests(-0.2, bubble_point, (5000., 1000.), (2000., 5000.))
x_ipr = x
y_ipr = [ipr.pressure(flow_rate) for flow_rate in x_ipr]
# y_ipr = range(0, 6000)
# x_ipr = [ipr.flow_rate(pressure) for pressure in y_ipr]

plt.plot(x, y, x_ipr, y_ipr)
plt.grid(True)
plt.xlabel('Flow rate')
plt.ylabel('Pressure')
plt.show()
