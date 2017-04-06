import numpy as np
import matplotlib.pyplot as plt

from oil import Oil
from water import Water
from gas import Gas
from flow import Flow
from tubing import Tubing
import helpers

well_head_pressure = 150.0
temperature = 170.0

x = []
y = []

for flow_rate in np.arange(1, 1000, 10):
    print(flow_rate)
    tubing = Tubing(10000, 1.995, 90.0, 0.00015)
    gas = Gas(0.7)
    oil = Oil(25.0, gas)
    water = Water(1.07, gas)
    flow = Flow(flow_rate, 600.0, 0, tubing)

    bubble_point = helpers.bubble_point(temperature, oil, water, flow)

    last_pressure = 0.0
    pressure = well_head_pressure
    last_depth = 0.0
    depth = 0.0

    delta_p = 0.5
    last_limit = 5.0
    limit = 10.0
    while depth < tubing.length:
        gas.update_conditions(pressure, temperature)
        oil.update_conditions(pressure, bubble_point, temperature, flow.water_cut)
        water.update_conditions(pressure, bubble_point, temperature, flow.water_cut)
        flow.update_conditions(pressure, bubble_point, gas, oil, water)

        dp_dz = flow.grav_pressure_gradient + flow.fric_pressure_gradient
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

plt.plot(x, y)
plt.grid(True)
plt.xlabel('Tubing length')
plt.ylabel('Pressure')
plt.show()
