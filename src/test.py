import numpy as np
import matplotlib.pyplot as plt

from oil import Oil
from water import Water
from gas import Gas
from flow import Flow
from tubing import Tubing
import helpers

well_head_pressure = 5.0
temperature = 175.0
delta_l = 1

x = []
y = []

for flow_rate in np.arange(1, 2000, 100):
    tubing = Tubing(8000, 1.995, 90.0, 0.000902255639097744)
    gas = Gas(0.65)
    oil = Oil(25.0, gas)
    water = Water(1.07, gas)
    flow = Flow(flow_rate, 10.0, 0, tubing)

    bubble_point = helpers.bubble_point(temperature, oil, water, flow)
    pressure = well_head_pressure
    for l in np.arange(0, tubing.length, delta_l):
        gas.update_conditions(pressure, temperature)
        oil.update_conditions(pressure, bubble_point, temperature, flow.water_cut)
        water.update_conditions(pressure, bubble_point, temperature, flow.water_cut)
        flow.update_conditions(pressure, bubble_point, gas, oil, water)

        delta_p = flow.grav_pressure_gradient + flow.fric_pressure_gradient
        pressure = pressure - delta_p

    x.append(flow_rate)
    y.append(pressure)
    print(flow_rate)
    print(pressure)

plt.plot(x, y)
plt.grid(True)
plt.xlabel('Tubing length')
plt.ylabel('Pressure')
plt.show()
