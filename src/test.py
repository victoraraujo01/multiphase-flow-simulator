import numpy as np
import matplotlib.pyplot as plt

from oil import Oil
from water import Water
from gas import Gas
from flow import Flow
from tubing import Tubing
from helpers import bubble_point

tubing = Tubing(8000, 1.995, 90.0, 0.000902255639097744)
gas = Gas(0.65)
oil = Oil(25.0, gas)
water = Water(1.07, gas)
flow = Flow(600.0, 10.0, 0.0, tubing)

well_head_pressure = 5.0
temperature = 175.0
delta_l = 1

bubble_point = bubble_point(temperature, oil, water, flow)

x = []
y = []

pressure = well_head_pressure
for l in np.arange(0, tubing.length, delta_l): 
    x.append(l)
    y.append(pressure)

    gas.update_conditions(pressure, temperature)
    oil.update_conditions(pressure, bubble_point, temperature, flow.water_cut)
    water.update_conditions(pressure, bubble_point, temperature, flow.water_cut)
    flow.update_conditions(pressure, bubble_point, gas, oil, water)

    delta_p = flow.grav_pressure_gradient + flow.fric_pressure_gradient
    pressure = pressure - delta_p

plt.plot(x, y)
plt.grid(True)
plt.xlabel('Tubing length')
plt.ylabel('Pressure')
plt.show()
