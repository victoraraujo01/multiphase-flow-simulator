from oil import Oil
from water import Water
from gas import Gas
from flow import Flow
from tubing import Tubing

bubble_point = 83.44862

tubing = Tubing(1.995, 90.0, 0.000902255639097744)
gas = Gas(0.65)
oil = Oil(25.0, bubble_point, gas)
water = Water(1.07, bubble_point, gas)
flow = Flow(600.0, 10.0, 0.0, tubing)

pressure = 5.0
temperature = 175.0

gas.update_conditions(pressure, temperature)
oil.update_conditions(pressure, temperature, flow.water_cut)
water.update_conditions(pressure, temperature, flow.water_cut)
flow.update_conditions(pressure, gas, oil, water)

print(flow.grav_pressure_gradient, flow. fric_pressure_gradient)
