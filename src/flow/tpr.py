import numpy as np

class TPR(object):
    def __init__(self, correlation, tubing, mixture, temperature, well_head_pressure):
        self.correlation = correlation
        self.tubing = tubing
        self.mixture = mixture
        self.well_head_pressure = well_head_pressure
        self.temperature = temperature

        self.step = 10.
        self.max_flow_rate = 6500.

    def simulate(self):
        flow_rates = np.arange(1., self.max_flow_rate, self.step)
        pressures = []

        bubble_point = self.mixture.bubble_point(self.temperature)

        for flow_rate in flow_rates:
            pressure = self.well_head_pressure
            depth = 0.0

            last_pressure = 0.0
            last_depth = 0.0

            delta_p = 0.5
            last_limit = 5.0
            limit = 10.0
            while depth < self.tubing.length:
                self.mixture.update_conditions(
                    pressure, bubble_point, self.temperature, self.tubing.diameter, flow_rate
                )
                self.correlation.update_conditions(self.mixture, self.tubing)

                dp_dz = (
                    self.correlation.grav_pressure_gradient +
                    self.correlation.fric_pressure_gradient
                )
                last_depth = depth
                depth = depth + delta_p / abs(dp_dz)

                last_pressure = pressure
                pressure = pressure + delta_p

                if pressure >= limit:
                    delta_p = limit / 10.0
                    tmp_limit = limit
                    limit = last_limit * 10.0
                    last_limit = tmp_limit

            factor = (self.tubing.length - last_depth) / (depth - last_depth)
            pressure = last_pressure + factor * (pressure - last_pressure)
            pressures.append(pressure)

        return (flow_rates, pressures)
