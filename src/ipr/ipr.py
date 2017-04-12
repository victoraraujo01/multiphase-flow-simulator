import math
from . import ipr_tests_analysis as analysis

class IPR(object):

    def __init__(self, b_param, avg_pressure, bubble_point, undersaturated_pi):
        self.avg_pressure = avg_pressure
        self.bubble_point = bubble_point
        self.undersaturated_pi = undersaturated_pi
        self.b_param = b_param

        self.flow_rate_at_bubble_point = self.calc_flow_rate_at_bubble_point()
        self.max_flow_rate = self.calc_max_flow_rate()

    @classmethod
    def from_tests(cls, b_param, bubble_point, first_test, secnd_test):
        high_pressure_test = first_test if first_test[0] > secnd_test[0] else secnd_test
        low_pressure_test = first_test if first_test[0] <= secnd_test[0] else secnd_test

        high_pressure, _ = high_pressure_test
        low_pressure, _ = low_pressure_test

        result = None
        if high_pressure > bubble_point and low_pressure > bubble_point:
            result = analysis.tests_above_bp(high_pressure_test, low_pressure_test)
        elif high_pressure > bubble_point and low_pressure < bubble_point:
            result = analysis.tests_on_opposite_sides_of_bp(
                b_param, bubble_point, high_pressure_test, low_pressure_test
            )
        else:
            result = analysis.tests_below_bp(
                b_param, bubble_point, high_pressure_test, low_pressure_test
            )

        avg_pressure, undersaturated_pi = result
        return cls(b_param, avg_pressure, bubble_point, undersaturated_pi)

    @classmethod
    def from_past_ipr(cls, past_ipr, new_avg_pressure):
        new_fetkovic_pi = past_ipr.undersaturated_pi * new_avg_pressure / past_ipr.avg_pressure
        return cls(past_ipr.b_param, new_avg_pressure, past_ipr.bubble_point, new_fetkovic_pi)

    @classmethod
    def saturated_ipr_from_test(cls, b_param, avg_pressure, first_test):
        well_pressure, flow_rate = first_test
        ipr = cls(b_param, avg_pressure, 1e10, 1)
        ipr.max_flow_rate = (
            flow_rate /
            (
                1.0 + b_param * (well_pressure / avg_pressure) -
                (1.0 + b_param) * (well_pressure / avg_pressure) ** 2
            )
        )
        return ipr

    def calc_flow_rate_at_bubble_point(self):
        return max(0.0, (self.avg_pressure - self.bubble_point) * self.undersaturated_pi)

    def calc_max_flow_rate(self):
        pressure = self.bubble_point
        if self.avg_pressure < self.bubble_point:
            pressure = self.avg_pressure

        qmax = (
            self.flow_rate_at_bubble_point +
            self.undersaturated_pi * pressure /
            (2.0 + self.b_param)
        )
        return qmax

    def flow_rate(self, well_pressure):
        if well_pressure >= self.bubble_point:
            return self.flow_rate_above_bubble_point(well_pressure)
        else:
            return self.flow_rate_below_bubble_point(well_pressure)

    def flow_rate_above_bubble_point(self, well_pressure):
        return max(0, self.undersaturated_pi * (self.avg_pressure - well_pressure))

    def flow_rate_below_bubble_point(self, well_pressure):
        pressure = self.bubble_point
        if self.avg_pressure < self.bubble_point:
            pressure = self.avg_pressure

        flow_rate = (
            (
                1.0 + self.b_param * (well_pressure / pressure) -
                (1.0 + self.b_param) * (well_pressure / pressure) ** 2
            ) * (self.max_flow_rate - self.flow_rate_at_bubble_point) +
            self.flow_rate_at_bubble_point
        )
        return flow_rate

    def pressure(self, flow_rate):
        if flow_rate <= self.flow_rate_at_bubble_point:
            return self.pressure_above_bubble_point(flow_rate)
        else:
            return self.pressure_below_bubble_point(flow_rate)

    def pressure_above_bubble_point(self, flow_rate):
        return self.avg_pressure - flow_rate / self.undersaturated_pi

    def pressure_below_bubble_point(self, flow_rate):
        pressure = self.bubble_point
        if self.avg_pressure < self.bubble_point:
            pressure = self.avg_pressure

        term_a = -(1 + self.b_param) / (pressure ** 2)
        term_b = self.b_param / pressure
        term_c = (
            1 - (flow_rate - self.flow_rate_at_bubble_point) /
            (self.max_flow_rate - self.flow_rate_at_bubble_point)
        )
        delta = term_b ** 2 - 4 * term_a * term_c

        pressure = 0
        if delta >= 0:
            pressure = max(0, (-term_b - math.sqrt(delta)) / (2 * term_a))
        return pressure
