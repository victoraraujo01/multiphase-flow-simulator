class IPR(object):

    def __init__(self, avg_pressure, bubble_point, undersaturated_ip):
        self.avg_pressure = avg_pressure
        self.bubble_point = bubble_point
        self.undersaturated_ip = undersaturated_ip
        self.b_param = -0.2

        self.flow_rate_at_bubble_point = self.calc_flow_rate_at_bubble_point()
        self.max_flow_rate = self.calc_max_flow_rate()

    def calc_flow_rate_at_bubble_point(self):
        return max(0.0, (self.avg_pressure - self.bubble_point) * self.undersaturated_ip)

    def calc_max_flow_rate(self):
        qmax = (
            self.flow_rate_at_bubble_point +
            self.undersaturated_ip * self.bubble_point /
            (2.0 + self.b_param)
        )
        return qmax

    def flow_rate(self, well_pressure):
        if well_pressure >= self.bubble_point:
            return self.flow_rate_above_bubble_point(well_pressure)
        else:
            return self.flow_rate_below_bubble_point(well_pressure)

    def flow_rate_above_bubble_point(self, well_pressure):
        return self.undersaturated_ip * (self.avg_pressure - well_pressure)

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
