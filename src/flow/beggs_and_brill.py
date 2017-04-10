import math
from .helpers import density_to_specific_gravity
from .helpers import estimate_fluid_property
from .flowpattern import FlowPattern

HORZ_LIQUID_HOLDUP_CONSTANTS = {
    FlowPattern.segregated:   (0.980, 0.4846, 0.0868),
    FlowPattern.intermittent: (0.845, 0.5351, 0.0173),
    FlowPattern.distributed:  (1.065, 0.5824, 0.0609)
}

LIQUID_HOLDUP_CONSTANTS = {
    FlowPattern.segregated:   (0.011, -3.7680, 3.5390, -1.6140),
    FlowPattern.intermittent: (2.960, 0.3050, -0.4473, 0.0978),
    FlowPattern.distributed:  (1.0, 1.0, 1.0, 1.0),
    FlowPattern.downward:     (4.700, -0.3692, 0.1244, -0.5056)
}


class BeggsAndBrill(object):
    """docstring for Flow"""
    def __init__(self):
        super(BeggsAndBrill, self).__init__()
        self.grav_pressure_gradient = 0
        self.fric_pressure_gradient = 0

    def update_conditions(self, mixture, tubing):
        froude_number = self._calc_froude_number(mixture.mixture_velocity, tubing.diameter)
        flow_pattern = self._calc_flow_pattern(froude_number, mixture.no_slip_liquid_fraction)
        liquid_velocity_number = self._calc_liquid_velocity_number(
            mixture.liquid_velocity,
            mixture.liquid_density,
            mixture.liquid_surface_tension
        )
        liquid_holdup = self._calc_liquid_holdup_with_incl(
            flow_pattern,
            froude_number,
            mixture.no_slip_liquid_fraction,
            liquid_velocity_number,
            tubing.inclination
        )
        mixture_density_holdup = estimate_fluid_property(
            mixture.liquid_density,
            mixture.gas.density,
            (1 - liquid_holdup)
        )
        mixture_specific_grav_holdup = density_to_specific_gravity(
            mixture_density_holdup
        )

        reynolds = self._calc_reynolds(
            mixture.mixture_density_no_slip,
            mixture.mixture_velocity,
            tubing.diameter,
            mixture.mixture_viscosity_no_slip
        )
        moody_friction = self._calc_moody_friction_factor(reynolds, tubing.rugosity)
        friction_factor = self._calc_friction_factor(
            mixture.no_slip_liquid_fraction,
            liquid_holdup,
            moody_friction
        )

        self.grav_pressure_gradient = self._calc_grav_pressure_gradient(
            mixture_specific_grav_holdup,
            tubing.inclination
        )
        self.fric_pressure_gradient = self._calc_fric_pressure_gradient(
            friction_factor,
            mixture.mixture_specific_grav_no_slip,
            mixture.mixture_velocity,
            tubing.diameter
        )

    def _calc_froude_number(self, mixture_velocity, tubing_diameter):
        """
        Calculates the mixture Froude number based on the mixture velocity (sum
        of all phases respective superficial velocities) and the tubing
        diameter.

        Args:
            mixture_velocity (double): Superficial mixture velocity
                (:math:`ft/s`)
            diameter (double): Tubing diameter (:math:`in`).

        Returns:
            THe Froude number.
        """
        return (mixture_velocity ** 2) / (32.174 * tubing_diameter / 12)

    def _calc_transition_froude_numbers(self, no_slip_liquid_fraction):
        """
        Calculates the Froude number limits used to determine flow pattern.

        Args:
            no_slip_liquid_fraction (double): The no slip liquid fraction

        Returns:
            A tuple with the four froude number limits in the format
            :math:`(Fr_1, Fr_2, Fr_3 and Fr_4)`.
        """
        return (
            316.0 * no_slip_liquid_fraction ** 0.302,  # Fr_1
            0.0009252 * no_slip_liquid_fraction ** -2.4684,  # Fr_2
            0.1 * no_slip_liquid_fraction ** -1.4516,  # Fr_3
            0.5 * no_slip_liquid_fraction ** -6.738  # Fr_4
        )

    def _calc_flow_pattern(self, froude_number, no_slip_liquid_fraction):
        """
        Returns the flow pattern based no the Froude number and the no slip
        liquid fraction.

        Args:
            froude_number (double): The mixture's froude number
            no_slip_liquid_fraction (double): The no slip liquid fraction

        Returns:
            The flow pattern using the `FlowPattern` enum.
        """
        fr1, fr2, fr3, fr4 = self._calc_transition_froude_numbers(
            no_slip_liquid_fraction
        )
        if froude_number >= fr1 or froude_number >= fr4:
            return FlowPattern.distributed
        elif froude_number > fr3:
            return FlowPattern.intermittent
        elif froude_number > fr2:
            return FlowPattern.transition
        else:
            return FlowPattern.segregated

    def _calc_horz_liquid_holdup(self,
                                 flow_pattern,
                                 froude_number,
                                 no_slip_liquid_fraction):
        """
        Returns the liquid fraction considering slippage for horizontal flow.

        Args:
            flow_pattern (FlowPattern): The flow pattern as determined by
                `flow_pattern`.
            froude_number (double): The mixture's froude number
            no_slip_liquid_fraction (double): The no slip liquid fraction

        Returns:
            The liquid fraction considering slippage for horizontal flow.
        """
        term_a, term_b, term_c = HORZ_LIQUID_HOLDUP_CONSTANTS[flow_pattern]
        answer = (
            term_a *
            no_slip_liquid_fraction ** term_b / (froude_number ** term_c)
        )
        return max(answer, no_slip_liquid_fraction)

    def _calc_liquid_velocity_number(self,
                                     liquid_superficial_velocity,
                                     liquid_density,
                                     superficial_tension):
        """
        Returns the liquid velocity number used in the correction of the liquid
        fraction for different inclinations.

        Args:
            _liquid_superficial_velocity (double): The superficial liquid velocity
                :math:`V_{so} + V_{sw}` in ft/s.
            _liquid_density (double): The liquid density in :math:`lbm/ft^3`.
            _superficial_tension (double): The superficial tension in
                :math:`dynes/cm`

        Returns:
            The liquid density in the same unit as the given oil and water density.
        """
        return (
            1.938 * liquid_superficial_velocity *
            (liquid_density / superficial_tension) ** (1 / 4)
        )

    def _calc_liquid_holdup_with_incl(self,
                                      flow_pattern,
                                      froude_number,
                                      no_slip_liquid_fraction,
                                      liquid_velocity_number,
                                      inclination):
        """
        Returns the liquid fraction considering slippage for any inclination.

        Args:
            flow_pattern (FlowPattern): The flow pattern as determined by
                `flow_pattern`.
            froude_number (double): The mixture's froude number.
            no_slip_liquid_fraction (double): The no slip liquid fraction.
            liquid_velocity_number (double): The liquid velocity number obtained
                using `liquid_velocity_number()`.
            inclination (double): The inclination angle with the horizontal in
                degrees.

        Returns:
            The liquid fraction considering slippage for any inclination.
        """
        if flow_pattern != FlowPattern.transition:
            horz_liquid_holdup = self._calc_horz_liquid_holdup(
                flow_pattern,
                froude_number,
                no_slip_liquid_fraction
            )

            term_d, term_e, term_f, term_g = LIQUID_HOLDUP_CONSTANTS[flow_pattern]

            inclination_rad = inclination * math.pi / 180

            c_parameter = max(0, (
                (1 - no_slip_liquid_fraction) *
                math.log(
                    term_d *
                    no_slip_liquid_fraction ** term_e *
                    liquid_velocity_number ** term_f *
                    froude_number ** term_g
                )
            ))
            if flow_pattern == FlowPattern.distributed:
                c_parameter = 0

            phi_parameter = (
                1 + c_parameter * (
                    math.sin(1.8 * inclination_rad) -
                    0.333 * (math.sin(1.8 * inclination_rad) ** 3)
                )
            )
            payne_factor = 0.924 if inclination > 0 else 0.685
            max_value = no_slip_liquid_fraction if inclination < 0 else 1
            min_value = no_slip_liquid_fraction if inclination >= 0 else 0
            return max(
                min_value,
                min(max_value, payne_factor * horz_liquid_holdup * phi_parameter)
            )
        else:
            _, fr2, fr3, _ = self._calc_transition_froude_numbers(
                no_slip_liquid_fraction
            )
            term_a = (fr3 - froude_number) / (fr3 - fr2)
            term_b = 1 - term_a
            answer = (
                (
                    term_a * self._calc_liquid_holdup_with_incl(
                        FlowPattern.segregated,
                        froude_number,
                        no_slip_liquid_fraction,
                        liquid_velocity_number,
                        inclination
                    )
                ) + (
                    term_b * self._calc_liquid_holdup_with_incl(
                        FlowPattern.intermittent,
                        froude_number,
                        no_slip_liquid_fraction,
                        liquid_velocity_number,
                        inclination
                    )
                )
            )
            max_value = no_slip_liquid_fraction if inclination < 0 else 1
            min_value = no_slip_liquid_fraction if inclination >= 0 else 0
            return max(
                min_value,
                min(max_value, answer)
            )

    def _calc_grav_pressure_gradient(self,
                                     mixture_specific_gravity,
                                     inclination):
        """
        Returns the pressure gradient due to gravity.

        Args:
            _liquid_density (double): The liquid density in :math:`lbm/ft^3`.
            _gas_density (double): The gas density in :math:`lbm/ft^3`.
            _liquid_fraction (double): The liquid fraction considering slippage
                and the given inclination.
            inclination (double): The inclination angle with the horizontal in
                degrees.

        Returns:
            The pressure gradient due to gravity :math:`\frac{dP}{dz}` in
            :math:`psi/ft`.
        """
        _inclination_rad = inclination * math.pi / 180
        return -0.433 * mixture_specific_gravity * math.sin(_inclination_rad)

    def _calc_reynolds(self, density, velocity, diameter, viscosity):
        return 124 * density * abs(velocity) * diameter / viscosity

    def _calc_moody_friction_factor(self, reynolds, rugosity):
        term_a = (
            2.457 * math.log(
                1 / (
                    (7 / reynolds) ** 0.9 + 0.27 * rugosity
                )
            )
        ) ** 16
        term_b = (37530 / reynolds) ** 16
        answer = (
            8 * (
                (8 / reynolds) ** 12 +
                1 / ((term_a + term_b) ** (3 / 2))
            ) ** (1 / 12)
        )
        return answer

    def _calc_friction_factor(self,
                              no_slip_liquid_fraction,
                              liquid_holdup,
                              moody_friction_factor):
        term_y = no_slip_liquid_fraction / (liquid_holdup ** 2)
        term_s = 0.0
        if term_y < 1.0 or term_y > 1.2:
            term_s = (
                math.log(term_y) /
                (
                    -0.0523 +
                    3.182 * math.log(term_y) -
                    0.8725 * (math.log(term_y) ** 2) +
                    0.01853 * (math.log(term_y) ** 4)
                )
            )
        else:
            math.log(2.2 * term_y - 1.2)
        return moody_friction_factor * math.exp(term_s)

    def _calc_fric_pressure_gradient(self,
                                     friction_factor,
                                     mixture_specific_gravity,
                                     mixture_velocity,
                                     diameter):
        _mixture_flow_rate = (
            math.pi * ((diameter / 12) ** 2) * mixture_velocity / 4
        ) * 86400 / 5.614583
        answer = (
            -1.1471e-5 * friction_factor * mixture_specific_gravity *
            _mixture_flow_rate ** 2 /
            (diameter ** 5)
        )
        return answer
