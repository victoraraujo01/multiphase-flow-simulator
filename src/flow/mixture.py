import math
from .helpers import free_gas_liquid_ratio
from .helpers import no_slip_fraction
from .helpers import density_to_specific_gravity
from .helpers import estimate_fluid_property


class Mixture(object):
    """docstring for Flow"""
    def __init__(self, gas, oil, water, water_cut, prod_gas_liquid_ratio):
        super(Mixture, self).__init__()
        self.gas = gas
        self.oil = oil
        self.water = water
        self.water_cut = water_cut
        self.prod_gas_liquid_ratio = prod_gas_liquid_ratio

        self.prod_flow_rate = None
        self.liquid_velocity = None
        self.mixture_velocity = None
        self.no_slip_liquid_fraction = None
        self.liquid_density = None
        self.liquid_surface_tension = None
        self.liquid_viscosity = None
        self.mixture_viscosity_no_slip = None
        self.mixture_density_no_slip = None
        self.mixture_specific_grav_no_slip = None

    def update_conditions(self, pressure, bubble_point, temperature, tubing_diameter, flow_rate):
        self.gas.update_conditions(pressure, temperature)
        self.oil.update_conditions(pressure, bubble_point, temperature, self.water_cut)
        self.water.update_conditions(pressure, bubble_point, temperature, self.water_cut)

        free_glr = free_gas_liquid_ratio(
            pressure,
            bubble_point,
            self.oil.gas_solubility,
            self.water.gas_solubility,
            self.water_cut,
            self.prod_gas_liquid_ratio
        )

        self.prod_flow_rate = flow_rate
        gas_flow_rate = self.in_situ_flow_rate(
            self.gas.formation_volume_factor, free_glr
        )
        oil_flow_rate = self.in_situ_flow_rate(
            self.oil.formation_volume_factor, (1 - self.water_cut)
        )
        water_flow_rate = self.in_situ_flow_rate(
            self.water.formation_volume_factor, self.water_cut
        )
        gas_velocity = self.superficial_velocity(gas_flow_rate, tubing_diameter)
        oil_velocity = self.superficial_velocity(oil_flow_rate, tubing_diameter)
        water_velocity = self.superficial_velocity(water_flow_rate, tubing_diameter)

        self.liquid_velocity = oil_velocity + water_velocity
        self.mixture_velocity = gas_velocity + self.liquid_velocity

        self.no_slip_liquid_fraction = no_slip_fraction(
            self.liquid_velocity,
            self.mixture_velocity
        )
        water_fraction = no_slip_fraction(
            water_velocity,
            self.liquid_velocity
        )

        self.liquid_density = estimate_fluid_property(
            self.oil.density,
            self.water.density,
            water_fraction
        )
        self.liquid_surface_tension = estimate_fluid_property(
            self.oil.oil_gas_surface_tension,
            self.water.water_gas_surface_tension,
            water_fraction
        )
        self.liquid_viscosity = estimate_fluid_property(
            self.oil.viscosity,
            self.water.viscosity,
            water_fraction
        )

        self.mixture_viscosity_no_slip = estimate_fluid_property(
            self.liquid_viscosity,
            self.gas.viscosity,
            (1 - self.no_slip_liquid_fraction)
        )
        self.mixture_density_no_slip = estimate_fluid_property(
            self.liquid_density,
            self.gas.density,
            (1 - self.no_slip_liquid_fraction)
        )
        self.mixture_specific_grav_no_slip = density_to_specific_gravity(
            self.mixture_density_no_slip
        )

    def bubble_point(self, temperature):
        """
        Calculates the mixture's bubble point based on the water cut and the
        prouction gas liquid ratio.

        Args:
            _water_cut: Water cut, WC.
            prod_gas_liquid_ratio: Production gas liquid ratio,
                :math:`GLR_p` (in the same unit as :math:`R_{so}` and
                :math:`R_{sw}`, suggestion: :math:`scf/stb`).

        Returns:
            The mixture's bubble point Pb (psi).
        """
        pressure_low = 0.0
        pressure_high = 100000.0
        bubble_point_ = 0.0
        error = 1.0
        while abs(error) > 1e-10:
            bubble_point_ = (pressure_low + pressure_high) / 2
            rso = self.oil.calc_gas_solubility(bubble_point_,
                                               pressure_high,
                                               temperature)
            rsw = self.water.calc_gas_solubility(bubble_point_,
                                                 pressure_high,
                                                 temperature)

            error = (self.prod_gas_liquid_ratio -
                     (1 - self.water_cut) * rso - self.water_cut * rsw)

            if error > 0.0:
                pressure_low = bubble_point_
            else:
                pressure_high = bubble_point_

        return bubble_point_

    def in_situ_flow_rate(self, formation_volume_factor, fraction):
        """
        Calculates the in-situ gas flow rate at the given conditions

        Args:
            _liquid_flow_rate (double): Total liquid flow rate (:math:`bpd`).
            _gas_formation_volume_factor (double): Gas formation volume factor,
                :math:`B_g` (:math:`bbl/scf`). Pay attention to the unit used.
            _free_gas_liquid_ratio (double): The free gas liquid ratio
                (:math:`scf/stb`).

        Returns:
            The in-situ gas flow rate in :math:`ft^3/s`.
        """
        answer = (
            fraction * self.prod_flow_rate * formation_volume_factor *
            0.000064984  # reduced 5.614583 / 86400 to improve speed
        )
        return answer

    def superficial_velocity(self, in_situ_flow_rate, tubing_diameter):
        """
        Transforms the passed in-situ flow rate into superficial velocity at
        the given diameter.

        Args:
            in_situ_flow_rate (double): In situ flow rate (:math:`ft^3/s`).
            diameter (double): The tubing diameter (:math:`in`).

        Returns:
            The superficial velocity in :math:`ft/s`.
        """
        return (
            4 * in_situ_flow_rate / (math.pi * (tubing_diameter / 12) ** 2)
        )
