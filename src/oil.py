import math
from helpers import specific_gravity_from_api


class Oil(object):
    """Class that represents the oil phase"""
    def __init__(self, api_gravity, dissolved_gas):
        super(Oil, self).__init__()
        self.api_gravity = api_gravity
        self.gas = dissolved_gas
        self.gas_solubility = None
        self.compressibility = None
        self.formation_volume_factor = None
        self.viscosity = None
        self.density = None
        self.oil_gas_surface_tension = None

    @property
    def specific_gravity(self):
        """
        Calculates the oil specific gravity based on it's API gravity.

        Returns:
            The oil specific gravity.
        """
        return specific_gravity_from_api(self.api_gravity)

    def update_conditions(self,
                          pressure,
                          bubble_point,
                          temperature,
                          water_cut):
        self.gas_solubility = self.calc_gas_solubility(
            pressure, bubble_point, temperature
        )
        if pressure >= bubble_point:
            self.compressibility = self.calc_compressibility(
                pressure, bubble_point, temperature
            )
        self.formation_volume_factor = self.calc_formation_volume_factor(
            pressure,
            bubble_point,
            temperature,
            self.gas_solubility,
            self.compressibility
        )
        self.viscosity = self.calc_viscosity(
            pressure, bubble_point, temperature, self.gas_solubility
        )
        self.density = self.calc_density(
            self.gas_solubility, self.formation_volume_factor, water_cut
        )
        self.oil_gas_surface_tension = self.calc_oil_gas_surf_tension(
            temperature, self.gas_solubility
        )

    def calc_gas_solubility(self, pressure, bubble_point, temperature):
        """
        Calculates gas solubility in oil (:math:`R_{so}`) using Standing's
        correlation. If pressure is higher than the mixture's bubble point,
        returns the Rso at bubble point.

        Args:
            pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure
                and must be above bubble point.
            temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The gas solubility in oil, :math:`R_{so}` (:math:`scf/stb`).
        """
        if pressure > bubble_point:
            pressure = bubble_point

        exponent = 0.0125 * self.api_gravity - 0.00091 * temperature
        first_term = (pressure + 14.7)/18.2 + 1.4

        answer = (
            self.gas.specific_gravity * (first_term * 10 ** exponent) ** 1.2048
        )
        return answer

    def calc_compressibility(self, pressure, bubble_point, temperature):
        """
        Calculates the isothermal oil compressibility using Vasquez
        correlation. This is the compressibility of the oil as a single-phase
        liquid with a certain amount of gas in solution. It is valid above the
        bubble point only (Vasquez-Beggs Correlation).

        Args:
            pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure
                and must be above bubble point.
            temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The compressibility of the oil as a single-phase liquid with a
            certain amount of gas in solution in :math:`psi^{-1}`.
        """
        if pressure < bubble_point:
            raise ValueError('Pressure must be below bubble point.')

        gas_solubility_in_oil_at_bp = self.calc_gas_solubility(
            bubble_point,
            bubble_point,
            temperature
        )

        numerator = (-1433 +
                     5 * gas_solubility_in_oil_at_bp +
                     17.2 * temperature -
                     1180 * self.gas.specific_gravity +
                     12.61 * self.api_gravity)
        denominator = (pressure + 14.7) * (10 ** 5)
        return numerator/denominator

    def calc_formation_volume_factor(self,
                                     pressure,
                                     bubble_point,
                                     temperature,
                                     gas_solubility_in_oil,
                                     oil_compressibility=0.0):
        """
        Calculates the oil formation volume factor (:math:`B_o`) using
        Standing's correlation. The bubble point is necessary because a
        different correlation must be used below and above it. Also, it is
        important to remember that if pressure is above bubble point, the gas
        solubility in oil supplied must be the one at the bubble point
        (:math:`R_{sob}`). Finally, the oil compressibility (:math:`c_o`) is
        only needed if pressure is above bubble point.

        Args:
            pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure.
            temperature (double): Temperature (fahrenheit degrees).
            gas_solubility_in_oil_ (double): Gas solubility in oil,
                :math:`R_{so}` (in :math:`scf/stb`). **If pressure is above
                bubble point, the gas solubility in oil supplied must be the
                one at the bubble point (:math:`R_{sob}`).**
            oil_compressibility (double, optional): Oil's compressibility
                (:math:`psi^{-1}`). Value can be omitted if pressure is below
                bubble point.

        Returns:
            The oil formation volume factor, in :math:`bbl/stb`.
        """
        result = (0.9759 + 12e-5 * (
            gas_solubility_in_oil *
            math.sqrt(self.gas.specific_gravity/self.specific_gravity) +
            1.25 * temperature
        ) ** 1.2)

        if pressure > bubble_point:
            result = (
                result *
                math.exp(oil_compressibility * (bubble_point - pressure))
            )

        return result

    def calc_dead_oil_viscosity(self, temperature):
        """
        Calculates the dead oil viscosity using the Beggs and Robinson
        correlation.

        Args:
            temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The dead oil viscosity in :math:`cp`.
        """
        term_x = (10 ** (3.0324 - 0.02023 * self.api_gravity) /
                  (temperature ** 1.163))
        dead_oil_viscosity = 10 ** term_x - 1
        return dead_oil_viscosity

    def calc_viscosity(self,
                       pressure,
                       bubble_point,
                       temperature,
                       gas_solubility_in_oil):
        """
        Calculates the live oil viscosity. If pressure is below bubble point,
        the Beggs and Robinson correlation will be used. Instead, if it is
        above bubble point, the Beggs and Velasquez correlation will be used.

        Args:
            pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure.
            temperature (double): Temperature (fahrenheit degrees).
            gas_solubility_in_oil_ (double): Gas solubility in oil,
                :math:`R_{so}` (in :math:`scf/stb`). **If pressure is above
                bubble point, the gas solubility in oil supplied must be the
                one at the bubble point (:math:`R_{sob}`)**.

        Returns:
            The live oil viscosity in :math:`cp`.
        """
        dead_oil_viscosity = self.calc_dead_oil_viscosity(temperature)

        live_oil_viscosity = (
            10.715 *
            (gas_solubility_in_oil + 100) ** (-0.515) *
            dead_oil_viscosity **
            (5.44 * (gas_solubility_in_oil + 150) ** (-0.338))
        )

        if pressure > bubble_point:
            live_oil_viscosity = (
                live_oil_viscosity *
                ((pressure + 14.7) / (bubble_point + 14.7)) ** (
                    2.6 * (pressure + 14.7) ** 1.187 * math.exp(
                        -11.513 - 8.98e-5 * (pressure + 14.7)
                    )
                )
            )

        return live_oil_viscosity

    def calc_dead_density(self, in_cubic_feet=False):
        """
        Calculates the dead oil density at standard conditions

        Args:
            in_cubic_feet (boolean, optional): if ``True`` density will be in
                :math:`lbm/scf`. Otherwise it will be in :math:`lbm/stb`.
                Defaults to ``False``.

        Returns:
            The oil density in :math:`lbm/scf` or :math:`lbm/stb`.
        """
        water_density = 62.4
        if not in_cubic_feet:
            water_density = 350.0
        return water_density * self.specific_gravity

    def calc_density(self,
                     gas_solubility_in_oil,
                     oil_formation_volume_factor,
                     water_cut):
        """
        Calculates the live oil density at the given conditions

        Args:
            gas_solubility_in_oil (double): Gas solubility in oil,
                :math:`R_{so}` (:math:`scf/stb`).
            oil_formation_volume_factor (double): Oil formation volume factor,
                :math:`B_o` (:math:`bbl/stb`)
            water_cut: Water cut, WC.

        Returns:
            The live oil density in :math:`lbm/ft^3`.
        """
        density = ((
            350 * self.specific_gravity +
            0.0764 * self.gas.specific_gravity * gas_solubility_in_oil *
            (1 - water_cut)
        ) / (
            5.615 * oil_formation_volume_factor
        ))
        return density

    def calc_dead_oil_gas_surf_tension(self, temperature):
        """
        Calculates the dead oil - gas surface tension using Abdul-Majeed
        correlation (an update to Baker and Swerdloff's correlation). Source:
        http://petrowiki.org/Interfacial_tension#Water-hydrocarbon_surface_tension

        Args:
            temperature (double): Temperature (fahrenheit degrees).
            oil_api_gravity: Oil's API gravity (API degrees).

        Returns:
            The dead oil - gas surface tension in :math:`dina/cm`
        """
        return (
            (1.17013 - 1.694e-3 * temperature) *
            (38.085 - 0.259 * self.api_gravity)
        )

    def calc_oil_gas_surf_tension(self,
                                  temperature,
                                  gas_solubility_in_oil):
        """
        Corrects the dead oil - gas surface tension using Abdul-Majeed proposed
        method to obtain the surface tension between live oil and gas. Source:
        http://petrowiki.org/Interfacial_tension#Water-hydrocarbon_surface_tension

        Args:
            temperature (double): Temperature (fahrenheit degrees).
            gas_solubility_in_oil_ (double): Gas solubility in oil,
                :math:`R_{so}` (in :math:`scf/stb`).
        Returns:
            The oil - gas surface tension in :math:`dina/cm`
        """
        dead_og_surf_tension = self.calc_dead_oil_gas_surf_tension(
            temperature
        )
        return dead_og_surf_tension * (
            0.056379 +
            0.94362 * math.exp(
                -3.8491e-3 * gas_solubility_in_oil
            )
        )
