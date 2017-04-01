import math


class Water(object):
    """docstring for Water"""
    def __init__(self, specific_gravity, bubble_point, dissolved_gas):
        super(Water, self).__init__()
        self.specific_gravity = specific_gravity
        self.bubble_point = bubble_point
        self.gas = dissolved_gas
        self.gas_solubility = None
        self.compressibility = None
        self.formation_volume_factor = None
        self.viscosity = None
        self.density = None
        self.water_gas_surface_tension = None

    def update_conditions(self, pressure, temperature, water_cut):
        self.gas_solubility = self._calc_gas_solubility(pressure, temperature)
        if pressure >= self.bubble_point:
            self.compressibility = self._calc_compressibility(
                pressure, temperature
            )
        self.formation_volume_factor = self._calc_formation_volume_factor(
            pressure, temperature, self.compressibility
        )
        self.viscosity = self._calc_viscosity(
            pressure, temperature
        )
        self.density = self._calc_density(
            self.gas_solubility, self.formation_volume_factor, water_cut
        )
        self.water_gas_surface_tension = self._calc_water_gas_surf_tension()

    def _calc_gas_solubility(self, pressure, temperature):
        """
        Calculates gas solubility in water (:math:`R_{sw}`) using Culberson and
        Maketta correlation. If pressure is higher than the mixture's bubble
        point, returns the :math:`R_{so}` at bubble point.

        Args:
            pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure
                and must be above bubble point.
            temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The gas solubility in water, :math:`R_{sw}` (:math:`scf/stb`).
        """
        term_a = (8.15839 -
                  6.12265e-2 * temperature +
                  1.91663e-4 * (temperature ** 2) -
                  2.1654e-7 * (temperature ** 3))

        term_b = (1.01021e-2 -
                  7.44241e-5 * temperature +
                  3.05553e-7 * (temperature ** 2) -
                  2.94883e-10 * (temperature ** 3))

        term_c = (-9.02505 +
                  0.130237 * temperature -
                  8.53425e-4 * (temperature ** 2) +
                  2.34122e-6 * (temperature ** 3) -
                  2.37049e-9 * (temperature ** 4)) * (10 ** -7)

        if pressure > self.bubble_point:
            pressure = self.bubble_point

        abs_pressure = pressure + 14.7
        return term_a + term_b * (abs_pressure) + term_c * (abs_pressure) ** 2

    def _calc_compressibility(self, pressure, temperature):
        """
        Calculates the isothermal water compressibility using Dodson and
        Standing correlation. This is the compressibility of the water as a
        liquid with a certain amount of gas in solution. It is valid above the
        bubble point only.

        Args:
            pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure
                and must be above bubble point.
            temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The compressibility of the water as a single-phase liquid with a
            certain amount of gas in solution in :math:`psi^{-1}`.
        """
        if pressure < self.bubble_point:
            raise ValueError('Pressure must be below bubble point.')

        gas_solubility_in_water_at_bp = self._calc_gas_solubility(
            self.bubble_point,
            temperature
        )

        term_a = 3.8546 - 1.34e-4 * (pressure + 14.7)
        term_b = -0.01052 + 4.77e-7 * (pressure + 14.7)
        term_c = 3.9267e-5 - 8.8e-10 * (pressure + 14.7)

        result = (
            (term_a + term_b * temperature + term_c * (temperature ** 2)) *
            (1 + 8.9e-3 * gas_solubility_in_water_at_bp) / 1e6
        )
        return result

    def _calc_formation_volume_factor(self,
                                      pressure,
                                      temperature,
                                      water_compressibility=0.0):
        """
        Calculates the water formation volume factor (:math:`B_w`) using
        Gould's correlation. The bubble point is necessary because a different
        correlation must be used below and above it and the water
        compressibility (:math:`c_w`) is only needed if pressure is above
        bubble point.

        Args:
            pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure
                and must be above bubble point.
            temperature (double): Temperature (fahrenheit degrees).
            water_compressibility: Waters compressibility (:math:`psi^{-1}`).
                Value can be omitted if pressure is below bubble point.

        Returns:
            The water formation volume factor, in :math:`bbl/stb`.
        """
        result = (1.0 +
                  1.2e-4 * (temperature - 60) +
                  1.0e-6 * (temperature - 60) ** 2)

        if pressure >= self.bubble_point:
            result = result - 3.33e-6 * (self.bubble_point + 14.7)
            result = (result * math.exp(water_compressibility *
                                        (self.bubble_point - pressure)))
        else:
            result = result - 3.33e-6 * (pressure + 14.7)

        return result

    def _calc_viscosity(self, pressure, temperature):
        """
        Calculates the water viscosity using the Kestin, Khalifa and Correa
        correlation.

        Args:
            _pressure (double): Pressure at which the oil is (:math:`psig`).
                Note that this value is relative to the atmospheric pressure.
            _temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The water viscosity in :math:`cp`.
        """

        viscosity = (
            109.574 * temperature ** (-1.12166) * (
                0.9994 + 4.0295e-5 * (pressure + 14.7) +
                3.1062e-9 * (pressure + 14.7) ** 2
            )
        )
        return viscosity

    def _calc_density(self,
                      gas_solubility_in_water,
                      water_formation_volume_factor,
                      water_cut):
        """
        Calculates the live water density at the given conditions

        Args:
            gas_solubility_in_water (double): Gas solubility in water,
                :math:`R_{sw}` (:math:`scf/stb`).
            water_formation_volume_factor (double): Water formation volume
                factor, :math:`B_w` (:math:`bbl/stb`)
            water_cut: Water cut, WC.

        Returns:
            The live water density in :math:`lbm/ft^3`.
        """
        density = (
            (
                350 * self.specific_gravity +
                0.0764 * self.gas.specific_gravity * gas_solubility_in_water *
                water_cut
            ) / (
                5.615 * water_formation_volume_factor
            )
        )
        return density

    def _calc_water_gas_surf_tension(self):
        """
        Returns the water - gas surface tension. Currently there is no
        correlation implemented for this property.

        Returns:
            The water - gas surface tension in :math:`dina/cm`
        """
        return 90.359630470339
