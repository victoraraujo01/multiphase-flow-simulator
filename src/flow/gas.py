import math


class Gas(object):
    """Class that represents the gas phase"""
    def __init__(self, specific_gravity):
        super(Gas, self).__init__()
        self.specific_gravity = specific_gravity
        self.deviation_factor = None
        self.formation_volume_factor = None
        self.density = None
        self.viscosity = None

    def update_conditions(self, pressure, temperature):
        self.deviation_factor = self.calc_deviation_factor(
            pressure, temperature
        )
        self.formation_volume_factor = self.calc_formation_volume_factor(
            pressure, temperature, self.deviation_factor, False
        )
        self.density = self.calc_density(self.formation_volume_factor, False)
        self.viscosity = self.calc_viscosity(temperature, self.density)

    def calc_deviation_factor(self, pressure, temperature):
        """
        Calculates the gas deviation factor or compressibility factor :math:`Z`
        using Papay's correlation.

        Args:
            pressure (double): Pressure at which the gas is (:math:`psig`).
                Note that this value is in psig, so it is relative to the
                atmospheric pressure.
            temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The gas deviation factor.
        """
        pseudo_critical_temperature = (168. +
                                       325. * self.specific_gravity -
                                       12.5 * self.specific_gravity ** 2)
        pseudo_critical_pressure = (677. +
                                    15.0 * self.specific_gravity -
                                    37.5 * self.specific_gravity ** 2)

        pseudo_reduced_temperature = ((temperature + 460) /
                                      pseudo_critical_temperature)
        pseudo_reduced_pressure = ((pressure + 14.7) /
                                   pseudo_critical_pressure)

        pseudo_reduced_ratio = (pseudo_reduced_pressure /
                                pseudo_reduced_temperature)

        deviation_factor = (1 - pseudo_reduced_ratio *
                            (0.3675 - 0.04188423 * pseudo_reduced_ratio))
        return deviation_factor

    def calc_formation_volume_factor(self,
                                     pressure,
                                     temperature,
                                     deviation_factor,
                                     in_cubic_feet=True):
        """
        Calculates the gas formation volume factor. This is a convenience
        method that uses the gas specific gravity to calculate the gas
        deviation factor under the supplied pressure and temperature and
        under standard conditions.

        Args:
            pressure (double): Pressure at which the gas is (:math:`psig`).
                Note that this value is in psig, so it is relative to the
                atmospheric pressure.
            temperature (double): Temperature (fahrenheit degrees).
            in_cubic_feet (boolean, optional): If ``true``, result will be in
                :math:`ft^3/scf`. If set to ``false``, result will be in
                :math:`bbl/scf`.

        Returns:
            The gas formation volume factor.
        """
        deviation_factor_std = self.calc_deviation_factor(0, 60.0)

        conversion_factor = 0.028269
        if not in_cubic_feet:
            conversion_factor = 0.00503475

        gas_formation_volume_factor = (
            conversion_factor *
            (temperature + 460) / (pressure + 14.7) *
            deviation_factor / deviation_factor_std
        )

        return gas_formation_volume_factor

    def calc_density(self,
                     gas_formation_volume_factor=1.0,
                     bg_in_cubic_feet=True):
        """
        Calculates the gas density at standard conditions (if no gas formation
        volume factor, :math:`B_g`, is given) or uses the given :math:`B_g` to
        calculate it at different conditions.

        Args:
            gas_formation_volume_factor (double, optional): Gas' formation
                volume factor :math:`B_g` (:math:`ft^3/scf` or :math:`bbl/scf`
                - use argument ``bg_in_cubic_feet`` accordingly).
            in_cubic_feet (boolean, optional): must be ``True`` if supplied gas
                formation volume factor is in :math:`ft^3/scf` or ``False`` if
                it's in :math:`bbl/scf`. Defaults to ``True``.

        Returns:
            The gas density in :math:`lbm/ft^3`.
        """
        # conversion_factor = 14.7 * 28.97 / (10.7316 * 520)
        conversion_factor = 0.0764106

        if not bg_in_cubic_feet:
            conversion_factor = conversion_factor / 5.614583333

        density = (conversion_factor * self.specific_gravity /
                   gas_formation_volume_factor)
        return density

    def calc_viscosity(self, temperature, gas_density):
        """
        Calculates the gas viscosity using the Lee et al. correlation.

        Args:
            _temperature (double): Temperature (fahrenheit degrees).

        Returns:
            The gas viscosity in :math:`cp`.
        """
        molecular_weight = 28.97 * self.specific_gravity
        x_exponent = 3.5 + 986 / (temperature + 460) + 0.01 * molecular_weight
        y_exponent = 2.4 - 0.2 * x_exponent
        _gas_viscosity = (
            (9.4 + 0.02 * molecular_weight) * ((temperature + 460.) ** 1.5) /
            (209. + 19. * molecular_weight + temperature + 460) *
            10 ** (-4) * math.exp(
                x_exponent * (gas_density / 62.4) ** y_exponent
            )
        )
        return _gas_viscosity
