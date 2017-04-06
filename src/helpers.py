def density_to_specific_gravity(_density):
    return _density / 62.4


def specific_gravity_from_api(api_gravity):
    """
    Converts the oil's API gravity to specific gravity.
    """
    return 141.5 / (api_gravity + 131.5)


def free_gas_liquid_ratio(pressure,
                          bubble_pnt,
                          gas_solubility_in_oil,
                          gas_solubility_in_water,
                          _water_cut,
                          prod_gas_liquid_ratio):
    """
    Calculates the free gas liquid ratio.

    Args:
        pressure: Pressure at which the gas is (:math:`psig`). Note that this
            value is in :math:`psig`, so it is relative to the atmospheric
            pressure.
        bubble_pnt: Mixture's bubble point (:math:`psig`). Note that this
            value is in :math:`psig`, so it is relative to the atmospheric
            pressure.
        gas_solubility_in_oil: Gas solubility in oil, :math:`R_{so}` (in the
            same unit as :math:`GLR_p` and :math:`R_{sw}`, suggestion:
            :math:`scf/stb`).
        gas_solubility_in_water: Gas solubility in water, :math:`R_{sw}` (in
            the same unit as :math:`GLR_p` and :math:`R_{so}`, suggestion:
            :math:`scf/stb`).
        _water_cut: Water cut, WC.
        prod_gas_liquid_ratio: Production gas liquid ratio,
            :math:`GLR_p` (in the same unit as :math:`R_{so}` and
            :math:`R_{sw}`, suggestion: :math:`scf/stb`).

    Returns:
        The free gas liquid ratio if pressure is below bubble point. Otherwise,
        returns zero.
    """
    if pressure >= bubble_pnt:
        return 0
    else:
        return (prod_gas_liquid_ratio -
                gas_solubility_in_oil * (1 - _water_cut) -
                gas_solubility_in_water * _water_cut)


def bubble_point(temperature, mixture):
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
        rso = mixture.oil.calc_gas_solubility(bubble_point_,
                                              pressure_high,
                                              temperature)
        rsw = mixture.water.calc_gas_solubility(bubble_point_,
                                                pressure_high,
                                                temperature)

        error = (mixture.prod_gas_liquid_ratio -
                 (1 - mixture.water_cut) * rso - mixture.water_cut * rsw)

        if error > 0.0:
            pressure_low = bubble_point_
        else:
            pressure_high = bubble_point_

    return bubble_point_


def no_slip_fraction(fluid_velocity, total_velocity):
    """
    Calculates the no slip gas fraction of the produced fluid based on the
    superficial velocity of each phase. Note that the suggested unit is
    :math:`ft/s`, but as long as all velocities are in the same unit, any unit
    can be used.

    Args:
        _oil_velocity (double): Superficial oil velocity (:math:`ft/s`).
        _gas_velocity (double): Superficial gas velocity (:math:`ft/s`).
        _water_velocity (double): Superficial water velocity (:math:`ft/s`).

    Returns:
        THe no slip gas fraction.
    """
    return fluid_velocity / total_velocity


def estimate_fluid_property(_first_fluid_property,
                            _second_fluid_property,
                            _fraction_of_second_fluid):
    """
    Returns the property of the liquid mixture by weighting the water and oil
    properties by their respective no slip fractions.

    Args:
        _oil_property (double): The property value for the oil fraction.
        _water_property (double): The property value for the water fraction.
        _water_fraction (double): The fraction of water in the liquid phase of
            the mixture.

    Returns:
        The liquid mixture property in the same unit as the given oil and water
        properties
    """
    return (
        _first_fluid_property * (1 - _fraction_of_second_fluid) +
        _second_fluid_property * _fraction_of_second_fluid
    )
