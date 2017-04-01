import math


def density_to_specific_gravity(_density):
    return _density / 62.4


def specific_gravity_from_api(api_gravity):
    """
    Converts the oil's API gravity to specific gravity.
    """
    return 141.5/(api_gravity + 131.5)


def water_cut(_oil_flow_rate, _water_flow_rate):
    """
    Calculates the water cut based on the passed oil and water flow rates.
    Supplied oil and water rates must be in the same unit.

    Args:
        oil_flow_rate: Current oil flow rate.
        water_flow_rate: Current water flow rate.

    Returns:
        THe water cut ratio, known as WC or BSW (between 0 and 1).
    """
    return _water_flow_rate / (_oil_flow_rate + _water_flow_rate)


def production_gas_liquid_ratio(_production_gas_oil_ratio, _water_cut):
    """
    Calculates the production gas liquid ratio based on the production gas oil
    ratio and water cut. The calculated production gas liquid ratio will be on
    the same unit as the gas oil ratio.

    Args:
        _production_gas_oil_ratio: Production gas oil ratio, :math:`RGO` (
            suggested unit: :math:`scf/stb`).
        _water_cut: Water cut.

    Returns:
        The production gas liquid ratio, :math:`RGL` (in the same unit as the
        supplied :math:`RGO`).
    """
    return _production_gas_oil_ratio * (1 - _water_cut)


def free_gas_liquid_ratio(_pressure,
                          _bubble_point,
                          _gas_solubility_in_oil,
                          _gas_solubility_in_water,
                          _water_cut,
                          _production_gas_liquid_ratio):
    """
    Calculates the free gas liquid ratio.

    Args:
        _pressure: Pressure at which the gas is (:math:`psig`). Note that this
            value is in :math:`psig`, so it is relative to the atmospheric
            pressure.
        _bubble_point: Mixture's bubble point (:math:`psig`). Note that this
            value is in :math:`psig`, so it is relative to the atmospheric
            pressure.
        _gas_solubility_in_oil: Gas solubility in oil, :math:`R_{so}` (in the
            same unit as :math:`GLR_p` and :math:`R_{sw}`, suggestion:
            :math:`scf/stb`).
        _gas_solubility_in_water: Gas solubility in water, :math:`R_{sw}` (in
            the same unit as :math:`GLR_p` and :math:`R_{so}`, suggestion:
            :math:`scf/stb`).
        _water_cut: Water cut, WC.
        _production_gas_liquid_ratio: Production gas liquid ratio,
            :math:`GLR_p` (in the same unit as :math:`R_{so}` and
            :math:`R_{sw}`, suggestion: :math:`scf/stb`).

    Returns:
        The free gas liquid ratio if pressure is below bubble point. Otherwise,
        returns zero.
    """
    if _pressure >= _bubble_point:
        return 0
    else:
        return (_production_gas_liquid_ratio -
                _gas_solubility_in_oil * (1 - _water_cut) -
                _gas_solubility_in_water * _water_cut)

# TODO: FIX
def bubble_point(_temperature,
                 _gas_specific_gravity,
                 _oil_api_gravity,
                 _water_cut,
                 _production_gas_liquid_ratio):
    """
    Calculates the mixture's bubble point based on the water cut and the
    prouction gas liquid ratio.

    Args:
        _water_cut: Water cut, WC.
        _production_gas_liquid_ratio: Production gas liquid ratio,
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
        bubble_point_ = (pressure_low + pressure_high)/2
        rso = self.gas.solubility_in_oil(bubble_point,
                                         pressure_high,
                                         _temperature,
                                         _gas_specific_gravity,
                                         _oil_api_gravity)
        rsw = self.gas.solubility_in_water(bubble_point,
                                           pressure_high,
                                           _temperature)

        error = (_production_gas_liquid_ratio -
                 (1 - _water_cut) * rso - _water_cut * rsw)

        if error > 0.0:
            pressure_low = bubble_point_
        else:
            pressure_high = bubble_point_

    return bubble_point


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
