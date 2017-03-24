"""
Formulas
"""


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
    Calculates the production gas oil ratio based on the production gas liquid
    ratio and water cut. The calculated production gas oil ratio will be on the
    same unit as the gas liquid ratio.

    Args:
        _production_gas_oil_ratio: Production gas oil ratio, RGO (suggested
                                  unit: scf/stb).
        _water_cut: Water cut.

    Returns:
        The production gas liquid ratio, RGL (in the same unit as the supplied
        RGO).
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
        _pressure: Pressure at which the gas is (psig). Note that this value is
                   in psig, so it is relative to the atmospheric pressure.
        _bubble_point: Mixture's bubble point (psig).  Note that this value is
                       in psig, so it is relative to the atmospheric pressure.
        _gas_solubility_in_oil: Gas solubility in oil, Rso (in the same unit
                                as GLR_p and Rsw, suggestion: scf/stb).
        _gas_solubility_in_water: Gas solubility in water, Rsw (in the same
                                  unit as GLR_p and Rso, suggestion: scf/stb).
        _water_cut: Water cut, WC.
        _production_gas_liquid_ratio: Production gas liquid ratio, GLR_p (in
                                      the same unit as Rso and Rsw, suggestion:
                                      scf/stb).

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


def specific_gravity_from_api(api_gravity):
    """
    Converts the oil's API gravity to specific gravity.
    """
    return 141.5/(api_gravity + 131.5)
