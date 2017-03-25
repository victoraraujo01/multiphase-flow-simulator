"""
Formulas
"""

import math


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


def gas_density(_gas_specific_gravity,
                _gas_formation_volume_factor=1.0,
                _bg_in_cubic_feet=True):
    """
    Calculates the gas density at standard conditions (if no gas formation
    volume factor, Bg, is given) or uses the given Bg to calculate it at
    different conditions.

    Args:
        _gas_specific_gravity (double):  Gas' specific gravity (no unit)
        _gas_formation_volume_factor (double, optional): Gas' formation volume
            factor (:math:`ft^3/scf` or :math:`bbl/scf` - use argument
            bg_in_cubic_feet accordingly).
        _in_cubic_feet (boolean, optional): must be `True` if supplied gas
            formation volume factor is in :math:`ft^3/scf` or `False` if it's
            in :math:`bbl/scf`. Defaults to True.

    Returns:
        The gas density in :math:`lbm/ft^3`.
    """
    # conversion_factor = 14.7 * 28.97 / (10.7316 * 520)
    conversion_factor = 0.0764106
    if not _bg_in_cubic_feet:
        conversion_factor = conversion_factor / 5.614

    density = (conversion_factor * _gas_specific_gravity /
               _gas_formation_volume_factor)
    return density


def dead_oil_density(_oil_specific_gravity, _in_cubic_feet=False):
    """
    Calculates the dead oil density at standard conditions

    Args:
        _oil_specific_gravity (double): Oil's specific gravity (no unit).
        _in_cubic_feet (boolean, optional): if `True` density will be in
            :math:`lbm/scf`. Otherwise it will be in :math:`lbm/stb`. Defaults
            to `False`.
    Returns:
        The oil density in :math:`lbm/scf` or :math:`lbm/stb`.
    """
    _water_density = 62.4
    if not _in_cubic_feet:
        _water_density = 350.0
    return _water_density * _oil_specific_gravity


def live_oil_density(_oil_specific_gravity,
                     _gas_specific_gravity,
                     _gas_solubility_in_oil,
                     _oil_formation_volume_factor,
                     _water_cut):
    """
    Calculates the live oil density at the given conditions

    Args:
        _oil_specific_gravity (double): Oil's specific gravity (no unit).
        _gas_specific_gravity (double): Gas' specific gravity (no unit).
        _gas_solubility_in_oil (double): Gas solubility in oil, Rso
            (:math:`scf/stb`).
        _oil_formation_volume_factor (double): Oil formation volume factor, Bo
            (:math:`bbl/stb`)
        _water_cut: Water cut, WC.

    Returns:
        The live oil density in :math:`lbm/ft^3`.
    """
    density = ((
        350 * _oil_specific_gravity +
        0.0764 * _gas_specific_gravity * _gas_solubility_in_oil *
        (1 - _water_cut)
    ) / (
        5.615 * _oil_formation_volume_factor
    ))
    return density


def live_water_density(_water_specific_gravity,
                       _gas_specific_gravity,
                       _gas_solubility_in_water,
                       _water_formation_volume_factor,
                       _water_cut):
    """
    Calculates the live water density at the given conditions

    Args:
        _water_specific_gravity (double): Water's specific gravity (no unit).
        _gas_specific_gravity (double): Gas' specific gravity (no unit).
        _gas_solubility_in_water (double): Gas solubility in water, Rsw
            (:math:`scf/stb`).
        _water_formation_volume_factor (double): Oil formation volume factor,
            Bw (:math:`bbl/stb`)
        _water_cut: Water cut, WC.

    Returns:
        The live water density in :math:`lbm/ft^3`.
    """
    density = (
        (
            350 * _water_specific_gravity +
            0.0764 * _gas_specific_gravity * _gas_solubility_in_water *
            _water_cut
        ) / (
            5.615 * _water_formation_volume_factor
        )
    )
    return density


def in_situ_oil_flow_rate(_liquid_flow_rate,
                          _oil_formation_volume_factor,
                          _water_cut):
    """
    Calculates the in-situ oil flow rate at the given conditions

    Args:
        _liquid_flow_rate (double): Total liquid flow rate (bpd).
        _oil_formation_volume_factor (double): Oil formation volume factor, Bo
            (:math:`bbl/stb`).
        _water_cut: Water cut, WC.

    Returns:
        The in-situ oil flow rate in :math:`ft^3/s`.
    """
    flow_rate = ((1 - _water_cut) *
                 _liquid_flow_rate * _oil_formation_volume_factor *
                 5.614583 / 86400)
    return flow_rate


def in_situ_water_flow_rate(_liquid_flow_rate,
                            _water_formation_volume_factor,
                            _water_cut):
    """
    Calculates the in-situ water flow rate at the given conditions

    Args:
        _liquid_flow_rate (double): Total liquid flow rate (bpd).
        _water_formation_volume_factor (double): Water formation volume factor,
            Bw (:math:`bbl/stb`).
        _water_cut: Water cut, WC.

    Returns:
        The in-situ water flow rate in :math:`ft^3/s`.
    """
    flow_rate = (_water_cut *
                 _liquid_flow_rate * _water_formation_volume_factor *
                 5.614583 / 86400)
    return flow_rate


def in_situ_gas_flow_rate(_liquid_flow_rate,
                          _gas_formation_volume_factor,
                          _free_gas_liquid_ratio):
    """
    Calculates the in-situ gas flow rate at the given conditions

    Args:
        _liquid_flow_rate (double): Total liquid flow rate (bpd).
        _gas_formation_volume_factor (double): Gas formation volume factor,
            Bg (:math:`bbl/scf`).
        _water_cut: Water cut, WC.

    Returns:
        The in-situ oil flow rate in :math:`ft^3/s`.
    """
    flow_rate = (_free_gas_liquid_ratio *
                 _liquid_flow_rate * _gas_formation_volume_factor *
                 5.614583 / 86400)
    return flow_rate


def superficial_velocity(_in_situ_flow_rate, _diameter):
    """
    Transforms the passed in-situ flow rate into superficial velocity at the
    given diameter.

    Args:
        _in_situ_flow_rate (double): Total liquid flow rate (:math:`ft^3/s`).
        _diameter (double): The tubing diameter (:math:`in`).

    Returns:
        The superficial velocity in :math:`ft/s`.
    """
    return 4 * _in_situ_flow_rate / (math.pi * (_diameter / 12) ** 2)
