"""
Formulas
"""
from enum import Enum
import math


class FlowPattern(Enum):
    distributed = 1
    intermittent = 2
    transition = 3
    segregated = 4


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
    volume factor, :math:`B_g`, is given) or uses the given :math:`B_g` to
    calculate it at different conditions.

    Args:
        _gas_specific_gravity (double):  Gas' specific gravity (no unit)
        _gas_formation_volume_factor (double, optional): Gas' formation volume
            factor :math:`B_g` (:math:`ft^3/scf` or :math:`bbl/scf` - use
            argument ``bg_in_cubic_feet`` accordingly).
        _in_cubic_feet (boolean, optional): must be ``True`` if supplied gas
            formation volume factor is in :math:`ft^3/scf` or ``False`` if it's
            in :math:`bbl/scf`. Defaults to ``True``.

    Returns:
        The gas density in :math:`lbm/ft^3`.
    """
    # conversion_factor = 14.7 * 28.97 / (10.7316 * 520)
    conversion_factor = 0.0764106
    if not _bg_in_cubic_feet:
        conversion_factor = conversion_factor / 5.614583333

    density = (conversion_factor * _gas_specific_gravity /
               _gas_formation_volume_factor)
    return density


def dead_oil_density(_oil_specific_gravity, _in_cubic_feet=False):
    """
    Calculates the dead oil density at standard conditions

    Args:
        _oil_specific_gravity (double): Oil's specific gravity (no unit).
        _in_cubic_feet (boolean, optional): if ``True`` density will be in
            :math:`lbm/scf`. Otherwise it will be in :math:`lbm/stb`. Defaults
            to ``False``.
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
        _gas_solubility_in_oil (double): Gas solubility in oil, :math:`R_{so}`
            (:math:`scf/stb`).
        _oil_formation_volume_factor (double): Oil formation volume factor,
            :math:`B_o` (:math:`bbl/stb`)
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
        _gas_solubility_in_water (double): Gas solubility in water,
            :math:`R_{sw}` (:math:`scf/stb`).
        _water_formation_volume_factor (double): Water formation volume factor,
            :math:`B_w` (:math:`bbl/stb`)
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
        _liquid_flow_rate (double): Total liquid flow rate (:math:`bpd`).
        _oil_formation_volume_factor (double): Oil formation volume factor,
            :math:`B_o` (:math:`bbl/stb`).
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
        _liquid_flow_rate (double): Total liquid flow rate (:math:`bpd`).
        _water_formation_volume_factor (double): Water formation volume factor,
            :math:`B_w` (:math:`bbl/stb`).
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
        _liquid_flow_rate (double): Total liquid flow rate (:math:`bpd`).
        _gas_formation_volume_factor (double): Gas formation volume factor,
            :math:`B_g` (:math:`bbl/scf`). Pay attention to the unit used.
        _free_gas_liquid_ratio (double): The free gas liquid ratio
            (:math:`scf/stb`).

    Returns:
        The in-situ gas flow rate in :math:`ft^3/s`.
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
        _in_situ_flow_rate (double): In situ flow rate (:math:`ft^3/s`).
        _diameter (double): The tubing diameter (:math:`in`).

    Returns:
        The superficial velocity in :math:`ft/s`.
    """
    return 4 * _in_situ_flow_rate / (math.pi * (_diameter / 12) ** 2)


def gas_fraction(_oil_velocity, _gas_velocity, _water_velocity):
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
        THe gas fraction.
    """
    total_velocity = _oil_velocity + _gas_velocity + _water_velocity
    return _gas_velocity / total_velocity


def liquid_fraction(_oil_velocity, _gas_velocity, _water_velocity):
    """
    Calculates the no slip liquid fraction (oil + water) of the produced fluid
    based on the superficial velocity of each phase. Note that the suggested
    unit is :math:`ft/s`, but as long as all velocities are in the same unit,
    any unit can be used.

    Args:
        _oil_velocity (double): Superficial oil velocity (:math:`ft/s`).
        _gas_velocity (double): Superficial gas velocity (:math:`ft/s`).
        _water_velocity (double): Superficial water velocity (:math:`ft/s`).

    Returns:
        THe liquid (oil + water) fraction.
    """
    total_velocity = _oil_velocity + _gas_velocity + _water_velocity
    return (_oil_velocity + _water_velocity) / total_velocity


def oil_fraction(_oil_velocity, _water_velocity):
    """
    Calculates the oil fraction of the produced fluid based on the superficial
    velocity of each **liquid** phase. Note that the suggested unit is
    :math:`ft/s`, but as long as all velocities are in the same unit, any unit
    can be used.

    Args:
        _oil_velocity (double): Superficial oil velocity (:math:`ft/s`).
        _water_velocity (double): Superficial water velocity (:math:`ft/s`).

    Returns:
        THe oil fraction.
    """
    total_velocity = _oil_velocity + _water_velocity
    return _oil_velocity / total_velocity


def water_fraction(_oil_velocity, _water_velocity):
    """
    Calculates the water fraction of the produced fluid based on the
    superficial velocity of each **liquid** phase. Note that the suggested unit
    is :math:`ft/s`, but as long as all velocities are in the same unit, any
    unit can be used.

    Args:
        _oil_velocity (double): Superficial oil velocity (:math:`ft/s`).
        _water_velocity (double): Superficial water velocity (:math:`ft/s`).

    Returns:
        THe oil fraction.
    """
    total_velocity = _oil_velocity + _water_velocity
    return _water_velocity / total_velocity


def mixture_froude_number(_mixture_velocity, _diameter):
    """
    Calculates the mixture Froude number based on the mixture velocity (sum of
    all phases respective superficial velocities) and the tubing diameter.

    Args:
        _mixture_velocity (double): Superficial mixture velocity (:math:`ft/s`)
        _diameter (double): Tubing diameter (:math:`in`).

    Returns:
        THe Froude number.
    """
    return 0.37267 * (_mixture_velocity ** 2) / _diameter


def transition_froude_numbers(_liquid_fraction):
    return (
        316.0 * _liquid_fraction ** 0.302,  # Fr_1
        0.0009252 * _liquid_fraction ** -2.4684,  # Fr_2
        0.1 * _liquid_fraction ** -1.4516,  # Fr_3
        0.5 * _liquid_fraction ** -6.738  # Fr_4
    )


def flow_pattern(_froude_number, _liquid_fraction):
    fr1, fr2, fr3, fr4 = transition_froude_numbers(_liquid_fraction)
    if _froude_number > fr1 or _froude_number > fr4:
        return FlowPattern.distributed
    elif _froude_number > fr3:
        return FlowPattern.intermittent
    elif _froude_number > fr2:
        return FlowPattern.transition
    else:
        return FlowPattern.segregated


def horizontal_liquid_phase_fraction(_flow_pattern,
                                     _froude_number,
                                     _no_slip_liquid_fraction):
    constants = {
        FlowPattern.segregated:   (0.980, 0.4846, 0.0868),
        FlowPattern.intermittent: (0.845, 0.5351, 0.0173),
        FlowPattern.distributed:  (1.065, 0.5824, 0.0609)
    }
    if _flow_pattern != FlowPattern.transition:
        term_a, term_b, term_c = constants[_flow_pattern]
        answer = (
            term_a *
            _no_slip_liquid_fraction ** term_b / _froude_number ** term_c
        )
        return max(answer, _no_slip_liquid_fraction)
    else:
        _, fr2, fr3, _ = transition_froude_numbers(_no_slip_liquid_fraction)
        term_a = (fr3 - _froude_number) / (fr3 - fr2)
        term_b = 1 - term_a
        answer = (
            (
                term_a * horizontal_liquid_phase_fraction(
                    FlowPattern.segregated,
                    _froude_number,
                    _no_slip_liquid_fraction
                )
            ) + (
                term_b * horizontal_liquid_phase_fraction(
                    FlowPattern.intermittent,
                    _froude_number,
                    _no_slip_liquid_fraction
                )
            )
        )
        return answer
