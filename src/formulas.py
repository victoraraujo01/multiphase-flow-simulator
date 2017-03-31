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
    downward = 5


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


def no_slip_gas_fraction(_oil_velocity, _gas_velocity, _water_velocity):
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
    total_velocity = _oil_velocity + _gas_velocity + _water_velocity
    return _gas_velocity / total_velocity


def no_slip_liquid_fraction(_oil_velocity, _gas_velocity, _water_velocity):
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
        THe no slip liquid (oil + water) fraction.
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
        THe water fraction.
    """
    total_velocity = _oil_velocity + _water_velocity
    return _water_velocity / total_velocity


def froude_number(_mixture_velocity, _diameter):
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


def transition_froude_numbers(_no_slip_liquid_fraction):
    """
    Calculates the Froude number limits used to determine flow pattern.

    Args:
        _no_slip_liquid_fraction (double): The no slip liquid fraction

    Returns:
        A tuple with the four froude number limits in the format
        :math:`(Fr_1, Fr_2, Fr_3 and Fr_4)`.
    """
    return (
        316.0 * _no_slip_liquid_fraction ** 0.302,  # Fr_1
        0.0009252 * _no_slip_liquid_fraction ** -2.4684,  # Fr_2
        0.1 * _no_slip_liquid_fraction ** -1.4516,  # Fr_3
        0.5 * _no_slip_liquid_fraction ** -6.738  # Fr_4
    )


def flow_pattern(_froude_number, _no_slip_liquid_fraction):
    """
    Returns the flow pattern based no the Froude number and the no slip liquid
    fraction.

    Args:
        _froude_number (double): The mixture's froude number
        _no_slip_liquid_fraction (double): The no slip liquid fraction

    Returns:
        The flow pattern using the `FlowPattern` enum.
    """
    fr1, fr2, fr3, fr4 = transition_froude_numbers(_no_slip_liquid_fraction)
    if _froude_number > fr1 or _froude_number > fr4:
        return FlowPattern.distributed
    elif _froude_number > fr3:
        return FlowPattern.intermittent
    elif _froude_number > fr2:
        return FlowPattern.transition
    else:
        return FlowPattern.segregated


def horz_liquid_holdup(_flow_pattern,
                       _froude_number,
                       _no_slip_liquid_fraction):
    """
    Returns the liquid fraction considering slippage for horizontal flow.

    Args:
        _flow_pattern (FlowPattern): The flow pattern as determined by
            `flow_pattern`.
        _froude_number (double): The mixture's froude number
        _no_slip_liquid_fraction (double): The no slip liquid fraction

    Returns:
        The liquid fraction considering slippage for horizontal flow.
    """
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
                term_a * horz_liquid_holdup(
                    FlowPattern.segregated,
                    _froude_number,
                    _no_slip_liquid_fraction
                )
            ) + (
                term_b * horz_liquid_holdup(
                    FlowPattern.intermittent,
                    _froude_number,
                    _no_slip_liquid_fraction
                )
            )
        )
        return max(answer, _no_slip_liquid_fraction)


def estimate_liquid_property(_oil_property, _water_property, _water_fraction):
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
        _oil_property * (1 - _water_fraction) +
        _water_property * _water_fraction
    )


def liquid_velocity_number(_liquid_superficial_velocity,
                           _liquid_density,
                           _superficial_tension):
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
        1.938 * _liquid_superficial_velocity *
        (_liquid_density / _superficial_tension) ** (1/4)
    )


def liquid_holdup_with_incl(_flow_pattern,
                            _froude_number,
                            _no_slip_liquid_fraction,
                            _liquid_velocity_number,
                            _inclination):
    """
    Returns the liquid fraction considering slippage for any inclination.

    Args:
        _flow_pattern (FlowPattern): The flow pattern as determined by
            `flow_pattern`.
        _froude_number (double): The mixture's froude number.
        _no_slip_liquid_fraction (double): The no slip liquid fraction.
        _liquid_velocity_number (double): The liquid velocity number obtained
            using `liquid_velocity_number()`.
        _inclination (double): The inclination angle with the horizontal in
            degrees.

    Returns:
        The liquid fraction considering slippage for any inclination.
    """
    _horz_liquid_holdup = horz_liquid_holdup(
        _flow_pattern,
        _froude_number,
        _no_slip_liquid_fraction
    )

    constants = {
        FlowPattern.segregated:   (0.011, -3.7680, 3.5390, -1.6140),
        FlowPattern.intermittent: (2.960, 0.3050, -0.4473, 0.0978),
        FlowPattern.distributed:  (0.0, 0.0, 0.0, 0.0),
        FlowPattern.downward:     (4.700, -0.3692, 0.1244, -0.5056)
    }

    term_d, term_e, term_f, term_g = constants[_flow_pattern]

    c_parameter = max(0, (
        (1 - _no_slip_liquid_fraction) *
        math.log(
            term_d *
            _no_slip_liquid_fraction ** term_e *
            _liquid_velocity_number ** term_f *
            _froude_number ** term_g
        )
    ))
    if _flow_pattern == FlowPattern.distributed:
        c_parameter = 0

    phi_parameter = (
        1 + c_parameter * (
            math.sin(1.8 * _inclination) -
            0.333 * math.sin(1.8 * _inclination) ** 3
        )
    )

    return _horz_liquid_holdup * phi_parameter


def gravitational_pressure_gradient(_liquid_density,
                                    _gas_density,
                                    _liquid_fraction,
                                    _inclination):
    """
    Returns the pressure gradient due to gravity.

    Args:
        _liquid_density (double): The liquid density in :math:`lbm/ft^3`.
        _gas_density (double): The gas density in :math:`lbm/ft^3`.
        _liquid_fraction (double): The liquid fraction considering slippage
            and the given inclination.
        _inclination (double): The inclination angle with the horizontal in
            degrees.

    Returns:
        The pressure gradient due to gravity :math:`\frac{dP}{dz}` in
        :math:`psi/ft`.
    """
    _mixture_specific_gravity = (
        (
            _liquid_density * _liquid_fraction +
            _gas_density * (1 - _liquid_fraction)
        ) / 62.4
    )
    _inclination_rad = _inclination * math.pi/180
    return -0.433 * _mixture_specific_gravity * math.sin(_inclination_rad)


def reynolds(_density, _velocity, _diameter, _viscosity):
    return 124 * _density * abs(_velocity) * _diameter / _viscosity


def moody_friction_factor(_reynolds, _rugosity):
    term_a = (
        2.457 * math.log(
            1 / (
                (7/_reynolds) ** 0.9 + 0.27 * _rugosity
            )
        )
    ) ** 16
    term_b = (37530 / _reynolds) ** 16
    answer = (
        8 * (
            (8 / _reynolds) ** 12 +
            1 / ((term_a + term_b) ** (3/2))
        ) ** (1/12)
    )
    return answer


def friction_factor(_no_slip_liquid_fraction,
                    _liquid_holdup,
                    _moody_friction_factor):
    term_y = _no_slip_liquid_fraction / (_liquid_holdup ** 2)
    term_s = math.log(2.2 * term_y - 1.2)
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
    return _moody_friction_factor * math.exp(term_s)
