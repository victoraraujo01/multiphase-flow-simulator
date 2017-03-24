"""
Correlations
"""
import math


def gas_solubility_in_oil(_pressure,
                          _bubble_point,
                          _temperature,
                          _gas_specific_gravity,
                          _oil_api_gravity):
    """
    Calculates gas solubility in oil (Rso) using Standing correlation. If
    pressure is higher than the mixture's bubble point, returns the Rso at
    bubble point.

    Args:
        _gas_specific_gravity: Gas' specific gravity (doesn't have an unit).
        _pressure: Pressure at which the gas is (psig). Note that this value is
                   in psig, so it is relative to the atmospheric pressure.
        _bubble_point: Mixture's bubble point (psig).  Note that this value is
                       in psig, so it is relative to the atmospheric pressure.
        _temperature: Temperature at which the gas is (fahrenheit degrees).
        _oil_api_gravity: Oil's API gravity (API degrees).

    Returns:
        The gas solubility in oil, Rso (scf/stb).
    """
    if _pressure > _bubble_point:
        _pressure = _bubble_point

    exponent = 0.0125 * _oil_api_gravity - 0.00091 * _temperature
    first_term = (_pressure + 14.7)/18.2 + 1.4
    return _gas_specific_gravity * (first_term * 10 ** exponent) ** 1.2048


def gas_solubility_in_water(_pressure, _bubble_point, _temperature):
    """
    Calculates gas solubility in water (Rsw) using Culberson and Maketta
    correlation. If pressure is higher than the mixture's bubble point, returns
    the Rso at bubble point.

    Args:
        _pressure: Pressure at which the gas is (psig). Note that this value is
                   in psig, so it is relative to the atmospheric pressure.
        _bubble_point: Mixture's bubble point (psig). Note that this value is
                       in psig, so it is relative to the atmospheric pressure.
        _temperature: Temperature at which the gas is (fahrenheit degrees).

    Returns:
        The gas solubility in water, Rsw (scf/stb).
    """
    term_a = (8.15839 -
              6.12265e-2 * _temperature +
              1.91663e-4 * (_temperature ** 2) -
              2.1654e-7 * (_temperature ** 3))

    term_b = (1.01021e-2 -
              7.44241e-5 * _temperature +
              3.05553e-7 * (_temperature ** 2) -
              2.94883e-10 * (_temperature ** 3))

    term_c = (-9.02505 +
              0.130237 * _temperature -
              8.53425e-4 * (_temperature ** 2) +
              2.34122e-6 * (_temperature ** 3) -
              2.37049e-9 * (_temperature ** 4)) * (10 ** -7)

    if _pressure > _bubble_point:
        _pressure = _bubble_point

    abs_pressure = _pressure + 14.7
    return term_a + term_b * (abs_pressure) + term_c * (abs_pressure) ** 2


def mixture_bubble_point(_temperature,
                         _gas_specific_gravity,
                         _oil_api_gravity,
                         _water_cut,
                         _production_gas_liquid_ratio):
    """
    Calculates the mixture's bubble point based on the water cut and the
    prouction gas liquid ratio.

    Args:
        _water_cut: Water cut, WC.
        _production_gas_liquid_ratio: Production gas liquid ratio, GLR_p (in
                                      the same unit as Rso and Rsw, suggestion:
                                      scf/stb).

    Returns:
        The mixture's bubble point Pb (psi).
    """
    pressure_low = 0.0
    pressure_high = 100000.0
    bubble_point = 0.0
    error = 1.0
    while abs(error) > 1e-10:
        bubble_point = (pressure_low + pressure_high)/2
        rso = gas_solubility_in_oil(bubble_point,
                                    pressure_high,
                                    _temperature,
                                    _gas_specific_gravity,
                                    _oil_api_gravity)
        rsw = gas_solubility_in_water(bubble_point,
                                      pressure_high,
                                      _temperature)

        error = (_production_gas_liquid_ratio -
                 (1 - _water_cut) * rso - _water_cut * rsw)

        if error > 0.0:
            pressure_low = bubble_point
        else:
            pressure_high = bubble_point

    return bubble_point


def oil_compressibility(_pressure,
                        _bubble_point,
                        _temperature,
                        _gas_solubility_in_oil_at_bp,
                        _gas_specific_gravity,
                        _oil_api_gravity):
    """
    Calculates the isothermal oil compressibility using Vasquez correlation.
    This is the compressibility of the oil as a single-phase liquid with a
    certain amount of gas in solution. It is valid above the bubble point
    only (Vasquez-Beggs Correlation).

    Args:
        _pressure: Pressure at which the oil is (psig). Note that this value is
                   in psig, so it is relative to the atmospheric pressure, and
                   must be above bubble point.
        _bubble_point: Mixture's bubble point (psig). Note that this value is
                       in psig, so it is relative to the atmospheric pressure.
        _temperature: Temperature (fahrenheit degrees).
        _gas_solubility_in_oil_at_bp: Gas solubility in oil at bubble point,
                                      Rsob (in scf/stb).
        _gas_specific_gravity: Gas' specific gravity (doesn't have an unit).
        _oil_api_gravity: Oil's API gravity (API degrees).

    Returns:
        The compressibility of the oil as a single-phase liquid with a certain
        amount of gas in solution in psi-1.
    """
    if _pressure < _bubble_point:
        raise ValueError('Pressure must be below bubble point.')

    numerator = (-1433 +
                 5 * _gas_solubility_in_oil_at_bp +
                 17.2 * _temperature -
                 1180 * _gas_specific_gravity +
                 12.61 * _oil_api_gravity)
    denominator = (_pressure + 14.7) * (10 ** 5)
    return numerator/denominator


def oil_formation_volume_factor(_pressure,
                                _bubble_point,
                                _temperature,
                                _gas_solubility_in_oil,
                                _gas_specific_gravity,
                                _oil_specific_gravity,
                                _oil_compressibility):
    """
    Calculates the oil formation volume factor (B_o) using Standing's
    correlation. The bubble point is necessary because a different correlation
    must be used below and above it. Also, it is important to remember that if
    pressure is above bubble point, the gas solubility in oil supplied must be
    the one at the bubble point (R_sob). Finally, the oil compressibility (c_o)
    is only needed if pressure is above bubble point.

    Args:
        _pressure: Pressure at which the oil is (psig). Note that this value is
                   in psig, so it is relative to the atmospheric pressure.
        _bubble_point: Mixture's bubble point (psig). Note that this value is
                       in psig, so it is relative to the atmospheric pressure.
        _temperature: Temperature (fahrenheit degrees).
        _gas_solubility_in_oil_: Gas solubility in oil, Rso (in scf/stb).
        _gas_specific_gravity: Gas' specific gravity (doesn't have an unit).
        _oil_specific_gravity: Oil's specific gravity (doesn't have an unit).
        _oil_api_gravity: Oil's API gravity (API degrees).
        _oil_compressibility: Oil's compressibility (psi-1). Value doesn't
                              matter if pressure is below bubble point.
    Returns:
        The oil formation volume factor, in bbl/stb.
    """
    result = (0.9759 + 12e-5 * (
        _gas_solubility_in_oil *
        math.sqrt(_gas_specific_gravity/_oil_specific_gravity) +
        1.25 * _temperature
    ) ** 1.2)

    if _pressure > _bubble_point:
        result = (result *
                  math.exp(_oil_compressibility * (_bubble_point - _pressure)))

    return result


def water_compressibility(_pressure,
                          _bubble_point,
                          _temperature,
                          _gas_solubility_in_water_at_bp):
    """
    Calculates the isothermal water compressibility using Dodson and Standing
    correlation. This is the compressibility of the water as a single-phase
    liquid with a certain amount of gas in solution. It is valid above the
    bubble point only.

    Args:
        _pressure: Pressure at which the water is (psig). Note that this value
                   is in psig, so it is relative to the atmospheric pressure,
                   and must be above bubble point.
        _bubble_point: Mixture's bubble point (psig). Note that this value is
                       in psig, so it is relative to the atmospheric pressure.
        _temperature: Temperature (fahrenheit degrees).
        _gas_solubility_in_water_at_bp: Gas solubility in water at bubble point
                                        Rswb (in scf/stb).

    Returns:
        The compressibility of the water as a single-phase liquid with a
        certain amount of gas in solution in psi-1.
    """
    if _pressure < _bubble_point:
        raise ValueError('Pressure must be below bubble point.')

    term_a = 3.8546 - 1.34e-4 * (_pressure + 14.7)
    term_b = -0.01052 + 4.77e-7 * (_pressure + 14.7)
    term_c = 3.9267e-5 - 8.8e-10 * (_pressure + 14.7)

    result = ((term_a + term_b * _temperature + term_c * (_temperature ** 2)) *
              (1 + 8.9e-3 * _gas_solubility_in_water_at_bp) / 1e6)
    return result


def water_formation_volume_factor(_pressure,
                                  _bubble_point,
                                  _temperature,
                                  _water_compressibility):
    """
    Calculates the water formation volume factor (B_w) using Gould's
    correlation. The bubble point is necessary because a different correlation
    must be used below and above it and the water compressibility (c_w)
    is only needed if pressure is above bubble point.

    Args:
        _pressure: Pressure at which the oil is (psig). Note that this value is
                   in psig, so it is relative to the atmospheric pressure.
        _bubble_point: Mixture's bubble point (psig). Note that this value is
                       in psig, so it is relative to the atmospheric pressure.
        _temperature: Temperature (fahrenheit degrees).
        _water_compressibility: Waters compressibility (psi-1). Value doesn't
                                matter if pressure is below bubble point.
    Returns:
        The water formation volume factor, in bbl/stb.
    """
    result = (1.0 +
              1.2e-4 * (_temperature - 60) +
              1.0e-6 * (_temperature - 60) ** 2)

    if _pressure >= _bubble_point:
        result = result - 3.33e-6 * (_bubble_point + 14.7)
        result = (result * math.exp(_water_compressibility *
                                    (_bubble_point - _pressure)))
    else:
        result = result - 3.33e-6 * (_pressure + 14.7)

    return result


def gas_deviation_factor(_pressure,
                         _temperature,
                         _gas_specific_gravity):
    """
    Calculates the gas deviation factor or compressibility factor Z using Papay
    correlation.

    Args:
        _pressure: Pressure at which the gas is (psig). Note that this value is
                   in psig, so it is relative to the atmospheric pressure.
        _temperature: Temperature (fahrenheit degrees).
        _gas_specific_gravity: Gas' specific gravity (doesn't have an unit).
    Returns:
        The gas deviation factor.
    """
    pseudo_critical_temperature = (168. +
                                   325. * _gas_specific_gravity -
                                   12.5 * _gas_specific_gravity ** 2)
    pseudo_critical_pressure = (677. +
                                15.0 * _gas_specific_gravity -
                                37.5 * _gas_specific_gravity ** 2)

    pseudo_reduced_temperature = ((_temperature + 460) /
                                  pseudo_critical_temperature)

    pseudo_reduced_pressure = ((_pressure + 14.7) /
                               pseudo_critical_pressure)

    pseudo_reduced_ratio = pseudo_reduced_pressure/pseudo_reduced_temperature
    deviation_factor = (1 - pseudo_reduced_ratio *
                        (0.3675 - 0.04188423 * pseudo_reduced_ratio))
    return deviation_factor
