"""
Correlations
TODO: should I check if pressure is above bubble point inside methods?
"""


def gas_solubility_in_oil(pressure,
                          temperature,
                          gas_specific_gravity,
                          oil_api_gravity):
    """
    Calculates gas solubility in oil (Rso) using Standing correlation.

    Args:
        gas_specific_gravity: Gas' specific gravity (doesn't have an unit)
        pressure: Pressure at which the gas is (psig).
        temperature: Temperature at which the gas is (fahrenheit degrees)
        oil_api_gravity: Oil's API gravity (API degrees)

    Returns:
        The gas solubility in oil, Rso (scf/stb)
    """
    exponent = 0.0125 * oil_api_gravity - 0.00091 * temperature
    first_term = (pressure + 14.7)/18.2 + 1.4
    return gas_specific_gravity * (first_term * 10 ** exponent) ** 1.2048


def gas_solubility_in_water(pressure, temperature):
    """
    Calculates gas solubility in water (Rsw) using Culberson and Maketta
    correlation.

    Args:
        pressure: Pressure at which the gas is (psig).
        temperature: Temperature at which the gas is (fahrenheit degrees)

    Returns:
        The gas solubility in water, Rsw (scf/stb)
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

    abs_pressure = pressure + 14.7
    return term_a + term_b * (abs_pressure) + term_c * (abs_pressure) ** 2


def isothermal_oil_compressibility(_pressure,
                                   _temperature,
                                   _gas_solubility_in_oil_at_bp,
                                   _gas_specific_gravity,
                                   _oil_api_gravity):
    """
    Calculates the isothermal oil compressibility using Vasquez correlation.
    This is the compressibility of the oil as a single-phase liquid with a
    certain amount of gas in solution. It is valid above the bubble point
    only (Vasquez- Beggs Correlation).

    Args:
        _pressure: Pressure (psig)
        _temperature: Temperature (fahrenheit degrees)
        _gas_solubility_in_oil_at_bp: Gas solubility in oil at bubble point,
                                      Rsob (in scf/stb)
        _gas_specific_gravity: Gas' specific gravity (doesn't have an unit)
        _oil_api_gravity: Oil's API gravity (API degrees)

    Returns:
        The compressibility of the oil as a single-phase liquid with a certain
        amount of gas in solution in psi-1.
    """
    numerator = (-1433 +
                 5 * _gas_solubility_in_oil_at_bp +
                 17.2 * _temperature -
                 1180 * _gas_specific_gravity +
                 12.61 * _oil_api_gravity)
    denominator = (_pressure + 14.7) * (10 ** 5)
    return numerator/denominator