"""
Correlations test
"""

from src import correlations
import pytest


class TestCorrelations(object):
    """
    Tests for correlations
    """
    # Physical conditions
    production_gas_liquid_ratio = 10  # scf/stb
    water_cut = 0
    pressures = [0., 5., 25.5, 100.5, 1000.5, 10088.07, 12515.47]
    temperature = 175  # fahrenheit

    # Fluid properties
    oil_api_gravity = 25
    gas_specific_gravity = 0.65
    bubble_point = 83.44862

    def test_gas_solubility_in_oil(self):
        answers = []
        expected_answers = [
            2.581821, 2.97368, 4.667221, 9.99999, 9.99999, 9.99999, 9.99999
        ]
        for pressure in self.pressures:
            answer = correlations.gas_solubility_in_oil(
                pressure,
                self.bubble_point,
                self.temperature,
                self.gas_specific_gravity,
                self.oil_api_gravity
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_gas_solubility_in_water(self):
        answers = []
        expected_answers = [
            2.224238, 2.248478, 2.347755, 2.627459,
            2.627459, 2.627459, 2.627459
        ]
        for pressure in self.pressures:
            answer = correlations.gas_solubility_in_water(
                pressure,
                self.bubble_point,
                self.temperature,
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_mixture_bubble_point(self):
        expected_answer = 83.44862
        answer = correlations.mixture_bubble_point(
            self.temperature,
            self.gas_specific_gravity,
            self.oil_api_gravity,
            self.water_cut,
            self.production_gas_liquid_ratio
        )
        assert answer == pytest.approx(expected_answer)

    def test_isothermal_oil_compressibility(self):
        expected_answer = 0.0001197419
        pressure = self.bubble_point
        rsob = correlations.gas_solubility_in_oil(
            pressure,
            self.bubble_point,
            self.temperature,
            self.gas_specific_gravity,
            self.oil_api_gravity
        )
        answer = correlations.isothermal_oil_compressibility(
            pressure,
            self.bubble_point,
            self.temperature,
            rsob,
            self.gas_specific_gravity,
            self.oil_api_gravity
        )
        assert answer == pytest.approx(expected_answer, 1e-6)
