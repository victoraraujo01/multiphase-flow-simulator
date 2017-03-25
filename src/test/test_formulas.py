"""
Correlations test
"""

import pytest
from src import formulas
from src import correlations


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
    oil_specific_gravity = formulas.specific_gravity_from_api(oil_api_gravity)
    gas_specific_gravity = 0.65
    water_specific_gravity = 1.07
    bubble_point = 83.44862

    def test_free_gas_liquid_ratio(self):
        answers = []
        expected_answers = [
            7.41817908307171, 7.02632049562826, 5.33277871673351,
            0.0, 0.0, 0.0, 0.0
        ]
        for pressure in self.pressures:
            gas_solubility_in_water = correlations.gas_solubility_in_water(
                pressure,
                self.bubble_point,
                self.temperature
            )
            gas_solubility_in_oil = correlations.gas_solubility_in_oil(
                pressure,
                self.bubble_point,
                self.temperature,
                self.gas_specific_gravity,
                self.oil_api_gravity
            )
            answer = formulas.free_gas_liquid_ratio(
                pressure,
                self.bubble_point,
                gas_solubility_in_oil,
                gas_solubility_in_water,
                self.water_cut,
                self.production_gas_liquid_ratio
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_gas_density(self):
        answers = []
        expected_answers = [
            0.0406297014274935, 0.0545374088375396, 0.112030674159871,
            0.32898148509036, 3.95691607568035, 26.8524150746719,
            17.0215642317133
        ]
        bg_in_cubic_feet = True
        for pressure in self.pressures:
            gas_form_volume_factor = correlations.gas_formation_volume_factor(
                pressure,
                self.temperature,
                self.gas_specific_gravity,
                bg_in_cubic_feet
            )
            answer = formulas.gas_density(
                self.gas_specific_gravity,
                gas_form_volume_factor,
                bg_in_cubic_feet
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_live_oil_density(self):
        answers = []
        expected_answers = [
            53.4959715812803, 53.4921112656614, 53.4754149044103,
            53.5157193774093, 53.9928785276566, 54.0480899080085,
            54.0492853935634
        ]
        for pressure in self.pressures:
            oil_compressibility = 0
            gas_solubility_in_oil = correlations.gas_solubility_in_oil(
                pressure,
                self.bubble_point,
                self.temperature,
                self.gas_specific_gravity,
                self.oil_api_gravity
            )
            if pressure >= self.bubble_point:
                oil_compressibility = correlations.oil_compressibility(
                    pressure,
                    self.bubble_point,
                    self.temperature,
                    gas_solubility_in_oil,
                    self.gas_specific_gravity,
                    self.oil_api_gravity
                )
            oil_form_volume_factor = correlations.oil_formation_volume_factor(
                pressure,
                self.bubble_point,
                self.temperature,
                gas_solubility_in_oil,
                self.gas_specific_gravity,
                self.oil_specific_gravity,
                oil_compressibility
            )
            answer = formulas.live_oil_density(
                self.oil_specific_gravity,
                self.gas_specific_gravity,
                gas_solubility_in_oil,
                oil_form_volume_factor,
                self.water_cut
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_live_water_density(self):
        answers = []
        expected_answers = [
            64.9444055973349, 64.9454585350703, 64.9497759367957,
            64.9656191042121, 65.1535477539374, 66.6008389877806,
            66.8445362992679
        ]
        for pressure in self.pressures:
            water_compressibility = 0
            gas_solubility_in_water = correlations.gas_solubility_in_water(
                pressure,
                self.bubble_point,
                self.temperature
            )
            if pressure >= self.bubble_point:
                water_compressibility = correlations.water_compressibility(
                    pressure,
                    self.bubble_point,
                    self.temperature,
                    gas_solubility_in_water
                )
            water_form_vol_factor = correlations.water_formation_volume_factor(
                pressure,
                self.bubble_point,
                self.temperature,
                water_compressibility
            )
            answer = formulas.live_water_density(
                self.water_specific_gravity,
                self.gas_specific_gravity,
                gas_solubility_in_water,
                water_form_vol_factor,
                self.water_cut
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)
