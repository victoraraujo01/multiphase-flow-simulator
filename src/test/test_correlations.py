"""
Correlations test
"""

from src import correlations
from src import formulas
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
    oil_specific_gravity = formulas.specific_gravity_from_api(oil_api_gravity)
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

    def test_oil_compressibility(self):
        answers = []
        expected_answers = [
            0., 0., 0., 0.000102018231614898, 0.0000115765369065,
            1.16329503591492e-06, 9.37936005636566e-07
        ]
        rsob = correlations.gas_solubility_in_oil(
            self.bubble_point,
            self.bubble_point,
            self.temperature,
            self.gas_specific_gravity,
            self.oil_api_gravity
        )
        for pressure in self.pressures:
            answer = 0.
            if pressure >= self.bubble_point:
                answer = correlations.oil_compressibility(
                    pressure,
                    self.bubble_point,
                    self.temperature,
                    rsob,
                    self.gas_specific_gravity,
                    self.oil_api_gravity
                )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_oil_formation_volume_factor(self):
        answers = []
        expected_answers = [
            1.05393838926934, 1.05407923616424, 1.05468843676191,
            1.05477542649677, 1.04545390558636, 1.04438595011707,
            1.04436284993531
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
            answer = correlations.oil_formation_volume_factor(
                pressure,
                self.bubble_point,
                self.temperature,
                gas_solubility_in_oil,
                self.gas_specific_gravity,
                self.oil_specific_gravity,
                oil_compressibility
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_water_compressibility(self):
        answers = []
        expected_answers = [
            0., 0., 0., 3.28222578859436e-06, 3.21086775559541e-06,
            2.49034442014193e-06, 2.29788351995951e-06
        ]
        rswb = correlations.gas_solubility_in_water(
            self.bubble_point,
            self.bubble_point,
            self.temperature,
        )
        for pressure in self.pressures:
            answer = 0.
            if pressure >= self.bubble_point:
                answer = correlations.water_compressibility(
                    pressure,
                    self.bubble_point,
                    self.temperature,
                    rswb
                )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_water_formation_volume_factor(self):
        answers = []
        expected_answers = [
            1.026976049, 1.02695939900076, 1.02689113400342,
            1.02664070603277, 1.02367946741586, 1.00143406717806,
            0.997783106257181
        ]
        for pressure in self.pressures:
            water_compressibility = 0
            rswb = correlations.gas_solubility_in_water(
                self.bubble_point,
                self.bubble_point,
                self.temperature,
            )
            if pressure >= self.bubble_point:
                water_compressibility = correlations.water_compressibility(
                    pressure,
                    self.bubble_point,
                    self.temperature,
                    rswb
                )
            answer = correlations.water_formation_volume_factor(
                pressure,
                self.bubble_point,
                self.temperature,
                water_compressibility
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_gas_deviation_factor(self):
        answers = []
        expected_answers = [
            0.995264833426288, 0.993657411592812, 0.987083853157865,
            0.963265409927398, 0.705764343983635, 1.03495623866993,
            2.0249883893761
        ]
        for pressure in self.pressures:
            answer = correlations.gas_deviation_factor(
                pressure,
                self.temperature,
                self.gas_specific_gravity
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_gas_formation_volume_factor(self):
        answers = []
        expected_answers = [
            0.217716165147761, 0.16219587571978, 0.0789582215069893,
            0.0268882693610031, 0.00223551437955921, 0.000329420752706747,
            0.000519678606823437
        ]
        for pressure in self.pressures:
            answer = correlations.gas_formation_volume_factor(
                pressure,
                self.temperature,
                self.gas_specific_gravity,
                False
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)

    def test_dead_oil_viscosity(self):
        expected_answer = 5.7287743853257
        answer = correlations.dead_oil_viscosity(
            self.temperature,
            self.oil_api_gravity,
        )
        assert answer == pytest.approx(expected_answer)

    def test_oil_viscosity(self):
        answers = []
        expected_answers = [
            5.59927133208733, 5.57988842602591, 5.49767846161432,
            5.26037907452966, 6.45257947952086, 82.6292113951839,
            104.890675713315
        ]
        for pressure in self.pressures:
            gas_solubility_in_oil = correlations.gas_solubility_in_oil(
                pressure,
                self.bubble_point,
                self.temperature,
                self.gas_specific_gravity,
                self.oil_api_gravity
            )
            answer = correlations.live_oil_viscosity(
                pressure,
                self.bubble_point,
                self.temperature,
                gas_solubility_in_oil,
                self.oil_api_gravity
            )
            answers.append(answer)
        assert answers == pytest.approx(expected_answers)
