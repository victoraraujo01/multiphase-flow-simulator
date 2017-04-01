"""
Correlations test
"""

from src import correlations
from src import formulas
import pytest


@pytest.fixture(scope="module")
def input():
    input_ = {}
    # Physical conditions
    input_["water_cut"] = 0
    input_["pressures"] = [0., 5., 25.5, 100.5, 1000.5, 10088.07, 12515.47]
    input_["temperature"] = 175  # fahrenheit

    # Fluid properties
    input_["oil_api_gravity"] = 25
    input_["oil_specific_gravity"] = formulas.specific_gravity_from_api(
        input_["oil_api_gravity"]
    )
    input_["gas_specific_gravity"] = 0.65
    input_["bubble_point"] = 83.44862

    # Production data
    input_["production_gas_liquid_ratio"] = 10  # scf/stb
    return input_


@pytest.fixture(scope="module")
def correlation_results(input):
    results = {}
    results["rsw"] = []
    results["rso"] = []
    results["bg"] = []
    results["bg_bbl"] = []
    results["co"] = []
    results["bo"] = []
    results["cw"] = []
    results["bw"] = []
    results["z"] = []
    results["uo"] = []
    results["ug"] = []
    results["uw"] = []

    for pressure in input["pressures"]:
        gas_solubility_in_water = correlations.gas_solubility_in_water(
            pressure,
            input["bubble_point"],
            input["temperature"]
        )
        gas_solubility_in_oil = correlations.gas_solubility_in_oil(
            pressure,
            input["bubble_point"],
            input["temperature"],
            input["gas_specific_gravity"],
            input["oil_api_gravity"]
        )
        gas_form_volume_factor = correlations.gas_formation_volume_factor(
            pressure,
            input["temperature"],
            input["gas_specific_gravity"],
            True
        )
        gas_form_volume_factor_bbl = correlations.gas_formation_volume_factor(
            pressure,
            input["temperature"],
            input["gas_specific_gravity"],
            False
        )
        gas_deviation_factor = correlations.gas_deviation_factor(
            pressure,
            input["temperature"],
            input["gas_specific_gravity"]
        )
        live_oil_viscosity = correlations.live_oil_viscosity(
            pressure,
            input["bubble_point"],
            input["temperature"],
            gas_solubility_in_oil,
            input["oil_api_gravity"]
        )
        gas_density = formulas.gas_density(
            input["gas_specific_gravity"], gas_form_volume_factor, True
        )
        gas_viscosity = correlations.gas_viscosity(
            input["temperature"], input["gas_specific_gravity"], gas_density
        )
        water_viscosity = correlations.water_viscosity(
            pressure, input["temperature"]
        )
        oil_compressibility = 0.
        water_compressibility = 0
        if pressure >= input["bubble_point"]:
            oil_compressibility = correlations.oil_compressibility(
                pressure,
                input["bubble_point"],
                input["temperature"],
                gas_solubility_in_oil,
                input["gas_specific_gravity"],
                input["oil_api_gravity"]
            )
            water_compressibility = correlations.water_compressibility(
                pressure,
                input["bubble_point"],
                input["temperature"],
                gas_solubility_in_water
            )
        oil_form_volume_factor = correlations.oil_formation_volume_factor(
            pressure,
            input["bubble_point"],
            input["temperature"],
            gas_solubility_in_oil,
            input["gas_specific_gravity"],
            input["oil_specific_gravity"],
            oil_compressibility
        )
        water_form_volume_factor = correlations.water_formation_volume_factor(
            pressure,
            input["bubble_point"],
            input["temperature"],
            water_compressibility
        )
        results["rsw"].append(gas_solubility_in_water)
        results["rso"].append(gas_solubility_in_oil)
        results["bg"].append(gas_form_volume_factor)
        results["bg_bbl"].append(gas_form_volume_factor_bbl)
        results["co"].append(oil_compressibility)
        results["bo"].append(oil_form_volume_factor)
        results["cw"].append(water_compressibility)
        results["bw"].append(water_form_volume_factor)
        results["z"].append(gas_deviation_factor)
        results["uo"].append(live_oil_viscosity)
        results["ug"].append(gas_viscosity)
        results["uw"].append(water_viscosity)
    return results


@pytest.fixture(scope="module")
def expected_answers():
    answers = {}
    answers["rso"] = [
        2.581821, 2.97368, 4.667221, 9.99999, 9.99999, 9.99999, 9.99999
    ]
    answers["rsw"] = [
        2.224238, 2.248478, 2.347755, 2.627459,
        2.627459, 2.627459, 2.627459
    ]
    answers["co"] = [
        0., 0., 0., 0.000102018231614898, 0.0000115765369065,
        1.16329503591492e-06, 9.37936005636566e-07
    ]
    answers["bo"] = [
        1.05393838926934, 1.05407923616424, 1.05468843676191,
        1.05477542649677, 1.04545390558636, 1.04438595011707,
        1.04436284993531
    ]
    answers["cw"] = [
        0., 0., 0., 3.28222578859436e-06, 3.21086775559541e-06,
        2.49034442014193e-06, 2.29788351995951e-06
    ]
    answers["bw"] = [
        1.026976049, 1.02695939900076, 1.02689113400342,
        1.02664070603277, 1.02367946741586, 1.00143406717806,
        0.997783106257181
    ]
    answers["z"] = [
        0.995264833426288, 0.993657411592812, 0.987083853157865,
        0.963265409927398, 0.705764343983635, 1.03495623866993,
        2.0249883893761
    ]
    answers["bg_bbl"] = [
        0.217716165147761, 0.16219587571978, 0.0789582215069893,
        0.0268882693610031, 0.00223551437955921, 0.000329420752706747,
        0.000519678606823437
    ]
    answers["uo"] = [
        5.59927133208733, 5.57988842602591, 5.49767846161432,
        5.26037907452966, 6.45257947952086, 82.6292113951839,
        104.890675713315
    ]
    answers["ug"] = [
        0.0130207829090467, 0.0130224390730264, 0.0130307467885985,
        0.0130747646161221, 0.0147783550021731, 0.0704594847727007,
        0.032399978006134
    ]
    answers["uw"] = [
        0.334024576050522, 0.334092052578831, 0.334369248728208,
        0.335390812475144, 0.348560031401639, 0.57570480492809,
        0.665379244980299
    ]
    return answers


def test_gas_solubility_in_oil(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["rso"]) == correlation_results["rso"]


def test_gas_solubility_in_water(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["rsw"]) == correlation_results["rsw"]


def test_mixture_bubble_point(input):
    expected_answer = input["bubble_point"]
    answer = correlations.mixture_bubble_point(
        input["temperature"],
        input["gas_specific_gravity"],
        input["oil_api_gravity"],
        input["water_cut"],
        input["production_gas_liquid_ratio"]
    )
    assert answer == pytest.approx(expected_answer)


def test_oil_compressibility(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["co"]) == correlation_results["co"]


def test_oil_formation_volume_factor(input, expected_answers,
                                     correlation_results):
    assert pytest.approx(expected_answers["bo"]) == correlation_results["bo"]


def test_water_compressibility(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["cw"]) == correlation_results["cw"]


def test_water_formation_volume_factor(input, expected_answers, 
                                       correlation_results):
    assert pytest.approx(expected_answers["bw"]) == correlation_results["bw"]


def test_gas_deviation_factor(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["z"]) == correlation_results["z"]


def test_gas_formation_volume_factor(input, expected_answers, 
                                     correlation_results):
    assert (pytest.approx(expected_answers["bg_bbl"]) ==
            correlation_results["bg_bbl"])


def test_dead_oil_viscosity(input):
    expected_answer = 5.7287743853257
    answer = correlations.dead_oil_viscosity(
        input["temperature"],
        input["oil_api_gravity"],
    )
    assert answer == pytest.approx(expected_answer)


def test_oil_viscosity(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["uo"]) == correlation_results["uo"]


def test_gas_viscosity(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["ug"], 1e-1) == correlation_results["ug"]

def test_water_viscosity(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["uw"]) == correlation_results["uw"]
