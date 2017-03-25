"""
Correlations test
"""

import pytest
from src import formulas
from src import correlations


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
    input_["water_specific_gravity"] = 1.07
    input_["bubble_point"] = 83.44862

    # Production data
    input_["production_gas_liquid_ratio"] = 10  # scf/stb
    input_["liquid_flow_rate"] = 600  # stb/day

    # Well data
    input_["diameter"] = 1.995  # in
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
    return results


@pytest.fixture(scope="module")
def expected_answers():
    answers = {}
    answers["expected_fglr"] = [
        7.41817908307171, 7.02632049562826,
        5.33277871673351, 0.0, 0.0, 0.0, 0.0
    ]
    answers["expected_rho_gas"] = [
        0.0406297014274935, 0.0545374088375396, 0.112030674159871,
        0.32898148509036, 3.95691607568035, 26.8524150746719,
        17.0215642317133
    ]
    answers["expected_rho_oil"] = [
        53.4959715812803, 53.4921112656614, 53.4754149044103,
        53.5157193774093, 53.9928785276566, 54.0480899080085,
        54.0492853935634
    ]
    answers["expected_rho_water"] = [
        64.9444055973349, 64.9454585350703, 64.9497759367957,
        64.9656191042121, 65.1535477539374, 66.6008389877806,
        66.8445362992679
    ]

    answers["expected_qo"] = [
        0.0410932261349932, 0.0410987177779217, 0.0411224706065278,
        0.0411258623501842, 0.0407624147610332, 0.0407207750067095,
        0.0407198743269328
    ]
    answers["expected_qg"] = [
        0.0629713499770325, 0.0444347536529818, 0.0164174587898616,
        0.0, 0.0, 0.0, 0.0
    ]
    return answers


def test_free_gas_liquid_ratio(input, expected_answers, correlation_results):
    answers = []
    expected = expected_answers["expected_fglr"]
    for (i, pressure) in enumerate(input["pressures"]):
        gas_solubility_in_water = correlation_results["rsw"][i]
        gas_solubility_in_oil = correlation_results["rso"][i]
        answer = formulas.free_gas_liquid_ratio(
            pressure,
            input["bubble_point"],
            gas_solubility_in_oil,
            gas_solubility_in_water,
            input["water_cut"],
            input["production_gas_liquid_ratio"]
        )
        answers.append(answer)
    assert answers == pytest.approx(expected)


def test_gas_density(input, expected_answers, correlation_results):
    answers = []
    expected = expected_answers["expected_rho_gas"]
    for (i, _) in enumerate(input["pressures"]):
        gas_form_volume_factor = correlation_results["bg"][i]
        answer = formulas.gas_density(
            input["gas_specific_gravity"],
            gas_form_volume_factor)
        answers.append(answer)
    assert answers == pytest.approx(expected)


def test_live_oil_density(input, expected_answers, correlation_results):
    answers = []
    expected = expected_answers["expected_rho_oil"]
    for (i, _) in enumerate(input["pressures"]):
        gas_solubility_in_oil = correlation_results["rso"][i]
        oil_form_volume_factor = correlation_results["bo"][i]
        answer = formulas.live_oil_density(
            input["oil_specific_gravity"],
            input["gas_specific_gravity"],
            gas_solubility_in_oil,
            oil_form_volume_factor,
            input["water_cut"]
        )
        answers.append(answer)
    assert answers == pytest.approx(expected)


def test_live_water_density(input, expected_answers, correlation_results):
    answers = []
    expected = expected_answers["expected_rho_water"]
    for (i, _) in enumerate(input["pressures"]):
        gas_solubility_in_water = correlation_results["rsw"][i]
        water_form_vol_factor = correlation_results["bw"][i]
        answer = formulas.live_water_density(
            input["water_specific_gravity"],
            input["gas_specific_gravity"],
            gas_solubility_in_water,
            water_form_vol_factor,
            input["water_cut"]
        )
        answers.append(answer)
    assert answers == pytest.approx(expected)


def test_in_situ_oil_flow_rate(input, expected_answers, correlation_results):
    answers = []
    expected = expected_answers["expected_qo"]
    for (i, _) in enumerate(input["pressures"]):
        oil_form_volume_factor = correlation_results["bo"][i]
        answer = formulas.in_situ_oil_flow_rate(
            input["liquid_flow_rate"],
            oil_form_volume_factor,
            input["water_cut"]
        )
        answers.append(answer)
    assert answers == pytest.approx(expected)


def test_in_situ_gas_flow_rate(input, expected_answers, correlation_results):
    answers = []
    expected = expected_answers["expected_qg"]
    for (i, pressure) in enumerate(input["pressures"]):
        gas_solubility_in_oil = correlation_results["rso"][i]
        gas_solubility_in_water = correlation_results["rsw"][i]
        gas_form_volume_factor_bbl = correlation_results["bg_bbl"][i]
        free_gas_liquid_ratio = formulas.free_gas_liquid_ratio(
            pressure,
            input["bubble_point"],
            gas_solubility_in_oil,
            gas_solubility_in_water,
            input["water_cut"],
            input["production_gas_liquid_ratio"]
        )
        answer = formulas.in_situ_gas_flow_rate(
            input["liquid_flow_rate"],
            gas_form_volume_factor_bbl,
            free_gas_liquid_ratio
        )
        answers.append(answer)
    assert answers == pytest.approx(expected)
