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
    results["fglr"] = []
    results["rho_gas_ft"] = []
    results["rho_gas_bbl"] = []
    results["rho_oil"] = []
    results["rho_water"] = []
    results["qg"] = []
    results["qo"] = []
    results["qw"] = []
    results["vsg"] = []
    results["vso"] = []
    results["vsw"] = []
    results["lambda_l"] = []
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
        free_gas_liquid_ratio = formulas.free_gas_liquid_ratio(
            pressure,
            input["bubble_point"],
            gas_solubility_in_oil,
            gas_solubility_in_water,
            input["water_cut"],
            input["production_gas_liquid_ratio"]
        )
        gas_density_bg_in_ft = formulas.gas_density(
            input["gas_specific_gravity"],
            gas_form_volume_factor
        )
        gas_density_bg_in_bbl = formulas.gas_density(
            input["gas_specific_gravity"],
            gas_form_volume_factor_bbl,
            False
        )
        oil_density = formulas.live_oil_density(
            input["oil_specific_gravity"],
            input["gas_specific_gravity"],
            gas_solubility_in_oil,
            oil_form_volume_factor,
            input["water_cut"]
        )
        water_density = formulas.live_water_density(
            input["water_specific_gravity"],
            input["gas_specific_gravity"],
            gas_solubility_in_water,
            water_form_volume_factor,
            input["water_cut"]
        )
        oil_flow_rate = formulas.in_situ_oil_flow_rate(
            input["liquid_flow_rate"],
            oil_form_volume_factor,
            input["water_cut"]
        )
        gas_flow_rate = formulas.in_situ_gas_flow_rate(
            input["liquid_flow_rate"],
            gas_form_volume_factor_bbl,
            free_gas_liquid_ratio
        )
        water_flow_rate = formulas.in_situ_water_flow_rate(
            input["liquid_flow_rate"],
            water_form_volume_factor,
            input["water_cut"]
        )
        gas_velocity = formulas.superficial_velocity(
            gas_flow_rate,
            input["diameter"]
        )
        oil_velocity = formulas.superficial_velocity(
            oil_flow_rate,
            input["diameter"]
        )
        water_velocity = formulas.superficial_velocity(
            water_flow_rate,
            input["diameter"]
        )
        no_slip_liquid_fraction = formulas.liquid_fraction(
            oil_velocity,
            gas_velocity,
            water_velocity
        )
        results["rsw"].append(gas_solubility_in_water)
        results["rso"].append(gas_solubility_in_oil)
        results["bg"].append(gas_form_volume_factor)
        results["bg_bbl"].append(gas_form_volume_factor_bbl)
        results["co"].append(oil_compressibility)
        results["bo"].append(oil_form_volume_factor)
        results["cw"].append(water_compressibility)
        results["bw"].append(water_form_volume_factor)
        results["fglr"].append(free_gas_liquid_ratio)
        results["rho_gas_ft"].append(gas_density_bg_in_ft)
        results["rho_gas_bbl"].append(gas_density_bg_in_bbl)
        results["rho_oil"].append(oil_density)
        results["rho_water"].append(water_density)
        results["qg"].append(gas_flow_rate)
        results["qo"].append(oil_flow_rate)
        results["qw"].append(water_flow_rate)
        results["vsg"].append(gas_velocity)
        results["vso"].append(oil_velocity)
        results["vsw"].append(water_velocity)
        results["lambda_l"].append(no_slip_liquid_fraction)
    return results


@pytest.fixture(scope="module")
def expected_answers():
    answers = {}
    answers["fglr"] = [
        7.41817908307171, 7.02632049562826,
        5.33277871673351, 0.0, 0.0, 0.0, 0.0
    ]
    answers["rho_gas"] = [
        0.0406297014274935, 0.0545374088375396, 0.112030674159871,
        0.32898148509036, 3.95691607568035, 26.8524150746719,
        17.0215642317133
    ]
    answers["rho_oil"] = [
        53.4959715812803, 53.4921112656614, 53.4754149044103,
        53.5157193774093, 53.9928785276566, 54.0480899080085,
        54.0492853935634
    ]
    answers["rho_water"] = [
        64.9444055973349, 64.9454585350703, 64.9497759367957,
        64.9656191042121, 65.1535477539374, 66.6008389877806,
        66.8445362992679
    ]

    answers["qo"] = [
        0.0410932261349932, 0.0410987177779217, 0.0411224706065278,
        0.0411258623501842, 0.0407624147610332, 0.0407207750067095,
        0.0407198743269328
    ]
    answers["qg"] = [
        0.0629713499770325, 0.0444347536529818, 0.0164174587898616,
        0.0, 0.0, 0.0, 0.0
    ]
    answers["qw"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    answers["vso"] = [
        1.89302804810466, 1.89328102980193, 1.89437524349607,
        1.89453148982855, 1.87778867002612, 1.87587046524483,
        1.87582897392764
    ]
    answers["vsg"] = [
        2.90088033832969, 2.04696109036049, 0.756297640533635,
        0.0, 0.0, 0.0, 0.0
    ]
    answers["vsw"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    answers["lambda_l"] = [0.394882, 0.480499, 0.714677, 1.0, 1.0, 1.0, 1.0]
    return answers


def test_free_gas_liquid_ratio(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["fglr"]) == correlation_results["fglr"]


def test_gas_density(input, expected_answers, correlation_results):
    expected = expected_answers["rho_gas"]
    assert (
        correlation_results["rho_gas_ft"] == pytest.approx(expected) and
        correlation_results["rho_gas_bbl"] == pytest.approx(expected)
    )


def test_live_oil_density(input, expected_answers, correlation_results):
    expected = expected_answers["rho_oil"]
    assert correlation_results["rho_oil"] == pytest.approx(expected)


def test_live_water_density(input, expected_answers, correlation_results):
    expected = expected_answers["rho_water"]
    assert correlation_results["rho_water"] == pytest.approx(expected)


def test_in_situ_oil_flow_rate(input, expected_answers, correlation_results):
    expected = expected_answers["qo"]
    assert correlation_results["qo"] == pytest.approx(expected)


def test_in_situ_gas_flow_rate(input, expected_answers, correlation_results):
    expected = expected_answers["qg"]
    assert correlation_results["qg"] == pytest.approx(expected)


def test_in_situ_water_flow_rate(input, expected_answers, correlation_results):
    expected = expected_answers["qw"]
    assert correlation_results["qw"] == pytest.approx(expected)


def test_superficial_velocity(input, expected_answers, correlation_results):
    assert (
        correlation_results["vsg"] == pytest.approx(expected_answers["vsg"]) and
        correlation_results["vso"] == pytest.approx(expected_answers["vso"]) and
        correlation_results["vsw"] == pytest.approx(expected_answers["vsw"])
    )


def test_no_slip_liquid_fraction(input, expected_answers, correlation_results):
    expected = expected_answers["lambda_l"]
    assert correlation_results["lambda_l"] == pytest.approx(expected)