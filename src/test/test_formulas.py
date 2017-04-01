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

    # Tubing data
    input_["diameter"] = 1.995  # in
    input_["inclination"] = 90  # degrees
    input_["rugosity"] = 0.000902255639097744
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
    results["fr"] = []
    results["fr_trans"] = []
    results["alpha_h_seg"] = []
    results["alpha_h_int"] = []
    results["alpha_h_dist"] = []
    results["alpha_h_tran"] = []
    results["pattern"] = []
    results["rho_liq"] = []
    results["sigma_l"] = []
    results["nlv"] = []
    results["alpha_seg"] = []
    results["alpha_int"] = []
    results["alpha_dist"] = []
    results["n_re"] = []
    results["f_n"] = []
    results["ftp_fn"] = []
    results["dp_grav"] = []
    results["dp_fric"] = []
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
        no_slip_liquid_fraction = formulas.no_slip_liquid_fraction(
            oil_velocity,
            gas_velocity,
            water_velocity
        )
        water_fraction = formulas.water_fraction(
            oil_velocity,
            water_velocity
        )
        froude = formulas.froude_number(
            oil_velocity + gas_velocity + water_velocity,
            input["diameter"]
        )
        froude_trans = formulas.transition_froude_numbers(
            no_slip_liquid_fraction
        )
        pattern = formulas.flow_pattern(froude, no_slip_liquid_fraction)
        alpha_h_seg = formulas.horz_liquid_holdup(
            formulas.FlowPattern.segregated,
            froude,
            no_slip_liquid_fraction
        )
        alpha_h_int = formulas.horz_liquid_holdup(
            formulas.FlowPattern.intermittent,
            froude,
            no_slip_liquid_fraction
        )
        alpha_h_dist = formulas.horz_liquid_holdup(
            formulas.FlowPattern.distributed,
            froude,
            no_slip_liquid_fraction
        )
        alpha_h_tran = formulas.horz_liquid_holdup(
            formulas.FlowPattern.transition,
            froude,
            no_slip_liquid_fraction
        )
        alpha_h = formulas.horz_liquid_holdup(
            pattern,
            froude,
            no_slip_liquid_fraction
        )
        rho_liq = formulas.estimate_fluid_property(
            oil_density,
            water_density,
            water_fraction
        )
        sigma_dog = correlations.dead_oil_gas_surface_tension(
            input["temperature"],
            input["oil_api_gravity"]
        )
        sigma_og = correlations.live_oil_gas_surface_tension(
            sigma_dog,
            gas_solubility_in_oil
        )
        sigma_wg = correlations.water_gas_surface_tension()
        liquid_surface_tension = formulas.estimate_fluid_property(
            sigma_og,
            sigma_wg,
            water_fraction
        )
        liquid_velocity_number = formulas.liquid_velocity_number(
            oil_velocity + water_velocity,
            rho_liq,
            liquid_surface_tension
        )
        alpha_seg = formulas.liquid_holdup_with_incl(
            alpha_h,
            formulas.FlowPattern.segregated,
            froude,
            no_slip_liquid_fraction,
            liquid_velocity_number,
            input["inclination"]
        )
        alpha_int = formulas.liquid_holdup_with_incl(
            alpha_h,
            formulas.FlowPattern.intermittent,
            froude,
            no_slip_liquid_fraction,
            liquid_velocity_number,
            input["inclination"]
        )
        alpha_dist = formulas.liquid_holdup_with_incl(
            alpha_h,
            formulas.FlowPattern.distributed,
            froude,
            no_slip_liquid_fraction,
            liquid_velocity_number,
            input["inclination"]
        )
        liquid_holdup = formulas.liquid_holdup_with_incl(
            alpha_h,
            pattern,
            froude,
            no_slip_liquid_fraction,
            liquid_velocity_number,
            input["inclination"]
        )
        oil_viscosity = correlations.live_oil_viscosity(
            pressure,
            input["bubble_point"],
            input["temperature"],
            gas_solubility_in_oil,
            input["oil_api_gravity"]
        )
        gas_viscosity = correlations.gas_viscosity(
            input["temperature"],
            input["gas_specific_gravity"],
            gas_density_bg_in_ft
        )
        water_viscosity = correlations.water_viscosity(
            pressure, input["temperature"]
        )
        liquid_viscosity = formulas.estimate_fluid_property(
            oil_viscosity,
            water_viscosity,
            water_fraction
        )
        mixture_viscosity_no_slip = formulas.estimate_fluid_property(
            liquid_viscosity,
            gas_viscosity,
            (1 - no_slip_liquid_fraction)
        )
        mixture_density_no_slip = formulas.estimate_fluid_property(
            rho_liq,
            gas_density_bg_in_ft,
            (1 - no_slip_liquid_fraction)
        )
        mixture_density_holdup = formulas.estimate_fluid_property(
            rho_liq,
            gas_density_bg_in_ft,
            (1 - liquid_holdup)
        )
        mixture_specific_grav_no_slip = formulas.density_to_specific_gravity(
            mixture_density_no_slip
        )
        mixture_specific_grav_holdup = formulas.density_to_specific_gravity(
            mixture_density_holdup
        )
        reynolds = formulas.reynolds(
            mixture_density_no_slip,
            oil_velocity + water_velocity + gas_velocity,
            input["diameter"],
            mixture_viscosity_no_slip
        )
        moody_friction = formulas.moody_friction_factor(
            reynolds,
            input["rugosity"]
        )
        friction_factor = formulas.friction_factor(
            no_slip_liquid_fraction,
            liquid_holdup,
            moody_friction
        )
        grav_pressure_gradient = formulas.gravitational_pressure_gradient(
            mixture_specific_grav_holdup,
            input["inclination"]
        )
        fric_pressure_gradient = formulas.frictional_pressure_gradient(
            friction_factor,
            mixture_specific_grav_no_slip,
            oil_velocity + water_velocity + gas_velocity,
            input["diameter"]
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
        results["fr"].append(froude)
        results["fr_trans"].append(froude_trans)
        results["alpha_h_seg"].append(alpha_h_seg)
        results["alpha_h_int"].append(alpha_h_int)
        results["alpha_h_dist"].append(alpha_h_dist)
        results["alpha_h_tran"].append(alpha_h_tran)
        results["pattern"].append(pattern)
        results["rho_liq"].append(rho_liq)
        results["sigma_l"].append(liquid_surface_tension)
        results["nlv"].append(liquid_velocity_number)
        results["alpha_seg"].append(alpha_seg)
        results["alpha_int"].append(alpha_int)
        results["alpha_dist"].append(alpha_dist)
        results["n_re"].append(reynolds)
        results["f_n"].append(moody_friction)
        results["ftp_fn"].append(friction_factor/moody_friction)
        results["dp_grav"].append(grav_pressure_gradient)
        results["dp_fric"].append(fric_pressure_gradient)
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
    answers["fr"] = [
        4.2930010412648, 2.9001960167938, 1.3124833540346,
        0.670479015410192, 0.658680735892434, 0.657335707862731,
        0.6573066297113
    ]
    answers["fr_trans"] = [
        (238.682299259474, 0.00916886746446259, 0.385272526093591, 261.80215203344),
        (253.255052983182, 0.00564867595258456, 0.289771480657344, 69.7795452024043),
        (285.514562859197, 0.00212007297799988, 0.162844879304547, 4.80814628152552),
        (316, 0.0009252, 0.1, 0.5),
        (316, 0.0009252, 0.1, 0.5),
        (316, 0.0009252, 0.1, 0.5),
        (316, 0.0009252, 0.1, 0.5)
    ]
    answers["alpha_h_seg"] = [
        0.550491330724074, 0.626376432808134, 0.813349452259821,
        1.01460229702376, 1.0161670044178, 1.01634731595193, 1.01635121853251
    ]
    answers["alpha_h_int"] = [
        0.501163698847922, 0.560441385674535,
        0.714677, 1.0, 1.0, 1.0, 1.0
        # 0.702664285046605, 0.850864188258459, 0.851125558383153,
        # 0.851155657072191, 0.851156308467212
    ]
    answers["alpha_h_dist"] = [
        0.567279421596702, 0.651336332312092, 0.861373062062183,
        1.09124621336544, 1.09242669094078, 1.09256269027836, 1.09256563370637
    ]
    answers["alpha_h_tran"] = [
        0.394882, 0.480499, 0.714677, 1.0, 1.0, 1.0, 1.0
        # -0.0113519048845774, -0.0453475970426878, -0.0890487660906993,
        # -0.0919502819555191, -0.0795397220250882, -0.0781140475828659,
        # -0.0780832012158941
    ]
    answers["pattern"] = [2, 2, 2, 1, 1, 1, 1]
    answers["rho_liq"] = [
        53.4959715812803, 53.4921112656614, 53.4754149044103,
        53.5157193774093, 53.9928785276566, 54.0480899080085,
        54.0492853935634
    ]
    answers["sigma_l"] = [
        27.2985627448439, 27.2504251653417, 27.0433185286994,
        26.4050356952703, 26.4050356952703, 26.4050356952703,
        26.4050356952703
    ]
    answers["nlv"] = [
        4.34066376295262, 4.3430814109979, 4.35354789930155,
        4.38080833941479, 4.35173977682921, 4.34840529840787,
        4.34833316327154
    ]
    answers["alpha_seg"] = [
        0.667661453417903, 0.711192036405742, 0.80758542712356, 1., 1., 1., 1.
    ]
    answers["alpha_int"] = [
        0.527262065160054, 0.587347992764985, 0.736112769123655, 1., 1., 1., 1.
    ]
    answers["alpha_dist"] = [
        0.501163698847922, 0.728120216442387, 0.928501833186807, 1., 1., 1., 1.
    ]
    answers["n_re"] = [
        11303.2719423257, 9331.16182613241, 6377.46183950293, 4767.93954013749,
        3886.99890720668, 303.538772330867, 239.117244836072
    ]
    answers["f_n"] = [
        0.0315981341515841, 0.0330847388475963, 0.0364849713769,
        0.0395497578694166, 0.0419685286859213, 0.210846210876276,
        0.267651126726036
    ]
    answers["ftp_fn"] = [
        1.44285834982483, 1.44128142170068, 1.43818299356602, 1, 1, 1, 1
    ]
    answers["dp_grav"] = [
        -0.195860358694611, -0.218172248441688, -0.273355531374995,
        -0.371351065551574, -0.37466212183454, -0.375045239265508,
        -0.375053534862387
    ]
    answers["dp_fric"] = [
        -0.0143839875404616, -0.0123652766188966, -0.00915353193836248,
        -0.00493117028077515, -0.00518650499949317, -0.0260299293127166,
        -0.033042025191692
    ]
    return answers


def test_free_gas_liquid_ratio(input, expected_answers, correlation_results):
    assert pytest.approx(expected_answers["fglr"]) == correlation_results["fglr"]


def test_gas_density(input, expected_answers, correlation_results):
    expected = expected_answers["rho_gas"]
    assert correlation_results["rho_gas_ft"] == pytest.approx(expected)
    assert correlation_results["rho_gas_bbl"] == pytest.approx(expected, 1e-4)


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
    assert correlation_results["vsg"] == pytest.approx(expected_answers["vsg"])
    assert correlation_results["vso"] == pytest.approx(expected_answers["vso"])
    assert correlation_results["vsw"] == pytest.approx(expected_answers["vsw"])


def test_no_slip_liquid_fraction(input, expected_answers, correlation_results):
    expected = expected_answers["lambda_l"]
    assert correlation_results["lambda_l"] == pytest.approx(expected)


def test_froude_number(input, expected_answers, correlation_results):
    expected = expected_answers["fr"]
    assert correlation_results["fr"] == pytest.approx(expected)


def test_trans_froude_number(input, expected_answers, correlation_results):
    expected = [tuple(map(pytest.approx, tup)) for tup in expected_answers["fr_trans"]]
    assert correlation_results["fr_trans"] == expected


def test_horz_liquid_holdup(input, expected_answers, correlation_results):
    assert correlation_results["alpha_h_seg"] == pytest.approx(expected_answers["alpha_h_seg"])
    assert correlation_results["alpha_h_int"] == pytest.approx(expected_answers["alpha_h_int"])
    assert correlation_results["alpha_h_dist"] == pytest.approx(expected_answers["alpha_h_dist"])
    assert correlation_results["alpha_h_tran"] == pytest.approx(expected_answers["alpha_h_tran"])


def test_pattern(input, expected_answers, correlation_results):
    results = [pat.value for pat in correlation_results["pattern"]]
    assert results == expected_answers["pattern"]


def test_liquid_density(input, expected_answers, correlation_results):
    assert correlation_results["rho_liq"] == pytest.approx(expected_answers["rho_liq"])


def test_liquid_surface_tension(input, expected_answers, correlation_results):
    assert correlation_results["sigma_l"] == pytest.approx(expected_answers["sigma_l"], 1e-2)


def test_liquid_velocity_number(input, expected_answers, correlation_results):
    assert correlation_results["nlv"] == pytest.approx(expected_answers["nlv"], 1e-2)


def test_liquid_holdup_with_incl(input, expected_answers, correlation_results):
    assert correlation_results["alpha_seg"] == pytest.approx(expected_answers["alpha_seg"], 1e-3)
    assert correlation_results["alpha_int"] == pytest.approx(expected_answers["alpha_int"], 1e-3)
    assert correlation_results["alpha_dist"] == pytest.approx(expected_answers["alpha_dist"], 1e-3)


def test_reynolds(input, expected_answers, correlation_results):
    assert correlation_results["n_re"] == pytest.approx(expected_answers["n_re"])


def test_moody_friction(input, expected_answers, correlation_results):
    assert correlation_results["f_n"] == pytest.approx(expected_answers["f_n"])


def test_friction_factor(input, expected_answers, correlation_results):
    assert correlation_results["ftp_fn"] == pytest.approx(expected_answers["ftp_fn"], 1e-5)


def test_grav_pressure_gradient(input, expected_answers, correlation_results):
    assert correlation_results["dp_grav"] == pytest.approx(expected_answers["dp_grav"], 1e-4)


def test_fric_pressure_gradient(input, expected_answers, correlation_results):
    assert correlation_results["dp_fric"] == pytest.approx(expected_answers["dp_fric"], 1e-5)
