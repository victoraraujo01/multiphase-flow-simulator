import math

def tests_above_bp(higher_pressure_test, lower_pressure_test):
    first_pressure, first_flow_rate = higher_pressure_test
    secnd_pressure, secnd_flow_rate = lower_pressure_test

    avg_pressure = (
        (secnd_flow_rate * first_pressure - first_flow_rate * secnd_pressure) /
        (secnd_flow_rate - first_flow_rate)
    )
    undersaturated_ip = (secnd_flow_rate - first_flow_rate) / (first_pressure - secnd_pressure)

    return (avg_pressure, undersaturated_ip)

def tests_on_opposite_sides_of_bp(b_param, bubble_point, higher_pressure_test, lower_pressure_test):
    first_pressure, first_flow_rate = higher_pressure_test
    secnd_pressure, secnd_flow_rate = lower_pressure_test

    undersaturated_ip = (
        (2 + b_param) * (secnd_flow_rate - first_flow_rate) /
        (
            (2 + b_param) * first_pressure +
            b_param * secnd_pressure -
            (1 + b_param) * (bubble_point + secnd_pressure ** 2 / bubble_point)
        )
    )
    avg_pressure = first_flow_rate / undersaturated_ip + first_pressure
    return (avg_pressure, undersaturated_ip)

def tests_below_bp(b_param, bubble_point, higher_pressure_test, lower_pressure_test):
    first_pressure, first_flow_rate = higher_pressure_test
    secnd_pressure, secnd_flow_rate = lower_pressure_test

    first_pressure_term = (
        1 + b_param * first_pressure / bubble_point -
        (1 + b_param) * (first_pressure / bubble_point) ** 2
    )
    secnd_pressure_term = (
        1 + b_param * secnd_pressure / bubble_point -
        (1 + b_param) * (secnd_pressure / bubble_point) ** 2
    )
    flow_rate_at_bp = (
        secnd_pressure_term * first_flow_rate - first_pressure_term * secnd_flow_rate
    ) / (secnd_pressure_term - first_pressure_term)

    avg_pressure = 0.
    undersaturated_ip = 0.
    if flow_rate_at_bp >= 0:
        undersaturated_ip = (
            (2 + b_param) * (secnd_flow_rate - flow_rate_at_bp) /
            (bubble_point * secnd_pressure_term)
        )
        avg_pressure = bubble_point + flow_rate_at_bp / undersaturated_ip
    else:
        first_avg_term = first_flow_rate - secnd_flow_rate
        secnd_avg_term = (
            b_param * (first_flow_rate * secnd_pressure - secnd_flow_rate * first_pressure)
        )
        thirt_avg_term = (
            (1 + b_param) *
            (first_flow_rate * secnd_pressure ** 2 - secnd_flow_rate * first_pressure ** 2)
        )
        avg_pressure = (
            (
                -secnd_avg_term -
                math.sqrt(secnd_avg_term ** 2 + 4 * first_avg_term * thirt_avg_term)
            ) / (2 * first_avg_term)
        )
        max_flow_rate = first_flow_rate / (
            1 + b_param * first_pressure / avg_pressure -
            (1 + b_param) * (first_pressure / avg_pressure) ** 2
        )
        undersaturated_ip = (2 + b_param) * max_flow_rate / avg_pressure
    return (avg_pressure, undersaturated_ip)
