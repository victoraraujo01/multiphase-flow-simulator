import matplotlib.pyplot as plt

from flow.oil import Oil
from flow.water import Water
from flow.gas import Gas
from flow.mixture import Mixture
from flow.beggs_and_brill import BeggsAndBrill
from flow.tubing import Tubing
from flow.tpr import TPR
from ipr.ipr import IPR

def main():
    well_head_pressure = 150.0
    temperature = 170.0

    prod_glr = 600.0
    water_cut = 0.0

    tubing = Tubing(10000, 1.995, 90.0, 0.00015)  # 2 3/8 EU
    gas = Gas(0.7)
    oil = Oil(25.0, gas)
    water = Water(1.07, gas)
    mixture = Mixture(gas, oil, water, water_cut, prod_glr)
    correlation = BeggsAndBrill()
    bubble_point = mixture.bubble_point(temperature)

    tpr = TPR(correlation, tubing, mixture, temperature, well_head_pressure)
    x, y = tpr.simulate()

    tpr.tubing = Tubing(10000, 2.992, 90.0, 0.00015)  # 3 1/2 EU
    _, y2 = tpr.simulate()

    tpr.tubing = Tubing(10000, 3.96, 90.0, 0.00015)  # 4 1/2 EU
    _, y3 = tpr.simulate()

    ipr = IPR.from_tests(-0.2, bubble_point, (5000., 1000.), (2000., 5000.))
    y_ipr = [ipr.pressure(flow_rate) for flow_rate in x]

    plt.plot(x, y, label="TPR 2 3/8\"")
    plt.plot(x, y2, label="TPR 3 1/2\"")
    plt.plot(x, y3, label="TPR 4 1/2\"")
    plt.plot(x, y_ipr, label="IPR")
    plt.grid(True)
    plt.xlabel('Flow rate')
    plt.ylabel('Pressure')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
