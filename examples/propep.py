import pypropep as ppp
from numpy import pi, sqrt


def main():
    ppp.init()
    n20 = ppp.PROPELLANTS["NITROUS OXIDE"]
    ipa = ppp.PROPELLANTS["ISOPROPYL ALCOHOL"]
    OF = 3.5

    sp = ppp.ShiftingPerformance()
    sp.add_propellants_by_mass([(ipa, 1.0), (n20, OF)])
    sp.set_state(15, 1)
    T = sp.properties[0].T
    Ve = sp.performance.Isp
    mdot = 10000 / Ve
    Ae = mdot * sp.performance.a_dotm / 101325
    At = Ae / sp.performance.ae_at
    De = sqrt(4 * Ae / pi) * 1000
    Dt = sqrt(4 * At / pi) * 1000
    print(f'exit velocity: {Ve:.0f}m/s\n'
          f'mass flow rate: {mdot:.4}kg/s\n'
          f'exit diameter: {De:.4}mm\n'
          f'throat diameter: {Dt:.4}mm\n'
          f'flame temp: {T:.0f}K')


if __name__ == "__main__":
    main()
