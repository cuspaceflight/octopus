from matplotlib import pyplot as plt
from numpy import array
from numpy import linspace

from octopus import Fluid,Orifice,STRAIGHT




def main():
    fluid = Fluid('N2O', T=250, P=18e5)
    fluid.calculate()

    orifice=Orifice(fluid,0,15e5,STRAIGHT,1e-2,1e-3)
    P_cc=linspace(0,17.9e5,100)
    mdot=[]
    for P in P_cc:
        orifice.P_cc=P
        mdot.append(orifice.m_dot_HEM())
    plt.plot(18e5-P_cc,mdot)
    plt.xlabel('Pressure Drop (Pa)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.show()


if __name__ == "__main__":
    main()
