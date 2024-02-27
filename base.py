import json
import numpy as np

from operatorTM import operatorTM
import calculationinfo


def performCalculation(calculation_info, temp_info):

    temporaryZforAzimuthalStudy = []
    temporarycontribZxforAzimuthalStudy = []
    temporarycontribZyforAzimuthalStudy = []

    for azimAngle in calculation_info.azimuthalAngles:

        temp_info.azimuthalAngles = azimAngle
        local_operator = operatorTM(temp_info)

        if float(azimAngle) == 0.0:

            if calculation_info.intransmission:

                X, Y, Zs, Zp = local_operator.calculateTsAndTpForAllAngles()

                if calculation_info.polarization == "TE":

                    Z = Zs

                if calculation_info.polarization == "TM":

                    Z = Zp

                # if np.size(calculationInfo.azimuthalAngles) == 1:

                #    sig_dataFit.emit(plotData(X,Y,Z), calculationInfo.polarization,os.path.splitext(calculationInfo.calculationOrderName)[0], True, calculationInfo.EMin, calculationInfo.EMax, calculationInfo.Zmin, calculationInfo.Zmax, False, 1.0, 'RdBu')

                # else:

                temporaryZforAzimuthalStudy.append(np.real(Z))
                temporarycontribZxforAzimuthalStudy.append(np.zeros_like(np.real(Z)))
                temporarycontribZyforAzimuthalStudy.append(np.real(Z))

            else:

                X, Y, Zs, Zp = local_operator.calculateRsAndRpForAllAngles()

                if calculation_info.polarization == "TE":

                    Z = Zs

                if calculation_info.polarization == "TM":

                    Z = Zp

                # if np.size(calculationInfo.azimuthalAngles) == 1:

                # self.sig_dataFit.emit(plotData(X,Y,Z), calculationInfo.polarization,os.path.splitext(calculationInfo.calculationOrderName)[0], True, calculationInfo.EMin, calculationInfo.EMax, calculationInfo.Zmin, calculationInfo.Zmax, False, 1.0, 'RdBu')

                # else:

                temporaryZforAzimuthalStudy.append(Z)

            # self.sig_msg.emit('{} at Phi = {} was succesfully computed'.format(tempCalculationInfo.calculationOrderName,str(azimAngle)))

        if float(azimAngle) > 0.0:

            # print('attempting OP with Phi {}'.format(azimAngle))

            if calculation_info.intransmission:

                X, Y, Zs, zss, zsp, zsx, zsy = (
                    local_operator.calculateTsAndTpForAllAngles()
                )

                # if np.size(calculationInfo.azimuthalAngles) == 1:

                # self.sig_dataFit.emit(plotData(X,Y,Zs), calculationInfo.polarization,os.path.splitext(calculationInfo.calculationOrderName)[0], True, calculationInfo.EMin, calculationInfo.EMax, calculationInfo.Zmin, calculationInfo.Zmax, False, 1.0, 'RdBu')

                # else:

                # print(zsx)
                # print(zsy)

                # plt.plot(Y,zsx)
                # plt.plot(Y,zsy)
                # plt.show()

                temporaryZforAzimuthalStudy.append(Zs)
                temporarycontribZxforAzimuthalStudy.append(zsx)
                temporarycontribZyforAzimuthalStudy.append(zsy)

            else:

                X, Y, Zs = local_operator.calculateRsAndRpForAllAngles()

                # if np.size(calculationInfo.azimuthalAngles) == 1:

                # self.sig_dataFit.emit(plotData(X,Y,Zs), calculationInfo.polarization,os.path.splitext(calculationInfo.calculationOrderName)[0], True, calculationInfo.EMin, calculationInfo.EMax, calculationInfo.Zmin, calculationInfo.Zmax, False, 1.0, 'RdBu')

                # else:

                temporaryZforAzimuthalStudy.append(Zs)

            # self.sig_msg.emit('{} at Phi = {} was succesfully computed'.format(tempCalculationInfo.calculationOrderName,str(azimAngle)))

    # if np.size(calculationInfo.azimuthalAngles) > 1:

    # if calculationInfo.intransmission is False:

    # self.sig_azimPlot.emit(plotData(Y,calculationInfo.azimuthalAngles,temporaryZforAzimuthalStudy), calculationInfo.EMin, calculationInfo.EMax, calculationInfo.Zmin, calculationInfo.Zmax, 1.0, 'RdBu')

    # else:

    # self.sig_azimPlotWithContrib.emit(plotData(Y,calculationInfo.azimuthalAngles,temporaryZforAzimuthalStudy), plotData(Y,calculationInfo.azimuthalAngles,temporarycontribZxforAzimuthalStudy), plotData(Y,calculationInfo.azimuthalAngles,temporarycontribZyforAzimuthalStudy), calculationInfo.EMin, calculationInfo.EMax, calculationInfo.Zmin, calculationInfo.Zmax, 1.0, 'RdBu')


if __name__ == "__main__":

    # Replace 'yourfile.json' with your actual file path
    with open("calculation_order.json", "r") as f:
        calculation_order = json.load(f)

    print(calculation_order)
    calculation_info = calculationinfo.CalculationInfo(
        calculation_order["intransmission"],
        calculation_order["angleMin"],
        calculation_order["angleMax"],
        calculation_order["angleStep"],
        calculation_order["wavelengthMin"],
        calculation_order["wavelengthMax"],
        calculation_order["wavelengthStep"],
        calculation_order["structure"],
        calculation_order["polarization"],
        calculation_order["azimuthalAngles"],
        calculation_order["inEnergy"],
        calculation_order["Zmin"],
        calculation_order["Zmax"],
    )
    temp_info = calculationinfo.CalculationInfo(
        calculation_order["intransmission"],
        calculation_order["angleMin"],
        calculation_order["angleMax"],
        calculation_order["angleStep"],
        calculation_order["wavelengthMin"],
        calculation_order["wavelengthMax"],
        calculation_order["wavelengthStep"],
        calculation_order["structure"],
        calculation_order["polarization"],
        calculation_order["azimuthalAngles"],
        calculation_order["inEnergy"],
        calculation_order["Zmin"],
        calculation_order["Zmax"],
    )

    performCalculation(calculation_info, temp_info)
