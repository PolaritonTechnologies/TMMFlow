# Author: Florian Le Roux
# -*- coding: utf-8 -*-

# libraries
import numpy as np  # numpy
import pickle
import math
import core
import mat
import utils
import csvparser

# from plotData import plotData
# from plotPopUp import PlotPopUp

import os
from os import walk
from calculationinfo import CalculationInfo

# matplotlib inline plots
import matplotlib.pylab as plt
import matplotlib.pyplot as pyplt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm


class operatorTM:

    def __init__(self, CalculationInfo):

        self.localParser = csvparser.csvparser()
        self.calculationInfo = CalculationInfo

    def calculateRsAndRpForAllAngles(self):

        ##check if the simulation has already been performed

        for root, dirs, files in walk(".", topdown=False):

            if self.calculationInfo.calculationOrderName in files:

                if self.calculationInfo.inEnergy:

                    with open(
                        str(
                            "ComputedSimulations\{}".format(
                                os.path.splitext(
                                    self.calculationInfo.calculationOrderName
                                )[0]
                            )
                            + "E"
                            + ".txt"
                        ),
                        "rb",
                    ) as f:
                        packed = pickle.load(f)

                        X = packed[0]
                        Ye = packed[1]
                        Zs = packed[2]
                        Zp = packed[3]
                        return X, Ye, Zs, Zp

                else:

                    with open(
                        str(
                            "ComputedSimulations\{}".format(
                                self.calculationInfo.calculationOrderName
                            )
                        ),
                        "rb",
                    ) as f:
                        packed = pickle.load(f)

                        X = packed[0]
                        Y = packed[1]
                        Zs = packed[2]
                        Zp = packed[3]
                        return X, Y, Zs, Zp

        e_listx_wvl = []
        e_listy_wvl = []
        e_listz_wvl = []

        # polar angle in radians - correction for Normal incidence to avoid
        # crash

        theta_min = float(self.calculationInfo.angleMin) * np.pi / 1.8e2
        theta_max = float(self.calculationInfo.angleMax) * np.pi / 1.8e2
        theta_stepInRadians = float(self.calculationInfo.angleStep) * np.pi / 1.8e2
        theta_steps = (theta_max - theta_min) // (theta_stepInRadians) + 2
        v_theta = np.linspace(theta_min, theta_max, int(theta_steps))

        if v_theta[0] == 0:

            v_theta[0] = 0.1 * np.pi / 1.8e2

        # wavelength in nm
        wvl_min = float(self.calculationInfo.wavelengthMin)
        wvl_max = float(self.calculationInfo.wavelengthMax)
        wvl_steps = (wvl_max - wvl_min) / (
            float(self.calculationInfo.wavelengthStep)
        ) + 1
        # print('nb steps:{}'.format(wvl_steps))
        wavelength = np.linspace(wvl_min, wvl_max, int(wvl_steps))

        # azimuthal angle radians
        self.phi_0 = float(self.calculationInfo.azimuthalAngles) * np.pi / 1.8e2

        currentStructure = list(reversed(self.calculationInfo.structure))

        # import the necessary indices, assign layer to corresponding index

        materialAtLayer = np.zeros(len(currentStructure))
        uniqueMaterials = []

        for ilayer, layer in enumerate(currentStructure):

            if layer[0] not in uniqueMaterials:

                uniqueMaterials.append(layer[0])

            materialAtLayer[ilayer] = self.findIndexInUniqueMaterials(
                layer[0], uniqueMaterials
            )

        # the necessary indices are imported inside the array indicesForCalculation
        # which contains as first dimension the name and second dimension
        # the result from csvparser with asked wavelength.

        indicesForCalculation = []

        for iUnique, uniqueMaterialName in enumerate(uniqueMaterials):

            splnx, splkx, splny, splky, splnz, splkz = self.localParser.importFromFile(
                uniqueMaterialName
            )

            importedMaterial = self.localParser.indexAccordingToPromptedWavelengths(
                wavelength, splnx, splkx, splny, splky, splnz, splkz
            )

            indicesForCalculation.append(importedMaterial)

        d_incidentMedium = 0.0  # incident medium has zero thickness
        d_substrate = 0.0  # substrate has zero thickness

        # multilayer thicknesses: incident medium and substrate have zero thickness
        # all thickness in nm

        d_list = np.zeros(len(currentStructure) + 2)
        d_list[0] = d_incidentMedium
        d_list[-1] = d_substrate

        indexXarray = np.zeros(len(currentStructure) + 2, dtype=complex)
        indexXarray[0] = 1.0  # air
        indexXarray[-1] = 1.5  # ~quartz/glass

        indexYarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexYarray[0] = 1.0  # air
        indexYarray[-1] = 1.5  # ~quartz/glass

        indexZarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexZarray[0] = 1.0  # air
        indexZarray[-1] = 1.5  # ~quartz/glass

        for ilayer, layer in enumerate(currentStructure):

            d_list[ilayer + 1] = currentStructure[ilayer][1]

        all_tensors = []

        for iwl, wvl in enumerate(wavelength):

            # multilayer composition

            for ilayer, layer in enumerate(currentStructure):

                currentIndexInIndices = int(materialAtLayer[ilayer])
                indexXarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][0][iwl]
                    + indicesForCalculation[currentIndexInIndices][1][iwl] * 1j
                )
                indexYarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][2][iwl]
                    + indicesForCalculation[currentIndexInIndices][3][iwl] * 1j
                )
                indexZarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][4][iwl]
                    + indicesForCalculation[currentIndexInIndices][5][iwl] * 1j
                )

            e_listx = indexXarray**2
            e_listy = indexYarray**2
            e_listz = indexZarray**2

            # optical constant tensor
            m_eps = np.zeros((len(e_listx), 3, 3), dtype=np.complex128)

            # filling dielectric tensor

            m_eps[:, 0, 0] = e_listx
            m_eps[:, 1, 1] = e_listy
            m_eps[:, 2, 2] = e_listz

            all_tensors.append(m_eps)

        total_R_s, total_R_p, total_r_s, total_r_p = self.calculate(
            wavelength, v_theta, 0, all_tensors, d_list
        )

        # all calculations linked to the azimuthal angle happen here
        # so far only the s polarization is handled
        # The calculation relies on exchanging y and x to calculate all the necessary ry and rx terms

        if self.phi_0 > 0:

            # print(total_r_s.shape)

            total_r_s_y = (
                total_r_s  # reflection coefficients calculated naturally above
            )

            # switched structure to calculate rx coefficients

            total_r_s_x, total_R_s_xTest = self.switchedCalculation()

            total_r_s_s = np.zeros_like(total_r_s_y)

            total_r_s_p = np.zeros_like(total_r_s_y)

            total_R_s = np.zeros_like(total_r_s_y)

            contribX = np.zeros_like(total_r_s_y)

            contribY = np.zeros_like(total_r_s_y)

            # print(total_r_s_x.shape)
            # print(total_r_s_s.shape)

            # theta form the columns, wvl form the lines

            cosphiSquared = np.cos(self.phi_0) ** 2
            sinphiSquared = np.sin(self.phi_0) ** 2
            sinphiS2D2 = np.sin(2 * self.phi_0) / 2

            for tIdx, angleTheta in enumerate(v_theta):

                for wIdx, wvlVal in enumerate(wavelength):

                    (
                        total_r_s_s[wIdx, tIdx],
                        total_r_s_p[wIdx, tIdx],
                        total_R_s[wIdx, tIdx],
                        contribX[wIdx, tIdx],
                        contribY[wIdx, tIdx],
                    ) = core.azimuthalCalculation(
                        total_r_s_y[wIdx, tIdx],
                        total_r_s_x[wIdx, tIdx],
                        cosphiSquared,
                        sinphiSquared,
                        sinphiS2D2,
                    )

        # convert angles back to degrees for display
        v_theta = v_theta * 1.8e2 / np.pi

        X = v_theta

        Y = wavelength

        # print('Operator: using Phi {}'.format(self.phi_0))

        if self.phi_0 == 0:

            Zs = total_R_s
            Zp = total_R_p

            packedData = X, Y, Zs, Zp

            with open(
                "ComputedSimulations\{}".format(
                    self.calculationInfo.calculationOrderName
                ),
                "wb",
            ) as f:
                pickle.dump(packedData, f)

            # save another set with energies instead of wavelength

            Ye = np.zeros_like(Y)

            for iwvl, wvl in enumerate(wavelength):

                Ye[iwvl] = 1239.82 / wvl

            packedDataInEnergy = X, Ye, Zs, Zp

            # save another set with energies instead of wavelength

            with open(
                str(
                    "ComputedSimulations\{}".format(
                        os.path.splitext(self.calculationInfo.calculationOrderName)[0]
                    )
                    + "E"
                    + ".txt"
                ),
                "wb",
            ) as f:
                pickle.dump(packedDataInEnergy, f)

            if self.calculationInfo.inEnergy:

                return X, Ye, Zs, Zp

            else:

                return X, Y, Zs, Zp

        else:

            Zs = total_R_s

            # only the s factors are available

            packedData = X, Y, np.real(Zs)

            with open(
                "ComputedSimulations\{}At{}".format(
                    self.calculationInfo.calculationOrderName, str(self.phi_0)
                ),
                "wb",
            ) as f:
                pickle.dump(packedData, f)

            # save another set with energies instead of wavelength

            Ye = np.zeros_like(Y)

            for iwvl, wvl in enumerate(wavelength):

                Ye[iwvl] = 1239.82 / wvl

            packedDataInEnergy = X, Ye, np.real(Zs)

            # plt.plot(Ye, np.real(Zs))
            # plt.show()

            # save another set with energies instead of wavelength

            with open(
                str(
                    "ComputedSimulations\{}".format(
                        os.path.splitext(self.calculationInfo.calculationOrderName)[0]
                    )
                    + "E"
                    + "{}".format(str(self.phi_0))
                    + ".txt"
                ),
                "wb",
            ) as f:
                pickle.dump(packedDataInEnergy, f)

            if self.calculationInfo.inEnergy:

                # return X, Ye, np.real(total_R_s_xTest), zss, zsp

                return X, Ye, np.real(Zs)

            else:

                return X, Y, np.real(Zs)

    def calculate(self, wavelength, v_theta, phi_0, all_tensors, d_list):

        total_R_s = np.zeros((len(wavelength), len(v_theta)))
        total_R_p = np.zeros_like(total_R_s)

        total_r_s = np.zeros_like(total_R_s, dtype=np.complex128)
        total_r_p = np.zeros_like(total_R_p, dtype=np.complex128)

        for i_w, wvl in enumerate(wavelength):

            # initializing reflectance output vector
            v_r_p = np.zeros_like(v_theta)
            v_r_s = np.zeros_like(v_theta)

            # initializing reflection coefficients output vector
            a_r_p = np.zeros_like(v_theta, dtype=np.complex128)
            a_r_s = np.zeros_like(v_theta, dtype=np.complex128)

            m_eps = all_tensors[i_w]

            ##angle loop
            for i_t, t in enumerate(v_theta):

                ##------Computing------
                m_r_ps = core.rt(wvl, t, phi_0, m_eps, d_list)[
                    "m_r_ps"
                ]  # reflection matrix

                # extracting values from the matrix
                r_pp = m_r_ps[0, 0]
                r_ps = m_r_ps[0, 1]
                r_sp = m_r_ps[1, 0]
                r_ss = m_r_ps[1, 1]

                r_s = r_sp + r_ss
                r_p = r_pp + r_ps

                a_r_p[i_t] = r_p  # reflection coefficient for p (TM)
                a_r_s[i_t] = r_s  # reflection coefficient for s (TE)

                v_r_p[i_t] = utils.R_ps_rl(m_r_ps)[
                    "R_p"
                ]  # getting p-polarized reflectance (TM)
                v_r_s[i_t] = utils.R_ps_rl(m_r_ps)[
                    "R_s"
                ]  # getting s-polarized reflectance (TE)

            total_r_s[i_w] = a_r_s
            total_r_p[i_w] = a_r_p
            total_R_s[i_w] = v_r_s
            total_R_p[i_w] = v_r_p

        return total_R_s, total_R_p, total_r_s, total_r_p

    def switchedCalculation(self):

        e_listx_wvl = []
        e_listy_wvl = []
        e_listz_wvl = []

        # polar angle in radians - correction for Normal incidence to avoid
        # crash

        theta_min = float(self.calculationInfo.angleMin) * np.pi / 1.8e2
        theta_max = float(self.calculationInfo.angleMax) * np.pi / 1.8e2
        theta_stepInRadians = float(self.calculationInfo.angleStep) * np.pi / 1.8e2
        theta_steps = (theta_max - theta_min) // (theta_stepInRadians) + 2
        v_theta = np.linspace(theta_min, theta_max, int(theta_steps))

        if v_theta[0] == 0:

            v_theta[0] = 0.16 * np.pi / 1.8e2

        # wavelength in nm
        wvl_min = float(self.calculationInfo.wavelengthMin)
        wvl_max = float(self.calculationInfo.wavelengthMax)
        wvl_steps = (wvl_max - wvl_min) / (
            float(self.calculationInfo.wavelengthStep)
        ) + 1
        wavelength = np.linspace(wvl_min, wvl_max, int(wvl_steps))

        currentStructure = list(reversed(self.calculationInfo.structure))

        # import the necessary indices, assign layer to corresponding index

        materialAtLayer = np.zeros(len(currentStructure))
        uniqueMaterials = []

        for ilayer, layer in enumerate(currentStructure):

            if layer[0] not in uniqueMaterials:

                uniqueMaterials.append(layer[0])

            materialAtLayer[ilayer] = self.findIndexInUniqueMaterials(
                layer[0], uniqueMaterials
            )

        # the necessary indices are imported inside the array indicesForCalculation
        # which contains as first dimension the name and second dimension
        # the result from csvparser with asked wavelength.

        indicesForCalculation = []

        for iUnique, uniqueMaterialName in enumerate(uniqueMaterials):

            # the switch is made here, y become x

            splny, splky, splnx, splkx, splnz, splkz = self.localParser.importFromFile(
                uniqueMaterialName
            )

            importedMaterial = self.localParser.indexAccordingToPromptedWavelengths(
                wavelength, splnx, splkx, splny, splky, splnz, splkz
            )

            indicesForCalculation.append(importedMaterial)

        d_incidentMedium = 0.0  # incident medium has zero thickness
        d_substrate = 0.0  # substrate has zero thickness

        # multilayer thicknesses: incident medium and substrate have zero thickness
        # all thickness in nm

        d_list = np.zeros(len(currentStructure) + 2)
        d_list[0] = d_incidentMedium
        d_list[-1] = d_substrate

        indexXarray = np.zeros(len(currentStructure) + 2, dtype=complex)
        indexXarray[0] = 1.0  # air
        indexXarray[-1] = 1.5  # ~quartz/glass

        indexYarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexYarray[0] = 1.0  # air
        indexYarray[-1] = 1.5  # ~quartz/glass

        indexZarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexZarray[0] = 1.0  # air
        indexZarray[-1] = 1.5  # ~quartz/glass

        for ilayer, layer in enumerate(currentStructure):

            d_list[ilayer + 1] = currentStructure[ilayer][1]

        all_tensors = []

        for iwl, wvl in enumerate(wavelength):

            # multilayer composition

            for ilayer, layer in enumerate(currentStructure):

                currentIndexInIndices = int(materialAtLayer[ilayer])
                indexXarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][0][iwl]
                    + indicesForCalculation[currentIndexInIndices][1][iwl] * 1j
                )
                indexYarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][2][iwl]
                    + indicesForCalculation[currentIndexInIndices][3][iwl] * 1j
                )
                indexZarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][4][iwl]
                    + indicesForCalculation[currentIndexInIndices][5][iwl] * 1j
                )

            e_listx = indexXarray**2
            e_listy = indexYarray**2
            e_listz = indexZarray**2

            # optical constant tensor
            m_eps = np.zeros((len(e_listx), 3, 3), dtype=np.complex128)

            # filling dielectric tensor

            m_eps[:, 0, 0] = e_listx
            m_eps[:, 1, 1] = e_listy
            m_eps[:, 2, 2] = e_listz

            all_tensors.append(m_eps)

        total_R_s, total_R_p, total_r_s, total_r_p = self.calculate(
            wavelength, v_theta, 0, all_tensors, d_list
        )

        return total_r_s, total_R_s

    def calculateTsAndTpForAllAngles(self):

        e_listx_wvl = []
        e_listy_wvl = []
        e_listz_wvl = []

        # polar angle in radians - correction for Normal incidence to avoid
        # crash

        theta_min = float(self.calculationInfo.angleMin) * np.pi / 1.8e2
        theta_max = float(self.calculationInfo.angleMax) * np.pi / 1.8e2
        theta_stepInRadians = float(self.calculationInfo.angleStep) * np.pi / 1.8e2
        theta_steps = (theta_max - theta_min) // (theta_stepInRadians) + 2
        v_theta = np.linspace(theta_min, theta_max, int(theta_steps))

        if v_theta[0] == 0:

            v_theta[0] = 0.1 * np.pi / 1.8e2

        # wavelength in nm
        wvl_min = float(self.calculationInfo.wavelengthMin)
        wvl_max = float(self.calculationInfo.wavelengthMax)
        wvl_steps = (wvl_max - wvl_min) / (
            float(self.calculationInfo.wavelengthStep)
        ) + 1
        # print('nb steps:{}'.format(wvl_steps))
        wavelength = np.linspace(wvl_min, wvl_max, int(wvl_steps))

        # azimuthal angle radians
        self.phi_0 = float(self.calculationInfo.azimuthalAngles) * np.pi / 1.8e2

        currentStructure = list(reversed(self.calculationInfo.structure))

        # import the necessary indices, assign layer to corresponding index

        materialAtLayer = np.zeros(len(currentStructure))
        uniqueMaterials = []

        for ilayer, layer in enumerate(currentStructure):

            if layer[0] not in uniqueMaterials:

                uniqueMaterials.append(layer[0])

            materialAtLayer[ilayer] = self.findIndexInUniqueMaterials(
                layer[0], uniqueMaterials
            )

        # the necessary indices are imported inside the array indicesForCalculation
        # which contains as first dimension the name and second dimension
        # the result from csvparser with asked wavelength.

        indicesForCalculation = []

        for iUnique, uniqueMaterialName in enumerate(uniqueMaterials):

            splnx, splkx, splny, splky, splnz, splkz = self.localParser.importFromFile(
                uniqueMaterialName
            )

            importedMaterial = self.localParser.indexAccordingToPromptedWavelengths(
                wavelength, splnx, splkx, splny, splky, splnz, splkz
            )

            indicesForCalculation.append(importedMaterial)

        d_incidentMedium = 0.0  # incident medium has zero thickness
        d_substrate = 0.0  # substrate has zero thickness

        # multilayer thicknesses: incident medium and substrate have zero thickness
        # all thickness in nm

        d_list = np.zeros(len(currentStructure) + 2)
        d_list[0] = d_incidentMedium
        d_list[-1] = d_substrate

        indexXarray = np.zeros(len(currentStructure) + 2, dtype=complex)
        indexXarray[0] = 1.0  # air
        indexXarray[-1] = 1.5  # ~quartz/glass

        indexYarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexYarray[0] = 1.0  # air
        indexYarray[-1] = 1.5  # ~quartz/glass

        indexZarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexZarray[0] = 1.0  # air
        indexZarray[-1] = 1.5  # ~quartz/glass

        for ilayer, layer in enumerate(currentStructure):

            d_list[ilayer + 1] = currentStructure[ilayer][1]

        all_tensors = []

        for iwl, wvl in enumerate(wavelength):

            # multilayer composition

            for ilayer, layer in enumerate(currentStructure):

                currentIndexInIndices = int(materialAtLayer[ilayer])
                indexXarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][0][iwl]
                    + indicesForCalculation[currentIndexInIndices][1][iwl] * 1j
                )
                indexYarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][2][iwl]
                    + indicesForCalculation[currentIndexInIndices][3][iwl] * 1j
                )
                indexZarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][4][iwl]
                    + indicesForCalculation[currentIndexInIndices][5][iwl] * 1j
                )

            e_listx = indexXarray**2
            e_listy = indexYarray**2
            e_listz = indexZarray**2

            # optical constant tensor
            m_eps = np.zeros((len(e_listx), 3, 3), dtype=np.complex128)

            # filling dielectric tensor

            m_eps[:, 0, 0] = e_listx
            m_eps[:, 1, 1] = e_listy
            m_eps[:, 2, 2] = e_listz

            all_tensors.append(m_eps)

        total_T_s, total_T_p, total_t_s, total_t_p = self.calculateT(
            wavelength, v_theta, 0, all_tensors, d_list
        )

        # all calculations linked to the azimuthal angle happen here
        # so far only the s polarization is handled
        # The calculation relies on exchanging y and x to calculate all the necessary ry and rx terms

        if self.phi_0 > 0:

            # print(total_t_s.shape)

            total_t_s_y = (
                total_t_s  # reflection coefficients calculated naturally above
            )

            # switched structure to calculate rx coefficients

            total_t_s_x, total_T_s_x = self.switchedCalculationT()

            total_t_s_s = np.zeros_like(total_t_s_y)

            total_t_s_p = np.zeros_like(total_t_s_y)

            total_T_s = np.zeros_like(total_t_s_y)

            contribX = np.zeros_like(total_t_s_y)

            contribY = np.zeros_like(total_t_s_y)

            # theta form the columns, wvl form the lines

            cosphiSquared = np.cos(self.phi_0) ** 2
            sinphiSquared = np.sin(self.phi_0) ** 2
            sinphiS2D2 = np.sin(2 * self.phi_0) / 2

            for tIdx, angleTheta in enumerate(v_theta):

                for wIdx, wvlVal in enumerate(wavelength):

                    (
                        total_t_s_s[wIdx, tIdx],
                        total_t_s_p[wIdx, tIdx],
                        total_T_s[wIdx, tIdx],
                        contribX[wIdx, tIdx],
                        contribY[wIdx, tIdx],
                    ) = core.azimuthalCalculation(
                        total_t_s_y[wIdx, tIdx],
                        total_t_s_x[wIdx, tIdx],
                        cosphiSquared,
                        sinphiSquared,
                        sinphiS2D2,
                    )

        # convert angles back to degrees for display
        v_theta = v_theta * 1.8e2 / np.pi

        X = v_theta

        Y = wavelength

        if self.phi_0 == 0:

            Zs = total_T_s
            Zp = total_T_p

            packedData = X, Y, Zs, Zp

            with open(
                "ComputedSimulations\T{}".format(
                    self.calculationInfo.calculationOrderName
                ),
                "wb",
            ) as f:
                pickle.dump(packedData, f)

            # save another set with energies instead of wavelength

            Ye = np.zeros_like(Y)

            for iwvl, wvl in enumerate(wavelength):

                Ye[iwvl] = 1239.82 / wvl

            packedDataInEnergy = X, Ye, Zs, Zp

            # save another set with energies instead of wavelength

            with open(
                str(
                    "ComputedSimulations\T{}".format(
                        os.path.splitext(self.calculationInfo.calculationOrderName)[0]
                    )
                    + "E"
                    + ".txt"
                ),
                "wb",
            ) as f:
                pickle.dump(packedDataInEnergy, f)

            if self.calculationInfo.inEnergy:

                return X, Ye, Zs, Zp

            else:

                return X, Y, Zs, Zp

        else:

            Zs = total_T_s
            zss = total_t_s_s
            zsp = total_t_s_p
            zsx = contribX
            zsy = contribY

            # only the s factors are available

            packedData = X, Y, np.real(Zs), zss, zsp

            with open(
                "ComputedSimulations\T{}At{}".format(
                    self.calculationInfo.calculationOrderName, str(self.phi_0)
                ),
                "wb",
            ) as f:
                pickle.dump(packedData, f)

            # save another set with energies instead of wavelength

            Ye = np.zeros_like(Y)

            for iwvl, wvl in enumerate(wavelength):

                Ye[iwvl] = 1239.82 / wvl

            # plt.plot(Ye,zsy)
            # plt.plot(Ye,zsx)
            # plt.show()

            # print('sending else')
            packedDataInEnergy = (
                X,
                Ye,
                np.real(Zs),
                zss,
                zsp,
                np.real(zsx),
                np.real(zsy),
            )

            # save another set with energies instead of wavelength

            with open(
                str(
                    "ComputedSimulations\T{}".format(
                        os.path.splitext(self.calculationInfo.calculationOrderName)[0]
                    )
                    + "E"
                    + "{}".format(str(self.phi_0))
                    + ".txt"
                ),
                "wb",
            ) as f:
                pickle.dump(packedDataInEnergy, f)

            if self.calculationInfo.inEnergy:

                return X, Ye, np.real(Zs), zss, zsp, np.real(zsx), np.real(zsy)

            else:

                return X, Y, np.real(Zs), zss, zsp, zsx, zsy

    def calculateT(self, wavelength, v_theta, phi_0, all_tensors, d_list):

        total_T_s = np.zeros((len(wavelength), len(v_theta)))
        total_T_p = np.zeros_like(total_T_s)

        total_t_s = np.zeros_like(total_T_s, dtype=np.complex128)
        total_t_p = np.zeros_like(total_T_p, dtype=np.complex128)

        for i_w, wvl in enumerate(wavelength):

            # initializing reflectance output vector
            v_t_p = np.zeros_like(v_theta)
            v_t_s = np.zeros_like(v_theta)

            # initializing reflection coefficients output vector
            a_t_p = np.zeros_like(v_theta, dtype=np.complex128)
            a_t_s = np.zeros_like(v_theta, dtype=np.complex128)

            m_eps = all_tensors[i_w]

            ##angle loop
            for i_t, t in enumerate(v_theta):

                ##------Computing------
                m_t_ps = core.rt(wvl, t, phi_0, m_eps, d_list)[
                    "m_t_ps"
                ]  # transmission matrix

                # extracting values from the matrix
                t_pp = m_t_ps[0, 0]
                t_ps = m_t_ps[0, 1]
                t_sp = m_t_ps[1, 0]
                t_ss = m_t_ps[1, 1]

                t_s = t_sp + t_ss
                t_p = t_pp + t_ps

                a_t_p[i_t] = t_p  # transmission coefficient for p (TM)
                a_t_s[i_t] = t_s  # transmission coefficient for s (TE)

                v_t_p[i_t] = utils.T_ps_rl(m_t_ps, t, 1, 1.5)[
                    "T_p"
                ]  # getting p-polarized transmittance (TM)
                v_t_s[i_t] = utils.T_ps_rl(m_t_ps, t, 1, 1.5)[
                    "T_s"
                ]  # getting s-polarized transmittance (TE)

            total_t_s[i_w] = a_t_s
            total_t_p[i_w] = a_t_p
            total_T_s[i_w] = v_t_s
            total_T_p[i_w] = v_t_p

        return total_T_s, total_T_p, total_t_s, total_t_p

    def switchedCalculationT(self):

        e_listx_wvl = []
        e_listy_wvl = []
        e_listz_wvl = []

        # polar angle in radians - correction for Normal incidence to avoid
        # crash

        theta_min = float(self.calculationInfo.angleMin) * np.pi / 1.8e2
        theta_max = float(self.calculationInfo.angleMax) * np.pi / 1.8e2
        theta_stepInRadians = float(self.calculationInfo.angleStep) * np.pi / 1.8e2
        theta_steps = (theta_max - theta_min) // (theta_stepInRadians) + 2
        v_theta = np.linspace(theta_min, theta_max, int(theta_steps))

        if v_theta[0] == 0:

            v_theta[0] = 0.16 * np.pi / 1.8e2

        # wavelength in nm
        wvl_min = float(self.calculationInfo.wavelengthMin)
        wvl_max = float(self.calculationInfo.wavelengthMax)
        wvl_steps = (wvl_max - wvl_min) / (
            float(self.calculationInfo.wavelengthStep)
        ) + 1
        wavelength = np.linspace(wvl_min, wvl_max, int(wvl_steps))

        currentStructure = list(reversed(self.calculationInfo.structure))

        # import the necessary indices, assign layer to corresponding index

        materialAtLayer = np.zeros(len(currentStructure))
        uniqueMaterials = []

        for ilayer, layer in enumerate(currentStructure):

            if layer[0] not in uniqueMaterials:

                uniqueMaterials.append(layer[0])

            materialAtLayer[ilayer] = self.findIndexInUniqueMaterials(
                layer[0], uniqueMaterials
            )

        # the necessary indices are imported inside the array indicesForCalculation
        # which contains as first dimension the name and second dimension
        # the result from csvparser with asked wavelength.

        indicesForCalculation = []

        for iUnique, uniqueMaterialName in enumerate(uniqueMaterials):

            # the switch is made here, y become x

            splny, splky, splnx, splkx, splnz, splkz = self.localParser.importFromFile(
                uniqueMaterialName
            )

            importedMaterial = self.localParser.indexAccordingToPromptedWavelengths(
                wavelength, splnx, splkx, splny, splky, splnz, splkz
            )

            indicesForCalculation.append(importedMaterial)

        d_incidentMedium = 0.0  # incident medium has zero thickness
        d_substrate = 0.0  # substrate has zero thickness

        # multilayer thicknesses: incident medium and substrate have zero thickness
        # all thickness in nm

        d_list = np.zeros(len(currentStructure) + 2)
        d_list[0] = d_incidentMedium
        d_list[-1] = d_substrate

        indexXarray = np.zeros(len(currentStructure) + 2, dtype=complex)
        indexXarray[0] = 1.0  # air
        indexXarray[-1] = 1.5  # ~quartz/glass

        indexYarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexYarray[0] = 1.0  # air
        indexYarray[-1] = 1.5  # ~quartz/glass

        indexZarray = np.zeros(len(currentStructure) + 2, dtype=np.complex128)
        indexZarray[0] = 1.0  # air
        indexZarray[-1] = 1.5  # ~quartz/glass

        for ilayer, layer in enumerate(currentStructure):

            d_list[ilayer + 1] = currentStructure[ilayer][1]

        all_tensors = []

        for iwl, wvl in enumerate(wavelength):

            # multilayer composition

            for ilayer, layer in enumerate(currentStructure):

                currentIndexInIndices = int(materialAtLayer[ilayer])
                indexXarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][0][iwl]
                    + indicesForCalculation[currentIndexInIndices][1][iwl] * 1j
                )
                indexYarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][2][iwl]
                    + indicesForCalculation[currentIndexInIndices][3][iwl] * 1j
                )
                indexZarray[ilayer + 1] = (
                    indicesForCalculation[currentIndexInIndices][4][iwl]
                    + indicesForCalculation[currentIndexInIndices][5][iwl] * 1j
                )

            e_listx = indexXarray**2
            e_listy = indexYarray**2
            e_listz = indexZarray**2

            # optical constant tensor
            m_eps = np.zeros((len(e_listx), 3, 3), dtype=np.complex128)

            # filling dielectric tensor

            m_eps[:, 0, 0] = e_listx
            m_eps[:, 1, 1] = e_listy
            m_eps[:, 2, 2] = e_listz

            all_tensors.append(m_eps)

        total_T_s, total_T_p, total_t_s, total_t_p = self.calculateT(
            wavelength, v_theta, 0, all_tensors, d_list
        )

        return total_t_s, total_T_s

    def findIndexInUniqueMaterials(self, layerName, UniqueMaterials):

        for iUnique, uniqueName in enumerate(UniqueMaterials):

            if uniqueName == layerName:

                return iUnique
