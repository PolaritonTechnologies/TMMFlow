import csv
import time
import numpy
from itertools import tee
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline as us


class csvparser:

    def importFromFile(self, filename):

        # this function is made specifically to import .csv files from COMPLETE EASE
        # with two header lines, it will not work otherwise

        # spline precision:

        sdef = 0.0001

        wavelength = []

        nx = []
        kx = []

        ny = []
        ky = []

        nz = []
        kz = []

        nameAndExtension = filename + ".csv"

        with open(str("Materials\{}".format(nameAndExtension)), "r") as csvfile:

            readerColumn, reader = tee(csv.reader(csvfile, delimiter=","))

            numberOfcolumns = len(next(readerColumn))

            del readerColumn

            if numberOfcolumns <= 3:

                # the material studied is isotropic, only one set of n and k

                for irow, row in enumerate(reader):

                    # ignore the extensions on the two first lines

                    if irow > 1:

                        wavelength.append(float(row[0]))

                        nx.append(float(row[1]))
                        kx.append(float(row[2]))

                wavelength = numpy.asarray(wavelength)

                splnx = us(wavelength, numpy.asarray(nx), s=sdef)
                splkx = us(wavelength, self.cleanForZeros(numpy.asarray(kx)), s=sdef)

                splny = splnx
                splky = splkx

                splnz = splnx
                splkz = splkx

                # this particular step look for zeros in kz to avoid crashing the algorithm,
                # it replaces them with 1e-12

                # plt.plot(wavelength,splnx)
                # plt.plot(wavelength,splkx)
                # plt.plot(wavelength,splny)
                # plt.plot(wavelength,splky)
                # plt.plot(wavelength,splnz)
                # plt.plot(wavelength,splkz)

                return splnx, splkx, splny, splky, splnz, splkz

            if numberOfcolumns <= 5:

                # the material studied is uniaxial, nxy, kxy and nz, kz

                for irow, row in enumerate(reader):

                    # ignore the extensions on the two first lines

                    if irow > 1:

                        wavelength.append(float(row[0]))

                        nx.append(float(row[1]))
                        kx.append(float(row[2]))

                        nz.append(float(row[3]))
                        kz.append(float(row[4]))

                wavelength = numpy.asarray(wavelength)

                splnx = us(wavelength, numpy.asarray(nx), s=sdef)
                splkx = us(wavelength, self.cleanForZeros(numpy.asarray(kx)), s=sdef)

                splny = splnx
                splky = splkx

                splnz = us(wavelength, numpy.asarray(nz), s=sdef)
                splkz = us(wavelength, self.cleanForZeros(numpy.asarray(kz)), s=sdef)

                # this particular step look for zeros in kz to avoid crashing the algorithm,
                # it replaces them with 1e-12

                # plt.plot(wavelength,splnx)
                # plt.plot(wavelength,splkx)
                # plt.plot(wavelength,splny)
                # plt.plot(wavelength,splky)
                # plt.plot(wavelength,splnz)
                # plt.plot(wavelength,splkz)

                return splnx, splkx, splny, splky, splnz, splkz

            else:

                # the material studied is biaxial, nx ny nz and kx ky kz

                for irow, row in enumerate(reader):

                    # ignore the extensions on the two first lines

                    if irow > 1:

                        wavelength.append(float(row[0]))

                        nx.append(float(row[1]))
                        kx.append(float(row[2]))

                        ny.append(float(row[3]))
                        ky.append(float(row[4]))

                        nz.append(float(row[5]))
                        kz.append(float(row[6]))

                wavelength = numpy.asarray(wavelength)

                splnx = us(wavelength, numpy.asarray(nx), s=sdef)
                splkx = us(wavelength, self.cleanForZeros(numpy.asarray(kx)), s=sdef)

                splny = us(wavelength, numpy.asarray(ny), s=sdef)
                splky = us(wavelength, self.cleanForZeros(numpy.asarray(ky)), s=sdef)

                splnz = us(wavelength, numpy.asarray(nz), s=sdef)
                splkz = us(wavelength, self.cleanForZeros(numpy.asarray(kz)), s=sdef)

                # this particular step look for zeros in kz to avoid crashing the algorithm,
                # it replaces them with 1e-12

                # plt.plot(wavelength,splnx)
                # plt.plot(wavelength,splkx)
                # plt.plot(wavelength,splny)
                # plt.plot(wavelength,splky)
                # plt.plot(wavelength,splnz)
                # plt.plot(wavelength,splkz)

                return splnx, splkx, splny, splky, splnz, splkz

    def cleanForZeros(self, arrayToClean: numpy.array):

        for irow, element in enumerate(arrayToClean):

            if element == 0:

                arrayToClean[irow] = 1e-12

        return arrayToClean

    def indexAccordingToPromptedWavelengths(
        self, wavelength, splnx, splkx, splny, splky, splnz, splkz
    ):

        result = []

        result.append(splnx(wavelength))
        result.append(splkx(wavelength))
        result.append(splny(wavelength))
        result.append(splky(wavelength))
        result.append(splnz(wavelength))
        result.append(splkz(wavelength))

        return result

    def write(
        self, filepath, wavelength, realX, complexX, realY, complexY, realZ, complexZ
    ):

        falseCounter = 2

        if realY is not False:

            falseCounter = 4

        if realZ is not False:

            falseCounter = 6

        write = numpy.zeros((numpy.size(wavelength), falseCounter + 1))

        for eltIdx, elt in enumerate(wavelength):

            if falseCounter == 6:

                write[eltIdx, 5] = realZ[eltIdx]

                write[eltIdx, 6] = complexZ[eltIdx]

            if falseCounter > 2:

                write[eltIdx, 3] = realY[eltIdx]

                write[eltIdx, 4] = complexY[eltIdx]

            if falseCounter > 0:

                write[eltIdx, 0] = wavelength[eltIdx]

                write[eltIdx, 1] = realX[eltIdx]

                write[eltIdx, 2] = complexX[eltIdx]

        numpy.savetxt(filepath, write, delimiter=",")

        temprows = []

        with open(filepath, "r") as csvfile:

            reader = csv.reader(csvfile)

            for rowRead in reader:

                temprows.append(rowRead)

        f = open(filepath, "w", newline="")

        with f:

            writer = csv.writer(f)

            if falseCounter <= 2:

                writer.writerow(["Opt. Consts vs nm", "", ""])
                writer.writerow(["nm", "n", "k"])

            if falseCounter > 2 and falseCounter <= 4:

                writer.writerow(["Opt. Consts vs nm", "", "", "", ""])
                writer.writerow(["nm", "nord", "kord", "nex", "kex"])

            if falseCounter > 4:

                writer.writerow(["Opt. Consts vs nm", "", "", "", "", "", ""])
                writer.writerow(["nm", "nx", "kx", "ny", "ky", "nz", "kz"])

            for row in temprows:

                writer.writerow(row)


if __name__ == "__main__":

    local_parser = csvparser()

    local_parser.importFromFile("LCE3.25neff1.6Go7e12Ga1e12")
