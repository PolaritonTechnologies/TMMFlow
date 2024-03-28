#include <fstream>
#include "json.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <complex>
#include "Eigen/Dense"
#include "spline.h"

using namespace Eigen;
using json = nlohmann::json;

class CSVParser
{
public:

    std::vector<tk::spline> permittivityFromIndex(const tk::spline& real_index_spline, const tk::spline& imag_index_spline, const std::vector<double>& wavelength)
    {
        std::vector<double> real_values, imag_values;

        for (double w : wavelength) {
            double n = real_index_spline(w);
            double k = imag_index_spline(w);

            real_values.push_back(n*n - k*k);
            imag_values.push_back(2*n*k);
        }

        // Create the new splines
        tk::spline real_spline(wavelength, real_values);
        tk::spline imag_spline(wavelength, imag_values);

        std::vector<tk::spline> permittivity_splines = {real_spline, imag_spline};
        return permittivity_splines;
    }

    std::vector<tk::spline> importIndexFromFile(const std::string &filename)
    {

        std::string fullfilename = filename + ".csv";
        std::filesystem::path fullpath = std::filesystem::current_path() / "Materials" / fullfilename;
        std::ifstream file(fullpath);

        std::vector<double> wavelength;
        std::vector<double> nx, kx;
        std::vector<double> ny, ky;
        std::vector<double> nz, kz;

        std::string line;
        std::getline(file, line);
        std::getline(file, line);

        int columnCount = 0;
        if (std::getline(file, line))
        {
            std::istringstream ss(line);
            std::string field;

            while (std::getline(ss, field, ','))
            {
                ++columnCount;
            }
        }

        file.clear();
        file.seekg(0);
        std::getline(file, line);
        std::getline(file, line);

        if (columnCount == 7)
        {
            while (std::getline(file, line))
            {
                std::istringstream ss(line);
                std::string field;
                std::getline(ss, field, ',');
                wavelength.push_back(std::stod(field));
                std::getline(ss, field, ',');
                nx.push_back(std::stod(field));
                std::getline(ss, field, ',');
                kx.push_back(std::stod(field));
                std::getline(ss, field, ',');
                ny.push_back(std::stod(field));
                std::getline(ss, field, ',');
                ky.push_back(std::stod(field));
                std::getline(ss, field, ',');
                nz.push_back(std::stod(field));
                std::getline(ss, field, ',');
                kz.push_back(std::stod(field));
            }

            tk::spline nx_spline(wavelength, nx);
            tk::spline ny_spline(wavelength, ny);
            tk::spline nz_spline(wavelength, nz);

            tk::spline kx_spline(wavelength, kx);
            tk::spline ky_spline(wavelength, ky);
            tk::spline kz_spline(wavelength, kz);

            std::vector<tk::spline> indices_splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};

            // Create the permittivity splines
            std::vector<tk::spline> permittivities_x_spline = permittivityFromIndex(nx_spline, kx_spline, wavelength);
            tk::spline real_x_spline = permittivities_x_spline[0];
            tk::spline imag_x_spline = permittivities_x_spline[1];

            std::vector<tk::spline> permittivities_y_spline = permittivityFromIndex(ny_spline, ky_spline, wavelength);
            tk::spline real_y_spline = permittivities_y_spline[0];
            tk::spline imag_y_spline = permittivities_y_spline[1];

            std::vector<tk::spline> permittivities_z_spline = permittivityFromIndex(nz_spline, kz_spline, wavelength);          
            tk::spline real_z_spline = permittivities_z_spline[0];
            tk::spline imag_z_spline = permittivities_z_spline[1];

            std::vector<tk::spline> permittivity_splines = {real_x_spline, imag_x_spline, real_y_spline, imag_y_spline, real_z_spline, imag_z_spline};

            //returns the permittivity tensors and not the refractive indices
            
            return permittivity_splines;
        }

        if (columnCount == 5)
        {
            while (std::getline(file, line))
            {
                std::istringstream ss(line);
                std::string field;
                std::getline(ss, field, ',');
                wavelength.push_back(std::stod(field));
                std::getline(ss, field, ',');
                nx.push_back(std::stod(field));
                std::getline(ss, field, ',');
                kx.push_back(std::stod(field));
                std::getline(ss, field, ',');
                nz.push_back(std::stod(field));
                std::getline(ss, field, ',');
                kz.push_back(std::stod(field));
            }

            tk::spline nx_spline(wavelength, nx); // ordinary extraordinary only
            tk::spline ny_spline(wavelength, nx);
            tk::spline nz_spline(wavelength, nz);

            tk::spline kx_spline(wavelength, kx); // ordinary extraordinary only
            tk::spline ky_spline(wavelength, kx);
            tk::spline kz_spline(wavelength, kz);

            std::vector<tk::spline> indices_splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};

            // Create the permittivity splines
            std::vector<tk::spline> permittivities_x_spline = permittivityFromIndex(nx_spline, kx_spline, wavelength);
            tk::spline real_x_spline = permittivities_x_spline[0];
            tk::spline imag_x_spline = permittivities_x_spline[1];

            std::vector<tk::spline> permittivities_z_spline = permittivityFromIndex(nz_spline, kz_spline, wavelength);          
            tk::spline real_z_spline = permittivities_z_spline[0];
            tk::spline imag_z_spline = permittivities_z_spline[1];

            std::vector<tk::spline> permittivity_splines = {real_x_spline, imag_x_spline, real_x_spline, imag_x_spline, real_z_spline, imag_z_spline};

            //returns the permittivity tensors and not the refractive indices
            
            return permittivity_splines;
        }

        if (columnCount == 3)
        {
            while (std::getline(file, line))
            {
                std::istringstream ss(line);
                std::string field;
                std::getline(ss, field, ',');
                wavelength.push_back(std::stod(field));
                std::getline(ss, field, ',');
                nx.push_back(std::stod(field));
                std::getline(ss, field, ',');
                kx.push_back(std::stod(field));
                std::getline(ss, field, ',');
                ny.push_back(std::stod(field));
                std::getline(ss, field, ',');
                ky.push_back(std::stod(field));
                std::getline(ss, field, ',');
                nz.push_back(std::stod(field));
                std::getline(ss, field, ',');
                kz.push_back(std::stod(field));
            }
            tk::spline nx_spline(wavelength, nx); // isotropic

            tk::spline kx_spline(wavelength, kx); // isotropic

            std::vector<tk::spline> indices_splines = {nx_spline, nx_spline, nx_spline, kx_spline, kx_spline, kx_spline};

            // Create the permittivity splines
            std::vector<tk::spline> permittivities_x_spline = permittivityFromIndex(nx_spline, kx_spline, wavelength);
            tk::spline real_x_spline = permittivities_x_spline[0];
            tk::spline imag_x_spline = permittivities_x_spline[1];

            tk::spline real_y_spline = permittivities_x_spline[0];
            tk::spline imag_y_spline = permittivities_x_spline[1];
          
            tk::spline real_z_spline = permittivities_x_spline[0];
            tk::spline imag_z_spline = permittivities_x_spline[1];

            std::vector<tk::spline> permittivity_splines = {real_x_spline, imag_x_spline, real_y_spline, imag_y_spline, real_z_spline, imag_z_spline};

            //returns the permittivity tensors and not the refractive indices
            
            return permittivity_splines;
        }

        tk::spline default_spline({0},{0});
        return {default_spline};
            
    }
        
};

struct CalculationInfo
{
    std::string calculation_type;
    double polarAngleMin;
    double polarAngleMax;
    double polarAngleStep;
    double azimAngleMin;
    double azimAngleMax;
    double azimAngleStep;
    double wavelengthMin;
    double wavelengthMax;
    double wavelengthStep;
    std::vector<std::string> structure_materials;
    std::vector<double> structure_thicknesses;
    std::string polarization;
    std::vector<double> azimuthalAngles;

    // Default constructor
    CalculationInfo() = default;

    CalculationInfo(
        std::string calculation_type,
        double polarAngleMin,
        double polarAngleMax,
        double polarAngleStep,
        double azimAngleMin,
        double azimAngleMax,
        double azimAngleStep,
        double wavelengthMin,
        double wavelengthMax,
        double wavelengthStep,
        std::vector<std::string> structure_materials,
        std::vector<double> structure_thicknesses,
        std::string polarization) : calculation_type(calculation_type),
                       polarAngleMin(polarAngleMin),
                       polarAngleMax(polarAngleMax),
                       polarAngleStep(polarAngleStep),
                       azimAngleMin(azimAngleMin),
                       azimAngleMax(azimAngleMax),
                       azimAngleStep(azimAngleStep),
                       wavelengthMin(wavelengthMin),
                       wavelengthMax(wavelengthMax),
                       wavelengthStep(wavelengthStep),
                       structure_materials(structure_materials),
                       structure_thicknesses(structure_thicknesses),
                       polarization(polarization) {}
};

CalculationInfo loadCalculationInfo(const std::filesystem::path filepath)
{
    std::ifstream file(filepath);
    json calculation_order;
    file >> calculation_order;

    // Reverse "structure_materials" and "structure_thicknesses"
    std::vector<std::string> structure_materials = calculation_order["structure_materials"];
    std::reverse(structure_materials.begin(), structure_materials.end());

    std::vector<double> structure_thicknesses = calculation_order["structure_thicknesses"];
    std::reverse(structure_thicknesses.begin(), structure_thicknesses.end());

    CalculationInfo calculation_info(
        calculation_order["calculation_type"],
        calculation_order["polarAngleMin"],
        calculation_order["polarAngleMax"],
        calculation_order["polarAngleStep"],
        calculation_order["azimAngleMin"],
        calculation_order["azimAngleMax"],
        calculation_order["azimAngleStep"],
        calculation_order["wavelengthMin"],
        calculation_order["wavelengthMax"],
        calculation_order["wavelengthStep"],
        structure_materials,
        structure_thicknesses,
        calculation_order["polarization"]
        );

    return calculation_info;
}

std::vector<std::string> getUniqueMembers(const std::vector<std::string> &inputArray)
{
    std::set<std::string> uniqueSet(inputArray.begin(), inputArray.end());
    std::vector<std::string> uniqueVector(uniqueSet.begin(), uniqueSet.end());
    return uniqueVector;
}

std::vector<Matrix3cd> reverseVector(const std::vector<Matrix3cd> &vec)
{
    std::vector<Matrix3cd> reversedVec(vec.rbegin(), vec.rend());
    return reversedVec;
}