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

    std::pair<std::vector<tk::spline>,bool> importIndexFromFile(const std::string &filename)
    {
        bool general_material;

        std::string fullfilename = filename + ".csv";
        std::filesystem::path fullpath = std::filesystem::current_path().parent_path() / "materials" / fullfilename;
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

            while (std::getline(ss, field, '\t'))
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
            general_material = true; // in-plane out of plane anisotropy + in-plane anisotropy
            while (std::getline(file, line))
            {
                std::istringstream ss(line);
                std::string field;
                std::getline(ss, field, '\t');
                wavelength.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                nx.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                kx.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                ny.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                ky.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                nz.push_back(std::stod(field));
                std::getline(ss, field, '\t');
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
            
            return std::make_pair(permittivity_splines,general_material);
        }

        if (columnCount == 5)
        {
            general_material = true; // in-plane out-of-plane anisotropy
            while (std::getline(file, line))
            {
                std::istringstream ss(line);
                std::string field;
                std::getline(ss, field, '\t');
                wavelength.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                nx.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                kx.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                nz.push_back(std::stod(field));
                std::getline(ss, field, '\t');
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
            
            return std::make_pair(permittivity_splines,general_material);
        }

        if (columnCount == 3)
        {
            general_material = false; // isotropic
            while (std::getline(file, line))
            {
                std::istringstream ss(line);
                std::string field;
                std::getline(ss, field, '\t');
                wavelength.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                nx.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                kx.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                ny.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                ky.push_back(std::stod(field));
                std::getline(ss, field, '\t');
                nz.push_back(std::stod(field));
                std::getline(ss, field, '\t');
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
            
            return std::make_pair(permittivity_splines,general_material);
        }

        tk::spline default_spline({0},{0});
        std::vector<tk::spline> default_splines(1,default_spline);
        return std::make_pair(default_splines,false);       
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
    std::vector<bool> incoherent;
    std::vector<double> structure_thicknesses;
    std::vector<double> azimuthalAngles;
    std::string polarization;
    std::string incidentMediumMaterial;
    std::string exitMediumMaterial;

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
        std::string polarization,
        std::vector<std::string> structure_materials,
        std::vector<bool> incoherent,
        std::vector<double> structure_thicknesses,
        std::string incidentMediumMaterial,
        std::string exitMediumMaterial) : calculation_type(calculation_type),
                       polarAngleMin(polarAngleMin),
                       polarAngleMax(polarAngleMax),
                       polarAngleStep(polarAngleStep),
                       azimAngleMin(azimAngleMin),
                       azimAngleMax(azimAngleMax),
                       azimAngleStep(azimAngleStep),
                       wavelengthMin(wavelengthMin),
                       wavelengthMax(wavelengthMax),
                       wavelengthStep(wavelengthStep),
                       polarization(polarization),
                       structure_materials(structure_materials),
                       incoherent(incoherent),
                       structure_thicknesses(structure_thicknesses),
                       incidentMediumMaterial(incidentMediumMaterial),
                       exitMediumMaterial(exitMediumMaterial) {}
};

CalculationInfo loadCalculationInfo(const std::filesystem::path filepath)
{
    std::ifstream file(filepath);
    json calculation_order;
    file >> calculation_order;

    std::vector<std::string> structure_materials = calculation_order["structure_materials"];
    std::vector<bool> incoherent = calculation_order["incoherent"];
    std::vector<double> structure_thicknesses = calculation_order["structure_thicknesses"];

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
        calculation_order["polarization"],
        structure_materials,
        incoherent,
        structure_thicknesses,
        calculation_order["incident_medium"],
        calculation_order["exit_medium"]
        );

    return calculation_info;
}

std::vector<std::string> getUniqueMembers(const std::vector<std::string> &inputArray)
{
    std::set<std::string> uniqueSet(inputArray.begin(), inputArray.end());
    std::vector<std::string> uniqueVector(uniqueSet.begin(), uniqueSet.end());
    return uniqueVector;
}