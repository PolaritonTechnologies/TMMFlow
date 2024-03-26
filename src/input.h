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

            std::ofstream filen("nz.csv");
            for(const auto& n : nz) {
                filen << n << "\n";
            }
            filen.close();

            std::ofstream filek("kz.csv");
            for(const auto& k : kz) {
                filek << k << "\n";
            }
            filek.close();

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
    double angleMin;
    double angleMax;
    double angleStep;
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
        double angleMin,
        double angleMax,
        double angleStep,
        double wavelengthMin,
        double wavelengthMax,
        double wavelengthStep,
        std::vector<std::string> structure_materials,
        std::vector<double> structure_thicknesses,
        std::string polarization,
        std::vector<double> azimuthalAngles) : calculation_type(calculation_type),
                       angleMin(angleMin),
                       angleMax(angleMax),
                       angleStep(angleStep),
                       wavelengthMin(wavelengthMin),
                       wavelengthMax(wavelengthMax),
                       wavelengthStep(wavelengthStep),
                       structure_materials(structure_materials),
                       structure_thicknesses(structure_thicknesses),
                       polarization(polarization),
                       azimuthalAngles(azimuthalAngles) {}
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
        calculation_order["angleMin"],
        calculation_order["angleMax"],
        calculation_order["angleStep"],
        calculation_order["wavelengthMin"],
        calculation_order["wavelengthMax"],
        calculation_order["wavelengthStep"],
        structure_materials,
        structure_thicknesses,
        calculation_order["polarization"],
        calculation_order["azimuthalAngles"]
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

/*
std::array<std::map<std::pair<double, double>, std::complex<double>>, 4> perform_calculation()
{
    std::string filename = "calculation_order.json";
    std::filesystem::path fullPath = std::filesystem::current_path() / filename;
    std::cout << fullPath << std::endl;
    CalculationInfo calculation_order = loadCalculationInfo(fullPath);

    std::vector<std::complex<double>> e_listx_wvl, e_listy_wvl, e_listz_wvl;

    double theta_min = calculation_order.angleMin * M_PI / 180.0;
    double theta_max = calculation_order.angleMax * M_PI / 180.0;
    double theta_stepInRadians = calculation_order.angleStep * M_PI / 180.0;
    int theta_steps = (theta_max - theta_min) / theta_stepInRadians + 1;
    VectorXd v_theta = VectorXd::LinSpaced(theta_steps, theta_min, theta_max);

    double wvl_min = calculation_order.wavelengthMin;
    double wvl_max = calculation_order.wavelengthMax;
    int wvl_steps = (wvl_max - wvl_min) / calculation_order.wavelengthStep + 1;
    VectorXd wavelength = VectorXd::LinSpaced(wvl_steps, wvl_min, wvl_max);

    // Phi initially kept at 0 for now
    double phi_0 = 0;

    std::vector<double> d_list = calculation_order.structure_thicknesses;

    std::vector<std::string> unique_materials = getUniqueMembers(calculation_order.structure_materials);
    CSVParser parser;
    std::map<std::string, std::vector<tk::spline>> materialSplines;

    for (const auto &material : unique_materials)
    {
        std::vector<tk::spline> splines = parser.importIndexFromFile(material);
        materialSplines[material] = splines;
    }

    std::map<std::pair<double, double>, std::complex<double>> resultDictionary_reflection_s;
    std::map<std::pair<double, double>, std::complex<double>> resultDictionary_reflection_p;
    std::map<std::pair<double, double>, std::complex<double>> resultDictionary_transmission_s;
    std::map<std::pair<double, double>, std::complex<double>> resultDictionary_transmission_p;

    for (int p = 0; p < v_theta.size(); p++)
    {
        for (int i = 0; i < wavelength.size(); i++)
        {

            std::pair<double, double> key = std::make_pair(v_theta[p] * 180.0 / M_PI, wavelength[i]);

            std::vector<Matrix3cd> e_list_3x3;

            // substrate: quartz/glass
            Matrix3cd substrate_layer_tensor =
                (Matrix3cd(3, 3) << std::complex<double>(1.5, 0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(1.5, 0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.5, 0))
                    .finished();
            e_list_3x3.push_back(substrate_layer_tensor);

            for (const auto &material : calculation_order.structure_materials)
            {
                Matrix3cd next_layer_tensor =
                    (Matrix3cd(3, 3) << std::complex<double>(materialSplines[material][0](wavelength[i]), materialSplines[material][1](wavelength[i])), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
                     std::complex<double>(0.0, 0.0), std::complex<double>(materialSplines[material][2](wavelength[i]), materialSplines[material][3](wavelength[i])), std::complex<double>(0.0, 0.0),
                     std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(materialSplines[material][4](wavelength[i]), materialSplines[material][5](wavelength[i])))
                        .finished();
                e_list_3x3.push_back(next_layer_tensor);
            }

            // out medium: air
            Matrix3cd incident_layer_tensor =
                (Matrix3cd(3, 3) << std::complex<double>(1, 0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(1, 0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1, 0))
                    .finished();
            e_list_3x3.push_back(incident_layer_tensor);

            auto [m_r_ps, m_t_ps] = calculate_tr(e_list_3x3, d_list, wavelength[i], v_theta[p], phi_0);

            resultDictionary_reflection_p[key] = m_r_ps.cwiseAbs2()(0, 0);
            resultDictionary_reflection_s[key] = m_r_ps.cwiseAbs2()(1, 1);
            resultDictionary_transmission_p[key] = m_t_ps.cwiseAbs2()(0, 0);
            resultDictionary_transmission_s[key] = m_t_ps.cwiseAbs2()(1, 1);

            // std::cout << "reflection matrix: " << m_r_ps << std::endl;
            // std::cout << "transmission matrix: " << m_t_ps << std::endl;
        }
    }

    std::array<std::map<std::pair<double, double>, std::complex<double>>, 4> resultArray = {resultDictionary_reflection_p, resultDictionary_reflection_s, resultDictionary_transmission_p, resultDictionary_transmission_s};

    return resultArray;
}
*/

void dump_to_file(const std::map<std::pair<double, double>, std::complex<double>> &dictionary, const std::string &filename)
{
    std::ofstream file(filename);

    // Create maps of unique angles and wavelengths to indices
    std::map<double, int> angleIndices;
    std::map<double, int> wavelengthIndices;
    for (const auto &[key, value] : dictionary)
    {
        angleIndices[key.first] = 0;
        wavelengthIndices[key.second] = 0;
    }

    // Assign indices to unique angles and wavelengths
    int index = 0;
    for (auto &[angle, angleIndex] : angleIndices)
    {
        angleIndex = index++;
    }
    index = 0;
    for (auto &[wavelength, wavelengthIndex] : wavelengthIndices)
    {
        wavelengthIndex = index++;
    }

    // Create 2D grid
    std::vector<std::vector<std::complex<double>>> grid(wavelengthIndices.size(), std::vector<std::complex<double>>(angleIndices.size()));

    // Fill grid with values from dictionary
    for (const auto &[key, value] : dictionary)
    {
        int rowIndex = wavelengthIndices[key.second];
        int columnIndex = angleIndices[key.first];
        grid[rowIndex][columnIndex] = value.real();
    }

    // Write the first line (list of unique angles)
    for (const auto &[angle, angleIndex] : angleIndices)
    {
        file << "," << angle;
    }
    file << "\n";

    // Write the remaining lines (wavelength and corresponding values)
    for (const auto &[wavelength, rowIndex] : wavelengthIndices)
    {
        file << wavelength;
        for (const auto &value : grid[rowIndex])
        {
            file << "," << value.real();
        }
        file << "\n";
    }

    file.close();
}

/*
int execute()
{
    std::array<std::map<std::pair<double, double>, std::complex<double>>, 4> resultArray = perform_calculation();
    dump_to_file(resultArray[0], "R_p.txt");
    dump_to_file(resultArray[1], "R_s.txt");
    dump_to_file(resultArray[2], "T_p.txt");
    dump_to_file(resultArray[3], "T_s.txt");
    return 0;
};
*/