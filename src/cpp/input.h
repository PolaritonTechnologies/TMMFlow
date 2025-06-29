#pragma once
#include <string>
#include <vector>
#include "spline.h"
#include <filesystem>
#include <fstream>
#include <sqlite3.h>
#include <iostream>
#include <nlohmann/json.hpp>
std::filesystem::path getBasePath();
std::filesystem::path getDatabasePath(const std::string &filename);

class CSVParser
{
public:
    CSVParser(const std::string &db_filename)
    {
        std::filesystem::path fullpath = getDatabasePath(db_filename);
        std::cout << "Database file: " << fullpath << std::endl;
        int rc = sqlite3_open(fullpath.c_str(), &db);
        if (rc)
        {
            std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
            db = nullptr;
        }
    }

    ~CSVParser()
    {
        if (db)
        {
            sqlite3_close(db);
        }
    }

    // The CSVParser methods are defined here instead of cpp due to difficulties in linking the spline header
    std::vector<tk::spline> permittivityFromIndex(const tk::spline &real_index_spline, const tk::spline &imag_index_spline, const std::vector<double> &wavelength)
    {
        std::vector<double> real_values, imag_values;
        for (double w : wavelength)
        {
            double n = real_index_spline(w);
            double k = imag_index_spline(w);

            real_values.push_back(n * n - k * k);
            imag_values.push_back(2 * n * k);
        }

        // Create the new splines
        tk::spline real_spline(wavelength, real_values);
        tk::spline imag_spline(wavelength, imag_values);

        std::vector<tk::spline> permittivity_splines = {real_spline, imag_spline};
        return permittivity_splines;
    }

    std::pair<std::vector<tk::spline>, bool> importIndexFromFile(const std::string &material_name)
    {
        if (!db)
        {
            std::cerr << "Database not initialized" << std::endl;
            return {{}, false};
        }

        // Prepare the SQL query
        const char *sql = "SELECT data FROM materials WHERE name = ? AND active = 1";
        sqlite3_stmt *stmt;
        int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK)
        {
            std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
            return {{}, false};
        }

        // Bind the material name parameter
        sqlite3_bind_text(stmt, 1, material_name.c_str(), -1, SQLITE_STATIC);

        // Execute the query and process the result
        std::string json_data;
        if ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
        {
            json_data = reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0));
        }
        else
        {
            std::cerr << "Failed to execute statement: " << sqlite3_errmsg(db) << std::endl;
            sqlite3_finalize(stmt);
            return {{}, false};
        }

        // Finalize the statement
        sqlite3_finalize(stmt);

        // Parse the JSON data
        nlohmann::json json;
        try
        {
            json = nlohmann::json::parse(json_data);
        }
        catch (const nlohmann::json::parse_error &e)
        {
            std::cerr << "JSON parse error: " << e.what() << std::endl;
            return {{}, false};
        }

        // Extract keys from the JSON object
        std::vector<std::string> keys;
        for (auto it = json.begin(); it != json.end(); ++it)
        {
            keys.push_back(it.key());
            // std::cout << it.key() << std::endl;
        }

        if (keys.size() == 3)
        {
            bool general_material = false;
            // std::cout << "0" << std::endl;
            std::vector<double> wavelength = json["wavelength"].get<std::vector<double>>();
            std::vector<double> nx = json["n"].get<std::vector<double>>();
            std::vector<double> kx = json["k"].get<std::vector<double>>();
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

            // returns the permittivity tensors and not the refractive indices

            return std::make_pair(permittivity_splines, general_material);
        }
        else if (keys.size() == 5)
        {
            // std::cout << "1" << std::endl;
            bool general_material = true;
            std::vector<double> wavelength = json["wavelength"].get<std::vector<double>>();
            std::vector<double> nx = json["n_ord"].get<std::vector<double>>();
            std::vector<double> kx = json["k_ord"].get<std::vector<double>>();
            std::vector<double> nz = json["n_exord"].get<std::vector<double>>();
            std::vector<double> kz = json["k_exord"].get<std::vector<double>>();

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

            // returns the permittivity tensors and not the refractive indices

            return std::make_pair(permittivity_splines, general_material);
        }
        else if (keys.size() == 7)
        {
            // std::cout << "2" << std::endl;
            bool general_material = true;
            std::vector<double> wavelength = json["wavelength"].get<std::vector<double>>();
            std::vector<double> nx = json["n_x"].get<std::vector<double>>();
            std::vector<double> kx = json["k_x"].get<std::vector<double>>();
            std::vector<double> ny = json["n_y"].get<std::vector<double>>();
            std::vector<double> ky = json["k_y"].get<std::vector<double>>();
            std::vector<double> nz = json["n_z"].get<std::vector<double>>();
            std::vector<double> kz = json["k_z"].get<std::vector<double>>();
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

            // returns the permittivity tensors and not the refractive indices
            return std::make_pair(permittivity_splines, general_material);
        }
        else
        {
            std::cerr << "Invalid number of keys for material: " << keys.size() << std::endl;
            return {{}, false};
            // tk::spline default_spline({0}, {0});
            // std::vector<tk::spline> default_splines(1, default_spline);
            // return std::make_pair(default_splines, false);
        }

        // // print wavelength, n and k values
    }

private:
    sqlite3 *db;
};

// std::pair<std::vector<tk::spline>, bool> importIndexFromFile(const std::string &filename)
// {
//     bool general_material;

//     std::filesystem::path fullpath = getMaterialsPath(filename);
//     std::ifstream file(fullpath);

//     std::vector<double> wavelength;
//     std::vector<double> nx, kx;
//     std::vector<double> ny, ky;
//     std::vector<double> nz, kz;

//     std::string line;
//     std::getline(file, line);
//     std::getline(file, line);

//     int columnCount = 0;
//     if (std::getline(file, line))
//     {
//         std::istringstream ss(line);
//         std::string field;

//         while (std::getline(ss, field, '\t'))
//         {
//             ++columnCount;
//         }
//     }

//     file.clear();
//     file.seekg(0);
//     std::getline(file, line);
//     std::getline(file, line);

//     if (columnCount == 7)
//     {
//         general_material = true; // in-plane out of plane anisotropy + in-plane anisotropy
//         while (std::getline(file, line))
//         {
//             std::istringstream ss(line);
//             std::string field;
//             std::getline(ss, field, '\t');
//             wavelength.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             nx.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             kx.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             ny.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             ky.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             nz.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             kz.push_back(std::stod(field));
//         }

//         tk::spline nx_spline(wavelength, nx);
//         tk::spline ny_spline(wavelength, ny);
//         tk::spline nz_spline(wavelength, nz);

//         tk::spline kx_spline(wavelength, kx);
//         tk::spline ky_spline(wavelength, ky);
//         tk::spline kz_spline(wavelength, kz);

//         std::vector<tk::spline> indices_splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};

//         // Create the permittivity splines
//         std::vector<tk::spline> permittivities_x_spline = permittivityFromIndex(nx_spline, kx_spline, wavelength);
//         tk::spline real_x_spline = permittivities_x_spline[0];
//         tk::spline imag_x_spline = permittivities_x_spline[1];

//         std::vector<tk::spline> permittivities_y_spline = permittivityFromIndex(ny_spline, ky_spline, wavelength);
//         tk::spline real_y_spline = permittivities_y_spline[0];
//         tk::spline imag_y_spline = permittivities_y_spline[1];

//         std::vector<tk::spline> permittivities_z_spline = permittivityFromIndex(nz_spline, kz_spline, wavelength);
//         tk::spline real_z_spline = permittivities_z_spline[0];
//         tk::spline imag_z_spline = permittivities_z_spline[1];

//         std::vector<tk::spline> permittivity_splines = {real_x_spline, imag_x_spline, real_y_spline, imag_y_spline, real_z_spline, imag_z_spline};

//         // returns the permittivity tensors and not the refractive indices

//         return std::make_pair(permittivity_splines, general_material);
//     }

//     if (columnCount == 5)
//     {
//         general_material = true; // in-plane out-of-plane anisotropy
//         while (std::getline(file, line))
//         {
//             std::istringstream ss(line);
//             std::string field;
//             std::getline(ss, field, '\t');
//             wavelength.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             nx.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             kx.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             nz.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             kz.push_back(std::stod(field));
//         }

//         tk::spline nx_spline(wavelength, nx); // ordinary extraordinary only
//         tk::spline ny_spline(wavelength, nx);
//         tk::spline nz_spline(wavelength, nz);

//         tk::spline kx_spline(wavelength, kx); // ordinary extraordinary only
//         tk::spline ky_spline(wavelength, kx);
//         tk::spline kz_spline(wavelength, kz);

//         std::vector<tk::spline> indices_splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};

//         // Create the permittivity splines
//         std::vector<tk::spline> permittivities_x_spline = permittivityFromIndex(nx_spline, kx_spline, wavelength);
//         tk::spline real_x_spline = permittivities_x_spline[0];
//         tk::spline imag_x_spline = permittivities_x_spline[1];

//         std::vector<tk::spline> permittivities_z_spline = permittivityFromIndex(nz_spline, kz_spline, wavelength);
//         tk::spline real_z_spline = permittivities_z_spline[0];
//         tk::spline imag_z_spline = permittivities_z_spline[1];

//         std::vector<tk::spline> permittivity_splines = {real_x_spline, imag_x_spline, real_x_spline, imag_x_spline, real_z_spline, imag_z_spline};

//         // returns the permittivity tensors and not the refractive indices

//         return std::make_pair(permittivity_splines, general_material);
//     }

//     if (columnCount == 3)
//     {
//         general_material = false; // isotropic
//         while (std::getline(file, line))
//         {
//             std::istringstream ss(line);
//             std::string field;
//             std::getline(ss, field, '\t');
//             wavelength.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             nx.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             kx.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             ny.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             ky.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             nz.push_back(std::stod(field));
//             std::getline(ss, field, '\t');
//             kz.push_back(std::stod(field));
//         }
//         tk::spline nx_spline(wavelength, nx); // isotropic

//         tk::spline kx_spline(wavelength, kx); // isotropic

//         std::vector<tk::spline> indices_splines = {nx_spline, nx_spline, nx_spline, kx_spline, kx_spline, kx_spline};

//         // Create the permittivity splines
//         std::vector<tk::spline> permittivities_x_spline = permittivityFromIndex(nx_spline, kx_spline, wavelength);
//         tk::spline real_x_spline = permittivities_x_spline[0];
//         tk::spline imag_x_spline = permittivities_x_spline[1];

//         tk::spline real_y_spline = permittivities_x_spline[0];
//         tk::spline imag_y_spline = permittivities_x_spline[1];

//         tk::spline real_z_spline = permittivities_x_spline[0];
//         tk::spline imag_z_spline = permittivities_x_spline[1];

//         std::vector<tk::spline> permittivity_splines = {real_x_spline, imag_x_spline, real_y_spline, imag_y_spline, real_z_spline, imag_z_spline};

//         // returns the permittivity tensors and not the refractive indices

//         return std::make_pair(permittivity_splines, general_material);
//     }
//
// tk::spline default_spline({0}, {0});
// std::vector<tk::spline> default_splines(1, default_spline);
// return std::make_pair(default_splines, false);

struct CalculationInfo
{
public:
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
    double polarization;
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
        double polarization,
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

CalculationInfo loadCalculationInfo(const std::string &json_text);
std::vector<std::string> getUniqueMembers(const std::vector<std::string> &inputArray);