#include <iostream>
#include <fstream>
#include "Eigen\Dense"
#include "spline.h"

using namespace Eigen;

class CSVParser
{
public:
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

            std::vector<tk::spline> splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};
            return splines;
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
                ny.push_back(std::stod(field));
                std::getline(ss, field, ',');
                ky.push_back(std::stod(field));
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

            std::vector<tk::spline> splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};
            return splines;
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
            tk::spline ny_spline(wavelength, nx);
            tk::spline nz_spline(wavelength, nx);

            tk::spline kx_spline(wavelength, kx); // isotropic
            tk::spline ky_spline(wavelength, kx);
            tk::spline kz_spline(wavelength, kx);

            std::vector<tk::spline> splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};
            return splines;
        }

        tk::spline nx_spline({0}, {0}); // isotropic
        tk::spline ny_spline({0}, {0});
        tk::spline nz_spline({0}, {0});

        tk::spline kx_spline({0}, {0}); // isotropic
        tk::spline ky_spline({0}, {0});
        tk::spline kz_spline({0}, {0});

        std::vector<tk::spline> splines = {nx_spline, ny_spline, nz_spline, kx_spline, ky_spline, kz_spline};
        return splines;
    }
};