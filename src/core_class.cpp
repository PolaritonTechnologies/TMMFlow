#include "input.h"
#include "core.h"

class FilterStack
{
public:
    std::filesystem::path full_path;
    std::map<std::string, std::vector<tk::spline>> material_splines;
    CalculationInfo calculation_order;

    // Default constructor
    FilterStack() = default;

    // Constructor
    FilterStack(std::string filename)
    {
        full_path = std::filesystem::current_path() / filename;
        material_splines = assemble_materials(filename);
        calculation_order = loadCalculationInfo(full_path);
    }

    void calculate_transmission_reflection_matrices(double wavelength, double theta_0, double phi_0)
    {
        std::vector<double> d_list = calculation_order.structure_thicknesses;

        auto [m_r_ps, m_t_ps] = calculate_tr(assemble_e_list_3x3(material_splines, wavelength), d_list, wavelength, theta_0, phi_0);

        double reflection = std::sqrt(std::pow(m_r_ps(0, 0).real(), 2) + std::pow(m_r_ps(1, 1).real(), 2));
        double transmission = std::sqrt(std::pow(m_t_ps(0, 0).real(), 2) + std::pow(m_t_ps(1, 1).real(), 2));
    }

    std::array<std::map<std::pair<double, double>, std::complex<double>>, 4> calculate_transmission_reflection_for_range(){};

private:
    std::string filename = "calculation_order.json";

    // Function to assemble materials from file
    std::map<std::string, std::vector<tk::spline>> assemble_materials(std::string full_path)
    {
        CalculationInfo calculation_order = loadCalculationInfo(full_path);

        std::vector<double> d_list = calculation_order.structure_thicknesses;

        std::vector<std::string> unique_materials = getUniqueMembers(calculation_order.structure_materials);
        CSVParser parser;
        std::map<std::string, std::vector<tk::spline>> material_splines;

        for (const auto &material : unique_materials)
        {
            std::vector<tk::spline> splines = parser.importIndexFromFile(material);
            material_splines[material] = splines;
        }

        return material_splines;
    }

    // Function to assemble e_list_3x3 for a given wavelength using the material_splines obtained from file
    std::vector<Matrix3cd> assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength)
    {
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
                (Matrix3cd(3, 3) << std::complex<double>(material_splines[material][0](wavelength), material_splines[material][1](wavelength)), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][2](wavelength), material_splines[material][3](wavelength)), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][4](wavelength), material_splines[material][5](wavelength)))
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

        return e_list_3x3;
    }
};