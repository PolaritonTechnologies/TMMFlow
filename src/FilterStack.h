#include "core.h"
#include "input.h"
#include <cstring>

/*
 * Class to represent a filter stack defined by a an input file
 */
class FilterStack
{
public:
    // Public variables
    std::filesystem::path full_path;
    std::map<std::string, std::vector<tk::spline>> material_splines;
    CalculationInfo calculation_order;

    // Save the current material order as int positions of the layers
    std::vector<int> material_order_int;

    std::vector<double> d_list_initial;
    std::vector<std::string> material_order_initial;

    // Public methods
    // std::pair<Matrix2cd, Matrix2cd> calculate_transmission_reflection_matrices(double wavelength, double theta_0, double phi_0);
    double calculate_reflection_transmission_absorption(const char *type, const char *polarization, double wavelength, double theta_0, double phi_0); //, std::vector<double> d_list);
    // void calculate_and_save_ar_reflection();
    void change_material_order(std::vector<int> new_material_order);
    void change_material_thickness(std::vector<double> material_thickness);

    // Default constructor
    FilterStack() = default;

    // Constructor
    FilterStack(const char *filename)
    {
        // First convert the filename to an std::string
        full_path = std::filesystem::current_path() / filename;
        material_splines = assemble_materials(filename);
        calculation_order = loadCalculationInfo(full_path);
        d_list_initial = calculation_order.structure_thicknesses;
        material_order_initial = calculation_order.structure_materials;

        material_order_int.resize(calculation_order.structure_materials.size());

        for (int i = 0; i < calculation_order.structure_materials.size(); ++i)
        {
            material_order_int[i] = i;
        }
    }

    // The destructor should probably be populated more
    ~FilterStack() = default;

private:
    // Private methods
    std::map<std::string, std::vector<tk::spline>> assemble_materials(std::string full_path);
    std::vector<Matrix3cd> assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength);
};

/*
std::pair<Matrix2cd, Matrix2cd> FilterStack::calculate_transmission_reflection_matrices(double wavelength, double theta_0, double phi_0)
{
    std::vector<double> d_list = calculation_order.structure_thicknesses;
    auto [m_r_ps, m_t_ps] = calculate_tr(assemble_e_list_3x3(material_splines, wavelength), d_list, wavelength, theta_0, phi_0);

    std::cout << m_r_ps << std::endl;
    std::cout << m_t_ps << std::endl;

    // double reflection = std::sqrt(std::pow(m_r_ps(0, 0).real(), 2) + std::pow(m_r_ps(1, 1).real(), 2));
    // double transmission = std::sqrt(std::pow(m_t_ps(0, 0).real(), 2) + std::pow(m_t_ps(1, 1).real(), 2));
    return {m_r_ps, m_t_ps};
}
*/

double FilterStack::calculate_reflection_transmission_absorption(const char *type, const char *polarization, double wavelength, double theta_0, double phi_0) //, std::vector<double> d_list = FilterStack::d_list_initial)
{
    std::vector<double> d_list = calculation_order.structure_thicknesses;

    // Add 0 to the beginning - substrate
    d_list.insert(d_list.begin(), 0.0);
    // Add 0 to the end - incident medium
    d_list.push_back(0.0);
    theta_0 = theta_0 * M_PI / 180.0;

    auto [m_r_ps, m_t_ps] = calculate_tr(assemble_e_list_3x3(material_splines, wavelength), d_list, wavelength, theta_0, phi_0);

    double reflectivity_s = m_r_ps.cwiseAbs2()(1, 1);
    double transmissivity_s = m_t_ps.cwiseAbs2()(1, 1);
    double reflectivity_p = m_r_ps.cwiseAbs2()(0, 0);
    double transmissivity_p = m_t_ps.cwiseAbs2()(0, 0);
    double absorption_s = 1 - reflectivity_s - transmissivity_s;
    double absorption_p = 1 - reflectivity_p - transmissivity_p;

    // To implement: azimuthal angle phi_0
    double reflectivity = (reflectivity_s + reflectivity_p) / 2;
    double transmissivity = (transmissivity_s + transmissivity_p) / 2;
    double absorption = 1 - reflectivity - transmissivity;

    // for(const auto& d : d_list) {
    //     std::cout << d << " ";
    // }
    // std::cout << theta_0 << " ";
    // std::cout << wavelength << " ";
    // std::cout << "reflectivity_s: " << reflectivity_s << std::endl;
    // std::cout << "transmissivity_s: " << transmissivity_s << std::endl;
    // std::cout << "reflectivity_p: " << reflectivity_p << std::endl;
    // std::cout << "transmissivity_p: " << transmissivity_p << std::endl;
    // std::cout << "absorption_s: " << absorption_s << std::endl;
    // std::cout << "absorption_p: " << absorption_p << std::endl;
    // std::cout << "reflectivity: " << reflectivity << std::endl;
    // std::cout << "transmissivity: " << transmissivity << std::endl;
    // std::cout << "absorption: " << absorption << std::endl;

    if (strcmp(type, "t") == 0)
    {
        if (strcmp(polarization, "s") == 0)
            return transmissivity_s;
        else if (strcmp(polarization, "p") == 0)
            return transmissivity_p;
        else
            return transmissivity;
    }
    else if (strcmp(type, "r") == 0)
    {
        if (strcmp(polarization, "s") == 0)
            return reflectivity_s;
        else if (strcmp(polarization, "p") == 0)
            return reflectivity_p;
        else
            return reflectivity;
    }
    else if (strcmp(type, "a") == 0)
    {
        if (strcmp(polarization, "s") == 0)
            return absorption_s;
        else if (strcmp(polarization, "p") == 0)
            return absorption_p;
        else
            return absorption;
    }
    else
    {
        throw std::invalid_argument("Invalid type argument. Please use 'r', 't', or 'a'.");
    }
}

/*
void FilterStack::calculate_and_save_ar_reflection()
{
    std::vector<std::complex<double>> e_listx_wvl, e_listy_wvl, e_listz_wvl;
    std::vector<double> d_list = calculation_order.structure_thicknesses;
    d_list.insert(d_list.begin(), 0.0); // Add 0 to the beginning - substrate
    d_list.push_back(0.0); // Add 0 to the end - incident medium
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

            // out medium: air
            Matrix3cd incident_layer_tensor =
                (Matrix3cd(3, 3) << std::complex<double>(1, 0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(1, 0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1, 0))
                    .finished();
            e_list_3x3.push_back(incident_layer_tensor);

            for (const auto &material : calculation_order.structure_materials)
            {
                Matrix3cd next_layer_tensor =
                    (Matrix3cd(3, 3) << std::complex<double>(materialSplines[material][0](wavelength[i]), materialSplines[material][1](wavelength[i])), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
                     std::complex<double>(0.0, 0.0), std::complex<double>(materialSplines[material][2](wavelength[i]), materialSplines[material][3](wavelength[i])), std::complex<double>(0.0, 0.0),
                     std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(materialSplines[material][4](wavelength[i]), materialSplines[material][5](wavelength[i])))
                        .finished();
                e_list_3x3.push_back(next_layer_tensor);
            }

            // substrate: quartz/glass
            Matrix3cd substrate_layer_tensor =
                (Matrix3cd(3, 3) << std::complex<double>(1.5, 0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(1.5, 0), std::complex<double>(0.0, 0.0),
                 std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.5, 0))
                    .finished();
            e_list_3x3.push_back(substrate_layer_tensor);

            auto [m_r_ps, m_t_ps] = calculate_tr(e_list_3x3, d_list, wavelength[i], v_theta[p], phi_0);

            resultDictionary_reflection_p[key] = m_r_ps.cwiseAbs2()(0, 0);
            resultDictionary_reflection_s[key] = m_r_ps.cwiseAbs2()(1, 1);
            resultDictionary_transmission_p[key] = m_t_ps.cwiseAbs2()(0, 0);
            resultDictionary_transmission_s[key] = m_t_ps.cwiseAbs2()(1, 1);

            // std::cout << "reflection matrix: " << m_r_ps << std::endl;
            // std::cout << "transmission matrix: " << m_t_ps << std::endl;
        }
    };

    std::array<std::map<std::pair<double, double>, std::complex<double>>, 4> resultArray = {resultDictionary_reflection_p, resultDictionary_reflection_s, resultDictionary_transmission_p, resultDictionary_transmission_s};

    dump_to_file(resultArray[0], "R_p.txt");
    dump_to_file(resultArray[1], "R_s.txt");
    dump_to_file(resultArray[2], "T_p.txt");
    dump_to_file(resultArray[3], "T_s.txt");
}
*/

// Function to assemble materials from file
std::map<std::string, std::vector<tk::spline>> FilterStack::assemble_materials(std::string full_path)
{
    CalculationInfo calculation_order = loadCalculationInfo(full_path);

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
std::vector<Matrix3cd> FilterStack::assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength)
{
    std::vector<Matrix3cd> e_list_3x3;

    // out medium: air
    Matrix3cd incident_layer_tensor =
        (Matrix3cd(3, 3) << std::complex<double>(1, 0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(1, 0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1, 0))
            .finished();
    e_list_3x3.push_back(incident_layer_tensor);

    for (const auto &material : calculation_order.structure_materials)
    {
        Matrix3cd next_layer_tensor =
            (Matrix3cd(3, 3) << std::complex<double>(material_splines[material][0](wavelength), material_splines[material][1](wavelength)), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][2](wavelength), material_splines[material][3](wavelength)), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][4](wavelength), (std::max(material_splines[material][5](wavelength), 1e-3))))
                .finished();
        e_list_3x3.push_back(next_layer_tensor);
    }

    // substrate: quartz/glass
    Matrix3cd substrate_layer_tensor =
        (Matrix3cd(3, 3) << std::complex<double>(2.25, 0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0))
            .finished();
    e_list_3x3.push_back(substrate_layer_tensor);

    // After the line where e_list_3x3 is defined or updated
    // for (const auto &matrix : e_list_3x3) {
    //    for (int i = 0; i < matrix.rows(); ++i) {
    //        for (int j = 0; j < matrix.cols(); ++j) {
    //            std::cout << matrix(i, j) << " ";
    //        }
    //        std::cout << std::endl;
    //    }
    //    std::cout << std::endl;
    //}

    return e_list_3x3;
}

void FilterStack::change_material_order(std::vector<int> new_material_order)
{
    // Make sure new_material_order has the same size as structure_materials
    assert(new_material_order.size() == calculation_order.structure_materials.size());

    // Create new vectors to hold the reordered materials and their corresponding values
    std::vector<std::string> reordered_material_list;
    std::vector<double> reordered_d_list;

    // Populate the new vector with the materials in the new order
    for (int i : new_material_order)
    {
        reordered_material_list.push_back(material_order_initial[i]);
        reordered_d_list.push_back(d_list_initial[i]);
    }

    // Replace the old vector with the new one
    calculation_order.structure_materials = reordered_material_list;
    calculation_order.structure_thicknesses = reordered_d_list;
    material_order_int = new_material_order;
}

void FilterStack::change_material_thickness(std::vector<double> new_material_thickness)
{
    // Make sure new_material_thickness has the same size as structure_thicknesses
    assert(new_material_thickness.size() == calculation_order.structure_thicknesses.size());

    // If the material order is the same as the initial one, just replace the
    // old vector with the new one otherwise reorder the d_list so that it works
    if (calculation_order.structure_materials != material_order_initial)
    {
        std::vector<double> reordered_d_list;

        // Populate the new vector with the materials in the new order
        for (int i : material_order_int)
        {
            reordered_d_list.push_back(new_material_thickness[i]);
        }
        calculation_order.structure_thicknesses = reordered_d_list;
    }
    else
    {
        // Replace the old vector with the new one
        calculation_order.structure_thicknesses = new_material_thickness;
    }
}