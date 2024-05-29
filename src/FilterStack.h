#include "core.h"
#include "input.h"
#include <cstring>
#include <omp.h>
#include <sstream>

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
    std::vector<double> d_list_in_initial_order;
    std::vector<bool> incoherent_initial;
    std::vector<std::string> material_order_initial;
    std::vector<double> unique_wavelengths_vector;
    std::map<int, std::vector<Matrix3cd>> dict_assembled_e_list_3x3;
    std::map<int, std::vector<Matrix3cd>> dict_optim_assembled_e_list_3x3;
    bool general_materials_in_stack = true;

    // Public methods

    double calculate_reflection_transmission_absorption(const char *type, const char *polarization, double wavelength, double theta_0, double phi_0);
    std::vector<std::vector<std::vector<double>>> calculate_reflection_transmission_absorption_para(const char *type, const char *polarization, std::vector<double> wavelengths, std::vector<double> thetas_0, std::vector<double> phis_0);
    double calculate_merit(std::vector<double> target_value_vector, std::vector<double> target_wavelength_vector, std::vector<double> target_polar_angle_vector, std::vector<double> target_azimuthal_angle_vector, std::vector<double> target_weights_vector, std::vector<char *> target_condition_vector, std::vector<double> target_tolerance_vector, std::vector<char *> target_type_vector, std::vector<char *> target_polarization_vector, std::vector<char *> target_arithmetic);
    bool check_general_materials();

    // void calculate_and_save_ar_reflection();
    void change_material_order(std::vector<int> new_material_order);
    void change_material_thickness(std::vector<double> material_thickness);

    void get_material_order();
    void get_thicknesses();
    void initialise_optimization(std::vector<double> target_wavelengths_vector);

    void setGeneralMaterialsInStack(bool value)
    {
        general_materials_in_stack = value;
    }
    bool getGeneralMaterialsInStack() const
    {
        return general_materials_in_stack;
    }

    void reset_filter();

    // Default constructor
    FilterStack() = default;

    // Constructor
    FilterStack(const char *filename)
    {
        // First convert the filename to an std::string
        full_path = std::filesystem::current_path() / filename;
        std::tie(material_splines, general_materials_in_stack) = assemble_materials(filename);
        setGeneralMaterialsInStack(general_materials_in_stack);
        calculation_order = loadCalculationInfo(full_path);
        d_list_initial = calculation_order.structure_thicknesses;
        d_list_in_initial_order = d_list_initial;
        incoherent_initial = calculation_order.incoherent;
        material_order_initial = calculation_order.structure_materials;
        material_order_int.resize(calculation_order.structure_materials.size());
        for (int i = 0; i < calculation_order.structure_materials.size(); ++i)
        {
            material_order_int[i] = i;
        }
        // Change the wavelength extrema if needed
        initialise_e_list_3x3(material_splines, 250, 1200, 0.1);
    }

    // The destructor should probably be populated more
    ~FilterStack() = default;

private:
    // Private methods
    std::pair<std::map<std::string, std::vector<tk::spline>>, bool> assemble_materials(std::string full_path);
    std::vector<Matrix3cd> assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength);
    void initialise_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength_min, double wavelength_max, double wavelength_step);
    void initialise_e_list_3x3_optim(std::map<std::string, std::vector<tk::spline>> material_splines, std::vector<double> wavelengths);
};

std::vector<std::vector<std::vector<double>>> FilterStack::calculate_reflection_transmission_absorption_para(const char *type, const char *polarization, std::vector<double> wavelengths, std::vector<double> thetas_0, std::vector<double> phis_0)
{
    std::vector<std::vector<std::vector<double>>> result(phis_0.size(), std::vector<std::vector<double>>(thetas_0.size(), std::vector<double>(wavelengths.size())));

    std::vector<double> d_list = calculation_order.structure_thicknesses;
    std::vector<bool> incoherent = calculation_order.incoherent;

    // Add 0 to the beginning and false - exit medium
    d_list.insert(d_list.begin(), 0.0);
    incoherent.insert(incoherent.begin(), false);
    // Add 0 to the end and false - incident medium
    d_list.push_back(0.0);
    incoherent.push_back(false);

#pragma omp parallel for collapse(3)
    for (int p = 0; p < phis_0.size(); p++)
    {
        for (int n = 0; n < thetas_0.size(); n++)
        {
            for (int i = 0; i < wavelengths.size(); i++)
            {
                double reflectivity, transmissivity;

                if (strcmp(polarization, "") == 0)
                {
                    std::tie(reflectivity, transmissivity) = calculate_rt_unpolarized(assemble_e_list_3x3(material_splines, wavelengths[i]), d_list, incoherent, wavelengths[i], thetas_0[n], phis_0[p]);
                }
                else
                {
                    std::tie(reflectivity, transmissivity) = calculate_rt(assemble_e_list_3x3(material_splines, wavelengths[i]), d_list, incoherent, polarization, wavelengths[i], thetas_0[n], phis_0[p]);
                }

                if (strcmp(type, "t") == 0)
                {
                    result[p][n][i] = transmissivity;
                }
                else if (strcmp(type, "r") == 0)
                {
                    result[p][n][i] = reflectivity;
                }
                else if (strcmp(type, "a") == 0)
                {
                    result[p][n][i] = 1 - reflectivity - transmissivity;
                }
                else
                {
                    // #pragma omp critical
                    {
                        std::cout << "Invalid type. Please use 't', 'r' or 'a'." << std::endl;
                    }
                }
            }
        }
    }
    return result;
}

double FilterStack::calculate_merit(std::vector<double> target_value_vector, std::vector<double> target_wavelength_vector, std::vector<double> target_polar_angle_vector, std::vector<double> target_azimuthal_angle_vector, std::vector<double> target_weights_vector, std::vector<char *> target_condition_vector, std::vector<double> target_tolerance_vector, std::vector<char *> target_type_vector, std::vector<char *> target_polarization_vector, std::vector<char *> target_arithmetic)
{
    std::vector<double> target_calculated_values;

    double merit = 0.0;

    // First calculate the actual values for the different targets
#pragma omp parallel for reduction(+ : merit)
    for (size_t i = 0; i < target_value_vector.size(); i++)
    {
        // std::cout << target_type_vector[i] << " " << target_polarization_vector[i] << " " << target_wavelength_vector[i] << " " << target_polar_angle_vector[i] << " " << target_azimuthal_angle_vector[i] << " " << target_condition_vector[i][0] << " ";

        // Calculate actual value
        double target_calculated = calculate_reflection_transmission_absorption(target_type_vector[i], target_polarization_vector[i], target_wavelength_vector[i], target_polar_angle_vector[i], target_azimuthal_angle_vector[i]);
        target_calculated_values.push_back(target_calculated);
        // std::cout << target_calculated << std::endl;0
    }

    // Iterate again over the target values again and potentially do arithmetic
    // operations between the targets to calculate the merit
    for (size_t i = 0; i < target_value_vector.size(); i++)
    {
        double target_calculated = target_calculated_values[i];
        // std::cout << target_calculated_values[0] << std::endl;
        // std::cout << target_calculated_values[1] << std::endl;

        // Perform arithmetic operation if target_arithmetic is not empty
        if (target_arithmetic[i][0] != '\0')
        {
            std::string str(target_arithmetic[i]);
            // std::cout << str << std::endl;
            size_t pos = 0;
            double target_calculated = target_calculated_values[std::stoi(str) - 1];
            // std::cout << target_calculated << std::endl;
            char operation;

            while ((pos = str.find_first_of("+-*/", pos)) != std::string::npos)
            {
                operation = str[pos];
                str = str.substr(pos + 1);
                pos = 0;
                int index = std::stoi(str.substr(0, str.find_first_of("+-*/", pos))) - 1;

                switch (operation)
                {
                case '+':
                    target_calculated += target_calculated_values[index];
                    break;
                case '-':
                    target_calculated -= target_calculated_values[index];
                    break;
                case '*':
                    target_calculated *= target_calculated_values[index];
                    break;
                case '/':
                    target_calculated /= target_calculated_values[index];
                    break;
                }
            }
            // std::cout << target_calculated << std::endl;
        }
        else
        {
            continue;
        }
        // Calculate merit for equal condition
        if (target_condition_vector[i][0] == '=' && target_calculated != target_value_vector[i])
        {
            merit += pow((target_calculated - target_value_vector[i]) / target_tolerance_vector[i], 2) * target_weights_vector[i];
        }

        // Calculate merit for larger condition
        else if (target_condition_vector[i][0] == '>' && target_calculated < target_value_vector[i])
        {
            merit += pow((target_calculated - target_value_vector[i]) / target_tolerance_vector[i], 2) * target_weights_vector[i];
        }

        // Calculate merit for smaller condition
        else if (target_condition_vector[i][0] == '<' && target_calculated > target_value_vector[i])
        {
            merit += pow((target_calculated - target_value_vector[i]) / target_tolerance_vector[i], 2) * target_weights_vector[i];
        }

        // std::cout << merit << std::endl;
    }

    return merit;
}

double FilterStack::calculate_reflection_transmission_absorption(const char *type, const char *polarization, double wavelength, double theta_0, double phi_0)
{
    std::vector<double> d_list = calculation_order.structure_thicknesses;
    std::vector<bool> incoherent = calculation_order.incoherent;

    // Add 0 to the beginning and false - exit medium
    d_list.insert(d_list.begin(), 0.0);
    incoherent.insert(incoherent.begin(), false);
    // Add 0 to the end and false - incident medium
    d_list.push_back(0.0);
    incoherent.push_back(false);
    // for (bool value : incoherent)
    // {
    //     std::cout << value << " ";
    // }
    // std::cout << std::endl;

    // Add 0 to the beginning - substrate
    // d_list.insert(d_list.begin(), 0.0);
    // Add 0 to the end - incident medium
    // d_list.push_back(0.0);

    theta_0 = theta_0 * M_PI / 180.0;

    double wavelength_key = static_cast<int>(wavelength * 10);
    double reflectivity, transmissivity;

    if (strcmp(polarization, "") == 0)
    {
        std::tie(reflectivity, transmissivity) = calculate_rt_unpolarized(dict_optim_assembled_e_list_3x3[wavelength_key], d_list, incoherent, wavelength, theta_0, phi_0);
    }
    else
    {
        std::tie(reflectivity, transmissivity) = calculate_rt(dict_optim_assembled_e_list_3x3[wavelength_key], d_list, incoherent, polarization, wavelength, theta_0, phi_0);
    }

    if (strcmp(type, "t") == 0)
    {
        return transmissivity;
    }
    else if (strcmp(type, "r") == 0)
    {
        return reflectivity;
    }
    else if (strcmp(type, "a") == 0)
    {
        return 1 - reflectivity - transmissivity;
    }
    else
    {
        std::cout << "Invalid type. Please use 't', 'r' or 'a'." << std::endl;
        return 0;
    }
}

// Function to assemble materials from file
std::pair<std::map<std::string, std::vector<tk::spline>>, bool> FilterStack::assemble_materials(std::string full_path)
{
    CalculationInfo calculation_order = loadCalculationInfo(full_path);

    std::vector<std::string> list_materials = calculation_order.structure_materials;

    // Append incident medium

    list_materials.push_back(calculation_order.incidentMediumMaterial);

    // Append exit medium

    list_materials.push_back(calculation_order.exitMediumMaterial);

    std::vector<std::string> unique_materials = getUniqueMembers(list_materials);

    CSVParser parser;
    std::map<std::string, std::vector<tk::spline>> material_splines;
    bool contains_general_material = false;

    for (const auto &material : unique_materials)
    {
        std::vector<tk::spline> splines;
        bool is_general_material;
        std::tie(splines, is_general_material) = parser.importIndexFromFile(material);
        if (is_general_material == true)
        {
            contains_general_material = true;
        }
        material_splines[material] = splines;
    }

    return std::make_pair(material_splines, contains_general_material);
}

void FilterStack::initialise_optimization(std::vector<double> target_wavelengths_vector)
{
    std::set<double> unique_wavelengths(target_wavelengths_vector.begin(), target_wavelengths_vector.end());
    std::vector<double> unique_wavelengths_vector_converted(unique_wavelengths.begin(), unique_wavelengths.end());
    unique_wavelengths_vector = unique_wavelengths_vector_converted;
    initialise_e_list_3x3_optim(material_splines, unique_wavelengths_vector);
}

// Function to assemble e_list_3x3 for a given wavelength using the material_splines obtained from file
void FilterStack::initialise_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength_min, double wavelength_max, double wavelength_step)
{
    for (double wavelength = wavelength_min; wavelength <= wavelength_max; wavelength += wavelength_step)
    {
        int wavelength_key = static_cast<int>(wavelength * 10);
        dict_assembled_e_list_3x3[wavelength_key] = assemble_e_list_3x3(material_splines, wavelength);
    }
}

void FilterStack::initialise_e_list_3x3_optim(std::map<std::string, std::vector<tk::spline>> material_splines, std::vector<double> wavelengths)
{
    for (double wavelength : wavelengths)
    {
        int wavelength_key = static_cast<int>(wavelength * 10);
        dict_optim_assembled_e_list_3x3[wavelength_key] = assemble_e_list_3x3(material_splines, wavelength);
    }
}

std::vector<Matrix3cd> FilterStack::assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength)
{
    std::vector<Matrix3cd> e_list_3x3;

    // exit medium
    Matrix3cd exit_layer_tensor =
        (Matrix3cd(3, 3) << std::complex<double>(material_splines[calculation_order.exitMediumMaterial][0](wavelength), std::max(material_splines[calculation_order.exitMediumMaterial][1](wavelength), 0.0)), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[calculation_order.exitMediumMaterial][2](wavelength), std::max(material_splines[calculation_order.exitMediumMaterial][3](wavelength), 0.0)), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[calculation_order.exitMediumMaterial][4](wavelength), std::max(material_splines[calculation_order.exitMediumMaterial][5](wavelength), 0.0)))
            .finished();
    e_list_3x3.push_back(exit_layer_tensor.conjugate());

    for (const auto &material : calculation_order.structure_materials)
    {
        Matrix3cd next_layer_tensor =
            (Matrix3cd(3, 3) << std::complex<double>(material_splines[material][0](wavelength), material_splines[material][1](wavelength)), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][2](wavelength), material_splines[material][3](wavelength)), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][4](wavelength), (std::max(material_splines[material][5](wavelength), 0.0))))
                .finished();
        e_list_3x3.push_back(next_layer_tensor.conjugate());
    }

    // incident medium
    Matrix3cd incident_layer_tensor =
        (Matrix3cd(3, 3) << std::complex<double>(material_splines[calculation_order.incidentMediumMaterial][0](wavelength), std::max(material_splines[calculation_order.incidentMediumMaterial][1](wavelength), 0.0)), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[calculation_order.incidentMediumMaterial][2](wavelength), std::max(material_splines[calculation_order.incidentMediumMaterial][3](wavelength), 0.0)), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[calculation_order.incidentMediumMaterial][4](wavelength), std::max(material_splines[calculation_order.incidentMediumMaterial][5](wavelength), 0.0)))
            .finished();
    e_list_3x3.push_back(incident_layer_tensor.conjugate());

    return e_list_3x3;
}

void FilterStack::change_material_order(std::vector<int> new_material_order)
{
    size_t size_arrays = new_material_order.size();

    // Make sure new_material_order has the same size as structure_materials
    assert(size_arrays == calculation_order.structure_materials.size());

    // Create new vectors to hold the reordered materials and their corresponding values
    std::vector<std::string> reordered_material_list(size_arrays);
    std::vector<double> reordered_d_list(size_arrays);
    std::map<int, std::vector<Matrix3cd>> reordered_dict_optim_assembled_e_list_3x3;
    std::vector<bool> reordered_incoherent(size_arrays);

    bool first_loop = true;
    std::vector<Matrix3cd> reordered_e_list_3x3(size_arrays + 2);
    size_t i;
    int wavelength_key;

    for (double wavelength : unique_wavelengths_vector)
    {

        wavelength_key = static_cast<int>(wavelength * 10);

        reordered_e_list_3x3[0] = dict_assembled_e_list_3x3[wavelength_key][0];

        i = 0;

        while (i < size_arrays)
        {

            if (first_loop)
            {

                reordered_material_list[i] = material_order_initial[new_material_order[i]];

                reordered_d_list[i] = d_list_in_initial_order[new_material_order[i]];
                reordered_incoherent[i] = incoherent_initial[new_material_order[i]];
            }

            reordered_e_list_3x3[i + 1] = dict_assembled_e_list_3x3[wavelength_key][new_material_order[i] + 1];

            i++;
        }

        first_loop = false;

        reordered_e_list_3x3[size_arrays + 1] = dict_assembled_e_list_3x3[wavelength_key][size_arrays + 1];

        reordered_dict_optim_assembled_e_list_3x3[wavelength_key] = reordered_e_list_3x3;
    }

    // Replace the old vector with the new one
    calculation_order.structure_materials = reordered_material_list;
    calculation_order.structure_thicknesses = reordered_d_list;
    calculation_order.incoherent = reordered_incoherent;
    dict_optim_assembled_e_list_3x3 = reordered_dict_optim_assembled_e_list_3x3;
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

    // Set the d_list_in_initial_order to the new d_list
    d_list_in_initial_order = new_material_thickness;
}

void FilterStack::reset_filter()
{
    // Function that resets the filter parameters to their initial state (as
    // read from the file)
    calculation_order.structure_materials = material_order_initial;
    calculation_order.structure_thicknesses = d_list_initial;
    calculation_order.incoherent = incoherent_initial;
    material_order_int.resize(calculation_order.structure_materials.size());

    for (int i = 0; i < calculation_order.structure_materials.size(); ++i)
    {
        material_order_int[i] = i;
    }
}

void FilterStack::get_material_order()
{
    for (int i = 0; i < calculation_order.structure_materials.size(); i++)
    {
        std::cout << calculation_order.structure_materials[i] << " ";
    }
}

void FilterStack::get_thicknesses()
{
    for (int i = 0; i < calculation_order.structure_thicknesses.size(); i++)
    {
        std::cout << calculation_order.structure_thicknesses[i] << " ";
    }
}