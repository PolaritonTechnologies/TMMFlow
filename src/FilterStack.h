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
    bool general_materials_in_stack = true;

    // Public methods
    // std::pair<Matrix2cd, Matrix2cd> calculate_transmission_reflection_matrices(double wavelength, double theta_0, double phi_0);
    double calculate_reflection_transmission_absorption(const char *type, const char *polarization, double wavelength, double theta_0, double phi_0, bool is_general_case); //, std::vector<double> d_list);
    bool check_general_materials();

    // void calculate_and_save_ar_reflection();
    void change_material_order(std::vector<int> new_material_order);
    void change_material_thickness(std::vector<double> material_thickness);

    void get_material_order();
    void get_thicknesses();

    void setGeneralMaterialsInStack(bool value) {
        general_materials_in_stack = value;
    }
    bool getGeneralMaterialsInStack() const {
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
    std::pair<std::map<std::string, std::vector<tk::spline>>, bool> assemble_materials(std::string full_path);
    std::vector<Matrix3cd> assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength);
};

double FilterStack::calculate_reflection_transmission_absorption(const char *type, const char *polarization, double wavelength, double theta_0, double phi_0, bool is_general_case) //, std::vector<double> d_list = FilterStack::d_list_initial)
{
    std::vector<double> d_list = calculation_order.structure_thicknesses;

    // Add 0 to the beginning - substrate
    d_list.insert(d_list.begin(), 0.0);
    // Add 0 to the end - incident medium
    d_list.push_back(0.0);
    theta_0 = theta_0 * M_PI / 180.0;

    double reflectivity_s, transmissivity_s, reflectivity_p, transmissivity_p;

    std::tie(reflectivity_s, reflectivity_p, transmissivity_s, transmissivity_p) = calculate_tr(assemble_e_list_3x3(material_splines, wavelength), d_list, wavelength, theta_0, phi_0, is_general_case);

    double absorption_s = 1 - reflectivity_s - transmissivity_s;
    double absorption_p = 1 - reflectivity_p - transmissivity_p;

    // To implement: azimuthal angle phi_0
    double reflectivity = (reflectivity_s + reflectivity_p) / 2;
    double transmissivity = (transmissivity_s + transmissivity_p) / 2;
    double absorption = 1 - reflectivity - transmissivity;

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

// Function to assemble materials from file
std::pair<std::map<std::string, std::vector<tk::spline>>, bool> FilterStack::assemble_materials(std::string full_path)
{
    CalculationInfo calculation_order = loadCalculationInfo(full_path);

    std::vector<std::string> unique_materials = getUniqueMembers(calculation_order.structure_materials);
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

// Function to assemble e_list_3x3 for a given wavelength using the material_splines obtained from file
std::vector<Matrix3cd> FilterStack::assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength)
{
    std::vector<Matrix3cd> e_list_3x3;

    // substrate: quartz/glass
    Matrix3cd substrate_layer_tensor =
        (Matrix3cd(3, 3) << std::complex<double>(2.25, 0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0))
            .finished();
    e_list_3x3.push_back(substrate_layer_tensor);

    // std::reverse(structure_materials.begin(), structure_materials.end());

    for (const auto &material : calculation_order.structure_materials)
    {
        Matrix3cd next_layer_tensor =
            (Matrix3cd(3, 3) << std::complex<double>(material_splines[material][0](wavelength), material_splines[material][1](wavelength)), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][2](wavelength), material_splines[material][3](wavelength)), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(material_splines[material][4](wavelength), (std::max(material_splines[material][5](wavelength), 1e-3))))
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
        reordered_d_list.push_back(calculation_order.structure_thicknesses[material_order_int[i]]);
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

void FilterStack::reset_filter()
{
    // Function that resets the filter parameters to their initial state (as
    // read from the file)
    calculation_order.structure_materials = material_order_initial;
    calculation_order.structure_thicknesses = d_list_initial;
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