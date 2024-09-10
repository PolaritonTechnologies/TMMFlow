/*
 * Header to represent a filter stack defined by a an input file
 */

#pragma once
#include <filesystem>
#include <map>
#include <vector>
#include <string>
#include "core.h"
#include "input.h"

class FilterStack{ 
    
public:

    // Constructor
    FilterStack(const char *json_text);

    // Public variables
    std::filesystem::path full_path;
    std::map<std::string, std::vector<tk::spline>> material_splines;
    CalculationInfo calculation_order;
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
    double calculate_reflection_transmission_absorption(const char *type, double s_polarization_percentage, double wavelength, double theta_0, double phi_0);
    std::vector<std::vector<std::vector<double>>> calculate_reflection_transmission_absorption_para(const char *type, const double s_polarization_percentage, std::vector<double> wavelengths, std::vector<double> thetas_0, std::vector<double> phis_0);
    double calculate_merit(std::vector<double> target_value_vector, std::vector<double> target_wavelength_vector, std::vector<double> target_polar_angle_vector, std::vector<double> target_azimuthal_angle_vector, std::vector<double> target_weights_vector, std::vector<char *> target_condition_vector, std::vector<double> target_tolerance_vector, std::vector<char *> target_type_vector, std::vector<double> target_polarization_vector, std::vector<char *> target_arithmetic);
    bool check_general_materials();
    void change_material_order(std::vector<int> new_material_order);
    void change_material_thickness(std::vector<double> material_thickness);
    void get_material_order();
    void get_thicknesses();
    void initialise_optimization(std::vector<double> target_wavelengths_vector);
    void setGeneralMaterialsInStack(bool value);
    bool getGeneralMaterialsInStack() const;
    void reset_filter();
    
private:
    // Private methods
    std::pair<std::map<std::string, std::vector<tk::spline>>, bool> assemble_materials(const std::string &json_text);
    std::vector<Matrix3cd> assemble_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength);
    void initialise_e_list_3x3(std::map<std::string, std::vector<tk::spline>> material_splines, double wavelength_min, double wavelength_max, double wavelength_step);
    void initialise_e_list_3x3_optim(std::map<std::string, std::vector<tk::spline>> material_splines, std::vector<double> wavelengths);

};