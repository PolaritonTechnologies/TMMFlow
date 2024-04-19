// g++ -shared -o run_filter_stack.so run_filter_stack.cpp -fPIC -fopenmp

#include "FilterStack.h"
#include <iostream>
#include <cmath>

extern "C"
{

	FilterStack *createFilterStack(const char *filename)
	{
		return new FilterStack(filename);
	}

	void destroyFilterStack(FilterStack *filter_stack)
	{
		delete filter_stack;
	}

	bool getGeneralMaterialsInStack(FilterStack *filter_stack)
	{
		return filter_stack->getGeneralMaterialsInStack();
	}

	double calculate_reflection_transmission_absorption(FilterStack *filter_stack, const char *type, const char *polarization, double wavelength, double theta_0, double phi_0, bool is_general_case)
	{

		double values = filter_stack->calculate_reflection_transmission_absorption(type, polarization, wavelength, theta_0, phi_0, is_general_case);

		while (std::isnan(values) || values > 1.0)
		{
			theta_0 = theta_0 + 0.0001;
			// std::cout << "Increasing theta_0 by 0.0001 to avoid kz crash. New theta_0: " << theta_0 << std::endl;
			values = filter_stack->calculate_reflection_transmission_absorption(type, polarization, wavelength, theta_0, phi_0, is_general_case);
		}

		return values;
	}

	char *calculate_reflection_transmission_absorption_para(FilterStack *filter_stack, const char *type, const char *polarization, double *wavelengths, size_t wavelengths_size, double *theta_0, size_t thetas_0_size, double *phis_0, size_t phis_0_size, bool is_general_case)
	{
		std::vector<double> wavelengths_vector(wavelengths, wavelengths + wavelengths_size);
		std::vector<double> thetas_0_vector(theta_0, theta_0 + thetas_0_size);
		std::vector<double> phis_0_vector(phis_0, phis_0 + phis_0_size);

		// Convert each value in thetas_0_vector from degrees to radians
		std::transform(thetas_0_vector.begin(), thetas_0_vector.end(), thetas_0_vector.begin(),
					   [](double theta)
					   { return theta * M_PI / 180.0; });
		// Convert each value in phis_0_vector from degrees to radians
		std::transform(phis_0_vector.begin(), phis_0_vector.end(), phis_0_vector.begin(),
					   [](double phi)
					   { return phi * M_PI / 180.0; });

		std::vector<std::vector<std::vector<double>>> result_array = filter_stack->calculate_reflection_transmission_absorption_para(type, polarization, wavelengths_vector, thetas_0_vector, phis_0_vector, is_general_case);
		if (!is_general_case)
		{
			for (size_t p = 0; p < result_array.size(); p++)
			{
				for (size_t n = 0; n < result_array.size(); n++)
				{
					for (size_t i = 0; i < result_array.size(); i++)
					{
						while (std::isnan(result_array[p][n][i]) || result_array[p][n][i] > 1.0)
						{
							double theta_0;
							theta_0 = thetas_0_vector[n] + 0.0001;
							// std::cout << "Increasing theta_0 by 0.0001 to avoid kz crash. New theta_0: " << theta_0 << std::endl;
							result_array[p][n][i] = filter_stack->calculate_reflection_transmission_absorption(type, polarization, wavelengths_vector[i], theta_0, phis_0_vector[p], is_general_case);
						}
					}
				}
			}
		}

		std::stringstream result_string;
		for (size_t p = 0; p < phis_0_size; p++)
		{
			if (p != 0)
			{
				result_string << "=";
			}

			for (size_t n = 0; n < thetas_0_size; n++)
			{
				if (n != 0)
				{
					result_string << "+";
				}

				for (size_t i = 0; i < wavelengths_size; i++)
				{
					if (i != 0)
					{
						result_string << "--";
					}

					result_string << result_array[p][n][i];
				}
			}
		}

		return std::strcpy(new char[result_string.str().length() + 1], result_string.str().c_str());
	}

	double calculate_merit(FilterStack *filter_stack, char **target_type, char **target_polarization, double *target_value, double *target_wavelength, double *target_polar_angle, double *target_azimuthal_angle, double *target_weights, char **target_condition, double *target_tolerance, size_t target_size, const char *core_selection)
	{
		std::vector<double> target_value_vector(target_value, target_value + target_size);
		std::vector<double> target_wavelength_vector(target_wavelength, target_wavelength + target_size);
		std::vector<double> target_polar_angle_vector(target_polar_angle, target_polar_angle + target_size);
		std::vector<double> target_azimuthal_angle_vector(target_azimuthal_angle, target_azimuthal_angle + target_size);
		std::vector<double> target_weights_vector(target_weights, target_weights + target_size);
		std::vector<double> target_tolerance_vector(target_tolerance, target_tolerance + target_size);
		std::vector<char *> target_condition_vector(target_condition, target_condition + target_size);
		std::vector<char *> target_type_vector(target_type, target_type + target_size);
		std::vector<char *> target_polarization_vector(target_polarization, target_polarization + target_size);

		double merit = filter_stack->calculate_merit(target_value_vector, target_wavelength_vector, target_polar_angle_vector, target_azimuthal_angle_vector, target_weights_vector, target_condition_vector, target_tolerance_vector, target_type_vector, target_polarization_vector, core_selection);
		return merit;
	}

	void change_material_thickness(FilterStack *filter_stack, double *d_list, size_t size_d_list)
	{
		std::vector<double> d_list_vector(d_list, d_list + size_d_list);
		filter_stack->change_material_thickness(d_list_vector);
	}

	void change_material_order(FilterStack *filter_stack, int *material_order, size_t size_material_order)
	{
		std::vector<int> material_order_vector(material_order, material_order + size_material_order);
		filter_stack->change_material_order(material_order_vector);
	}

	void reset_filter(FilterStack *filter_stack)
	{
		filter_stack->reset_filter();
	}

	void get_material_order(FilterStack *filter_stack)
	{
		filter_stack->get_material_order();
	}

	void get_thicknesses(FilterStack *filter_stack)
	{
		filter_stack->get_thicknesses();
	}
};