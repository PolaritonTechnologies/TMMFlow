// g++ -shared -o run_filter_stack.so run_filter_stack.cpp -fPIC

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

	double calculate_reflection_transmission_absorption(FilterStack *filter_stack, const char *type, const char *polarization, double wavelength, double theta_0, double phi_0)
	{

		double values = filter_stack->calculate_reflection_transmission_absorption(type, polarization, wavelength, theta_0, phi_0);

		while (std::isnan(values) || values > 1.0)
		{
			theta_0 = theta_0 + 0.0001;
			// std::cout << "Increasing theta_0 by 0.0001 to avoid kz crash. New theta_0: " << theta_0 << std::endl;
			values = filter_stack->calculate_reflection_transmission_absorption(type, polarization, wavelength, theta_0, phi_0);
		}

		return values;
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

// int main()
// {
//     FilterStack filter_stack = FilterStack("calculation_order.json");

//     auto [m_r_ps, m_t_ps] = filter_stack.calculate_transmission_reflection_matrices(400, 0.001, 0);

//     std::cout << "reflection matrix: " << m_r_ps << std::endl;
//     std::cout << "transmission matrix: " << m_t_ps << std::endl;

//     filter_stack.calculate_and_save_ar_reflection();

//     return 0;
// }