// g++ -shared -o run_filter_stack.so run_filter_stack.cpp -fPIC

#include "FilterStack.h"

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

	double calculate_reflection(FilterStack *filter_stack, double wavelength, double theta_0, double phi_0, double *d_list, size_t size_d_list)
	{
		// Convert double* to std::vector<double>
		std::vector<double> d_list_vector(d_list, d_list + size_d_list);

		double reflection = filter_stack->calculate_reflection(wavelength, theta_0, phi_0, d_list_vector);
		// double reflection = 1.0;
		return reflection;
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