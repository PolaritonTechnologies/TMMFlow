// g++ -shared -o run.so run.cpp -fPIC
#include "core.h"

#include <vector>
#include <complex>
#include <cmath>
#include "Eigen/Dense"

extern "C"
{
    double run_function(double wavelength, double theta_0, double phi_0)
    {

        std::vector<Matrix3cd> e_list_3x3 = {
            (Matrix3cd(3, 3) << std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0))
                .finished(),

            (Matrix3cd(3, 3) << std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931))
                .finished(),

            (Matrix3cd(3, 3) << std::complex<double>(2.90627602, 0.84588284), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(4.3259341, 4.83066447), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(2.98285307, 1.06542853))
                .finished(),

            (Matrix3cd(3, 3) << std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931))
                .finished(),

            (Matrix3cd(3, 3) << std::complex<double>(2.25, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0.0), std::complex<double>(0.0, 0.0),
             std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0.0))
                .finished()};

        std::vector<double> d_list = {0, 20, 105, 100, 0};

        auto [m_r_ps, m_t_ps] = calculate_tr(e_list_3x3, d_list, wavelength, theta_0, phi_0);
        std::cout << "reflection matrix: " << m_r_ps << std::endl;
        std::cout << "transmission matrix: " << m_t_ps << std::endl;

        double reflection = std::sqrt(std::pow(m_r_ps(0, 0).real(), 2) + std::pow(m_r_ps(1, 1).real(), 2));
        double transmission = std::sqrt(std::pow(m_t_ps(0, 0).real(), 2) + std::pow(m_t_ps(1, 1).real(), 2));

        return reflection;
        // return {m_r_ps, m_t_ps};
    }
}