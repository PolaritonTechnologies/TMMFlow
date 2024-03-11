/* Parameters
----------
'wl'= vacuum incident wavelength in nm
'theta_0,phi_0'= polar and azimuth angles as defined in Mansuripur
                    JAP 67(10)
'e_list_3x3'= [n_layer+2,3,3] numpy array: it contains n_layers+2 3x3
                dielectric tensors:
e_list_3x3[0]= 3x3 incident medium dielectric tensor: must be real,diagonal
                and isotropic,
e_list_3x3[n_layers+1] = 3x3 substrate dielectric tensor: must be real,
                            diagonal and isotropic,
e_list_3x3[n]=3x3 dielectric tensor of the n_th layers: arbitrary
'd_list'= n_layers+2 numpy array: contains layer thinknesses:
d_list[0]=d_list[n_layers+1]=0: for the incident medium and substrate
d_list[n]=d_n n_th layer thickness in nm

Returns
-------
'a dictionary'= {'m_r_ps':m_r_ps,  #reflection matrix
                    'm_t_ps':m_t_ps,  #transmission matrix
                    'm_Kn':m_Kn, # wavevectors (n_layers,n_k,n_xyz)
                    'm_En':m_En, # electric field amplitudes (n_layers,n_k,n_xyz,n_pol)
                    'm_Hn':m_Hn, # magnetic field amplitudes (n_layers,n_k,n_xyz,n_pol)
                    'wl': wl,'theta_0': theta_0,'phi_0': phi_0,  #inputs
                    'e_list_3x3': e_list_3x3,'d_list': d_list}   #inputs
""" */

#include <iostream>
#include <map>
#include <vector>
#include <complex>
#include "Eigen/Dense"
#include <unordered_map>
#include <variant>
// #include "nullspace.cpp"
#include "kz_eigenvectors.cpp"
#include "kz_eigenvalues.cpp"
#include "m_abc.cpp"
// #include "json.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Eigen;

std::pair<Matrix2cd, Matrix2cd> calculate_tr_per_layer(
    Matrix2cd m_a12,
    Matrix2cd m_a34,
    Matrix2cd m_b12,
    Matrix2cd m_b34,
    Matrix2cd m_a12_np1,
    Matrix2cd m_a34_np1,
    Matrix2cd m_b12_np1,
    Matrix2cd m_b34_np1,
    Matrix2cd m_c12_np1,
    Matrix2cd m_c34_np1,
    Matrix2cd m_R_np1,
    Matrix2cd m_T)
{
    /*
    This function is supposed to be used recursively to calculate the reflection
    and transmission matrices for each layer. The iteration over each layer
    happens in another function, however.
    */
    // Define matrix R and T for this iteration
    Matrix2cd m_R = Matrix2cd::Zero();
    Matrix2cd m_Tn = Matrix2cd::Zero();

    Matrix2cd f1 = m_b12_np1 * m_c12_np1 + (m_b34_np1 * m_c34_np1) * m_R_np1;
    Matrix2cd f2_inv = m_a12_np1 * m_c12_np1 + (m_a34_np1 * m_c34_np1) * m_R_np1;
    Matrix2cd f2 = f2_inv.inverse();

    // There is a sign issue in f1 for the bottom left entry
    std::cout << m_b12_np1 << std::endl;
    std::cout << m_c12_np1 << std::endl;
    std::cout << m_b34_np1 << std::endl;
    std::cout << m_c34_np1 << std::endl;
    std::cout << m_R_np1 << std::endl;
    // std::cout << f2 << std::endl;

    Matrix2cd f_np1 = f1 * f2;
    Matrix2cd r1_inv = f_np1 * m_a34 - m_b34;
    Matrix2cd r1 = r1_inv.inverse();
    Matrix2cd r2 = m_b12 - f_np1 * m_a12;

    // std::cout << f_np1 << std::endl;
    // std::cout << m_a34 << std::endl;
    // std::cout << m_b34 << std::endl;
    // std::cout << r1_inv << std::endl;

    // std::cout << r1 << std::endl;
    // std::cout << r2 << std::endl;

    // Calculate the reflection for this iteration
    m_R = r1 * r2;

    // Repeat the process for the transmission
    Matrix2cd f1_trans_inv = m_a12_np1 * m_c12_np1 + (m_a34_np1 * m_c34_np1) * m_R_np1;
    Matrix2cd f1_trans = f1_trans_inv.inverse();
    Matrix2cd f2_trans = m_a12 + m_a34 * m_R;

    // Calculate overall transmission and
    m_Tn = f1_trans * f2_trans;
    m_T = m_T * m_Tn;

    return {m_R, m_T};
}

// Function to replace the numpy real_if_close
double real_if_close(std::complex<double> c, double tol = 1e-7)
{
    if (std::abs(c.imag()) < tol)
    {
        return c.real();
    }
    else
    {
        throw std::runtime_error("Imaginary part is not close to zero");
    }
};

std::pair<Matrix2cd, Matrix2cd> calculate_tr(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0)
// void calculate_tr(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0)
{
    // Assuming e_list_3x3 is a std::vector<Eigen::Matrix3cd>
    std::complex<double> n_0 = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_s = std::sqrt(e_list_3x3[e_list_3x3.size() - 1](0, 0));

    // wavevector modulus and in plane components
    std::complex<double> k0 = 2.0 * M_PI / wavelength;
    std::complex<double> kx = -k0 * n_0.real() * sin(theta_0) * cos(phi_0);
    std::complex<double> ky = -k0 * n_0.real() * sin(theta_0) * sin(phi_0);

    // m_a_np1, m_b_np1, m_a12_np1, m_a34_np1, m_b12_np1, m_b34_np1, m_c12_np1, m_c34_np1, m_K_np1 = chunk_1(wavelength, theta_0, phi_0, e_list_3x3[-1], d_list[-1], n_0, n_s)
    Vector4cd v_kz1 = kz_eigenvalues(k0, kx, ky, e_list_3x3.back());
    auto [v_e, v_kz] = kz_eigenvectors(k0, kx, ky, v_kz1, e_list_3x3.back());

    // MatrixXcd v_e = result.first;
    // tor4cd v_kz = result.second;

    auto [m_a_np1, m_b_np1, m_a12_np1, m_a34_np1, m_b12_np1, m_b34_np1, m_c12_np1, m_c34_np1] = m_abc(k0, kx, ky, v_kz, v_e, d_list.back());

    Matrix2cd m_T = Matrix2cd::Identity();
    Matrix2cd m_R_np1 = Matrix2cd::Zero();
    Matrix2cd m_R_0 = Matrix2cd::Zero();

    // Worked up to here
    // Now iterate over all layers starting from the second last going backwards
    for (int i = d_list.size() - 2; i >= 0; --i)
    {
        std::cout << "i: " << i << std::endl;

        v_kz1 = kz_eigenvalues(k0, kx, ky, e_list_3x3[i]);

        // std::cout << "v_kz1: " << v_kz1 << std::endl;
        // std::cout << "e_list_3x3: " << e_list_3x3[i] << std::endl;

        auto [v_e, v_kz] = kz_eigenvectors(k0, kx, ky, v_kz1, e_list_3x3[i]);

        // v_e = result.first;
        // v_kz = result.second;

        auto [m_a, m_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34] = m_abc(k0, kx, ky, v_kz, v_e, d_list[i]);

        // std::cout << v_a << std::endl;
        // std::cout << v_b << std::endl;
        // std::cout << m_a12 << std::endl;
        // std::cout << m_a34 << std::endl;
        // std::cout << m_b12 << std::endl;
        // std::cout << m_b34 << std::endl;
        // std::cout << m_c12 << std::endl;
        // std::cout << m_c34 << std::endl;

        std::pair<Matrix2cd, Matrix2cd> result3 = calculate_tr_per_layer(m_a12, m_a34, m_b12, m_b34, m_a12_np1, m_a34_np1, m_b12_np1, m_b34_np1, m_c12_np1, m_c34_np1, m_R_np1, m_T);

        Matrix2cd m_R = result3.first;
        m_T = result3.second;

        if (i == 0)
        {
            m_R_0 = m_R;
        }
        // std::cout << "m_R" << m_R << std::endl;
        // std::cout << "m_T " << m_T << std::endl;

        // In the next iteration m_a12 --> m_a12_np1, similarly m_R --> m_R_np1.
        m_a12_np1 = m_a12;
        m_a34_np1 = m_a34;
        m_b12_np1 = m_b12;
        m_b34_np1 = m_b34;
        m_c12_np1 = m_c12;
        m_c34_np1 = m_c34;
        m_R_np1 = m_R;
    }

    // This has to be calculated outside the loop
    // rotating m_R to the s,p states
    Matrix2cd p_inc = Matrix2cd::Zero();
    p_inc(0, 0) = cos(theta_0) * cos(phi_0);
    p_inc(0, 1) = -sin(phi_0);
    p_inc(1, 0) = cos(theta_0) * sin(phi_0);
    p_inc(1, 1) = cos(phi_0);
    Matrix2cd p_inc_inv = p_inc.inverse();

    // Finally the  R matrix output...
    Matrix2cd m_r_ps = p_inc_inv * m_R_0 * p_inc;

    // rotating m_T to the s,p states
    double theta_s = asin(real_if_close(sin(theta_0) * n_0.real() / n_s.real()));
    Matrix2cd p_sub = Matrix2cd::Zero();
    p_sub(0, 0) = cos(theta_s) * cos(phi_0);
    p_sub(0, 1) = -sin(phi_0);
    p_sub(1, 0) = cos(theta_s) * sin(phi_0);
    p_sub(1, 1) = cos(phi_0);
    Matrix2cd p_sub_inv = p_sub.inverse();

    // Finally the  T matrix output...
    Matrix2cd m_t_ps = p_sub_inv * m_T * p_inc;

    return {m_r_ps, m_t_ps};
};

int main()
{
    /*  Example input  */
    double wavelength = 400;
    double theta_0 = 0.0017453292519943296;
    double phi_0 = 0;

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

    // calculate_tr(e_list_3x3, d_list, wavelength, theta_0, phi_0);
    auto [m_r_ps, m_t_ps] = calculate_tr(e_list_3x3, d_list, wavelength, theta_0, phi_0);
    std::cout << "reflection matrix: " << m_r_ps << std::endl;
    std::cout << "transmission matrix: " << m_t_ps << std::endl;

    return 0;
}
