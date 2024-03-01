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
#include "nullspace.cpp"
#include "kz_eigenvectors.cpp"
#include "kz_eigenvalues.cpp"
#include "m_abc.cpp"

using namespace Eigen;

std::map<std::string, int> rt(float wl, float theta_0, float phi_0, const std::vector<Eigen::MatrixXcf> &e_list_3x3, const std::vector<float> &d_list)
{
    std::map<std::string, int> result_dictionary;

    // Can implement incident medium and input checks here

    // Calculate n_0
    std::complex<float> n_0 = std::sqrt(e_list_3x3.front()(0, 0));
    std::cout << n_0 << std::endl;

    // Calculate n_s
    std::complex<float> n_s = std::sqrt(e_list_3x3.back()(0, 0));
    std::cout << n_s << std::endl;

    // Calculate wavevector modulus and in plane components
    std::complex<float> k0 = 2.0 * M_PI / wl;
    std::complex<float> kx = -k0 * n_0 * (sin(theta_0) * cos(phi_0));
    std::complex<float> ky = -k0 * n_0 * (sin(theta_0) * sin(phi_0));

    std::cout << k0 << " " << kx << " " << ky << std::endl;

    // kz,v_e and boundary and propagation matrix for R and T
    std::vector<Eigen::MatrixXcd> m_a12(e_list_3x3.size(), Eigen::MatrixXcd::Zero(2, 2));
    std::vector<Eigen::MatrixXcd> m_a34 = m_a12;
    std::vector<Eigen::MatrixXcd> m_b12 = m_a12;
    std::vector<Eigen::MatrixXcd> m_b34 = m_a12;
    std::vector<Eigen::MatrixXcd> m_c12 = m_a12;
    std::vector<Eigen::MatrixXcd> m_c34 = m_a12;
    std::vector<Eigen::VectorXcd> m_a(e_list_3x3.size(), Eigen::VectorXcd::Zero(4));
    std::vector<Eigen::VectorXcd> m_b = m_a;
    std::vector<Eigen::MatrixXcd> m_Kn(e_list_3x3.size(), Eigen::MatrixXcd::Zero(4, 3));

    for (size_t n = 0; n < e_list_3x3.size(); ++n)
    {
        Eigen::VectorXcd v_kz = kz_eigenvalues(k0, kx, ky, e_list_3x3[n]);
        Eigen::VectorXcd v_e;
        std::tie(v_e, v_kz) = kz_eigenvectors(k0, kx, ky, v_kz, e_list_3x3[n]);

        // Storing the wavevectors
        m_Kn[n].col(0).setConstant(kx); // kx
        m_Kn[n].col(1).setConstant(ky); // ky
        m_Kn[n].col(2) = v_kz;          // kz

        std::tie(m_a[n], m_b[n], m_a12[n], m_a34[n], m_b12[n], m_b34[n], m_c12[n], m_c34[n]) = m_abc(k0, kx, ky, v_kz, v_e, d_list[n]);
    }

    std::vector<Matrix2cd> m_R(m_c12.size());

    for (int n = e_list_3x3.size() - 3; n >= 0; --n)
    {
        // building the first factor for the F_n+1 matrix
        Matrix2cd f1 = m_b12[n + 1] * m_c12[n + 1] + m_b34[n + 1] * m_c34[n + 1] * m_R[n + 1];

        // building the second factor for the F_n+1 matrix
        Matrix2cd f2_inv = m_a12[n + 1] * m_c12[n + 1] + m_a34[n + 1] * m_c34[n + 1] * m_R[n + 1];
        Matrix2cd f2 = f2_inv.inverse();

        // F_n+1 matrix
        Matrix2cd f_np1 = f1 * f2;

        // R_n
        Matrix2cd r1_inv = f_np1 * m_a34[n] - m_b34[n];
        Matrix2cd r1 = r1_inv.inverse();
        Matrix2cd r2 = m_b12[n] - f_np1 * m_a12[n];
        m_R[n] = r1 * r2;
    }

    Matrix2cd p_inc = Matrix2cd::Zero();
    p_inc(0, 0) = std::cos(theta_0) * std::cos(phi_0);
    p_inc(0, 1) = -std::sin(phi_0);
    p_inc(1, 0) = std::cos(theta_0) * std::sin(phi_0);
    p_inc(1, 1) = std::cos(phi_0);
    Matrix2cd p_inc_inv = p_inc.inverse();

    // Finally the R matrix output...
    Matrix2cd m_r_ps = p_inc_inv * m_R[0] * p_inc;

    // looping for T over the layers
    std::vector<Matrix2cd> m_Tn(m_c12.size());
    Matrix2cd m_T = Matrix2cd::Identity();
    for (int n = e_list_3x3.size() - 3; n >= 0; --n)
    {

        // building the first factor for the T_n
        Matrix2cd f1_inv = m_a12[n + 1] * m_c12[n + 1] + m_a34[n + 1] * m_c34[n + 1] * m_R[n + 1];
        Matrix2cd f1 = f1_inv.inverse();

        // building the second factor for the T_n
        Matrix2cd f2 = m_a12[n] + m_a34[n] * m_R[n];

        // T_n
        m_Tn[n] = f1 * f2;

        // T
        m_T = m_T * m_Tn[n];
    }

    // rotating m_T to the s,p states
    double theta_s = std::asin(std::sin(theta_0) * n_0 / n_s);
    Matrix2cd p_sub = Matrix2cd::Zero();
    p_sub(0, 0) = std::cos(theta_s) * std::cos(phi_0);
    p_sub(0, 1) = -std::sin(phi_0);
    p_sub(1, 0) = std::cos(theta_s) * std::sin(phi_0);
    p_sub(1, 1) = std::cos(phi_0);
    Matrix2cd p_sub_inv = p_sub.inverse();

    // Finally the T matrix output...
    Matrix2cd m_t_ps = p_sub_inv * m_T * p_inc;

    // initializing the fields
    // Assuming Matrix4x3x2cd is a 3D matrix type with complex double elements
    std::vector<Matrix4x3x2cd> m_En(e_list_3x3.size(), Matrix4x3x2cd::Zero());
    std::vector<Matrix4x3x2cd> m_Hn = m_En;
    m_En[0](1, 1, 0) = 1.0;                // TE polarization
    m_En[0](0, 0, 1) = -std::cos(theta_0); // TM polarization

    for (int n = 0; n < e_list_3x3.size(); ++n)
    {

        // forward electric and magnetic fields in the i_th layer

        // Ey1 Ez1
        m_En[n](0, 1, 0) = m_a[n](0) * m_En[n](0, 0, 0); // TE
        m_En[n](0, 2, 0) = m_b[n](0) * m_En[n](0, 0, 0);
        m_En[n](0, 1, 1) = m_a[n](0) * m_En[n](0, 0, 1); // TM
        m_En[n](0, 2, 1) = m_b[n](0) * m_En[n](0, 0, 1);

        // Ex2 Ez2
        m_En[n](1, 0, 0) = m_a[n](1) * m_En[n](1, 1, 0); // TE
        m_En[n](1, 2, 0) = m_b[n](1) * m_En[n](1, 1, 0);
        m_En[n](1, 0, 1) = m_a[n](1) * m_En[n](1, 1, 1); // TM
        m_En[n](1, 2, 1) = m_b[n](1) * m_En[n](1, 1, 1);

        // Hx1 Hy1 Hz1
        m_Hn[n](0, 0, 0) = m_b12[n](0, 0) * m_En[n](0, 0, 0) / k0; // TE
        m_Hn[n](0, 1, 0) = m_b12[n](1, 0) * m_En[n](0, 0, 0) / k0;
        m_Hn[n](0, 2, 0) = (-ky + kx * m_a[n](0)) * m_En[n](0, 0, 0) / k0;
        m_Hn[n](0, 0, 1) = m_b12[n](0, 0) * m_En[n](0, 0, 1) / k0; // TM
        m_Hn[n](0, 1, 1) = m_b12[n](1, 0) * m_En[n](0, 0, 1) / k0;
        m_Hn[n](0, 2, 1) = (-ky + kx * m_a[n](0)) * m_En[n](0, 0, 1) / k0;

        // Hx2 Hy2 Hz2
        m_Hn[n](1, 0, 0) = m_b12[n](0, 1) * m_En[n](1, 1, 0) / k0; // TE
        m_Hn[n](1, 1, 0) = m_b12[n](1, 1) * m_En[n](1, 1, 0) / k0;
        m_Hn[n](1, 2, 0) = (-ky * m_a[n](1) + kx) * m_En[n](1, 1, 0) / k0;
        m_Hn[n](1, 0, 1) = m_b12[n](0, 1) * m_En[n](1, 1, 1) / k0; // TM
        m_Hn[n](1, 1, 1) = m_b12[n](1, 1) * m_En[n](1, 1, 1) / k0;
        m_Hn[n](1, 2, 1) = (-ky * m_a[n](1) + kx) * m_En[n](1, 1, 1) / k0;

        // exiting one before the last, because then I have no backpropagation
        if (n == e_list_3x3.size() - 1)
        {
            break;
        }

        // backward electric and magnetic fields in the i_th layer

        // Ex3 Ey4
        m_En[n](2, 0, 0) = m_R[n](0, 0) * m_En[n](0, 0, 0) + m_R[n](0, 1) * m_En[n](1, 1, 0); // TE
        m_En[n](3, 1, 0) = m_R[n](1, 0) * m_En[n](0, 0, 0) + m_R[n](1, 1) * m_En[n](1, 1, 0);
        m_En[n](2, 0, 1) = m_R[n](0, 0) * m_En[n](0, 0, 1) + m_R[n](0, 1) * m_En[n](1, 1, 1); // TM
        m_En[n](3, 1, 1) = m_R[n](1, 0) * m_En[n](0, 0, 1) + m_R[n](1, 1) * m_En[n](1, 1, 1);

        // Ey3 Ez3
        m_En[n](2, 1, 0) = m_a[n](2) * m_En[n](2, 0, 0);
        m_En[n](2, 2, 0) = m_b[n](2) * m_En[n](2, 0, 0);
        m_En[n](2, 1, 1) = m_a[n](2) * m_En[n](2, 0, 1);
        m_En[n](2, 2, 1) = m_b[n](2) * m_En[n](2, 0, 1);

        // Ex4 Ez4
        m_En[n](3, 0, 0) = m_a[n](3) * m_En[n](3, 1, 0);
        m_En[n](3, 2, 0) = m_b[n](3) * m_En[n](3, 1, 0);
        m_En[n](3, 0, 1) = m_a[n](3) * m_En[n](3, 1, 1);
        m_En[n](3, 2, 1) = m_b[n](3) * m_En[n](3, 1, 1);

        // Hx3 Hy3 Hz3
        m_Hn[n](2, 0, 0) = m_b34[n](0, 0) * m_En[n](2, 0, 0) / k0; // TE
        m_Hn[n](2, 1, 0) = m_b34[n](1, 0) * m_En[n](2, 0, 0) / k0;
        m_Hn[n](2, 2, 0) = (-ky + kx * m_a[n](2)) * m_En[n](2, 0, 0) / k0;
        m_Hn[n](2, 0, 1) = m_b34[n](0, 0) * m_En[n](2, 0, 1) / k0; // TM
        m_Hn[n](2, 1, 1) = m_b34[n](1, 0) * m_En[n](2, 0, 1) / k0;
        m_Hn[n](2, 2, 1) = (-ky + kx * m_a[n](2)) * m_En[n](2, 0, 1) / k0;

        // Hx4 Hy4 Hz4
        m_Hn[n](3, 0, 0) = m_b34[n](0, 1) * m_En[n](3, 1, 0) / k0; // TE
        m_Hn[n](3, 1, 0) = m_b34[n](1, 1) * m_En[n](3, 1, 0) / k0;
        m_Hn[n](3, 2, 0) = (-ky * m_a[n](3) + kx) * m_En[n](3, 1, 0) / k0;
        m_Hn[n](3, 0, 1) = m_b34[n](0, 1) * m_En[n](3, 1, 1) / k0; // TM
        m_Hn[n](3, 1, 1) = m_b34[n](1, 1) * m_En[n](3, 1, 1) / k0;
        m_Hn[n](3, 2, 1) = (-ky * m_a[n](3) + kx) * m_En[n](3, 1, 1) / k0;

        // Ex1 Ey2 n_th+1 layer
        m_En[n + 1](0, 0, 0) = m_Tn[n](0, 0) * m_En[n](0, 0, 0) + m_Tn[n](0, 1) * m_En[n](1, 1, 0); // TE
        m_En[n + 1](1, 1, 0) = m_Tn[n](1, 0) * m_En[n](0, 0, 0) + m_Tn[n](1, 1) * m_En[n](1, 1, 0);
        m_En[n + 1](0, 0, 1) = m_Tn[n](0, 0) * m_En[n](0, 0, 1) + m_Tn[n](0, 1) * m_En[n](1, 1, 1); // TM
        m_En[n + 1](1, 1, 1) = m_Tn[n](1, 0) * m_En[n](0, 0, 1) + m_Tn[n](1, 1) * m_En[n](1, 1, 1);
    }

    for (auto &matrix : m_Kn)
    {
        matrix.col(0) *= -1;
        matrix.col(2) *= -1;
    }

    for (auto &matrix : m_En)
    {
        matrix.col(0) *= -1;
        matrix.col(2) *= -1;
    }

    for (auto &matrix : m_Hn)
    {
        matrix.col(0) *= -1;
        matrix.col(2) *= -1;
    }

    using FloatVector = std::vector<double>;

    // Define the variant type that can hold any of the types we need
    using Value = std::variant<Eigen::MatrixXcd, double, FloatVector>;

    std::unordered_map<std::string, Value> dictionary;

    // Assuming m_r_ps, m_t_ps, m_Kn, m_En, m_Hn are Eigen::MatrixXcd
    // wl, theta_0, phi_0 are double
    // e_list_3x3 is a 3x3 Eigen::MatrixXcd
    // d_list is FloatVector

    dictionary["m_r_ps"] = m_r_ps;
    dictionary["m_t_ps"] = m_t_ps;
    dictionary["m_Kn"] = m_Kn;
    dictionary["m_En"] = m_En;
    dictionary["m_Hn"] = m_Hn;
    dictionary["wl"] = wl;
    dictionary["theta_0"] = theta_0;
    dictionary["phi_0"] = phi_0;
    dictionary["e_list_3x3"] = e_list_3x3;
    dictionary["d_list"] = d_list;

    return result_dictionary;
}

int main()
{
    /*  Example input  */
    float wl = 400;
    float theta_0 = 0.0017453292519943296;
    float phi_0 = 0;

    std::vector<Eigen::MatrixXcf> e_list_3x3 = {
        (Eigen::MatrixXcf(3, 3) << std::complex<float>(-23.38459267, 4.76594931), std::complex<float>(0, 0), std::complex<float>(0, 0),
         std::complex<float>(0, 0), std::complex<float>(-23.38459267, 4.76594931), std::complex<float>(0, 0),
         std::complex<float>(0, 0), std::complex<float>(0, 0), std::complex<float>(-23.38459267, 4.76594931))
            .finished(),
        (Eigen::MatrixXcf(3, 3) << std::complex<float>(2.25, 0), std::complex<float>(0, 0), std::complex<float>(0, 0),
         std::complex<float>(0, 0), std::complex<float>(2.25, 0), std::complex<float>(0, 0),
         std::complex<float>(0, 0), std::complex<float>(0, 0), std::complex<float>(2.25, 0))
            .finished()};

    std::vector<float> d_list = {0, 20, 105, 100, 0};

    rt(wl, theta_0, phi_0, e_list_3x3, d_list);

    return 0;
}
