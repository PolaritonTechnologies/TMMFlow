#include <iostream>
#include <map>
#include <vector>
#include <complex>
#include "Eigen/Dense"
#include <unordered_map>
#include <variant>
// #include "nullspace.cpp"
// #include "kz_eigenvectors.cpp"
// #include "kz_eigenvalues.cpp"
// #include "m_abc.cpp"
// #include "json.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Eigen;

MatrixXcd nullspace(MatrixXcd A, double atol = 1e-5)
{
    // Compute the singular value decomposition
    JacobiSVD<MatrixXcd> svd(A, ComputeThinU | ComputeThinV);

    // Get the singular values
    // Singular values are fine, however, the tolerances seem to be different in C++ than in Python
    VectorXd singular_values = svd.singularValues();

    // std::cout << "singular values: " << std::endl
    //           << singular_values << std::endl;
    // std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl
    //           << svd.matrixU() << std::endl;
    // std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl
    //           << svd.matrixV() << std::endl;

    // Define a mask that cuts off all values that are smaller than atol
    std::vector<int> mask;

    // Iterate over all singular values and check if they are smaller than atol
    for (int i = 0; i < singular_values.size(); i++)
    {
        if (singular_values(i) <= atol)
        {
            // mask[i] = i;
            mask.insert(std::end(mask), i);
        }
    }
    for (int i = 0; i < mask.size(); i++)
    {
        std::cout << mask[i] << std::endl;
    }

    // Take care here! In the python code it says that we only want to select
    // specific columns of V.
    MatrixXcd null_space = svd.matrixV()(Eigen::all, mask);

    std::cout << "Nullspace: " << std::endl
              << null_space << std::endl;

    return null_space;
}

Vector4cd kz_eigenvalues(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Matrix3cd m_eps)
{

    Vector4cd v_kz;
    // Define accuracy of floating point comparisons
    double epsilon = 1e-9;

    bool diag_flag = m_eps.isApprox(m_eps.diagonal().asDiagonal().toDenseMatrix());
    bool iso_flag = (m_eps(0, 0) == m_eps(1, 1)) && (m_eps(1, 1) == m_eps(2, 2));

    // Are we diagonal and isotropic?

    /*
    // Create a diagonal matrix from the diagonal of the input matrix
    Eigen::MatrixXcd diagonalMatrix = m_eps.diagonal().asDiagonal();

    bool diag_flag = true;
    // Check if the input matrix is equal to its diagonal matrix
    // The norm of their difference should be close to zero for them to be considered equal
    if ((m_eps - diagonalMatrix).norm() > epsilon)
    {
        diag_flag = false;
    }

    std::complex<double> firstElement = m_eps.diagonal()(0);
    bool iso_flag = true;

    // Check if all elements of the diagonal of m_eps are close to the first element
    for (int i = 1; i < m_eps.diagonal().size(); ++i)
    {
        if (std::abs(m_eps.diagonal()(i) - firstElement) > epsilon)
        {
            iso_flag = false;
            break; // No need to check further if we found a non-equal element
        }
    }
    */

    // diagonal isotropic material
    if (diag_flag && iso_flag)
    {
        std::complex<double> kz = sqrt(pow(k0, 2) * m_eps(0, 0) - pow(kx, 2) - pow(ky, 2));
        std::fill(v_kz.begin(), v_kz.begin() + 2, -kz);
        std::fill(v_kz.begin() + 2, v_kz.end(), kz);
    }

    // general material
    else
    {
        // coefficients for the quartic equation
        std::complex<double> A = (kx / k0) * (((m_eps(0, 2) + m_eps(2, 0)) / m_eps(2, 2)) + (ky / k0) * ((m_eps(1, 2) + m_eps(2, 1)) / m_eps(2, 2)));

        std::complex<double> B = (std::pow((kx / k0), 2) * (1.0 + m_eps(0, 0) / m_eps(2, 2)) + std::pow((ky / k0), 2) * (1.0 + m_eps(1, 1) / m_eps(2, 2)) + ((kx * ky) / std::pow(k0, 2)) * (m_eps(0, 1) + m_eps(1, 0)) / m_eps(2, 2) + ((m_eps(0, 2) * m_eps(2, 0) + m_eps(1, 2) * m_eps(2, 1)) / m_eps(2, 2) - m_eps(0, 0) - m_eps(1, 1)));

        std::complex<double> C = ((std::pow(kx, 2) + std::pow(ky, 2)) / std::pow(k0, 2) * ((kx / k0) * (m_eps(0, 2) + m_eps(2, 0)) / m_eps(2, 2) + (ky / k0) * (m_eps(1, 2) + m_eps(2, 1)) / m_eps(2, 2)) + (kx / k0) * ((m_eps(0, 1) * m_eps(1, 2) + m_eps(1, 0) * m_eps(2, 1)) / m_eps(2, 2) - (m_eps(1, 1) / m_eps(2, 2)) * (m_eps(0, 2) + m_eps(2, 0))) + (ky / k0) * ((m_eps(0, 1) * m_eps(2, 0) + m_eps(1, 0) * m_eps(0, 2)) / m_eps(2, 2) - (m_eps(0, 0) / m_eps(2, 2)) * (m_eps(1, 2) + m_eps(2, 1))));

        std::complex<double> D1 = ((std::pow(kx, 2) + std::pow(ky, 2)) / std::pow(k0, 2)) * (std::pow((kx / k0), 2) * m_eps(0, 0) / m_eps(2, 2) + std::pow((ky / k0), 2) * m_eps(1, 1) / m_eps(2, 2) + ((kx * ky) / std::pow(k0, 2)) * (m_eps(0, 1) + m_eps(1, 0)) / m_eps(2, 2) - m_eps(0, 0) * m_eps(1, 1) / m_eps(2, 2));

        std::complex<double> D2 = std::pow((kx / k0), 2) * ((m_eps(0, 1) * m_eps(1, 0) + m_eps(0, 2) * m_eps(2, 0)) / m_eps(2, 2) - m_eps(0, 0));

        std::complex<double> D3 = std::pow((ky / k0), 2) * ((m_eps(0, 1) * m_eps(1, 0) + m_eps(1, 2) * m_eps(2, 1)) / m_eps(2, 2) - m_eps(1, 1));

        std::complex<double> D4 = ((kx * ky) / std::pow(k0, 2)) * ((m_eps(0, 2) * m_eps(2, 1) + m_eps(2, 0) * m_eps(1, 2)) / m_eps(2, 2) - m_eps(0, 1) - m_eps(1, 0));

        std::complex<double> D5 = (m_eps(0, 0) * m_eps(1, 1) + (m_eps(0, 1) * m_eps(1, 2) * m_eps(2, 0) + m_eps(1, 0) * m_eps(2, 1) * m_eps(0, 2)) / m_eps(2, 2) - m_eps(0, 1) * m_eps(1, 0) - (m_eps(0, 0) / m_eps(2, 2)) * m_eps(1, 2) * m_eps(2, 1) - (m_eps(1, 1) / m_eps(2, 2)) * m_eps(0, 2) * m_eps(2, 0));

        std::complex<double> D = D1 + D2 + D3 + D4 + D5;
        // companion matrix
        MatrixXcd m_comp = MatrixXcd::Zero(4, 4);
        m_comp(1, 0) = 1.0;
        m_comp(2, 1) = 1.0;
        m_comp(3, 2) = 1.0;
        m_comp(0, 3) = -D;
        m_comp(1, 3) = -C;
        m_comp(2, 3) = -B;
        m_comp(3, 3) = -A;

        // eigenvalues
        ComplexEigenSolver<MatrixXcd> ces;
        ces.compute(m_comp);
        VectorXcd eigvals = ces.eigenvalues();
        for (int i = 0; i < 4; ++i)
        {
            v_kz[i] = k0 * eigvals[i];
        }
    }

    // We should first round the numbers as otherwise the sorting is
    // meaningless. If the numbers are close to zero, they shouldn't be sorted
    // (in this case I chose 1e-6 as a threshold)
    // Round the real and imaginary parts of each complex number
    for (auto &kz : v_kz)
    {
        kz = std::complex<double>(std::round(kz.real() / epsilon) * epsilon, std::round(kz.imag() / epsilon) * epsilon);
    }

    // Sort the complex numbers by their imaginary parts
    std::sort(v_kz.begin(), v_kz.end(), [](std::complex<double> a, std::complex<double> b)
              { return std::imag(a) < std::imag(b); });

    // output sorted by imaginary part
    // std::sort(v_kz.begin(), v_kz.end(), [](std::complex<double> a, std::complex<double> b)
    //   { return std::imag(a) < std::imag(b); });

    // for (const auto &element : v_kz)
    // {
    // std::cout << element << std::endl;
    // }

    return v_kz;
}

std::pair<MatrixXcd, Vector4cd> kz_eigenvectors(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Vector4cd v_kz, Matrix3cd m_eps)
{
    // Initializing vector and matrix
    MatrixXcd v_e = MatrixXcd::Zero(4, 3);
    Matrix3cd m_k = Matrix3cd::Zero();
    MatrixXcd m_char = MatrixXcd::Zero(3, 3);

    // std::cout << v_e << std::endl;

    // Are we diagonal and isotropic?
    bool diag_flag = m_eps.isApprox(m_eps.diagonal().asDiagonal().toDenseMatrix());
    bool iso_flag = (m_eps(0, 0) == m_eps(1, 1)) && (m_eps(1, 1) == m_eps(2, 2));
    // std::cout << "Here-1" << std::endl;

    // Diagonal isotropic material
    if (diag_flag && iso_flag)
    {
        // std::cout << "Here-2" << std::endl;
        if (kx == std::complex<double>(0.0) && ky == std::complex<double>(0.0))
        {
            v_e.row(0) = Vector3cd(1.0, 0.0, 0.0);
            v_e.row(1) = Vector3cd(0.0, 1.0, 0.0);
            v_e.row(2) = Vector3cd(1.0, 0.0, 0.0);
            v_e.row(3) = Vector3cd(0.0, 1.0, 0.0);
        }
        else if (kx == std::complex<double>(0.0))
        {

            v_e.row(0) = Vector3cd(1.0, 0.0, 0.0);
            v_e.row(1) = Vector3cd(0.0, -v_kz[1], ky);
            v_e.row(2) = Vector3cd(1.0, 0.0, 0.0);
            v_e.row(3) = Vector3cd(0.0, -v_kz[3], ky);
        }
        else if (ky == std::complex<double>(0.0))
        {

            v_e.row(0) = Vector3cd(-v_kz[1], 0.0, kx);
            v_e.row(1) = Vector3cd(0.0, 1.0, 0.0);
            v_e.row(2) = Vector3cd(-v_kz[3], 0.0, kx);
            v_e.row(3) = Vector3cd(0.0, 1.0, 0.0);
        }
        else
        {

            v_e.row(0) = Vector3cd(-v_kz[1], 0.0, kx);
            v_e.row(1) = Vector3cd(-ky, kx, 0.0);
            v_e.row(2) = Vector3cd(-v_kz[3], 0.0, kx);
            v_e.row(3) = Vector3cd(-ky, kx, 0.0);
        }
    }

    // General material
    else
    {
        // std::cout << "Here-3" << std::endl;
        for (int m = 0; m < 4; ++m)
        {
            // std::cout << m << std::endl;

            // k matrix
            m_k(0, 0) = 0.0;
            m_k(0, 1) = -v_kz[m];
            m_k(0, 2) = ky;
            m_k(1, 0) = v_kz[m];
            m_k(1, 1) = 0.0;
            m_k(1, 2) = -kx;
            m_k(2, 0) = -ky;
            m_k(2, 1) = kx;
            m_k(2, 2) = 0.0;
            // std::cout << "Here-3a" << std::endl;

            // Characteristic matrix
            m_char = m_k * m_k / (k0 * k0);
            m_char = m_char + m_eps;
            // std::cout << "Here-3b" << std::endl;
            std::cout << m_char << std::endl;

            // Calculating the null space
            // This here is the problem as it becomes empty when it shouldnt.
            // The tolerances in C++ somehow seem to be lower for the nullspace
            // stuff
            MatrixXcd null_space = nullspace(m_char, 1e-5);
            // std::cout << "Here-3c" << std::endl;

            // std::cout << v_e.row(m) << std::endl;

            v_e.row(m) = null_space.col(0);
            // std::cout << "Here-3d" << std::endl;
        }

        // std::cout << "Here-4" << std::endl;
        // c    leaning small elements from the eigenvectors
        for (int m = 0; m < 4; ++m)
        {
            std::complex<double> max_e = v_e.row(m).cwiseAbs().maxCoeff();
            Vector3cd v_e_rel = v_e.row(m).cwiseAbs() / max_e;
            for (int i = 0; i < 3; ++i)
            {
                if (std::abs(v_e_rel[i]) < 1.0e-12)
                {
                    v_e(m, i) = 0.0;
                }
            }
        }

        // e    igenvector swapping to get appropriate polarization states
        if (std::abs(v_e(0, 0)) == 0.0)
        {

            // std::cout << "Here-5" << std::endl;
            Vector3cd swap_e = v_e.row(0);
            v_e.row(0) = v_e.row(1);
            v_e.row(1) = swap_e;
            std::complex<double> swap_kz = v_kz[0];
            v_kz[0] = v_kz[1];
            v_kz[1] = swap_kz;
        }

        if (std::abs(v_e(2, 0)) == 0.0)
        {
            // std::cout << "Here-6" << std::endl;
            Vector3cd swap_e = v_e.row(2);
            v_e.row(2) = v_e.row(3);
            v_e.row(3) = swap_e;
            std::complex<double> swap_kz = v_kz[2];
            v_kz[2] = v_kz[3];
            v_kz[3] = swap_kz;
        }
    }

    // Normalizing eigenvectors
    for (int m = 0; m < 4; ++m)
    {
        v_e.row(m) = v_e.row(m) / v_e.row(m).norm();
    }

    return {v_e, v_kz};
}

std::tuple<Vector4cd, Vector4cd, Matrix2cd, Matrix2cd, Matrix2cd, Matrix2cd, Matrix2cd, Matrix2cd> m_abc(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Vector4cd v_kz, MatrixXcd v_e, double d)
{
    // matrix allocation
    Matrix2cd m_a12 = Matrix2cd::Identity();
    Matrix2cd m_a34 = Matrix2cd::Identity();
    Matrix2cd m_b12 = Matrix2cd::Zero();
    Matrix2cd m_b34 = Matrix2cd::Zero();
    Matrix2cd m_c12 = Matrix2cd::Zero();
    Matrix2cd m_c34 = Matrix2cd::Zero();

    for (int i = 0; i < 4; i++)
    {
        std::cout << v_e.row(i) << std::endl;
    };

    std::complex<double> a1 = v_e(0, 1) / v_e(0, 0);
    std::complex<double> a2 = v_e(1, 0) / v_e(1, 1);

    // a34 matrix
    std::complex<double> a3 = v_e(2, 1) / v_e(2, 0);
    std::complex<double> a4 = v_e(3, 0) / v_e(3, 1);

    m_a34(0, 1) = a4;
    m_a34(1, 0) = a3;

    // b12 matrix
    std::complex<double> b1 = v_e(0, 2) / v_e(0, 0);
    std::complex<double> b2 = v_e(1, 2) / v_e(1, 1);

    m_b12(0, 0) = -v_kz[0] * a1 + ky * b1;
    m_b12(0, 1) = -v_kz[1] + ky * b2;
    m_b12(1, 0) = v_kz[0] - kx * b1;
    m_b12(1, 1) = v_kz[1] * a2 - kx * b2;

    // b34 matrix
    std::complex<double> b3 = v_e(2, 2) / v_e(2, 0);
    std::complex<double> b4 = v_e(3, 2) / v_e(3, 1);

    // Something is wrong here! It might be that the order of indices is wrong (0,1) --> (1,0) or something is odd with the minus sign
    m_b34(0, 0) = -v_kz[2] * a3 + ky * b3;
    m_b34(0, 1) = -v_kz[3] + ky * b4;
    m_b34(1, 0) = v_kz[2] - kx * b3;
    m_b34(1, 1) = v_kz[3] * a4 - kx * b4;

    std::cout << "m_b34: " << std::endl
              << m_b34 << std::endl;

    // c12 matrix
    m_c12(0, 0) = std::exp(std::complex<double>(0, 1) * v_kz[0] * d);
    m_c12(1, 1) = std::exp(std::complex<double>(0, 1) * v_kz[1] * d);

    // c34 matrix
    m_c34(0, 0) = std::exp(std::complex<double>(0, 1) * v_kz[2] * d);
    m_c34(1, 1) = std::exp(std::complex<double>(0, 1) * v_kz[3] * d);

    // pure coefficient vectors
    Vector4cd v_a;
    v_a << a1, a2, a3, a4;

    Vector4cd v_b;
    v_b << b1, b2, b3, b4;

    return {v_a, v_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34};
}

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