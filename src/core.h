#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Dense"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Eigen;

/**
 * @brief Function that returns the kernel (or nullspace) that is the linear
 * maps that are mapped on the null vector
 *
 * @param A Matrix we want to calculate the kernel for
 * @param atol Tolerance to what is regarded as zero from the singular value
 * decomposition. Here the results from Eigen in C++ are very different from the
 * results from scipy and python
 * @return Returns the kernel of the matrix
 */
MatrixXcd nullspace(MatrixXcd A, double atol = 1e-5)
{
    // Compute the singular value decomposition
    JacobiSVD<MatrixXcd> svd(A, ComputeThinU | ComputeThinV);

    // Get the singular values
    // Singular values are fine, however, the tolerances seem to be different in C++ than in Python
    VectorXd singular_values = svd.singularValues();

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

    // Take care here! In the python code it says that we only want to select
    // specific columns of V.
    MatrixXcd null_space = svd.matrixV()(Eigen::all, mask);

    return null_space;
}

/**
 * @brief Function that returns the eigenvalues of the wavevector in the z
 *
 * @param k0, kx, ky The wavevector modulus and in plane components
 * @param m_eps The dielectric tensor of the material
 * @return Returns eigenvalues
 */
Vector4cd kz_eigenvalues(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Matrix3cd m_eps)
{

    Vector4cd v_kz;
    // Define accuracy of floating point comparisons
    double epsilon = 1e-9;

    // Are we diagonal and isotropic?
    bool diag_flag = m_eps.isApprox(m_eps.diagonal().asDiagonal().toDenseMatrix());
    bool iso_flag = (m_eps(0, 0) == m_eps(1, 1)) && (m_eps(1, 1) == m_eps(2, 2));

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

    return v_kz;
}
/**
 * @brief Function that returns the eigenvectors of the wavevector in the z
 *
 * @param k0, kx, ky The wavevector modulus and in plane components
 * @param v_kz The eigenvalues of the wavevector in the z direction
 * @param m_eps The dielectric tensor of the material
 * @return Returns eigenvectors and again the eigenvalues
 */
std::pair<MatrixXcd, Vector4cd> kz_eigenvectors(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Vector4cd v_kz, Matrix3cd m_eps)
{
    // Initializing vector and matrix
    MatrixXcd v_e = MatrixXcd::Zero(4, 3);
    Matrix3cd m_k = Matrix3cd::Zero();
    MatrixXcd m_char = MatrixXcd::Zero(3, 3);

    // Are we diagonal and isotropic?
    bool diag_flag = m_eps.isApprox(m_eps.diagonal().asDiagonal().toDenseMatrix());
    bool iso_flag = (m_eps(0, 0) == m_eps(1, 1)) && (m_eps(1, 1) == m_eps(2, 2));

    // Diagonal isotropic material
    if (diag_flag && iso_flag)
    {
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
        for (int m = 0; m < 4; ++m)
        {
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

            // Characteristic matrix
            m_char = m_k * m_k / (k0 * k0);
            m_char = m_char + m_eps;

            // Calculating the null space. The tolerances in C++ somehow
            // seem to be lower for the nullspace stuff
            MatrixXcd null_space = nullspace(m_char);

            v_e.row(m) = null_space.col(0);
        }

        // cleaning small elements from the eigenvectors
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

        // eigenvector swapping to get appropriate polarization states
        if (std::abs(v_e(0, 0)) == 0.0)
        {

            Vector3cd swap_e = v_e.row(0);
            v_e.row(0) = v_e.row(1);
            v_e.row(1) = swap_e;
            std::complex<double> swap_kz = v_kz[0];
            v_kz[0] = v_kz[1];
            v_kz[1] = swap_kz;
        }

        if (std::abs(v_e(2, 0)) == 0.0)
        {
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

/**
 * @brief Function to calculate matrix elements of the current layer from
 * incident wave vectors and eigenvalues and vectors
 *
 * @param k0, kx, ky The wavevector modulus and in plane components
 * @param v_kz The eigenvalues of the wavevector in the z direction
 * @param v_e The eigenvectors of the wavevector in the z direction
 * @param d Thickness of the layer
 * @return Returns a whole number of matrices for the TMM
 */
std::tuple<Vector4cd, Vector4cd, Matrix2cd, Matrix2cd, Matrix2cd, Matrix2cd, Matrix2cd, Matrix2cd> m_abc(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Vector4cd v_kz, MatrixXcd v_e, double d)
{
    // matrix allocation
    Matrix2cd m_a12 = Matrix2cd::Identity();
    Matrix2cd m_a34 = Matrix2cd::Identity();
    Matrix2cd m_b12 = Matrix2cd::Zero();
    Matrix2cd m_b34 = Matrix2cd::Zero();
    Matrix2cd m_c12 = Matrix2cd::Zero();
    Matrix2cd m_c34 = Matrix2cd::Zero();

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

    m_b34(0, 0) = -v_kz[2] * a3 + ky * b3;
    m_b34(0, 1) = -v_kz[3] + ky * b4;
    m_b34(1, 0) = v_kz[2] - kx * b3;
    m_b34(1, 1) = v_kz[3] * a4 - kx * b4;

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

/**
 * @brief Function to calculate transmisison and reflection matrices for each
 * layer (which needs the next layers matrices as well). It can be recursively
 * called to calculate the total transmission and reflection matrices of layered
 * stacks. The function is a bit messy but allows for recursive calculations and
 * thereby circumvents the need of weird data structures.
 *
 * @param m_a12, m_a34, m_b12, m_b34 Matrices for the current layer
 * @param m_a12_np1, m_a34_np1, m_b12_np1, m_b34_np1 Matrices for the next layer
 * @param m_R_np1 The reflection matrix for the next layer
 * @param m_T The total transmission matrix
 * @return Total reflection and transmission matrices
 */
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

    // Define matrix R and T for this iteration
    Matrix2cd m_R = Matrix2cd::Zero();
    Matrix2cd m_Tn = Matrix2cd::Zero();

    Matrix2cd f1 = m_b12_np1 * m_c12_np1 + (m_b34_np1 * m_c34_np1) * m_R_np1;
    Matrix2cd f2_inv = m_a12_np1 * m_c12_np1 + (m_a34_np1 * m_c34_np1) * m_R_np1;
    Matrix2cd f2 = f2_inv.inverse();

    Matrix2cd f_np1 = f1 * f2;
    Matrix2cd r1_inv = f_np1 * m_a34 - m_b34;
    Matrix2cd r1 = r1_inv.inverse();
    Matrix2cd r2 = m_b12 - f_np1 * m_a12;

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

/**
 * @brief Function to replace the numpy real_if_close
 *
 * @param c The complex number to check
 * @param tol Tolerance to check for
 * @return The real part of the complex number
 */
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

/**
 * @brief Checks if a matrix is real, diagonal, and isotropic.
 *
 * @param matrix The matrix to check.
 * @return true if the matrix is real, diagonal, and isotropic; false otherwise.
 */
bool isRealDiagonalIsotropic(const MatrixXcd &matrix)
{
    bool diag_flag = matrix.isApprox(matrix.diagonal().asDiagonal().toDenseMatrix());
    bool iso_flag = (matrix(0, 0) == matrix(1, 1)) && (matrix(0, 0) == matrix(2, 2));
    bool real_flag = std::abs(matrix(0, 0)) == std::real(matrix(0, 0));

    return diag_flag && iso_flag && real_flag;
}

std::pair<Matrix2cd, Matrix2cd> calculate_tr(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0)
{
    // Incident medium and substrate hva to be real, diagonal and isotropic
    // incident medium check
    MatrixXcd firstMatrix = e_list_3x3[0];
    if (!isRealDiagonalIsotropic(firstMatrix))
    {
        throw std::runtime_error("Incident medium must be real, diagonal and isotropic");
    }

    std::complex<double> n_0 = std::sqrt(firstMatrix(0, 0));

    // substrate check
    MatrixXcd lastMatrix = e_list_3x3[e_list_3x3.size() - 1];
    if (!isRealDiagonalIsotropic(lastMatrix))
    {
        throw std::runtime_error("Substrate must be real, diagonal and isotropic");
    }

    std::complex<double> n_s = std::sqrt(lastMatrix(0, 0));

    // wavevector modulus and in plane components
    std::complex<double> k0 = 2.0 * M_PI / wavelength;
    std::complex<double> kx = -k0 * n_0.real() * sin(theta_0) * cos(phi_0);
    std::complex<double> ky = -k0 * n_0.real() * sin(theta_0) * sin(phi_0);

    Vector4cd v_kz1 = kz_eigenvalues(k0, kx, ky, e_list_3x3.back());
    auto [v_e, v_kz] = kz_eigenvectors(k0, kx, ky, v_kz1, e_list_3x3.back());

    auto [m_a_np1, m_b_np1, m_a12_np1, m_a34_np1, m_b12_np1, m_b34_np1, m_c12_np1, m_c34_np1] = m_abc(k0, kx, ky, v_kz, v_e, d_list.back());

    Matrix2cd m_T = Matrix2cd::Identity();
    Matrix2cd m_R_np1 = Matrix2cd::Zero();
    Matrix2cd m_R_0 = Matrix2cd::Zero();

    // Worked up to here
    // Now iterate over all layers starting from the second last going backwards
    for (int i = d_list.size() - 2; i >= 0; --i)
    {
        v_kz1 = kz_eigenvalues(k0, kx, ky, e_list_3x3[i]);

        auto [v_e, v_kz] = kz_eigenvectors(k0, kx, ky, v_kz1, e_list_3x3[i]);

        auto [m_a, m_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34] = m_abc(k0, kx, ky, v_kz, v_e, d_list[i]);

        std::pair<Matrix2cd, Matrix2cd> result3 = calculate_tr_per_layer(m_a12, m_a34, m_b12, m_b34, m_a12_np1, m_a34_np1, m_b12_np1, m_b34_np1, m_c12_np1, m_c34_np1, m_R_np1, m_T);

        Matrix2cd m_R = result3.first;
        m_T = result3.second;

        if (i == 0)
        {
            m_R_0 = m_R;
        }

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

    // The reflection matrix
    Matrix2cd m_r_ps = p_inc_inv * m_R_0 * p_inc;

    // rotating m_T to the s,p states
    double theta_s = asin(real_if_close(sin(theta_0) * n_0.real() / n_s.real()));
    Matrix2cd p_sub = Matrix2cd::Zero();
    p_sub(0, 0) = cos(theta_s) * cos(phi_0);
    p_sub(0, 1) = -sin(phi_0);
    p_sub(1, 0) = cos(theta_s) * sin(phi_0);
    p_sub(1, 1) = cos(phi_0);
    Matrix2cd p_sub_inv = p_sub.inverse();

    // The transmission matrix
    Matrix2cd m_t_ps = p_sub_inv * m_T * p_inc;

    return {m_r_ps, m_t_ps};
};

// int main()
// {
//     double wavelength = 400;
//     double theta_0 = 0.0017453292519943296;
//     double phi_0 = 0;

//     std::vector<Matrix3cd> e_list_3x3 = {
//         (Matrix3cd(3, 3) << std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0))
//             .finished(),

//         (Matrix3cd(3, 3) << std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931))
//             .finished(),

//         (Matrix3cd(3, 3) << std::complex<double>(2.90627602, 0.84588284), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(4.3259341, 4.83066447), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(2.98285307, 1.06542853))
//             .finished(),

//         (Matrix3cd(3, 3) << std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(-23.38459267, 4.76594931))
//             .finished(),

//         (Matrix3cd(3, 3) << std::complex<double>(2.25, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0.0), std::complex<double>(0.0, 0.0),
//          std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(2.25, 0.0))
//             .finished()};

//     std::vector<double> d_list = {0, 20, 105, 100, 0};

//     auto [m_r_ps, m_t_ps] = calculate_tr(e_list_3x3, d_list, wavelength, theta_0, phi_0);
//     std::cout << "reflection matrix: " << m_r_ps << std::endl;
//     std::cout << "transmission matrix: " << m_t_ps << std::endl;

//     return 0;
// }