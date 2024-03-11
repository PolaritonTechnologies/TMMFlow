#include "Eigen/Dense"
#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>

using namespace Eigen;

Eigen::Vector4cd kz_eigenvalues(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Matrix3cd m_eps)
{

    Eigen::Vector4cd v_kz;
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

/*
int main()
{
    // Example input
    double k0 = 0.015707963267948967;
    std::complex<double> kx(-2.741555386204003e-05, 0);
    std::complex<double> ky(0, 0);

    Eigen::Matrix3cd m_eps(3, 3);
    m_eps << std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0, 0), std::complex<double>(0, 0),
        std::complex<double>(0, 0), std::complex<double>(-23.38459267, 4.76594931), std::complex<double>(0, 0),
        std::complex<double>(0, 0), std::complex<double>(0, 0), std::complex<double>(-23.38459267, 4.76594931);

    kz_eigenvalues(k0, kx, ky, m_eps);

    return 0;
}
*/