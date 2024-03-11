#include "Eigen/Dense"
#include <array>
#include <complex>
#include <iostream>
#include <vector>
#include "nullspace.cpp"

// using namespace std;
using namespace Eigen;

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

/*
// Assuming the nullspace function is defined here...
int main()
{
    // Define test values
    std::complex<double> k0(0.015707963267948967, 0.0);
    std::complex<double> kx(-2.741555386204003e-05, 0.0);
    std::complex<double> ky(0.0, 0.0);

    Eigen::Vector4cd v_kz{std::complex<double>(-0.01570794, 0.0), std::complex<double>(-0.01570794, 0.0), std::complex<double>(0.01570794, 0.0), std::complex<double>(0.01570794, 0.0)};

    Eigen::Matrix3cd m_eps{{std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0)},
                           {std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)},
                           {std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)}};

    auto [v_e, v_kz2] = kz_eigenvectors(k0, kx, ky, v_kz, m_eps);

    // Print the results
    std::cout << "Nullspace: " << std::endl
              << v_e << std::endl;

    // v_e
    // (0.999998,-0)           (0,0) (-0.00174533,0)
    // (0,0)           (1,0)           (0,0)
    // (-0.999998,-0)           (0,0) (-0.00174533,0)
    // (0,0)           (1,0)           (0,0)

    return 0;
}
*/