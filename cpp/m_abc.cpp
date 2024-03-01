#include "Eigen/Dense"
#include <array>
#include <complex>
#include <iostream>
#include <vector>

std::tuple<Eigen::Vector4cd, Eigen::Vector4cd, Eigen::Matrix2cd, Eigen::Matrix2cd, Eigen::Matrix2cd, Eigen::Matrix2cd, Eigen::Matrix2cd, Eigen::Matrix2cd> m_abc(std::complex<double> k0, std::complex<double> kx, std::complex<double> ky, Eigen::Vector4cd v_kz, Eigen::Matrix4cd v_e, double d)
{
    // matrix allocation
    Eigen::Matrix2cd m_a12 = Eigen::Matrix2cd::Identity();
    Eigen::Matrix2cd m_a34 = Eigen::Matrix2cd::Identity();
    Eigen::Matrix2cd m_b12 = Eigen::Matrix2cd::Zero();
    Eigen::Matrix2cd m_b34 = Eigen::Matrix2cd::Zero();
    Eigen::Matrix2cd m_c12 = Eigen::Matrix2cd::Zero();
    Eigen::Matrix2cd m_c34 = Eigen::Matrix2cd::Zero();

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
    Eigen::Vector4cd v_a;
    v_a << a1, a2, a3, a4;

    Eigen::Vector4cd v_b;
    v_b << b1, b2, b3, b4;

    return {v_a, v_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34};
}

/*
int main()
{
    // Define test values
    std::complex<double> k0(0.015707963267948967, 0.0);
    std::complex<double> kx(-0.00002741555386204, 0.0);
    std::complex<double> ky(0.0, 0.0);

    Eigen::Vector4cd v_kz{std::complex<double>(-0.01570794, 0.0), std::complex<double>(-0.01570794, 0.0), std::complex<double>(0.01570794, 0.0), std::complex<double>(0.01570794, 0.0)};

    Eigen::Matrix4cd v_e{{std::complex<double>(0.99999848, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(-0.00174533, 0.0), std::complex<double>(0.0, 0.0)},
                         {std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0)},
                         {std::complex<double>(-0.99999848, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(-0.00174533, 0.0), std::complex<double>(0.0, 0.0)},
                         {std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0)}};

    double d = 0.0;

    auto [v_a, v_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34] = m_abc(k0, kx, ky, v_kz, v_e, d);

    // Print the results
    std::cout << "v_a: " << std::endl
              << v_a << std::endl;
    std::cout << "v_b: " << std::endl
              << v_b << std::endl;
    std::cout << "m_a12: " << std::endl
              << m_a12 << std::endl;
    std::cout << "m_a34: " << std::endl
              << m_a34 << std::endl;
    std::cout << "m_b12: " << std::endl
              << m_b12 << std::endl;
    std::cout << "m_b34: " << std::endl
              << m_b34 << std::endl;
    std::cout << "m_c12: " << std::endl
              << m_c12 << std::endl;
    std::cout << "m_c34: " << std::endl
              << m_c34 << std::endl;

    // v_a:
    //   (0,0)
    //   (0,0)
    // (-0,-0)
    //   (0,0)
    // v_b:
    // (-0.00174533,0)
    //   (0,0)
    // (0.00174533,-0)
    //   (0,0)
    // m_a12:
    // (1,0) (0,0)
    // (0,0) (1,0)
    // m_a34:
    //   (1,0)   (0,0)
    // (-0,-0)   (1,0)
    // m_b12:
    // (0,0) (0.0157079,0)
    // (-0.015708,0)         (0,0)
    // m_b34:
    //  (0,0) (-0.0157079,0)
    //   (0.015708,0)          (0,0)
    // m_c12:
    // (1,-0)  (0,0)
    //  (0,0) (1,-0)
    // m_c34:
    // (1,0) (0,0)
    // (0,0) (1,0)

    return 0;
}
*/