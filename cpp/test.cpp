#include "Eigen/Dense"
#include <complex>
#include <array>
#include <iostream>

Eigen::MatrixXd nullspace(Eigen::MatrixXcd A, double atol = 1e-9)
{
    // Compute the singular value decomposition
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Get the singular values
    Eigen::VectorXd singular_values = svd.singularValues();

    // Find the rank of the matrix A
    int rank = 0;
    for (int i = 0; i < singular_values.size(); i++)
    {
        if (singular_values(i) > atol)
        {
            rank++;
        }
    }

    // The null space is formed by the last (n-rank) columns of V
    Eigen::MatrixXcd null_space = svd.matrixV().rightCols(A.cols() - rank);

    return null_space;
}

// Now call the function nullspace with the following parameters
int main()
{
    double atol = 1e-07;

    // Define the complex numbers
    std::array<std::array<std::complex<double>, 3>, 3> data = {{{{std::complex<double>(-1.41965504, -3.98478163), std::complex<double>(0.0, 0.0), std::complex<double>(0.00405773850, 0.00181320587)}},
                                                                {{std::complex<double>(0.0, 0.0), std::complex<double>(-3.55271368e-15, 1.77635684e-15), std::complex<double>(0.0, 0.0)}},
                                                                {{std::complex<double>(0.00405773850, 0.00181320587), std::complex<double>(0.0, 0.0), std::complex<double>(2.98285002, 1.06542853)}}}};

    // Now call nullspace with atol and data
    Eigen::MatrixXcd null_space = nullspace(Eigen::Map<Eigen::Matrix<std::complex<double>, 3, 3, Eigen::RowMajor>>(data[0].data()), atol);

    // Print the result
    std::cout << "The null space of the matrix is: " << std::endl
              << null_space << std::endl;

    return 0;
}