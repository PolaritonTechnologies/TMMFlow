#include "Eigen/Dense"
#include <complex>
#include <iostream>
#include <vector>

// using namespace std;
using namespace Eigen;

MatrixXcd nullspace(MatrixXcd A, double atol = 1e-9)
{
    // Compute the singular value decomposition
    JacobiSVD<MatrixXcd> svd(A, ComputeThinU | ComputeThinV);

    // Get the singular values
    VectorXd singular_values = svd.singularValues();

    std::cout << "singular values: " << std::endl
              << singular_values << std::endl;
    std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl
              << svd.matrixU() << std::endl;
    std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl
              << svd.matrixV() << std::endl;
    // singular values: 4.23012 3.16742 3.97205e-15
    // Its right singular vectors are the columns of the thin V matrix:
    // (0.910642,-0.413193)     (-0.000661001,0.00145366)       (-0,0)
    // (0,0)                    (0,-0)                          (-0.894427,0.447214)
    // (0.000779099,0.00139393) (0.941729,0.336369)             (-0,0)
    // Its right singular vectors are the columns of the thin V matrix:
    // (0.0836123,0.996497)     (-0.000133519,-0.00159129)      (0,0)
    // (0,0)                    (0,0)                           (1,0)
    // (0.00159688,0)               (0.999999,0)                      (0,0)

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
    // for (int i = 0; i < mask.size(); i++)
    // {
    //     std::cout << "mask: " << std::endl
    //               << mask[i] << std::endl;
    // }

    // Only get the columns of V that are not masked
    MatrixXcd null_space = svd.matrixV()(Eigen::all, mask);

    return null_space.transpose();

    /*
    Initial stuff that copilot suggested
    // Find the rank of the matrix A
    int rank = 0;
    for (int i = 0; i < singular_values.size(); i++)
    {
        if (singular_values(i) > atol)
        {
            rank++;
        }
    }
    std::cout << "Rank: " << std::endl
              << rank << std::endl;
    std::cout << "A cols: " << std::endl
              << A.cols() << std::endl;

    // The null space is formed by the last (n-rank) columns of V
    MatrixXcd null_space = svd.matrixV().rightCols(A.cols() - rank);
    */
}

// Assuming the nullspace function is defined here...
int main()
{
    // Define test values
    // Matrix3cd A{{std::complex<double>(-1.41965504, -3.98478163), std::complex<double>(0, 0), std::complex<double>(0.00405773850, 0.00181320587)},
    // {std::complex<double>(0, 0), std::complex<double>(-3.55271368e-15, 1.77635684e-15), std::complex<double>(0, 0)},
    // {std::complex<double>(0.00405773850, 0.00181320587), std::complex<double>(0, 0), std::complex<double>(2.98285002, 1.06542853)}};
    Matrix3cd A{{std::complex<double>(-1.41965504, -3.98478163), std::complex<double>(0, 0), std::complex<double>(0.00405773850, 0.00181320587)},
                {std::complex<double>(0, 0), std::complex<double>(-3.55271368e-15, 1.77635684e-15), std::complex<double>(0, 0)},
                {std::complex<double>(0.00405773850, 0.00181320587), std::complex<double>(0, 0), std::complex<double>(2.98285002, 1.06542853)}};
    double atol = 1e-07;

    MatrixXcd result = nullspace(A, atol);

    // Print the results
    std::cout << "Nullspace: " << std::endl
              << result << std::endl;

    return 0;
}