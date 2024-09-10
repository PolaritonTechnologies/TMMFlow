#include <complex>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace Eigen;

Matrix2cd coherent_block(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const char &polarization, const double &wavelength, const double &theta_0, const std::complex<double> &n_0);
Matrix2cd incoherent_block(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const char &polarization, const double &wavelength, const double &theta_0, const double &n_0, const bool &closing_layer);
std::tuple<double, double> calculate_rt(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const std::vector<bool> &incoherent, const char &polarization, const double &wavelength, const double &theta_0);
std::tuple<double, double> calculate_rt_azim_pola(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const std::vector<bool> &incoherent, const double &wavelength, const double &theta_0, const double &phi_0, const double &s_polarization_percentage);