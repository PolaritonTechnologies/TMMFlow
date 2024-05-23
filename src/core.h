#include <complex>
#include "Eigen/Dense"
using namespace Eigen;
//Troparevsky et al. Optics Express Vol. 18, Issue 24, pp. 24715-24721 (2010)
//Formalism by Al-Ghezi et al. Optics Express 35770 (2020)
//Reflectivity and Transmissivity coefficients in "Generalized matrix method for calculation of Internal light energy flux in mixed coherent and incoherent multilayers" by Emanuele Centurioni (2005)
std::tuple<double, double, Matrix2cd> calculate_rt_coherent(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, const char* polarization, double wavelength, double theta_0, double phi_0)
{
    std::complex<double> n_s = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_0 = std::sqrt(e_list_3x3[e_list_3x3.size() - 1](0, 0));
    // Wavevector modulus and in plane components
    std::complex<double> k0 = 2.0 * M_PI * n_0 / wavelength;
    std::complex<double> kx = k0 * n_0 * sin(theta_0);
    std::complex<double> kz_i;
    std::complex<double> kz_im1;
    std::complex<double> pre_factor_dynamical;
    std::complex<double> quotient_dynamical;
    // Here total_Matrix represents the electric displacement
    Matrix2cd total_Matrix = Matrix2cd::Identity();
    Matrix2cd dynamical;
    Matrix2cd propa;
    double transmissivity;  
    for (int i = 1; i < d_list.size() ; ++i){
        propa = Matrix2cd::Identity();
        if (strcmp(polarization, "s") == 0){
            // Calculate kz value according to Berreman formalism - s polarization
            kz_i = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i](1, 1)/(wavelength*wavelength) - (kx * kx));
            kz_im1 = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i-1](1, 1)/(wavelength*wavelength) - (kx * kx));
            pre_factor_dynamical = 0.5*(e_list_3x3[i](1, 1))/(e_list_3x3[i-1](1, 1));       
            // Anisotropic dynamical matrix
            dynamical(0,0) = std::complex<double>(1,0) + (kz_im1/kz_i);    
            dynamical(0,1) = std::complex<double>(1,0) - (kz_im1/kz_i);
        }
        if (strcmp(polarization, "p") == 0){
            // Calculate kz value according to Berreman formalism
            kz_i = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i](1, 1)/(wavelength*wavelength) - ((kx * kx) * e_list_3x3[i](1, 1))/(e_list_3x3[i](2, 2)));
            kz_im1 = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i-1](1, 1)/(wavelength*wavelength) - ((kx * kx) * e_list_3x3[i-1](1, 1))/(e_list_3x3[i-1](2, 2)));
            quotient_dynamical = (e_list_3x3[i-1](1, 1) * kz_i)/(e_list_3x3[i](1, 1) * kz_im1);
            pre_factor_dynamical = (0.5/quotient_dynamical)*sqrt((kz_i*kz_i + kx*kx)/(kz_im1*kz_im1 + kx*kx));     
            dynamical(0,0) = std::complex<double>(1,0) + quotient_dynamical; // Anisotropic dynamical matrix
            dynamical(0,1) = std::complex<double>(1,0) - quotient_dynamical;
        }
        dynamical(1,0) = dynamical(0,1);
        dynamical(1,1) = dynamical(0,0);
        if (i < d_list.size()-1){ //Propagation Matrix
            propa(0,0) = std::exp(std::complex<double>(0,1) * kz_i * d_list[i]);
            propa(1,1) = std::exp(std::complex<double>(0,-1) * kz_i * d_list[i]);        
        }
        total_Matrix = propa * pre_factor_dynamical * dynamical * total_Matrix;
    }
    std::complex<double> reflection_coeff = (total_Matrix(1,0)/total_Matrix(0,0)); 
    std::complex<double> transmission_coeff = pow(n_0,2)/pow(n_s,2) * std::complex<double>(1,0) / total_Matrix(0,0);    
    if (strcmp(polarization, "s") == 0){ 
    transmissivity = ((n_s*std::sqrt(1-((n_0.real()*n_0.real())/(n_s.real()*n_s.real())*sin(theta_0)*sin(theta_0)))).real())/(std::sqrt(n_0.real()*(1-sin(theta_0)*sin(theta_0)))) * std::norm(transmission_coeff);
    }
    if (strcmp(polarization, "p") == 0){
    transmissivity = ((std::conj(n_s)*std::sqrt(1-((n_0.real()*n_0.real())/(n_s.real()*n_s.real())*sin(theta_0)*sin(theta_0)))).real())/(std::sqrt(n_0.real()*(1-sin(theta_0)*sin(theta_0)))) * std::norm(transmission_coeff);
    }
    return {std::norm(reflection_coeff), transmissivity, total_Matrix};
};

std::tuple<double, double> calculate_rt(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, const char* polarization, double wavelength, double theta_0, double phi_0, double n_exit_medium)
{
    std::complex<double> n_s = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_0 = std::sqrt(e_list_3x3[d_list.size()-1](0, 0));  
    auto[reflectivity_coherent, transmissivity_coherent, total_Matrix] = calculate_rt_coherent(e_list_3x3, d_list, polarization, wavelength, theta_0, phi_0);
    std::complex<double> reflection_coeff_reverse = -(total_Matrix(0,1)/total_Matrix(0,0)); 
    std::complex<double> transmission_coeff_reverse = pow(n_s,2)/pow(n_0,2) * (total_Matrix.determinant() / total_Matrix(0,0));
    double reflectivity_no_backside_reverse = (reflection_coeff_reverse * std::conj(reflection_coeff_reverse)).real();
    // Propagation angle in the substrate
    double cos_theta_sub = std::sqrt(1 - (n_0.real() * n_0.real() * sin(theta_0)*sin(theta_0))/(n_s.real()*n_s.real()));
    // Attenuation through one pass in the incoherent layer
    double attenuation = std::exp((d_list[0] / cos_theta_sub) * sqrt(4 * M_PI * e_list_3x3[0](1, 1) / (wavelength * wavelength)).imag());
    double r_subs_exit;
    double transmissivity_no_backside_reverse;
    if (strcmp(polarization, "s") == 0){
        transmissivity_no_backside_reverse = ((n_0.real() * std::sqrt(1 - sin(theta_0) * sin(theta_0)))) / (n_s.real() * std::sqrt(1 - ((n_0.real() * n_0.real()) / (n_s.real() * n_s.real()) * sin(theta_0) * sin(theta_0)))) * std::norm(transmission_coeff_reverse);
        r_subs_exit = (n_s.real() * cos_theta_sub - n_exit_medium * std::cos(theta_0)) / (n_s.real() * cos_theta_sub + n_exit_medium * std::cos(theta_0));
    }
    if (strcmp(polarization, "p") == 0){
        transmissivity_no_backside_reverse = ((n_0.real() * std::sqrt(1 - sin(theta_0) * sin(theta_0))) / (n_s.real() * std::sqrt(1 - ((n_0.real() * n_0.real()) / (n_s.real() * n_s.real()) * sin(theta_0) * sin(theta_0))))) * std::norm(transmission_coeff_reverse);
        r_subs_exit = (n_exit_medium * cos_theta_sub - n_s.real() * std::cos(theta_0)) / (n_exit_medium * cos_theta_sub + n_s.real() * std::cos(theta_0));
    }
    double R_subs_exit = std::norm(r_subs_exit);
    return {reflectivity_coherent + (transmissivity_coherent * transmissivity_no_backside_reverse * R_subs_exit * std::pow(attenuation, 4)) / (1 - reflectivity_no_backside_reverse * R_subs_exit * std::pow(attenuation, 4)), (transmissivity_coherent * (1 - R_subs_exit) * std::pow(attenuation, 2)) / (1 - reflectivity_no_backside_reverse * R_subs_exit * std::pow(attenuation, 4))};
};
std::tuple<double, double> calculate_rt_unpolarized(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0, double n_exit_medium)
{
    auto[reflectivity_s, transmissivity_s] = calculate_rt(e_list_3x3, d_list, "s", wavelength, theta_0, phi_0, n_exit_medium);
    auto[reflectivity_p, transmissivity_p] = calculate_rt(e_list_3x3, d_list, "p", wavelength, theta_0, phi_0, n_exit_medium);
    return {0.5 * (reflectivity_s + reflectivity_p), 0.5 * (transmissivity_p + transmissivity_s)};
};