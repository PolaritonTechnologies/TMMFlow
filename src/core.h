#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Dense"
#include <chrono>
#include <thread>

using namespace Eigen;

std::tuple<double, double> calculate_rt_s_coherent(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0)
{
    //Troparevsky et al. Optics Express Vol. 18, Issue 24, pp. 24715-24721 (2010)
    //Formalism by Al-Ghezi et al. Optics Express 35770 (2020)
    std::complex<double> n_s = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_0 = std::sqrt(e_list_3x3[e_list_3x3.size() - 1](0, 0));

    // Wavevector modulus and in plane components
    std::complex<double> k0 = 2.0 * M_PI / wavelength;
    std::complex<double> kx = k0 * n_0.real() * sin(theta_0);
    std::complex<double> kz_i = 0;
    std::complex<double> kz_im1 = 0;
    std::complex<double> pre_factor_dynamic = 0;

    // Here total_Matrix_s represents the electric displacement
    Matrix2cd total_Matrix_s = Matrix2cd::Identity();
    Matrix2cd propa_and_dynamical = Matrix2cd::Zero();
    Matrix2cd propa = Matrix2cd::Zero();

    for (int i = 1; i < d_list.size() ; ++i)
    {
        propa_and_dynamical = Matrix2cd::Zero();
        propa = Matrix2cd::Zero();
        e_list_3x3[i](1, 1) = std::conj(e_list_3x3[i](1, 1));
        e_list_3x3[i](2, 2) = std::conj(e_list_3x3[i](2, 2));
        // Calculate kz value according to Berreman formalism - s polarization
        kz_i = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i](1, 1)/(wavelength*wavelength) - (kx * kx));
        kz_im1 = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i-1](1, 1)/(wavelength*wavelength) - (kx * kx));
        pre_factor_dynamic = 0.5*(e_list_3x3[i](1, 1))/(e_list_3x3[i-1](1, 1));       
        // Anisotropic dynamical matrix
        propa_and_dynamical(0,0) = std::complex<double>(1,0) + (kz_im1/kz_i);
        propa_and_dynamical(0,1) = std::complex<double>(1,0) - (kz_im1/kz_i);
        propa_and_dynamical(1,0) = propa_and_dynamical(0,1);
        propa_and_dynamical(1,1) = propa_and_dynamical(0,0);
        if (i<d_list.size()){
        // Propagation matrix
        propa(0,0) = std::exp(std::complex<double>(0,1) * kz_i * d_list[i]);
        propa(1,1) = std::exp(std::complex<double>(0,-1) * kz_i * d_list[i]);        
        propa_and_dynamical = propa * pre_factor_dynamic * propa_and_dynamical;
        }
        else{
        propa_and_dynamical = pre_factor_dynamic * propa_and_dynamical;    
        }
        total_Matrix_s = propa_and_dynamical * total_Matrix_s;
    }
    //std::cout << "Total matrix s: " << total_Matrix_s << std::endl;
    std::complex<double> reflection_s_coeff = 0;
    std::complex<double> transmission_s_coeff = 0;
    reflection_s_coeff =  (total_Matrix_s(1,0)/total_Matrix_s(0,0)); 
    transmission_s_coeff = (n_0 * n_0)/(n_s * n_s * total_Matrix_s(0,0));
    double reflectivity_s_coherent, transmissivity_s_coherent;
    reflectivity_s_coherent = (reflection_s_coeff * std::conj(reflection_s_coeff)).real();
    transmissivity_s_coherent = (n_s.real() / n_0.real()) * (transmission_s_coeff * std::conj(transmission_s_coeff)).real();
    
    //std::cout << "Reflectivity s coherent: " << reflectivity_s_coherent << std::endl;
    //std::cout << "Transmissivity s coherent: " << transmissivity_s_coherent << std::endl;

    return {reflectivity_s_coherent, transmissivity_s_coherent};
};

std::tuple<double, double> calculate_rt_s(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0, double n_exit_medium)
{
    double reflectivity_s, transmissivity_s;

    auto [reflectivity_s_no_backside, transmissivity_s_no_backside] = calculate_rt_s_coherent(e_list_3x3, d_list, wavelength, theta_0, phi_0);

    //std::cout << "Reflectivity s no correction: " << reflectivity_s_no_backside << std::endl;  
    //std::cout << "Transmissivity s no correction: " << transmissivity_s_no_backside << std::endl;

    // Incoherent layer

    // Flip the coherent structure to obtain the reverse quantities
    std::vector<Matrix3cd> e_list_3x3_reverse = e_list_3x3;
    std::reverse(e_list_3x3_reverse.begin(), e_list_3x3_reverse.end());
    std::vector<double> d_list_reverse = d_list;
    std::reverse(d_list_reverse.begin(), d_list_reverse.end());

    // Convert the wavelength of interest based on the Incoherent layer refractive index
    std::complex<double> n_s = std::sqrt(e_list_3x3[0](0, 0));
    double wavelength_reverse = wavelength / n_s.real();
    double theta_0_reverse = asin(sin(theta_0) / n_s.real());

    //std::cout << "Wavelength reverse: " << wavelength_reverse << std::endl;
    //std::cout << "Theta 0 reverse: " << theta_0_reverse << std::endl;

    auto [reflectivity_s_no_backside_reverse, transmissivity_s_no_backside_reverse] = calculate_rt_s_coherent(e_list_3x3_reverse, d_list_reverse, wavelength_reverse, theta_0_reverse, phi_0);
    //std::cout << "Transmissivity s coherent reverse: " << transmissivity_s_no_backside_reverse << std::endl;
    //std::cout << "Reflectivity s coherent reverse: " << reflectivity_s_no_backside_reverse << std::endl;

    // propagation angle in the substrate
    double cos_theta_sub = std::sqrt(1 - (sin(theta_0)*sin(theta_0)/ n_s.real()*n_s.real()));
    //std::cout << "Cos theta sub: " << cos_theta_sub << std::endl;

    // Attenuation through one pass in the incoherent layer
    double attenuation = std::exp((d_list[0] / cos_theta_sub) * sqrt(4 * M_PI * e_list_3x3[0](1, 1) / (wavelength * wavelength)).imag());
    //std::cout << "Attenuation: " << attenuation << std::endl;

    // coefficient from reflectivity and transmissivity to the exit medium from the substrate
    double r_s_subs_exit = (n_s.real() * cos_theta_sub - n_exit_medium * std::cos(theta_0)) / (n_s.real() * cos_theta_sub + n_exit_medium * std::cos(theta_0));
    double R_s_subs_exit = r_s_subs_exit * r_s_subs_exit;
    double T_s_subs_exit = 1 - R_s_subs_exit;
    
    //std::cout << "R_s_subs_exit: " << R_s_subs_exit << std::endl;
    //std::cout << "T_s_subs_exit: " << T_s_subs_exit << std::endl;

    reflectivity_s = reflectivity_s_no_backside + (transmissivity_s_no_backside * transmissivity_s_no_backside_reverse * R_s_subs_exit * std::pow(attenuation, 4)) / (1 - reflectivity_s_no_backside_reverse * R_s_subs_exit * std::pow(attenuation, 4));
    transmissivity_s = (transmissivity_s_no_backside * T_s_subs_exit * std::pow(attenuation, 2)) / (1 - reflectivity_s_no_backside_reverse * R_s_subs_exit * std::pow(attenuation, 4));
    
    //std::cout << "Reflectivity s: " << reflectivity_s << std::endl;
    //std::cout << "Transmissivity s: " << transmissivity_s << std::endl;

    return {reflectivity_s, transmissivity_s};
};


std::tuple<double, double> calculate_rt_p_coherent(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0)
{

    std::complex<double> n_s = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_0 = std::sqrt(e_list_3x3[e_list_3x3.size() - 1](0, 0));

    // Wavevector modulus and in plane components
    std::complex<double> k0 = 2.0 * M_PI / wavelength;
    std::complex<double> kx = k0 * n_0.real() * sin(theta_0);
    std::complex<double> kz_i = 0;
    std::complex<double> kz_im1 = 0;
    std::complex<double> pre_factor_dynamic = 0;
    std::complex<double> quotient_dynamical = 0;

    // Here total_Matrix_p represents the electric displacement
    Matrix2cd total_Matrix_p = Matrix2cd::Identity();
    Matrix2cd propa_and_dynamical = Matrix2cd::Zero();
    Matrix2cd propa = Matrix2cd::Zero();

    for (int i = 1; i < d_list.size() ; ++i)
    {
        propa_and_dynamical = Matrix2cd::Zero();
        propa = Matrix2cd::Zero();
        e_list_3x3[i](1, 1) = std::conj(e_list_3x3[i](1, 1));
        e_list_3x3[i](2, 2) = std::conj(e_list_3x3[i](2, 2));
        // Calculate kz value according to Berreman formalism
        kz_i = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i](1, 1)/(wavelength*wavelength) - ((kx * kx) * e_list_3x3[i](1, 1))/(e_list_3x3[i](2, 2)));
        kz_im1 = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i-1](1, 1)/(wavelength*wavelength) - ((kx * kx) * e_list_3x3[i-1](1, 1))/(e_list_3x3[i-1](2, 2)));
        pre_factor_dynamic = 0;
        quotient_dynamical = (e_list_3x3[i-1](1, 1) * kz_i)/(e_list_3x3[i](1, 1) * kz_im1);
        pre_factor_dynamic = (0.5/quotient_dynamical)*sqrt((kz_i*kz_i + kx*kx)/(kz_im1*kz_im1 + kx*kx));     
        // Anisotropic dynamical matrix
        propa_and_dynamical(0,0) = std::complex<double>(1,0) + quotient_dynamical;
        propa_and_dynamical(0,1) = std::complex<double>(1,0) - quotient_dynamical;
        propa_and_dynamical(1,0) = propa_and_dynamical(0,1);
        propa_and_dynamical(1,1) = propa_and_dynamical(0,0);
        if (i<d_list.size()){
        // Propagation matrix
        propa(0,0) = std::exp(std::complex<double>(0,1) * kz_i * d_list[i]);
        propa(1,1) = std::exp(std::complex<double>(0,-1) * kz_i * d_list[i]);        
        propa_and_dynamical = propa * pre_factor_dynamic * propa_and_dynamical;
        }
        else{
        propa_and_dynamical = pre_factor_dynamic * propa_and_dynamical;    
        }
        total_Matrix_p = propa_and_dynamical * total_Matrix_p;
    }
    std::complex<double> reflection_p_coeff = 0;
    std::complex<double> transmission_p_coeff = 0;
    reflection_p_coeff =  (total_Matrix_p(1,0)/total_Matrix_p(0,0));
    transmission_p_coeff = (n_0 * n_0)/(n_s * n_s * total_Matrix_p(0,0));
    double reflectivity_p_coherent, transmissivity_p_coherent;
    reflectivity_p_coherent = (reflection_p_coeff * std::conj(reflection_p_coeff)).real();
    transmissivity_p_coherent = (n_s.real() / n_0.real()) * (transmission_p_coeff * std::conj(transmission_p_coeff)).real();
    
    return {reflectivity_p_coherent, transmissivity_p_coherent};
};

std::tuple<double, double> calculate_rt_p(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0, double n_exit_medium)
{
    double reflectivity_p, transmissivity_p;

     auto [reflectivity_p_no_backside, transmissivity_p_no_backside] = calculate_rt_p_coherent(e_list_3x3, d_list, wavelength, theta_0, phi_0);

    // Incoherent layer

    // Flip the coherent structure to obtain the reverse quantities
    std::vector<Matrix3cd> e_list_3x3_reverse = e_list_3x3;
    std::reverse(e_list_3x3_reverse.begin(), e_list_3x3_reverse.end());
    std::vector<double> d_list_reverse = d_list;
    std::reverse(d_list_reverse.begin(), d_list_reverse.end());

    // Convert the wavelength of interest based on the Incoherent layer refractive index
    std::complex<double> n_s = std::sqrt(e_list_3x3[0](0, 0));
    double wavelength_reverse = wavelength / n_s.real();
    double theta_0_reverse = asin(sin(theta_0) / n_s.real());

    auto [reflectivity_p_no_backside_reverse, transmissivity_p_no_backside_reverse] = calculate_rt_p_coherent(e_list_3x3_reverse, d_list_reverse, wavelength_reverse, theta_0_reverse, phi_0);

    // propagation angle in the substrate
    double cos_theta_sub = std::sqrt(1 - (sin(theta_0)*sin(theta_0)/ n_s.real()*n_s.real()));

    // Attenuation through one pass in the incoherent layer
    double attenuation = std::exp((d_list[0] / cos_theta_sub) * sqrt(4 * M_PI * e_list_3x3[0](1, 1) / (wavelength * wavelength)).imag());

    // coefficient from reflectivity and transmissivity to the exit medium from the substrate
    double r_p_subs_exit = (n_s.real() * cos_theta_sub - n_exit_medium * std::cos(theta_0)) / (n_s.real() * cos_theta_sub + n_exit_medium * std::cos(theta_0));
    double R_p_subs_exit = r_p_subs_exit * r_p_subs_exit;
    double T_p_subs_exit = 1 - R_p_subs_exit;

    reflectivity_p = reflectivity_p_no_backside + (transmissivity_p_no_backside * transmissivity_p_no_backside_reverse * R_p_subs_exit * std::pow(attenuation, 4)) / (1 - reflectivity_p_no_backside_reverse * R_p_subs_exit * std::pow(attenuation, 4));
    transmissivity_p = (transmissivity_p_no_backside * T_p_subs_exit * std::pow(attenuation, 2)) / (1 - reflectivity_p_no_backside_reverse * R_p_subs_exit * std::pow(attenuation, 4));

    return {reflectivity_p, transmissivity_p};
};   

std::tuple<double, double> calculate_rt(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, double wavelength, double theta_0, double phi_0, double n_exit_medium)
{
    double reflectivity_s, transmissivity_s, reflectivity_p, transmissivity_p;
    std::tie(reflectivity_s, transmissivity_s) = calculate_rt_s(e_list_3x3, d_list, wavelength, theta_0, phi_0, n_exit_medium);
    std::tie(reflectivity_p, transmissivity_p) = calculate_rt_p(e_list_3x3, d_list, wavelength, theta_0, phi_0, n_exit_medium);
    return {0.5 * (reflectivity_s + reflectivity_p), 0.5 * (transmissivity_p + transmissivity_s)};
};