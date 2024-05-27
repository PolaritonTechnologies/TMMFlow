#include <complex>
#include "Eigen/Dense"
#include <iostream>
using namespace Eigen;
//Troparevsky et al. Optics Express Vol. 18, Issue 24, pp. 24715-24721 (2010)
//Formalism by Al-Ghezi et al. Optics Express 35770 (2020)
//Reflectivity and Transmissivity coefficients in "Generalized matrix method for calculation of Internal light energy flux in mixed coherent and incoherent multilayers" by Emanuele Centurioni (2005)
Matrix2cd coherent_block(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, const char* polarization, double wavelength, double theta_0, double phi_0)
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
    double transmissivity_coherent;  
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
    std::complex<double> reflection_coherent_coeff = (total_Matrix(1,0)/total_Matrix(0,0)); 
    std::complex<double> transmission_coherent_coeff = pow(n_0,2)/pow(n_s,2) * std::complex<double>(1,0) / total_Matrix(0,0);    
    if (strcmp(polarization, "s") == 0){ 
    transmissivity_coherent = ((n_s*std::sqrt(1-((n_0.real()*n_0.real())/(n_s.real()*n_s.real())*sin(theta_0)*sin(theta_0)))).real())/(std::sqrt(n_0.real()*(1-sin(theta_0)*sin(theta_0)))) * std::norm(transmission_coherent_coeff);
    }
    if (strcmp(polarization, "p") == 0){
    transmissivity_coherent = ((std::conj(n_s)*std::sqrt(1-((n_0.real()*n_0.real())/(n_s.real()*n_s.real())*sin(theta_0)*sin(theta_0)))).real())/(std::sqrt(n_0.real()*(1-sin(theta_0)*sin(theta_0)))) * std::norm(transmission_coherent_coeff);
    }      
 
    std::complex<double> reflection_coeff_reverse = -(total_Matrix(0,1)/total_Matrix(0,0)); 
    std::complex<double> transmission_coeff_reverse = pow(n_s,2)/pow(n_0,2) * (total_Matrix.determinant() / total_Matrix(0,0));
    double reflectivity_no_backside_reverse = std::norm(reflection_coeff_reverse);

    //Coherent Packet as Incoherent Interface Equivalent
    Matrix2cd coherent_as_incoherent = Matrix2cd::Identity();
    std::complex<double> pre_factor = 1/std::norm(transmission_coherent_coeff);
    coherent_as_incoherent(0,1) = -reflectivity_no_backside_reverse;
    coherent_as_incoherent(1,0) = std::norm(reflection_coherent_coeff);
    coherent_as_incoherent(1,1) = std::norm(transmission_coherent_coeff) * std::norm(transmission_coeff_reverse) - std::norm(reflection_coherent_coeff) * reflectivity_no_backside_reverse;
    coherent_as_incoherent = pre_factor * coherent_as_incoherent;

    return coherent_as_incoherent;
};

Matrix2cd incoherent_block(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, const char* polarization, double wavelength, double theta_0, double phi_0, double n_0)
{
    Matrix2cd interface = Matrix2cd::Identity();
    Matrix2cd propa = Matrix2cd::Identity();

    std::complex<double> n_s = std::sqrt(e_list_3x3[0](0, 0));

    //Careful here with the angle that will come out in the first coherent layer (check for absorption effects in the future)
    double n_exit_medium = std::sqrt(e_list_3x3[d_list.size() - 1](0, 0)).real();

    //Propagation angle in the incoherent layer
    double cos_theta_sub = std::sqrt(1 - (n_0 * n_0 * sin(theta_0)*sin(theta_0))/(n_s.real()*n_s.real()));

    double r_forward;
    double t_forward;
    double r_backward;
    double t_backward;

    if (strcmp(polarization, "s") == 0)
    {
        //for light moving forward
        r_forward = (n_s.real() * cos_theta_sub - n_exit_medium * std::cos(theta_0)) / (n_s.real() * cos_theta_sub + n_exit_medium * std::cos(theta_0));
        t_forward = 2 * n_s.real() * cos_theta_sub / (n_s.real() * cos_theta_sub + n_exit_medium * std::cos(theta_0));
        
        //for light moving backward
        r_backward = (n_exit_medium * std::cos(theta_0) - n_s.real() * cos_theta_sub) / (n_exit_medium * std::cos(theta_0) + n_s.real() * cos_theta_sub);
        t_backward = 2 * n_exit_medium * std::cos(theta_0) / (n_exit_medium * std::cos(theta_0) + n_s.real() * cos_theta_sub);
    }

    if (strcmp(polarization, "p") == 0)
    {
        //for light moving forward
        r_forward = (n_exit_medium * cos_theta_sub - n_s.real() * std::cos(theta_0)) / (n_exit_medium * cos_theta_sub + n_s.real() * std::cos(theta_0));
        t_forward = 2 * n_s.real() * cos_theta_sub / (n_exit_medium * cos_theta_sub + n_s.real() * std::cos(theta_0));
        
        //for light moving backward
        r_backward = (n_s.real() * std::cos(theta_0) - n_exit_medium * cos_theta_sub) / (n_s.real() * std::cos(theta_0) + n_exit_medium * cos_theta_sub);
        t_backward = 2 * n_exit_medium * cos_theta_sub / (n_s.real() * std::cos(theta_0) + n_exit_medium * cos_theta_sub);
    }

    std::complex<double> pre_factor_substrate = 1/std::norm(t_forward);

    interface(0,1) = -std::norm(r_backward);
    interface(1,0) = std::norm(r_forward);
    interface(1,1) = std::norm(t_forward) * std::norm(t_backward) - std::norm(r_forward) * std::norm(r_backward);
   
    //Substrate Layer Propagation Matrix
    propa(0,0) = std::exp(std::complex<double> (0,-4 * M_PI * n_s.real() * cos_theta_sub * d_list[0] / wavelength));
    propa(1,1) = std::exp(std::complex<double> (0,4 * M_PI * n_s.real() * cos_theta_sub * d_list[0] / wavelength));
    
    Matrix2cd incoherent_matrix = propa * pre_factor_substrate * interface;

    return incoherent_matrix;
};

std::tuple<double, double> calculate_rt(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, std::vector<bool> incoherent, const char* polarization, double wavelength, double theta_0, double phi_0)
{
    //Formalism from "Generalized matrix method for calculation of Internal light energy flux in mixed coherent and incoherent multilayers" by Emanuele Centurioni (2005)
    std::complex<double> n_exit_medium = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_0 = std::sqrt(e_list_3x3[d_list.size()-1](0, 0));
    //Loop through the e_list to find the first incoherent layer, assemble the sub_e_list with coherent layers
    int i = 0;
    std::vector<Matrix3cd> sub_e_list;
    std::vector<double> sub_d_list;
    std::vector<std::tuple<std::vector<Matrix3cd>, std::vector<double>, bool>> calculation_block;
    bool currently_incoherent = false;
    //There should exist at least one incoherent layer - a substrate
    while (i < d_list.size()){
        i++;
        sub_d_list.push_back(d_list[i]);
        sub_e_list.push_back(e_list_3x3[i]);
        if (currently_incoherent != incoherent[i]){
            calculation_block.push_back({sub_e_list, sub_d_list, currently_incoherent});
            sub_e_list.clear();
            sub_d_list.clear();
            sub_e_list.push_back(e_list_3x3[i]);
            sub_d_list.push_back(d_list[i]);
            currently_incoherent = !currently_incoherent;
        }
    }
    Matrix2cd general_Matrix = Matrix2cd::Identity();
    for (int i = 0; i < calculation_block.size(); i++){
        // the block is coherent
        if (std::get<2>(calculation_block[i]) == false)
        {
            general_Matrix = incoherent_block(std::get<0>(calculation_block[i]), std::get<1>(calculation_block[i]), polarization, wavelength, theta_0, phi_0, n_0.real()) * general_Matrix;
        }
        //the block is incoherent
        if (std::get<2>(calculation_block[i]) == true)
        {
            general_Matrix = coherent_block(std::get<0>(calculation_block[i]), std::get<1>(calculation_block[i]), polarization, wavelength, theta_0, phi_0) * general_Matrix;
        }
    }
    std::complex<double> transmission_coeff = pow(n_0,2)/pow(n_exit_medium,2) * std::complex<double>(1,0) / general_Matrix(0,0);
    double transmissivity;
    if (strcmp(polarization, "s") == 0){
        transmissivity = ((n_exit_medium*std::sqrt(1-((n_0.real()*n_0.real())/(n_exit_medium.real()*n_exit_medium.real())*sin(theta_0)*sin(theta_0)))).real())/(n_0.real()*std::sqrt((1-sin(theta_0)*sin(theta_0)))) * transmission_coeff.real();
    }
    if (strcmp(polarization, "p") == 0){
        transmissivity = ((std::conj(n_exit_medium)*std::sqrt(1-((n_0.real()*n_0.real())/(n_exit_medium.real()*n_exit_medium.real())*sin(theta_0)*sin(theta_0)))).real())/(n_0.real()*std::sqrt((1-sin(theta_0)*sin(theta_0)))) * transmission_coeff.real();
    }
    return {((general_Matrix(1,0)/general_Matrix(0,0))).real(), transmissivity};      
};

std::tuple<double, double> calculate_rt_unpolarized(std::vector<Matrix3cd> e_list_3x3, std::vector<double> d_list, std::vector<bool> incoherent, double wavelength, double theta_0, double phi_0)
{
    auto[reflectivity_s, transmissivity_s] = calculate_rt(e_list_3x3, d_list, incoherent, "s", wavelength, theta_0, phi_0);
    auto[reflectivity_p, transmissivity_p] = calculate_rt(e_list_3x3, d_list, incoherent, "p", wavelength, theta_0, phi_0);
    return {0.5 * (reflectivity_s + reflectivity_p), 0.5 * (transmissivity_p + transmissivity_s)};
};
