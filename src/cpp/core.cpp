#include <complex>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace Eigen;
//Troparevsky et al. Optics Express Vol. 18, Issue 24, pp. 24715-24721 (2010)
//Formalism by Al-Ghezi et al. Optics Express 35770 (2020)
//Reflectivity and Transmissivity coefficients in "Generalized matrix method for calculation of Internal light energy flux in mixed coherent and incoherent multilayers" by Emanuele Centurioni (2005)
Matrix2cd coherent_block(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const char* &polarization, const double &wavelength, const double &theta_0, const double &phi_0, const std::complex<double> &n_0)
{ 
    std::complex<double> n_exit = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_inc = std::sqrt(e_list_3x3[e_list_3x3.size()-1](0, 0));
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
    std::complex<double> transmission_coherent_coeff =  pow(n_inc,2)/pow(n_exit,2) * std::complex<double>(1,0) / total_Matrix(0,0);    
    if (strcmp(polarization, "s") == 0){ 
    transmissivity_coherent = ((n_exit*std::sqrt(1-((n_0.real()*n_0.real())/(n_exit.real()*n_exit.real())*sin(theta_0)*sin(theta_0)))).real())/(std::sqrt(n_0.real()*(1-sin(theta_0)*sin(theta_0)))) * std::norm(transmission_coherent_coeff);
    }
    if (strcmp(polarization, "p") == 0){
    transmissivity_coherent = ((std::conj(n_exit)*std::sqrt(1-((n_0.real()*n_0.real())/(n_exit.real()*n_exit.real())*sin(theta_0)*sin(theta_0)))).real())/(std::sqrt(n_0.real()*(1-sin(theta_0)*sin(theta_0)))) * std::norm(transmission_coherent_coeff);
    }      
    std::complex<double> reflection_coeff_reverse = -(total_Matrix(0,1)/total_Matrix(0,0)); 
    std::complex<double> transmission_coeff_reverse =  pow(n_exit,2)/pow(n_0,2) * (total_Matrix.determinant() / total_Matrix(0,0));
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

Matrix2cd incoherent_block(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const char* &polarization, const double &wavelength, const double &theta_0, const double &phi_0, const double &n_0, const bool &closing_layer)
{
    Matrix2cd interface = Matrix2cd::Identity();
    Matrix2cd propa = Matrix2cd::Identity();
    std::complex<double> n_s = std::sqrt(e_list_3x3[d_list.size() - 1](0, 0));
    //Careful here with the angle that will come out in the first coherent layer (check for absorption effects in the future)
    double n_exit_medium = std::sqrt(e_list_3x3[0](0, 0)).real();
    //Propagation angle in the incoherent layer and the exit layer
    double cos_theta_sub = std::sqrt(1 - (n_0 * n_0 * sin(theta_0)*sin(theta_0))/(n_s.real()*n_s.real()));
    double cos_theta_exit = std::sqrt(1 - (n_0 * n_0 * sin(theta_0)*sin(theta_0))/(n_exit_medium*n_exit_medium));
    double r_forward;
    double t_forward;
    double r_backward;
    double t_backward;
    double kx = 2.0 * M_PI * n_0 / wavelength * sin(theta_0);
    double kz_inc;
    double kz_exit;
    double k_inc;
    double k_exit;
    if (strcmp(polarization, "s") == 0)
    {
        kz_inc = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[1](1, 1).real()/(wavelength*wavelength) - (kx * kx));
        kz_exit = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[0](1, 1).real()/(wavelength*wavelength) - (kx * kx));

        r_forward = (kz_inc-kz_exit)/(kz_inc+kz_exit);
        t_forward =  2 * kz_inc / (kz_inc + kz_exit);
        t_backward = 2 * kz_exit / (kz_inc + kz_exit);
    }

    if (strcmp(polarization, "p") == 0)
    {
        k_inc = std::sqrt(e_list_3x3[1](1, 1).real()) * 2 * M_PI / wavelength;
        k_exit = std::sqrt(e_list_3x3[0](1, 1).real()) * 2 * M_PI / wavelength;
        r_forward = (k_exit * cos_theta_sub - k_inc * cos_theta_exit)/(k_exit * cos_theta_sub + k_inc * cos_theta_exit);
        t_forward = 2 * k_inc * cos_theta_sub / (k_inc * cos_theta_exit + k_exit * cos_theta_sub);
        t_backward = 2 * k_exit * cos_theta_exit / (k_inc * cos_theta_exit + k_exit * cos_theta_sub);
    }
    r_backward = -r_forward;
    std::complex<double> pre_factor_substrate;
    if(closing_layer){
        pre_factor_substrate = 1/std::norm(t_forward);
        interface(0,1) = -std::norm(r_backward);
        interface(1,0) = std::norm(r_forward);
        interface(1,1) = std::norm(t_forward) * std::norm(t_backward) - std::norm(r_forward) * std::norm(r_backward);
    }
    else{
        pre_factor_substrate = 1;
        interface = Matrix2cd::Identity();
    }   
    //Substrate Layer Propagation Matrix
    propa(0,0) = std::exp(std::complex<double> (0,-4 * M_PI * n_s.real() * cos_theta_sub * d_list[0] / wavelength));
    propa(1,1) = std::exp(std::complex<double> (0,4 * M_PI * n_s.real() * cos_theta_sub * d_list[0] / wavelength));
    return propa * pre_factor_substrate * interface;
};

std::tuple<double, double> calculate_rt(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const std::vector<bool> &incoherent, const char* &polarization, const double &wavelength, const double &theta_0, const double &phi_0)
{ 
    //Formalism from "Generalized matrix method for calculation of Internal light energy flux in mixed coherent and incoherent multilayers" by Emanuele Centurioni (2005)
    std::complex<double> n_exit_medium = std::sqrt(e_list_3x3[0](0, 0));
    std::complex<double> n_0 = std::sqrt(e_list_3x3[d_list.size()-1](0, 0));
    std::vector<Matrix3cd> sub_e_list;
    sub_e_list.push_back(e_list_3x3[0]);
    std::vector<double> sub_d_list;
    sub_d_list.push_back(d_list[0]);
    std::vector<std::tuple<std::vector<Matrix3cd>, std::vector<double>, bool>> calculation_blocks;
    bool currently_incoherent = incoherent[1];
    for (int i = 1; i < d_list.size(); i++) {
        if (currently_incoherent != incoherent[i]){
            if (!currently_incoherent){
                 sub_e_list.push_back(e_list_3x3[i]);
                 sub_d_list.push_back(0);
            }
            calculation_blocks.push_back({sub_e_list, sub_d_list, currently_incoherent});
            sub_e_list.clear();
            sub_d_list.clear();
            sub_e_list.push_back(e_list_3x3[i-1]);
            sub_d_list.push_back(0.0);
            sub_e_list.push_back(e_list_3x3[i]);
            sub_d_list.push_back(d_list[i]);
            currently_incoherent = !currently_incoherent;
        }
        else{
            sub_e_list.push_back(e_list_3x3[i]);
            sub_d_list.push_back(d_list[i]);
        }
    }
    calculation_blocks.push_back({sub_e_list, sub_d_list, currently_incoherent});
    Matrix2cd general_Matrix = Matrix2cd::Identity();
    for (int i = 0; i < calculation_blocks.size(); i++){
        if (std::get<2>(calculation_blocks[i]) == false)
        {
            general_Matrix =  coherent_block(std::get<0>(calculation_blocks[i]), std::get<1>(calculation_blocks[i]), polarization, wavelength, theta_0, phi_0, n_0.real()) * general_Matrix;
        }
        if (std::get<2>(calculation_blocks[i]) == true)
        {
            bool closing_layer;
            if (i == calculation_blocks.size()-1 || i == 0){
                closing_layer = true;
            }
            else{
                closing_layer = false;
            }
            general_Matrix = incoherent_block(std::get<0>(calculation_blocks[i]), std::get<1>(calculation_blocks[i]), polarization, wavelength, theta_0, phi_0, n_0.real(), closing_layer) * general_Matrix;
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

std::tuple<double, double> calculate_rt_unpolarized(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const std::vector<bool> &incoherent, const double &wavelength, const double &theta_0, const double &phi_0)
{
    const char* temp_s = "s";
    const char* temp_p = "p";
    auto[reflectivity_s, transmissivity_s] = calculate_rt(e_list_3x3, d_list, incoherent, temp_s, wavelength, theta_0, phi_0);
    auto[reflectivity_p, transmissivity_p] = calculate_rt(e_list_3x3, d_list, incoherent, temp_p, wavelength, theta_0, phi_0);
    return {0.5 * (reflectivity_s + reflectivity_p), 0.5 * (transmissivity_p + transmissivity_s)};
};