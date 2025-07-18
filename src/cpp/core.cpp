#include <complex>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace Eigen;
//Troparevsky et al. Optics Express Vol. 18, Issue 24, pp. 24715-24721 (2010)
//Formalism by Al-Ghezi et al. Optics Express 35770 (2020)
//Reflectivity and Transmissivity coefficients in "Generalized matrix method for calculation of Internal light energy flux in mixed coherent and incoherent multilayers" by Emanuele Centurioni (2005)
Matrix2cd coherent_block(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const char &polarization, const double &wavelength, const double &theta_0, const std::complex<double> &n_0)
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
        if (polarization == 's'){
            // Calculate kz value according to Berreman formalism - s polarization
            kz_i = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i](1, 1)/(wavelength*wavelength) - (kx * kx));
            kz_im1 = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[i-1](1, 1)/(wavelength*wavelength) - (kx * kx));
            pre_factor_dynamical = 0.5*(e_list_3x3[i](1, 1))/(e_list_3x3[i-1](1, 1));       
            // Anisotropic dynamical matrix
            dynamical(0,0) = std::complex<double>(1,0) + (kz_im1/kz_i);    
            dynamical(0,1) = std::complex<double>(1,0) - (kz_im1/kz_i);
        }
        if (polarization == 'p'){
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
    if (polarization == 's'){
    transmissivity_coherent = ((n_exit*std::sqrt(1-((n_0.real()*n_0.real())/(n_exit.real()*n_exit.real())*sin(theta_0)*sin(theta_0)))).real())/(std::sqrt(n_0.real()*(1-sin(theta_0)*sin(theta_0)))) * std::norm(transmission_coherent_coeff);
    }
    if (polarization == 'p'){
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

Matrix2cd incoherent_block(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const char &polarization, const double &wavelength, const double &theta_0, const double &n_0, const bool &closing_layer)
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
    if (polarization == 's')
    {
        kz_inc = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[1](1, 1).real()/(wavelength*wavelength) - (kx * kx));
        kz_exit = std::sqrt(4.0 * M_PI * M_PI * e_list_3x3[0](1, 1).real()/(wavelength*wavelength) - (kx * kx));

        r_forward = (kz_inc-kz_exit)/(kz_inc+kz_exit);
        t_forward =  2 * kz_inc / (kz_inc + kz_exit);
        t_backward = 2 * kz_exit / (kz_inc + kz_exit);
    }

    if (polarization == 'p')
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

std::tuple<double, double> calculate_rt(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const std::vector<bool> &incoherent, const char &polarization, const double &wavelength, const double &theta_0)
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
            general_Matrix = coherent_block(std::get<0>(calculation_blocks[i]), std::get<1>(calculation_blocks[i]), polarization, wavelength, theta_0, n_0.real()) * general_Matrix;
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
            general_Matrix = incoherent_block(std::get<0>(calculation_blocks[i]), std::get<1>(calculation_blocks[i]), polarization, wavelength, theta_0, n_0.real(), closing_layer) * general_Matrix;
        }
    }
    std::complex<double> transmission_coeff = pow(n_0,2)/pow(n_exit_medium,2) * std::complex<double>(1,0) / general_Matrix(0,0);
    double transmissivity;
    if (polarization == 's'){
        transmissivity = ((n_exit_medium*std::sqrt(1-((n_0.real()*n_0.real())/(n_exit_medium.real()*n_exit_medium.real())*sin(theta_0)*sin(theta_0)))).real())/(n_0.real()*std::sqrt((1-sin(theta_0)*sin(theta_0)))) * transmission_coeff.real();
    }
    if (polarization == 'p'){
        transmissivity = ((std::conj(n_exit_medium)*std::sqrt(1-((n_0.real()*n_0.real())/(n_exit_medium.real()*n_exit_medium.real())*sin(theta_0)*sin(theta_0)))).real())/(n_0.real()*std::sqrt((1-sin(theta_0)*sin(theta_0)))) * transmission_coeff.real();
    }
    return {((general_Matrix(1,0)/general_Matrix(0,0))).real(), transmissivity};      
};

std::tuple<double, double> calculate_rt_azim_pola(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const std::vector<bool> &incoherent, const double &wavelength, const double &theta_0, const double &phi_0, const double &s_polarization_percentage)
{
    const char temp_s = 's';
    const char temp_p = 'p';

    std::tuple<double, double> s_pol = calculate_rt(e_list_3x3, d_list, incoherent, temp_s, wavelength, theta_0);

    // Create a modified e_list_3x3 for p-polarization where the 0,0 diagonal term and 1,1 diagonal term are swapped - this is not
    // computationally efficient but necessary without other flags on the tensor.

    std::vector<Matrix3cd> e_list_3x3_swapped = e_list_3x3;
    for (auto &matrix : e_list_3x3_swapped) {
        std::swap(matrix(0, 0), matrix(1, 1));
    }
    // Calculate p-polarization with the modified e_list_3x3
    std::tuple<double, double> p_pol = calculate_rt(e_list_3x3_swapped, d_list, incoherent, temp_p, wavelength, theta_0);

    //Azimuthal angle 
    if (phi_0 == 0){
        return {s_polarization_percentage * std::get<0>(s_pol) +  (1-s_polarization_percentage) * std::get<0>(p_pol), s_polarization_percentage * std::get<1>(s_pol) + (1-s_polarization_percentage) * std::get<1>(p_pol)};
    }
    else{
        double cos_phi_squared = cos(phi_0) * cos(phi_0);
        double sin_phi_squared = sin(phi_0) * sin(phi_0);
        // The transmissivity and reflectivity are rotated according to the azimuthal angle
        // at theta 0 degrees, s- matches the y- and p- matches the x-axis
        return {(cos_phi_squared * (1-s_polarization_percentage) + sin_phi_squared * (s_polarization_percentage)) * std::get<0>(p_pol) + (cos_phi_squared * s_polarization_percentage + sin_phi_squared * (1-s_polarization_percentage)) * std::get<0>(s_pol) , (cos_phi_squared * (1-s_polarization_percentage) + sin_phi_squared * (s_polarization_percentage)) * std::get<1>(p_pol) + (cos_phi_squared * s_polarization_percentage + sin_phi_squared * (1-s_polarization_percentage)) * std::get<1>(s_pol)};
    }
};

double calculate_electric_field_at_z(const std::vector<Matrix3cd> &e_list_3x3, const std::vector<double> &d_list, const std::vector<bool> &incoherent, const char &polarization, const double &wavelength, const double &theta_0, const double &z)
{ 
    //We first need to find the layer where the height z is located
    double z_sum = 0;
    int layer = 0;
    while (z_sum < z){
        z_sum += d_list[layer];
        layer++;
    }
    //We now know the index of the last layer, we need to perform the calculation below up to there by creating a copy of the e_list_3x3 and d_list truncated up to that point
    std::vector<Matrix3cd> truncated_e_list_3x3(e_list_3x3.begin(), e_list_3x3.begin() + layer);
    std::vector<double> truncated_d_list(d_list.begin(), d_list.begin() + layer);
    std::vector<bool> truncated_incoherent(incoherent.begin(), incoherent.begin() + layer);
    //We need to replace the last layer with a layer of thickness z - z_sum
    truncated_d_list[layer-1] = z - z_sum;
    //Formalism from "Generalized matrix method for calculation of Internal light energy flux in mixed coherent and incoherent multilayers" by Emanuele Centurioni (2005)
    std::complex<double> n_exit_medium = std::sqrt(truncated_e_list_3x3[0](0, 0));
    std::complex<double> n_0 = std::sqrt(truncated_e_list_3x3[truncated_d_list.size()-1](0, 0));
    std::vector<Matrix3cd> sub_e_list;
    sub_e_list.push_back(truncated_e_list_3x3[0]);
    std::vector<double> sub_d_list;
    sub_d_list.push_back(truncated_d_list[0]);
    std::vector<std::tuple<std::vector<Matrix3cd>, std::vector<double>, bool>> calculation_blocks;
    bool currently_incoherent = incoherent[1];
    for (int i = 1; i < truncated_d_list.size(); i++) {
        if (currently_incoherent != incoherent[i]){
            if (!currently_incoherent){
                 sub_e_list.push_back(truncated_e_list_3x3[i]);
                 sub_d_list.push_back(0);
            }
            calculation_blocks.push_back({sub_e_list, sub_d_list, currently_incoherent});
            sub_e_list.clear();
            sub_d_list.clear();
            sub_e_list.push_back(truncated_e_list_3x3[i-1]);
            sub_d_list.push_back(0.0);
            sub_e_list.push_back(truncated_e_list_3x3[i]);
            sub_d_list.push_back(truncated_d_list[i]);
            currently_incoherent = !currently_incoherent;
        }
        else{
            sub_e_list.push_back(truncated_e_list_3x3[i]);
            sub_d_list.push_back(truncated_d_list[i]);
        }
    }
    calculation_blocks.push_back({sub_e_list, sub_d_list, currently_incoherent});
    Matrix2cd general_Matrix = Matrix2cd::Identity();
    for (int i = 0; i < calculation_blocks.size(); i++){
        if (std::get<2>(calculation_blocks[i]) == false)
        {
            general_Matrix = coherent_block(std::get<0>(calculation_blocks[i]), std::get<1>(calculation_blocks[i]), polarization, wavelength, theta_0, n_0.real()) * general_Matrix;
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
            general_Matrix = incoherent_block(std::get<0>(calculation_blocks[i]), std::get<1>(calculation_blocks[i]), polarization, wavelength, theta_0, n_0.real(), closing_layer) * general_Matrix;
        }
    }
    // we need to multiply this by [1,general_Matrix(1,0)/general_Matrix(0,0)] and extract the two components - as defined in Eq. 43 of the paper
    std::complex<double> electric_field_plus = 1; // with reference to Eo+R
    std::complex<double> electric_field_minus = general_Matrix(1,0)/general_Matrix(0,0);
    // with Eq. 52 and 53 we determine Ef and Eb for backward and forward propagating waves (all components)
    // depending on the polarisation
    std::complex<double> electric_field_forward_x;
    std::complex<double> electric_field_forward_y;
    std::complex<double> electric_field_forward_z;

    std::complex<double> electric_field_backward_x;
    std::complex<double> electric_field_backward_y;
    std::complex<double> electric_field_backward_z;

    //for s-polarisation
    if (polarization == 's'){
        electric_field_forward_x = electric_field_plus;
        electric_field_backward_x = electric_field_minus;
        electric_field_forward_y = 0;
        electric_field_backward_y = 0;
        electric_field_forward_z = 0;
        electric_field_backward_z = 0;
    }
    //for p-polarisation
    if (polarization == 'p'){
        electric_field_forward_x = 0;
        electric_field_backward_x = 0;
        electric_field_forward_y = cos(theta_0) * electric_field_plus;
        electric_field_backward_y = -sin(theta_0) * electric_field_minus;
        electric_field_forward_z = cos(theta_0) * electric_field_plus;
        electric_field_backward_z = sin(theta_0) * electric_field_minus;
    }
    //we now need to calculate the electric field at the height z
    //we sum the components of the forward and backward propagating waves in the different direction
    //components are in phase already
    //the norm of the electric field is the square of the sum of the components
    //we might need to double check the choice of x,y,z components vs the polarisation of wave that was chosen in the software
    std::complex<double> electric_field_x = electric_field_forward_x + electric_field_backward_x;
    std::complex<double> electric_field_y = electric_field_forward_y + electric_field_backward_y;
    std::complex<double> electric_field_z = electric_field_forward_z + electric_field_backward_z;
    
    double electric_field_norm = std::norm(electric_field_x) + std::norm(electric_field_y) + std::norm(electric_field_z);
    return {electric_field_norm};
};