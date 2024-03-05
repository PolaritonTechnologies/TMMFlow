import numpy as np
import scipy as sp

from core import kz_eigenvalues, kz_eigenvectors, m_abc


def calculate_tr_per_layer(
    m_a12,
    m_a34,
    m_b12,
    m_b34,
    m_a12_np1,
    m_a34_np1,
    m_b12_np1,
    m_b34_np1,
    m_c12_np1,
    m_c34_np1,
    m_R_np1,
    m_T,
):
    """
    This function should allow the calculation of the reflection and
    transmission without any weird data structures
    """
    # looping for R over the layers
    m_R = np.zeros((2, 2), dtype=np.complex128)

    ## Used to be the for loop block that was iterating from the second last
    # layer to the first one (in our example 3 to 0).

    # building the first factor for the F_n+1 matrix-----
    f1 = np.dot(m_b12_np1, m_c12_np1) + np.dot(np.dot(m_b34_np1, m_c34_np1), m_R_np1)

    # building the second factor for the F_n+1 matrix-----
    f2_inv = np.dot(m_a12_np1, m_c12_np1) + np.dot(
        np.dot(m_a34_np1, m_c34_np1), m_R_np1
    )
    f2 = np.linalg.inv(f2_inv)

    # F_n+1 matrix-----
    f_np1 = np.dot(f1, f2)

    # R_n
    r1_inv = np.dot(f_np1, m_a34) - m_b34
    r1 = np.linalg.inv(r1_inv)
    r2 = m_b12 - np.dot(f_np1, m_a12)
    m_R = np.dot(r1, r2)

    # looping for T over the layers
    m_Tn = np.zeros_like(m_R)

    # for n in range(len(e_list_3x3) - 2, -1, -1):
    ## Again this was the former for loop block

    # building the first factor for the T_n
    f1_trans_inv = np.dot(m_a12_np1, m_c12_np1) + np.dot(
        np.dot(m_a34_np1, m_c34_np1), m_R_np1
    )
    f1_trans = np.linalg.inv(f1_trans_inv)

    # building the second factor for the T_n
    f2_trans = m_a12 + np.dot(m_a34, m_R)

    # T_n
    m_Tn = np.dot(f1_trans, f2_trans)

    # T
    m_T = np.dot(m_T, m_Tn)

    ## End of for loop block

    # In the next iteration m_a12 --> m_a12_np1, similarly m_R --> m_R_np1.
    # However, I think the m_a12 -> m_a12_np1 should be done via chunk_1
    # I am not yet 100% sure if we need the m_Tn (but I think so from what we need for chunk 3)
    return m_R, m_T


def calculate_tr(e_list_3x3, d_list, wavelength, theta_0, phi_0):
    # refractive indices
    n_0 = np.sqrt(e_list_3x3[0, 0, 0])
    n_s = np.sqrt(e_list_3x3[-1, 0, 0])

    # wavevector modulus and in plane components
    k0 = 2.0 * np.pi / wavelength
    kx = -k0 * n_0 * np.sin(theta_0) * np.cos(phi_0)
    ky = -k0 * n_0 * np.sin(theta_0) * np.sin(phi_0)

    # m_a_np1, m_b_np1, m_a12_np1, m_a34_np1, m_b12_np1, m_b34_np1, m_c12_np1, m_c34_np1, m_K_np1 = chunk_1(wavelength, theta_0, phi_0, e_list_3x3[-1], d_list[-1], n_0, n_s)
    v_kz = kz_eigenvalues(k0, kx, ky, e_list_3x3[-1])
    v_e, v_kz = kz_eigenvectors(k0, kx, ky, v_kz, e_list_3x3[-1])
    (
        m_a_np1,
        m_b_np1,
        m_a12_np1,
        m_a34_np1,
        m_b12_np1,
        m_b34_np1,
        m_c12_np1,
        m_c34_np1,
    ) = m_abc(k0, kx, ky, v_kz, v_e, d_list[-1])

    m_T = np.identity(2, dtype=np.complex128)
    m_R_np1 = np.zeros((2, 2), dtype=np.complex128)

    # Now iterate over all layers starting from the second last going backwards
    for i in range(len(d_list) - 2, -1, -1):
        v_kz = kz_eigenvalues(k0, kx, ky, e_list_3x3[i])
        v_e, v_kz = kz_eigenvectors(k0, kx, ky, v_kz, e_list_3x3[i])
        m_a, m_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34 = m_abc(
            k0, kx, ky, v_kz, v_e, d_list[i]
        )

        m_R, m_T = calculate_tr_per_layer(
            m_a12,
            m_a34,
            m_b12,
            m_b34,
            m_a12_np1,
            m_a34_np1,
            m_b12_np1,
            m_b34_np1,
            m_c12_np1,
            m_c34_np1,
            m_R_np1,
            m_T,
        )

        if i == 0:
            m_R_0 = m_R

        # In the next iteration m_a12 --> m_a12_np1, similarly m_R --> m_R_np1.
        (
            m_a12_np1,
            m_a34_np1,
            m_b12_np1,
            m_b34_np1,
            m_c12_np1,
            m_c34_np1,
        ) = (m_a12, m_a34, m_b12, m_b34, m_c12, m_c34)
        m_R_np1 = m_R

    ## This has to be calculated outside the loop
    # rotating m_R to the s,p states
    p_inc = np.zeros((2, 2))
    p_inc[0, 0] = np.cos(theta_0) * np.cos(phi_0)
    p_inc[0, 1] = -np.sin(phi_0)
    p_inc[1, 0] = np.cos(theta_0) * np.sin(phi_0)
    p_inc[1, 1] = np.cos(phi_0)
    p_inc_inv = np.linalg.inv(p_inc)

    # Finally the  R matrix output...
    m_r_ps = np.dot(np.dot(p_inc_inv, m_R_0), p_inc)

    # rotating m_T to the s,p states
    theta_s = sp.arcsin(np.real_if_close(np.sin(theta_0) * n_0 / n_s))
    p_sub = np.zeros((2, 2), dtype=np.complex128)
    p_sub[0, 0] = np.cos(theta_s) * np.cos(phi_0)
    p_sub[0, 1] = -np.sin(phi_0)
    p_sub[1, 0] = np.cos(theta_s) * np.sin(phi_0)
    p_sub[1, 1] = np.cos(phi_0)
    p_sub_inv = np.linalg.inv(p_sub)

    # Finally the  T matrix output...
    m_t_ps = np.dot(np.dot(p_sub_inv, m_T), p_inc)


    return m_r_ps, m_t_ps

## Select values for the calculation
wavelength = 400
theta_0 = 0.0017453292519943296
phi_0 = 0

e_list_3x3 = np.array(
    [
        [
            [1.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, 1.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j],
        ],
        [
            [-23.38459267 + 4.76594931j, 0.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, -23.38459267 + 4.76594931j, 0.0 + 0.0j],
            [0.0 + 0.0j, 0.0 + 0.0j, -23.38459267 + 4.76594931j],
        ],
        [
            [2.90627602 + 0.84588284j, 0.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, 4.3259341 + 4.83066447j, 0.0 + 0.0j],
            [0.0 + 0.0j, 0.0 + 0.0j, 2.98285307 + 1.06542853j],
        ],
        [
            [-23.38459267 + 4.76594931j, 0.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, -23.38459267 + 4.76594931j, 0.0 + 0.0j],
            [0.0 + 0.0j, 0.0 + 0.0j, -23.38459267 + 4.76594931j],
        ],
        [
            [2.25 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, 2.25 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, 0.0 + 0.0j, 2.25 + 0.0j],
        ],
    ]
)
d_list = np.array([0.0, 20.0, 105.0, 100.0, 0.0])


m_r_ps, m_t_ps = calculate_tr(e_list_3x3, d_list, wavelength, theta_0, phi_0)
print("Reflection matrix: ", m_r_ps)
print("Transmission matrix: ", m_t_ps)
