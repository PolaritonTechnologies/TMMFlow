# #added modification for normal incidence case - Florian Le Roux
## based on Analysis of multilayer thin-film structures containing magneto
## optic and anisotropic media at oblique incidence using 2x2 matrices -
## Mansuripur JAP 67

""" core.py This module contains the core routine for the solution of the
transfer matrix problem of a multilayer with arbitrary dielectric tensors. For
information about each function see the docstrings.

    The most two important functions are:

    rt(...) -- Calculates reflection and transmission quantities for a
    multilayer structure with a general dielectric tensors

    mo_rt(...) -- Calculates reflection and transmission matrix for an
    isotropic multilayer with off diagonal components corresponding to:
    e_xy=-e_yx!=0        Polar Kerr or Faraday effect (mo_flag='pp')
    e_xz=-e_zx!=0        Transverse Kerr or Faraday effect (mo_flag='tt')
    e_yz=-e_zy!=0        Longitudinal  Kerr or Faraday effect (mo_flag='ll')
"""

import numpy as np
import scipy as sp
import scipy.linalg as sp_la


def nullspace(A, atol=1e-9):
    """Compute an approximate basis for the nullspace of A using the singular
    value decomposition of `A`.
    A = array([[-1.41965504e+00-3.98478163e+00j,  0.00000000e+00+0.00000000e+00j,
         4.05773850e-03+1.81320587e-03j],
       [ 0.00000000e+00+0.00000000e+00j, -3.55271368e-15+1.77635684e-15j,
         0.00000000e+00+0.00000000e+00j],
       [ 4.05773850e-03+1.81320587e-03j,  0.00000000e+00+0.00000000e+00j,
         2.98285002e+00+1.06542853e+00j]])

    atol = 1e-07

    return array([[ 0.00000000e+00+0.00000000e+00j],
       [-9.92738180e-01-1.20295079e-01j],
       [ 1.38777878e-17-3.33066907e-16j]])

    Parameters
    ----------
    'A' = ndarray;  A should be at most 2-D.  A 1-D array with length k will be
     treated  as a 2-D with shape (1, k)
    'atol' = float; The absolute tolerance for a zero singular value.  Singular
     values smaller than `atol` are considered to be zero.
    'rtol' = float; The relative tolerance.  Singular values less than
     rtol*smax are considered to be zero, where smax is the largest singular
     value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
    tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Returns
    -------
    'ns' = ndarray; If `A` is an array with shape (m, k), then `ns` will be an
    array with shape (k, n), where n is the estimated dimension of the
    nullspace of `A`.  The columns of `ns` are a basis for the
    nullspace; each element in numpy.dot(A, ns) will be approximately
    zero.
    """

    # singular value decomposition
    u, s, vh = sp_la.svd(A)
    # u = array([[-3.35604811e-01-9.42001519e-01j,  1.39329749e-03+7.80229655e-04j, -1.28415381e-19+2.03315489e-19j], 
    # [ 0.00000000e+00+0.00000000e+00j, 0.00000000e+00+0.00000000e+00j, 8.34134427e-01-5.51561200e-01j], 
    # [ 1.45419123e-03-6.59821242e-04j,  4.13931197e-01-9.10306769e-01j,  5.55112070e-17+2.22044461e-16j]])
    # vh = array([[ 9.99998725e-01+0.00000000e+00j, 0.00000000e+00+0.00000000e+00j, 1.33519208e-04+1.59129156e-03j],
    # [-1.59688328e-03+0.00000000e+00j,  0.00000000e+00+0.00000000e+00j, 8.36122711e-02+9.96497084e-01j], 
    # [ 0.00000000e+00+0.00000000e+00j, -9.92738180e-01-1.20295079e-01j, 1.38777878e-17-3.33066907e-16j]]) 
    # --> Slightly different and in different order in C++
    # s = array([4.23012351e+00, 3.16741723e+00, 3.97205465e-15])  --> that is correct
    null_mask = s <= atol
    # null_mask: array([False, False,  True]) --> that is correct
    null_space = sp.compress(null_mask, vh, axis=0)
    # Only get the last column in this case of the matrix vh
    # array([[ 0.00000000e+00+0.00000000e+00j, -9.92738180e-01-1.20295079e-01j, 1.38777878e-17-3.33066907e-16j]]) 
    return sp.transpose(null_space)


def kz_eigenvalues(k0, kx, ky, m_eps):
    """Calculates the kz from the characteristic equation using a companion
    matrix method.
    k0 = 0.015707963267948967 
    kx = (-2.741555386204003e-05+0j) 
    ky = (-0+0j)
    m_eps = [[-23.38459267+4.76594931j,   0.        +0.j        ,
          0.        +0.j        ],
       [  0.        +0.j        , -23.38459267+4.76594931j,
          0.        +0.j        ],
       [  0.        +0.j        ,   0.        +0.j        ,
        -23.38459267+4.76594931j]]

    Parameters
    ----------
    'ko'= vacuum wavevector
    'kx,ky'= in plane wavevector components
    'm_eps'= 3x3 complex dielectric tensor

    Returns
    -------
    'v_kz'= kz wavevector components
    """

    # output
    v_kz = np.zeros(4, dtype=np.complex128)

    # are we diagonal and isotropic?
    diag_flag = (m_eps == np.diag(np.diagonal(m_eps))).all()
    iso_flag = m_eps[0, 0] == m_eps[1, 1] == m_eps[2, 2]

    # diagonal isotropic material
    if diag_flag and iso_flag:
        kz = np.sqrt((k0**2) * m_eps[0, 0] - kx**2 - ky**2 + 0j)
        v_kz[0:2] = -kz
        v_kz[2:4] = kz

    # general material
    else:

        # coefficients for the quartic equation
        A = (kx / k0) * (
            ((m_eps[0, 2] + m_eps[2, 0]) / m_eps[2, 2])
            + (ky / k0) * ((m_eps[1, 2] + m_eps[2, 1]) / m_eps[2, 2])
        )

        B = (
            ((kx / k0) ** 2) * (1.0 + m_eps[0, 0] / m_eps[2, 2])
            + ((ky / k0) ** 2) * (1.0 + m_eps[1, 1] / m_eps[2, 2])
            + ((kx * ky) / ((k0) ** 2)) * (m_eps[0, 1] + m_eps[1, 0]) / m_eps[2, 2]
            + (
                (m_eps[0, 2] * m_eps[2, 0] + m_eps[1, 2] * m_eps[2, 1]) / m_eps[2, 2]
                - m_eps[0, 0]
                - m_eps[1, 1]
            )
        )

        C = (
            ((kx**2 + ky**2) / (k0**2))
            * (
                (kx / k0) * (m_eps[0, 2] + m_eps[2, 0]) / m_eps[2, 2]
                + (ky / k0) * (m_eps[1, 2] + m_eps[2, 1]) / m_eps[2, 2]
            )
            + (kx / k0)
            * (
                (m_eps[0, 1] * m_eps[1, 2] + m_eps[1, 0] * m_eps[2, 1]) / m_eps[2, 2]
                - (m_eps[1, 1] / m_eps[2, 2]) * (m_eps[0, 2] + m_eps[2, 0])
            )
            + (ky / k0)
            * (
                (m_eps[0, 1] * m_eps[2, 0] + m_eps[1, 0] * m_eps[0, 2]) / m_eps[2, 2]
                - (m_eps[0, 0] / m_eps[2, 2]) * (m_eps[1, 2] + m_eps[2, 1])
            )
        )

        D1 = ((kx**2 + ky**2) / (k0**2)) * (
            ((kx / k0) ** 2) * m_eps[0, 0] / m_eps[2, 2]
            + ((ky / k0) ** 2) * m_eps[1, 1] / m_eps[2, 2]
            + ((kx * ky) / (k0**2)) * (m_eps[0, 1] + m_eps[1, 0]) / m_eps[2, 2]
            - m_eps[0, 0] * m_eps[1, 1] / m_eps[2, 2]
        )

        D2 = ((kx / k0) ** 2) * (
            (m_eps[0, 1] * m_eps[1, 0] + m_eps[0, 2] * m_eps[2, 0]) / m_eps[2, 2]
            - m_eps[0, 0]
        )
        D3 = ((ky / k0) ** 2) * (
            (m_eps[0, 1] * m_eps[1, 0] + m_eps[1, 2] * m_eps[2, 1]) / m_eps[2, 2]
            - m_eps[1, 1]
        )
        D4 = ((kx * ky) / (k0**2)) * (
            (m_eps[0, 2] * m_eps[2, 1] + m_eps[2, 0] * m_eps[1, 2]) / m_eps[2, 2]
            - m_eps[0, 1]
            - m_eps[1, 0]
        )
        D5 = (
            m_eps[0, 0] * m_eps[1, 1]
            + (
                m_eps[0, 1] * m_eps[1, 2] * m_eps[2, 0]
                + m_eps[1, 0] * m_eps[2, 1] * m_eps[0, 2]
            )
            / m_eps[2, 2]
            - m_eps[0, 1] * m_eps[1, 0]
            - (m_eps[0, 0] / m_eps[2, 2]) * m_eps[1, 2] * m_eps[2, 1]
            - (m_eps[1, 1] / m_eps[2, 2]) * m_eps[0, 2] * m_eps[2, 0]
        )
        D = D1 + D2 + D3 + D4 + D5

        # companion matrix
        m_comp = np.zeros((4, 4), dtype=np.complex128)
        m_comp[1, 0] = 1.0
        m_comp[2, 1] = 1.0
        m_comp[3, 2] = 1.0
        m_comp[0, 3] = -D
        m_comp[1, 3] = -C
        m_comp[2, 3] = -B
        m_comp[3, 3] = -A

        # eigenvalues
        v_kz = k0 * np.linalg.eigvals(m_comp)

    # output sorted by imaginary part
    return v_kz[np.argsort(np.imag(v_kz))]


def kz_eigenvectors(k0, kx, ky, v_kz, m_eps):
    """Calculates the kz field eigenvectors from the characteristic equation.

    k0 = 0.015707963267948967
    kx = (-2.741555386204003e-05+0j)
    ky = (-0+0j)
    v_kz = array([-0.01570794-0.j, -0.01570794-0.j,  0.01570794+0.j,  0.01570794+0.j])
    m_eps = array([[1.+0.j, 0.+0.j, 0.+0.j],
       [0.+0.j, 1.+0.j, 0.+0.j],
       [0.+0.j, 0.+0.j, 1.+0.j]])
    
    return parameters:
    v_e = array([[ 0.99999848+0.j,  0.        +0.j, -0.00174533+0.j],
       [ 0.        +0.j,  1.        +0.j,  0.        +0.j],
       [-0.99999848+0.j,  0.        +0.j, -0.00174533+0.j],
       [ 0.        +0.j,  1.        +0.j,  0.        +0.j]])
    
    v_kz = array([-0.01570794-0.j, -0.01570794-0.j,  0.01570794+0.j,  0.01570794+0.j]) 

    Parameters
    ----------
    'ko'= vacuum wavevector
    'kx,ky'= in plane wavevector components
    'v_kz'= off plane wavevector components
    'm_eps'= 3x3 complex dielectric tensor

    Returns
    -------
    'v_e'= kz field eigenvectors
    """

    # initializing vector and matrix
    v_e = np.zeros((4, 3), dtype=np.complex128)
    m_k = np.zeros_like(m_eps)
    m_char = np.zeros_like(m_eps)

    # are we diagonal and isotropic?
    diag_flag = (m_eps == np.diag(np.diagonal(m_eps))).all()
    iso_flag = m_eps[0, 0] == m_eps[1, 1] == m_eps[2, 2]

    # diagonal isotropic material
    if diag_flag and iso_flag:

        if kx == 0.0 and ky == 0.0:

            v_e[0, :] = np.array([1.0, 0.0, 0.0])
            v_e[1, :] = np.array([0.0, 1.0, 0.0])
            v_e[2, :] = np.array([1.0, 0.0, 0.0])
            v_e[3, :] = np.array([0.0, 1.0, 0.0])

        elif kx == 0.0:

            v_e[0, :] = np.array([1.0, 0.0, 0.0])
            v_e[1, :] = np.array([0.0, -v_kz[1], ky])
            v_e[2, :] = np.array([1.0, 0.0, 0.0])
            v_e[3, :] = np.array([0.0, -v_kz[3], ky])

        elif ky == 0.0:

            v_e[0, :] = np.array([-v_kz[1], 0.0, kx])
            v_e[1, :] = np.array([0.0, 1.0, 0.0])
            v_e[2, :] = np.array([-v_kz[3], 0.0, kx])
            v_e[3, :] = np.array([0.0, 1.0, 0.0])

        else:

            v_e[0, :] = np.array([-v_kz[1], 0.0, kx])
            v_e[1, :] = np.array([-ky, kx, 0.0])
            v_e[2, :] = np.array([-v_kz[3], 0.0, kx])
            v_e[3, :] = np.array([-ky, kx, 0.0])

    # general material
    else:

        # normal incidence id taken into account by relying
        # on continuity

        for m in range(4):

            # k matrix
            m_k[0, 0] = 0.0
            m_k[0, 1] = -v_kz[m]
            m_k[0, 2] = ky
            m_k[1, 0] = v_kz[m]
            m_k[1, 1] = 0.0
            m_k[1, 2] = -kx
            m_k[2, 0] = -ky
            m_k[2, 1] = kx
            m_k[2, 2] = 0.0

            # Characteristic matrix
            m_char = np.dot(m_k, m_k) / (k0**2)
            m_char = m_char + m_eps

            # Calculating the null space
            null_space = nullspace(m_char, atol=1e-7)

            v_e[m, :] = null_space[:, 0]

        # cleaning small elements from the eigenvectors
        for m in range(4):
            max_e = np.abs(v_e[m, :]).max()
            v_e_rel = np.abs(v_e[m, :]) / max_e
            v_e[m, v_e_rel < 1.0e-12] = 0.0

        # eigenvector swapping to get appropriate polarization states
        if np.abs(v_e[0, 0]) == 0.0:
            swap_e = v_e[0, :].copy()
            v_e[0, :] = v_e[1, :].copy()
            v_e[1, :] = swap_e.copy()
            swap_kz = v_kz[0].copy()
            v_kz[0] = v_kz[1].copy()
            v_kz[1] = swap_kz.copy()

        if np.abs(v_e[2, 0]) == 0.0:
            swap_e = v_e[2, :].copy()
            v_e[2, :] = v_e[3, :].copy()
            v_e[3, :] = swap_e.copy()
            swap_kz = v_kz[2].copy()
            v_kz[2] = v_kz[3].copy()
            v_kz[3] = swap_kz.copy()

    # normalizing eigenvectors
    for m in range(4):
        v_e[m, :] = v_e[m, :] / np.abs(np.sqrt(np.dot(v_e[m, :], np.conj(v_e[m, :]))))

    return v_e, v_kz


def m_abc(k0, kx, ky, v_kz, v_e, d):
    """Calculates layer by layer boundary and propagation matrixes
        to solve the transfer matrix problem
    
    k0 = 0.015707963267948967
    kx = (-2.741555386204003e-05+0j)
    ky = (-0+0j)
    v_kz = array([-0.01570794-0.j, -0.01570794-0.j,  0.01570794+0.j,  0.01570794+0.j])
    v_e =array([[ 0.99999848+0.j,  0.        +0.j, -0.00174533+0.j],
       [ 0.        +0.j,  1.        +0.j,  0.        +0.j],
       [-0.99999848+0.j,  0.        +0.j, -0.00174533+0.j],
       [ 0.        +0.j,  1.        +0.j,  0.        +0.j]]) 
    d = 0.0

    results:
    for v_a, v_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34
    (array([ 0.+0.j,  0.+0.j, -0.-0.j,  0.+0.j]), 
    array([-0.00174533+0.j,  0.        +0.j,  0.00174533-0.j,  0.        +0.j]), 
    array([[1.+0.j, 0.+0.j], [0.+0.j, 1.+0.j]]), 
    array([[ 1.+0.j,  0.+0.j], [-0.-0.j,  1.+0.j]]), 
    array([[ 0.        +0.j,  0.01570794+0.j], [-0.01570799+0.j,  0.        -0.j]]), 
    array([[ 0.        +0.j, -0.01570794+0.j], [ 0.01570799+0.j,  0.        +0.j]]), 
    array([[1.+0.j, 0.+0.j], [0.+0.j, 1.+0.j]]), 
    array([[1.+0.j, 0.+0.j], [0.+0.j, 1.+0.j]]))

    Parameters
    ----------
    'ko'= vacuum wavevector
    'kx,ky'= in plane wavevector components
    'v_kz'= off plane wavevector components
    'v_e'= kz field eigenvectors
    'd'= layer thickness in nm

    Returns
    -------
    'm_a12,m_a34,m_b12,m_b34,m_c12,m_c34'= boundary condition and propagation
                                           matrixes
    """

    # matrix allocation
    m_a12 = np.identity(2, dtype=np.complex128)
    m_a34 = np.identity(2, dtype=np.complex128)
    m_b12 = np.zeros((2, 2), dtype=np.complex128)
    m_b34 = np.zeros_like(m_b12)
    m_c12 = np.zeros((2, 2), dtype=np.complex128)
    m_c34 = np.zeros_like(m_c12)

    a1 = v_e[0, 1] / v_e[0, 0]
    a2 = v_e[1, 0] / v_e[1, 1]

    # a34 matrix

    a3 = v_e[2, 1] / v_e[2, 0]
    a4 = v_e[3, 0] / v_e[3, 1]

    m_a34[0, 1] = a4
    m_a34[1, 0] = a3

    # b12 matrix
    b1 = v_e[0, 2] / v_e[0, 0]

    b2 = v_e[1, 2] / v_e[1, 1]

    m_b12[0, 0] = -v_kz[0] * a1 + ky * b1
    m_b12[0, 1] = -v_kz[1] + ky * b2
    m_b12[1, 0] = v_kz[0] - kx * b1
    m_b12[1, 1] = v_kz[1] * a2 - kx * b2

    # b34 matrix
    b3 = v_e[2, 2] / v_e[2, 0]

    b4 = v_e[3, 2] / v_e[3, 1]

    m_b34[0, 0] = -v_kz[2] * a3 + ky * b3
    m_b34[0, 1] = -v_kz[3] + ky * b4
    m_b34[1, 0] = v_kz[2] - kx * b3
    m_b34[1, 1] = v_kz[3] * a4 - kx * b4

    # c12 matrix
    m_c12[0, 0] = np.exp(1j * v_kz[0] * d)
    m_c12[1, 1] = np.exp(1j * v_kz[1] * d)

    # c34 matrix
    m_c34[0, 0] = np.exp(1j * v_kz[2] * d)
    m_c34[1, 1] = np.exp(1j * v_kz[3] * d)

    # pure coefficient vectors
    v_a = np.array([a1, a2, a3, a4])
    v_b = np.array([b1, b2, b3, b4])

    return v_a, v_b, m_a12, m_a34, m_b12, m_b34, m_c12, m_c34


def rt(wl, theta_0, phi_0, e_list_3x3, d_list):
    """
    wl = 400
    theta_0 = 0.0017453292519943296
    phi_0 = 0
    e_list_3x3 = array([[[  1.        +0.j        ,   0.        +0.j        ,
           0.        +0.j        ],
        [  0.        +0.j        ,   1.        +0.j        ,
           0.        +0.j        ],
        [  0.        +0.j        ,   0.        +0.j        ,
           1.        +0.j        ]],

       [[-23.38459267+4.76594931j,   0.        +0.j        ,
           0.        +0.j        ],
        [  0.        +0.j        , -23.38459267+4.76594931j,
           0.        +0.j        ],
        [  0.        +0.j        ,   0.        +0.j        ,
         -23.38459267+4.76594931j]],

       [[  2.90627602+0.84588284j,   0.        +0.j        ,
           0.        +0.j        ],
        [  0.        +0.j        ,   4.3259341 +4.83066447j,
           0.        +0.j        ],
        [  0.        +0.j        ,   0.        +0.j        ,
           2.98285307+1.06542853j]],

       [[-23.38459267+4.76594931j,   0.        +0.j        ,
           0.        +0.j        ],
        [  0.        +0.j        , -23.38459267+4.76594931j,
           0.        +0.j        ],
        [  0.        +0.j        ,   0.        +0.j        ,
         -23.38459267+4.76594931j]],

       [[  2.25      +0.j        ,   0.        +0.j        ,
           0.        +0.j        ],
        [  0.        +0.j        ,   2.25      +0.j        ,
           0.        +0.j        ],
        [  0.        +0.j        ,   0.        +0.j        ,
           2.25      +0.j        ]]])
   
    d_list = array([  0.,  20., 105., 100.,   0.])

    results:
    # Does not agree
    m_r_ps =  array([[-0.86308751-0.3615758j ,  0.        +0.j        ],
       [ 0.        +0.j        , -0.86216249-0.37338276j]])

    # Does not agree either
    m_t_ps = array([[-5.08256779e-05+4.8402379e-05j,  0.00000000e+00+0.0000000e+00j],
       [ 0.00000000e+00+0.0000000e+00j, -1.51900814e-05-4.6948622e-06j]])

    # Matches
    m_Kn = array([[[ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          1.57079393e-02+0.j        ],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          1.57079393e-02+0.j        ],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -1.57079393e-02-0.j        ],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -1.57079393e-02-0.j        ]],

       [[ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          7.70111787e-03+0.07634936j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          7.70111787e-03+0.07634936j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -7.70111787e-03-0.07634936j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -7.70111787e-03-0.07634936j]],

       [[ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          2.70549844e-02+0.00385721j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          3.65196650e-02+0.01631886j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -2.70549844e-02-0.00385721j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -3.65196650e-02-0.01631886j]],

       [[ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          7.70111787e-03+0.07634936j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          7.70111787e-03+0.07634936j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -7.70111787e-03-0.07634936j],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -7.70111787e-03-0.07634936j]],

       [[ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          2.35619290e-02+0.j        ],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
          2.35619290e-02+0.j        ],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -2.35619290e-02-0.j        ],
        [ 2.74155539e-05-0.j        , -0.00000000e+00+0.j        ,
         -2.35619290e-02-0.j        ]]])


    # Does not match
    m_En =array([[[[-0.00000000e+00-0.00000000e+00j,
           9.99998477e-01-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
          -1.74532837e-03+0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 1.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -8.63086194e-01-3.61575249e-01j],
         [ 0.00000000e+00-0.00000000e+00j,
           0.00000000e+00-0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -1.50637111e-03-6.31068500e-04j]],

        [[-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [-8.62162486e-01-3.73382762e-01j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
           3.89229624e-02-7.45758588e-02j],
         [ 0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
           2.51134546e-05+1.65095847e-05j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 3.94677534e-02-7.57414448e-02j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
           6.05563053e-02+1.49306219e-02j],
         [ 0.00000000e+00-0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
           7.47851201e-06-2.09902444e-05j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 5.29118045e-02-1.23961774e-02j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
          -8.70853662e-02+8.94763510e-02j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00+0.00000000e+00j,
           9.89608868e-05-6.73391873e-05j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j],
         [-2.26668512e-02+1.20574744e-03j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
           1.06280715e-01-1.38016863e-02j],
         [ 0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
           1.02475079e-04+7.42274016e-06j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 9.82924551e-03+1.09584541e-02j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
          -1.88024759e-05+3.27105181e-05j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
          -1.09532516e-08-7.85641902e-09j]],

        [[-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [-8.54621629e-06-9.94922552e-08j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -3.20231676e-05+1.56918281e-05j],
         [ 0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
           4.42970202e-09+1.19456996e-08j]],

        [[-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [-6.64386515e-06-4.59536995e-06j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
          -5.08256435e-05+4.84023462e-05j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
           5.91383315e-08-5.63186967e-08j]],

        [[-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [-1.51900814e-05-4.69486220e-06j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00+0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]]]]) 

    # Does not match
    m_Hn = array([[[[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           1.00000000e+00-0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-9.99998477e-01-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [ 1.74532837e-03-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           8.63087508e-01+3.61575800e-01j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-8.62161173e-01-3.73382194e-01j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-1.50475664e-03-6.51675526e-04j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
           3.81562441e-01+1.52624840e-01j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-3.87494969e-01-1.54701403e-01j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [ 6.88841895e-05-1.32193692e-04j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           4.28821979e-02-3.01657003e-01j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00+0.00000000e+00j]],

        [[ 8.61932414e-02+2.51103086e-01j,
          -0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [ 9.23484734e-05-2.16354001e-05j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
          -1.71965300e-01+1.32727391e-01j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[ 5.39511218e-02+2.07451278e-02j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-3.95610984e-05+2.10442521e-06j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
          -1.86444404e-01-2.32635842e-03j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[ 1.14675124e-02+3.56889786e-02j,
          -0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [ 1.71552610e-05+1.91261008e-05j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
          -1.68209391e-04-7.53534516e-05j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[ 3.70635249e-06+4.15881005e-05j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-1.49159537e-08-1.73646655e-10j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           9.19708745e-05+1.47957020e-04j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[ 1.90787542e-05-3.45458119e-05j,
          -0.00000000e+00+0.00000000e+00j],
         [ 0.00000000e+00-0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-1.15957263e-08-8.02042952e-09j,
          -0.00000000e+00-0.00000000e+00j]]],


       [[[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
          -7.62385168e-05+7.26035685e-05j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[ 2.27851067e-05+7.04228854e-06j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-2.65116800e-08-8.19407618e-09j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]],

        [[-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j],
         [ 0.00000000e+00+0.00000000e+00j,
           0.00000000e+00+0.00000000e+00j],
         [-0.00000000e+00-0.00000000e+00j,
          -0.00000000e+00-0.00000000e+00j]]]]) 

    Calculates reflection matrix, transmission matrix and field amplitudes for a
       multilayer structure with a general dielectric tensors

    Parameters
    ----------
    'wl'= vacuum incident wavelength in nm
    'theta_0,phi_0'= polar and azimuth angles as defined in Mansuripur
                     JAP 67(10)
    'e_list_3x3'= [n_layer+2,3,3] numpy array: it contains n_layers+2 3x3
                  dielectric tensors:
    e_list_3x3[0]= 3x3 incident medium dielectric tensor: must be real,diagonal
                   and isotropic,
    e_list_3x3[n_layers+1] = 3x3 substrate dielectric tensor: must be real,
                             diagonal and isotropic,
    e_list_3x3[n]=3x3 dielectric tensor of the n_th layers: arbitrary
    'd_list'= n_layers+2 numpy array: contains layer thinknesses:
    d_list[0]=d_list[n_layers+1]=0: for the incident medium and substrate
    d_list[n]=d_n n_th layer thickness in nm

    Returns
    -------
    'a dictionary'= {'m_r_ps':m_r_ps,  #reflection matrix
                     'm_t_ps':m_t_ps,  #transmission matrix
                     'm_Kn':m_Kn, # wavevectors (n_layers,n_k,n_xyz)
                     'm_En':m_En, # electric field amplitudes (n_layers,n_k,n_xyz,n_pol)
                     'm_Hn':m_Hn, # magnetic field amplitudes (n_layers,n_k,n_xyz,n_pol)
                     'wl': wl,'theta_0': theta_0,'phi_0': phi_0,  #inputs
                     'e_list_3x3': e_list_3x3,'d_list': d_list}   #inputs
    """

    # incident medium check
    diag_flag = (e_list_3x3[0] == np.diag(np.diagonal(e_list_3x3[0]))).all()
    iso_flag = e_list_3x3[0, 0, 0] == e_list_3x3[0, 1, 1] == e_list_3x3[0, 2, 2]
    real_flag = np.abs(e_list_3x3[0, 0, 0]) == np.real(e_list_3x3[0, 0, 0])
    if not (diag_flag and iso_flag and real_flag):
        raise Exception("Incident medium must be real, diagonal and isotropic")

    n_0 = np.sqrt(e_list_3x3[0, 0, 0])

    # substrate check
    diag_flag = (e_list_3x3[-1] == np.diag(np.diagonal(e_list_3x3[-1]))).all()
    iso_flag = e_list_3x3[-1, 0, 0] == e_list_3x3[-1, 1, 1] == e_list_3x3[-1, 2, 2]
    real_flag = np.abs(e_list_3x3[-1, 0, 0]) == np.real(e_list_3x3[-1, 0, 0])
    if not (diag_flag and iso_flag and real_flag):
        raise Exception("Substrate must be real, diagonal and isotropic")

    n_s = np.sqrt(e_list_3x3[-1, 0, 0])

    # wavevector modulus and in plane components
    k0 = 2.0 * np.pi / wl
    kx = -k0 * n_0 * np.sin(theta_0) * np.cos(phi_0)
    ky = -k0 * n_0 * np.sin(theta_0) * np.sin(phi_0)

    # kz,v_e and boundary and propagation matrix for R and T
    m_a12 = np.zeros((len(e_list_3x3), 2, 2), dtype=np.complex128)
    m_a34 = np.zeros_like(m_a12)
    m_b12 = np.zeros_like(m_a12)
    m_b34 = np.zeros_like(m_a12)
    m_c12 = np.zeros_like(m_a12)
    m_c34 = np.zeros_like(m_a12)
    m_a = np.zeros((len(e_list_3x3), 4), dtype=np.complex128)
    m_b = np.zeros_like(m_a)
    m_Kn = np.zeros((len(e_list_3x3), 4, 3), dtype=np.complex128)

    #! This first iterates over all available layers and calculates the eigenvalues and vectors
    for n in range(len(e_list_3x3)):
        v_kz = kz_eigenvalues(k0, kx, ky, e_list_3x3[n])
        v_e, v_kz = kz_eigenvectors(k0, kx, ky, v_kz, e_list_3x3[n])

        # storing the wavevectors
        m_Kn[n, :, 0] = kx  # kx
        m_Kn[n, :, 1] = ky  # ky
        m_Kn[n, :, 2] = v_kz  # kz

        # print(n,v_e)
        m_a[n], m_b[n], m_a12[n], m_a34[n], m_b12[n], m_b34[n], m_c12[n], m_c34[n] = (
            m_abc(k0, kx, ky, v_kz, v_e, d_list[n])
        )

    # looping for R over the layers
    m_R = np.zeros_like(m_c12)

    #! Now the iteration goes backwards from the last layer to the first layer
    #! to calculate the reflection. m_R[n+1] is thus the value that was
    #! calculated in the iteration before
    for n in range(len(e_list_3x3) - 2, -1, -1):

        # building the first factor for the F_n+1 matrix-----
        f1 = np.dot(m_b12[n + 1], m_c12[n + 1]) + np.dot(
            np.dot(m_b34[n + 1], m_c34[n + 1]), m_R[n + 1]
        )

        # building the second factor for the F_n+1 matrix-----
        f2_inv = np.dot(m_a12[n + 1], m_c12[n + 1]) + np.dot(
            np.dot(m_a34[n + 1], m_c34[n + 1]), m_R[n + 1]
        )
        f2 = np.linalg.inv(f2_inv)

        # F_n+1 matrix-----
        f_np1 = np.dot(f1, f2)

        # R_n
        r1_inv = np.dot(f_np1, m_a34[n]) - m_b34[n]
        r1 = np.linalg.inv(r1_inv)
        r2 = m_b12[n] - np.dot(f_np1, m_a12[n])
        m_R[n] = np.dot(r1, r2)
        print(m_R[n])

    # rotating m_R to the s,p states
    p_inc = np.zeros((2, 2))
    p_inc[0, 0] = np.cos(theta_0) * np.cos(phi_0)
    p_inc[0, 1] = -np.sin(phi_0)
    p_inc[1, 0] = np.cos(theta_0) * np.sin(phi_0)
    p_inc[1, 1] = np.cos(phi_0)
    p_inc_inv = np.linalg.inv(p_inc)

    # Finally the  R matrix output...
    m_r_ps = np.dot(np.dot(p_inc_inv, m_R[0]), p_inc)

    #! Now again the iteration is done for the transmission (at least this
    #! should be merged with the previous loop)
    # looping for T over the layers
    m_Tn = np.zeros_like(m_c12)
    m_T = np.identity(2, dtype=np.complex128)
    for n in range(len(e_list_3x3) - 2, -1, -1):

        # building the first factor for the T_n
        f1_inv = np.dot(m_a12[n + 1], m_c12[n + 1]) + np.dot(
            np.dot(m_a34[n + 1], m_c34[n + 1]), m_R[n + 1]
        )
        f1 = np.linalg.inv(f1_inv)

        # building the second factor for the T_n
        f2 = m_a12[n] + np.dot(m_a34[n], m_R[n])

        # T_n
        m_Tn[n] = np.dot(f1, f2)

        # T
        m_T = np.dot(m_T, m_Tn[n])
        print(m_T)

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

    # initializing the fields
    m_En = np.zeros((len(e_list_3x3), 4, 3, 2), dtype=np.complex128)
    m_Hn = np.zeros_like(m_En)
    m_En[0, 1, 1, 0] = 1.0  # TE polarization
    m_En[0, 0, 0, 1] = -np.cos(theta_0)  # TM polarization

    # loop over all layers
    for n in range(len(e_list_3x3)):

        # forward electric and magnetic fields in the i_th layer

        # Ey1 Ez1
        m_En[n, 0, 1, 0] = m_a[n, 0] * m_En[n, 0, 0, 0]  # TE
        m_En[n, 0, 2, 0] = m_b[n, 0] * m_En[n, 0, 0, 0]
        m_En[n, 0, 1, 1] = m_a[n, 0] * m_En[n, 0, 0, 1]  # TM
        m_En[n, 0, 2, 1] = m_b[n, 0] * m_En[n, 0, 0, 1]

        # Ex2 Ez2
        m_En[n, 1, 0, 0] = m_a[n, 1] * m_En[n, 1, 1, 0]  # TE
        m_En[n, 1, 2, 0] = m_b[n, 1] * m_En[n, 1, 1, 0]
        m_En[n, 1, 0, 1] = m_a[n, 1] * m_En[n, 1, 1, 1]  # TM
        m_En[n, 1, 2, 1] = m_b[n, 1] * m_En[n, 1, 1, 1]

        # Hx1 Hy1 Hz1
        m_Hn[n, 0, 0, 0] = m_b12[n, 0, 0] * m_En[n, 0, 0, 0] / k0  # TE
        m_Hn[n, 0, 1, 0] = m_b12[n, 1, 0] * m_En[n, 0, 0, 0] / k0
        m_Hn[n, 0, 2, 0] = (-ky + kx * m_a[n, 0]) * m_En[n, 0, 0, 0] / k0
        m_Hn[n, 0, 0, 1] = m_b12[n, 0, 0] * m_En[n, 0, 0, 1] / k0  # TM
        m_Hn[n, 0, 1, 1] = m_b12[n, 1, 0] * m_En[n, 0, 0, 1] / k0
        m_Hn[n, 0, 2, 1] = (-ky + kx * m_a[n, 0]) * m_En[n, 0, 0, 1] / k0

        # Hx2 Hy2 Hz2
        m_Hn[n, 1, 0, 0] = m_b12[n, 0, 1] * m_En[n, 1, 1, 0] / k0  # TE
        m_Hn[n, 1, 1, 0] = m_b12[n, 1, 1] * m_En[n, 1, 1, 0] / k0
        m_Hn[n, 1, 2, 0] = (-ky * m_a[n, 1] + kx) * m_En[n, 1, 1, 0] / k0
        m_Hn[n, 1, 0, 1] = m_b12[n, 0, 1] * m_En[n, 1, 1, 1] / k0  # TM
        m_Hn[n, 1, 1, 1] = m_b12[n, 1, 1] * m_En[n, 1, 1, 1] / k0
        m_Hn[n, 1, 2, 1] = (-ky * m_a[n, 1] + kx) * m_En[n, 1, 1, 1] / k0

        # exiting one before the last, because then I have no backpropagation
        if n == len(e_list_3x3) - 1:
            break

        # backward electric and magnetic fields in the i_th layer

        # Ex3 Ey4
        m_En[n, 2, 0, 0] = (
            m_R[n, 0, 0] * m_En[n, 0, 0, 0] + m_R[n, 0, 1] * m_En[n, 1, 1, 0]
        )  # TE
        m_En[n, 3, 1, 0] = (
            m_R[n, 1, 0] * m_En[n, 0, 0, 0] + m_R[n, 1, 1] * m_En[n, 1, 1, 0]
        )
        m_En[n, 2, 0, 1] = (
            m_R[n, 0, 0] * m_En[n, 0, 0, 1] + m_R[n, 0, 1] * m_En[n, 1, 1, 1]
        )  # TM
        m_En[n, 3, 1, 1] = (
            m_R[n, 1, 0] * m_En[n, 0, 0, 1] + m_R[n, 1, 1] * m_En[n, 1, 1, 1]
        )

        # Ey3 Ez3
        m_En[n, 2, 1, 0] = m_a[n, 2] * m_En[n, 2, 0, 0]
        m_En[n, 2, 2, 0] = m_b[n, 2] * m_En[n, 2, 0, 0]
        m_En[n, 2, 1, 1] = m_a[n, 2] * m_En[n, 2, 0, 1]
        m_En[n, 2, 2, 1] = m_b[n, 2] * m_En[n, 2, 0, 1]

        # Ex4 Ez4
        m_En[n, 3, 0, 0] = m_a[n, 3] * m_En[n, 3, 1, 0]
        m_En[n, 3, 2, 0] = m_b[n, 3] * m_En[n, 3, 1, 0]
        m_En[n, 3, 0, 1] = m_a[n, 3] * m_En[n, 3, 1, 1]
        m_En[n, 3, 2, 1] = m_b[n, 3] * m_En[n, 3, 1, 1]

        # Hx3 Hy3 Hz3
        m_Hn[n, 2, 0, 0] = m_b34[n, 0, 0] * m_En[n, 2, 0, 0] / k0  # TE
        m_Hn[n, 2, 1, 0] = m_b34[n, 1, 0] * m_En[n, 2, 0, 0] / k0
        m_Hn[n, 2, 2, 0] = (-ky + kx * m_a[n, 2]) * m_En[n, 2, 0, 0] / k0
        m_Hn[n, 2, 0, 1] = m_b34[n, 0, 0] * m_En[n, 2, 0, 1] / k0  # TM
        m_Hn[n, 2, 1, 1] = m_b34[n, 1, 0] * m_En[n, 2, 0, 1] / k0
        m_Hn[n, 2, 2, 1] = (-ky + kx * m_a[n, 2]) * m_En[n, 2, 0, 1] / k0

        # Hx4 Hy4 Hz4
        m_Hn[n, 3, 0, 0] = m_b34[n, 0, 1] * m_En[n, 3, 1, 0] / k0  # TE
        m_Hn[n, 3, 1, 0] = m_b34[n, 1, 1] * m_En[n, 3, 1, 0] / k0
        m_Hn[n, 3, 2, 0] = (-ky * m_a[n, 3] + kx) * m_En[n, 3, 1, 0] / k0
        m_Hn[n, 3, 0, 1] = m_b34[n, 0, 1] * m_En[n, 3, 1, 1] / k0  # TM
        m_Hn[n, 3, 1, 1] = m_b34[n, 1, 1] * m_En[n, 3, 1, 1] / k0
        m_Hn[n, 3, 2, 1] = (-ky * m_a[n, 3] + kx) * m_En[n, 3, 1, 1] / k0

        # Ex1 Ey2 n_th+1 layer
        m_En[n + 1, 0, 0, 0] = (
            m_Tn[n, 0, 0] * m_En[n, 0, 0, 0] + m_Tn[n, 0, 1] * m_En[n, 1, 1, 0]
        )  # TE
        m_En[n + 1, 1, 1, 0] = (
            m_Tn[n, 1, 0] * m_En[n, 0, 0, 0] + m_Tn[n, 1, 1] * m_En[n, 1, 1, 0]
        )
        m_En[n + 1, 0, 0, 1] = (
            m_Tn[n, 0, 0] * m_En[n, 0, 0, 1] + m_Tn[n, 0, 1] * m_En[n, 1, 1, 1]
        )  # TM
        m_En[n + 1, 1, 1, 1] = (
            m_Tn[n, 1, 0] * m_En[n, 0, 0, 1] + m_Tn[n, 1, 1] * m_En[n, 1, 1, 1]
        )

    # flipping the x and z coordinates for the right reference frame
    m_Kn[:, :, 0] = -m_Kn[:, :, 0]
    m_Kn[:, :, 2] = -m_Kn[:, :, 2]
    m_En[:, :, 0] = -m_En[:, :, 0]
    m_En[:, :, 2] = -m_En[:, :, 2]
    m_Hn[:, :, 0] = -m_Hn[:, :, 0]
    m_Hn[:, :, 2] = -m_Hn[:, :, 2]

    # Output in a for of a dictionary-
    return {
        "m_r_ps": m_r_ps,
        "m_t_ps": m_t_ps,
        "m_Kn": m_Kn,
        "m_En": m_En,
        "m_Hn": m_Hn,
        "wl": wl,
        "theta_0": theta_0,
        "phi_0": phi_0,
        "e_list_3x3": e_list_3x3,
        "d_list": d_list,
    }


"""
# commented out because only needed for off-diagonal components 
def mo_rt(wl, theta_0, phi_0, e_list, e_list_off, d_list, mo_flag):
    # Calculates reflection and transmission matrix for an isotropic multilayer
    #    with off diagonal components corresponding to:
    #    e_xy=-e_yx!=0        Polar Kerr or Faraday effect (mo_flag='pp')
    #    e_xz=-e_zx!=0        Transverse Kerr or Faraday effect (mo_flag='tt')
    #    e_yz=-e_zy!=0        Longitudinal  Kerr or Faraday effect (mo_flag='ll')

    # Parameters
    # ----------
    # 'wl'= vacuum incident wavelength in nm
    # 'theta_0,phi_0'= polar and azimuth angles as defined in Mansuripur
    #                  JAP 67(10)
    # 'e_list' =
    #         [n_layer+2] numpy array: contains n_layers+2 dielectric constants:
    #         e_list[0] = incident medium dielectric constant: must be real
    #         e_list[n_layers+1]=substrate dielectric constant: must be real
    #         e_list[n]=dielectric tensor of the n_th layer: can be complex
    # 'e_list_off' =
    #         [n_layer+2] numpy array: it contains n_layers+2 off diagonal
    #                     dielectric constants:
    #         e_list_off[0] = must be 0: incident medium is diagonal
    #         e_list_off[n_layers+1] = must be 0: substrate is diagonal
    #         e_list_off[n] = off diagonal dielectric constant of the n_th layer:
    #                         can be complex
    # 'd_list' =
    #         n_layers+2 numpy array: contains layer thinknesses:
    #         d_list[0]=d_list[n_layers+1]=0: incident medium and substrate
    #         d_list[n]=d_n n_th layer thickness in nm
    # 'mo_flag' =
    #         'pp' e_xy=-e_yx!=0  Polar Kerr or Faraday effect
    #         'tt' e_xz=-e_zx!=0  Transverse Kerr or Faraday effect
    #         'll' e_yz=-e_zy!=0  Longitudinal  Kerr or Faraday effect

    # Returns
    # -------
    # 'a dictionary'= {'m_r_ps':m_r_ps,          #reflection matrix
    #                  'm_t_ps':m_t_ps,          #transmission matrix
    #                  'wl': wl,'theta_0': theta_0,'phi_0': phi_0,#inputs
    #                  'e_list': e_list,'e_list_off': e_list_off, #inputs
    #                  'e_list_3x3': e_list_3x3, #full dielectric tensor
    #                  'd_list': d_list}         #layer thicknesses

    # filling the dielectric tensor depending on mo_flag
    e_list_3x3 = np.zeros((len(e_list), 3, 3), dtype=np.complex128)
    e_list_3x3[:, 0, 0] = e_list
    e_list_3x3[:, 1, 1] = e_list
    e_list_3x3[:, 2, 2] = e_list
    if mo_flag == "pp":
        e_list_3x3[:, 0, 1] = e_list_off
        e_list_3x3[:, 1, 0] = -e_list_off
    elif mo_flag == "tt":
        e_list_3x3[:, 0, 2] = e_list_off
        e_list_3x3[:, 2, 0] = -e_list_off
    elif mo_flag == "ll":
        e_list_3x3[:, 1, 2] = e_list_off
        e_list_3x3[:, 2, 1] = -e_list_off
    else:
        raise Exception("mo_flag must be either 'pp', 'tt' or 'll'...")

    # computing reflection and transmission matrix
    rt_out = rt(wl, theta_0, phi_0, e_list_3x3, d_list)
    m_r_ps = rt_out["m_r_ps"]
    m_t_ps = rt_out["m_t_ps"]

    # Output in the for of a dictionary
    return {
        "m_r_ps": m_r_ps,
        "m_t_ps": m_t_ps,
        "wl": wl,
        "theta_0": theta_0,
        "phi_0": phi_0,
        "e_list": e_list,
        "e_list_off": e_list_off,
        "e_list_3x3": e_list_3x3,
        "d_list": d_list,
    }

"""


def azimuthalCalculation(e_s_y, e_s_x, cosphiSquared, sinphiSquared, sinphiS2D2):
    """
    e_s_y = (-1.5190081436892855e-05-4.6948622023882075e-06j)
    e_s_x = (-5.082514751739978e-05+4.840250961008275e-05j)
    cosphiSquared = 0.8830222215594891 
    sinphiSquared = 0.11697777844051097
    sinphiS2D2 = 0.3213938048432696
    """

    # see Experimental notes for calculation, values are passed to reduce calculation time

    e_s_s = e_s_y * cosphiSquared + e_s_x * sinphiSquared

    e_s_p = (e_s_y - e_s_x) * sinphiS2D2

    c_X = (e_s_x**2) * sinphiSquared

    c_Y = (e_s_y**2) * cosphiSquared

    contribX = np.real(np.abs(c_X))

    contribY = np.real(np.abs(c_Y))

    E_s = np.abs(e_s_s) ** 2 + np.abs(e_s_p) ** 2

    return e_s_s, e_s_p, E_s, contribX, contribY
