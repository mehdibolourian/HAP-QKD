#%% 
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import gamma, hyp2f1
import pandas as pd
from matplotlib import pyplot as plt
from scipy.integrate import dblquad, quad
from scipy.linalg import cholesky, inv
from scipy.special import factorial as fac
import multiprocess as mlp
import os
import functools as fnct
import time

# np.seterr(all='raise')

def get_zernikes_index_range(n) -> list:
    """
    Returns a list of Zernike indexes

    Parameters
    ----------
    n : int
        The maximum Zernike radial index to be returned

    Returns
    -------
    list
        a list of indexes in the format [[n, m, j]], where
        n is the radial, m is the azimuth and j is the incremental
    """
    out = []
    j = 0
    for n in range(n + 1):
        for m in np.arange(-n, n + 1, 1):
            if (n - np.abs(m)) % 2 == 0:
                out.append([n, m, j])
                j += 1         
    return out

def calculate_j_Noll(n, m):

    j = n * (n + 1) // 2 + abs(m)
    
    # Determine the appropriate case for the piecewise function
    case1 = (m > 0) & (n % 4 == 0) | (m > 0) & (n % 4 == 1)
    case2 = (m < 0) & (n % 4 == 2) | (m < 0) & (n % 4 == 3)
    case3 = (m >= 0) & (n % 4 == 2) | (m >= 0) & (n % 4 == 3)
    case4 = (m <= 0) & (n % 4 == 0) | (m <= 0) & (n % 4 == 1)
    
    # Apply conditions
    j += np.where(case1 | case2, 0, 0)
    j += np.where(case3 | case4, 1, 0)
    
    return j

def eta_0(alpha, beta) -> float:
    """
    Compute the smf coupling efficiency without turbulence, Eq. 19 Scriminich22

    Parameters
    ----------
    alpha : float
        alpha parameter
    beta : float
        beta parameter

    Returns
    -------
    float
        value of eta_0
    """

    eta_0 = 2*((np.exp(-beta**2)-np.exp(-beta**2*alpha**2))/(beta*np.sqrt(1-alpha**2)))**2
    return eta_0

def noll_to_zernike(j):
    ''' Convert from Noll index to Zernike
    
    This function converts from the Noll to the Zernike index notation using the 
    code defined in https://sourceforge.net/p/octave/optics/ci/default/tree/inst/zernike_noll_to_mn.m#l50
    '''
        
    n = int( np.fix(np.sqrt(2*j-1) + 0.5) - 1)
    s = np.mod(n,2)
    me = 2 * np.fix( (2*j + 1 - n*(n+1)) / 4)
    mo = 2 * np.fix( (2*(j+1) - n*(n+1)) / 4) - 1
    m = int( (mo*s + me*(1-s)) * (1 - 2*np.mod(j,2)) )
    return (m,n)

def zernike_N(m : int,n : int) -> float:
    ''' Calculate the normalization factor of the Zernike polynomial
    
    This function calculates the normalization factor N of the Zernike polynomials.
    '''
    
    return np.sqrt((2*(n+1)) / (1 + (m==0)))

def zernike_R(m: int, n: int, ro: float) -> float:
    ''' Calculate the Zernike radial function '''
    
    R = 0.
    
    for s in range(0, int( (n-np.abs(m))/2 ) + 1):
        num = ((-1)**s) * fac(n-s)
        den = fac(s) * fac(0.5*(n + m)-s) * fac(0.5*(n - m)-s)
        R += num/den * (ro**(n-2*s))
        
    return R

def compute_zernike(n_max):
    ''' Create function for Nmax Zernike polynomials
    
    This function creates a list of Nmax functions continuing for the first
    Nmax Zernike polynomials ordered according to the Noll indices '''
    
    Z = []
    
    for ii in range(n_max):
        (m,n) = noll_to_zernike(ii+1)
        # print(m,n)
        if m>= 0:
            jj = lambda theta, ro, m=m, n=n: zernike_N(m,n) * zernike_R(np.abs(m),n,ro) * np.cos(m*theta)
            Z.append(jj)
        else:
            jj = lambda theta,ro, m=m, n=n: - zernike_N(m,n) * zernike_R(np.abs(m),n,ro) * np.sin(m*theta)
            Z.append(jj)

    return Z

def circular_pupil(r, obs):
    ''' Calculate the circular pupil function '''
    
    if (obs <= r) and (r <= 1):
        return 1.
    else:
        return 0.
    
def subfun (pupil_function,i,j,Z,pupil_area,obs):
    fun = lambda theta, r: pupil_function(r) * Z[i](theta, r) * Z[j](theta, r)
    result = 1/pupil_area * dblquad(fun, obs, 1, 0, 2*np.pi)[0]
    return result

def calculate_CZZ(W_fbp, R_pupil, Z, n_max, obs):
    ''' Calculate the CZZ matrix
    
    This function calculates the CZZ matrix for the pupil P and the back-propagated
    fiber mode present in the Fibre function '''

    CZZ = np.zeros((n_max, n_max))
    
    backprop_fiber = lambda r : np.sqrt(2/(np.pi * (W_fbp/R_pupil)**2)) * \
            np.exp( -(r * R_pupil / W_fbp)**2 )
    
    pupil_function = lambda r: backprop_fiber(r) * circular_pupil(r, obs) * r

    pupil_area = 2*np.pi*quad(pupil_function, obs, 1)[0]
    
    for i in range(n_max):
        t0 = time.time()
        fun2 = lambda j: subfun(pupil_function,i,j,Z,pupil_area,obs) 
        queue = range(i+1)
        pool_threads = os.cpu_count() - 1
        pool = mlp.Pool(pool_threads)
        result = pool.map(fnct.partial(fun2), queue)
        pool.close()
        pool.join()
        t1 = time.time()
        total = t1-t0
        print('row %d done' % i)
        print('in %f seconds' % total)
        for j in range(i+1):
            CZZ[i, j] = result[j]
            if i != j:
                CZZ[j, i] = result[j]
        
        # for j in range(i + 1):
        #     fun = lambda theta, r: pupil_function(r) * Z[i](theta, r) * Z[j](theta, r)
        #     result = 1/pupil_area * dblquad(fun, obs, 1, 0, 2*np.pi)[0]
        #     CZZ[i, j] = result
        #     if i != j:
        #         CZZ[j, i] = result
        # print('row %d done' % i)

    return CZZ



def compute_conversion_matrix(CZZ):
    ''' Compute the conversion matrix
    
    This function computes the conversion matrix for the CZZ matrix present
    in self.CZZ '''
        
    Q = cholesky(CZZ)
    M = inv(Q.conj().T)
    return M

def beta(D_rx, MFD, lmbd, f) -> float:
    """
    Calculate beta parameter for the computation of eta_0, Eq. 20 Scriminich22

    Parameters
    ----------
    D_rx : float
        size of the receiver aperture at the lens before the fiber
    MFD : float
        mode field diameter of the fiber 
    lmbd : float
        wavelength of the beam
    f : float
        focal length before the fiber

    Returns
    -------
    float
        value of beta
    """

    beta = np.pi*D_rx*MFD/4/lmbd/f
    return beta

def beta_opt(alpha) -> float:
    """
    Compute the optimal beta parameter given alpha

    Parameters
    ----------
    alpha : float
        alpha parameter

    Returns
    -------
    float
        value of beta
    """

    return 1.22*np.exp(-0.55*alpha) - 0.1*np.exp(-8*alpha)

def G(n, beta = 11/3) -> float:
    """
    Compute the geometrical factor to be used for b_n in Eq. 22 Scriminich22

    Parameters
    ----------
    n : int
        Zernike radial index
    beta : float
        beta parameter of the turbulence spectrum (default Kolmogorov = 11/3)

    Returns
    -------
    float
        value of G
    """
    return (n + 1)/np.pi*gamma((beta + 4)/2)*gamma(beta/2)*gamma((2*n + 2-beta)/2)*np.sin(np.pi*(beta - 2)/2)/gamma((2*n + 4 + beta)/2)

def bn2_zernike(D_rx, r_0, n) -> float:
    """
    Compute the Zernike coefficient of order n, Eq. 22 Scriminich22.

    Parameters
    ----------
    D_rx : float
        size of the receiver aperture
    r_0 : float
        value of fried parameter
    n : int or np.ndarray
        Zernike radial index

    Returns
    -------
    float
        value of bn2
    """

    return (D_rx/r_0)**(5/3)*G(n)

def bn2(D_rx, r_0, n, obs):
    """
    Compute the annular coefficient of order n [Dai and Mahajan 2007, eq. 39].

    Parameters
    ----------
    D_rx : float
        size of the receiver aperture
    r_0 : float
        value of fried parameter
    n : int or np.ndarray
        Zernike radial index
    obs : float
        Obstruction ratio of receiver aperture.

    Returns
    -------
    float
        value of bn2
    """

    pi = np.pi
    constant_term = 0.023 * (pi ** (8 / 3)) / (2 ** (5 / 3) * gamma(17 / 6))
    
    part1 = (n + 1) * gamma(n - 5 / 6) / ((1 - obs**2) * (1 - obs**(2 * (n + 1))))
    part1 *= (D_rx / r_0) ** (5 / 3)
    
    term1 = (1 + obs ** (2 * n + 17 / 3)) * gamma(14 / 3)
    term1 /= (gamma(17 / 6) * gamma(n + 23 / 6))
    
    term2 = (2 * obs ** (2 * (n + 1))) / fac(n + 1)
    term2 *= hyp2f1(n - 5 / 6, -11 / 6, n + 2, obs ** 2)
    
    part2 = term1 - term2
    
    result = constant_term * part1 * part2
    
    return result

def integrand(x, xi, bj2):
        kappa = 1
        integrand = (np.cos(np.sum((0.5*np.arctan(2*bj2*kappa*x))) - xi*x*kappa)/(np.prod(1 + (4*(x*kappa)**2*bj2**2))**0.25))/kappa
        return integrand

def eta_ao(bj2) -> float:
    """
    Compute the smf coupling efficiency with turbulence, Eq. 24 Scriminich22

    Parameters
    ----------
    bj2 : list
        list or numpy array of Zernike coefficients (without order 0)

    Returns
    -------
    float
        value of eta_ao
    """
    return np.prod(1/np.sqrt(1 + 2*np.array(bj2)))

def eta_s(scint_index):
    return (1 + scint_index)**(-1/4)

def compute_eta_xi_probability_distribution(xi, bj2) -> float:
    """
    Compute the probability distribution of xi, Eq. 33 Scriminich22

    Parameters
    ----------
    xi : float
        input parameter for the probability distribution
    bj2 : numpy.array
        Zernike coefficients squared (without order 0)
        
    Returns
    -------
    float
        value of p_xi(xi)
    """
    
    integrand = lambda x : (np.cos(np.sum((0.5*np.arctan(2*bj2*x))) - xi*x)/(np.prod(1 + (4*(x)**2*bj2**2))**0.25))
    integral = quad_vec(integrand, 0, np.inf)[0]

    return integral/np.pi

def compute_eta_smf_probability_distribution(eta_smf, eta_max, bj2) -> float:
    """
    Compute the probability distribution of eta_smf, Eq. 34 Scriminich22

    Parameters
    ----------
    eta_smf : float
        input parameter for the probability distribution
    eta_max : float
        maximum normalized coupled flux computed as eta_0*eta_S
    bj2 : numpy.array
        Zernike coefficients squared (without order 0)
        
    Returns
    -------
    float
        value of p_smf(eta_smf)
    """

    return compute_eta_xi_probability_distribution(np.log(eta_max/eta_smf), bj2)/eta_smf

# %%
if __name__ == "__main__":
    # print(beta_opt(0.4195))
    # Zernike indices look-up table
    array_of_zernike_index = get_zernikes_index_range(30) 
    lut_zernike_index_pd = pd.DataFrame(array_of_zernike_index[1:], columns = ["n", "m", "j"])
    n = np.array(lut_zernike_index_pd["n"].values)
    j = np.array(lut_zernike_index_pd["j"].values)
    print("j_max", np.max(j))
    D_rx = 0.4
    r0 = 0.05
    obs = 0
    bj2 = bn2(D_rx, r0, n, obs)
    # bj2 = bj2 / (np.sum(bj2) / (D_rx/r0)**(5/3))
    n_max = 6
    residual_phase_max = 3.7313031113440047
    residual_phase = np.sum(bj2) - np.sum(bj2[n <= n_max])
    print("Residual phase", residual_phase)
    print("Fitting error", (residual_phase_max - residual_phase)*100/ residual_phase_max)
    bj2[n <= n_max] = 0
    bj_wvln = np.sqrt(bj2)/(2*np.pi)
    print("Max variance in wavelength units", np.max(bj_wvln))
    print("Rayleigh criterion", 0.05)

    eta_max = 0.81
    eta = 10**(-np.linspace(8, 0, 1000))
    # eta = np.linspace(1e-8, 1, 1000)
    smf_pdf = compute_eta_smf_probability_distribution(eta, eta_max, bj2)
    smf_pdf = smf_pdf/np.sum(smf_pdf)
    mean_eta = np.sum(eta*smf_pdf)
    mean_eta_theor = eta_ao(bj2)*eta_max
    print("mean eta ", mean_eta)
    print("mean eta theor ", mean_eta_theor)
    print("relative error", 100*np.abs(mean_eta_theor - mean_eta)/mean_eta)

    plt.figure()
    plt.plot(eta, smf_pdf)
    plt.show()

    # # Test conversion matrix
    # n_max = 50
    # Z = compute_zernike(n_max)
    # CZZ = calculate_CZZ(0.41/2, 0.41/2, Z, n_max, 0)
    # plt.figure()
    # plt.imshow(CZZ)
    # plt.colorbar()
    # plt.show()
    # print("CZZ?", np.all(np.linalg.eigvals(CZZ) > 0))
    # M = compute_conversion_matrix(CZZ)
    # plt.figure()
    # plt.imshow(M)
    # plt.colorbar()
    # plt.show()

# %%