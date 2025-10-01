import numpy as np
from netsquid.components.models.qerrormodels import QuantumErrorModel
import netsquid.util.simtools as simtools
from scipy.special import i0, i1
from scipy.integrate import quad
from matplotlib import pyplot as plt
from scipy.special import erf
from netsquid.util.simlog import warn_deprecated
from balloon_qnet import cn2
from scipy.integrate import quad_vec
from balloon_qnet import smf_coupling as smf
import pandas as pd
from math import fsum
import random
import warnings
from scipy.special import gamma
from scipy.special import kv
from balloon_qnet import transmittance as transmittance
from scipy.linalg import inv

""" This file contains the main functions of the free-space loss model for downlink, horizontal link and uplink."""

#np.seterr(all='raise')

RE = 6371 # Earth's radius [km]
c = 299792458 

# Zernike indices look-up table
max_n_wfs = 150 # Maximum radial index n returned by WFS 
array_of_zernike_index = smf.get_zernikes_index_range(max_n_wfs)
lut_zernike_index_pd = pd.DataFrame(array_of_zernike_index[1:], columns = ["n", "m", "j"])
lut_zernike_index_pd["j_Noll"] = smf.calculate_j_Noll(lut_zernike_index_pd["n"], lut_zernike_index_pd["m"])

def compute_channel_length(ground_station_alt, aerial_platform_alt, zenith_angle):
    """Compute channel length that corresponds to a particular ground station altitude, aerial 
    platform altitude and zenith angle.

    ## Parameters
        `ground_station_alt` : float 
            Altitude of the ground station [km].
        `aerial_platform_alt` : float 
            Altitude of the aerial platform [km].
        `zenith_angle` : float
            Zenith angle of aerial platform [degrees].
    ## Returns
    `length` : float
        Length of the channel [km].
    """
    zenith_angle = np.deg2rad(zenith_angle)
    RA = RE + aerial_platform_alt
    RG = RE + ground_station_alt
    length = np.sqrt(RA**2 + RG**2*(np.cos(zenith_angle)**2 - 1)) - RG*np.cos(zenith_angle)
    return length

def compute_height_min_horiz (length, height):
    """Compute minimal height of a horizontal channel between two ballons at the same height.
    
    ## Parameters 
        `length` : float
            length of the horizontal channel [km]
        `height̀` : float
            height of the balloons [km]
    ## Returns 
    `hmin` : float
        Minimal height of the channel [km] 
    """
    RS = RE + height
    theta = np.arcsin((length)/(2*RS))
    hmin = np.cos(theta)*RS -RE
    return hmin
    
def sec(theta):
    """Compute secant of angle theta.

    ## Parameters
    `theta` : float
        Angle for which secant will be calculated [degrees].
    ## Returns
    `sec` : float
        Secant result.
    """
    theta = np.deg2rad(theta)
    sec = 1/np.cos(theta)
    return sec

def lognormal_pdf(eta, mu, sigma):
    """Compute lognormal distribution probability density function (PDF).

    ## Parameters
    `eta` : np.ndarray
        Input random variable values to calculate PDF for.
    `mu` : float
        Mean value of lognormal distribution.
    `sigma` : float
        Standard deviation of lognormal distribution.
    ## Returns
    `pdf` : np.ndarray
        PDF of lognormal distribution for values of eta.
    """

    pdf = np.exp(-(np.log(eta) + mu)**2/(2*sigma**2))/(eta*sigma*np.sqrt(2*np.pi))
    return pdf

def lognormal_cdf(eta, mu, sigma):
    """Compute lognormal distribution cumulative density function (CDF).

    ## Parameters
    `eta` : np.ndarray
        Input random variable values to calculate CDF for.
    `mu` : float
        Mean value of lognormal distribution.
    `sigma` : float
        Standard deviation of lognormal distribution.
    ## Returns
    `cdf` : np.ndarray
        CDF of lognormal distribution for values of eta.
    """
    cdf = (1 + erf((np.log(eta) + mu)/(sigma*np.sqrt(2))))/2
    return cdf

def truncated_lognormal_pdf(eta, mu, sigma):
    """Compute truncated lognormal distribution probability density function (PDF) according to [Vasylyev et al., 2018].

    ## Parameters
    `eta` : np.ndarray
        Input random variable values to calculate PDF for.
    `mu` : float
        Mean value of truncated lognormal distribution.
    `sigma` : float
        Standard deviation of truncated lognormal distribution.
    ## Returns
    `pdf` : np.ndarray
        PDF of truncated lognormal distribution for values of eta.
    """
    lognormal_cdf_dif = lognormal_cdf(1, mu, sigma)
    if np.size(eta) == 1:
        if eta < 0 or eta > 1:
            pdf = 0
        else: 
            pdf = lognormal_pdf(eta, mu, sigma)/lognormal_cdf_dif
    else:
        pdf = np.zeros(np.size(eta))
        eta_domain = (eta >= 0) & (eta <= 1)
        pdf = lognormal_pdf(eta[eta_domain], mu, sigma)/lognormal_cdf_dif
    return pdf

def lognegative_weibull_pdf(eta, eta_0, wandering_variance, R, l):
    """Compute log-negative Weiibull distribution probability density function (PDF) according to [Vasylyev et al., 2018].

    ## Parameters
    `eta` : np.ndarray
        Input random variable values to calculate PDF for.
    `eta_0` : float
        Maximal transmittance of the Gaussian beam.
    `wandering_variance` : float
        Wandering variance of the Gaussian beam.
    `R` : float
        Scale parameter of distribution.
    `l` : float
        Shape parameter of distribution.
    ## Returns
    `pdf` : np.ndarray
        PDF of log-negative Weibull distribution for values of eta.
    """
    if np.size(eta) == 1:
        if eta < 0 or eta > eta_0:
            pdf = 0
        else: 
            pdf = (R**2/(wandering_variance*eta*l))*((np.log(eta_0/eta))**(2/l - 1))*np.exp(-(R**2/(2*wandering_variance))*(np.log(eta_0/eta))**(2/l))
    else:
        pdf = np.zeros(np.size(eta))
        eta_domain = (eta >= 0) & (eta <= eta_0)
        pdf[eta_domain] = (R**2/(wandering_variance*eta[eta_domain]*l))*((np.log(eta_0/eta[eta_domain]))**(2/l - 1))*np.exp(-(R**2/(2*wandering_variance))*(np.log(eta_0/eta[eta_domain]))**(2/l))
    return pdf

class HorizontalChannel(QuantumErrorModel):
    """Model for photon loss on a horizontal free-space channel.

    Uses probability density of atmospheric transmittance (PDT) from [Vasylyev et al., 2018] to
    sample the loss probability of the photon.

    ## Parameters
    ----------
    `W0` : float
        Waist radius of the beam at the transmitter [m].
    `rx_aperture` : float
        Diameter of the receiving telescope [m].
    `obs_ratio` : float
        Obscuration ratio of the receiving telescope.
    `Cn2` : float
        Index of refraction structure constant [m**(-2/3)].
    `wavelength` : float
        Wavelength of the radiation [m].
    `pointing_error` : float
        Pointing error [rad].
    `tracking_efficiency` : float
        Efficiency of the coarse tracking mechanism.
    `Tatm` : float
        Atmospheric transmittance (square of the transmission coefficient).
    `rng` : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    """
    def __init__(self, W0, rx_aperture, obs_ratio, Cn2, wavelength, pointing_error = 0, tracking_efficiency = 0, Tatm = 1, rng = None):
        super().__init__()
        self.rng = rng if rng else simtools.get_random_state()
        self.W0 = W0
        self.rx_aperture = rx_aperture
        self.obs_ratio = obs_ratio
        self.Cn2 = Cn2
        self.wavelength = wavelength
        self.pointing_error = pointing_error
        self.tracking_efficiency = tracking_efficiency
        self.Tatm = Tatm
        self.required_properties = ['length']

    @property
    def rng(self):
        """ :obj:`~numpy.random.RandomState`: Random number generator."""
        return self.properties['rng']

    @rng.setter
    def rng(self, value):
        if not isinstance(value, np.random.RandomState):
            raise TypeError("{} is not a valid numpy RandomState".format(value))
        self.properties['rng'] = value
        
    @property
    def Tatm(self):
        """ :float: atmosphere transmittance. """
        return self.properties['Tatm']

    @Tatm.setter
    def Tatm(self, value):
        if (value < 0) or (value > 1):
            raise ValueError
        self.properties['Tatm'] = value

    @property
    def pointing_error(self):
        """ :float: pointing error variance. """
        return self.properties['pointing_error']

    @pointing_error.setter
    def pointing_error(self, value):
        if (value < 0):
            raise ValueError
        self.properties['pointing_error'] = value

    @property
    def tracking_efficiency(self):
        """ :float: efficiency of the coarse tracking mechanism. """
        return self.properties['tracking_efficiency']

    @tracking_efficiency.setter
    def tracking_efficiency(self, value):
        if (value < 0) or (value > 1):
            raise ValueError
        self.properties['tracking_efficiency'] = value

    @property
    def W0(self):
        """float: beam waist radius at the transmitter [m]."""
        return self.properties['W0']

    @W0.setter
    def W0(self, value):
        if value < 0:
            raise ValueError
        self.properties['W0'] = value

    @property
    def rx_aperture(self):
        """float: diameter of the receiving telescope [m]."""
        return self.properties['rx_aperture']

    @rx_aperture.setter
    def rx_aperture(self, value):
        if value < 0:
            raise ValueError
        self.properties['rx_aperture'] = value

    @property
    def obs_ratio(self):
        """float: obscuration ratio of the receiving telescope."""
        return self.properties['obs_ratio']

    @obs_ratio.setter
    def obs_ratio(self, value):
        if value < 0 or (value > 1):
            raise ValueError
        self.properties['obs_ratio'] = value

    @property
    def Cn2(self):
        """float: index of refraction structure constant [m**(-2/3)]."""
        return self.properties['Cn2']

    @Cn2.setter
    def Cn2(self, value):
        if value < 0:
            raise ValueError
        self.properties['Cn2'] = value

    @property
    def wavelength(self):
        """float: wavelength of the radiation [m]."""
        return self.properties['wavelength']

    @wavelength.setter
    def wavelength(self, value):
        if value < 0:
            raise ValueError
        self.properties['wavelength'] = value

    def _compute_rytov_variance(self, length):
        """Compute rytov variance for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        ## Returns
        `rytov_var` : float
            Rytov variance for given length.
        """
        k = 2*np.pi/self.wavelength
        rytov_var = 1.23*self.Cn2*k**(7/6)*length**(11/6)
        return rytov_var
    
    def _compute_wandering_variance(self, length):
        """Compute beam wandering variance for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        ## Returns
        `wandering_var` : float
            Beam wandering variance for given length [m^2].
        """
        k = 2*np.pi/self.wavelength
        Lambda_0 = 2*length/(k*self.W0**2)
        Theta_0 = 1
        rytov_var = self._compute_rytov_variance(length)
        f = lambda xi: (Theta_0 + (1 - Theta_0)*xi)**2 + 1.63*(rytov_var)**(6/5)*Lambda_0*(1 - xi)**(16/5)
        integrand = lambda xi: xi**2/f(xi)**(1/6)
        wandering_var = 7.25*self.Cn2*self.W0**(-1/3)*length**3*quad(integrand, 0, 1)[0]
        return wandering_var
    
    def _compute_scintillation_index_plane(self, rytov_var, length):
        """Compute aperture-averaged scintillation index of plane wave for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        `rytov_var` : float 
            Rytov variance.
        ## Returns
        `scint_index` : float
            Scintillation index for requested input parameters.
        """
        k = 2*np.pi/self.wavelength
        d = np.sqrt(k*self.rx_aperture**2/(4*length))
        first_term = 0.49*rytov_var/(1 + 0.65*d**2 + 1.11*rytov_var**(6/5))**(7/6)
        second_term = 0.51*rytov_var*(1 + 0.69*rytov_var**(6/5))**(-5/6)/(1 + 0.9*d**2 + 0.62*d**2*rytov_var**(6/5))
        scint_index = np.exp(first_term + second_term) - 1
        return scint_index
    
    def _compute_scintillation_index_spherical(self, rytov_var, length):
        """Compute aperture-averaged scintillation index of spherical wave for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        `rytov_var` : float 
            Rytov variance.
        ## Returns
        `scint_index` : float
            Scintillation index for requested input parameters.
        """
        k = 2*np.pi/self.wavelength
        d = np.sqrt(k*self.rx_aperture**2/(4*length))
        beta_0_sq = 0.4065*rytov_var
        first_term = 0.49*beta_0_sq/(1 + 0.18*d**2 + 0.56*beta_0_sq**(6/5))**(7/6)
        second_term = 0.51*beta_0_sq*(1 + 0.69*beta_0_sq**(6/5))**(-5/6)/(1 + 0.9*d**2 + 0.62*d**2*beta_0_sq**(6/5))
        return np.exp(first_term + second_term) - 1

    def _compute_coherence_width_plane(self, length):
        """Compute coherence width of plane wave for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        ## Returns
        `coherence_width` : float
            Coherence width for requested input parameters.
        """ 
        k = 2*np.pi/self.wavelength
        coherence_width = (0.42*length*self.Cn2*k**2)**(-3/5) 
        return coherence_width
    
    def _compute_coherence_width_spherical(self, length):
        """Compute coherence width of spherical wave for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        ## Returns
        `coherence_width` : float
            Coherence width for requested input parameters.
        """ 
        k = 2*np.pi/self.wavelength
        coherence_width = (0.16*length*self.Cn2*k**2)**(-3/5) 
        return coherence_width

    def _compute_coherence_width_gaussian(self, length):
        """Compute coherence width of gaussian wave for a horizontal channel (valid also for the 
        strong tubulence regime) [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        ## Returns
        `coherence_width` : float
            Coherence width for requested input parameters.
        """
        k = 2*np.pi/self.wavelength
        Lambda_0 = 2*length/(k*self.W0**2)
        Lambda = Lambda_0/(1 + Lambda_0**2)
        Theta = 1/(1 + Lambda_0**2)
        rho_plane = (1.46*self.Cn2*length*k**2)**(-3/5)
        q = length/(k*rho_plane**2)
        Theta_e = (Theta - 2*q*Lambda/3)/(1 + 4*q*Lambda/3)
        Lambda_e = Lambda/(1 + 4*q*Lambda/3)
        if Theta_e >= 0:
            a_e = (1 - Theta_e**(8/3))/(1 - Theta_e)
        else:
            a_e = (1 + np.abs(Theta_e)**(8/3))/(1 - Theta_e)
        coherence_width = 2.1*rho_plane*(8/(3*(a_e + 0.618*Lambda_e**(11/6))))**(3/5)
        return coherence_width
    
    def _compute_long_term_beam_size_at_receiver(self, rytov_var, length):
        """Compute long-term beamsize at receiver for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        `rytov_var` : float 
            Rytov variance.
        ## Returns
        `W_LT` : float
            Long-term beamsize at receiver for requested input parameters [m].
        """
        k = 2*np.pi/self.wavelength
        W_LT = self.W0*np.sqrt(1 + (self.wavelength*length/(np.pi*self.W0**2))**2 + 1.63*rytov_var**(6/5)*2*length/(k*self.W0**2))
        return W_LT
    
    def _compute_short_term_beam_size_at_receiver(self, long_term_beamsize, wandering_var):
        """Compute short-term beamsize at receiver for a horizontal channel [Andrews/Phillips, 2005].

        ## Parameters
        `long_term_beamsize` : float 
            Long-term beamsize at the receiver [m].
        `wandering_var` : float 
            Beam wandering variance at receiver [m^2].
        ## Returns
        `W_ST` : float
            Short-term beamsize at receiver for requested input parameters [m].
        """
        W_ST = np.sqrt(long_term_beamsize**2 - wandering_var)
        return W_ST
    
    def _compute_pdt(self, eta, length):
        """Compute probability distribution of atmospheric transmittance (PDT) [Vasylyev et al., 2018].

        ## Parameters
        `eta` : np.ndarray
            Input random variable values to calculate PDT for.
        `length` : float 
            Length of the channel [km].
        ## Returns
        `integral` : np.ndarray
            PDT function for input eta.
        """
        z = length*1e3
        rx_radius = self.rx_aperture/2
        rytov_var = self._compute_rytov_variance(z)
        pointing_var = (self.pointing_error*z)**2
        wandering_var = (self._compute_wandering_variance(z) + pointing_var)*(1 - self.tracking_efficiency)
        wandering_percent = 100*np.sqrt(wandering_var)/rx_radius
        if wandering_percent > 100:
            print("Warning ! The total wandering is larger than the aperture of the receiver. Use smaller values of pointing error.")

        W_LT = self._compute_long_term_beam_size_at_receiver(rytov_var, z)
        W_ST = self._compute_short_term_beam_size_at_receiver(W_LT, wandering_var)

        X = (rx_radius/W_ST)**2
        T0 = np.sqrt(1 - np.exp(-2*X))
        l = 8 * X * np.exp(-4*X) * i1(4*X) / (1 - np.exp(-4*X)*i0(4*X))/np.log(2*T0**2/(1 - np.exp(-4*X)*i0(4*X)))
        R = rx_radius * np.log(2*T0**2/(1 - np.exp(-4*X)*i0(4*X)))**(-1./l)

        if wandering_var >= 1e-7:
            pdt = lognegative_weibull_pdf(eta, T0, wandering_var, R, l)
        else: 
            pdt = np.zeros(np.size(eta))
            delta_eta = np.abs(eta[1] - eta[0])
            pdt[np.abs(eta - T0)  < delta_eta] = 1
        return pdt
    
    def _compute_channel_pdf(self, eta_ch, length):
        """Compute probability density function (PDF) of free-space channel efficiency.

        ## Parameters
        `eta_ch` : np.ndarray
            Input random variable values to calculate pdf for.
        `length` : float 
            Length of the channel [km].
        ## Returns
        `ch_pdf` : np.ndarray
            Channel PDF for input eta.
        """
        pdt = self._compute_pdt(eta_ch, length)
        pdt = pdt/np.sum(pdt)

        z = length*1e3
        n = lut_zernike_index_pd["n"]
        n = np.array(lut_zernike_index_pd["n"].values)
        rytov_var = self._compute_rytov_variance(z)
        r0 = self._compute_coherence_width_gaussian(z)
        eta_s = (1 + rytov_var)**(-1/4)
        bj2 = smf.bn2(self.rx_aperture, r0, n, self.obs_ratio)
        beta_opt = smf.beta_opt(self.obs_ratio)
        eta_smf_max = smf.eta_0(self.obs_ratio, beta_opt)
        ch_pdf = pdt*self.Tatm*smf.eta_ao(bj2)*eta_s*eta_smf_max
        return ch_pdf
    
    def _compute_mean_channel_efficiency(self, eta_ch, length, detector_efficiency = 1):
        """Compute mean channel efficiency, including losses at the detector.

        ## Parameters
        `eta_ch` : np.ndarray
            Input random variable values to calculate pdf for.
        `length` : float 
            Length of the channel [km].
        `detector_efficiency` : float
            Efficiency of detector at receiver (default 1).
        ## Returns
        `ch_pdf` : np.ndarray
            Channel PDF for input eta.
        """
        pdt = self._compute_pdt(eta_ch, length)
        pdt = pdt/np.sum(pdt)
        z = length*1e3
        n = lut_zernike_index_pd["n"]
        n = np.array(lut_zernike_index_pd["n"].values)
        rytov_var = self._compute_rytov_variance(z)
        r0 = self._compute_coherence_width_gaussian(z)
        eta_s = (1 + rytov_var)**(-1/4)
        bj2 = smf.bn2(self.rx_aperture, r0, n,self.obs_ratio)
        beta_opt = smf.beta_opt(self.obs_ratio)
        eta_smf_max = smf.eta_0(self.obs_ratio, beta_opt)
        mean_transmittance = np.sum(eta_ch*pdt)*self.Tatm*smf.eta_ao(bj2)*eta_s*eta_smf_max*detector_efficiency
        return mean_transmittance
    
    def _draw_pdt_sample(self, length, n_samples):
        """Draw random sample from probability distribution of atmospheric transmittance (PDT).

        ## Parameters
        `length` : float 
            Length of the channel [km].
        `n_samples` : int
            Number of samples to return.
        ## Returns
        `samples` : float
            Random samples of PDT.
        """
        eta = np.linspace(1e-7, 1, 1000)
        pdt = self._compute_pdt(eta, length)
        pdt = np.abs(pdt/np.sum(pdt))
        samples = np.random.choice(eta, n_samples, p = pdt)
        return samples

    def _draw_channel_pdf_sample(self, length, n_samples):
        """Draw random sample from free-space channel probability distribution.
    
        ## Parameters
        `length` : float 
            Length of the channel [km].
        `n_samples` : int
            Number of samples to return.
        ## Returns
        `samples` : float
            Random samples of channel PDF.
        """
        z = length*1e3
        eta = np.linspace(1e-7, 1, 1000)
        ch_pdf = self._compute_channel_pdf(eta, length)
        ch_pdf = np.abs(ch_pdf/np.sum(ch_pdf))
        ch_pdf_samples = np.random.choice(eta, n_samples, p = ch_pdf)
        rytov_var = self._compute_rytov_variance(z)
        scint_index = self._compute_scintillation_index_spherical(rytov_var, z)
        eta_s = (1 + scint_index)**(-1/4)
        beta_opt = smf.beta_opt(self.obs_ratio)
        eta_smf_max = smf.eta_0(self.obs_ratio, beta_opt)

        return self.Tatm*ch_pdf_samples*eta_smf_max*eta_s
    
    def _compute_loss_probability(self, length, n_samples):
        """Compute loss probability of photon in horizontal channel, taking all losses into account.

        ## Parameters
        `length` : float 
            Length of the channel [km].
        `n_samples` : int
            Number of samples to return.
        ## Returns
        `prob_loss` : float
            Probability that a photon is lost in the channel.
        """
        T = self._draw_channel_pdf_sample(length, n_samples)
        prob_loss = 1 - T
        return prob_loss
    
    def error_operation(self, qubits, **kwargs):
        """Error operation to apply to qubits.

        ## Parameters
        qubits : tuple of :obj:`~netsquid.qubits.qubit.Qubit`
            Qubits to apply noise to.
        """
        if 'channel' in kwargs:
            warn_deprecated("channel parameter is deprecated. "
                            "Pass length parameter directly instead.",
                            key="FreeSpaceLossModel.compute_model.channel")
            kwargs['length'] = kwargs['channel'].properties["length"]
            del kwargs['channel']

        prob_loss = self._compute_loss_probability(length = kwargs['length'], n_samples = len(qubits))
        for idx, qubit in enumerate(qubits):
            if qubit is None:
                continue
            self.lose_qubit(qubits, idx, prob_loss[idx], rng = self.properties['rng'])

class DownlinkChannel(QuantumErrorModel):
    """Model for photon loss on a downlink free-space channel.

    Uses probability density of atmospheric transmittance (PDT) from [Vasylyev et al., 2018] to
    sample the loss probability of the photon.

    ## Parameters
    ----------
    `W0` : float
        Waist radius of the beam at the transmitter [m].
    `rx_aperture` : float
        Diameter of the receiving telescope [m].
    `obs_ratio` : float
        Obscuration ratio of the receiving telescope.
    `n_max` : int
        Maximum radial index of correction of AO system.
    `Cn0` : float
        Reference index of refraction structure constant at ground level [m**(-2/3)].
    `wind_speed` : float
        Rms speed of the wind [m/s].  
    `wavelength` : float
        Wavelength of the radiation [m].
    `ground_station_alt` : float 
        Altitude of the ground station [km].
    `aerial_platform_alt` : float 
        Altitude of the aerial platform [km].
    `zenith_angle` : float
        Zenith angle of aerial platform [degrees].
    `pointing_error` : float
        Pointing error [rad].
    `tracking_efficiency` : float
        Efficiency of the coarse tracking mechanism.
    `Tatm` : float
        Atmospheric transmittance (square of the transmission coefficient).
    `integral_gain: float`
        Integral gain of the AO system integral controller.
    `control_delay: float`
        Delay of the AO system loop [s].
    `integration_time: float`
        Integration time of the AO system integral controller [s].
    `rng` : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    """
    def __init__(self, W0, rx_aperture, obs_ratio, n_max, Cn0, wind_speed, wavelength, 
                 ground_station_alt, aerial_platform_alt, zenith_angle, pointing_error = 0, 
                 tracking_efficiency = 0, Tatm = 1, integral_gain = 1, control_delay = 13.32e-4, integration_time = 6.66e-4, rng = None):
        super().__init__()
        self.rng = rng if rng else simtools.get_random_state()
        self.W0 = W0
        self.rx_aperture = rx_aperture
        self.obs_ratio = obs_ratio
        self.n_max = n_max
        self.Cn2 = Cn0
        self.wind_speed = wind_speed
        self.wavelength = wavelength
        self.ground_station_alt = ground_station_alt
        self.aerial_platform_alt = aerial_platform_alt
        self.zenith_angle = zenith_angle
        self.pointing_error = pointing_error
        self.integral_gain = integral_gain
        self.control_delay = control_delay
        self.integration_time = integration_time
        self.tracking_efficiency = tracking_efficiency
        self.Tatm = Tatm
        self.required_properties = ['length']

    @property
    def rng(self):
        """ :obj:`~numpy.random.RandomState`: Random number generator."""
        return self.properties['rng']

    @rng.setter
    def rng(self, value):
        if not isinstance(value, np.random.RandomState):
            raise TypeError("{} is not a valid numpy RandomState".format(value))
        self.properties['rng'] = value
    
    @property
    def Tatm(self):
        """ :float: atmosphere transmittance. """
        return self.properties['Tatm']

    @Tatm.setter
    def Tatm(self, value):
        if (value < 0) or (value > 1):
            raise ValueError
        self.properties['Tatm'] = value

    @property
    def pointing_error(self):
        """ :float: pointing error variance. """
        return self.properties['pointing_error']

    @pointing_error.setter
    def pointing_error(self, value):
        if (value < 0):
            raise ValueError
        self.properties['pointing_error'] = value

    @property
    def tracking_efficiency(self):
        """ :float: efficiency of the coarse tracking mechanism. """
        return self.properties['tracking_efficiency']

    @tracking_efficiency.setter
    def tracking_efficiency(self, value):
        if (value < 0) or (value > 1):
            raise ValueError
        self.properties['tracking_efficiency'] = value

    @property
    def W0(self):
        """float: beam waist radius at the transmitter [m]."""
        return self.properties['W0']

    @W0.setter
    def W0(self, value):
        if value < 0:
            raise ValueError
        self.properties['W0'] = value

    @property
    def rx_aperture(self):
        """float: diameter of the receiving telescope [m]."""
        return self.properties['rx_aperture']

    @rx_aperture.setter
    def rx_aperture(self, value):
        if value < 0:
            raise ValueError
        self.properties['rx_aperture'] = value

    @property
    def obs_ratio(self):
        """float: obscuration ratio of the receiving telescope."""
        return self.properties['obs_ratio']

    @obs_ratio.setter
    def obs_ratio(self, value):
        if value < 0 or (value > 1):
            raise ValueError
        self.properties['obs_ratio'] = value

    @property
    def n_max(self):
        """float: maximum radial index of correction of AO system."""
        return self.properties['n_max']

    @n_max.setter
    def n_max(self, value):
        if value < 0:
            raise ValueError
        self.properties['n_max'] = value


    @property
    def integral_gain(self):
        """float: integral gain of the AO system integral controller."""
        return self.properties['integral_gain']

    @integral_gain.setter
    def integral_gain(self, value):
        if value < 0:
            raise ValueError
        self.properties['integral_gain'] = value

    @property
    def control_delay(self):
        """float: delay of the AO system loop [s]."""
        return self.properties['control_delay']

    @control_delay.setter
    def control_delay(self, value):
        if value < 0:
            raise ValueError
        self.properties['control_delay'] = value

    @property
    def integration_time(self):
        """float: integration time of the AO system integral controller [s]."""
        return self.properties['integration_time']

    @integration_time.setter
    def integration_time(self, value):
        if value < 0:
            raise ValueError
        self.properties['integration_time'] = value
    
    @property
    def Cn0(self):
        """float: index of refraction structure constant [m**(-2/3)]."""
        return self.properties['Cn2']

    @Cn0.setter
    def Cn2(self, value):
        if value < 0:
            raise ValueError
        self.properties['Cn2'] = value

    @property
    def wavelength(self):
        """float: wavelength of the radiation [m]."""
        return self.properties['wavelength']

    @wavelength.setter
    def wavelength(self, value):
        if value < 0:
            raise ValueError
        self.properties['wavelength'] = value

    @property
    def ground_station_alt(self):
        """float: Altitude of the ground station [km]."""
        return self.properties['ground_station_alt']

    @ground_station_alt.setter
    def ground_station_alt(self, value):
        if value < 0:
            raise ValueError
        self.properties['ground_station_alt'] = value

    @property
    def aerial_platform_alt(self):
        """float: Altitude of the aerial platform [km]."""
        return self.properties['aerial_platform_alt']

    @aerial_platform_alt.setter
    def aerial_platform_alt(self, value):
        if value < 0:
            raise ValueError
        self.properties['aerial_platform_alt'] = value

    @property
    def zenith_angle(self):
        """float: Zenith angle of aerial platform [degrees]."""
        return self.properties['zenith_angle']

    @zenith_angle.setter
    def zenith_angle(self, value):
        if value < 0:
            raise ValueError
        self.properties['zenith_angle'] = value
    
    def _compute_Cn2(self, h):
        """Compute index of refraction structure constant [Andrews/Phillips, 2005].
        Uses the Hufnagel-Valley (HV) model.

        ## Parameters
        `h` : np.ndarray
            Values of h corresponding to slant path of the channel to integrate over [m].
        ## Returns
        `Cn2` : float
            Index of refraction structure constant [m^(-2/3)].
        """
        Cn2 = cn2.hufnagel_valley(h, self.wind_speed, self.Cn0)
        return Cn2
    
    def _compute_rytov_variance_plane(self):
        """Compute rytov variance of a plane wave for a downlink channel [Andrews/Phillips, 2005].

        ## Returns
        `rytov_var` : float
            Rytov variance for given length.
        """
        ground_station_alt = self.ground_station_alt*1e3
        aerial_platform_alt = self.aerial_platform_alt*1e3
        k = 2*np.pi/self.wavelength
        integrand = lambda h : self._compute_Cn2(h)*(h - ground_station_alt)**(5/6)
        rytov_var = 2.25*k**(7/6)*sec(self.zenith_angle)**(11/6)*quad(integrand, ground_station_alt, aerial_platform_alt)[0]
        return rytov_var

    def _compute_rytov_variance_spherical(self):
        """Compute rytov variance of a spherical wave for a downlink channel [Andrews/Phillips, 2005].

        ## Returns
        `rytov_var` :float
            Rytov variance for given length.
        """
        ground_station_alt = self.ground_station_alt*1e3
        aerial_platform_alt = self.aerial_platform_alt*1e3
        k = 2*np.pi/self.wavelength
        integrand = lambda h : self._compute_Cn2(h)*(h - ground_station_alt)**(5/6)*((aerial_platform_alt - h)/(aerial_platform_alt - ground_station_alt))**(5/6)
        rytov_var = 2.25*k**(7/6)*sec(self.zenith_angle)**(11/6)*quad(integrand, ground_station_alt, aerial_platform_alt)[0]
        return rytov_var
    
    def _compute_wandering_variance(self):
        """Compute beam wandering variance for a downlink channel [Andrews/Phillips, 2005].

        ## Returns
        `wandering_var` : float
            Beam wandering variance for given length [m^2].
        """
        ground_station_alt = self.ground_station_alt*1e3
        aerial_platform_alt = self.aerial_platform_alt*1e3
        k = 2*np.pi/self.wavelength
        length = 1e3*compute_channel_length(self.ground_station_alt, self.aerial_platform_alt, self.zenith_angle)
        Lambda_0 = 2*length/(k*self.W0**2)
        Theta_0 = 1
        rytov_var = self._compute_rytov_variance_spherical()
        f = lambda h: (Theta_0 + (1 - Theta_0)*(h - ground_station_alt)/(aerial_platform_alt - ground_station_alt))**2 + 1.63*(rytov_var)**(6/5)*Lambda_0*((aerial_platform_alt - h)/(aerial_platform_alt - ground_station_alt))**(16/5)
        integrand = lambda h: self._compute_Cn2(h)*(h - ground_station_alt)**2/f(h)**(1/6)
        wandering_var = 7.25*sec(self.zenith_angle)**3*self.W0**(-1/3)*quad(integrand, ground_station_alt, aerial_platform_alt)[0]
        return wandering_var
    
    def _compute_scintillation_index_plane(self, rytov_var, length):
        """Compute aperture-averaged scintillation index of plane wave for a downlink channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        `rytov_var` : float 
            Rytov variance.
        ## Returns
        `scint_index` : float
            Scintillation index for requested input parameters.
        """
        k = 2*np.pi/self.wavelength
        d = np.sqrt(k*self.rx_aperture**2/(4*length))
        first_term = 0.49*rytov_var/(1 + 0.65*d**2 + 1.11*rytov_var**(6/5))**(7/6)
        second_term = 0.51*rytov_var*(1 + 0.69*rytov_var**(6/5))**(-5/6)/(1 + 0.9*d**2 + 0.62*d**2*rytov_var**(6/5))
        return np.exp(first_term + second_term) - 1
    
    def _compute_scintillation_index_spherical(self, rytov_var, length):
        """Compute aperture-averaged scintillation index of spherical wave for a downlink channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        `rytov_var` : float 
            Rytov variance.
        ## Returns
        `scint_index` : float
            Scintillation index for requested input parameters.
        """
        k = 2*np.pi/self.wavelength
        d = np.sqrt(k*self.rx_aperture**2/(4*length))
        beta_0_sq = 0.4065*rytov_var
        first_term = 0.49*beta_0_sq/(1 + 0.18*d**2 + 0.56*beta_0_sq**(6/5))**(7/6)
        second_term = 0.51*beta_0_sq*(1 + 0.69*beta_0_sq**(6/5))**(-5/6)/(1 + 0.9*d**2 + 0.62*d**2*beta_0_sq**(6/5))
        return np.exp(first_term + second_term) - 1
    
    def _compute_coherence_width_plane(self):
        """Compute coherence width of plane wave for a downlink channel [Andrews/Phillips, 2005].

        ## Returns
        `coherence_width` : float
            Coherence width for requested input parameters.
        """ 
        ground_station_alt = self.ground_station_alt*1e3
        aerial_platform_alt = self.aerial_platform_alt*1e3
        k = 2*np.pi/self.wavelength
        integrand = lambda h : self._compute_Cn2(h)
        coherence_width = (0.42*k**2*sec(self.zenith_angle)*quad(integrand, ground_station_alt, aerial_platform_alt)[0])**(-3/5)
        return coherence_width
    
    def _compute_coherence_width_spherical(self):
        """Compute coherence width of spherical wave for a downlink channel [Andrews/Phillips, 2005].

        ## Returns
        `coherence_width` : float
            Coherence width for requested input parameters.
        """ 
        ground_station_alt = self.ground_station_alt*1e3
        aerial_platform_alt = self.aerial_platform_alt*1e3
        k = 2*np.pi/self.wavelength
        integrand = lambda h : self._compute_Cn2(h)*((aerial_platform_alt - h)/(aerial_platform_alt - ground_station_alt))**(5/3)
        coherence_width = (0.42*k**2*sec(self.zenith_angle)*quad(integrand, ground_station_alt, aerial_platform_alt)[0])**(-3/5)
        return coherence_width
    
    def _compute_coherence_width_gaussian(self, length):
        """Compute coherence width of gaussian wave for an downlink channel [Andrews/Phillips, 2005].

        ## Parameters 
        `length` : float
            Length of the channel [km].
        ## Returns
        `coherence_width` : float
            Coherence width for requested input parameters.
        """ 
        ground_station_alt = self.ground_station_alt*1e3
        aerial_platform_alt = self.aerial_platform_alt*1e3
        k = 2*np.pi/self.wavelength
        z = length*1e3
        Lambda_0 = 2*z/(k*self.W0**2)
        Lambda = Lambda_0/(1 + Lambda_0**2)
        Theta = 1/(1 + Lambda_0**2)
        Theta_bar = 1 - Theta
        integrand_1 = lambda h : self._compute_Cn2(h)*(Theta + Theta_bar*(aerial_platform_alt - h)/(aerial_platform_alt - ground_station_alt))**(5/3)
        mu_1d = quad(integrand_1, ground_station_alt, aerial_platform_alt)[0]
        integrand_2 = lambda h : self._compute_Cn2(h)*((h - ground_station_alt)/(aerial_platform_alt - ground_station_alt))**(5/3)
        mu_2d = quad(integrand_2, ground_station_alt, aerial_platform_alt)[0]
        coherence_width = (np.cos(np.deg2rad(self.zenith_angle))/(0.423*k**2*(mu_1d + 0.622*mu_2d*Lambda**(11/6))))**(3/5)
        return coherence_width
    
    def _compute_long_term_beam_size_at_receiver(self, rytov_var, length):
        """Compute long-term beamsize at receiver for a downlink channel [Andrews/Phillips, 2005].

        ## Parameters
        `length` : float 
            Length of the channel [m].
        `rytov_var` : float 
            Rytov variance.
        ## Returns
        `W_LT` : float
            Long-term beamsize at receiver for requested input parameters [m].
        """
        k = 2*np.pi/self.wavelength
        return self.W0*np.sqrt(1 + (self.wavelength*length/(np.pi*self.W0**2))**2 + 1.63*rytov_var**(6/5)*2*length/(k*self.W0**2))
    
    def _compute_short_term_beam_size_at_receiver(self, long_term_beamsize, wandering_var):
        """Compute short-term beamsize at receiver for a downlink channel [Andrews/Phillips, 2005].

        ## Parameters
        `long_term_beamsize` : float 
            Long-term beamsize at the receiver [m].
        `wandering_var` : float 
            Beam wandering variance at receiver [m^2].
        ## Returns
        `W_ST` : float
            Short-term beamsize at receiver for requested input parameters [m].
        """
        return np.sqrt(long_term_beamsize**2 - wandering_var)
    
    def _compute_lognormal_parameters(self, r, R, l, short_term_beamsize, scint_index):
        """Compute mean and standard deviation of lognormal distribution [Vasylyev et al., 2018].

        ## Parameters
        `r` : float 
            Deflection radius from center of receiver aperture [m].
        `R` : float 
            Weibull distribution R parameter.
        `l` : float 
            Weibull distribution l parameter.
        `short_term_beamsize` : float 
            Short-term beamsize at receiver [m].
        `scint_index` : float
            Scintillation index of horizontal channel.
        ## Returns
        `mu, sigma` : tuple (float, float)
            Mean value (mu) and standard deviation (sigma) of lognormal distribution.
        """
        rx_radius = self.rx_aperture/2
        eta_0 = 1 - np.exp(-2*rx_radius**2/short_term_beamsize**2)
        eta_mean = eta_0*np.exp(-(r/R)**l)
        eta_var = (1 + scint_index)*eta_mean**2
        mu = -np.log(eta_mean**2/np.sqrt(eta_var))
        sigma = np.sqrt(np.log(eta_var/eta_mean**2))
        return mu, sigma

    def _compute_pdt_parameters(self, length):
        """Compute parameters useful for the calculation of the probability distribution 
        of atmospheric transmittance [Vasylyev et al., 2018].

        ## Parameters
        `length` : float
            Length of the channel [km].
        ## Returns
        `lognormal_params` : function
            Output of _compute_lognormal_parameters. When evaluated at specific deflection radius r, returns
            mean value (mu) and standard deviation (sigma) of lognormal distribution at r.
        `wandering_var` : float 
            Beam wandering variance at receiver [m^2]. 
        `W_LT` : float
            Long-term beamsize at receiver for requested input parameters [m].
        """
        z = length*1e3
        rx_radius = self.rx_aperture/2
        pointing_var = (self.pointing_error*z)**2
        rytov_var = self._compute_rytov_variance_spherical()
        scint_index = self._compute_scintillation_index_spherical(rytov_var, z)
        W_LT = self._compute_long_term_beam_size_at_receiver(rytov_var, z)
        wandering_var = (self._compute_wandering_variance() + pointing_var)*(1 - self.tracking_efficiency)
        wandering_percent = 100*np.sqrt(wandering_var)/rx_radius

        if wandering_percent > 100:
            print('Warning ! The total wandering is larger than the aperture of the receiver. Use smaller values of pointing error.')

        W_ST = self._compute_short_term_beam_size_at_receiver(W_LT, wandering_var)

        X = (rx_radius/W_ST)**2
        T0 = np.sqrt(1 - np.exp(-2*X))
        l = 8 * X * np.exp(-4*X) * i1(4*X) / (1 - np.exp(-4*X)*i0(4*X))/np.log(2*T0**2/(1 - np.exp(-4*X)*i0(4*X)))
        R = rx_radius * np.log(2*T0**2/(1 - np.exp(-4*X)*i0(4*X)))**(-1./l)

        lognormal_params = lambda r : self._compute_lognormal_parameters(r, R, l, W_ST, scint_index)

        return lognormal_params, wandering_var, W_LT
    
    def _compute_pdt(self, eta, length):
        """Compute probability distribution of atmospheric transmittance (PDT) [Vasylyev et al., 2018].

        ## Parameters
        `eta` : np.ndarray
            Input random variable values to calculate PDT for.
        `length` : float
            Length of the channel [km].
        ## Returns
        `integral` : np.ndarray
            PDT function for input eta.
        """
        lognormal_params, wandering_var, W_LT = self._compute_pdt_parameters(length)
        if wandering_var == 0: 
            pdt = truncated_lognormal_pdf(eta, lognormal_params(0)[0], lognormal_params(0)[1])
        else:   
            integrand = lambda r: r*truncated_lognormal_pdf(eta, lognormal_params(r)[0], lognormal_params(r)[1])*np.exp(-r**2/(2*wandering_var))/wandering_var
            pdt = quad_vec(integrand, 0, self.rx_aperture/2 + W_LT)[0]
        return pdt

    def _compute_conversion_matrix(self, j_max):
        """Compute conversion matrix [Canuet et al., 2019].
        """
        Z = smf.compute_zernike(j_max)
        CZZ = smf.calculate_CZZ(self.rx_aperture/2, self.rx_aperture/2, Z, j_max, self.obs_ratio)
        M = smf.compute_conversion_matrix(CZZ)
        return M

    def _compute_attenuation_factors(self):
        """Compute attenuation factors of turbulent phase mode variances up to maximum order of correction n_max [Roddier, 1999].

        ## Returns
        `gamma_j` : np.ndarray
            Attenuation factors.
        """

        n = lut_zernike_index_pd["n"]
        n = np.array(lut_zernike_index_pd["n"].values)
        n_corrected = n[n <= self.n_max]
        open_loop_tf = lambda v: self.integral_gain*np.exp(-self.control_delay*v)*(1 - np.exp(-self.integration_time*v))/(self.integration_time*v)**2
        e_error = lambda v: 1/(1 + open_loop_tf(v))
        gamma_j = np.ones_like(n, dtype = float)
        cutoff_freq = 0.3*(n_corrected + 1)*self.wind_speed/self.rx_aperture
        for index in range(0, np.size(n_corrected)):
            if n_corrected[index] == 1:
                PSD_turbulence = lambda v: v**(-2/3) if v <= cutoff_freq[index] else v**(-17/3)
            else:
                PSD_turbulence = lambda v: v**(0) if v <= cutoff_freq[index] else v**(-17/3)
            gamma_j[index] = quad(lambda v: e_error(v)**2*PSD_turbulence(v), 1e-2, np.inf)[0]/quad(PSD_turbulence, 1e-2, np.inf)[0]

        return gamma_j
    
    def _compute_smf_coupling_pdf(self, eta_smf, eta_max, length):
        """Compute probability density function (PDF) of single mode fiber (SMF) coupling efficiency [Canuet et al., 2018].

        ## Parameters
        `eta_smf` : np.ndarray
            Input random variable values to calculate pdf for.
        `eta_max` : float
            Theoretical maximum coupling efficiency.
        `length` : float
            Length of the channel [km].
        ## Returns
        `smf_pdf` : np.ndarray
            SMF PDF for input eta.
        """
        z = length*1e3
        n = np.array(lut_zernike_index_pd["n"].values)
        j_Noll_as_index = np.array(lut_zernike_index_pd["j_Noll"].values) - 2
        rytov_var = self._compute_rytov_variance_spherical()

        #Check of the condition for aperture averaging
        if rytov_var <1:
                check = np.sqrt(self.wavelength*length*1e3)
        else:
                check = 0.36*np.sqrt(self.wavelength*length*1e3)* (rytov_var**(-3/5))
        # if self.rx_aperture < check:
        #     print("Warning ! The aperture averaging hypothesis is not valid for this set of parameters. Use bigger values of receiving aperture size")
        
        scint_index = self._compute_scintillation_index_spherical(rytov_var, z)
        r0 = self._compute_coherence_width_gaussian(z)
        eta_s = np.exp(-np.log(1 + scint_index))
        bj2 = smf.bn2(self.rx_aperture, r0, n, self.obs_ratio)
        gamma_j = self._compute_attenuation_factors()
        bj2 = bj2*gamma_j

        # Check if we are below the Rayleigh criterion
        bj_wvln = np.sqrt(bj2)/(2*np.pi)
        bj_wlvn_max = np.max(bj_wvln)
        # if bj_wlvn_max > 0.05:
        #     print(f" Warning ! The maximum Zernike coefficient std in wavelenghts is {bj_wlvn_max}. The SMF PDF is accurate below the Rayleigh criterion (0.05). You may need to use higher order of correction or smaller integration time of the AO system.")

        beta_opt = smf.beta_opt(self.obs_ratio)
        eta_smf_max = smf.eta_0(self.obs_ratio, beta_opt)
        eta_max = eta_max*eta_s*eta_smf_max
        smf_pdf = smf.compute_eta_smf_probability_distribution(eta_smf, eta_max, bj2)
        return smf_pdf
    
    def _compute_channel_pdf(self, eta_ch, length):
        """Compute probability density function (PDF) of free-space channel efficiency [Scriminich et al., 2022].

        ## Parameters
        `eta_ch` : np.ndarray
            Input random variable values to calculate pdf for.
        `length` : float
            Length of the channel [km].
        ## Returns
        `ch_pdf` : np.ndarray
            Channel PDF for input eta.
        """
        N = 10
        pdt = self._compute_pdt(eta_ch, length)
        pdt = pdt/np.sum(pdt)
        eta_rx = np.random.choice(eta_ch, N, p = pdt)
        integral = 0
        for index in range(0, N):
            integral = integral + self._compute_smf_coupling_pdf(eta_ch, eta_rx[index], length)
        ch_pdf = integral
        ch_pdf = self.Tatm*ch_pdf/np.sum(ch_pdf)
        return ch_pdf
    
    def _compute_mean_channel_efficiency(self, eta_ch, length, detector_efficiency = 1):
        """Compute mean channel efficiency, including losses at the detector.

        ## Parameters
        `eta_ch` : np.ndarray
            Input random variable values to calculate pdf for.
        `length` : float 
            Length of the channel [km].
        `detector_efficiency` : float
            Efficiency of detector at receiver (default 1).
        ## Returns
        `ch_pdf` : np.ndarray
            Channel PDF for input eta.
        """
        pdt = self._compute_pdt(eta_ch, length)
        pdt = pdt/np.sum(pdt)
        
        z = length*1e3
        n = np.array(lut_zernike_index_pd["n"].values)
        j_Noll_as_index = np.array(lut_zernike_index_pd["j_Noll"].values) - 2

        rytov_var = self._compute_rytov_variance_spherical()

        #Check of the condition for aperture averaging
        if rytov_var <1:
                check = np.sqrt(self.wavelength*length*1e3)
        else:
                check = 0.36*np.sqrt(self.wavelength*length*1e3)* (rytov_var**(-3/5))
        # if self.rx_aperture < check:
        #     print("Warning ! The aperture averaging hypothesis is not valid for this set of parameters. Use bigger values of receiving aperture size")

        scint_index = self._compute_scintillation_index_spherical(rytov_var, z)
        r0 = self._compute_coherence_width_gaussian(z)
        eta_s = np.exp(-np.log(1 + scint_index))
        bj2 = smf.bn2(self.rx_aperture, r0, n,self.obs_ratio)

        gamma_j = self._compute_attenuation_factors()
        bj2 = bj2*gamma_j

        # Check if we are below the Rayleigh criterion
        bj_wvln = np.sqrt(bj2)/(2*np.pi)
        bj_wlvn_max = np.max(bj_wvln)
        # if bj_wlvn_max > 0.05:
        #     print(f" Warning ! The maximum Zernike coefficient std in wavelenghts is {bj_wlvn_max}. The SMF PDF is accurate below the Rayleigh criterion (0.05). You may need to use higher order of correction or smaller integration time of the AO system.")

        beta_opt = smf.beta_opt(self.obs_ratio)
        eta_smf_max = smf.eta_0(self.obs_ratio, beta_opt)

        mean_transmittance = np.sum(eta_ch*pdt)*self.Tatm*eta_s*eta_smf_max*detector_efficiency *smf.eta_ao(bj2)
        return mean_transmittance
    
    def _draw_pdt_sample(self, length):
        """Draw random sample from probability distribution of atmospheric transmittance (PDT).

        ## Parameters
        `length` : float
            Length of the channel [km].
        ## Returns
        `sample` : float
            Random sample of PDT.
        """
        eta = np.linspace(1e-7, 1, 1000)
        pdt = self._compute_pdt(eta, length)
        pdt = np.abs(pdt/np.sum(pdt))
        sample = np.random.choice(eta, 1, p = pdt)
        return sample
    
    def _draw_smf_pdf_sample(self, length):
        """Draw random sample from probability distribution of single-mode fiber (SMF) coupling efficiency [Canuet et al., 2018].

        ## Parameters
        `length` : float
            Length of the channel [km].
        ## Returns
        `sample` : float
            Random sample of PDT.
        """
        eta = np.linspace(1e-7, 1, 1000)
        smf_pdf = self._compute_smf_coupling_pdf(eta, 1, length)
        smf_pdf = np.abs(smf_pdf/np.sum(smf_pdf))
        plt.figure()
        plt.plot(eta, smf_pdf)
        plt.show()
        sample = np.random.choice(eta, 1, p = smf_pdf)
        return sample
    
    def _draw_channel_pdf_sample(self, length, n_samples):
        """Draw random sample from free-space channel probability distribution [Scriminich et al., 2022].
        To be more efficient, the sample is calculated as the product of a sample from the PDT and the SMF coupling efficiency PDF,
        instead of the function of the channel PDF. 

        ## Parameters
        `length` : float
            Length of the channel [km].
        `n_samples` : int
            Number of samples to return.
        ## Returns
        `sample` : float
            Random sample of channel PDF.
        """
        eta = np.linspace(1e-4, 1, 1000)
        pdt = self._compute_pdt(eta, length)
        pdt = np.abs(pdt/np.sum(pdt))
        pdt /= pdt.sum()
        smf_pdf = self._compute_smf_coupling_pdf(eta, 1, length)
        smf_pdf = np.abs(smf_pdf/np.sum(smf_pdf))
        smf_pdf /= smf_pdf.sum()
        pdt_sample = np.random.choice(eta, n_samples, p = pdt)
        smf_sample = np.random.choice(eta, n_samples, p = smf_pdf)
        sample = self.Tatm*pdt_sample*smf_sample
        return sample
    
    def _compute_loss_probability(self, length, n_samples):
        """Compute loss probability of photon in downlink channel, taking all losses into account.

        ## Parameters
        `length` : float
            Length of the channel [km].
        `n_samples` : int
            Number of samples to return.
        ## Returns
        `prob_loss` : float
            Probability that a photon is lost in the channel.
        """
        T = self._draw_channel_pdf_sample(length, n_samples)
        prob_loss = 1 - T
        return prob_loss
    
    def error_operation(self, qubits, **kwargs):
        """Error operation to apply to qubits.

        Parameters
        ----------
        qubits : tuple of :obj:`~netsquid.qubits.qubit.Qubit`
            Qubits to apply noise to.

        """
        if 'channel' in kwargs:
            warn_deprecated("channel parameter is deprecated. "
                            "Pass length parameter directly instead.",
                            key="FreeSpaceLossModel.compute_model.channel")
            kwargs['length'] = kwargs['channel'].properties["length"]
            del kwargs['channel']

        prob_loss = self._compute_loss_probability(length = kwargs['length'], n_samples = len(qubits))
        for idx, qubit in enumerate(qubits):
            if qubit is None:
                continue
            self.lose_qubit(qubits, idx, prob_loss[idx], rng = self.properties['rng'])
        
class CachedChannel(QuantumErrorModel):
    """Class that performs error operation on qubits from precalculated probability of loss samples saved in an array
    to speed up execution time.

    ## Parameters
    ----------
    `loss_array` : np.ndarray
        Probability of loss samples array.
    `rng` : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    """
    
    def __init__(self, loss_array, rng = None):
        super().__init__()
        self.loss_array = loss_array
        self.rng = rng if rng else simtools.get_random_state()
        
    @property
    def rng(self):
        """ :obj:`~numpy.random.RandomState`: Random number generator."""
        return self.properties['rng']
    
    @rng.setter
    def rng(self, value):
        if not isinstance(value, np.random.RandomState):
            raise TypeError("{} is not a valid numpy RandomState".format(value))
        self.properties['rng'] = value      

    @property
    def loss_array(self):
        """ :np.ndarray: probability of loss samples array. """
        return self.properties['loss_array']

    def loss_array(self, value):
        if np.any((value < 0) | (value > 1)):  # Properly use '|' for element-wise logical OR
            raise ValueError("Values in loss_array must be in the range [0, 1]")
        self.properties['loss_array'] = value
    
    def error_operation(self, qubits, **kwargs):
        """Error operation to apply to qubits.

        Parameters
        ----------
        qubits : tuple of :obj:`~netsquid.qubits.qubit.Qubit`
        Qubits to apply noise to.

        """
        if 'channel' in kwargs:
            warn_deprecated("channel parameter is deprecated. "
                            "Pass length parameter directly instead.",
                            key="FreeSpaceLossModel.compute_model.channel")
            kwargs['length'] = kwargs['channel'].properties["length"]
            del kwargs['channel']

        for idx, qubit in enumerate(qubits):
            if qubit is None:
                continue
            prob_loss = random.choice(self.loss_array)
            # print(prob_loss)
            self.lose_qubit(qubits, idx, prob_loss, rng=self.properties['rng'])

class UplinkChannel(DownlinkChannel):
    """Model for photon loss on an uplink free-space channel.

    In the current implementation and considering the principle of reciprocity, the same model as
    the downlink is used (which is inherited), with the difference that the anisoplanatic error is 
    taken into account when calculating the residual zernike coefficient variances after correction.

    ## Parameters
    ----------
    `R_rx` : float
        Radius of the receiver aperture on the balloon [m].
    `D_tx` : float
        Diameter of the transmitting telescope [m].
    `obs_ratio` : float
        Obscuration ratio of the transmitting telescope.
    `n_max` : int
        Maximum radial index of correction of AO system.
    `Cn0` : float
        Reference index of refraction structure constant at ground level [m**(-2/3)].
    `wind_speed` : float
        Rms speed of the wind [m/s].  
    `wavelength` : float
        Wavelength of the radiation [m].
    `ground_station_alt` : float 
        Altitude of the ground station [km].
    `aerial_platform_alt` : float 
        Altitude of the aerial platform [km].
    `zenith_angle` : float
        Zenith angle of aerial platform [degrees].
    `pointing_error` : float
        Pointing error [rad].
    `tracking_efficiency` : float
        Efficiency of the coarse tracking mechanism.
    `Tatm` : float
        Atmospheric transmittance (square of the transmission coefficient).
    `integral_gain: float`
        Integral gain of the AO system integral controller.
    `control_delay: float`
        Delay of the AO system loop [s].
    `integration_time: float`
        Integration time of the AO system integral controller [s].
    `rng` : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    """

    def __init__(self, R_rx, D_tx, obs_ratio, n_max, Cn0, wind_speed, wavelength, 
                 ground_station_alt, aerial_platform_alt, zenith_angle, pointing_error = 0, 
                 tracking_efficiency = 0, Tatm = 1, integral_gain = 1, control_delay = 13.32e-4, integration_time = 6.66e-4, rng = None):
        super().__init__(R_rx, D_tx, obs_ratio, n_max, Cn0, wind_speed, wavelength, 
                         ground_station_alt, aerial_platform_alt, zenith_angle, pointing_error, 
                         tracking_efficiency, Tatm, integral_gain, control_delay, integration_time, rng)
        self.D_tx = D_tx
        self.R_rx = R_rx

    def _compute_anisoplanatic_error(self, length):
        """Compute anisoplanatic error

        ## Parameters
        `length` : float
            Length of the channel [m]
        ## Returns
        `var_aniso` : float
            Wavefront error variance attributed to anisoplanatism.
        """
        ground_station_alt = self.ground_station_alt*1e3
        aerial_platform_alt = self.aerial_platform_alt*1e3
        k = 2*np.pi/self.wavelength

        Lambda_0 = 2*length/(k*self.R_rx**2)
        Lambda = Lambda_0/(1 + Lambda_0**2)
        Theta = 1/(1 + Lambda_0**2)
        Theta_bar = 1 - Theta
        integrand_1 = lambda h : self._compute_Cn2(h)*(Theta + Theta_bar*(h - ground_station_alt)/(aerial_platform_alt - ground_station_alt))**(5/3)
        mu_1u = quad(integrand_1, ground_station_alt, aerial_platform_alt)[0]
        integrand_2 = lambda h : self._compute_Cn2(h)*((aerial_platform_alt - h)/(aerial_platform_alt - ground_station_alt))**(5/3)
        mu_2u = quad(integrand_2, ground_station_alt, aerial_platform_alt)[0]
        isoplanatic = (np.cos(np.deg2rad(self.zenith_angle))**(8/5))/(self.aerial_platform_alt - self.ground_station_alt)/(2.91*k**2*(mu_1u + 0.62*mu_2u*Lambda**(11/6)))**(3/5)

        var_isoplanatic = (self.pointing_error/isoplanatic)**(5/3)

        return var_isoplanatic

    def _compute_smf_coupling_pdf(self, eta_smf, eta_max, length):
        """Compute probability density function (PDF) of single mode fiber (SMF) coupling efficiency [Canuet et al., 2018].
        Overwrites the parent DownlinkChannel method to include isoplanatic error.

        ## Parameters
        `eta_smf` : np.ndarray
            Input random variable values to calculate pdf for.
        `eta_max` : float
            Theoretical maximum coupling efficiency.
        `length` : float
            Length of the channel [km].
        ## Returns
        `smf_pdf` : np.ndarray
            SMF PDF for input eta.
        """
        z = length*1e3
        n = np.array(lut_zernike_index_pd["n"].values)
        rytov_var = self._compute_rytov_variance_spherical()

        #Check of the condition for aperture averaging
        if rytov_var <1:
                check = np.sqrt(self.wavelength*length*1e3)
        else:
                check = 0.36*np.sqrt(self.wavelength*length*1e3)* (rytov_var**(-3/5))
        if self.rx_aperture < check:
            print("Problem aperture averaging: AO level : {} , Aperture : {} , check: {}".format(self.n_max,self.rx_aperture,check))
            # raise ValueError('The aperture averaging hypothesis is not valid for this set of parameters. Use bigger values of receiving aperture size')
        scint_index = self._compute_scintillation_index_spherical(rytov_var, z)
        r0 = self._compute_coherence_width_gaussian(z)
        eta_s = np.exp(-np.log(1 + scint_index))
        bj2 = smf.bn2(self.D_tx, r0, n,self.obs_ratio)

        gamma_j = self._compute_attenuation_factors()
        bj2 = bj2*gamma_j

        # Check if we are below the Rayleigh criterion
        bj_wvln = np.sqrt(bj2)/(2*np.pi)
        bj_wlvn_max = np.max(bj_wvln)
        if bj_wlvn_max > 0.05:
            print("Problem Rayleigh: AO level : {} , check: {}".format(self.n_max,bj_wlvn_max))
            # raise ValueError(f"The maximum Zernike coefficient std in wavelenghts is {bj_wlvn_max}. The SMF PDF is accurate below the Rayleigh criterion (0.05). You may need to use higher order of correction or smaller integration time of the AO system.")
        
        # Compute anisoplanatic error
        var_aniso = self._compute_anisoplanatic_error(z)

        beta_opt = smf.beta_opt(self.obs_ratio)
        eta_smf_max = smf.eta_0(self.obs_ratio, beta_opt)
        eta_max = eta_max*eta_s*eta_smf_max*np.exp(-var_aniso)
        smf_pdf = smf.compute_eta_smf_probability_distribution(eta_smf, eta_max, bj2)
        return smf_pdf

    def _compute_mean_channel_efficiency(self, eta_ch, length, detector_efficiency = 1):
        """Compute mean channel efficiency, including losses at the detector.
        Overwrites the parent DownlinkChannel method to include isoplanatic error.

        ## Parameters
        `eta_ch` : np.ndarray
            Input random variable values to calculate pdf for.
        `length` : float 
            Length of the channel [km].
        `detector_efficiency` : float
            Efficiency of detector at receiver (default 1).
        ## Returns
        `ch_pdf` : np.ndarray
            Channel PDF for input eta.
        """
        pdt = self._compute_pdt(eta_ch, length)
        pdt = pdt/np.sum(pdt)
        
        z = length*1e3
        n = np.array(lut_zernike_index_pd["n"].values)
        rytov_var = self._compute_rytov_variance_spherical()

        #Check of the condition for aperture averaging
        if rytov_var <1:
                check = np.sqrt(self.wavelength*length*1e3)
        else:
                check = 0.36*np.sqrt(self.wavelength*length*1e3)* (rytov_var**(-3/5))
        if self.rx_aperture < check:
            print("Problem aperture averaging: AO level : {} , Aperture : {} , check: {}".format(self.n_max,self.rx_aperture,check))
            # raise ValueError('The aperture averaging hypothesis is not valid for this set of parameters. Use bigger values of receiving aperture size')
        
        scint_index = self._compute_scintillation_index_spherical(rytov_var, z)
        r0 = self._compute_coherence_width_gaussian(z)
        eta_s = np.exp(-np.log(1 + scint_index))
        bj2 = smf.bn2(self.rx_aperture, r0, n,self.obs_ratio)

        gamma_j = self._compute_attenuation_factors()
        bj2 = bj2*gamma_j

        # Check if we are below the Rayleigh criterion
        bj_wvln = np.sqrt(bj2)/(2*np.pi)
        bj_wlvn_max = np.max(bj_wvln)
        if bj_wlvn_max > 0.05:
            print("Problem Rayleigh: AO level : {} , check: {}".format(self.n_max,bj_wlvn_max))
            # raise ValueError(f"The maximum Zernike coefficient std in wavelenghts is {bj_wlvn_max}. The SMF PDF is accurate below the Rayleigh criterion (0.05). You may need to use higher order of correction or smaller integration time of the AO system.")

        # Compute anisoplanatic error
        var_aniso = self._compute_anisoplanatic_error(z)

        beta_opt = smf.beta_opt(self.obs_ratio)
        eta_smf_max = smf.eta_0(self.obs_ratio, beta_opt)

        mean_transmittance = np.sum(eta_ch*pdt)*self.Tatm*smf.eta_ao(bj2)*eta_s*eta_smf_max*np.exp(-var_aniso)*detector_efficiency
        return mean_transmittance

if __name__ == "__main__":
    obs_ratio_ground = 0.4195
    obs_ratio_drone = 0
    n_max_ground = 1
    n_max_drone = 3
    n_samples = 10
    w0 = 0.15
    D_rx = 0.41
    Cn0 = 9.6*1e-14
    u_rms = 10
    wvln = 1550e-9
    ground_station_alt = 0.02
    aerial_platform_alt = 20
    zenith_angle = 0
    pointing_error = 0e-6
    tracking_efficiency = 0

    RE = 6371       
    H = 30
    RS = RE + H
    RG = RE
    theta_z = 70
    slant_length = np.sqrt(RS**2 + RG**2*(np.cos(np.deg2rad(theta_z))**2 - 1)) - RG*np.cos(np.deg2rad(theta_z))
    a = np.arccos((RG**2 + RS**2 - slant_length**2)/(2*RG*RS))
    length_horizontal = 2*np.sin(a)*(RE + H)
    transmittance_horiz = transmittance.horizontal(H, length_horizontal, wvln*1e9)
    Cn2_horizontal = cn2.hufnagel_valley(H*10**3, u_rms, Cn0)

    # Test horizontal channel
    obj = HorizontalChannel(w0, 0.2, obs_ratio_drone, Cn2_horizontal, wvln, pointing_error, tracking_efficiency, transmittance_horiz)  
    print(obj._draw_pdt_sample(length_horizontal, n_samples))
    print("rytov_var", obj._compute_rytov_variance(length_horizontal*1e3))
    print(obj._draw_channel_pdf_sample(length_horizontal, n_samples))

    # # Plot collection efficiency pdf
    # plt.figure()
    # D_rx = np.round(np.linspace(0.2, 1, 9), 2)
    # eta = np.linspace(1e-7, 1, 1000) 
    # for index in range(0, np.size(D_rx)):
    #     obj = DownlinkChannel(w0, D_rx[index], obs_ratio_ground, n_max_ground, Cn0, u_rms, wvln, ground_station_alt, aerial_platform_alt, zenith_angle, pointing_error, tracking_efficiency)
    #     length = compute_channel_length(ground_station_alt, aerial_platform_alt, zenith_angle)
    #     pdt = obj._compute_pdt(eta, length)
    #     pdt = np.abs(pdt/np.sum(pdt))
    #     plt.plot(eta, pdt, label = f"$D_{{\mathrm{{Rx}}}}$ = {D_rx[index]} m")
    # plt.legend()
    # plt.xlim([eta[0], eta[-1]])
    # plt.xlabel("Collection efficiency")
    # plt.ylabel("Counts")
    # plt.grid()
    # plt.savefig(f"collection_zenith{zenith_angle}.pdf")
    # plt.show()

    # Test downlink
    obj = DownlinkChannel(w0, D_rx, obs_ratio_ground, n_max_ground, Cn0, u_rms, wvln, ground_station_alt, aerial_platform_alt, zenith_angle, pointing_error)
    length = compute_channel_length(ground_station_alt, aerial_platform_alt, zenith_angle)
    print(obj._draw_pdt_sample(length))
    print(obj._draw_smf_pdf_sample(length))
    print(obj._draw_channel_pdf_sample(length, n_samples))
    
    # # print(obj._compute_loss_probability(length, n_samples))
    # eta_ch = np.linspace(1e-4, 1, 1000)
    # scint_pdf = obj._compute_scintillation_pdf(eta_ch)
    # print("mean scint", np.sum(scint_pdf*eta_ch))
    # plt.figure()
    # plt.plot(eta_ch, scint_pdf)
    # plt.show()
    # channel_pdf = obj._compute_channel_pdf(eta_ch, length)
    # plt.figure()
    # plt.semilogx(eta_ch, channel_pdf)
    # plt.show()

    #obj = UplinkChannel(0.25, 0.4, 9.6*1e-14, 10, 850e-9)
    #print(obj._draw_pdt_sample(0.010, 20, 80))
    
