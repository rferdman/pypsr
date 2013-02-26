import numpy as np

# Actually output is cos(phi0) where phi0 is the half width
# Based on equation (17) in Rafikov and Lai (2006), ApJ, 641, 438.

#  Precession perios is in DAYS, to match the units of time here, expressed as MJD (days).

def width_rl(t, alpha, delta, T1, rho, inclination, prec_period):

    phi_SO = 2.*np.pi*(t-T1)/prec_period  

    cos_zeta = -np.cos(delta)*np.cos(inclination) + np.sin(delta)*np.sin(inclination)*np.cos(phi_SO)
    sin_zeta = np.sqrt(1.0 - cos_zeta*cos_zeta)

    cos_phi0 = (np.cos(rho) - cos_zeta*np.cos(alpha))/(sin_zeta*np.sin(alpha))

    return cos_phi0


# These are based on equation (3) in Kramer (1998), ApJ, 509, 856.
def width_mk(t, alpha, lambd, T0, rho, inclination, prec_period):

    phi = 2.*np.pi*(T0 - t)/prec_period  


    cos_beta = np.cos(lambd)*np.cos(inclination) + \
        np.sin(lambd)*np.sin(inclination)*np.cos(phi)

    beta = np.arccos(cos_beta)
    sigma = beta - alpha

    sqrt_arg =  ( (np.sin(rho/2.))**2. - (np.sin(sigma/2.))**2. ) / \
                (np.sin(alpha) * np.sin(beta))   

    W = 4.0*np.arcsin( np.sqrt( sqrt_arg ))

    return W
