import numpy as np

def compute_x_channeledTPST_positive(beta_tpst, beta_x_at_teca, delta_mu_xtpst, E, dEX, TECA_jaw, dx_at_teca, deltaP_P, dispersion_tpst):
    return np.sqrt(beta_tpst / beta_x_at_teca) * (
        np.cos(delta_mu_xtpst) + 
        np.sqrt((beta_x_at_teca * (E + dEX)) / ((TECA_jaw - dx_at_teca * deltaP_P) * (TECA_jaw - dx_at_teca * deltaP_P)) - 1) * 
        np.sin(delta_mu_xtpst)
    ) * (TECA_jaw - dx_at_teca * deltaP_P) + dispersion_tpst * deltaP_P

