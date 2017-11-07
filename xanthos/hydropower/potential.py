'''
Created on April 25, 2017

Original:
@author: Sean Turner (sean.turner@pnnl.gov)

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
'''

import numpy as np
import pandas as pd
import os

def HydropowerPotential(settings, q_grids):

    sww = 9810 # specific weight of water N/m^3
    av_hours_in_month = 730.5 # hours in a month (average)
    watthr_to_twh = 10 ** -12 # convert watt-hours to Tw-hr
    twh_to_exajoule = 0.0036 # convert Tw-hr to Exajoules

    hyd_grid_data = pd.read_csv(settings.GridData)  # Contains elevations for energy head
    q_grids_ = np.apply_along_axis(constrain_q, 1, q_grids, settings.q_ex)

    # COMPUTE ENERGY PRODUCTION FOR GCAM REGIONS (ExaJoules)
    q_ex_h = settings.ef * sww * q_grids_ * av_hours_in_month * watthr_to_twh  # Energy before head is multiplied
    e_grids = q_ex_h * hyd_grid_data["elevD"][:, np.newaxis]  # Multiply by elevation difference (head) to get energy (TWh)
    e_grids_ = pd.DataFrame(e_grids.T)
    e_grids__ = e_grids_.set_index(pd.period_range(settings.hpot_start_date, periods=len(e_grids_), freq="M"))  # Time-stamp the data
    e_grids_annual = e_grids__.resample("A").sum() * twh_to_exajoule  # Sum to annual energy production, then convert to EJ

    # TECHNICAL POTENTIAL
    techpot = e_grids_annual.groupby(hyd_grid_data["regID"], axis=1).sum()
    techpot_ = techpot.drop(techpot.columns[0], axis=1)
    filename_techpot = os.path.join(settings.OutputFolder, "Tech_hyd_pot_by_GCAM_region_EJ_{}.csv".format(settings.ProjectName))
    pd.DataFrame.to_csv(techpot_, filename_techpot)

    # TECHNICAL, EXPLOITABLE POTENTIAL
    techpot_expl = e_grids_annual.groupby(hyd_grid_data["regID"] * hyd_grid_data["inGrandELEC"], axis=1).sum()
    techpot_expl_ = techpot_expl.drop(techpot_expl.columns[0], axis=1)

    filename_techpot_expl = os.path.join(settings.OutputFolder, "Tech_expl_hyd_pot_by_GCAM_region_EJ_{}.csv".format(settings.ProjectName))
    pd.DataFrame.to_csv(techpot_expl_, filename_techpot_expl)


# CONSTRAIN EACH GRID OF INFLOW TO REPRESENT MAXIMUM TURBINE CAPACITY, DEFINED BY q_ex
def constrain_q(q, ex):
    q_max = np.percentile(q, ex * 100)
    return np.clip(q, 0, q_max)