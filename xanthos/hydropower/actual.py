"""
Calculate country and GCAM-region level hydropower production from Xanthos streamflow.

Created on April 25, 2017

@author: Sean Turner (sean.turner@pnnl.gov)
@Project: Xanthos V2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import math
import numpy as np
import pandas as pd
import os


class HydropowerActual:
    """
    Compute country and GCAM-region level hydropower production time series based on xanthos streamflow output.

    Gridded streamflow is used to drive dam simulations for 1593 large dams (~54% global hydropower installed capacity).
    Each dam is has been pre-assigned an optimized look-up table, which assigns a turbine release at each time step
    using current reservoir storage, inflow and month of year. Power output time series are summed across dams for each
    country and scaled up to account for unsimulated capacity (i.e., for dams not represented in the model).
    """

    secs_in_month = 2629800  # number of seconds in an average month
    cumecs_to_Mm3permonth = 2.6298  # m3/s to Mm3/month
    sww = 9810  # specific weight of water (N/m^3)
    hours_in_year = 8766  # number of hours in a year
    mwh_to_exajoule = 3.6 * (10 ** -9)  # megawatts to exajoules

    def __init__(self, settings, q_grids):
        """Load inputs, run simulation, and output results."""
        self.settings = settings
        self.q_grids = q_grids

        # read input data
        self.res_data = pd.read_csv(settings.HydroDamData)  # Read dam data
        self.grid_data = pd.read_csv(settings.GridData)  # Read grid data
        self.drainage_area = np.loadtxt(settings.DrainArea)  # Read upstream drainage area for all grids
        self.missing_cap = pd.read_csv(settings.MissingCap)  # Read missing installed capacity (i.e., dams not represented)
        self.rule_curves = np.load(settings.rule_curves)  # rule curves for all 1593 dams

        # assign from inputs
        fname = "actual_hydro_by_gcam_region_EJperyr_{}.csv".format(settings.ProjectName)
        self.filename_hydro = os.path.join(settings.OutputFolder, fname)
        self.start_date = settings.hact_start_date  # Get start date for simulation "M/YYYY"
        self.loc_refs = self.grid_data[["ID", "long", "lati"]]  # Get latitiude and longitude for all grid squares
        self.grid_ids = self.res_data.iloc[:, 0:2].apply(self.get_grid_id, 1)  # Get grid indices for all dams
        self.q = np.transpose(self.q_grids[self.grid_ids - 1, :])  # get flows for all grids containing dams
        self.dr_ar_assumed = np.apply_along_axis(self.get_drain_area,
                                                 1,
                                                 np.array(self.loc_refs)[self.grid_ids - 1, 1:3],
                                                 self.drainage_area)  # Get assumed drainage areas (grid)
        self.dr_ar_actual = np.array(self.res_data["CATCH"])  # Get actual drainage areas
        self.q_ = self.q * self.dr_ar_actual / self.dr_ar_assumed  # Correct inflow by drainage
        self.q_Mm3 = self.q_ * HydropowerActual.cumecs_to_Mm3permonth  # Get inflow in Mm^3 / month
        self.power_all_dams = np.empty([len(self.q_Mm3), len(self.res_data)])  # array to hold power data

        # assigned during run
        self.rc = None
        self.inflow = None
        self.cap = None
        self.cap_live = None
        self.head = None
        self.q_max = None
        self.efficiency = None
        self.l_inflow = None
        self.s = None
        self.power = None
        self.s_states = None
        self.q_states = None
        self.env_flow = None
        self.hydro_gcam_regions_EJ = None

        # run simulation
        self.hydro_sim()

        # convert power production time series to GCAM region enery production
        self.to_region()

        # write output
        self.write_output()

    @staticmethod
    def find_nearest(array, value):
        """Get value from within an array closest to a value."""
        idx = (np.abs(array - value)).idxmin()  # idxmin instead of argmin
        return array[idx]

    @staticmethod
    def find_nearest_idx(array, value):
        """Get index of an array closest to a value."""
        return (np.abs(array - value)).idxmin()

    def get_grid_id(self, longlat):
        """Get the grid location of a dam based on longitude/latitude."""
        lon = self.find_nearest(self.loc_refs["long"], longlat["LONG_DD"])
        lat = self.find_nearest(self.loc_refs["lati"], longlat["LAT_DD"])
        return int(self.loc_refs[(self.loc_refs["long"] == lon) & (self.loc_refs["lati"] == lat)]["ID"])

    def get_drain_area(self, x, drainage_area):
        """Get drainage area implied by the routing network."""
        lonseq = np.unique(self.loc_refs["long"])
        latseq = np.unique(self.loc_refs["lati"])[::-1]
        return drainage_area[latseq.tolist().index(x[1]), lonseq.tolist().index(x[0])]

    def env_flow_constraint(self):
        """Apply environmental flow constraints."""
        inflow_mmf = self.inflow.groupby(self.inflow.index.month).mean()
        inflow_maf = self.inflow.mean()
        inflow_efr_perc = ((inflow_mmf < 0.4 * inflow_maf) * 1 * 0.6 +
                           ((inflow_mmf >= 0.4 * inflow_maf) & (inflow_mmf <= 0.8 * inflow_maf)) * 1 * 0.45 +
                           (inflow_mmf > 0.8 * inflow_maf) * 1 * 0.3 * [not i for i in (inflow_mmf < 1)])
        self.env_flow = inflow_efr_perc * inflow_mmf

    def init_sim(self):
        """Initialize simulation."""
        self.s = [self.cap] * (self.l_inflow + 1)  # Declare storage variable to be simulated
        self.power = [0] * self.l_inflow
        self.q_states = self.inflow.groupby(self.inflow.index.month).quantile((0, 0.2375, 0.4750, 0.7125, 0.95, 1))

    def get_power(self, res):
        """
        MAIN SIMULATION FUNCTION.

        NEXT VERSION TO INCLUDE ev, area, max_depth, installed_cap
        """
        # start simulation (method for getting q state should be improved--currently inaccurate)
        for t in range(0, self.l_inflow):
            mth = self.inflow.index[t].month
            s_state = self.s[t] / self.cap
            env = self.env_flow[mth]
            active = self.s[t] + self.inflow[t] - (self.cap - self.cap_live)
            r_breaks = np.linspace(0, 1, 5)
            release_options = self.rc[:, mth - 1]
            release_options[(np.isnan(release_options))] = 1.1
            release = r_breaks[(release_options <= s_state)][-1] * self.q_max
            r = min(min(max(release, env), active), self.q_max)
            self.s[t + 1] = max(min(self.s[t] + self.inflow[t] - r, self.cap), 0)
            h = (np.mean([self.s[t], self.s[t + 1]]) / self.cap) * self.head
            self.power[t] = max(self.efficiency * HydropowerActual.sww * h * (r / HydropowerActual.secs_in_month), 0)

        self.power_all_dams[:, res] = self.power  # MegaWatts

    def sim_vars(self, idx, res):
        """Calculate simulation variable for a target reservoir."""
        # get rule curves
        self.rc = self.rule_curves[:, :, res]

        self.inflow = \
            pd.DataFrame(self.q_Mm3).set_index(pd.period_range(self.start_date, periods=len(self.q_Mm3), freq="M"))[res]
        self.l_inflow = len(self.inflow)
        self.cap = self.res_data["CAP"][res]
        self.cap_live = self.res_data["CAPLIVE"][res]

        if math.isnan(self.cap_live) is True:
            self.cap_live = self.cap

        self.installed_cap = self.res_data["ECAP"][res]
        self.q_max = self.res_data["FLOW_M3S"][res] * HydropowerActual.cumecs_to_Mm3permonth
        self.efficiency = self.res_data["EFF"][res]
        self.head = self.res_data["HEAD"][res]

        if math.isnan(self.head) is True:
            self.head = self.installed_cap / (
                self.efficiency * HydropowerActual.sww * (self.q_max / HydropowerActual.secs_in_month))

    def hydro_sim(self):
        """Simulate hydropower operations and add power to array."""
        for idx, res in enumerate(self.res_data.index.tolist()):
            # calculate sim variables
            self.sim_vars(idx, res)

            # environmental flow constraint
            self.env_flow_constraint()

            # initialize simulation
            self.init_sim()

            # get power
            self.get_power(res)

    def to_region(self):
        """Convert power production time series to GCAM region enery production."""
        power_all_dams_monthly = pd.DataFrame(self.power_all_dams).set_index(
            pd.period_range(self.start_date, periods=len(self.power_all_dams), freq="M")
        )
        energy_all_dams = power_all_dams_monthly.resample("A").mean() * (
            HydropowerActual.hours_in_year * HydropowerActual.mwh_to_exajoule
        )
        energy_all_countries = energy_all_dams.groupby(self.res_data["COUNTRY"], axis=1).sum()
        energy_all_countries_total = energy_all_countries.multiply(list(self.missing_cap["factor"]))
        self.hydro_gcam_regions_EJ = energy_all_countries_total.groupby(list(self.missing_cap["GCAM_ID"]), axis=1).sum()

    def write_output(self):
        """Write results to CSV."""
#        df = pd.DataFrame(self.hydro_gcam_regions_EJ)
#        cols = ['region_{}'.format(i) for i in df.columns]
#        cols.insert(0, 'year')
#        df.columns = cols
        odf = self.hydro_gcam_regions_EJ.T
        odf.reset_index(inplace=True)
        odf.rename(columns={'index': 'region'}, inplace=True)
        pd.DataFrame.to_csv(odf, self.filename_hydro, index=False)
