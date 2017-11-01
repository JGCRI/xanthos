from scipy.io import loadmat


class ReadParams:
    """
    Read parameters from calibrated input file.
    Calibrated using VIC run via MATLAB GA function under WATCH

    1971 Jan - 1990 Dec = 240 months.

    @:param basin_id            basin ids
    @:param pars                A, B, C, D, M parameters
    @:param kge                 error metric
    @:param precipitation       precipitation
    @:param pet                 Potential Evapotranspiration
    @:param maxstor             maximum soil water storage
    @:param rsim                simulated runoff (Q)
    @:param robs                observed runoff from VIC model
    @:param ea                  actual evapotranspiration (ET)
    @:param gw                  groundwater storage
    @:param sm                  soil moisture storage
    @:param recharge            groundwater recharge
    @:param dr                  direct runoff
    @:param baseflow            baseflow
    @:param snowpack            snowpack
    @:param rainfall            rainfall
    @:param snowfall            snowfall
    @:param snowmelt            snowmelt
    @:param obsbfi              observed base flow index
    @:param simbfi              simulated base flow index
    """

    def __init__(self, f, v):

        # read in .MAT file
        m = self.get_data(f, v)

        # assign attributes
        self.basin_id = m[0]
        self.pars = m[1]
        self.kge = m[2]
        self.precipitation = m[3]
        self.pet = m[4]
        self.maxstor = m[5]
        self.rsim = m[6]
        self.robs = m[7]
        self.ea = m[8]
        self.gw = m[9]
        self.sm = m[10]
        self.recharge = m[11]
        self.dr = m[12]
        self.baseflow = m[13]
        self.snowpack = m[14]
        self.rainfall = m[15]
        self.snowfall = m[16]
        self.snowmelt = m[17]
        self.obsbfi = m[18]
        self.simbfi = m[19]

    @staticmethod
    def get_data(f, v):
        """
        Read input .MAT calibration file.  Returns 20 variables.
        """
        return loadmat(f)[v][0][0]


class ReadWatch:

    def __init__(self, f, v):

        # read in .MAT file
        m = self.get_data(f, v)

        """
        Watch plus from ISI-MIP.
        
        1971 jan - 1990 dec = 240 months
        
        Use parameters from 1990 to spin-up (usually 5 to 10 years)
        Spin-up time to let model reach steady-state 
        
        Yaling used 1991 - 2010 as validation
        """
        # assign attributes
        self.basin_id = m[0] # 0.5 basin IDs
        self.grids = m[1] # coords
        self.pet = m[2] # PET hargreaves (inputs:
        self.p = m[3] # watch+wfdi data
        self.tmin = m[4] # watch degree Celcius (from watch Kelvin) min air temp per month
        self.robs = m[5] # VIC (mm/month)
        self.maxstor = m[6] # not used!
        self.thickness = m[7] # soil thickness for each grid cell; not used!
        self.psi = m[8] # not used
        self.lmb = m[9] # not used
        self.phi = m[10] # not used
        self.bfi = m[11] # not used in model only GA obj function
        self.vicsurf = m[12] # derived from vic runoff and BFI (mm/month)
        self.vicbase = m[13] # derived from vic

    @staticmethod
    def get_data(f, v):
        """
        Read input .MAT calibration file.  Returns 20 variables.
        """
        return loadmat(f)[v][0][0]


if __name__ == "__main__":

    f = '/users/ladmin/repos/github/xanthos/example/input/runoff/abcd/climate/watch/WATCH_basin_1_grid.mat'
    fc = '/users/ladmin/repos/github/xanthos/example/input/runoff/abcd/climate/calibrated_dist/Cal_dist_noconstrain_basin_1.mat'


    b = ReadWatch(f, 'WATCH_basin_grid')
    c = ReadParams(fc, 'dist_nocstrn_cal')

    # 40 years - 480 months
    print b.robs.shape

    # 20 years - 240 months
    print c.robs.shape