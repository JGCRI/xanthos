"""
Functions for computing drought statistics from hydrological model
output.

Created on January 21, 2019

@author: Robert Link (robert.link@pnnl.gov), Caleb Braun (caleb.braun@pnnl.gov)
@Project: Xanthos V2.3

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2019, Battelle Memorial Institute
"""

import numpy as np
import pandas as pd
import logging
import os
from xanthos.data_writer.out_writer import OutWriter


class DroughtStats:
    """
    Analyze drought impacts based on runoff or soil moisture.

    Drought metrics based on work done by Sheffield and Wood (2008).

    Sheffield, J. and Eric F. Wood (2008) "Projected changes in drought
      occurrence under future global warming from multi-model,
      multi-scenario, IPCC AR4 simulations", _Clim Dyn_ 31: 79--105,
      doi:10.1007/s00382-007-0340-z
    """
    MONTHS_IN_YEAR = 12

    def __init__(self, settings, runoff, soil_moisture):
        """Run drought statistics based on given configuration settings."""
        # The calculations expect the output variable to be (ntime x ngrid), so
        # we need to transpose
        if settings.drought_var.lower() == 'q':
            hydroout = runoff.T
        elif settings.drought_var.lower() == 'soilmoisture':
            hydroout = soil_moisture.T
        else:
            raise ValueError("Invalid drought variable specified (must be 'q' or 'soilmoisture')")

        # Create a template for drought output files
        output_path = os.path.join(settings.OutputFolder, "drought_{}_{}".format("{}", settings.OutputNameStr))

        # Initialize an empty OutWriter so that outputs will be in correct format
        out_writer = OutWriter(settings, 0, {})

        # Calculate thresholds if they're not provided, otherwise compute statistics
        if settings.drought_thresholds is None:
            logging.info("\tCalculating drought thresholds")
            thresholds = self.calculate_thresholds(hydroout, settings)
            np.save(output_path.format("thresholds"), thresholds)
        else:
            logging.info("\tCalculating drought statistics")
            threshvals = np.load(settings.drought_thresholds)
            severity, intensity, duration = self.droughtstats(hydroout, threshvals)

            # Convert to the standard Xanthos output format and write out
            for varname, arr in zip(["severity", "intensity", "duration"], [severity, intensity, duration]):
                df = pd.DataFrame(arr.T, columns=[str(x) for x in range(arr.shape[0])])
                out_writer.write_data(output_path.format(varname), varname, df)

    @classmethod
    def calculate_thresholds(cls, histout, settings):
        """Calculate historical drought thresholds.

        :param histout:     matrix[ntime x ngrid] of hydrological outputs.
        :param settings:    main configuration settings
        :return:            array of quantiles
        """
        # Get subset of data over which we're calculating the thresholds
        syear = settings.threshold_start_year
        eyear = settings.threshold_end_year
        smonth = (syear - settings.StartYear) * cls.MONTHS_IN_YEAR
        emonth = (eyear + 1 - syear) * cls.MONTHS_IN_YEAR

        histout = histout[smonth:emonth, :]

        return cls.getthresh(histout, settings.threshold_nper)

    def droughtstats(self, hydroout, threshvals):
        """Compute Severity, Intensity, and Duration statistics for a matrix of hydrological output.

        :param hydroout: matrix[ntime x ngrid] of outputs for a hydrological
                         variable, such as runoff or soil moisture.
        :param threshvals: matrix[K x ngrid] of drought thresholds.  The K different
                         sets of thresholds will be recycled as necessary over all
                         the ntime time values.  Typically, either K==1 (single set
                         of thresholds for all times), or K==12 (different set of
                         thresholds for each month of the year).
        :return: (S, I, D) matrices of severity, intensity, and duration, each with
                         the same dimensions as hydroout.

        Duration is defined as the number of consecutive months that the
        hydrological variable has been under the drought threshold.  By definition,
        the duration of a grid cell that is not currently under drought conditions
        is zero.

        Severity is the cumulative sum since the start of the current drought
        episode of the amount by which the hydrological variable fell short of the
        drought threshold.  If q_0(t) is the (possibly time variable) threshold, and
        t_0 is the start of the current drought episode, then
        $$
           S(t) = \sum_{t'=t_0}^t q_0(t') - q(t').
        $$
        Severity is also defined to be zero for a grid cell that is not currently
        under drought conditions.

        Intensity is the average severity over a drought period.  I = S/D.  It is
        defined to be (you guessed it) zero for a grid cell not currently under
        drought conditions.
        """
        S = np.empty_like(hydroout)
        I = np.empty_like(hydroout)
        D = np.empty_like(hydroout)

        # Handle the first case specially, since there won't be a previous time step.

        isdrought = hydroout[0, :] < threshvals[0, :]
        D[0, :] = np.where(isdrought, 1.0, 0.0)
        I[0, :] = S[0, :] = np.where(isdrought, (threshvals[0, :] - hydroout[0, :]) / threshvals[0, :], 0.0)

        ntime = hydroout.shape[0]
        nthresh = threshvals.shape[0]

        for t in range(1, ntime):
            tm1 = t-1
            m = t % nthresh         # "month"
            thresh = threshvals[m, :]
            hydro = hydroout[t, :]

            dtm1 = D[tm1, :]
            stm1 = S[tm1, :]

            dt = D[t, :]
            st = S[t, :]
            it = I[t, :]

            isdrought = hydro < thresh
            dt[:] = np.where(isdrought, dtm1+1, 0.0)
            st[:] = np.where(isdrought, stm1 + (thresh - hydro) / thresh, 0.0)
            it[:] = np.where(isdrought, st / dt, 0.0)

        return (S, I, D)

    @staticmethod
    def getthresh(histout, nper, quantile=0.1):
        """
        Compute quantile-based thresholds using data from a reference period.

        :param histout: Matrix[ntime x ngrid] of output variable during the reference period.
        :param nper: Number of periods to calculate thresholds for.  Generally either
                     1 (single threshold for all periods), or 12 (thresholds by month)
        :param quantile: Quantile to use for determining the threshold.
        :return: Matrix[nper x ngrid] of threshold values.

        The array passed as histout should contain only values from the reference
        period, from which drought thresholds will be calculated.
        """
        (ntime, ngrid) = histout.shape
        # Check that ntime is a multiple of nper
        nyear = int(ntime / nper)

        histout = np.reshape(histout, (nyear, nper, ngrid))

        # Equivalent to np.quantile(histout, quantile, axis=0) in numpy versions >= 1.15.0
        return np.percentile(histout, quantile * 100, axis=0)
