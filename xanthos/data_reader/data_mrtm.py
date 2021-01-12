import numpy as np

from xanthos.data_reader.data_utils import DataUtils
from xanthos.utils.math import sub2ind


class DataMrtm(DataUtils):

    def __init__(self, config_obj, coords):

        self.config = config_obj

        map_index = sub2ind([self.config.ngridrow, self.config.ngridcol], coords[:, 4].astype(int) - 1,
                            coords[:, 3].astype(int) - 1)

        self.flow_dist = self.load_routing_data(self.config.flow_distance, self.config.ngridrow, self.config.ngridcol,
                                                map_index, rep_val=1000)

        self.flow_dir = self.load_routing_data(self.config.flow_direction, self.config.ngridrow, self.config.ngridcol,
                                               map_index)

        self.str_velocity = self.load_routing_data(self.config.strm_veloc, self.config.ngridrow, self.config.ngridcol,
                                                   map_index, rep_val=0)

        self.instream_flow = np.zeros((self.config.ncell, ), dtype=float)

        self.chs_prev = self.load_chs_data()

        super().__init__(nmonths=config_obj.nmonths)

    def load_chs_data(self):
        """Load channel velocity file into array if in future mode, else stage zeros array."""
        try:
            # Initialize channel storage/soil moisture.
            if self.config.HistFlag == "True":
                return np.zeros((self.config.ncell,), dtype=float)

            # For future runs, initialize with the last value of the historical channel storage/soil moisture
            else:
                return self.load_data(self.config.ChStorageFile, 0, self.config.ChStorageVarName)[:, -1]

        except AttributeError:
            return np.zeros((self.config.ncell,), dtype=float)

    def load_routing_data(self, fn, ngridrow, ngridcol, map_index, skip=68, rep_val=None):
        """Load routing data.

        DRT data, 280 x 720, -9999 for missing values, convert to 67420 X 1

        @:param fle             file to load
        @:param ngridrow        number of grids per row
        @:param ngridcol        number of grids per column
        @:param map_index
        @:param skip
        @:param rep_val         value to replace with when less than value
        """
        fd = self.load_data(fn)
        v = self.vectorize(fd, ngridrow, ngridcol, map_index, skip=skip)

        if rep_val is None:
            return v

        else:
            v[np.where(v < rep_val)[0]] = rep_val
            return v

    @staticmethod
    def vectorize(data, ngridrow, ngridcol, map_index, skip):
        """Convert 2D Map (360 x 720) Matrix to 1D Map(67420)."""
        new = np.zeros((ngridrow, ngridcol), dtype=float) - 9999

        for i in range(0, data.shape[0]):
            new[i + skip, :] = data[data.shape[0] - 1 - i, :]

        new = new.reshape((ngridrow * ngridcol,), order='F')

        return new[map_index]
