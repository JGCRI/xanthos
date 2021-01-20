import numpy as np

from xanthos.data_reader.data_utils import DataUtils
from xanthos.data_reader.data_reference import DataReference
from xanthos.utils.math import sub2ind


class DataMrtmManaged(DataUtils):

    NCELL = 67420
    NGRIDROW = 360
    NGRIDCOL = 720

    def __init__(self, config_obj=None, start_year=1971, end_year=2001, flow_distance_file=None,
                 flow_direction_file=None, stream_velocity_file=None, historical_mode="True",
                 hist_channel_storage_file=None, hist_channel_storage_varname=None):

        if config_obj is None:
            self.start_year = start_year
            self.end_year = end_year
            self.nmonths = (self.end_year - self.start_year + 1) * 12
            self.historical_mode = historical_mode

            # load reference data
            self.reference_data = DataReference(nmonths=self.nmonths)

            # routing files
            self.flow_distance_file = flow_distance_file
            self.flow_direction_file = flow_direction_file
            self.stream_velocity_file = stream_velocity_file

            # only used when historical mode is False to initialize channel storage
            self.hist_channel_storage_file = hist_channel_storage_file
            self.hist_channel_storage_varname = hist_channel_storage_varname

        else:
            self.start_year = config_obj.StartYear
            self.end_year = config_obj.EndYear
            self.nmonths = config_obj.nmonths
            self.reference_data = DataReference(config=config_obj)
            self.flow_distance_file = config_obj.flow_distance
            self.flow_direction_file = config_obj.flow_direction
            self.stream_velocity_file = config_obj.strm_veloc
            self.historical_mode = config_obj.HistFlag

            if self.historical_mode != "True":
                self.hist_channel_storage_file = config_obj.ChStorageFile
                self.hist_channel_storage_varname = config_obj.ChStorageVarName

        map_index = sub2ind([DataMrtmManaged.NGRIDROW, DataMrtmManaged.NGRIDCOL], self.reference_data.coords[:, 4].astype(int) - 1,
                            self.reference_data.coords[:, 3].astype(int) - 1)

        self.flow_dist = self.load_routing_data(self.flow_distance_file,
                                                DataMrtmManaged.NGRIDROW,
                                                DataMrtmManaged.NGRIDCOL,
                                                map_index, rep_val=1000)

        self.flow_dir = self.load_routing_data(self.flow_direction_file,
                                               DataMrtmManaged.NGRIDROW,
                                               DataMrtmManaged.NGRIDCOL,
                                               map_index)

        self.str_velocity = self.load_routing_data(self.stream_velocity_file,
                                                   DataMrtmManaged.NGRIDROW,
                                                   DataMrtmManaged.NGRIDCOL,
                                                   map_index, rep_val=0)

        self.instream_flow = np.zeros((DataMrtmManaged.NCELL, ), dtype=float)

        self.chs_prev = self.load_chs_data()

        super().__init__(nmonths=self.nmonths)

    def load_chs_data(self):
        """Load channel velocity file into array if in future mode, else stage zeros array."""
        try:
            # Initialize channel storage/soil moisture.
            if self.historical_mode == "True":
                return np.zeros((DataMrtmManaged.NCELL,), dtype=float)

            # For future runs, initialize with the last value of the historical channel storage/soil moisture
            else:
                return self.load_data(self.hist_channel_storage_file, 0, self.hist_channel_storage_varname)[:, -1]

        except AttributeError:
            return np.zeros((DataMrtmManaged.NCELL,), dtype=float)

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
