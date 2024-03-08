from read_pattern_file import *
import scipy.io as sio
from coordinate_transformations import *
from project_pattern import *

def main():
    threeD_pat = read_antenna_pattern_file("./", "aw3939_3600_T0.msi", "msi")

    # Load map structure from .mat file
    map_struct = sio.loadmat("SLCmap_5May2022.mat")['SLC']

    # Extract SLC data and create a column map for easy access
    SLC = map_struct[0][0]
    column_map = dict(zip([name for name in SLC.dtype.names], [i for i in range(len(SLC.dtype.names))]))

    # Calculate the map data
    map_ = SLC[column_map['dem']] + 0.3048 * SLC[column_map['hybrid_bldg']]

    # Define additional variables for the add_rad_patt function
    map_res = SLC[column_map["cellsize"]]
    tx_antenna_height = 3
    rx_antenna_height = 1.5
    tirem_preds = np.zeros(map_.shape)
    antenna_0_az_bearing_angle = 120
    antenna_0_el_deviation_angle_from_zenith = 90
    antenna_inclined_tow_bearing_angle=90
    antenna_threeD_gain = threeD_pat
    endpoint_coords = [40.76895, -111.84167] # example basestation coordinate
    [BS_x, BS_y, indx] = lon_lat_to_grid_xy(endpoint_coords[1], endpoint_coords[0], SLC,
                                        column_map)  # establish basestation pixel location.

    tx_antenna_raster_idx = [BS_y, BS_x]

    add_rad_patt(tirem_preds, map_, map_res, SLC[column_map["axis"]], tx_antenna_raster_idx, tx_antenna_height, rx_antenna_height,
                 antenna_0_az_bearing_angle, antenna_0_el_deviation_angle_from_zenith, antenna_threeD_gain,
                 antenna_inclined_tow_bearing_angle=90)

if __name__ == "__main__":
    main()