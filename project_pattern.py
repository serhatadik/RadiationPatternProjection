from utils import *
from plotter import *
from tqdm import tqdm

def add_rad_patt(initial_preds, map_, map_res, UTM_boundaries, tx_antenna_raster_idx, tx_antenna_height, rx_antenna_height,
                 antenna_0_az_bearing_angle, antenna_0_el_deviation_angle_from_zenith, antenna_threeD_gain,
                 antenna_inclined_tow_bearing_angle=90):
    (dim1, dim2) = initial_preds.shape
    elevs = np.zeros((dim1, dim2))
    bearings = np.zeros((dim1, dim2))
    tx_elevation = map_[tx_antenna_raster_idx[0], tx_antenna_raster_idx[1]]
    tx_total_height = tx_elevation + tx_antenna_height
    assert antenna_0_el_deviation_angle_from_zenith >= 90 and antenna_0_el_deviation_angle_from_zenith <= 180

    beta = antenna_0_el_deviation_angle_from_zenith - 90
    gamma = antenna_inclined_tow_bearing_angle - 90
    ant_dir = np.matmul(rotz(gamma), np.matmul(roty(beta), np.array([0, 0, 1])))

    print("Now calculating angles for gain projection. \n")
    for i in tqdm(range(dim1)):
        for j in range(dim2):
            rx_total_height = map_[i, j] + rx_antenna_height

            slope = (rx_total_height - tx_total_height) / (
                        map_res * (np.linalg.norm(np.array([i, j]) - np.array(tx_antenna_raster_idx))))

            # Elevation angle

            vec3_tx_to_rx = np.array(
                [j - tx_antenna_raster_idx[1], i - tx_antenna_raster_idx[0], rx_total_height - tx_total_height])

            el_angle = find_angle_bw_vecs(vec3_tx_to_rx, ant_dir)

            elevs[i, j] = el_angle

            vec_tx_to_rx = np.array([i, j]) - np.array(tx_antenna_raster_idx)

            norm_vec_tx_to_rx = np.linalg.norm(vec_tx_to_rx)

            dot_product = np.dot(vec_tx_to_rx, np.array([1, 0]))

            angle_radians = np.arccos(dot_product / norm_vec_tx_to_rx)
            angle_degrees = np.degrees(angle_radians)

            if vec_tx_to_rx[1] < 0:
                angle_degrees = 360 - angle_degrees

            bearings[i, j] = angle_degrees

    bearings = np.nan_to_num(bearings, nan=0)
    elevs = np.nan_to_num(elevs, nan=0)
    #plotter(elevs, "Elevation before Subtraction")
    final_preds = np.zeros_like(initial_preds)
    print("Now applying the projection. \n")
    for i in tqdm(range(dim1)):
        for j in range(dim2):
            if bearings[i, j] - antenna_0_az_bearing_angle >= 0:
                bearings[i, j] = bearings[i, j] - antenna_0_az_bearing_angle
            else:
                bearings[i, j] = 360 + bearings[i, j] - antenna_0_az_bearing_angle

            if elevs[i, j] - 90 >= 0:
                elevs[i, j] = elevs[i, j] - 90
            else:
                elevs[i, j] = 360 + elevs[i, j] - 90

            bearings[i, j] = round(bearings[i, j])
            if bearings[i, j] == 360:
                bearings[i, j] = 359

            elevs[i, j] = round(elevs[i, j])
            if elevs[i, j] == 360:
                elevs[i, j] = 359
            final_preds[i, j] = initial_preds[i, j] + antenna_threeD_gain[(bearings[i, j], elevs[i, j])]
    N1 = map_.shape[0]
    N2 = map_.shape[1]
    en = UTM_boundaries

    UTM_long = np.linspace(en[0, 2], en[0, 3] - map_res, N1)
    UTM_lat = np.linspace(en[0, 0], en[0, 1] - map_res, N2)
    plotter(final_preds, "Projected Antenna Gain", UTM_long, UTM_lat, "dBi")
    return final_preds




