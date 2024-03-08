import os
import warnings
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

def read_antenna_pattern_file(pattern_folder_dir, pattern_file_name, extension):
    if not os.path.isdir(pattern_folder_dir):
        warnings.warn(f'Error: The following folder does not exist:\n{pattern_folder_dir}')
        return

    if pattern_file_name not in os.listdir(pattern_folder_dir):
        warnings.warn(f'Error: The file does not exist in the folder:\n {pattern_folder_dir}')
        return

    pattern_file = os.path.join(pattern_folder_dir, pattern_file_name)
    threeD_pat = {}

    with open(pattern_file, 'r') as file:
        lines = file.readlines()
        if extension in ['ana', '.ana', 'msi', '.msi']:
            azimuth_gain, elevation_gain = {}, {}
            is_horizontal, is_vertical = False, False
            for line in lines:
                line = line.strip()
                if line.startswith(('Horizontal', 'HORIZONTAL')):
                    is_horizontal, is_vertical = True, False
                elif line.startswith(('Vertical', 'VERTICAL')):
                    is_horizontal, is_vertical = False, True
                elif line:
                    if is_horizontal or is_vertical:
                        angle, gain = map(float, line.split())
                        if extension in ['msi', '.msi']:
                            gain = -gain  # specific adjustment for msi files
                        if is_horizontal:
                            azimuth_gain[angle] = gain
                        else:
                            elevation_gain[angle] = gain
            if len(azimuth_gain) == 1:  # Spread single azimuth gain across all angles if only one is present
                azimuth_gain = {i: list(azimuth_gain.values())[0] for i in np.arange(1.0, 360.0)}
            for az in azimuth_gain:
                for el in elevation_gain:
                    threeD_pat[(az, el)] = azimuth_gain[az] + elevation_gain[el]

        elif extension in ['pat', '.pat']:
            azimuth_gain, elevation_gain = {}, defaultdict(list)
            is_horizontal, is_vertical, flag, cnt_999 = False, False, 0, None
            for cnt, line in enumerate(lines, 1):
                line = line.strip()
                if cnt == 1:
                    KYPAT = line[-1]
                if KYPAT == '2' and cnt > 1:
                    if line == "999":
                        flag, cnt_999 = 1, cnt
                        continue
                    if flag == 1 and cnt == cnt_999 + 1:
                        num_splits, num_elv = line.split(',')
                    if flag == 0:
                        angle, gain = map(float, line.split(','))
                        azimuth_gain[angle] = gain
                    elif flag == 1 and cnt > cnt_999 + 2 and cnt not in [cnt_999 + int(num_elv) + 3 + n * 181 for n in
                                                                         range(int(num_splits))]:
                        angle, gain = map(float, line.split(','))
                        elevation_gain[angle].append(gain)
            if len(azimuth_gain) == 1:  # Spread single azimuth gain across all angles if only one is present
                azimuth_gain = {i: list(azimuth_gain.values())[0] for i in range(1, 360)}
            for az in azimuth_gain:
                for el in elevation_gain:
                    threeD_pat[(az, el)] = azimuth_gain[az] + np.mean(elevation_gain[el])

    fig = plt.figure(layout='constrained')
    ax1 = fig.add_subplot(1, 1, 1, projection='polar')
    ax1.set_theta_direction(-1)
    theta = np.array(list(elevation_gain.keys())) * np.pi / 180
    ax1.plot(theta, np.array(list(elevation_gain.values())))
    plt.title("Elevation Angles vs Gain (Polar Plot)")
    plt.show()

    fig = plt.figure(layout='constrained')
    ax1 = fig.add_subplot(1, 1, 1, projection='polar', theta_offset=np.pi / 2)
    ax1.set_theta_direction(-1)
    ax1.plot(np.array(list(azimuth_gain.keys())) * np.pi / 180, np.array(list(azimuth_gain.values())))
    plt.title("Azimuth Angles vs Gain (Polar Plot)")
    plt.show()

    return threeD_pat
