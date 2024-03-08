import numpy as np


def rotz(gamma):
    rotmat = np.array([[np.cos(np.deg2rad(gamma)), -np.sin(np.deg2rad(gamma)), 0],
                       [np.sin(np.deg2rad(gamma)), np.cos(np.deg2rad(gamma)), 0], [0, 0, 1]])
    return rotmat

def rotx(alpha):
    rotmat = np.array([[1, 0, 0], [0, np.cos(np.deg2rad(alpha)), -np.sin(np.deg2rad(alpha))],
                       [0, np.sin(np.deg2rad(alpha)), np.cos(np.deg2rad(alpha))]])
    return rotmat

def roty(beta):
    rotmat = np.array([[np.cos(np.deg2rad(beta)), 0, np.sin(np.deg2rad(beta))], [0, 1, 0],
                       [-np.sin(np.deg2rad(beta)), 0, np.cos(np.deg2rad(beta))]])
    return rotmat

def find_angle_bw_vecs(vec1, vec2):
    dot_ = np.dot(vec1, vec2)
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    angle_deg = np.rad2deg(np.arccos(dot_ / (norm1 * norm2)))
    return angle_deg

