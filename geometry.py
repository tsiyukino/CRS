# geometry.py

import math

def distance(atom1, atom2):
    # Calculate distance between two atoms
    dx = atom2.x - atom1.x
    dy = atom2.y - atom1.y
    dz = atom2.z - atom1.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def calculate_angle(atom1, atom2, atom3):
    """
    Calculate angle between three atoms (atom2 center)
    Returns angle in radians
    """
    # Vector from atom2 to atom1
    ba_x = atom1.x - atom2.x
    ba_y = atom1.y - atom2.y
    ba_z = atom1.z - atom2.z
    
    # Vector from atom2 to atom3
    bc_x = atom3.x - atom2.x
    bc_y = atom3.y - atom2.y
    bc_z = atom3.z - atom2.z
    
    ba_length = math.sqrt(ba_x*ba_x + ba_y*ba_y + ba_z*ba_z)
    bc_length = math.sqrt(bc_x*bc_x + bc_y*bc_y + bc_z*bc_z)
    
    dot_product = ba_x*bc_x + ba_y*bc_y + ba_z*bc_z
    
    cos_theta = dot_product / (ba_length * bc_length)
    
    # numerical errors
    if cos_theta > 1.0:
        cos_theta = 1.0
    if cos_theta < -1.0:
        cos_theta = -1.0
    
    # Return angle radians
    return math.acos(cos_theta)

def normalize_vector(x, y, z):

    # Return unit vector (length = 1)
    length = math.sqrt(x*x + y*y + z*z)

    if length < 1e-10:  # avoid division by zero
        return 0.0, 0.0, 0.0
    return x/length, y/length, z/length

def dot_product(x1, y1, z1, x2, y2, z2):
    return x1*x2 + y1*y2 + z1*z2

def vector_length(x, y, z):
    return math.sqrt(x*x + y*y + z*z)