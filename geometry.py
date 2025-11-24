# geometry.py

import math

# Constants
EPSILON = 1e-10  # Small value to avoid numerical issues

def distance(atom1, atom2):
    """
    Calculate Euclidean distance between two atoms.
    
    Args:
        atom1: Atom object with x, y, z coordinates
        atom2: Atom object with x, y, z coordinates
    
    Returns:
        float: Distance in Angstroms
    """
    dx = atom2.x - atom1.x
    dy = atom2.y - atom1.y
    dz = atom2.z - atom1.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def calculate_angle(atom1, atom2, atom3):
    """
    Calculate angle between three atoms (atom2 is the vertex/center).
    
    Args:
        atom1: First atom
        atom2: Center atom (vertex of angle)
        atom3: Third atom
    
    Returns:
        float: Angle in radians (0 to π)
    
    Raises:
        ValueError: If atoms are collinear or overlapping
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
    
    # Check for degenerate case (atoms on top of each other)
    if ba_length < EPSILON or bc_length < EPSILON:
        raise ValueError("Cannot calculate angle: atoms are too close together")
    
    dot_product = ba_x*bc_x + ba_y*bc_y + ba_z*bc_z
    
    cos_theta = dot_product / (ba_length * bc_length)
    
    # Clamp to valid range to handle numerical errors
    # (but could optionally warn if significantly out of range)
    if cos_theta > 1.0:
        if cos_theta > 1.01:  # More than rounding error
            print(f"Warning: cos(theta) = {cos_theta} > 1.0, clamping to 1.0")
        cos_theta = 1.0
    if cos_theta < -1.0:
        if cos_theta < -1.01:
            print(f"Warning: cos(theta) = {cos_theta} < -1.0, clamping to -1.0")
        cos_theta = -1.0
    
    # Return angle in radians
    return math.acos(cos_theta)

def normalize_vector(x, y, z):
    """
    Normalize a vector to unit length.
    
    Args:
        x, y, z: Vector components
    
    Returns:
        tuple: (x_norm, y_norm, z_norm) - unit vector
    
    Note: Returns (0, 0, 0) for zero-length vectors to avoid division by zero.
    Callers should check for this case if it matters.
    """
    length = math.sqrt(x*x + y*y + z*z)
    
    if length < EPSILON:  # Avoid division by zero
        # Could optionally raise an error instead
        return 0.0, 0.0, 0.0
    
    return x/length, y/length, z/length

def dot_product(x1, y1, z1, x2, y2, z2):
    """
    Calculate dot product of two vectors.
    
    Args:
        x1, y1, z1: First vector components
        x2, y2, z2: Second vector components
    
    Returns:
        float: Dot product
    """
    return x1*x2 + y1*y2 + z1*z2

def vector_length(x, y, z):
    """
    Calculate length (magnitude) of a vector.
    
    Args:
        x, y, z: Vector components
    
    Returns:
        float: Length of vector
    """
    return math.sqrt(x*x + y*y + z*z)

def cross_product(x1, y1, z1, x2, y2, z2):
    """
    Calculate cross product of two vectors.
    
    Args:
        x1, y1, z1: First vector components
        x2, y2, z2: Second vector components
    
    Returns:
        tuple: (x, y, z) components of cross product vector
    """
    return (
        y1*z2 - z1*y2,
        z1*x2 - x1*z2,
        x1*y2 - y1*x2
    )

def angle_between_vectors(x1, y1, z1, x2, y2, z2):
    """
    Calculate angle between two vectors.
    
    Args:
        x1, y1, z1: First vector
        x2, y2, z2: Second vector
    
    Returns:
        float: Angle in radians
    """
    dot = dot_product(x1, y1, z1, x2, y2, z2)
    len1 = vector_length(x1, y1, z1)
    len2 = vector_length(x2, y2, z2)
    
    if len1 < EPSILON or len2 < EPSILON:
        raise ValueError("Cannot calculate angle between zero-length vectors")
    
    cos_theta = dot / (len1 * len2)
    cos_theta = max(-1.0, min(1.0, cos_theta))  # Clamp to [-1, 1]
    
    return math.acos(cos_theta)