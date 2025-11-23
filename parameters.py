# parameters.py

import math

# Bond parameters: {bond_type: (r0, kb)}
BOND_PARAMS = {
    'C-C': (1.54, 310.0),
    'C=C': (1.34, 600.0),
    'C-H': (1.09, 340.0),
    'C-O': (1.43, 320.0),
    'C=O': (1.23, 570.0),
    'O-H': (0.96, 450.0),
    'N-H': (1.01, 410.0),
    'N-C': (1.47, 320.0),
    'C-N': (1.47, 320.0),
    'H-H': (0.74, 340.0),
}

# Angle parameters: {angle_type: (theta0_degrees, ka)}
ANGLE_PARAMS = {
    'C-C-C': (109.5, 60.0),
    'C-C-H': (109.5, 50.0),
    'H-C-H': (109.5, 35.0),
    'H-O-H': (104.5, 55.0),
    'C-O-H': (108.5, 55.0),
    'O-C-O': (120.0, 80.0),
    'H-N-H': (106.7, 40.0),
    'C-N-H': (109.5, 50.0),
    'N-C-C': (109.5, 60.0),
}

# VDW parameters: {element: (sigma, epsilon)}
VDW_PARAMS = {
    'H': (2.5, 0.015),
    'C': (3.4, 0.086),
    'N': (3.25, 0.086),
    'O': (3.12, 0.120),
}

def get_bond_params(atom1_element, atom2_element):
    """
    Get bond parameters for two atoms
    Returns: (r0, kb) or None if not found
    """
    # Try both orders
    key1 = f"{atom1_element}-{atom2_element}"
    key2 = f"{atom2_element}-{atom1_element}"
    
    if key1 in BOND_PARAMS:
        return BOND_PARAMS[key1]
    elif key2 in BOND_PARAMS:
        return BOND_PARAMS[key2]
    else:
        return None

def get_angle_params(atom1_element, atom2_element, atom3_element):
    """
    Get angle parameters for three atoms (atom2 is center)
    Returns: (theta0_radians, ka) or None if not found
    """
    # Try both orders (reverse matters for asymmetric angles)
    key1 = f"{atom1_element}-{atom2_element}-{atom3_element}"
    key2 = f"{atom3_element}-{atom2_element}-{atom1_element}"
    
    if key1 in ANGLE_PARAMS:
        theta_deg, ka = ANGLE_PARAMS[key1]
        return (math.radians(theta_deg), ka)
    elif key2 in ANGLE_PARAMS:
        theta_deg, ka = ANGLE_PARAMS[key2]
        return (math.radians(theta_deg), ka)
    else:
        return None

def get_vdw_params(element):
    """
    Get VDW parameters for an element
    Returns: (sigma, epsilon) or default values if not found
    """
    if element in VDW_PARAMS:
        return VDW_PARAMS[element]
    else:
        # Default for unknown elements
        return (3.0, 0.05)

def combine_vdw_params(element1, element2):
    """
    Get combined VDW parameters for two elements
    Uses Lorentz-Berthelot combining rules
    Returns: (sigma_combined, epsilon_combined)
    """
    sigma1, epsilon1 = get_vdw_params(element1)
    sigma2, epsilon2 = get_vdw_params(element2)
    
    # Combining rules
    sigma_combined = (sigma1 + sigma2) / 2.0
    epsilon_combined = math.sqrt(epsilon1 * epsilon2)
    
    return (sigma_combined, epsilon_combined)

def add_bond_params(atom1_element, atom2_element, r0, kb):
    """
    Add custom bond parameters to database
    """
    key = f"{atom1_element}-{atom2_element}"
    BOND_PARAMS[key] = (r0, kb)

def add_angle_params(atom1_element, atom2_element, atom3_element, theta0_deg, ka):
    """
    Add custom angle parameters to database
    """
    key = f"{atom1_element}-{atom2_element}-{atom3_element}"
    ANGLE_PARAMS[key] = (theta0_deg, ka)

def add_vdw_params(element, sigma, epsilon):
    """
    Add custom VDW parameters to database
    """
    VDW_PARAMS[element] = (sigma, epsilon)

def print_available_params():
    """
    Print all available parameters
    """
    print("=== Bond Parameters ===")
    for bond_type, (r0, kb) in BOND_PARAMS.items():
        print(f"{bond_type:8s}: r0={r0:.2f} Å, kb={kb:.1f} kcal/(mol·Å²)")
    
    print("\n=== Angle Parameters ===")
    for angle_type, (theta, ka) in ANGLE_PARAMS.items():
        print(f"{angle_type:10s}: θ0={theta:.1f}°, ka={ka:.1f} kcal/(mol·rad²)")
    
    print("\n=== VDW Parameters ===")
    for element, (sigma, epsilon) in VDW_PARAMS.items():
        print(f"{element:2s}: σ={sigma:.2f} Å, ε={epsilon:.3f} kcal/mol")