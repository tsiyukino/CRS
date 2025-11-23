# energy.py

import math
from geometry import distance, calculate_angle, normalize_vector

def calculate_bond_energy_and_force(molecule):
    """
    Calculate bond stretching energy and add forces to atoms
    Returns total bond energy
    """
    total_energy = 0.0
    
    for bond in molecule.bonds:
        atom1 = molecule.atoms[bond.atom1_index]
        atom2 = molecule.atoms[bond.atom2_index]
        
        # Current bond length
        r = distance(atom1, atom2)
        
        # Energy: E = 0.5 * kb * (r - r0)^2
        dr = r - bond.r0
        energy = 0.5 * bond.kb * dr * dr
        total_energy += energy
        
        # Force magnitude: F = -kb * (r - r0)
        force_magnitude = bond.kb * dr
        
        # Direction from atom1 to atom2
        dx = atom2.x - atom1.x
        dy = atom2.y - atom1.y
        dz = atom2.z - atom1.z
        
        # Normalize direction
        dir_x, dir_y, dir_z = normalize_vector(dx, dy, dz)
        
        # Apply force to atom1 (in direction of atom2)
        atom1.fx += force_magnitude * dir_x
        atom1.fy += force_magnitude * dir_y
        atom1.fz += force_magnitude * dir_z
        
        # Apply opposite force to atom2
        atom2.fx -= force_magnitude * dir_x
        atom2.fy -= force_magnitude * dir_y
        atom2.fz -= force_magnitude * dir_z
    
    return total_energy

def calculate_angle_energy_and_force(molecule):
    """
    Calculate angle bending energy and add forces to atoms
    Returns total angle energy
    """
    total_energy = 0.0
    
    for angle in molecule.angles:
        atom1 = molecule.atoms[angle.atom1_index]
        atom2 = molecule.atoms[angle.atom2_index]  # center atom
        atom3 = molecule.atoms[angle.atom3_index]
        
        # Current angle
        theta = calculate_angle(atom1, atom2, atom3)
        
        # Energy: E = 0.5 * ka * (theta - theta0)^2
        dtheta = theta - angle.theta0
        energy = 0.5 * angle.ka * dtheta * dtheta
        total_energy += energy
        
        # Simplified angle force using numerical gradient
        # Calculate force by seeing how energy changes with small movements
        delta = 0.0001  # Small displacement
        
        # Save current positions
        x1, y1, z1 = atom1.x, atom1.y, atom1.z
        x2, y2, z2 = atom2.x, atom2.y, atom2.z
        x3, y3, z3 = atom3.x, atom3.y, atom3.z
        
        # Numerical gradient for atom1
        atom1.x += delta
        theta_plus = calculate_angle(atom1, atom2, atom3)
        atom1.x = x1 - delta
        theta_minus = calculate_angle(atom1, atom2, atom3)
        atom1.x = x1
        dtheta_dx1 = (theta_plus - theta_minus) / (2 * delta)
        fx1 = -angle.ka * dtheta * dtheta_dx1
        
        atom1.y += delta
        theta_plus = calculate_angle(atom1, atom2, atom3)
        atom1.y = y1 - delta
        theta_minus = calculate_angle(atom1, atom2, atom3)
        atom1.y = y1
        dtheta_dy1 = (theta_plus - theta_minus) / (2 * delta)
        fy1 = -angle.ka * dtheta * dtheta_dy1
        
        atom1.z += delta
        theta_plus = calculate_angle(atom1, atom2, atom3)
        atom1.z = z1 - delta
        theta_minus = calculate_angle(atom1, atom2, atom3)
        atom1.z = z1
        dtheta_dz1 = (theta_plus - theta_minus) / (2 * delta)
        fz1 = -angle.ka * dtheta * dtheta_dz1
        
        # Numerical gradient for atom3
        atom3.x += delta
        theta_plus = calculate_angle(atom1, atom2, atom3)
        atom3.x = x3 - delta
        theta_minus = calculate_angle(atom1, atom2, atom3)
        atom3.x = x3
        dtheta_dx3 = (theta_plus - theta_minus) / (2 * delta)
        fx3 = -angle.ka * dtheta * dtheta_dx3
        
        atom3.y += delta
        theta_plus = calculate_angle(atom1, atom2, atom3)
        atom3.y = y3 - delta
        theta_minus = calculate_angle(atom1, atom2, atom3)
        atom3.y = y3
        dtheta_dy3 = (theta_plus - theta_minus) / (2 * delta)
        fy3 = -angle.ka * dtheta * dtheta_dy3
        
        atom3.z += delta
        theta_plus = calculate_angle(atom1, atom2, atom3)
        atom3.z = z3 - delta
        theta_minus = calculate_angle(atom1, atom2, atom3)
        atom3.z = z3
        dtheta_dz3 = (theta_plus - theta_minus) / (2 * delta)
        fz3 = -angle.ka * dtheta * dtheta_dz3
        
        # Apply forces
        atom1.fx += fx1
        atom1.fy += fy1
        atom1.fz += fz1
        
        atom3.fx += fx3
        atom3.fy += fy3
        atom3.fz += fz3
        
        # Equal and opposite on center atom (Newton's 3rd law)
        atom2.fx -= (fx1 + fx3)
        atom2.fy -= (fy1 + fy3)
        atom2.fz -= (fz1 + fz3)
    
    return total_energy

def calculate_vdw_energy_and_force(molecule, cutoff=10.0):
    """
    Calculate Van der Waals energy between non-bonded atoms
    Returns total VDW energy
    cutoff: ignore pairs beyond this distance (Angstroms)
    """
    total_energy = 0.0
    
    # Check all pairs of atoms
    for i in range(len(molecule.atoms)):
        for j in range(i + 1, len(molecule.atoms)):
            # Skip if atoms are bonded
            if is_bonded(molecule, i, j):
                continue
            
            atom1 = molecule.atoms[i]
            atom2 = molecule.atoms[j]
            
            # Distance between atoms
            r = distance(atom1, atom2)
            
            # Skip if beyond cutoff
            if r > cutoff:
                continue
            
            # Combine rules for VDW parameters
            epsilon = math.sqrt(atom1.epsilon * atom2.epsilon)
            sigma = (atom1.sigma + atom2.sigma) / 2.0
            
            # Lennard-Jones: E = epsilon * [(sigma/r)^12 - 2*(sigma/r)^6]
            sr = sigma / r
            sr6 = sr * sr * sr * sr * sr * sr
            sr12 = sr6 * sr6
            
            energy = epsilon * (sr12 - 2.0 * sr6)
            total_energy += energy
            
            # Force: F = -dE/dr = -24*epsilon * [(2*sigma^12/r^13) - (sigma^6/r^7)]
            force_magnitude = -24.0 * epsilon * (2.0 * sr12 / r - sr6 / r)
            
            # Direction from atom1 to atom2
            dx = atom2.x - atom1.x
            dy = atom2.y - atom1.y
            dz = atom2.z - atom1.z
            
            # Normalize direction
            dir_x, dir_y, dir_z = normalize_vector(dx, dy, dz)
            
            # Apply force
            atom1.fx += force_magnitude * dir_x
            atom1.fy += force_magnitude * dir_y
            atom1.fz += force_magnitude * dir_z
            
            atom2.fx -= force_magnitude * dir_x
            atom2.fy -= force_magnitude * dir_y
            atom2.fz -= force_magnitude * dir_z
    
    return total_energy

def is_bonded(molecule, atom1_index, atom2_index):
    # Check if two atoms are directly bonded
    for bond in molecule.bonds:
        if (bond.atom1_index == atom1_index and bond.atom2_index == atom2_index) or \
           (bond.atom1_index == atom2_index and bond.atom2_index == atom1_index):
            return True
    return False

def calculate_total_energy(molecule, include_vdw=True, vdw_cutoff=10.0):
    ##Calculate total energy and forces for entire molecule

    # Zero all forces
    molecule.zero_all_forces()
    
    bond_energy = calculate_bond_energy_and_force(molecule)
    angle_energy = calculate_angle_energy_and_force(molecule)
    
    vdw_energy = 0.0
    if include_vdw:
        vdw_energy = calculate_vdw_energy_and_force(molecule, vdw_cutoff)
    
    total = bond_energy + angle_energy + vdw_energy
    
    return total, bond_energy, angle_energy, vdw_energy