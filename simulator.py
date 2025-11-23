# simulator.py

import math
from molecule import Molecule
from parameters import get_bond_params, get_angle_params
from optimizer import optimize_geometry
from geometry import distance, calculate_angle

def create_molecule_from_atoms(atom_list, bond_list, angle_list=None):
    """
    Create molecule from lists of atoms, bonds, and angles
    
    Args:
        atom_list: [(element, x, y, z), ...]
        bond_list: [(atom1_idx, atom2_idx), ...] - gets params from database
        angle_list: [(atom1_idx, atom2_idx, atom3_idx), ...] or None
    
    Returns:
        Molecule object
    """
    mol = Molecule()
    
    # Add atoms
    for element, x, y, z in atom_list:
        mol.add_atom(element, x, y, z)
    
    # Add bonds with parameters from database
    for atom1_idx, atom2_idx in bond_list:
        elem1 = mol.atoms[atom1_idx].element
        elem2 = mol.atoms[atom2_idx].element
        
        params = get_bond_params(elem1, elem2)
        if params is None:
            print(f"Warning: No parameters for {elem1}-{elem2} bond, skipping")
            continue
        
        r0, kb = params
        mol.add_bond(atom1_idx, atom2_idx, r0, kb)
    
    # Add angles with parameters from database
    if angle_list is not None:
        for atom1_idx, atom2_idx, atom3_idx in angle_list:
            elem1 = mol.atoms[atom1_idx].element
            elem2 = mol.atoms[atom2_idx].element
            elem3 = mol.atoms[atom3_idx].element
            
            params = get_angle_params(elem1, elem2, elem3)
            if params is None:
                print(f"Warning: No parameters for {elem1}-{elem2}-{elem3} angle, skipping")
                continue
            
            theta0, ka = params
            mol.add_angle(atom1_idx, atom2_idx, atom3_idx, theta0, ka)
    
    return mol

def analyze_geometry(molecule, name="Molecule"):
    """
    Print analysis of optimized geometry
    """
    print(f"\n=== {name} Geometry Analysis ===")
    
    # Print all bond lengths
    print("\nBond Lengths:")
    for i, bond in enumerate(molecule.bonds):
        atom1 = molecule.atoms[bond.atom1_index]
        atom2 = molecule.atoms[bond.atom2_index]
        d = distance(atom1, atom2)
        print(f"  {atom1.element}-{atom2.element} ({bond.atom1_index}-{bond.atom2_index}): {d:.4f} Å (equilibrium: {bond.r0:.4f} Å)")
    
    # Print all angles
    if len(molecule.angles) > 0:
        print("\nAngles:")
        for i, angle in enumerate(molecule.angles):
            atom1 = molecule.atoms[angle.atom1_index]
            atom2 = molecule.atoms[angle.atom2_index]
            atom3 = molecule.atoms[angle.atom3_index]
            theta = calculate_angle(atom1, atom2, atom3)
            theta_deg = math.degrees(theta)
            theta0_deg = math.degrees(angle.theta0)
            print(f"  {atom1.element}-{atom2.element}-{atom3.element} ({angle.atom1_index}-{angle.atom2_index}-{angle.atom3_index}): {theta_deg:.2f}° (equilibrium: {theta0_deg:.2f}°)")

def run_simulation(molecule, name="Molecule", step_size=0.001, max_steps=2000, include_vdw=False, print_every=500):
    """
    Run complete simulation: optimize and analyze
    
    Args:
        molecule: Molecule object to simulate
        name: Name for printing
        step_size: Optimization step size
        max_steps: Maximum optimization steps
        include_vdw: Include Van der Waals interactions
        print_every: Print frequency (0 for silent)
    
    Returns:
        (converged, final_energy, steps)
    """
    print(f"\n{'='*60}")
    print(f"Simulating {name}")
    print(f"{'='*60}")
    
    print("\nInitial geometry:")
    molecule.print_geometry()
    
    print(f"\nOptimizing (step_size={step_size}, max_steps={max_steps}, VDW={include_vdw})...")
    converged, final_energy, steps = optimize_geometry(
        molecule,
        step_size=step_size,
        max_steps=max_steps,
        include_vdw=include_vdw,
        print_every=print_every
    )
    
    print("\nFinal geometry:")
    molecule.print_geometry()
    
    analyze_geometry(molecule, name)
    
    print(f"\nSummary:")
    print(f"  Converged: {converged}")
    print(f"  Steps: {steps}")
    print(f"  Final energy: {final_energy:.4f} kcal/mol")
    
    return converged, final_energy, steps