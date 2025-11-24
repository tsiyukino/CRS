# auto_angles.py - Automatically generate angles from bond topology

def find_angles(atoms, bonds):
    """
    Automatically find all valid angles from bond connectivity
    
    Args:
        atoms: List of (element, x, y, z) tuples
        bonds: List of (i, j) tuples (atom indices)
    
    Returns:
        List of (i, j, k) tuples representing angles (j is center)
    """
    
    # Build adjacency list
    adjacency = {}
    for i in range(len(atoms)):
        adjacency[i] = []
    
    for bond in bonds:
        i, j = bond[0], bond[1]
        adjacency[i].append(j)
        adjacency[j].append(i)
    
    angles = []
    
    # For each atom, find angles where it's the center
    for center in range(len(atoms)):
        neighbors = adjacency[center]
        
        # Need at least 2 neighbors to form angle
        if len(neighbors) < 2:
            continue
        
        # All pairs of neighbors form angles
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                atom1 = neighbors[i]
                atom2 = neighbors[j]
                angles.append((atom1, center, atom2))
    
    return angles


def export_molecule_code(atoms, bonds, molecule_name="Molecule"):
    """
    Generate complete Python code for molecule
    
    Args:
        atoms: List of (element, x, y, z) tuples
        bonds: List of (i, j) tuples
        molecule_name: Name for the molecule
    
    Returns:
        String containing Python code
    """
    
    angles = find_angles(atoms, bonds)
    
    code = f"""# {molecule_name} - Generated from Molecular Editor

from simulator import create_molecule_from_atoms, run_simulation

atoms = [
"""
    
    for i, atom in enumerate(atoms):
        elem, x, y, z = atom
        code += f"    ('{elem}', {x:.2f}, {y:.2f}, {z:.2f}),  # {i}\n"
    
    code += "]\n\nbonds = [\n"
    
    for bond in bonds:
        code += f"    ({bond[0]}, {bond[1]}),\n"
    
    code += "]\n\nangles = [\n"
    
    for angle in angles:
        code += f"    ({angle[0]}, {angle[1]}, {angle[2]}),\n"
    
    code += f"""]

# Create and optimize molecule
print("Creating {molecule_name}...")
mol = create_molecule_from_atoms(atoms, bonds, angles)
print(f"Atoms: {{len(mol.atoms)}}")
print(f"Bonds: {{len(mol.bonds)}}")
print(f"Angles: {{len(mol.angles)}}")

print("\\nOptimizing...")
converged, energy, steps = run_simulation(
    mol, "{molecule_name}", 
    step_size=0.001, 
    max_steps=2000,
    include_vdw=False,
    print_every=500
)

print(f"\\nFinal: Energy={{energy:.4f}} kcal/mol, Steps={{steps}}")

# Optional: Create 3D visualization
try:
    from visualize_3d import create_3d_html
    create_3d_html(mol, "{molecule_name}", energy, "{molecule_name.lower()}_3d.html")
    print(f"✓ 3D visualization saved to {molecule_name.lower()}_3d.html")
except ImportError:
    print("Note: visualize_3d.py not found, skipping 3D visualization")
"""
    
    return code


if __name__ == "__main__":
    # Example: Ethanol drawn by hand
    print("Example: Auto-generating angles for ethanol")
    print("="*60)
    
    # Atoms from drawing (2D coords converted to 3D)
    atoms = [
        ('C', 0.00, 0.00, 0.00),  # 0
        ('C', 1.54, 0.00, 0.00),  # 1
        ('O', 2.04, 1.26, 0.00),  # 2
        ('H', -0.40, -0.52, -0.90),  # 3
        ('H', -0.40, -0.52, 0.90),   # 4
        ('H', -0.40, 1.04, 0.00),    # 5
        ('H', 1.94, -0.52, -0.90),   # 6
        ('H', 1.94, -0.52, 0.90),    # 7
        ('H', 2.99, 1.26, 0.00),     # 8
    ]
    
    # Bonds from drawing
    bonds = [
        (0, 1),  # C-C
        (1, 2),  # C-O
        (0, 3), (0, 4), (0, 5),  # C-H
        (1, 6), (1, 7),  # C-H
        (2, 8),  # O-H
    ]
    
    # Auto-generate angles
    angles = find_angles(atoms, bonds)
    
    print(f"Found {len(angles)} angles:")
    for angle in angles:
        elem1 = atoms[angle[0]][0]
        elem2 = atoms[angle[1]][0]
        elem3 = atoms[angle[2]][0]
        print(f"  {elem1}{angle[0]}-{elem2}{angle[1]}-{elem3}{angle[2]}")
    
    print("\n" + "="*60)
    print("Generated Python code:")
    print("="*60)
    code = export_molecule_code(atoms, bonds, "Ethanol")
    print(code[:500] + "...")
    
    # Save to file
    with open('/mnt/user-data/outputs/generated_ethanol.py', 'w') as f:
        f.write(code)
    print("\n✓ Saved to generated_ethanol.py")
