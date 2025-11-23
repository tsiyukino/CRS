# test_simulator.py

from simulator import create_molecule_from_atoms, analyze_geometry, run_simulation

# Test 1: Create water molecule
print("=== Test 1: Create Water ===")
atoms = [
    ('O', 0.0, 0.0, 0.0),
    ('H', 1.2, 0.0, 0.0),
    ('H', 0.0, 1.2, 0.0),
]
bonds = [(0, 1), (0, 2)]
angles = [(1, 0, 2)]

water = create_molecule_from_atoms(atoms, bonds, angles)
print(f"Created molecule with {len(water.atoms)} atoms")
print(f"Number of bonds: {len(water.bonds)}")
print(f"Number of angles: {len(water.angles)}")
print(f"Expected: 3 atoms, 2 bonds, 1 angle")

# Test 2: Run water simulation
print("\n=== Test 2: Run Water Simulation ===")
converged, energy, steps = run_simulation(water, "Water", step_size=0.001, max_steps=1000, print_every=200)
print(f"\nExpected: O-H bonds ~0.96 Å, H-O-H angle ~104.5°")

# Test 3: Create methane
print("\n=== Test 3: Create Methane ===")
atoms = [
    ('C', 0.0, 0.0, 0.0),
    ('H', 1.2, 0.0, 0.0),
    ('H', 0.0, 1.2, 0.0),
    ('H', 0.0, 0.0, 1.2),
    ('H', 0.8, 0.8, 0.8),
]
bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
angles = [(1, 0, 2), (1, 0, 3), (1, 0, 4), (2, 0, 3), (2, 0, 4), (3, 0, 4)]

methane = create_molecule_from_atoms(atoms, bonds, angles)
print(f"Created molecule with {len(methane.atoms)} atoms")
print(f"Number of bonds: {len(methane.bonds)}")
print(f"Number of angles: {len(methane.angles)}")

# Test 4: Run methane simulation
print("\n=== Test 4: Run Methane Simulation ===")
converged, energy, steps = run_simulation(methane, "Methane", step_size=0.001, max_steps=2000, print_every=500)
print(f"\nExpected: All C-H bonds ~1.09 Å, all H-C-H angles ~109.5°")

# Test 5: Create H2 (simple)
print("\n=== Test 5: Simple H2 ===")
atoms = [('H', 0.0, 0.0, 0.0), ('H', 1.0, 0.0, 0.0)]
bonds = [(0, 1)]
h2 = create_molecule_from_atoms(atoms, bonds, None)
converged, energy, steps = run_simulation(h2, "H2", step_size=0.001, max_steps=500, print_every=100)
print(f"\nExpected: H-H bond ~0.74 Å")

# Test 6: Ammonia
print("\n=== Test 6: Ammonia ===")
atoms = [
    ('N', 0.0, 0.0, 0.0),
    ('H', 1.0, 0.0, 0.0),
    ('H', 0.0, 1.0, 0.0),
    ('H', 0.0, 0.0, 1.0),
]
bonds = [(0, 1), (0, 2), (0, 3)]
angles = [(1, 0, 2), (1, 0, 3), (2, 0, 3)]

ammonia = create_molecule_from_atoms(atoms, bonds, angles)
converged, energy, steps = run_simulation(ammonia, "Ammonia", step_size=0.001, max_steps=1000, print_every=200)
print(f"\nExpected: N-H bonds ~1.01 Å, H-N-H angles ~106.7°")

# Test 7: Missing parameters warning
print("\n=== Test 7: Missing Parameters ===")
atoms = [('C', 0.0, 0.0, 0.0), ('F', 1.5, 0.0, 0.0)]
bonds = [(0, 1)]
try:
    mol = create_molecule_from_atoms(atoms, bonds, None)
    print(f"Bonds created: {len(mol.bonds)}")
    print(f"Expected: Warning message and 0 bonds (C-F not in database)")
except:
    print("Error occurred")

# Test 8: Silent simulation
print("\n=== Test 8: Silent Simulation ===")
atoms = [('H', 0.0, 0.0, 0.0), ('H', 1.5, 0.0, 0.0)]
bonds = [(0, 1)]
h2_silent = create_molecule_from_atoms(atoms, bonds, None)
converged, energy, steps = run_simulation(h2_silent, "H2", step_size=0.001, max_steps=500, print_every=0)
print(f"Ran silently: converged={converged}, steps={steps}")

# Test 9: Just analyze without running
print("\n=== Test 9: Analyze Only ===")
atoms = [('O', 0.0, 0.0, 0.0), ('H', 0.96, 0.0, 0.0), ('H', -0.23, 0.94, 0.0)]
bonds = [(0, 1), (0, 2)]
angles = [(1, 0, 2)]
perfect_water = create_molecule_from_atoms(atoms, bonds, angles)
analyze_geometry(perfect_water, "Perfect Water")
print("Expected: Bonds and angles very close to equilibrium (no optimization ran)")

print("\n=== All Tests Complete ===")