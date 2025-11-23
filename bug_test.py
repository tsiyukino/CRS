#!/usr/bin/env python3
"""
Targeted tests to identify potential bugs in the molecular mechanics simulator
"""

import math
from molecule import Molecule, Atom
from geometry import calculate_angle, distance, normalize_vector
from energy import calculate_vdw_energy_and_force, calculate_total_energy
from optimizer import line_search_optimize

print("=" * 70)
print("BUG TEST SUITE")
print("=" * 70)

# Test 1: Division by zero in calculate_angle when atoms are coincident
print("\n=== Test 1: Coincident atoms in angle calculation ===")
try:
    atom1 = Atom('H', 0.0, 0.0, 0.0)
    atom2 = Atom('O', 0.0, 0.0, 0.0)  # Same position as atom1!
    atom3 = Atom('H', 1.0, 0.0, 0.0)
    angle = calculate_angle(atom1, atom2, atom3)
    print(f"Result: angle = {angle} radians")
    print("BUG: Should have raised an error or handled division by zero!")
except ZeroDivisionError as e:
    print(f"ERROR CAUGHT: {e}")
    print("BUG CONFIRMED: Division by zero not properly handled")

# Test 2: VDW force direction check
print("\n=== Test 2: VDW force direction ===")
mol = Molecule()
# Two atoms far apart (should attract)
mol.add_atom('C', 0.0, 0.0, 0.0)
mol.add_atom('C', 5.0, 0.0, 0.0)  # 5 Angstroms apart

mol.zero_all_forces()
vdw_energy = calculate_vdw_energy_and_force(mol, cutoff=10.0)

atom1 = mol.atoms[0]
atom2 = mol.atoms[1]

print(f"Distance: {distance(atom1, atom2):.2f} Å")
print(f"VDW Energy: {vdw_energy:.6f} kcal/mol")
print(f"Force on atom1: fx={atom1.fx:.6f}, fy={atom1.fy:.6f}, fz={atom1.fz:.6f}")
print(f"Force on atom2: fx={atom2.fx:.6f}, fy={atom2.fy:.6f}, fz={atom2.fz:.6f}")

# At this distance, atoms should attract
# Energy should be negative (attractive well)
# Force on atom1 should be positive in x (pulling toward atom2)
if vdw_energy < 0:
    print("✓ Energy is negative (attractive) - CORRECT")
else:
    print("✗ Energy is positive - WRONG (should be negative)")

if atom1.fx > 0:
    print("✓ Force on atom1 is positive x-direction (toward atom2) - CORRECT")
else:
    print("✗ Force on atom1 is negative - WRONG (should attract)")

# Test at very close distance (should repel strongly)
print("\n=== Test 2b: VDW at very close distance (repulsive) ===")
mol2 = Molecule()
mol2.add_atom('C', 0.0, 0.0, 0.0)
mol2.add_atom('C', 2.0, 0.0, 0.0)  # 2 Angstroms apart (very close)

mol2.zero_all_forces()
vdw_energy2 = calculate_vdw_energy_and_force(mol2, cutoff=10.0)

atom1 = mol2.atoms[0]
atom2 = mol2.atoms[1]

print(f"Distance: {distance(atom1, atom2):.2f} Å")
print(f"VDW Energy: {vdw_energy2:.6f} kcal/mol")
print(f"Force on atom1: fx={atom1.fx:.6f}, fy={atom1.fy:.6f}, fz={atom1.fz:.6f}")

# At close distance, should repel
# Energy should be very positive
# Force should push atoms apart (negative for atom1)
if vdw_energy2 > 0:
    print("✓ Energy is positive (repulsive) - CORRECT")
else:
    print("✗ Energy is negative - WRONG (should be repulsive)")

if atom1.fx < 0:
    print("✓ Force on atom1 is negative x-direction (pushing away) - CORRECT")
else:
    print("✗ Force on atom1 is positive - WRONG (should repel)")

# Test 3: line_search_optimize returning correct energy
print("\n=== Test 3: line_search_optimize energy return value ===")
mol3 = Molecule()
mol3.add_atom('H', 0.0, 0.0, 0.0)
mol3.add_atom('H', 2.0, 0.0, 0.0)
mol3.add_bond(0, 1, 0.74, 340.0)

# Get initial energy
initial_energy, _, _, _ = calculate_total_energy(mol3, include_vdw=False)
print(f"Initial energy: {initial_energy:.6f} kcal/mol")

# Run optimization
converged, final_energy_returned, steps = line_search_optimize(
    mol3, max_steps=100, force_threshold=0.01, include_vdw=False, print_every=0
)

# Calculate actual final energy
actual_final_energy, _, _, _ = calculate_total_energy(mol3, include_vdw=False)

print(f"Energy returned by optimizer: {final_energy_returned:.6f} kcal/mol")
print(f"Actual final energy: {actual_final_energy:.6f} kcal/mol")
print(f"Difference: {abs(final_energy_returned - actual_final_energy):.6f} kcal/mol")

if abs(final_energy_returned - actual_final_energy) < 0.0001:
    print("✓ Returned energy matches actual final energy - CORRECT")
else:
    print("✗ BUG: Returned energy doesn't match actual final energy!")
    print("   This suggests the optimizer is returning stale energy values")

# Test 4: Normalize vector with zero length
print("\n=== Test 4: Normalize zero-length vector ===")
result = normalize_vector(0.0, 0.0, 0.0)
print(f"normalize_vector(0, 0, 0) = {result}")
print("Note: Returns (0,0,0) which could cause issues if used for direction")

# Test 5: Check if angle force conserves momentum
print("\n=== Test 5: Force conservation check ===")
from simulator import create_molecule_from_atoms
atoms = [
    ('O', 0.0, 0.0, 0.0),
    ('H', 1.0, 0.0, 0.0),
    ('H', 0.0, 1.0, 0.0),
]
bonds = [(0, 1), (0, 2)]
angles = [(1, 0, 2)]

water = create_molecule_from_atoms(atoms, bonds, angles)
calculate_total_energy(water, include_vdw=False)

# Sum of all forces should be zero (Newton's 3rd law)
total_fx = sum(atom.fx for atom in water.atoms)
total_fy = sum(atom.fy for atom in water.atoms)
total_fz = sum(atom.fz for atom in water.atoms)

print(f"Total force: fx={total_fx:.10f}, fy={total_fy:.10f}, fz={total_fz:.10f}")

if abs(total_fx) < 1e-8 and abs(total_fy) < 1e-8 and abs(total_fz) < 1e-8:
    print("✓ Forces sum to zero - CORRECT (Newton's 3rd law holds)")
else:
    print("✗ WARNING: Forces don't sum to zero - potential issue")

print("\n" + "=" * 70)
print("BUG TEST SUITE COMPLETE")
print("=" * 70)
