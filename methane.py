from simulator import create_molecule_from_atoms, run_simulation
from auto_angles import find_angles

# Generated with wedge-and-dash notation
# Wedge bonds: coming out of plane (+z)
# Dash bonds: going behind plane (-z)
# Normal bonds: in plane (z=0)

atoms = [
    ('C', 0.000, 0.000, 0.000),  # 0
    ('H', -1.301, 0.746, 0.800),  # 1
    ('H', -0.044, -1.499, 0.000),  # 2
    ('C', 1.306, 0.738, 0.800),  # 3
    ('C', 2.639, 0.050, 0.000),  # 4
    ('H', 0.697, 2.108, 1.600),  # 5
    ('H', 2.024, 2.054, 0.000),  # 6
    ('H', 1.470, -0.298, -0.800),  # 7
    ('H', 2.639, -1.450, 0.000),  # 8
    ('H', 3.914, 0.841, 0.800),  # 9
    ('H', 4.111, -0.240, -0.800),  # 10
]

bonds = [
    (0, 2),
    (0, 1),  # wedge
    (0, 3),  # wedge
    (3, 4),  # dash
    (3, 5),  # wedge
    (3, 6),  # dash
    (0, 7),  # dash
    (4, 10),  # dash
    (4, 9),  # wedge
    (4, 8),
]

angles = find_angles(atoms, bonds)
mol = create_molecule_from_atoms(atoms, bonds, angles)

# Optimize geometry
converged, energy, steps = run_simulation(
    mol, "MyMolecule",
    step_size=0.001,
    max_steps=2000,
    include_vdw=False,
    print_every=500
)

print(f"\nFinal: Energy={energy:.4f} kcal/mol, Steps={steps}")
