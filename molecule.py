class Atom:
    def __init__(self, element, x, y, z):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.fx = 0.0  # force in x 
        self.fy = 0.0  # force in y 
        self.fz = 0.0  # force in z 
        
        # VDW parameters - set defaults based on element
        if element == 'H':
            self.sigma = 2.5
            self.epsilon = 0.015
        elif element == 'C':
            self.sigma = 3.4
            self.epsilon = 0.086
        elif element == 'O':
            self.sigma = 3.12
            self.epsilon = 0.120
        elif element == 'N':
            self.sigma = 3.25
            self.epsilon = 0.086
        else:
            self.sigma = 3.0
            self.epsilon = 0.05
    
    def zero_forces(self):
        # Reset forces to zero before calculating new forces
        self.fx = 0.0
        self.fy = 0.0
        self.fz = 0.0
    
    def get_position(self):
        # Return position
        return (self.x, self.y, self.z)
    
    def update_position(self, step_size):
        # Move atom based on current forces
        self.x += step_size * self.fx
        self.y += step_size * self.fy
        self.z += step_size * self.fz


class Bond:
    def __init__(self, atom1_index, atom2_index, r0, kb):
        self.atom1_index = atom1_index  # index of first atom in molecule atom list
        self.atom2_index = atom2_index  # index of second atom
        self.r0 = r0  # equilibrium bond length (Angstroms)
        self.kb = kb  # force constant (kcal/mol·Å²)


class Angle:
    def __init__(self, atom1_index, atom2_index, atom3_index, theta0, ka):
        self.atom1_index = atom1_index 
        self.atom2_index = atom2_index  
        self.atom3_index = atom3_index  
        self.theta0 = theta0  # equilibrium angle (radians)
        self.ka = ka  # force constant (kcal/mol·rad²)


class Molecule:
    def __init__(self):
        self.atoms = []   # Atom objects
        self.bonds = []   # Bond objects
        self.angles = []  # Angle objects
    
    def add_atom(self, element, x, y, z):
        # Add atom and return index
        atom = Atom(element, x, y, z)
        self.atoms.append(atom)
        return len(self.atoms) - 1  # return index of newly added atom
    
    def add_bond(self, atom1_index, atom2_index, r0, kb):
        # Add bond between two atoms
        bond = Bond(atom1_index, atom2_index, r0, kb)
        self.bonds.append(bond)
    
    def add_angle(self, atom1_index, atom2_index, atom3_index, theta0, ka):
        # Add angle (atom2 center)
        angle = Angle(atom1_index, atom2_index, atom3_index, theta0, ka)
        self.angles.append(angle)
    
    def zero_all_forces(self):
        # Reset forces on all atoms to zero
        for atom in self.atoms:
            atom.zero_forces()
    
    def get_max_force(self):
        # Return the maximum force component on any atom
        max_f = 0.0
        for atom in self.atoms:
            max_f = max(max_f, abs(atom.fx), abs(atom.fy), abs(atom.fz))
        return max_f
    
    def update_all_positions(self, step_size):
        # Move all atoms based on their forces
        for atom in self.atoms:
            atom.update_position(step_size)
    
    def print_geometry(self):
        # Print current atomic positions
        print("\nCurrent Geometry:")
        for i, atom in enumerate(self.atoms):
            print(f"Atom {i} ({atom.element}): x={atom.x:.4f}, y={atom.y:.4f}, z={atom.z:.4f}")