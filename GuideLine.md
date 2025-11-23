# Building a Molecular Mechanics Simulator from Scratch

## Part 1: Theory Foundation

### What are we doing?
We treat molecules as balls (atoms) connected by springs (bonds). We calculate the total energy based on:
1. How stretched/compressed bonds are
2. How bent angles are
3. How close non-bonded atoms are

Then we move atoms to minimize this energy, finding the stable molecular geometry.

### Why does this work?
Molecules naturally move toward their lowest energy state. By calculating energy and following the energy gradient downhill, we simulate this process.

---

## Part 2: Required Formulas

### 2.1 Bond Stretching Energy
When a bond is stretched or compressed from its natural length:

```
E_bond = (1/2) × k_b × (r - r_0)²
```

Where:
- `r` = current bond length
- `r_0` = equilibrium bond length (natural length)
- `k_b` = bond force constant (bond stiffness)

**Typical values:**
- C-C single bond: r_0 = 1.54 Å, k_b = 300 kcal/(mol·Å²)
- C=C double bond: r_0 = 1.34 Å, k_b = 600 kcal/(mol·Å²)
- C-H bond: r_0 = 1.09 Å, k_b = 340 kcal/(mol·Å²)
- O-H bond: r_0 = 0.96 Å, k_b = 450 kcal/(mol·Å²)

**Force from this energy:**
```
F = -dE/dr = -k_b × (r - r_0)
```

---

### 2.2 Angle Bending Energy
When a bond angle deviates from its natural angle:

```
E_angle = (1/2) × k_a × (θ - θ_0)²
```

Where:
- `θ` = current angle in radians
- `θ_0` = equilibrium angle
- `k_a` = angle force constant

**Typical values:**
- C-C-C: θ_0 = 109.5° (tetrahedral), k_a = 60 kcal/(mol·rad²)
- H-O-H: θ_0 = 104.5°, k_a = 55 kcal/(mol·rad²)
- C=C-C: θ_0 = 120° (trigonal), k_a = 70 kcal/(mol·rad²)

---

### 2.3 Van der Waals (Non-bonded) Energy
Between atoms not directly bonded (Lennard-Jones potential):

```
E_vdw = ε × [(σ/r)¹² - 2(σ/r)⁶]
```

Where:
- `r` = distance between atoms
- `ε` = well depth (how strongly they attract)
- `σ` = equilibrium distance

**Typical values (for carbon-carbon):**
- ε = 0.1 kcal/mol
- σ = 3.4 Å

**Combining rules for different atoms:**
```
ε_ij = sqrt(ε_i × ε_j)
σ_ij = (σ_i + σ_j) / 2
```

---

### 2.4 Distance Calculation
Between two atoms at positions (x₁,y₁,z₁) and (x₂,y₂,z₂):

```
r = sqrt((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²)
```

---

### 2.5 Angle Calculation
For three atoms A-B-C (B is the central atom):

```
Vector BA = (x_a - x_b, y_a - y_b, z_a - z_b)
Vector BC = (x_c - x_b, y_c - y_b, z_c - z_b)

cos(θ) = (BA · BC) / (|BA| × |BC|)
θ = arccos(cos(θ))
```

Dot product: BA · BC = BA_x×BC_x + BA_y×BC_y + BA_z×BC_z

---

## Part 3: Algorithm

### Step 1: Data Structure Setup

```
For each atom, store:
- Element type (C, H, O, etc.)
- Position (x, y, z)
- Force (fx, fy, fz) - starts at zero each step
- Van der Waals parameters (ε, σ)

For each bond, store:
- Atom 1 index
- Atom 2 index
- Equilibrium length r_0
- Force constant k_b

For each angle, store:
- Atom 1 index (first)
- Atom 2 index (center)
- Atom 3 index (last)
- Equilibrium angle θ_0
- Force constant k_a
```

### Step 2: Energy Minimization Loop

```
Repeat until converged:
    1. Set all atom forces to zero
    
    2. For each bond:
        - Calculate current length r
        - Calculate energy contribution
        - Calculate force on each atom
        - Add force to atoms
    
    3. For each angle:
        - Calculate current angle θ
        - Calculate energy contribution
        - Calculate force on each atom
        - Add force to atoms
    
    4. For each pair of non-bonded atoms:
        - Calculate distance r
        - If r < cutoff (typically 10 Å):
            - Calculate VDW energy
            - Calculate force on each atom
            - Add force to atoms
    
    5. Move atoms:
        For each atom:
            x_new = x_old + step_size × fx
            y_new = y_old + step_size × fy
            z_new = z_old + step_size × fz
    
    6. Check convergence:
        If max_force < threshold (e.g., 0.01):
            Stop - we found the minimum
```

### Step 3: Force Calculations (Detailed)

**Bond force:**
```
r = distance(atom1, atom2)
F_magnitude = -k_b × (r - r_0)

direction_x = (x2 - x1) / r
direction_y = (y2 - y1) / r
direction_z = (z2 - z1) / r

atom1.fx += F_magnitude × direction_x
atom1.fy += F_magnitude × direction_y
atom1.fz += F_magnitude × direction_z

atom2.fx -= F_magnitude × direction_x  (opposite direction)
atom2.fy -= F_magnitude × direction_y
atom2.fz -= F_magnitude × direction_z
```

**Angle force (simplified):**
For angle A-B-C, the force pushes/pulls atoms to change the angle. This is complex, so start by:
1. Calculate current angle θ
2. Calculate desired change: Δθ = θ - θ_0
3. Apply small forces perpendicular to bonds to change angle

(Full derivation requires vector calculus - can use numerical gradients instead for simplicity)

**VDW force:**
```
r = distance(atom1, atom2)
F_magnitude = 24 × ε × [(2×σ¹²/r¹³) - (σ⁶/r⁷)]

direction_x = (x2 - x1) / r
direction_y = (y2 - y1) / r
direction_z = (z2 - z1) / r

atom1.fx += F_magnitude × direction_x
atom1.fy += F_magnitude × direction_y
atom1.fz += F_magnitude × direction_z

(opposite on atom2)
```

---

## Part 4: Implementation Strategy

### Phase 1: Bond-only simulator (simplest)
1. Create atom and bond classes
2. Implement distance calculation
3. Implement bond energy and force
4. Implement simple gradient descent
5. Test on H₂ molecule (2 atoms, 1 bond)

**Expected result:** Bond length converges to r_0

### Phase 2: Add angles
1. Implement angle calculation (dot product)
2. Implement angle energy
3. Implement angle forces (start with numerical gradient)
4. Test on H₂O molecule (3 atoms, 2 bonds, 1 angle)

**Expected result:** H-O-H angle converges to ~104.5°

### Phase 3: Add VDW
1. Implement Lennard-Jones potential
2. Add all-pairs check (with cutoff)
3. Test on multiple molecules or noble gases

**Expected result:** Molecules maintain reasonable distances

### Phase 4: Test real molecules
- Methane (CH₄): should be tetrahedral
- Ethane (C₂H₆): should have staggered conformation
- Benzene (C₆H₆): should be planar hexagon

---

## Part 5: Numerical Methods

### Gradient Descent
```
step_size = 0.01 (start small!)

For each atom:
    x += step_size × fx
    y += step_size × fy
    z += step_size × fz
```

**Problems:**
- Can overshoot
- Can oscillate
- Slow convergence

**Solutions:**
- Adaptive step size: reduce if energy increases
- Add damping: velocity = 0.9 × velocity + force
- Use better optimizers (steepest descent, conjugate gradient)

### Convergence Criteria
```
max_force = max(|fx|, |fy|, |fz|) for all atoms

If max_force < 0.01 kcal/(mol·Å):
    Converged!
```

Or check energy change:
```
If |E_new - E_old| < 0.001 kcal/mol:
    Converged!
```

---

## Part 6: Common Issues and Solutions

### Issue 1: Atoms fly apart
- Step size too large → reduce to 0.001
- Forces too strong → check force constant units
- No attractive forces → add VDW

### Issue 2: Never converges
- Stuck in local minimum → try different starting geometry
- Oscillating → add damping or reduce step size
- Angles not working → check angle force implementation

### Issue 3: Wrong geometry
- Wrong equilibrium values → check r_0, θ_0 from literature
- Missing angle terms → make sure all angles are defined
- VDW too strong/weak → adjust ε values

---

## Part 7: Parameter Tables

### Bond Parameters
| Bond | r_0 (Å) | k_b (kcal/mol·Å²) |
|------|---------|-------------------|
| C-C  | 1.54    | 310              |
| C=C  | 1.34    | 600              |
| C-H  | 1.09    | 340              |
| C-O  | 1.43    | 320              |
| C=O  | 1.23    | 570              |
| O-H  | 0.96    | 450              |
| N-H  | 1.01    | 410              |

### Angle Parameters
| Angle   | θ_0 (degrees) | k_a (kcal/mol·rad²) |
|---------|---------------|---------------------|
| C-C-C   | 109.5         | 60                  |
| C-C-H   | 109.5         | 50                  |
| H-C-H   | 109.5         | 35                  |
| H-O-H   | 104.5         | 55                  |
| C-O-H   | 108.5         | 55                  |
| O-C-O   | 120.0         | 80                  |

### VDW Parameters
| Atom | σ (Å) | ε (kcal/mol) |
|------|-------|--------------|
| H    | 2.5   | 0.015        |
| C    | 3.4   | 0.086        |
| N    | 3.25  | 0.086        |
| O    | 3.12  | 0.120        |

---

## Part 8: Testing Strategy

### Test 1: H₂ molecule
Input: 2 H atoms at distance 1.0 Å
Expected: Converge to r ≈ 0.74 Å

### Test 2: Water (H₂O)
Input: O at origin, H atoms at (1,0,0) and (0,1,0)
Expected: 
- O-H distances ≈ 0.96 Å
- H-O-H angle ≈ 104.5°

### Test 3: Methane (CH₄)
Input: C at origin, 4 H atoms randomly around it
Expected:
- All C-H distances ≈ 1.09 Å
- All H-C-H angles ≈ 109.5° (tetrahedral)

---

## Part 9: Extensions (After Basic Version Works)

1. **Dihedral angles** - for rotation around bonds
2. **Partial charges** - for electrostatics (Coulomb)
3. **Molecular dynamics** - add velocities and simulate motion over time
4. **Bond breaking** - gradually reduce k_b to zero
5. **Constraints** - fix certain atoms in place
6. **Different force fields** - AMBER, CHARMM, etc.

---

## Part 10: Key Equations Summary

```
Total Energy = Σ E_bonds + Σ E_angles + Σ E_vdw

E_bond = (1/2) × k_b × (r - r_0)²

E_angle = (1/2) × k_a × (θ - θ_0)²

E_vdw = ε × [(σ/r)¹² - 2(σ/r)⁶]

F = -dE/dr

Update: position_new = position_old + step_size × force
```

---

This is everything you need to build a working molecular mechanics simulator. Start with Phase 1 (bonds only), get it working perfectly, then add complexity step by step.