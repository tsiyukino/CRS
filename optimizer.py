# optimizer.py

from energy import calculate_total_energy

# Optimization constants
DEFAULT_MAX_STEPS = 1000
DEFAULT_STEP_SIZE = 0.001
DEFAULT_FORCE_THRESHOLD = 0.01
DEFAULT_ENERGY_THRESHOLD = 0.001
MIN_STEP_SIZE = 1e-8

def optimize_geometry(molecule, max_steps=DEFAULT_MAX_STEPS, step_size=DEFAULT_STEP_SIZE, 
                     force_threshold=DEFAULT_FORCE_THRESHOLD, 
                     energy_threshold=DEFAULT_ENERGY_THRESHOLD,
                     include_vdw=True, print_every=100):
    """
    Optimize molecular geometry using steepest descent with adaptive step size.
    
    Args:
        molecule: Molecule object to optimize
        max_steps: Maximum number of optimization steps
        step_size: Initial step size for position updates
        force_threshold: Convergence criterion for max force (kcal/mol·Å)
        energy_threshold: Convergence criterion for energy change (kcal/mol)
        include_vdw: Whether to include VDW interactions
        print_every: Print status every N steps (0 = no printing)
    
    Returns:
        tuple: (converged, final_energy, num_steps)
    """
    
    prev_energy = None
    current_step_size = step_size
    steps_since_energy_increase = 0
    
    for step in range(max_steps):
        # Calculate energy and forces
        total_energy, bond_energy, angle_energy, vdw_energy = \
            calculate_total_energy(molecule, include_vdw=include_vdw)
        
        # Get maximum force on any atom
        max_force = molecule.get_max_force()
        
        # Print progress
        if print_every > 0 and step % print_every == 0:
            print(f"Step {step}: E={total_energy:.4f}, MaxF={max_force:.4f}, StepSize={current_step_size:.6f}")
        
        # Check force convergence
        if max_force < force_threshold:
            if print_every > 0:
                print(f"\nConverged! Max force {max_force:.6f} < {force_threshold}")
                print(f"Final energy: {total_energy:.4f} kcal/mol")
                print(f"  Bond: {bond_energy:.4f}")
                print(f"  Angle: {angle_energy:.4f}")
                print(f"  VDW: {vdw_energy:.4f}")
            return True, total_energy, step
        
        # Check energy change convergence
        if prev_energy is not None:
            energy_change = abs(total_energy - prev_energy)
            if energy_change < energy_threshold and steps_since_energy_increase > 5:
                if print_every > 0:
                    print(f"\nConverged! Energy change {energy_change:.6f} < {energy_threshold}")
                    print(f"Final energy: {total_energy:.4f} kcal/mol")
                return True, total_energy, step
        
        # Adaptive step size: reduce if energy increased, increase if improving
        if prev_energy is not None:
            if total_energy > prev_energy:
                current_step_size *= 0.5  # Reduce step size
                steps_since_energy_increase = 0
                if print_every > 0 and step % print_every == 0:
                    print(f"  Energy increased, reducing step size to {current_step_size:.6f}")
            else:
                steps_since_energy_increase += 1
                # Gradually increase step size if making good progress
                if steps_since_energy_increase > 10 and current_step_size < step_size:
                    current_step_size = min(current_step_size * 1.1, step_size)
        
        # Check if step size became too small
        if current_step_size < MIN_STEP_SIZE:
            if print_every > 0:
                print(f"\nStopping: step size {current_step_size} < minimum {MIN_STEP_SIZE}")
                print(f"This may indicate the optimization is stuck.")
            return False, total_energy, step
        
        # Update positions
        molecule.update_all_positions(current_step_size)
        
        prev_energy = total_energy
    
    # Did not converge
    if print_every > 0:
        print(f"\nDid not converge after {max_steps} steps")
        print(f"Final energy: {total_energy:.4f}, Max force: {max_force:.4f}")
    
    return False, total_energy, max_steps

def optimize_with_damping(molecule, max_steps=DEFAULT_MAX_STEPS, 
                         step_size=0.01, damping=0.9, 
                         force_threshold=DEFAULT_FORCE_THRESHOLD,
                         include_vdw=True, print_every=100):
    """
    Optimize geometry using velocity damping (momentum method).
    Helps avoid oscillations and can converge faster than steepest descent.
    
    Args:
        molecule: Molecule object to optimize
        max_steps: Maximum optimization steps
        step_size: Step size for velocity updates
        damping: Velocity damping factor (0.9 = keep 90% of previous velocity)
        force_threshold: Convergence criterion for max force
        include_vdw: Whether to include VDW interactions
        print_every: Print frequency (0 for silent)
    
    Returns:
        tuple: (converged, final_energy, num_steps)
    """
    
    # Initialize velocities for all atoms
    velocities = []
    for atom in molecule.atoms:
        velocities.append({'vx': 0.0, 'vy': 0.0, 'vz': 0.0})
    
    prev_energy = None
    
    for step in range(max_steps):
        # Calculate energy and forces
        total_energy, bond_energy, angle_energy, vdw_energy = \
            calculate_total_energy(molecule, include_vdw=include_vdw)
        
        max_force = molecule.get_max_force()
        
        # Print progress
        if print_every > 0 and step % print_every == 0:
            print(f"Step {step}: E={total_energy:.4f}, MaxF={max_force:.4f}")
        
        # Check convergence
        if max_force < force_threshold:
            if print_every > 0:
                print(f"\nConverged! Max force {max_force:.6f} < {force_threshold}")
                print(f"Final energy: {total_energy:.4f} kcal/mol")
                print(f"  Bond: {bond_energy:.4f}")
                print(f"  Angle: {angle_energy:.4f}")
                print(f"  VDW: {vdw_energy:.4f}")
            return True, total_energy, step
        
        # Update velocities with damping and current forces
        for i, atom in enumerate(molecule.atoms):
            velocities[i]['vx'] = damping * velocities[i]['vx'] + step_size * atom.fx
            velocities[i]['vy'] = damping * velocities[i]['vy'] + step_size * atom.fy
            velocities[i]['vz'] = damping * velocities[i]['vz'] + step_size * atom.fz
        
        # Update positions using velocities
        for i, atom in enumerate(molecule.atoms):
            atom.x += velocities[i]['vx']
            atom.y += velocities[i]['vy']
            atom.z += velocities[i]['vz']
        
        prev_energy = total_energy
    
    if print_every > 0:
        print(f"\nDid not converge after {max_steps} steps")
        print(f"Final energy: {total_energy:.4f}, Max force: {max_force:.4f}")
    
    return False, total_energy, max_steps

def line_search_optimize(molecule, max_steps=500, initial_step=0.1,
                        force_threshold=DEFAULT_FORCE_THRESHOLD, 
                        include_vdw=True, print_every=50):
    """
    Line search optimization: finds optimal step size along force direction.
    More robust than simple steepest descent but slower per iteration.
    
    Args:
        molecule: Molecule object to optimize
        max_steps: Maximum optimization steps
        initial_step: Initial step size to try
        force_threshold: Convergence criterion for max force
        include_vdw: Whether to include VDW interactions
        print_every: Print frequency (0 for silent)
    
    Returns:
        tuple: (converged, final_energy, num_steps)
    """
    
    for step in range(max_steps):
        # Calculate energy and forces
        total_energy, bond_energy, angle_energy, vdw_energy = \
            calculate_total_energy(molecule, include_vdw=include_vdw)
        max_force = molecule.get_max_force()
        
        if print_every > 0 and step % print_every == 0:
            print(f"Step {step}: E={total_energy:.4f}, MaxF={max_force:.4f}")
        
        # Check convergence
        if max_force < force_threshold:
            if print_every > 0:
                print(f"\nConverged at step {step}")
                print(f"Final energy: {total_energy:.4f} kcal/mol")
                print(f"  Bond: {bond_energy:.4f}")
                print(f"  Angle: {angle_energy:.4f}")
                print(f"  VDW: {vdw_energy:.4f}")
            return True, total_energy, step
        
        # Store current positions and forces
        saved_positions = [(atom.x, atom.y, atom.z) for atom in molecule.atoms]
        saved_forces = [(atom.fx, atom.fy, atom.fz) for atom in molecule.atoms]
        
        # Try different step sizes
        best_step = initial_step
        best_energy = total_energy
        
        for test_step in [initial_step * 0.1, initial_step * 0.5, initial_step, initial_step * 2]:
            # Restore positions
            for i, atom in enumerate(molecule.atoms):
                atom.x, atom.y, atom.z = saved_positions[i]
            
            # Move with test step size
            for i, atom in enumerate(molecule.atoms):
                atom.x += test_step * saved_forces[i][0]
                atom.y += test_step * saved_forces[i][1]
                atom.z += test_step * saved_forces[i][2]
            
            # Calculate energy at new position
            test_energy, _, _, _ = calculate_total_energy(molecule, include_vdw=include_vdw)
            
            if test_energy < best_energy:
                best_energy = test_energy
                best_step = test_step
        
        # Restore and move with best step
        for i, atom in enumerate(molecule.atoms):
            atom.x, atom.y, atom.z = saved_positions[i]
            atom.x += best_step * saved_forces[i][0]
            atom.y += best_step * saved_forces[i][1]
            atom.z += best_step * saved_forces[i][2]
    
    if print_every > 0:
        print(f"\nDid not converge after {max_steps} steps")
        total_energy, bond_energy, angle_energy, vdw_energy = \
            calculate_total_energy(molecule, include_vdw=include_vdw)
        max_force = molecule.get_max_force()
        print(f"Final energy: {total_energy:.4f}, Max force: {max_force:.4f}")
    
    return False, total_energy, max_steps