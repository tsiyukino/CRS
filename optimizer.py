# optimizer.py

from energy import calculate_total_energy

def optimize_geometry(molecule, max_steps=1000, step_size=0.001, 
                     force_threshold=0.01, energy_threshold=0.001,
                     include_vdw=True, print_every=100):
    """
    Optimize molecular geometry using steepest descent
    
    Args:
        molecule: Molecule object to optimize
        max_steps: Maximum number of optimization steps
        step_size: Initial step size for position updates
        force_threshold: Convergence criterion for max force (kcal/mol·Å)
        energy_threshold: Convergence criterion for energy change (kcal/mol)
        include_vdw: Whether to include VDW interactions
        print_every: Print status every N steps (0 = no printing)
    
    Returns:
        (converged, final_energy, num_steps)
    """
    
    prev_energy = None
    current_step_size = step_size
    
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
            if energy_change < energy_threshold:
                if print_every > 0:
                    print(f"\nConverged! Energy change {energy_change:.6f} < {energy_threshold}")
                    print(f"Final energy: {total_energy:.4f} kcal/mol")
                return True, total_energy, step
        
        # Adaptive step size: reduce if energy increased
        if prev_energy is not None and total_energy > prev_energy:
            current_step_size *= 0.5  # Reduce step size
            if print_every > 0 and step % print_every == 0:
                print(f"  Energy increased, reducing step size to {current_step_size:.6f}")
        
        # Update positions
        molecule.update_all_positions(current_step_size)
        
        prev_energy = total_energy
    
    # Did not converge
    if print_every > 0:
        print(f"\nDid not converge after {max_steps} steps")
        print(f"Final energy: {total_energy:.4f}, Max force: {max_force:.4f}")
    
    return False, total_energy, max_steps

def optimize_with_damping(molecule, max_steps=1000, step_size=0.01,
                         damping=0.9, force_threshold=0.01,
                         include_vdw=True, print_every=100):
    """
    Optimize geometry using velocity damping (momentum method)
    Helps avoid oscillations
    
    Args:
        damping: Velocity damping factor (0.9 = keep 90% of previous velocity)
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
    
    return False, total_energy, max_steps

def line_search_optimize(molecule, max_steps=500, initial_step=0.1,
                        force_threshold=0.01, include_vdw=True, 
                        print_every=50):
    """
    Simple line search optimization
    Finds optimal step size along force direction each iteration
    """
    
    for step in range(max_steps):
        # Calculate energy and forces
        total_energy, _, _, _ = calculate_total_energy(molecule, include_vdw=include_vdw)
        max_force = molecule.get_max_force()
        
        if print_every > 0 and step % print_every == 0:
            print(f"Step {step}: E={total_energy:.4f}, MaxF={max_force:.4f}")
        
        # Check convergence
        if max_force < force_threshold:
            if print_every > 0:
                print(f"\nConverged at step {step}")
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
    
    return False, total_energy, max_steps