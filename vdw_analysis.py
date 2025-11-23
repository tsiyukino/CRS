#!/usr/bin/env python3
"""
Detailed analysis of VDW force calculation
"""

import math

# Manual calculation of VDW force
def manual_vdw_force(r, sigma, epsilon):
    """
    Calculate VDW force manually using correct formula

    E_vdw = ε × [(σ/r)¹² - 2(σ/r)⁶]

    dE/dr = ε × [12(σ¹²)(-r⁻¹³) - 2×6(σ⁶)(-r⁻⁷)]
          = ε × [-12σ¹²/r¹³ + 12σ⁶/r⁷]

    F = -dE/dr = ε × [12σ¹²/r¹³ - 12σ⁶/r⁷]
              = 12ε/r × [(σ/r)¹² - (σ/r)⁶]
              = 24ε/r × [2(σ/r)¹² - (σ/r)⁶]
    """
    sr = sigma / r
    sr6 = sr ** 6
    sr12 = sr ** 12

    # Correct formula
    force = 24.0 * epsilon / r * (2.0 * sr12 - sr6)

    return force

# Code's calculation
def code_vdw_force(r, sigma, epsilon):
    """
    What the code calculates (from energy.py:181)
    """
    sr = sigma / r
    sr6 = sr ** 6
    sr12 = sr ** 12

    # From the code: force_magnitude = -24.0 * epsilon * (2.0 * sr12 / r - sr6 / r)
    force = -24.0 * epsilon * (2.0 * sr12 / r - sr6 / r)

    return force

# Test at different distances
sigma = 3.4  # Carbon
epsilon = 0.086

print("=" * 70)
print("VDW Force Analysis: Comparing Manual vs Code Implementation")
print("=" * 70)
print(f"\nParameters: σ={sigma} Å, ε={epsilon} kcal/mol")
print(f"\nMinimum at r = 2^(1/6) × σ = {2**(1/6) * sigma:.3f} Å")
print("\n" + "-" * 70)
print(f"{'Distance (Å)':<15} {'Manual F':<20} {'Code F':<20} {'Match?':<10}")
print("-" * 70)

for r in [2.0, 2.5, 3.0, 3.5, 3.8174, 4.0, 5.0, 6.0]:
    manual_f = manual_vdw_force(r, sigma, epsilon)
    code_f = code_vdw_force(r, sigma, epsilon)
    match = "YES" if abs(manual_f - code_f) < 1e-6 else "NO"
    print(f"{r:<15.2f} {manual_f:<20.6f} {code_f:<20.6f} {match:<10}")

print("\n" + "=" * 70)
print("Analysis:")
print("=" * 70)

# Check the sign
r_test = 5.0
manual = manual_vdw_force(r_test, sigma, epsilon)
code = code_vdw_force(r_test, sigma, epsilon)

print(f"\nAt r={r_test} Å (beyond minimum, attractive region):")
print(f"  Manual formula: F = {manual:.6f} (positive = attractive)")
print(f"  Code formula:   F = {code:.6f}")

if manual > 0 and code > 0:
    print("  ✓ Both give attractive force (correct)")
elif manual > 0 and code < 0:
    print("  ✗ BUG: Code gives wrong sign!")
elif abs(manual - code) > 1e-6:
    print("  ✗ BUG: Different magnitudes!")
else:
    print("  ✓ Formulas match")

r_test2 = 2.0
manual2 = manual_vdw_force(r_test2, sigma, epsilon)
code2 = code_vdw_force(r_test2, sigma, epsilon)

print(f"\nAt r={r_test2} Å (inside minimum, repulsive region):")
print(f"  Manual formula: F = {manual2:.6f} (negative = repulsive)")
print(f"  Code formula:   F = {code2:.6f}")

if manual2 < 0 and code2 < 0:
    print("  ✓ Both give repulsive force (correct)")
elif manual2 * code2 < 0:
    print("  ✗ BUG: Opposite signs!")
else:
    print("  ? Check needed")
