import numpy as np

def CLT(materials, laminate_orientations, ply_thicknesses):
    total_thickness = sum(ply_thicknesses)
    Z = [-total_thickness / 2.0]
    for t in ply_thicknesses:
        Z.append(Z[-1] + t)

    # Initialize A matrix components
    A11, A12, A22, A66 = 0, 0, 0, 0

    for ply_idx in range(len(ply_thicknesses)):
        # Get material properties for current ply
        myE11, myE22, myG12, myNu12 = materials[ply_idx]
        
        # Orientation of current ply
        angle_deg = laminate_orientations[ply_idx]
        angle_rad = np.deg2rad(angle_deg)
        c = np.cos(angle_rad)
        s = np.sin(angle_rad)

        # Calculate reduced stiffnesses
        denom = myE11 - myNu12**2 * myE22
        myQ11 = (myE11**2) / denom
        myQ12 = (myNu12 * myE11 * myE22) / denom
        myQ22 = (myE11 * myE22) / denom
        myQ66 = myG12

        # Transform stiffnesses to laminate coordinates
        Q11_S = myQ11 * c**4 + 2 * (myQ12 + 2 * myQ66) * c**2 * s**2 + myQ22 * s**4
        Q12_S = (myQ11 + myQ22 - 4 * myQ66) * c**2 * s**2 + myQ12 * (c**4 + s**4)
        Q22_S = myQ11 * s**4 + 2 * (myQ12 + 2 * myQ66) * c**2 * s**2 + myQ22 * c**4
        Q66_S = (myQ11 + myQ22 - 2 * myQ12 - 2 * myQ66) * c**2 * s**2 + myQ66 * (c**4 + s**4)

        # Calculate thickness of current ply
        thickness = Z[ply_idx + 1] - Z[ply_idx]

        # Accumulate A matrix components
        A11 += Q11_S * thickness
        A12 += Q12_S * thickness
        A22 += Q22_S * thickness
        A66 += Q66_S * thickness

    return A11, A12, A22, A66


def smeared_engineering_constants(A11, A12, A22, A66, total_thickness):
    Ex_bar = A11 / total_thickness
    Ey_bar = A22 / total_thickness
    Gxy_bar = A66 / total_thickness
    nu_xy_bar = abs(A12 / A11)
    return Ex_bar, Ey_bar, Gxy_bar, nu_xy_bar


# Input parameters for hybrid laminate with individual ply thicknesses
# Material properties for each ply
# Example: 8 plies of CFRP and 2 plies of Metal
CFRP = (125744, 10030, 5555, 0.271)  # (E11, E22, G12, Nu12) for CFRP
Metal = (141000, 141000, 36000, 0.3)  # (E11, E22, G12, Nu12) for Metal

# Define materials per ply
materials = [CFRP,CFRP,CFRP,CFRP, Metal ,CFRP,CFRP,CFRP, CFRP]

# Ply orientations for each ply
laminate_orientations = [
    60, -60, 0, 0,                   # Orientations for CFRP plies
    0,                               # Orientations for Metal plies     
    60, -60, 0, 0                    # Orientations for CFRP plies                             
]

# Ply thicknesses for each ply (in mm)
ply_thicknesses = [
    0.125, 0.125, 0.125, 0.125,                              # Thicknesses for CFRP plies
    9.0,                                                     # Thicknesses for Metal plies
    0.125, 0.125, 0.125, 0.125                               # Thicknesses for CFRP plies
]

# Calculate total thickness
total_thickness = sum(ply_thicknesses)

# Calculate A matrix components
A11, A12, A22, A66 = CLT(materials, laminate_orientations, ply_thicknesses)

# Calculate smeared engineering constants
Ex_bar, Ey_bar, Gxy_bar, nu_xy_bar = smeared_engineering_constants(A11, A12, A22, A66, total_thickness)

# Output results
print(f"Smeared Engineering Constants for Hybrid Laminate with Individual Ply Thicknesses:")
print(f"Ex_bar = {Ex_bar:.2f} MPa")
print(f"Ey_bar = {Ey_bar:.2f} MPa")
print(f"Gxy_bar = {Gxy_bar:.2f} MPa")
print(f"nu_xy_bar = {nu_xy_bar:.4f}")
