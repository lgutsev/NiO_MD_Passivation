# Import necessary libraries
import re
import math

# Constants for unit conversion
EV_TO_KCAL = 23.0605
METER_TO_ANGSTROM = 1e10

def format_buckingham_parameters(params):
    """
    Format Buckingham parameters for LAMMPS.
    
    Args:
        params (dict): Dictionary of Buckingham parameters.
        
    Returns:
        list: Formatted Buckingham parameters for LAMMPS.
    """
    buckingham_lines = []
    for (atom_type1, atom_type2), values in params.items():
        A = values['A']  # eV
        rho = values['rho']  # Å no conversion necessary
        C = values['C']   #  eV * Å^6  no conversion necessary
        buckingham_line = f"pair_coeff {atom_type1} {atom_type2} buck/coul/long {A:.6e} {rho:.6e} {C:.6e} # {atom_type1}-{atom_type2}"
        buckingham_lines.append(buckingham_line)
    return buckingham_lines

def read_opls_aa_parameters(filepath):
    """
    Read OPLS-AA parameters from a file.
    
    Args:
        filepath (str): Path to the OPLS-AA parameter file.
        
    Returns:
        dict: Dictionary of OPLS-AA parameters.
    """
    opls_data = {
        'pair': [],
        'bond': [],
        'angle': [],
        'dihedral': [],
        'improper': [],
    }
    
    with open(filepath, 'r') as file:
        lines = file.readlines()
    
    section = None
    improper_start = False
    
    for line in lines:
        line = line.split('#')[0].strip()  # Remove comments and extra whitespace
        if not line:
            continue
        
        if line.startswith('Pair Coeffs'):
            section = 'pair'
            continue
        elif line.startswith('Bond Coeffs'):
            section = 'bond'
            continue
        elif line.startswith('Angle Coeffs'):
            section = 'angle'
            continue
        elif line.startswith('Dihedral Coeffs'):
            section = 'dihedral'
            continue
        elif line.startswith('Improper Coeffs'):
            section = 'improper'
            improper_start = True
            continue
        
        if section == 'pair':
            parts = line.split()
            if len(parts) >= 3:
                try:
                    atom_type = int(parts[0])
                    epsilon = float(parts[1])
                    sigma = float(parts[2])    
                    opls_data[section].append({'atom_type': atom_type, 'epsilon': epsilon, 'sigma': sigma})
                except ValueError as e:
                    print(f"Error parsing line: {line.strip()}, error: {e}")
        
        elif section == 'bond':
            parts = line.split()
            if len(parts) >= 3:
                try:
                    bond_type = int(parts[0])
                    k = float(parts[1])
                    r0 = float(parts[2])  # nm
                    opls_data[section].append((bond_type, k, r0))
                except ValueError as e:
                    print(f"Error parsing line: {line.strip()}, error: {e}")
        
        elif section == 'angle':
            parts = line.split()
            if len(parts) >= 3:
                try:
                    angle_type = int(parts[0])
                    k = float(parts[1])
                    theta0 = float(parts[2])  # degrees
                    opls_data[section].append((angle_type, k, theta0))
                except ValueError as e:
                    print(f"Error parsing line: {line.strip()}, error: {e}")
        
        elif section == 'dihedral':
            parts = line.split()
            if len(parts) >= 5:
                try:
                    dihedral_type = int(parts[0])
                    v1 = float(parts[1])
                    v2 = float(parts[2])
                    v3 = float(parts[3])
                    v4 = float(parts[4])
                    opls_data[section].append((dihedral_type, v1, v2, v3, v4))
                except ValueError as e:
                    print(f"Error parsing line: {line.strip()}, error: {e}")
        
        elif section == 'improper' and improper_start:
            parts = line.split()
            if len(parts) == 4:
                try:
                    improper_type = int(parts[0])
                    k = float(parts[1])
                    sign = int(parts[2])
                    multiplicity = int(parts[3])
                    opls_data[section].append((improper_type, k, sign, multiplicity))
                except ValueError as e:
                    print(f"Error parsing line: {line.strip()}, error: {e}")
            else:
                improper_start = False
    
    print("Finished parsing OPLS-AA parameters.")
    return opls_data
    
def format_lj_parameters(opls_data):
    """
    Format Lennard-Jones parameters for LAMMPS.
    
    Args:
        opls_data (dict): Dictionary of OPLS-AA parameters.
        
    Returns:
        str: Formatted Lennard-Jones parameters for LAMMPS.
    """
    lj_params = []
    for params in opls_data['pair']:
        lj_params.append(f"pair_coeff {params['atom_type']} {params['atom_type']} lj/cut/coul/long {params['epsilon']:.3f} {params['sigma']:.3f} # {params['atom_type']}-{params['atom_type']}")
    
    return '\n'.join(lj_params)

def format_bond_parameters(opls_data):
    """
    Format bond parameters for LAMMPS.
    
    Args:
        opls_data (dict): Dictionary of OPLS-AA parameters.
        
    Returns:
        str: Formatted bond parameters for LAMMPS.
    """
    bond_params = []
    for bond_type, k, r0 in opls_data['bond']:
        bond_params.append(f"bond_coeff {bond_type} {k:.3f} {r0:.3f} # bond_type {bond_type}")
    
    return '\n'.join(bond_params)

def format_angle_parameters(opls_data):
    """
    Format angle parameters for LAMMPS.
    
    Args:
        opls_data (dict): Dictionary of OPLS-AA parameters.
        
    Returns:
        str: Formatted angle parameters for LAMMPS.
    """
    angle_params = []
    for angle_type, k, theta0 in opls_data['angle']:
        angle_params.append(f"angle_coeff {angle_type} {k:.3f} {theta0:.3f} # angle_type {angle_type}")
    
    return '\n'.join(angle_params)

def format_dihedral_parameters(opls_data):
    """
    Format dihedral parameters for LAMMPS.
    
    Args:
        opls_data (dict): Dictionary of OPLS-AA parameters.
        
    Returns:
        str: Formatted dihedral parameters for LAMMPS.
    """
    dihedral_params = []
    for dihedral_type, v1, v2, v3, v4 in opls_data['dihedral']:
        dihedral_params.append(f"dihedral_coeff {dihedral_type} opls {v1:.3f} {v2:.3f} {v3:.3f} {v4:.3f} # dihedral_type {dihedral_type}")
    
    return '\n'.join(dihedral_params)

def format_improper_parameters(opls_data):
    """
    Format improper parameters for LAMMPS.
    
    Args:
        opls_data (dict): Dictionary of OPLS-AA parameters.
        
    Returns:
        str: Formatted improper parameters for LAMMPS.
    """
    improper_params = []
    for improper_type, k, sign, multiplicity in opls_data['improper']:
        improper_params.append(f"improper_coeff {improper_type} {k:.3f} {sign} {multiplicity} # improper_type {improper_type}")
    
    return '\n'.join(improper_params)
    

def apply_combining_rules(opls_data, manual_lj_dict, ligand_atoms, mixing_method='LB'):
    """
    Apply combining rules to generate pair coefficients for LAMMPS.

    Parameters:
    opls_data (dict): OPLS force field data.
    manual_lj_dict (dict): Manually defined Lennard-Jones parameters for surface atoms.
    ligand_atoms (list): List of ligand atom types.
    mixing_method (str): The mixing method to use ('LB' for Lorentz-Berthelot, 'geom' for geometric).

    Returns:
    list: List of combined pair coefficient lines for LAMMPS input.
    """
    import math
    combined_lines = []
    surface_atoms = ['Ni', 'O']

    for atom_type1 in surface_atoms:
        if (atom_type1, atom_type1) not in manual_lj_dict:
            continue
        
        epsilon1 = manual_lj_dict[(atom_type1, atom_type1)]['epsilon']
        sigma1 = manual_lj_dict[(atom_type1, atom_type1)]['sigma']
        
        for params in opls_data['pair']:
            atom_type2 = params['atom_type']
            if atom_type2 not in ligand_atoms:
                continue
            
            epsilon2 = params['epsilon']
            sigma2 = params['sigma']
                
            if mixing_method == 'LB':
                # Lorentz-Berthelot mixing rules
                combined_epsilon = math.sqrt(epsilon1 * epsilon2)
                combined_sigma = (sigma1 + sigma2) / 2
            elif mixing_method == 'geom':
                # Geometric mixing rules
                combined_epsilon = math.sqrt(epsilon1 * epsilon2)
                combined_sigma = math.sqrt(sigma1 * sigma2)
            else:
                raise ValueError("Invalid mixing method. Use 'LB' for Lorentz-Berthelot or 'geom' for geometric.")
                
            combined_line = f"pair_coeff {atom_type1} {atom_type2} lj/cut/coul/long {combined_epsilon:.3f} {combined_sigma:.3f} # {atom_type1}-{atom_type2}"
            combined_lines.append(combined_line)

    return combined_lines
    

def combine_parameters(lj_params, bond_params, angle_params, dihedral_params, improper_params, combined_lj_lines, ni_o_lines, mixing_method):
    """
    Combine all parameter types into a single string.
    
    Args:
        lj_params (str): Lennard-Jones parameters.
        bond_params (str): Bond parameters.
        angle_params (str): Angle parameters.
        dihedral_params (str): Dihedral parameters.
        improper_params (str): Improper parameters.
        combined_lj_lines (list): Combined Lennard-Jones parameters.
        ni_o_lines (list): Ni-O Lennard-Jones parameters.
        mixing_method (str): The mixing method used ('LB' or 'geom').
        
    Returns:
        str: Combined parameters.
    """
    # Debugging: Print the mixing method
    print(f"Debug: Using mixing method: {mixing_method}")

    if mixing_method == 'LB':
        mixing_method_label = "Combined Lorentz-Berthelot Parameters"
    elif mixing_method == 'geom':
        mixing_method_label = "Geometric Mixing Rules"
    else:
        raise ValueError("Invalid mixing method. Use 'LB' for Lorentz-Berthelot or 'geom' for geometric.")

    return '\n\n'.join([
        "# Ligand OPLS-AA Pair Coeffs\n" + lj_params,
        "# Ligand OPLS-AA Bond Coeffs\n" + bond_params,
        "# Ligand OPLS-AA Angle Coeffs\n" + angle_params,
        "# Ligand OPLS-AA Dihedral Coeffs\n" + dihedral_params,
        "# Ligand OPLS-AA Improper Coeffs\n" + improper_params,
        "\n# Ni-O Lennard-Jones Parameters\n" + '\n'.join(ni_o_lines),
        f"\n# {mixing_method_label} between Ni, O, and Ligand Atoms\n" + '\n'.join(combined_lj_lines),
    ])
    

def save_combined_parameters(filepath, combined_data):
    """
    Save combined parameters to a file.
    
    Args:
        filepath (str): Path to the output file.
        combined_data (str): Combined parameters.
    """
    with open(filepath, 'w') as file:
        file.write(combined_data)

###
##
#ANALYSIS TOOLS
##
###

# Function to parse LJ parameters 
def parse_lj_parameters(lines):
    lj_params = {}
    for line in lines:
        if line.startswith('pair_coeff'):
            parts = line.split()
            if 'lj/cut' in parts or 'lj/cut/coul/long' in parts:
                atom_type = parts[2]
                epsilon = float(parts[4])
                sigma = float(parts[5])
                lj_params[atom_type] = {'epsilon': epsilon, 'sigma': sigma}
    return lj_params

# Function to parse Buckingham parameters
def parse_buckingham_parameters(lines):
    buckingham_params = {}
    for line in lines:
        if line.startswith('pair_coeff'):
            parts = line.split()
            if 'buck/coul/long' in parts:
                atom_type1 = parts[1]
                atom_type2 = parts[2]
                A = float(parts[4]) 
                rho = float(parts[5]) 
                C = float(parts[6]) 
                buckingham_params[(atom_type1, atom_type2)] = {'A': A, 'rho': rho, 'C': C}
    return buckingham_params

# Function to verify OPLS-AA parameters
def verify_opls_aa_params(lj_params):
    for atom_type, params in lj_params.items():
        epsilon = params['epsilon']
        sigma = params['sigma']
        if not (0 < epsilon < 10) or not (3 < sigma < 5):
            print(f"Unusual OPLS-AA parameters for atom type {atom_type}: epsilon = {epsilon}, sigma = {sigma}")
        else:
            print(f"OPLS-AA parameters for atom type {atom_type} are within expected range.")

# Function to verify Buckingham parameters
def verify_buckingham_params(buckingham_params):
    for pair, params in buckingham_params.items():
        A = params['A']
        rho = params['rho']
        C = params['C']
        print(f"Buckingham parameters for pair {pair}: A = {A:.2e} eV, rho = {rho:.2e} Å, C = {C:.2e} eV*Å^6")

# Function to verify combined parameters
def verify_combined_params(lj_params, lines):
    mismatches = []
    for line in lines:
        if line.startswith('pair_coeff'):
            parts = line.split()
            if 'lj/cut/coul/long' in parts:
                atom_type1 = parts[1]
                atom_type2 = parts[2]
                epsilon_comb = float(parts[4])
                sigma_comb = float(parts[5])
                if atom_type1 in lj_params and atom_type2 in lj_params:
                    epsilon1 = lj_params[atom_type1]['epsilon']
                    sigma1 = lj_params[atom_type1]['sigma']
                    epsilon2 = lj_params[atom_type2]['epsilon']
                    sigma2 = lj_params[atom_type2]['sigma']
                    expected_epsilon = math.sqrt(epsilon1 * epsilon2)
                    expected_sigma = (sigma1 + sigma2) / 2
                    if not (math.isclose(epsilon_comb, expected_epsilon, rel_tol=1e-2) and math.isclose(sigma_comb, expected_sigma, rel_tol=1e-2)):
                        mismatches.append((atom_type1, atom_type2, expected_epsilon, epsilon_comb, expected_sigma, sigma_comb))
    for mismatch in mismatches:
        atom_type1, atom_type2, expected_epsilon, epsilon_comb, expected_sigma, sigma_comb = mismatch
        print(f"Mismatch in combined parameters for {atom_type1}-{atom_type2}: expected epsilon = {expected_epsilon:.2e}, sigma = {expected_sigma:.2e}, got epsilon = {epsilon_comb:.2e}, sigma = {sigma_comb:.2e}")

###
##
#CONVERSION TOOLS
##
###
def lammps_data_convert_to_real_numbers(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('Atoms # full'):
                file.write(line)
                break
            file.write(line)
        
        # Write the rest of the file with converted numbers
        for line in lines[lines.index('Atoms # full\n') + 1:]:
            parts = line.split()
            if len(parts) >= 6:
                parts[4] = f"{float(parts[4]):.6f}"
                parts[5] = f"{float(parts[5]):.6f}"
                parts[6] = f"{float(parts[6]):.6f}"
                file.write(' '.join(parts) + '\n')
            else:
                file.write(line)
