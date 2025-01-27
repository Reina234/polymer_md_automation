def replace_itp_names(itp_path, output_path, new_name="SOL"):
    """
    Replaces all residue/molecule names in a GROMACS .itp file with a specified new name.

    Args:
        itp_path (str): Path to the input .itp file.
        output_path (str): Path to save the modified .itp file.
        new_name (str, optional): The new name to replace existing residue/molecule names. Default is "SOL".

    Returns:
        None: Writes the modified content to the specified output file.
    """
    import re

    with open(itp_path, "r") as file:
        lines = file.readlines()

    updated_lines = []
    in_moleculetype = False
    in_atoms_section = False

    for line in lines:
        stripped_line = line.strip()

        # Detect [ moleculetype ] section
        if stripped_line.startswith("[ moleculetype ]"):
            in_moleculetype = True
            updated_lines.append(line)
            continue

        # Detect [ atoms ] section
        if stripped_line.startswith("[ atoms ]"):
            in_atoms_section = True
            updated_lines.append(line)
            continue

        # Replace name in [ moleculetype ] section (2nd column after header)
        if in_moleculetype and stripped_line and not stripped_line.startswith(";"):
            parts = line.split()
            if len(parts) > 1:
                parts[0] = new_name  # Replace the molecule name
                line = " ".join(parts) + "\n"
            in_moleculetype = False  # Only applies to the next line after header

        # Replace residue name in [ atoms ] section (5th column)
        if in_atoms_section and stripped_line and not stripped_line.startswith(";"):
            parts = re.split(r"(\s+)", line)  # Preserve whitespace
            if len(parts) > 9:
                parts[8] = new_name  # Replace residue name (column 5)
                line = "".join(parts)

        updated_lines.append(line)

    # Write modified content to the output file
    with open(output_path, "w") as file:
        file.writelines(updated_lines)

    print(f"✅ Successfully replaced names with '{new_name}' in {output_path}")


replace_itp_names("temp/solvent.itp", "itp_manual.itp")
from modules.utils.atomistic.file_utils import rename_specific_residue_name_from_gro

rename_specific_residue_name_from_gro(
    "temp/polymer_in_solvent.gro", "TMZK", "SOL", output_name="manual_gro"
)
