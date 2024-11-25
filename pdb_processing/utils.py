from typing import Optional
import os 

def validate_pdb_lightweight(pdb_path: str) -> None:
    """
    Validates a PDB file to ensure compatibility with Packmol.
    Checks for required fields and formatting.
    """
    required_columns = [6, 7, 8]  # Columns for x, y, z coordinates
    try:
        with open(pdb_path, "r") as file:
            for line in file:
                if line.startswith(("ATOM", "HETATM")):
                    fields = line.split()
                    if len(fields) < max(required_columns):
                        raise ValueError(f"Malformed ATOM/HETATM line: {line.strip()}")
                    # Check if coordinates are valid floats
                    for col in required_columns:
                        try:
                            float(fields[col])
                        except ValueError:
                            raise ValueError(f"Invalid coordinate in line: {line.strip()}")
                elif line.startswith("TER"):
                    continue
                else:
                    pass  # Ignore other records for Packmol compatibility
        print("[+] PDB validated for Packmol.")
    except Exception as e:
        raise ValueError(f"PDB validation for Packmol failed: {e}")



def handle_existing_file(file_path: str) -> Optional[str]:
    """
    Check if the file exists and prompt the user to decide whether to overwrite it.
    :param file_path: Path to the file to check.
    :return: The file path if it should not be overwritten, otherwise None.
    """
    if os.path.exists(file_path):
        print(f"[!] File '{file_path}' already exists.")
        overwrite = input("Do you want to overwrite the existing file? [y/n]: ").strip().lower()
        if overwrite == "y":
            print("[+] Overwriting existing file. Proceeding with file generation.")
            return file_path
        elif overwrite != "n":
            raise ValueError("Invalid input. Please enter 'y' for yes or 'n' for no.")
        print("[!] File generation aborted. Returning existing file path.")
    return None
