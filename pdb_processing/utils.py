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
