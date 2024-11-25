import os

def run_packmol_with_system(input_script_path: str):
    """
    Runs the Packmol input script using os.system.
    """
    command = f"packmol < {input_script_path}"
    exit_code = os.system(command)
    if exit_code == 0:
        print("[INFO] Packmol executed successfully.")
    else:
        print(f"[ERROR] Packmol execution failed with exit code {exit_code}.")

run_packmol_with_system("output/packmol_files/packmol_input.inp")