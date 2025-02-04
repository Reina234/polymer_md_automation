import subprocess
import numpy as np
import csv
import os


def run_simulation(solvent_smiles, compressibility, monomer_smiles_list, N, T):
    """
    Dummy simulation runner. Replace this with your actual simulation call.
    Returns an object with paths to the simulation output files.
    """

    class GromacsOutputs:
        def __init__(self):
            self.tpr = "simulation.tpr"  # Input/control file
            self.xtc = "trajectory.xtc"  # Trajectory file
            self.gro = "final.gro"  # Final configuration file

    return GromacsOutputs()


def extract_radius_of_gyration(outputs):
    # Run gmx gyrate to get Rg vs. time.
    cmd = ["gmx", "gyrate", "-f", outputs.xtc, "-s", outputs.tpr, "-o", "gyrate.xvg"]
    subprocess.run(cmd, check=True)
    data = np.loadtxt("gyrate.xvg", comments=("#", "@"))
    Rg_mean = np.mean(data[:, 1])
    Rg_std = np.std(data[:, 1])
    return Rg_mean, Rg_std


def extract_diffusion_coefficient(outputs):
    # Use gmx msd to compute mean-squared displacement.
    cmd = ["gmx", "msd", "-f", outputs.xtc, "-s", outputs.tpr, "-o", "msd.xvg"]
    subprocess.run(cmd, check=True)
    data = np.loadtxt("msd.xvg", comments=("#", "@"))
    t = data[:, 0]
    msd = data[:, 1]
    # Fit the last 20% of the data (the diffusive regime)
    start = int(0.8 * len(t))
    slope, _ = np.polyfit(t[start:], msd[start:], 1)
    # For 3D diffusion, D = slope/6.
    D = slope / 6.0
    return D


def extract_sasa(outputs):
    # Run gmx sasa.
    cmd = ["gmx", "sasa", "-f", outputs.xtc, "-s", outputs.tpr, "-o", "sasa.xvg"]
    subprocess.run(cmd, check=True)
    data = np.loadtxt("sasa.xvg", comments=("#", "@"))
    SASA_mean = np.mean(data[:, 1])
    SASA_std = np.std(data[:, 1])
    return SASA_mean, SASA_std


def extract_end_to_end_distance(outputs):
    # Use gmx distance to compute the end-to-end distance.
    cmd = [
        "gmx",
        "distance",
        "-f",
        outputs.xtc,
        "-s",
        outputs.tpr,
        "-o",
        "end_to_end.xvg",
    ]
    subprocess.run(cmd, check=True)
    data = np.loadtxt("end_to_end.xvg", comments=("#", "@"))
    E2E_mean = np.mean(data[:, 1])
    E2E_std = np.std(data[:, 1])
    return E2E_mean, E2E_std


def extract_shape_descriptors(outputs):
    """
    Compute shape descriptors (asphericity and acylindricity) from the eigenvalues
    of the gyration tensor. We assume that 'gmx gyrate -dg' outputs a file 'eigen.xvg'
    with columns: time, eigen1, eigen2, eigen3.
    """
    cmd = [
        "gmx",
        "gyrate",
        "-f",
        outputs.xtc,
        "-s",
        outputs.tpr,
        "-o",
        "eigen.xvg",
        "-dg",
    ]
    subprocess.run(cmd, check=True)
    data = np.loadtxt("eigen.xvg", comments=("#", "@"))
    asphericities = []
    acylindricities = []
    for row in data:
        # Sort eigenvalues in descending order
        eig = np.sort(row[1:4])[::-1]
        asph = eig[0] - 0.5 * (eig[1] + eig[2])
        acyl = eig[1] - eig[2]
        asphericities.append(asph)
        acylindricities.append(acyl)
    asph_mean = np.mean(asphericities)
    asph_std = np.std(asphericities)
    acyl_mean = np.mean(acylindricities)
    acyl_std = np.std(acylindricities)
    return asph_mean, asph_std, acyl_mean, acyl_std


def extract_persistence_length(outputs):
    """
    Placeholder for persistence length extraction. This might involve computing
    the bond vector correlation along the polymer chain and fitting an exponential.
    Here we assume a command that outputs 'persistence.xvg' exists.
    """
    try:
        cmd = [
            "gmx",
            "persistence",
            "-f",
            outputs.xtc,
            "-s",
            outputs.tpr,
            "-o",
            "persistence.xvg",
        ]
        subprocess.run(cmd, check=True)
        data = np.loadtxt("persistence.xvg", comments=("#", "@"))
        distances = data[:, 0]
        correlation = data[:, 1]
        mask = correlation > 0  # Avoid log(0)
        distances = distances[mask]
        log_corr = np.log(correlation[mask])
        slope, _ = np.polyfit(distances, log_corr, 1)
        lp = -1.0 / slope if slope != 0 else np.nan
    except Exception:
        lp = np.nan  # If analysis is not available, return NaN.
    return lp


def extract_pair_correlation(outputs):
    """
    Compute the radial distribution function (RDF) and extract the location
    (r_peak) and height (g_peak) of the first peak.
    """
    cmd = ["gmx", "rdf", "-f", outputs.xtc, "-s", outputs.tpr, "-o", "rdf.xvg"]
    subprocess.run(cmd, check=True)
    data = np.loadtxt("rdf.xvg", comments=("#", "@"))
    r = data[:, 0]
    g_r = data[:, 1]
    peak_index = np.argmax(g_r)
    r_peak = r[peak_index]
    g_peak = g_r[peak_index]
    return r_peak, g_peak


def extract_relaxation_time(outputs):
    """
    Placeholder for viscoelastic analysis. Assuming we can extract a time series
    (e.g. from pressure or stress autocorrelation) that decays exponentially:
       autocorr ~ A * exp(-t/tau)
    We fit log(autocorr) vs. time to extract the relaxation time tau.
    """
    try:
        cmd = [
            "gmx",
            "energy",
            "-f",
            outputs.xtc,
            "-s",
            outputs.tpr,
            "-o",
            "viscoelastic.xvg",
        ]
        subprocess.run(cmd, check=True)
        data = np.loadtxt("viscoelastic.xvg", comments=("#", "@"))
        t = data[:, 0]
        autocorr = data[:, 1]
        mask = autocorr > 0
        t_fit = t[mask]
        log_autocorr = np.log(autocorr[mask])
        slope, _ = np.polyfit(t_fit, log_autocorr, 1)
        tau = -1.0 / slope if slope != 0 else np.nan
    except Exception:
        tau = np.nan
    return tau


def process_simulation(solvent_smiles, compressibility, monomer, N, T, csv_filename):
    # Run the simulation and obtain Gromacs output paths.
    outputs = run_simulation(solvent_smiles, compressibility, [monomer], N, T)

    # Extract mesoscale parameters.
    Rg_mean, Rg_std = extract_radius_of_gyration(outputs)
    D = extract_diffusion_coefficient(outputs)
    SASA_mean, SASA_std = extract_sasa(outputs)
    E2E_mean, E2E_std = extract_end_to_end_distance(outputs)
    asph_mean, asph_std, acyl_mean, acyl_std = extract_shape_descriptors(outputs)
    lp = extract_persistence_length(outputs)
    r_peak, g_peak = extract_pair_correlation(outputs)
    tau = extract_relaxation_time(outputs)

    row = {
        "monomer": monomer,
        "N": N,
        "T": T,
        "solvent_smiles": solvent_smiles,
        "compressibility": compressibility,
        "Rg_mean": Rg_mean,
        "Rg_std": Rg_std,
        "Diffusion_Coefficient": D,
        "SASA_mean": SASA_mean,
        "SASA_std": SASA_std,
        "End_to_End_mean": E2E_mean,
        "End_to_End_std": E2E_std,
        "Asphericity_mean": asph_mean,
        "Asphericity_std": asph_std,
        "Acylindricity_mean": acyl_mean,
        "Acylindricity_std": acyl_std,
        "Persistence_Length": lp,
        "RDF_r_peak": r_peak,
        "RDF_g_peak": g_peak,
        "Relaxation_Time": tau,
    }

    # Append the new row to the CSV file (create header if file is new).
    file_exists = os.path.isfile(csv_filename)
    with open(csv_filename, mode="a", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

    print("Simulation processed and results appended to", csv_filename)


if __name__ == "__main__":
    # Example input parameters.
    solvent_smiles = "O"  # e.g., water
    compressibility = 4.5e-10  # example compressibility value
    monomer = "C=C"  # example monomer SMILES
    N = 100  # polymer chain length
    T = 300  # temperature in Kelvin
    csv_filename = "sim_results.csv"

    process_simulation(solvent_smiles, compressibility, monomer, N, T, csv_filename)
