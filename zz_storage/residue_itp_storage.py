import sqlite3
import json


class ResidueITPStorage:
    """Handles `.itp` storage & retrieval per **residue SMILES**."""

    def __init__(self, db_path="residue_itp.db"):
        self.db_path = db_path
        self._initialize_db()

    def _initialize_db(self):
        """Creates SQLite storage if it doesn't exist."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS itp_storage (
                    residue_smiles TEXT PRIMARY KEY,
                    atoms TEXT,
                    bonds TEXT,
                    angles TEXT,
                    dihedrals TEXT,
                    constraints TEXT
                )
                """
            )
            conn.commit()

    def store_residue_itp(self, segmented_itp):
        """Stores `.itp` data per residue SMILES."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            for residue_smiles, sections in segmented_itp.items():
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO itp_storage (residue_smiles, atoms, bonds, angles, dihedrals, constraints)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (
                        residue_smiles,
                        json.dumps(sections.get("atoms", [])),
                        json.dumps(sections.get("bonds", [])),
                        json.dumps(sections.get("angles", [])),
                        json.dumps(sections.get("dihedrals", [])),
                        json.dumps(sections.get("constraints", [])),
                    ),
                )
            conn.commit()

    def load_residue_itp(self, residue_smiles):
        """Retrieves `.itp` data for a given residue SMILES."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(
                "SELECT atoms, bonds, angles, dihedrals, constraints FROM itp_storage WHERE residue_smiles = ?",
                (residue_smiles,),
            )
            row = cursor.fetchone()
            if row:
                return {
                    "atoms": json.loads(row[0]) if row[0] else [],
                    "bonds": json.loads(row[1]) if row[1] else [],
                    "angles": json.loads(row[2]) if row[2] else [],
                    "dihedrals": json.loads(row[3]) if row[3] else [],
                    "constraints": json.loads(row[4]) if row[4] else [],
                }
            return None


itp_storage = ResidueITPStorage()
