import xml.etree.ElementTree as ET
from typing import List, Dict
import xml.etree.ElementTree as ET
from typing import List, Tuple, Dict, Optional
import xml.dom.minidom


class VOTCAMappingGenerator:
    """
    Generates a VOTCA-compatible XML mapping file for coarse-grained simulations.
    """

    def __init__(
        self,
        molecule_name: str,
        bead_mappings: List[Dict[str, any]],
        bonds: Optional[List[Tuple[str, str]]] = None,
        angles: Optional[List[Tuple[str, str, str]]] = None,
        verify_itp: bool = False,
        itp_data: Optional[List[str]] = None,
    ):
        """
        Initializes the VOTCA XML generator with mapping and bonded interaction data.

        :param molecule_name: Name of the molecule in CG representation.
        :param bead_mappings: List of dicts, each containing:
                              - "unique_name": Bead identifier.
                              - "bead_type": Bead type.
                              - "atom_indices": List of atom indices in this bead.
        :param bonds: Optional list of bead-pair tuples defining bonded interactions.
        :param angles: Optional list of bead-triples defining angle interactions.
        :param verify_itp: If True, validate bond consistency against an .itp file.
        :param itp_data: Optional .itp file data as a list of strings.
        """
        self.molecule_name = molecule_name
        self.bead_mappings = bead_mappings
        self.bonds = bonds if bonds else []
        self.angles = angles if angles else []

        if verify_itp and itp_data:
            self._verify_bonds_against_itp(itp_data)

    def _verify_bonds_against_itp(self, itp_data: List[str]):
        """
        Validates computed CG bead bonds against those found in an .itp file.
        """
        bonds_section = False
        itp_bonds = []

        for line in itp_data:
            if "[ bonds ]" in line:
                bonds_section = True
                continue
            if bonds_section and line.strip() == "":
                break  # End of section
            if bonds_section:
                parts = line.split()
                if len(parts) >= 2:
                    atom1, atom2 = int(parts[0]), int(parts[1])
                    bead1 = next(
                        (
                            b["unique_name"]
                            for b in self.bead_mappings
                            if atom1 in b["atom_indices"]
                        ),
                        None,
                    )
                    bead2 = next(
                        (
                            b["unique_name"]
                            for b in self.bead_mappings
                            if atom2 in b["atom_indices"]
                        ),
                        None,
                    )
                    if bead1 and bead2:
                        itp_bonds.append((bead1, bead2))

        # Compare with existing bonds
        missing_bonds = [b for b in itp_bonds if b not in self.bonds]
        if missing_bonds:
            print(f"[WARNING] Missing bonds in VOTCA mapping: {missing_bonds}")

    def save_to_xml(self, filename: str):
        """
        Saves the mapping data to a VOTCA-compatible XML file.
        """
        root = ET.Element("cg_molecule")
        ET.SubElement(root, "name").text = self.molecule_name
        ET.SubElement(root, "ident").text = self.molecule_name

        topology_elem = ET.SubElement(root, "topology")
        cg_beads_elem = ET.SubElement(topology_elem, "cg_beads")

        # Add bead mappings
        for bead in self.bead_mappings:
            bead_elem = ET.SubElement(cg_beads_elem, "cg_bead")
            ET.SubElement(bead_elem, "name").text = bead["unique_name"]
            ET.SubElement(bead_elem, "type").text = bead["bead_type"]
            ET.SubElement(bead_elem, "mapping").text = (
                f"M{self.bead_mappings.index(bead) + 1}"
            )
            ET.SubElement(bead_elem, "beads").text = " ".join(
                map(str, bead["atom_indices"])
            )

        # Add bonded interactions
        if self.bonds or self.angles:
            cg_bonded_elem = ET.SubElement(topology_elem, "cg_bonded")

            if self.bonds:
                bond_elem = ET.SubElement(cg_bonded_elem, "bond")
                ET.SubElement(bond_elem, "name").text = "bond"
                ET.SubElement(bond_elem, "beads").text = "\n".join(
                    f"{b1} {b2}" for b1, b2 in self.bonds
                )

            if self.angles:
                angle_elem = ET.SubElement(cg_bonded_elem, "angle")
                ET.SubElement(angle_elem, "name").text = "angle"
                ET.SubElement(angle_elem, "beads").text = "\n".join(
                    f"{b1} {b2} {b3}" for b1, b2, b3 in self.angles
                )

        # Add maps
        maps_elem = ET.SubElement(root, "maps")
        for i, bead in enumerate(self.bead_mappings):
            map_elem = ET.SubElement(maps_elem, "map")
            ET.SubElement(map_elem, "name").text = f"M{i + 1}"
            ET.SubElement(map_elem, "weights").text = " ".join(
                ["1.0"] * len(bead["atom_indices"])
            )

        xml_str = ET.tostring(root, encoding="utf-8").decode("utf-8")
        formatted_xml = self._prettify_xml(xml_str)

        with open(filename, "w") as f:
            f.write(formatted_xml)

        print(f"[INFO] VOTCA XML mapping saved as {filename}")

    def _prettify_xml(self, xml_str: str) -> str:
        """Formats XML with proper indentation."""
        dom = xml.dom.minidom.parseString(xml_str)
        return dom.toprettyxml(indent="  ")  # Two-space indentation
