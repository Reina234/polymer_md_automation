import xml.etree.ElementTree as ET
import xml.dom.minidom
import re


class VOTCASingleBeadMapping:
    """
    Generates a VOTCA XML mapping for solvents, mapping the entire molecule to a single bead.
    """

    def __init__(self, molecule_name, itp_file):
        """
        :param molecule_name: The name of the solvent molecule in the topology (.tpr).
        :param itp_file: Path to the .itp file for extracting atomic masses.
        :param output_xml: Path to save the generated VOTCA XML.
        """
        self.molecule_name = molecule_name
        self.itp_file = itp_file
        self.atom_masses = self._parse_itp_mass()

    def _parse_itp_mass(self):
        """
        Extracts atomic masses from the .itp file.
        :return: List of atomic masses.
        """
        masses = []
        reading_atoms = False

        with open(self.itp_file, "r") as file:
            for line in file:
                if "[ atoms ]" in line:
                    reading_atoms = True
                    continue
                if reading_atoms and line.strip() == "":
                    break  # End of section

                if reading_atoms:
                    parts = re.split(r"\s+", line.strip())
                    if len(parts) > 7:
                        try:
                            masses.append(float(parts[7]))  # Mass is in the 8th column
                        except ValueError:
                            pass  # Ignore lines without valid masses

        return masses

    def generate_xml(self, output_file: str):
        """
        Generates and saves the VOTCA XML mapping file.
        """
        root = ET.Element("cg_molecule")

        # Add molecule name and identifier
        ET.SubElement(root, "name").text = self.molecule_name
        ET.SubElement(root, "ident").text = self.molecule_name

        topology_elem = ET.SubElement(root, "topology")
        cg_beads_elem = ET.SubElement(topology_elem, "cg_beads")

        # Create a single bead mapping
        bead_elem = ET.SubElement(cg_beads_elem, "cg_bead")
        ET.SubElement(bead_elem, "name").text = "S1"
        ET.SubElement(bead_elem, "type").text = "SOL"
        ET.SubElement(bead_elem, "mapping").text = "Solvent"
        ET.SubElement(bead_elem, "beads").text = f"1:{self.molecule_name}:*"

        weights = " ".join([f"{round(m)}" for m in self.atom_masses])

        maps_elem = ET.SubElement(root, "maps")
        map_elem = ET.SubElement(maps_elem, "map")
        ET.SubElement(map_elem, "name").text = "Solvent"
        ET.SubElement(map_elem, "weights").text = weights

        # Format XML output
        xml_str = ET.tostring(root, encoding="utf-8").decode("utf-8")
        formatted_xml = self._prettify_xml(xml_str)

        # Save XML
        with open(output_file, "w") as f:
            f.write(formatted_xml)

        print(f"[INFO] VOTCA XML mapping saved as {output_file}")

    def _prettify_xml(self, xml_str):
        """
        Formats XML with proper indentation.
        """
        dom = xml.dom.minidom.parseString(xml_str)
        return dom.toprettyxml(indent="  ")
