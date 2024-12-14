class AtomtypesManager:
    """
    Manages solvent atomtypes with persistent storage in a JSON file.
    """

    STORAGE_FILE = SOLVENT_JSON_PATH

    def __init__(self):
        # Ensure the JSON file exists
        if not os.path.exists(self.STORAGE_FILE):
            with open(self.STORAGE_FILE, "w") as file:
                json.dump({}, file)
        logger.info(f"[+] Using storage file: {self.STORAGE_FILE}")

    def extract_and_store_atomtypes(self, solvent_file: str, solvent_name: str) -> None:
        """
        Extracts the [ atomtypes ] section from a solvent file and stores it in JSON.

        Args:
            solvent_file (str): Path to the solvent .itp file.
            solvent_name (str): Name of the solvent to use as a key.
        """
        content = self._read_file(solvent_file)
        atomtypes_section = self._extract_section(content, "atomtypes")
        if not atomtypes_section:
            raise ValueError(f"No [ atomtypes ] section found in {solvent_file}.")

        # Load existing data from JSON
        storage_data = self._load_storage()

        # Store the atomtypes section
        storage_data[solvent_name] = atomtypes_section
        self._save_storage(storage_data)
        logger.info(f"[+] Stored atomtypes for {solvent_name}.")

    def retrieve_atomtypes(self, solvent_name: str) -> List[str]:
        """
        Retrieves the [ atomtypes ] section for a specified solvent.

        Args:
            solvent_name (str): Name of the solvent.

        Returns:
            List[str]: Lines of the [ atomtypes ] section.
        """
        storage_data = self._load_storage()
        if solvent_name not in storage_data:
            raise ValueError(f"No atomtypes data found for solvent: {solvent_name}")
        logger.info(f"[+] Retrieved atomtypes for {solvent_name}.")
        return storage_data[solvent_name]

    def _read_file(self, file_path: str) -> List[str]:
        """Helper to read a file."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        with open(file_path, "r") as file:
            return file.readlines()

    def _extract_section(self, content: List[str], section_name: str) -> List[str]:
        """Helper to extract a specific section."""
        in_section = False
        section_lines = []
        for line in content:
            if line.strip().startswith(f"[ {section_name} ]"):
                in_section = True
                section_lines.append(line)
                continue
            if in_section:
                if line.strip().startswith("[") and not line.strip().startswith(
                    f"[ {section_name} ]"
                ):
                    break
                section_lines.append(line)
        return section_lines

    def _load_storage(self) -> Dict[str, List[str]]:
        """Helper to load data from the JSON storage."""
        with open(self.STORAGE_FILE, "r") as file:
            return json.load(file)

    def _save_storage(self, data: Dict[str, List[str]]) -> None:
        """Helper to save data to the JSON storage."""
        with open(self.STORAGE_FILE, "w") as file:
            json.dump(data, file, indent=4)
        logger.info(f"[+] Updated storage file: {self.STORAGE_FILE}")
