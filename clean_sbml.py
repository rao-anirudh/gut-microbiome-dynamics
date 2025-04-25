from libsbml import readSBML, writeSBMLToFile


def truncate_id(original_id, max_len=250):
    """
    Truncates an identifier to ensure it does not exceed the specified length.

    Parameters:
    original_id (str): The original identifier to be truncated.
    max_len (int): The maximum allowed length for the identifier. Default is 250 characters.

    Returns:
    str: A truncated identifier, with a hash added to ensure uniqueness.
    """
    # Truncate the ID to half of the max_len, then append a hash-based suffix
    return original_id[:max_len // 2] + "_" + str(abs(hash(original_id)) % 10 ** 6)


def fix_long_ids(sbml_path, output_path, max_len=256):
    """
    Reads an SBML file, checks the length of all IDs, and truncates those
    that exceed the specified maximum length. Saves the cleaned model to a new file.

    Parameters:
    sbml_path (str): Path to the original SBML file to be processed.
    output_path (str): Path where the cleaned SBML file will be saved.
    max_len (int): The maximum allowed length for any identifier in the SBML model. Default is 256 characters.

    This function processes the following collections in the SBML model:
    - Reactions
    - Species
    - Parameters
    - Compartments

    It updates the IDs and optionally the names of any elements whose ID exceeds the maximum length.
    """
    # Read the SBML model from the input file
    doc = readSBML(sbml_path)
    model = doc.getModel()

    # Process the key collections: reactions, species, parameters, and compartments
    for collection in [
        model.getListOfReactions(),
        model.getListOfSpecies(),
        model.getListOfParameters(),
        model.getListOfCompartments(),
    ]:
        # Iterate over each object in the collection
        for obj in collection:
            # Check if the object's ID exceeds the maximum length
            if len(obj.getId()) > max_len:
                # Truncate the ID and create a new, unique identifier
                new_id = truncate_id(obj.getId(), max_len)

                # Print the truncation action for clarity
                print(f"Truncating ID:\n   Old: {obj.getId()}\n   New: {new_id}")

                # Set the new truncated ID and optionally update the name
                obj.setId(new_id)
                obj.setName(new_id)  # Optional: update the name too

    # Write the cleaned SBML model to the output file
    writeSBMLToFile(doc, output_path)
    print(f"\nClean model saved to: {output_path}")


fix_long_ids("MODEL1310110043_url_large.xml", "MODEL1310110043_url_large_cleaned.xml")
