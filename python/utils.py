"""Utility functions for genome data analysis."""

def parse_region(region):
    """Parse a region string into its chromosome and start and end positions."""
    region_parts = region.split(":")
    chromosome = str(region_parts[0])

    try:
        start_position = int(region_parts[1].split("-")[0])
        end_position = int(region_parts[1].split("-")[1])
    except IndexError:
        start_position, end_position = None, None

    return chromosome, start_position, end_position

def get_info_field_column(vcf_data):
    """Return the column index of the INFO field in a VCF file."""
    index = vcf_data.apply(lambda col: col.astype(str).str.contains("SVTYPE=").any(), axis=0).idxmax()
    return index

def get_info_field_value(info_field, field_name):
    """
    Get the value of a field in the INFO field of a VCF file.

    Args:
        info_field (str): The INFO field.
        field_name (str): The name of the field to get the value of.

    Returns:
        str: The value of the field.
    """

    # Split the INFO field into its parts.
    info_field_parts = info_field.split(";")

    # Get the field value.
    field_value = ""
    for info_field_part in info_field_parts:
        if info_field_part.startswith("{}=".format(field_name)):
            field_value = info_field_part.split("=")[1]
            break

    # Return the field value.
    return field_value
