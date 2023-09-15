# Utility functions for genome data analysis.

def parse_region(region):
    # Parse a region string into its chromosome and start and end positions.
    region_parts = region.split(":")
    chromosome = region_parts[0]
    start_position = int(region_parts[1].split("-")[0])
    end_position = int(region_parts[1].split("-")[1])

    return chromosome, start_position, end_position
