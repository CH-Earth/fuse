#!/bin/bash
#
# clean_camels_spat.sh
#
# Remove the unnecessary time dimension from the latitude, longitude,
# and hruId variables in one or more CAMELS-SPAT NetCDF forcing files.
#
# The shell expands any wildcard before passing the filenames to this script.
# The final command-line argument is interpreted as the output directory.
#
# Usage:
#   ./clean_camels_spat.sh <input_file(s)> <output_directory>
#
# Examples:
#
#   Process one file:
#
#   ./clean_camels_spat.sh \
#       ~/data/CAMELS-SPAT/raw/CAN_05BB001_daymet_distributed.nc \
#       ~/data/CAMELS-SPAT/clean
#
#   Process all matching files:
#
#   ./clean_camels_spat.sh \
#       ~/data/CAMELS-SPAT/raw/CAN_*_daymet_distributed.nc \
#       ~/data/CAMELS-SPAT/clean
#

set -euo pipefail

# Coordinate variables that should not retain the time dimension.
COORD_VARS="latitude,longitude,hruId"

# Check that at least one input file and one output directory were provided.
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input_file(s)> <output_directory>" >&2
    exit 1
fi

# Check that the required NCO commands are available.
for command_name in ncwa ncks; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: Required command '${command_name}' was not found." >&2
        echo "Install NCO before running this script." >&2
        exit 1
    fi
done

# The final argument is the output directory.
OUTPUT_DIR="${!#}"

# Create the output directory if it does not already exist.
mkdir -p "${OUTPUT_DIR}"

# Create a temporary working directory and remove it automatically on exit.
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

# Process every argument except the final output-directory argument.
for ((i = 1; i < $#; i++)); do
    INPUT="${!i}"

    if [[ ! -f "${INPUT}" ]]; then
        echo "ERROR: Input file does not exist: ${INPUT}" >&2
        exit 1
    fi

    BASENAME=$(basename "${INPUT}")
    OUTPUT="${OUTPUT_DIR}/${BASENAME}"

    WORK_FILE="${TMP_DIR}/${BASENAME}"
    COORD_FILE="${TMP_DIR}/coords_${i}.nc"

    echo "Processing: ${INPUT}"
    echo "Writing:    ${OUTPUT}"

    # Work on a temporary copy so the original file is never modified.
    cp "${INPUT}" "${WORK_FILE}"

    # Extract the coordinate variables from the first time step and average
    # over the singleton time selection. This removes the time dimension.
    ncwa -O \
        -a time \
        -v "${COORD_VARS}" \
        -d time,0,0 \
        "${WORK_FILE}" \
        "${COORD_FILE}"

    # Remove the original time-dependent coordinate variables.
    ncks -O \
        -x \
        -v "${COORD_VARS}" \
        "${WORK_FILE}" \
        "${WORK_FILE}"

    # Append the cleaned coordinate variables back into the working file.
    ncks -A \
        -v "${COORD_VARS}" \
        "${COORD_FILE}" \
        "${WORK_FILE}"

    # Move the completed file into the output directory.
    mv "${WORK_FILE}" "${OUTPUT}"

    echo "Completed:  ${OUTPUT}"
    echo
done

echo "All files processed successfully."
