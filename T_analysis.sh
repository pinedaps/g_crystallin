#!/usr/bin/env bash

# -----------------------------
#   T_analysis.sh
#   Temperature workflow:
#     • Parse T inputs
#     • Run topology generation
#     • Run duello scan
#     • Run plotting
# -----------------------------

set -euo pipefail

#######################################
# Help message
#######################################

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Temperature input options (choose one):
  --pH <value>
	Provide a pH value
  --tmin <value> --tmax <value> --tstep <value>
        Generate a temperature array from min to max with given step.
  --temps <comma-separated list>
        Provide an explicit list, e.g. --temps 300,310,320
  --pdb <path>        
        Pdb path

Other options:
  --outdir <path>     Output directory (default: ./T_analysis_out)  
  -h, --help          Show this help message

Example:
  $0 --pH 7.1 --tmin 280 --tmax 320 --tstep 10 --pdb ../raw_data/1AMM
  $0 --pH 7.1 --temps 290,300,310 --pdb ../raw_data/1AMM
EOF
}

#######################################
# Default values
#######################################

OUTDIR="./T_analysis_out"
T_ARRAY=()

#######################################
# Parse arguments
#######################################

# To support both long options and getopts,
# we manually parse long options:
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tmin)
            TMIN="$2"
            shift 2
            ;;
        --tmax)
            TMAX="$2"
            shift 2
            ;;
        --tstep)
            TSTEP="$2"
            shift 2
            ;;
        --temps)
            USER_TEMPS="$2"
            shift 2
            ;;
        --pH)
            PH="$2"
            shift 2
            ;;
        --pdb)
            PDB="$2"
            shift 2
            ;;
	--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

#######################################
# Build temperature array
#######################################

if [[ -n "${USER_TEMPS:-}" ]]; then
    IFS=',' read -ra T_ARRAY <<< "$USER_TEMPS"
elif [[ -n "${TMIN:-}" && -n "${TMAX:-}" && -n "${TSTEP:-}" ]]; then
    T_ARRAY=()
    T="$TMIN"
    while (( $(echo "$T <= $TMAX" | bc -l) )); do
        T_ARRAY+=("$T")
        T=$(echo "$T + $TSTEP" | bc -l)
    done
else
    echo "ERROR: You must provide either --temps or --tmin/--tmax/--tstep." >&2
    exit 1
fi

#########################################
# Create output directory and filenames
#########################################

mkdir -p "$OUTDIR"
TOPO_DIR="$OUTDIR/topologies"
SCAN_DIR="$OUTDIR/scans"
PLOT_DIR="$OUTDIR/plots"

mkdir -p "$TOPO_DIR" "$SCAN_DIR" "$PLOT_DIR"

FILE="${PDB##*/}"
XYZ_OUT="${OUTDIR}/${FILE}.xyz"

echo "pH: $PH"
echo "Temperatures: ${T_ARRAY[*]}"
echo "Output directory: $OUTDIR"
echo "$XYZ_OUT"
echo

#######################################
# Step 1: Generate topology files
#######################################

echo "=== Generating topology files ==="

for T in "${T_ARRAY[@]}"; do
    TOPO_OUT="${TOPO_DIR}/topology_${FILE}_T${T}.yaml"
    echo "  Running topology for pdb = $FILE at T = $T → $TOPO_OUT"

    python3 pdb2xyz/__init__2.py \
        -i "$PDB" \
	-o "$XYZ_OUT" \
	-t "$TOPO_OUT" \
	--pH "$PH" \
        --T "$T"
done

echo "Topology generation complete."
echo

#######################################
# Step 2: Run duello scan
#######################################

echo "=== Running duello scan ==="

for T in "${T_ARRAY[@]}"; do
    TOPO_IN="${TOPO_DIR}/topology_${FILE}_T${T}.yaml"
    SCAN_OUT="${SCAN_DIR}/scan_T${T}.dat"
    echo "  duello scan for T = $T → $SCAN_OUT"
    
    duello scan --mol1 "$XYZ_OUT" \
                --mol2 "$XYZ_OUT" \
                --rmin 20 \
                --rmax 80 \
                --dr 0.5 \
		--resolution 0.7 \
		--cutoff 80  \
                --top "$TOPO_IN"  \
		--molarity 0.115  \
		--temperature "$T" \
		--pmf "$SCAN_OUT"
done

echo "Duello scans complete."
echo

exit

#######################################
# Step 3: Plot results
#######################################
echo "=== Plotting results ==="

# Replace with your actual plotting script:
python3 plot_duello_results.py \
    --input_dir "$SCAN_DIR" \
    --output_dir "$PLOT_DIR"

echo "Plots generated in: $PLOT_DIR"
echo
echo "=== Done! ==="

