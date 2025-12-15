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
  --T <value>
        Provide a temperature value
  --emin <value> --emax <value> --estep <value>
        Generate a epsilon array from min to max with given step.
  --es <comma-separated list>
        Provide an explicit list, e.g. --es 0.4,0.5,0.6
  --pdb <path>        
        Pdb path

Other options:
  --outdir <path>     Output directory (default: ./<input_pdb>)  
  -h, --help          Show this help message

Example:
  $0 --pH 7.1 --T 293 --emin 0.400 --emax 0.8368 --estep 0.02 --pdb pdbs/XXXX --outdir XXXX
  $0 --pH 7.1 --T 293 --es 0.4,0.5,0.6 --pdb pdbs/XXXX --outdir XXXX
EOF
}

#######################################
# Parse arguments
#######################################

# To support both long options and getopts,
# we manually parse long options:
while [[ $# -gt 0 ]]; do
    case "$1" in
        --emin)
            EMIN="$2"
            shift 2
            ;;
        --emax)
            EMAX="$2"
            shift 2
            ;;
        --estep)
            ESTEP="$2"
            shift 2
            ;;
        --es)
            USER_ES="$2"
            shift 2
            ;;
        --pH)
            PH="$2"
            shift 2
            ;;
	--T)
            T="$2"
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
# Default values
#######################################

WD="$(pwd)"
FILE="${PDB##*/}"
OUTDIR="$FILE"
E_ARRAY=()

#######################################
# Build epsilon array
#######################################

if [[ -n "${USER_ES:-}" ]]; then
    IFS=',' read -ra E_ARRAY <<< "$USER_ES"
elif [[ -n "${EMIN:-}" && -n "${EMAX:-}" && -n "${ESTEP:-}" ]]; then
    E_ARRAY=()
    E="$EMIN"
    while (( $(echo "$E <= $EMAX" | bc -l) )); do
        E_ARRAY+=("$E")
        E=$(echo "$E + $ESTEP" | bc -l)
    done
else
    echo "ERROR: You must provide either --es or --emin/--emax/--estep." >&2
    exit 1
fi

#########################################
# Create output directory and filenames
#########################################

mkdir -p "$OUTDIR"

FILE="${PDB##*/}"
XYZ_OUT="${WD}/${OUTDIR}/${FILE}.xyz"

echo "pH: $PH"
echo "T: $T"
echo "Epsilon values: ${E_ARRAY[*]}"
echo "Output directory: $OUTDIR"
echo "$XYZ_OUT"
echo

#######################################
# Step 1: Generate topology files
#######################################

for E in "${E_ARRAY[@]}"; do
    echo "=== Generating topology file ==="
    echo 
    E_DIR="$OUTDIR/epsilon_$E"
    mkdir -p "$E_DIR"
    cd "$E_DIR"
    TOPO_OUT="topology_${FILE}_epsilon${E}.yaml"
    echo "  Running topology for pdb = $FILE at epsilon = $E → $TOPO_OUT"
    
    python3 ${WD}/pdb2xyz/__init__AH_Hakan.py \
        -i "${WD}/$PDB" \
	-o "$XYZ_OUT" \
	-t "$TOPO_OUT" \
	--pH "$PH" \
        --T "$T" \
        --epsilon "$E"
    cd $WD

    echo "Topology generation complete."
    echo
done

#######################################
# Step 2: Run duello scan
#######################################

for E in "${E_ARRAY[@]}"; do
    echo "=== Running duello scan ==="
    echo
    E_DIR="$OUTDIR/epsilon_$E"
    cd "$E_DIR"
    TOPO_IN="topology_${FILE}_epsilon${E}.yaml"
    SCAN_OUT="scan_epsilon${E}.dat"
    echo "  duello scan for epsilon = $E → $SCAN_OUT"
    
    duello scan --mol1 "$XYZ_OUT" \
                --mol2 "$XYZ_OUT" \
                --rmin 23 \
                --rmax 80 \
                --dr 0.5 \
		--resolution 0.7 \
		--cutoff 1000  \
                --top "$TOPO_IN"  \
		--molarity 0.115  \
		--temperature "$T" \
		--pmf "$SCAN_OUT"
    cd $WD

    echo "Duello scans complete."
    echo
done

echo "=== Done! ==="

