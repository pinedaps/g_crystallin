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
# Step 2: Resolution sensitivity analysis
#######################################

echo "=== Running resolution sensitivity analysis ==="

# Resolution range
RES_START=100.0      # initial coarse resolution
RES_MIN=0.1        # lowest allowed resolution
RES_STEP=1       # decrement
TOLERANCE=0.10     # 10% relative change threshold for stopping

for T in "${T_ARRAY[@]}"; do
    echo "  Sensitivity analysis for T = $T"

    TOPO_IN="${TOPO_DIR}/topology_${FILE}_T${T}.yaml"

    prev_B2=""
    res=$RES_START

    while (( $(echo "$res >= $RES_MIN" | bc -l) )); do
        OUT_PREFIX="${SCAN_DIR}/res_sens_T${T}_res${res}"
        SCAN_OUT="${OUT_PREFIX}.dat"
        JSON_OUT="${OUT_PREFIX}.json"
	TIME_LOG="${OUT_PREFIX}_time.log"

        echo "    → Testing resolution = $res"
	echo "    duello scan --mol1 "$XYZ_OUT" \
                            --mol2 "$XYZ_OUT" \
                            --rmin 20 \
                            --rmax 40 \
                            --dr 1 \
                            --resolution "$res" \
                            --cutoff 80 \
                            --top "$TOPO_IN" \
                            --molarity 0.115 \
                            --temperature "$T" \
                            --pmf "$SCAN_OUT" "	
	/usr/bin/time -f "elapsed=%e user=%U sys=%S" \
    		-o "$TIME_LOG" \
        	duello scan --mol1 "$XYZ_OUT" \
                    	    --mol2 "$XYZ_OUT" \
                    	    --rmin 20 \
                    	    --rmax 80 \
                    	    --dr 0.5 \
                     	    --resolution "$res" \
                    	    --cutoff 80 \
                    	    --top "$TOPO_IN" \
                    	    --molarity 0.115 \
                    	    --temperature "$T" \
                    	    --pmf "$SCAN_OUT"

	# Extract timing info
	
	ELAPSED=$(grep "elapsed" "$TIME_LOG" | sed 's/elapsed=\([0-9.]*\).*/\1/')
	USER=$(grep "user" "$TIME_LOG" | sed 's/.*user=\([0-9.]*\).*/\1/')
	SYS=$(grep "sys" "$TIME_LOG" | sed 's/.*sys=\([0-9.]*\)/\1/')

	echo "      CPU time: elapsed=${ELAPSED}s user=${USER}s sys=${SYS}s"
        
	# Extract B2 from JSON
        if [[ ! -f "$JSON_OUT" ]]; then
            echo "      WARNING: $JSON_OUT not produced — skipping"
            break
        fi

        B2=$(jq -r '.B2' "$JSON_OUT")
	echo "$JSON_OUT"	
	if ! [[ "$B2" =~ ^-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$ ]]; then
    		echo "      ERROR: Non-numeric B2 value '$B2'. Cannot compute relative change."
    		break
	else
		echo " $B2 "
	fi
        
	if [[ -z "$prev_B2" ]]; then
            echo "      B2 = $B2 (initial reference)"
            prev_B2="$B2"
        else
            # Compute relative change |ΔB2| / |B2_prev|
            diff=$(echo "($B2 - $prev_B2)" | bc -l)
            rel_change=$(echo "scale=6; sqrt(($diff * $diff)) / ($prev_B2 == 0 ? 1 : $prev_B2)" | bc -l)

            echo "      B2 = $B2; relative change = $rel_change"
            
	    # stop if variation < 10%
            below_tol=$(echo "$rel_change < $TOLERANCE" | bc -l)
            if [[ "$below_tol" -eq 1 ]]; then
                echo "      ✔ Converged: variation < 10%. Stopping."
                break
            fi

            prev_B2="$B2"
        fi

        # decrement resolution
	res=$(printf "%.2f" "$(echo "$res - $RES_STEP" | bc -l)")
    done
done

echo "Resolution sensitivity analysis complete."

#######################################
# Step 3: Plot results
#######################################

echo "=== Plotting results ==="

python3 plots/plot_potential.py "${SCAN_DIR}/"

echo "Plots generated in: $PLOT_DIR"
echo
echo "=== Done! ==="

