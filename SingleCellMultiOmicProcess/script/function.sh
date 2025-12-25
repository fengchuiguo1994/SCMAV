function odir {
    local outputDir=$1
    if [ ! -d $outputDir ]; then
        echo "Create $outputDir"
        mkdir $outputDir
    fi
}

function submitjobOptionalSLURM {
    local slurmfile="$1"
    shift
    if [ ! -f "$slurmfile" ]; then
        echo "Error: SLURM file '$slurmfile' does not exist." >&2
        return 1
    fi
    local dependency=""
    if [ $# -gt 0 ]; then
        dependency="--dependency=afterany:$(IFS=:; echo "$*")"
    fi
    local jidnow
    if [ -z "$dependency" ]; then
        jidnow=$(sbatch "$slurmfile")
    else
        jidnow=$(sbatch "$dependency" "$slurmfile")
    fi
    echo $jidnow | awk '{print $NF}'
}

function submitjobOptionalLSF {
    local lsf_script="$1"
    shift 
    if [ ! -f "$lsf_script" ]; then
        echo "Error: LSF script '$lsf_script' does not exist." >&2
        return 1
    fi
    local dependency=""
    if [ $# -gt 0 ]; then
        local dep_conditions=()
        for jobid in "$@"; do
            dep_conditions+=("done($jobid)")
        done
        dependency="-w \"$(IFS=' && '; echo "${dep_conditions[*]}")\""
    fi
    local jobid
    if [ -z "$dependency" ]; then
        jobid=$(bsub < "$lsf_script" | awk '{print $2}' | sed 's/[<>]//g')
    else
        jobid=$(eval "bsub $dependency < \"$lsf_script\"" | awk '{print $2}' | sed 's/[<>]//g')
    fi
    echo $jobid
}

function submitjobOptionalSGE {
    local sge_script="$1"
    shift
    if [ ! -f "$sge_script" ]; then
        echo "Error: SGE script '$sge_script' does not exist." >&2
        return 1
    fi
    local dependency=""
    if [ $# -gt 0 ]; then
        dependency="-hold_jid $(IFS=,; echo "$*")"
    fi
    local jobid
    if [ -z "$dependency" ]; then
        jobid=$(qsub "$sge_script")
    else
        jobid=$(qsub $dependency "$sge_script")
    fi
    echo $jobid
}

function submitjobOptionalTORQUE {
    local pbs_script="$1"
    shift 
    if [ ! -f "$pbs_script" ]; then
        echo "Error: PBS script '$pbs_script' does not exist." >&2
        return 1
    fi
    local dependency=""
    if [ $# -gt 0 ]; then
        dependency="-W depend=afterany:$(IFS=:; echo "$*")"
    fi
    local jobid
    if [ -z "$dependency" ]; then
        jobid=$(qsub "$pbs_script")
    else
        jobid=$(qsub $dependency "$pbs_script")
    fi
    echo $jobid | awk -F. '{print $1}'
}

function submitjobOptionalShell {
    local script="$1"
    shift
    if [ ! -f "$script" ]; then
        echo "Error: Script '$script' not found!" >&2
        return 1
    fi
    if [ $# -gt 0 ]; then
        echo "Waiting for dependencies: $*"
        for pid in "$@"; do
            wait "$pid" 2>/dev/null || echo "Warning: PID $pid not found"
        done
    fi
    bash "$script" &
    local current_pid=$!
    echo "$current_pid"
}

function universal_submit {
    local script=$1
    shift
    local deps=("$@")
    if which sbatch >/dev/null; then
        submitjobOptionalSLURM "$script" "${deps[@]}"
    elif which bsub >/dev/null; then
        submitjobOptionalLSF "$script" "${deps[@]}"
    elif which qsub >/dev/null; then
        if qsub --version 2>&1 | grep -q PBS; then
            submitjobOptionalTORQUE "$script" "${deps[@]}"
        else
            submitjobOptionalSGE "$script" "${deps[@]}"
        fi
    else
        echo "Error: No supported job scheduler found" >&2
        return 1
    fi
}
