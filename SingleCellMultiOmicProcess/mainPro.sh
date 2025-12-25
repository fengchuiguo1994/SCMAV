#!/usr/bin/bash

set -e
config=$1
inout=$2
source $config

SCRIPT_DIR=$(dirname "$0")
source $SCRIPT_DIR/script/function.sh

case $type in 
scRNA)
    echo 11
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Rbc=$(realpath $Rbc)
        Rrd=$(realpath $Rrd)
        outdir=$(realpath $outdir)
        echo $task
        case $task in
        local)
            echo 22
            bash $SCRIPT_DIR/script/process_scRNA_local.sh $Rbc $Rrd $outdir $outprefix $config
        ;;
        slurm)
            echo 33
            bash $SCRIPT_DIR/script/process_scRNA_slurm.sh $Rbc $Rrd $outdir $outprefix $config
        ;;
        *)
            echo "Invalid parameter task: $task"
            exit 1 
        ;;
        esac
    done
;;

scATAC)
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Dbc=$(realpath $Dbc)
        Drd1=$(realpath $Drd1)
        Drd2=$(realpath $Drd2)
        outdir=$(realpath $outdir)
        case $task in 
        local)
            echo 33
            bash $SCRIPT_DIR/script/process_scATAC_local.sh $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        slurm)
            echo 44
            bash $SCRIPT_DIR/script/process_scATAC_slurm.sh $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        *)
            echo "Invalid parameter task: $task"
            exit 1 
        ;;
        esac
    done
;;

scHi-C)
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        echo 22
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Dbc=$(realpath $Dbc)
        Drd1=$(realpath $Drd1)
        Drd2=$(realpath $Drd2)
        outdir=$(realpath $outdir)
        case $task in 
        local)
            echo 33
            bash $SCRIPT_DIR/script/process_scHi-C_local.sh $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        slurm)
            bash $SCRIPT_DIR/script/process_scHi-C_slurm.sh $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        *)
        echo "Invalid parameter task: $task"
        exit 1 
        ;;
        esac
    done
;;

scChIATAC)
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        echo 22
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Dbc=$(realpath $Dbc)
        Drd1=$(realpath $Drd1)
        Drd2=$(realpath $Drd2)
        outdir=$(realpath $outdir)
        case $task in 
        local)
            echo 33
            bash $SCRIPT_DIR/script/process_scChIATAC_local.sh $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        slurm)
            echo 33
            bash $SCRIPT_DIR/script/process_scChIATAC_slurm.sh $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        *)
        echo "Invalid parameter task: $task"
        exit 1 
        ;;
        esac
    done
;;

scRNA+scATAC)
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        echo 22
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Rbc=$(realpath $Rbc)
        Rrd=$(realpath $Rrd)
        Dbc=$(realpath $Dbc)
        Drd1=$(realpath $Drd1)
        Drd2=$(realpath $Drd2)
        outdir=$(realpath $outdir)
        case $task in 
        local)
            echo 33
            bash $SCRIPT_DIR/script/process_scRNA_scATAC_local.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        slurm)
            bash $SCRIPT_DIR/script/process_scRNA_scATAC_slurm.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        *)
        echo "Invalid parameter task: $task"
        exit 1 
        ;;
        esac
    done
;;

scRNA+scHi-C)
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        echo 22
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Rbc=$(realpath $Rbc)
        Rrd=$(realpath $Rrd)
        Dbc=$(realpath $Dbc)
        Drd1=$(realpath $Drd1)
        Drd2=$(realpath $Drd2)
        outdir=$(realpath $outdir)
        case $task in 
        local)
            echo 33
            bash $SCRIPT_DIR/script/process_scRNA_scHi-C_local.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        slurm)
            bash $SCRIPT_DIR/script/process_scRNA_scHi-C_slurm.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        *)
        echo "Invalid parameter task: $task"
        exit 1 
        ;;
        esac
    done
;;

scRNA+scChIATAC)
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        echo 22
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Rbc=$(realpath $Rbc)
        Rrd=$(realpath $Rrd)
        Dbc=$(realpath $Dbc)
        Drd1=$(realpath $Drd1)
        Drd2=$(realpath $Drd2)
        outdir=$(realpath $outdir)
        case $task in 
        local)
            echo 33
            bash $SCRIPT_DIR/script/process_scRNA_scChIATAC_local.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        slurm)
            bash $SCRIPT_DIR/script/process_scRNA_scChIATAC_slurm.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        *)
        echo "Invalid parameter task: $task"
        exit 1 
        ;;
        esac
    done
;;

ChAIR)
    cat $inout | while read Rbc Rrd Dbc Drd1 Drd2 outdir outprefix; do 
        odir $outdir && odir $outdir/jobdir
        jobdir=$outdir/jobdir
        Rbc=$(realpath $Rbc)
        Rrd=$(realpath $Rrd)
        Dbc=$(realpath $Dbc)
        Drd1=$(realpath $Drd1)
        Drd2=$(realpath $Drd2)
        outdir=$(realpath $outdir)
        case $task in 
        local)
            echo 33
            bash $SCRIPT_DIR/script/process_scRNA_scChIATAC_local.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        slurm)
            bash $SCRIPT_DIR/script/process_scRNA_scChIATAC_slurm.sh $Rbc $Rrd $Dbc $Drd1 $Drd2 $outdir $outprefix $config
        ;;
        *)
        echo "Invalid parameter task: $task"
        exit 1 
        ;;
        esac
    done
;;

*)
echo "Invalid parameter type: $type"
exit 1  # 失败退出
;;
esac
