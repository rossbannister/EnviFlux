# ================================================
# Script to plot Enviflux assimilation experiments
# Ross Bannister, August 2025
# ================================================

PYTHON_DIR=$HOME/DataAssim/SourceSink/ToyModel/3D/graphics
ENVIFLUX_DIR=/media/ross/banny/EnviFlux

# Main experiments with realistic fluxes
BACKGROUND_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Background
TRUTH_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Truth
CVT_FILE=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Covs/CVT.nc
GROUP_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Assim

# Blob experiments
#BACKGROUND_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Background
#TRUTH_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Truth
#CVT_FILE=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Covs/CVT.nc
#GROUP_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Assim

#MIN_TRACER=-13.0
#MAX_TRACER=13.0
#MIN_FLUX=-0.00055
#MAX_FLUX=0.00055
MIN_TRACER=0.0
MAX_TRACER=0.0
MIN_FLUX=0.0
MAX_FLUX=0.0


#EXPERIMENTS="wFac0.6_Bias0.0 wFac0.6_Bias0.0_z0 wFac0.7_Bias0.0 wFac0.7_Bias0.0_z0 wFac0.8_Bias0.0 wFac0.8_Bias0.0_z0 wFac0.9_Bias0.0 wFac0.9_Bias0.0_z0 wFac0.95_Bias0.0 wFac0.95_Bias0.0_z0 wFac1.05_Bias0.0 wFac1.05_Bias0.0_z0 wFac1.1_Bias0.0 wFac1.1_Bias0.0_z0 wFac1.2_Bias0.0 wFac1.2_Bias0.0_z0 wFac1.3_Bias0.0 wFac1.3_Bias0.0_z0 wFac1.4_Bias0.0 wFac1.4_Bias0.0_z0"
#OBS_DIR="Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0 Bias_0.0_TotalCol Bias_0.0_z0"

EXPERIMENTS="wFac1.3_Bias0.0 wFac1.3_Bias0.0_z0"
OBS_DIR="Bias_0.0_TotalCol Bias_0.0_z0"

OUTPUT_TYPE=png

# Convert lists into arrays
EXPERIMENTS_ARR=($EXPERIMENTS)
OBS_DIR_ARR=($OBS_DIR)

NItemsE=${#EXPERIMENTS_ARR[@]}
NItemsO=${#OBS_DIR_ARR[@]}
echo $NItemsE $NItemsO

if [ $NItemsE -eq $NItemsO ]; then

  let N=$NItemsE-1

  touch ProcessingDirs.txt

  for NUM in $(seq 0 $N)
  do
    # Extract the NUMth items
    echo $NUM
    EXP=${EXPERIMENTS_ARR[$NUM]}
    OB=${OBS_DIR_ARR[$NUM]}
    echo $EXP $OB >> ProcessingDirs.txt
    ASSIM_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Assim/$EXP
    OBS_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Observations/$OB
    BG_OBS_DIR=$ASSIM_DIR

    python3 $PYTHON_DIR/Plot_PriorPostInc.py \
      $ASSIM_DIR \
      $BACKGROUND_DIR \
      $OBS_DIR \
      $TRUTH_DIR \
      $BG_OBS_DIR \
      $CVT_FILE \
      $MIN_TRACER \
      $MAX_TRACER \
      $MIN_FLUX \
      $MAX_FLUX \
      $OUTPUT_TYPE

  done
fi


echo "===== End of script =============================="
