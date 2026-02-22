# ================================================
# Script to plot Enviflux assimilation experiments
# Ross Bannister, August 2025
# ================================================

PYTHON_DIR=$HOME/DataAssim/SourceSink/ToyModel/3D/graphics
ENVIFLUX_DIR=/media/ross/banny/EnviFlux

# Main experiments with realistic fluxes
#BACKGROUND_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Background
#OBS_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Observations/Bias_0.0_TotalCol
#TRUTH_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Truth
#CVT_FILE=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Covs/CVT.nc
#GROUP_DIR=$ENVIFLUX_DIR/Bias/Sept2025_experiments/Assim

# Blob experiments
BACKGROUND_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Background
OBS_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Observations/Unbiased_TC
#OBS_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Observations/Bias-12000_TC
TRUTH_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Truth
CVT_FILE=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Covs/CVT.nc
GROUP_DIR=$ENVIFLUX_DIR/RecoverBlobs/Sept2025_experiments/Assim


MIN_TRACER=-13.0
MAX_TRACER=13.0
MIN_FLUX=-0.00055
MAX_FLUX=0.00055

OUTPUT_TYPE=png


#TRACER_FACTORS="0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8"
#FLUX_FACTORS="0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8"
TRACER_FACTORS="1.0"
FLUX_FACTORS="100.0 10.0 1.0 0.1 0.01 0.001"

#EXPERIMENTS="wFac0.7_Bias0.0 wFac0.8_Bias0.0 wFac0.95_Bias0.0 wFac0.95_Bias0.0_homo wFac0.9_Bias0.0 wFac1.05_Bias0.0 wFac1.0_Bias0.0 wFac1.0_Bias0.0_homo wFac1.1_Bias0.0 wFac1.2_Bias0.0 wFac1.3_Bias0.0 wFac1.4_Bias0.0"
#EXPERIMENTS="wFac1.0_Bias0.0_ConstFluxStd_AssimAlso"
#EXPERIMENTS="wFac1.0_Bias0.0_SigmaFlux10 wFac1.0_Bias0.0_SigmaFlux100"
EXPERIMENTS="wFac1.0_Bias0.0"

for EXPERIMENT in $EXPERIMENTS
do
  for FAC_TRACER in $TRACER_FACTORS
  do
    for FAC_FLUX in $FLUX_FACTORS
    do

      echo $EXPERIMENT $FAC_TRACER $FAC_FLUX

      # Main experiments with realistic fluxes
      #EXP=TracerFac_${FAC_TRACER}__FluxFac_$FAC_FLUX
      #ASSIM_DIR=$GROUP_DIR/$EXPERIMENT/$EXP
      #BG_OBS_DIR=$GROUP_DIR/$EXPERIMENT/TracerFac_1.0__FluxFac_1.0

      # Blob experiments
      EXP=TracerFac${FAC_TRACER}_FluxFac${FAC_FLUX}
      ASSIM_DIR=$GROUP_DIR/${EXP}_${EXPERIMENT}
      BG_OBS_DIR=$ASSIM_DIR

      #echo $ASSIM_DIR

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
  done
done
echo "===== End of script =============================="
