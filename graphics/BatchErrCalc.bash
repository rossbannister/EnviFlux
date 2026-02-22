BIN_DIR=/home/ross/DataAssim/SourceSink/ToyModel/3D/src

BASE_DIR=/media/ross/banny/EnviFlux/Bias/Sept2025_experiments
BACKGROUND_DIR=$BASE_DIR/Background
TRUTH_DIR=$BASE_DIR/Truth
CVT_DIR=$BASE_DIR/Covs

EXPERIMENTS="wFac0.6_Bias0.0 wFac1.4_Bias0.0" # wFac0.6_Bias0.0_z0 wFac0.7_Bias0.0 wFac0.7_Bias0.0_z0 wFac0.8_Bias0.0 wFac0.8_Bias0.0_z0 wFac0.9_Bias0.0 wFac0.9_Bias0.0_z0 wFac0.95_Bias0.0 wFac0.95_Bias0.0_z0 wFac1.05_Bias0.0 wFac1.05_Bias0.0_z0 wFac1.1_Bias0.0 wFac1.1_Bias0.0_z0 wFac1.2_Bias0.0 wFac1.2_Bias0.0_z0 wFac1.4_Bias0.0_z0 wFac1.3_Bias0.0_z0 wFac1.3_Bias0.0 wFac1.0_Bias0.00"

echo "" > ErrorsOfExps.dat

for EXP in $EXPERIMENTS; do
  $BIN_DIR/Master_CalcError.out \
    $BASE_DIR/Assim/$EXP/Anal.nc \
    $TRUTH_DIR/Truth.nc \
    $BACKGROUND_DIR/Background.nc \
    $CVT_DIR/CVT.nc \
    $BASE_DIR/Assim/$EXP/ErrStats.dat
  echo $EXP >> ErrorsOfExps.dat
  cat $BASE_DIR/Assim/$EXP/ErrStats.dat >> ErrorsOfExps.dat
  echo "==========================================================" >> ErrorsOfExps.dat
done
