CC=g++
CFLAGS_CVT=-lfftw3 -lm
#CFLAGS_INOUT=-lnetcdf_c++ -lnetcdf
CFLAGS_INOUT=-lnetcdf

# Following is for Ubuntu
LIB=-L/usr/lib/x86_64-linux-gnu
# Following is for UoR cluster (and can't use the ARRBOUNDCHECK)
#LIB=-L/usr/lib/gcc/x86_64-redhat-linux/4.8.2:/home/users/sws98rnb/Lib/lib:/usr/lib64

#ARRBOUNDCHECK=-fsanitize=address
OPTIMISE=-O3


# Linking stage
Master_AdjointTests_CVT.out: Master_AdjointTests_CVT.o cvt.o cvt_basic.o misc.o inout.o
	$(CC) -o Master_AdjointTests_CVT.out Master_AdjointTests_CVT.o cvt.o cvt_basic.o misc.o inout.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_AdjointTests_SemiLagrangian.out: Master_AdjointTests_SemiLagrangian.o interpolate.o misc.o inout.o cvt_basic.o trajectory.o transport_model.o
	$(CC) -o Master_AdjointTests_SemiLagrangian.out Master_AdjointTests_SemiLagrangian.o interpolate.o misc.o inout.o cvt_basic.o trajectory.o transport_model.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_AdjointTests_SL+ObsOp.out: Master_AdjointTests_SL+ObsOp.o interpolate.o misc.o inout.o cvt_basic.o trajectory.o transport_model.o ObsOperators.o
	$(CC) -o Master_AdjointTests_SL+ObsOp.out Master_AdjointTests_SL+ObsOp.o interpolate.o misc.o inout.o cvt_basic.o trajectory.o transport_model.o ObsOperators.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_InverseTests_CVT.out: Master_InverseTests_CVT.o cvt.o cvt_basic.o misc.o inout.o
	$(CC) -o Master_InverseTests_CVT.out Master_InverseTests_CVT.o cvt.o cvt_basic.o misc.o inout.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_ImpliedCov.out: Master_ImpliedCov.o cvt.o cvt_basic.o misc.o inout.o
	$(CC) -o Master_ImpliedCov.out Master_ImpliedCov.o cvt.o cvt_basic.o misc.o inout.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_MakeWinds.out: Master_MakeWinds.o cvt_basic.o misc.o inout.o program_arguments.o interpolate.o
	$(CC) -o Master_MakeWinds.out Master_MakeWinds.o cvt_basic.o misc.o inout.o program_arguments.o interpolate.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_TestReadWrite.out: Master_TestReadWrite.o cvt_basic.o misc.o inout.o
	$(CC) -o Master_TestReadWrite.out Master_TestReadWrite.o cvt_basic.o misc.o inout.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_ForwardTraj.out: Master_ForwardTraj.o trajectory.o interpolate.o program_arguments.o misc.o cvt_basic.o inout.o
	$(CC) -o Master_ForwardTraj.out Master_ForwardTraj.o trajectory.o interpolate.o program_arguments.o misc.o cvt_basic.o inout.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_SemiLagrangian.out: Master_SemiLagrangian.o program_arguments.o misc.o inout.o transport_model.o trajectory.o interpolate.o cvt_basic.o
	$(CC) -o Master_SemiLagrangian.out Master_SemiLagrangian.o program_arguments.o misc.o inout.o transport_model.o trajectory.o interpolate.o cvt_basic.o  $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_MakeFields.out: Master_MakeFields.o program_arguments.o misc.o inout.o cvt_basic.o
	$(CC) -o Master_MakeFields.out Master_MakeFields.o program_arguments.o misc.o inout.o cvt_basic.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_ReplaceFirstField.out: Master_ReplaceFirstField.o program_arguments.o misc.o inout.o cvt_basic.o
	$(CC) -o Master_ReplaceFirstField.out Master_ReplaceFirstField.o program_arguments.o misc.o inout.o cvt_basic.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_GenerateObs.out: Master_GenerateObs.o program_arguments.o misc.o inout.o transport_model.o ObsOperators.o cvt_basic.o trajectory.o interpolate.o
	$(CC) -o Master_GenerateObs.out Master_GenerateObs.o program_arguments.o misc.o inout.o transport_model.o ObsOperators.o cvt_basic.o trajectory.o interpolate.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_AdjointTestsObOps.out: Master_AdjointTestsObOps.o inout.o cvt_basic.o misc.o cvt.o ObsOperators.o
	$(CC) -o Master_AdjointTestsObOps.out Master_AdjointTestsObOps.o inout.o cvt_basic.o misc.o cvt.o ObsOperators.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_MakeBG.out: Master_MakeBG.o program_arguments.o inout.o cvt_basic.o misc.o cvt.o
	$(CC) -o Master_MakeBG.out Master_MakeBG.o program_arguments.o inout.o cvt_basic.o misc.o cvt.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_GradTest.out: Master_GradTest.o misc.o inout.o transport_model.o ObsOperators.o cvt_basic.o cvt.o trajectory.o interpolate.o PenAndGrad.o
	$(CC) -o Master_GradTest.out Master_GradTest.o misc.o inout.o transport_model.o ObsOperators.o cvt_basic.o cvt.o trajectory.o interpolate.o PenAndGrad.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_Assimilate.out: Master_Assimilate.o program_arguments.o misc.o inout.o transport_model.o ObsOperators.o cvt_basic.o cvt.o trajectory.o interpolate.o PenAndGrad.o
	$(CC) -o Master_Assimilate.out Master_Assimilate.o program_arguments.o misc.o inout.o transport_model.o ObsOperators.o cvt_basic.o cvt.o trajectory.o interpolate.o PenAndGrad.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_CalcError.out: Master_CalcError.o program_arguments.o misc.o inout.o cvt_basic.o
	$(CC) -o Master_CalcError.out Master_CalcError.o program_arguments.o misc.o inout.o cvt_basic.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)

Master_Invicat2Enviflux.out: Master_Invicat2Enviflux.o program_arguments.o misc.o inout.o cvt_basic.o interpolate.o
	$(CC) -o Master_Invicat2Enviflux.out Master_Invicat2Enviflux.o program_arguments.o misc.o inout.o cvt_basic.o interpolate.o $(LIB) $(CFLAGS_CVT) $(CFLAGS_INOUT) $(ARRBOUNDCHECK)


# Compilation stage
Master_AdjointTests_CVT.o: source.h Master_AdjointTests_CVT.cpp
	$(CC) -c Master_AdjointTests_CVT.cpp -o Master_AdjointTests_CVT.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_AdjointTests_SemiLagrangian.o: source.h Master_AdjointTests_SemiLagrangian.cpp
	$(CC) -c Master_AdjointTests_SemiLagrangian.cpp -o Master_AdjointTests_SemiLagrangian.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_AdjointTests_SL+ObsOp.o: source.h Master_AdjointTests_SL+ObsOp.cpp
	$(CC) -c Master_AdjointTests_SL+ObsOp.cpp -o Master_AdjointTests_SL+ObsOp.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_InverseTests_CVT.o: source.h Master_InverseTests_CVT.cpp
	$(CC) -c Master_InverseTests_CVT.cpp -o Master_InverseTests_CVT.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_ImpliedCov.o: source.h Master_ImpliedCov.cpp
	$(CC) -c Master_ImpliedCov.cpp -o Master_ImpliedCov.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_MakeWinds.o: source.h Master_MakeWinds.cpp
	$(CC) -c Master_MakeWinds.cpp -o Master_MakeWinds.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_TestReadWrite.o: source.h Master_TestReadWrite.cpp
	$(CC) -c Master_TestReadWrite.cpp -o Master_TestReadWrite.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_ForwardTraj.o: source.h Master_ForwardTraj.cpp
	$(CC) -c Master_ForwardTraj.cpp -o Master_ForwardTraj.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_SemiLagrangian.o: source.h Master_SemiLagrangian.cpp
	$(CC) -c Master_SemiLagrangian.cpp -o Master_SemiLagrangian.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_MakeFields.o: source.h Master_MakeFields.cpp
	$(CC) -c Master_MakeFields.cpp -o Master_MakeFields.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_ReplaceFirstField.o: source.h Master_ReplaceFirstField.cpp
	$(CC) -c Master_ReplaceFirstField.cpp -o Master_ReplaceFirstField.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_GenerateObs.o: source.h Master_GenerateObs.cpp
	$(CC) -c Master_GenerateObs.cpp -o Master_GenerateObs.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_AdjointTestsObOps.o: source.h Master_AdjointTestsObOps.cpp
	$(CC) -c Master_AdjointTestsObOps.cpp -o Master_AdjointTestsObOps.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_MakeBG.o: source.h Master_MakeBG.cpp
	$(CC) -c Master_MakeBG.cpp -o Master_MakeBG.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_GradTest.o: source.h Master_GradTest.cpp
	$(CC) -c Master_GradTest.cpp -o Master_GradTest.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_Assimilate.o: source.h Master_Assimilate.cpp
	$(CC) -c Master_Assimilate.cpp -o Master_Assimilate.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_CalcError.o: source.h Master_CalcError.cpp
	$(CC) -c Master_CalcError.cpp -o Master_CalcError.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

Master_Invicat2Enviflux.o: source.h Master_Invicat2Enviflux.cpp
	$(CC) -c Master_Invicat2Enviflux.cpp -o Master_Invicat2Enviflux.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

cvt.o: source.h cvt.cpp
	$(CC) -c cvt.cpp -o cvt.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

cvt_basic.o: source.h cvt_basic.cpp
	$(CC) -c cvt_basic.cpp -o cvt_basic.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

misc.o: source.h misc.cpp
	$(CC) -c misc.cpp -o misc.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

inout.o: source.h inout.cpp
	$(CC) -c inout.cpp -o inout.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

transport_model.o: source.h transport_model.cpp
	$(CC) -c transport_model.cpp -o transport_model.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

trajectory.o: source.h trajectory.cpp
	$(CC) -c trajectory.cpp -o trajectory.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

interpolate.o: source.h interpolate.cpp
	$(CC) -c interpolate.cpp -o interpolate.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

ObsOperators.o: source.h ObsOperators.cpp
	$(CC) -c ObsOperators.cpp -o ObsOperators.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

PenAndGrad.o: source.h PenAndGrad.cpp
	$(CC) -c PenAndGrad.cpp -o PenAndGrad.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

program_arguments.o: source.h program_arguments.cpp
	$(CC) -c program_arguments.cpp -o program_arguments.o -I. $(ARRBOUNDCHECK) $(OPTIMISE)

# Tidy up
clean:
	rm *.o *.out
