#================================================================================
# Compile and link options for ubuntu
#===================================================================================
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
   bits=-m64
else
   bits=-m32
endif


# Have been using these for ages, but probably don't need all options
#LinkOpts=-fcheck=bounds -I/usr/local/include -L/usr/local/lib -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 $(bits) -lnetcdff -lrefblas -fcheck=bounds
#CompileOpts=-fcheck=bounds -I/usr/local/include -I/usr/include -L/usr/local/lib -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 $(bits) -lnetcdff -lrefblas -fdefault-real-8 -fcheck=bounds

# For code using shtools
LinkOpts_shtools=-fcheck=bounds -lSHTOOLS -O3 $(bits) -fdefault-real-8 -L/usr/local/lib
CompileOpts_shtools=-fcheck=bounds -I/usr/local/include -O3 $(bits) -fdefault-real-8

# For basic compilation
CompileOpts_modules=-fcheck=bounds -O3 $(bits) -fdefault-real-8

# For compilation using netcdf and fft
CompileOpts_netcdffft=-fcheck=bounds -O3 $(bits) -fdefault-real-8 -I/usr/include

# Main linking options
LinkOpts=-fcheck=bounds -O3 $(bits) -fdefault-real-8 -L/usr/local/lib -lnetcdff -lfftpack5.1d -llapack -lrefblas


# Linking stage
AdjointTestsMain.out: main_data.o cparam.o onedvar_data.o AdjointTestsMain.o cvt.o Readwrite4cvt.o
	gfortran -o $@ $^ $(LinkOpts)

Calibration.out: main_data.o onedvar_data.o cparam.o Calibration.o cvt.o Readwrite4cvt.o
	gfortran -o $@ $^ $(LinkOpts)

ImpliedCov.out: main_data.o cparam.o onedvar_data.o ImpliedCov.o cvt.o Readwrite4cvt.o
	gfortran -o $@ $^ $(LinkOpts)

genCVT_ALP_data.out: main_data.o onedvar_data.o genCVT_ALP_data.o
	gfortran -o $@ $^ $(LinkOpts_shtools)


# Compiling stage
main_data.o: main_data.f90
	gfortran -c $< $(CompileOpts_modules)

cparam.o: cparam.f90
	gfortran -c $< $(CompileOpts_modules)

onedvar_data.o: onedvar_data.f90 main_data.o
	gfortran -c $< $(CompileOpts_modules)

AdjointTestsMain.o: AdjointTestsMain.f90 main_data.o cparam.o onedvar_data.o
	gfortran -c $< $(CompileOpts_netcdffft)

cvt.o: cvt.f90 main_data.o cparam.o onedvar_data.o
	gfortran -c $< $(CompileOpts_netcdffft)

Calibration.o: Calibration.f90 main_data.o onedvar_data.o
	gfortran -c $< $(CompileOpts_netcdffft)

Readwrite4cvt.o: Readwrite4cvt.f90 main_data.o onedvar_data.o
	gfortran -c $< $(CompileOpts_netcdffft)

ImpliedCov.o: ImpliedCov.f90 main_data.o onedvar_data.o
	gfortran -c $< $(CompileOpts_netcdffft)

genCVT_ALP_data.o: genCVT_ALP_data.f90 main_data.o onedvar_data.o
	gfortran -c $< $(CompileOpts_shtools)

clean:
	rm *.o *.mod *.out
