PROGRAM Calibrate

! This program reads in methane forecasts in order to calibrate part of the CVT used in INVICAT.
! The output is a cov file (in nc format).

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields

USE onedvar_data, ONLY: nforecasts_calibtracer,        &
                        forecasts_calib_tracer_dir,    &
                        L,                             &
                        ALPsize,                       &
                        CVT_std,                       &
                        CVT_HorizCors,                 &
                        CVT_VertCors,                  &
                        CVT_TemporalCors,              &
                        n1d,                           &
                        GLatEq


USE cparam,       ONLY: deg2rad, Re, RCH4, Rdry


IMPLICIT none

INCLUDE "Readwrite4cvt.interface"
INCLUDE "cvt.interface"

! Local variables
INTEGER           :: forecastno, lon, lat, lev, lev1, lev2, mode, err, ll, mm
INTEGER           :: phase, month, sep, it
REAL              :: tracer(1:ylat, 1:ylon, 1:ylev, 1:nforecasts_calibtracer)
REAL              :: tracer1(1:ylat, 1:ylon, 1:ylev, 1:nforecasts_calibtracer)
REAL              :: tracer_verttrans(1:ylat, 1:ylon, 1:ylev)
REAL              :: tracer_verttrans_sh(1:L+1, 1:2*L+1, 1:ylev)
REAL              :: tracer_hortrans(1:2, 0:L, 0:L, 1:ylev, 1:nforecasts_calibtracer)
CHARACTER(len=48) :: filenames(1:nforecasts_calibtracer), fname
REAL              :: alt(1:ylat, 1:ylon, 0:ylev+1)
REAL              :: p(1:ylat, 1:ylon, 1:ylev)
REAL              :: T(1:ylat, 1:ylon, 1:ylev)
REAL              :: mass_ch4(1:nforecasts_calibtracer)
REAL              :: ppbv_ch4(1:nforecasts_calibtracer)
REAL              :: volume_air(1:nforecasts_calibtracer)
REAL              :: volume_ch4(1:nforecasts_calibtracer)
REAL              :: mass_air(1:nforecasts_calibtracer)
REAL              :: level_mean_ppbv_ch4(1:ylev)
REAL              :: totalmass_ch4, totalvolume_ch4, totalvolume, totalmass, level_mean, level_volume
REAL              :: dx, dy, dz, thisalt, total
REAL              :: thisvolume, thisvolume_ch4, thismass_ch4, thismass, thisdensity
REAL              :: mean_t, mean_ch4, mean_tpert_ch4pert, mean_tpert_tpert, global_tendency
REAL              :: nforecasts_real
REAL              :: temporal_mean(1:ylat, 1:ylon, 1:ylev)
REAL              :: vert_evals(1:ylev), Workv(1:3*ylev), psi, lambda, Workt(1:3*nmonth)
INTEGER           :: time1, time2, tdiff, field, indexll
REAL              :: dist(1:L+1), tracer_cov_impl(1:L+1, 1:ylev), dlen, num, den, corfn, derivL, dL, deltaL
REAL              :: sDt

! Function
REAL              :: CorrelationFn


! ==================================================================================

! Allocate and initialise CVT
CALL AllocateCVT (.TRUE.)
CALL InitializeHorizTrans ()


! ==================================================================================

! ----------------------------------------------------------------------------------
! Deal with the tracer covariance model
! ----------------------------------------------------------------------------------

nforecasts_real = REAL(nforecasts_calibtracer)

! ##### USER INPUT #####
! Set the filenames
filenames(1)  = 'ptom_t021_1995010206_prior.nc'
filenames(2)  = 'ptom_t021_1996010206_prior.nc'
filenames(3)  = 'ptom_t021_1997010206_prior.nc'
filenames(4)  = 'ptom_t021_1998010206_prior.nc'
filenames(5)  = 'ptom_t021_1999010206_prior.nc'
filenames(6)  = 'ptom_t021_2000010206_prior.nc'
filenames(7)  = 'ptom_t021_2001010206_prior.nc'
filenames(8)  = 'ptom_t021_2002010206_prior.nc'
filenames(9)  = 'ptom_t021_2003010206_prior.nc'
filenames(10) = 'ptom_t021_2004010206_prior.nc'

! Set the horizontal lengthscales of the tracer errors if required (in metres).
! (Note setting them all to zero here will result in these being found from data.
!  Check they are realistic afterwards!)
CVT_HorizCors % lengthscale_tracer(1:ylev) = 400000.0 ! in m
! ######################

! Note that level 1 is the top
! Set the altitudes of the very bottomost level
alt(1:ylat, 1:ylon, ylev+1) = 0.0

! Initialise the level means to zero
level_mean_ppbv_ch4(1:ylev) = 0.0


! Open the output file
OPEN (12, file='GlobalMethane.dat')
WRITE (12,*) '#    Year No, total mass,               mass(methane),            total volume,' // &
             '             volume(methane),          mean ppbv'


! ----------------------------------------------------------------------------------
! Read in the forecasts
! ----------------------------------------------------------------------------------
DO forecastno = 1, nforecasts_calibtracer

  PRINT *, 'Reading methane tracer number ', forecastno
  CALL read_concentration_fc ( TRIM(forecasts_calib_tracer_dir) // '/' // filenames(forecastno), &
                               tracer(1:ylat, 1:ylon, 1:ylev, forecastno),                       &
                               CVT_HorizCors % TOMCAT_longs(1:ylon),                             &
                               CVT_HorizCors % TOMCAT_lats(1:ylat),                              &
                               alt(1:ylat, 1:ylon, 1:ylev),                                      &
                               p(1:ylat, 1:ylon, 1:ylev),                                        &
                               T(1:ylat, 1:ylon, 1:ylev))

  ! Set the altitudes of the very topmost level
  DO lat = 1, ylat
    DO lon = 1, ylon
      dz = alt(lat,lon,1) - alt(lat,lon,2)
      alt(1:ylat, 1:ylon, 0) = alt(1:ylat, 1:ylon, 1) + dz
    END DO
  END DO

  ! Set the boundary longitudes and latitudes
  IF (forecastno == 1) THEN
    CVT_HorizCors % TOMCAT_longs(0)      = CVT_HorizCors % TOMCAT_longs(1) -    &
                                          (CVT_HorizCors % TOMCAT_longs(2) - CVT_HorizCors % TOMCAT_longs(1))
    CVT_HorizCors % TOMCAT_longs(ylon+1) = CVT_HorizCors % TOMCAT_longs(ylon) + &
                                           CVT_HorizCors % TOMCAT_longs(ylon) - CVT_HorizCors % TOMCAT_longs(ylon-1)
    CVT_HorizCors % TOMCAT_lats(0)       = CVT_HorizCors % TOMCAT_lats(1) -    &
                                          (CVT_HorizCors % TOMCAT_lats(2) - CVT_HorizCors % TOMCAT_lats(1))
    CVT_HorizCors % TOMCAT_lats(ylat+1)  = CVT_HorizCors % TOMCAT_lats(ylat) + &
                                           CVT_HorizCors % TOMCAT_lats(ylat) - CVT_HorizCors % TOMCAT_lats(ylat-1)
    !print *, 'Here are the latitudes'
    !print *, CVT_HorizCors % TOMCAT_lats(0:ylat+1)
  END IF
  !PRINT *, 'Methane at 28, 60, 58     ', tracer(28,60,58,forecastno)
  !PRINT *, 'Temperature at 28, 60, 58 ', T(28,60,58)
  !PRINT *, 'Pressure at 28, 60, 58    ', p(28,60,58)
  !PRINT *, 'Altitude at 28, 60, 58    ', alt(28,60,58)


  ! Compute the total methane mass for this forecast (see formula in notes)
  totalmass_ch4   = 0.0
  totalvolume_ch4 = 0.0
  totalvolume     = 0.0
  totalmass       = 0.0
  DO lev = 1, ylev
    level_mean   = 0.0
    level_volume = 0.0
    DO lat = 1, ylat
      dy = deg2rad * (CVT_HorizCors % TOMCAT_lats(lat-1) - CVT_HorizCors % TOMCAT_lats(lat+1)) / 2.0
      DO lon = 1, ylon
        dx = deg2rad * (CVT_HorizCors % TOMCAT_longs(lon+1) - CVT_HorizCors % TOMCAT_longs(lon-1)) / 2.0
        dz              = (alt(lat,lon,lev-1) - alt(lat,lon,lev+1)) / 2.0
        thisalt         = Re + alt(lat,lon,lev)
        thisvolume      = thisalt * thisalt * cos(CVT_HorizCors % TOMCAT_lats(lat)*deg2rad) * dx * dy * dz
        thisdensity     = p(lat,lon,lev) / (T(lat,lon,lev) * Rdry)
        thismass        = thisvolume * thisdensity
        thisvolume_ch4  = thisvolume * tracer(lat,lon,lev,forecastno) / 1000000000.0
        thismass_ch4    = p(lat,lon,lev) * thisvolume_ch4 / (T(lat,lon,lev) * RCH4)

        totalmass_ch4   = totalmass_ch4 + thismass_ch4
        totalvolume     = totalvolume + thisvolume
        totalvolume_ch4 = totalvolume_ch4 + thisvolume_ch4
        totalmass       = totalmass + thismass

        level_mean      = level_mean + thisvolume_ch4
        level_volume    = level_volume + thisvolume

        IF ((lat == 22) .AND. (lon == 30) .AND. (lev == ylev)) PRINT *, 'Example ground air density ', thisdensity

      END DO
    END DO
    level_mean_ppbv_ch4(lev) = 1000000000.0 * level_mean / level_volume
  END DO

  mass_ch4(forecastno)   = totalmass_ch4
  ppbv_ch4(forecastno)   = 1000000000.0 * totalvolume_ch4 / totalvolume
  volume_air(forecastno) = totalvolume
  mass_air(forecastno)   = totalmass
  volume_ch4(forecastno) = totalvolume_ch4

  PRINT *, 'Mass of CH4 for this forecast = ', mass_ch4(forecastno)
  PRINT *, 'Total volume of atmosphere    = ', volume_air(forecastno)
  PRINT *, 'Mean vol mixing ratio         = ', ppbv_ch4(forecastno)
  PRINT *, 'Total mass of the air         = ', mass_air(forecastno)
  PRINT *, 'Volume of CH4                 = ', volume_ch4(forecastno)
  !PRINT *, 'Time and level averaged ppbv  = ', level_mean_ppbv_ch4(1:ylev)


  WRITE (12,*) forecastno, mass_air(forecastno), mass_ch4(forecastno), volume_air(forecastno),  &
               volume_ch4(forecastno), ppbv_ch4(forecastno)

END DO

CLOSE (12)


! ----------------------------------------------------------------------------------
! Do a linear regression to find the best fit rate of change of global mean methane
! ----------------------------------------------------------------------------------

PRINT *, 'Doing linear regression'
! Treat time as in years
mean_t   = 0.0
mean_ch4 = 0.0
DO forecastno = 1, nforecasts_calibtracer
  mean_t   = mean_t   + real(forecastno)
  mean_ch4 = mean_ch4 + ppbv_ch4(forecastno)
END DO
mean_t   = mean_t   / nforecasts_real
mean_ch4 = mean_ch4 / nforecasts_real
PRINT *, 'Mean t   = ', mean_t
PRINT *, 'Mean ch4 = ', mean_ch4

mean_tpert_ch4pert = 0.0
mean_tpert_tpert   = 0.0
DO forecastno = 1, nforecasts_calibtracer
  mean_tpert_ch4pert = mean_tpert_ch4pert + (real(forecastno) - mean_t) * (ppbv_ch4(forecastno) - mean_ch4)
  mean_tpert_tpert   = mean_tpert_tpert   + (real(forecastno) - mean_t) * (real(forecastno) - mean_t)
END DO
mean_tpert_ch4pert = mean_tpert_ch4pert / nforecasts_real
mean_tpert_tpert   = mean_tpert_tpert  / nforecasts_real
global_tendency    = mean_tpert_ch4pert / mean_tpert_tpert
PRINT *, 'Global tendency = ', global_tendency, ' ppbv/year'


! ----------------------------------------------------------------------------------
! Detrend the data according to the formula in the notes (use option 2)
! ----------------------------------------------------------------------------------

PRINT *, 'Detrending the data'

DO forecastno = 1, nforecasts_calibtracer
  DO lat = 1, ylat
    DO lon = 1, ylon
      DO lev = 1, ylev
        tracer(lat,lon,lev,forecastno) = tracer(lat,lon,lev,forecastno) -      &
                                         global_tendency * level_mean_ppbv_ch4(lev) / mean_ch4
      END DO
    END DO
  END DO
END DO


! ----------------------------------------------------------------------------------
! Output the detrended data (for inspection)
! ----------------------------------------------------------------------------------

PRINT *, 'Writing time evolving (detrended) fields'
CALL Write_time_evolving_fields (TRIM(forecasts_calib_tracer_dir) // '/Detrended.nc',        &
                                 ylat, ylon, ylev, nforecasts_calibtracer,                   &
                                 CVT_HorizCors % TOMCAT_lats(1:ylat),                        &
                                 CVT_HorizCors % TOMCAT_longs(1:ylon),                       &
                                 tracer(1:ylat, 1:ylon, 1:ylev, 1:nforecasts_calibtracer),   &
                                 'CH4_detrend')


! ----------------------------------------------------------------------------------
! Part I of the calibration: compute the temporal mean
! ----------------------------------------------------------------------------------
PRINT *, 'Computing temporal mean'
DO lat = 1, ylat
  DO lon = 1, ylon
    DO lev = 1, ylev
      temporal_mean(lat,lon,lev) = SUM(tracer(lat,lon,lev,1:nforecasts_calibtracer)) / nforecasts_real
    END DO
  END DO
END DO

! ----------------------------------------------------------------------------------
! Part II of the calibration: remove the mean
! ----------------------------------------------------------------------------------
PRINT *, 'Removing mean'
DO lat = 1, ylat
  DO lon = 1, ylon
    DO lev = 1, ylev
      tracer1(lat,lon,lev,1:nforecasts_calibtracer) = tracer(lat,lon,lev,1:nforecasts_calibtracer) - temporal_mean(lat,lon,lev)
    END DO
  END DO
END DO

! ----------------------------------------------------------------------------------
! Part III of the calibration: compute the standard deviations
! Note these are used to normalise the data to find the vert cov matrix
! Then the standard deviations are overwritten with user's choice
! ----------------------------------------------------------------------------------
PRINT *, 'Computing standard deviations'
DO lat = 1, ylat
  DO lon = 1, ylon
    DO lev = 1, ylev
      CVT_std % std_tracer(lat,lon,lev) = SQRT( SUM(tracer1(lat,lon,lev,1:nforecasts_calibtracer) *    &
                                                    tracer1(lat,lon,lev,1:nforecasts_calibtracer)) /   &
                                          nforecasts_real )
    END DO
  END DO
END DO


! ----------------------------------------------------------------------------------
! Part IV of the calibration: normalize by the standard deviations
! ----------------------------------------------------------------------------------
PRINT *, 'Normalizing'
DO lat = 1, ylat
  DO lon = 1, ylon
    DO lev = 1, ylev
      tracer1(lat,lon,lev,1:nforecasts_calibtracer) = tracer1(lat,lon,lev,1:nforecasts_calibtracer) /  &
                                                      CVT_std % std_tracer(lat,lon,lev)
    END DO
  END DO
END DO

! ----------------------------------------------------------------------------------
! Overwrite the standard deviations with user's choice
! ----------------------------------------------------------------------------------
PRINT *, 'Overwriting standard deviations'
DO lat = 1, ylat
  DO lon = 1, ylon
    DO lev = 1, ylev
      CVT_std % std_tracer(lat,lon,lev) = 0.1
    END DO
  END DO
END DO

! ----------------------------------------------------------------------------------
! Part V of the calibration: find the global average vertical covariance matrix
! ----------------------------------------------------------------------------------
PRINT *, 'Find global average cov matrix'
DO lev1 = 1, ylev
  DO lev2 = lev1, ylev
    CVT_VertCors % vert_eigenvec_tracer(lev1,lev2) = SUM(tracer1(1:ylat,1:ylon,lev1,1:nforecasts_calibtracer) *  &
                                                         tracer1(1:ylat,1:ylon,lev2,1:nforecasts_calibtracer)) / &
                                                     (real(ylat) * real(ylon) * nforecasts_real)
    CVT_VertCors % glob_av_vert_cov(lev1,lev2)     = CVT_VertCors % vert_eigenvec_tracer(lev1,lev2)
  END DO
END DO

! Fill -in the lower triangular part of CVT_VertCors % glob_av_vert_cov (using symmetry)
DO lev1 = 1, ylev-1
  DO lev2 = lev1+1, ylev
    CVT_VertCors % glob_av_vert_cov(lev2,lev1) = CVT_VertCors % glob_av_vert_cov(lev1,lev2)
  END DO
END DO

! Compute the eigenvectors and eigenvalues
PRINT *, 'Find the eigen modes'
CALL DSYEV('V',                                       & ! Eigenvalues and vectors to be computed
           'U',                                       & ! Upper triangular of matrix specified
            ylev,                                     & ! Order of the matrix to be diagonalised
            CVT_VertCors % vert_eigenvec_tracer(1:ylev, 1:ylev), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
            ylev,                                     & ! Leading dimension of the matrix
            vert_evals(1:ylev),                       & ! Eigenvalues
            Workv(1:3*ylev),                          & ! Work array
            3*ylev,                                   & ! Size of work array
            err)
PRINT *, 'error flag after eigen solver for tracer ', err



PRINT *, 'Eigenvalues of global vertical error covariance matrix (before ridge regression): '
PRINT *, vert_evals(1:ylev)
! Do a ridge regression to remove negative and zero eigenvalues
CALL RidgeRegression (ylev,                             &
                      vert_evals(1:ylev),               &
                      0.01)

PRINT *, 'Eigenvalues of global vertical error covariance matrix (after ridge regression): '
PRINT *, vert_evals(1:ylev)



! ----------------------------------------------------------------------------------
! Part VI of the calibration: Compute local variances according to option
! ----------------------------------------------------------------------------------
PRINT *, 'Computing local variances'

SELECT CASE (CVT_VertCors % vert_covs)

CASE (1)
  ! The vertical variance spectrum is a function of lat and long
  CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1:ylon) = 0.0
  DO lat = 1, ylat
    DO lon = 1, ylon
      ! Project tracer ensemble onto this long/lat for each ensemble member
      DO forecastno = 1, nforecasts_calibtracer
        DO mode = 1, ylev
          psi = SUM(CVT_VertCors % vert_eigenvec_tracer(1:ylev, mode) * &
                    tracer1(lat, lon, 1:ylev, forecastno))
          CVT_VertCors % var_vspec_tracer(mode, lat, lon) = CVT_VertCors % var_vspec_tracer(mode, lat, lon) + &
                                                            psi * psi
        END DO
      END DO
    END DO
  END DO
  ! Compute the variance over the forecasts
  CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1:ylon) = CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1:ylon) / &
                                                 real(nforecasts_calibtracer)
  CVT_VertCors % std_vspec_tracer(1:ylev, 1:ylat, 1:ylon) = SQRT(CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1:ylon))

CASE (2)
  ! The vertical variance spectrum is a function of lat
  CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1) = 0.0
  DO lat = 1, ylat
    ! Project tracer ensemble onto this lat for each ensemble member and long
    DO lon = 1, ylon
      DO forecastno = 1, nforecasts_calibtracer
        DO mode = 1, ylev
          psi = SUM(CVT_VertCors % vert_eigenvec_tracer(1:ylev, mode) * tracer1(lat, lon, 1:ylev, forecastno))
          CVT_VertCors % var_vspec_tracer(mode, lat, 1) = CVT_VertCors % var_vspec_tracer(mode, lat, 1) + psi * psi
        END DO
      END DO
    END DO
  END DO
  ! Compute the variance over the forecasts and longitudes
  CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1) = CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1) / &
                                            real(nforecasts_calibtracer * ylon)
  CVT_VertCors % std_vspec_tracer(1:ylev, 1:ylat, 1) = SQRT(CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1))

  ! Copy this to each longitude
  DO lon = 2, ylon
    CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, lon) = CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1)
    CVT_VertCors % std_vspec_tracer(1:ylev, 1:ylat, lon) = CVT_VertCors % std_vspec_tracer(1:ylev, 1:ylat, 1)
  END DO

CASE (3)
  ! The variances are the eigenvalues
  ! Copy to each longitude and latitude
  DO lat = 1, ylat
    DO lon = 1, ylon
      CVT_VertCors % var_vspec_tracer(1:ylev, lat, lon) = vert_evals(1:ylev)
      CVT_VertCors % std_vspec_tracer(1:ylev, lat, lon) = SQRT(vert_evals(1:ylev))
    END DO
  END DO

END SELECT




! ----------------------------------------------------------------------------------
! Part VII of the calibration: Compute adjustments to ensure implied correlation
! ----------------------------------------------------------------------------------
PRINT *, 'Compute adjustment for correlation (if requested)'


IF (CVT_VertCors % ForceVCorr) THEN

  ! Compute the adjustment factors
  SELECT CASE (CVT_VertCors % vert_covs)

  CASE (1)
    ! The vertical variance spectrum is a function of lat and long
    DO lat = 1, ylat
      DO lon = 1, ylon
        DO lev = 1, ylev
          CVT_VertCors % vert_adjust_tracer(lev, lat, lon) = SQRT(SUM(CVT_VertCors % vert_eigenvec_tracer(lev,1:ylev) *  &
                                                                      CVT_VertCors % vert_eigenvec_tracer(lev,1:ylev) *  &
                                                                      CVT_VertCors % var_vspec_tracer(1:ylev, lat, lon)))
        END DO
      END DO
    END DO

  CASE (2)
    ! The vertical variance spectrum is a function of lat
    DO lat = 1, ylat
      DO lev = 1, ylev
        CVT_VertCors % vert_adjust_tracer(lev, lat, 1) = SQRT(SUM(CVT_VertCors % vert_eigenvec_tracer(lev,1:ylev) *  &
                                                                  CVT_VertCors % vert_eigenvec_tracer(lev,1:ylev) *  &
                                                                  CVT_VertCors % var_vspec_tracer(1:ylev, lat, 1)))
      END DO
      ! Copy to all longitudes
      DO lon = 2, ylon
        CVT_VertCors % vert_adjust_tracer(1:ylev, 1:ylat, lon) = CVT_VertCors % vert_adjust_tracer(1:ylev, 1:ylat, 1)
      END DO
    END DO

  CASE (3)
    ! The variances are the eigenvalues
    ! The adjustment should be unity, but still compute to check
    DO lev = 1, ylev
      CVT_VertCors % vert_adjust_tracer(lev, 1, 1) = SQRT(SUM(CVT_VertCors % vert_eigenvec_tracer(lev,1:ylev) *  &
                                                              CVT_VertCors % vert_eigenvec_tracer(lev,1:ylev) *  &
                                                              CVT_VertCors % var_vspec_tracer(1:ylev, 1, 1)))
    END DO
    ! Copy to all longitudes and latitudes
    DO lat = 1, ylat
      DO lon = 1, ylon
        CVT_VertCors % vert_adjust_tracer(1:ylev, lat, lon) = CVT_VertCors % vert_adjust_tracer(1:ylev, 1, 1)
      END DO
    END DO

  END SELECT

ELSE

  ! Set the adjustment factors to unity as no adjustment is to be done
  CVT_VertCors % vert_adjust_tracer(1:ylev, 1:ylat, 1:ylon) = 1.0

END IF


IF (PRODUCT(CVT_HorizCors % lengthscale_tracer(1:ylev)) == 0.0) THEN
  PRINT *, 'Horizontal lengthscales for tracer not set manually, so determining spectra from data'
  PRINT *, 'This may not be a reliable process, so please check the results afterwards'

  ! ----------------------------------------------------------------------------------
  ! Part VIII of the calibration: Perform inverse vertical and spectral transforms
  ! ----------------------------------------------------------------------------------
  PRINT *, 'Performing inverse vertical and spectral transforms'

  DO forecastno = 1, nforecasts_calibtracer

    PRINT *, 'Inverse transform on fc ', forecastno
    ! Inverse vertical transform
    CALL cvt_v_inv ( tracer_verttrans(1:ylat, 1:ylon, 1:ylev),         & ! out
                     tracer(1:ylat, 1:ylon, 1:ylev, forecastno) )        ! in

    ! Change horizontal grid (TOMCAT TO SHtools)
    CALL ChangeHorizGrid (ylev,                                        & ! No of levels
                          ylon,                                        & ! No of longitudes tomcat grid
                          ylat,                                        & ! No of latitudes tomcat grid
                          tracer_verttrans(1:ylat, 1:ylon, 1:ylev),    & ! Data on tomcat grid
                          CVT_HorizCors % TOMCAT_longs(1:ylon),        & ! Longs on tomcat grid
                          CVT_HorizCors % TOMCAT_lats(1:ylat),         & ! Lats on tomcat grid
                          2*L+1,                                       & ! No of longitudes shtools grid
                          L+1,                                         & ! No of latitudes shtools grid
                          tracer_verttrans_sh(1:L+1, 1:2*L+1, 1:ylev), & ! Data on shtools grid
                          CVT_HorizCors % SHtools_longs(1:2*L+1),      & ! Longitudes of shtools grid
                          CVT_HorizCors % SHtools_lats(1:L+1),         & ! Latitudes of shtools grid
                          CVT_HorizCors % TOMCAT_to_SH)                  ! Details of grid change

    ! Do spherical transform
    DO lev = 1, ylev
      CALL Spherical_inv (tracer_hortrans(1:2, 0:L, 0:L, lev, forecastno),&! OUT
                          tracer_verttrans_sh(1:L+1, 1:2*L+1, lev))        ! IN
    END DO
  END DO




  ! ----------------------------------------------------------------------------------
  ! Part IX of the tracer calibration: Find provisional horizontal spectrum
  ! ----------------------------------------------------------------------------------

  PRINT *, 'Finding horizontal spectrum for tracer'
  CVT_HorizCors % var_hspec_tracer(0:L, 1:ylev) = 0.0

  DO forecastno = 1, nforecasts_calibtracer
    PRINT *, 'Computing contribution to hspec from forecast ', forecastno
    DO lev = 1, ylev
      DO ll = 0, L
        !mm = 0 contribution (cosine contribution only)
        CVT_HorizCors % var_hspec_tracer(ll, lev) = CVT_HorizCors % var_hspec_tracer(ll, lev) +        &
                                                    tracer_hortrans(1,ll,0,lev,forecastno) *           &
                                                    tracer_hortrans(1,ll,0,lev,forecastno)
        IF (ll > 0) THEN
          DO mm = 1, ll
            DO phase = 1, 2
              CVT_HorizCors % var_hspec_tracer(ll, lev) = CVT_HorizCors % var_hspec_tracer(ll, lev) +  &
                                                          tracer_hortrans(phase,ll,mm,lev,forecastno) *&
                                                          tracer_hortrans(phase,ll,mm,lev,forecastno)
            END DO
          END DO
        END IF
      END DO
    END DO
  END DO

  ! Normalize
  PRINT *, 'Normalizing horizontal variance spectrum'
  DO ll = 0, L
    CVT_HorizCors % var_hspec_tracer(ll, 1:ylev) = CVT_HorizCors % var_hspec_tracer(ll, 1:ylev) /   &
                                                   (REAL(nforecasts_calibtracer) * REAL(2*ll + 1))
  END DO



  ! ----------------------------------------------------------------------------------
  ! Part X of the calibration: Compute adjustments to ensure implied correlation
  ! ----------------------------------------------------------------------------------
  PRINT *, 'Modifying horizontal variance spectrum to give correlation (tracer)'

  DO lev = 1, ylev
    lambda = 0.0
    DO ll = 0, L
      indexll = CVT_HorizCors % Plm_index(ll, 0)
      lambda = lambda + CVT_HorizCors % var_hspec_tracer(ll, lev) * SQRT(REAL(2*ll + 1)) *  &
                        CVT_HorizCors % assocLegPoly(1, indexll)
    END DO
    lambda = 1.0 / lambda
    CVT_HorizCors % var_hspec_tracer(0:L, lev) = CVT_HorizCors % var_hspec_tracer(0:L, lev) * lambda
  END DO



  ! ----------------------------------------------------------------------------------
  ! Part XI of the calibration: Compute the provisional implied covariances at each level
  ! ----------------------------------------------------------------------------------

  ! Compute separation function (metres)
  DO sep = 1, L+1
    ! What is the separation function (dictated by the Gaussian grid)?
    dist(sep) = CVT_HorizCors % GaussianCoLats(sep) * deg2rad * Re
  END DO

  ! Compute the implied covariances at each of these separations (separately for each level)
  DO lev = 1, ylev
    DO sep = 1, L+1
      total = 0.0
      DO ll = 0, L
        indexll = CVT_HorizCors % Plm_index(ll, 0)
        total   = total + CVT_HorizCors % var_hspec_tracer(ll, lev) * SQRT(REAL(2*ll+1)) * &
                          CVT_HorizCors % assocLegPoly(sep, indexll)
      END DO
      tracer_cov_impl(sep, lev) = total
    END DO
  END DO


  ! ----------------------------------------------------------------------------------
  ! Part XII of the calibration: Fit lengthscales to these implied correlations
  ! ----------------------------------------------------------------------------------

  dL = 1.0   ! To compute finite difference derivatives
  DO lev = 1, ylev
    PRINT *, '******************** level', lev
    ! First guess - scale 1000km
    CVT_HorizCors % lengthscale_tracer(lev) = 400000.0 ! in m
    DO it = 1, 10
      num = 0.0
      den = 0.0
      DO sep = 1, (L+1) / 6   ! Only do local fit (over only a subset of points)
        corfn  = CorrelationFn(CVT_HorizCors % GaussianCoLats(sep-1) * deg2rad * Re,  &
                               CVT_HorizCors % horiztracer_cor_tpe,                   &
                               CVT_HorizCors % lengthscale_tracer(lev))               ! in metres
        derivL = ( CorrelationFn(CVT_HorizCors % GaussianCoLats(sep-1) * deg2rad * Re,&
                                 CVT_HorizCors % horiztracer_cor_tpe,                 &
                                 CVT_HorizCors % lengthscale_tracer(lev) + dL)        & ! in metres
                   - corfn ) / dL
        num = num + (corfn - tracer_cov_impl(sep, lev)) * derivL
        den = den + derivL * derivL
      END DO
      deltaL = -1.0 * num / den
      CVT_HorizCors % lengthscale_tracer(lev) = CVT_HorizCors % lengthscale_tracer(lev) + deltaL
      PRINT *, 'It', it, CVT_HorizCors % lengthscale_tracer(lev)
    END DO
  END DO

  ! Output the correlation functions (implied and fitted)
  DO lev = 1, ylev
    ! Create a filename
    WRITE (fname, '(A,I0.3,A)') 'Tracer_horiz_cor_lev', lev, '.dat'
    OPEN (12, file=fname)
    DO sep = 1, L+1
      WRITE (12,*) sep, CVT_HorizCors % GaussianCoLats(sep) * deg2rad * Re,  &
                   tracer_cov_impl(sep, lev),                                & ! The target corr fn
                   CorrelationFn(CVT_HorizCors % GaussianCoLats(sep-1) * deg2rad * Re,  &
                               CVT_HorizCors % horiztracer_cor_tpe,                     &
                               CVT_HorizCors % lengthscale_tracer(lev))        ! The fitted function
    END DO
    CLOSE (12)
  END DO

ELSE

  PRINT *, 'Horizontal lengthscales for tracer set manually (earlier in Calibration code).'

END IF



! ----------------------------------------------------------------------------------
! Part XIII of the calibration: Derive the horizontal spectra for these lengthscales
! ----------------------------------------------------------------------------------
DO lev = 1, ylev

  ! Set the correlation functions given the lengthscales
  DO sep = 1, L+1
    CVT_HorizCors % horizcors_tracer(sep,lev) = CorrelationFn(CVT_HorizCors % GaussianCoLats(sep-1) * deg2rad * Re, &
                                                              CVT_HorizCors % horiztracer_cor_tpe,                  &
                                                              CVT_HorizCors % lengthscale_tracer(lev)) ! in metres
  END DO

  ! Derive the horizontal variance spectrum for this level
  PRINT *, 'Deriving horizontal spectrum for tracer'
  CVT_HorizCors % var_hspec_tracer(0:L, lev) = 0.0
  DO ll = 0, L
    ! Index for this total wn in CVT_assocLegPoly array
    indexll = CVT_HorizCors % Plm_index(ll, 0)
    DO lat = 1, L+1
      CVT_HorizCors % var_hspec_tracer(ll,lev) = CVT_HorizCors % var_hspec_tracer(ll, lev) +   &
                                                 CVT_HorizCors % horizcors_tracer(lat, lev) *  &
                                                 CVT_HorizCors % assocLegPoly(lat, indexll) *  &
                                                 CVT_HorizCors % GaussianWts(lat)
    END DO
    CVT_HorizCors % var_hspec_tracer(ll,lev) = CVT_HorizCors % var_hspec_tracer(ll,lev) /      &
                                               ( 2.0 * SQRT(REAL(2*ll+1)) )
  END DO

  ! Do a ridge regression to eliminate zero and negative values
  CALL RidgeRegression (L+1, CVT_HorizCors % var_hspec_tracer(0:L,lev), 0.01)

  ! Modify to ensure correlation
  PRINT *, 'Modifying horizontal variance spectrum to give correlation (flux)'
  lambda = 0.0
  DO ll = 0, L
    indexll = CVT_HorizCors % Plm_index(ll, 0)
    lambda  = lambda + CVT_HorizCors % var_hspec_tracer(ll, lev) * SQRT(REAL(2*ll + 1)) *      &
                       CVT_HorizCors % assocLegPoly(1, indexll)
  END DO
  lambda = 1.0 / lambda
  CVT_HorizCors % var_hspec_tracer(0:L, lev) = CVT_HorizCors % var_hspec_tracer(0:L, lev) * lambda
  PRINT *, 'lev = ', lev, '    Lambda = ', lambda

END DO


! Compute standard deviations
PRINT *, 'Computing horizontal standard deviation spectrum'
CVT_HorizCors % std_hspec_tracer(0:L, 1:ylev) = SQRT(CVT_HorizCors % var_hspec_tracer(0:L, 1:ylev))



PRINT *, 'data around the equator ...'
PRINT *, CVT_HorizCors % SHtools_lats(GLatEq-1), CVT_HorizCors % GaussianWts(GLatEq-1)
PRINT *, CVT_HorizCors % SHtools_lats(GLatEq),   CVT_HorizCors % GaussianWts(GLatEq)
PRINT *, CVT_HorizCors % SHtools_lats(GLatEq+1), CVT_HorizCors % GaussianWts(GLatEq+1)





! ==================================================================================
! ----------------------------------------------------------------------------------
! Deal with the surface flux covariance model
! ----------------------------------------------------------------------------------

! Repeat for each field
DO field = 1, nfields

  ! ----------------------------------------------------------------------------------
  ! Part I of the flux calibration: specify the standard deviations
  ! ----------------------------------------------------------------------------------
  ! ##### USER INPUT #####
  PRINT *, 'Setting variance of flux field'
  CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, field) = 0.00001
  ! ######################


  ! ----------------------------------------------------------------------------------
  ! Part II of the flux calibration: specify the temporal correlation function
  !                                  and its spectrum
  ! ----------------------------------------------------------------------------------

  ! ##### USER INPUT #####
  CVT_TemporalCors % temporal_covs(field) = 1
  ! 0 = No temporal covariances for flux
  ! 1 = Temporal covariances, spatial scale a fn of timescale
  ! 2 = Temporal covariances, spatial scale a fn of time
  ! ######################

  IF (CVT_TemporalCors % temporal_covs(field) > 0) THEN

    ! For these options need to specify a correlation shape and timescale for the temporal flux correlations

    ! ##### USER INPUT #####
    PRINT *, 'Setting temporal correlation function for flux, field', field
    CVT_TemporalCors % temporal_cor_tpe(field) = 3
    ! 1 = Lorentzian
    ! 2 = Gaussian
    ! 3 = SOAR
    ! 4 = exponential

    CVT_TemporalCors % timescale_flux(field)   = 6.0
    ! ######################

    ! Set the correlation function as a function of separation (may also be a function of field)
    DO time1 = 1, nmonth
      CVT_TemporalCors % temporalcors_flux(time1,field) = CorrelationFn(REAL(time1-1),                             &
                                                                        CVT_TemporalCors % temporal_cor_tpe(field),&
                                                                        CVT_TemporalCors % timescale_flux(field)) ! in months
    END DO

    ! Use this to generate the upper triangular part of a correlation matrix
    PRINT *, 'Generating temporal correlation function'
    DO time1 = 1, nmonth
      DO time2 = time1, nmonth
        tdiff = ABS(time1 - time2)
        CVT_TemporalCors % temp_eigenvec_flux(time1,time2,field) = CVT_TemporalCors % temporalcors_flux(tdiff+1,field)
        CVT_TemporalCors % temp_cor_matrix(time1,time2,field)    = CVT_TemporalCors % temp_eigenvec_flux(time1,time2,field)
      END DO
    END DO

    ! Fill-in the lower part of the correlation matrix (use symmetry)
    DO time1 = 1, nmonth-1
      DO time2 = time1+1, nmonth
        CVT_TemporalCors % temp_cor_matrix(time2,time1,field) = CVT_TemporalCors % temp_cor_matrix(time1,time2,field)
      END DO
    END DO

    ! Compute the eigenvectors and eigenvalues
    PRINT *, 'Finding eigenvalues and eigenvectors'
    CALL DSYEV('V',                                       & ! Eigenvalues and vectors to be computed
               'U',                                       & ! Upper triangular of matrix specified
                nmonth,                                   & ! Order of the matrix to be diagonalised
                CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, 1:nmonth, field), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
                nmonth,                                   & ! Leading dimension of the matrix
                CVT_TemporalCors % var_tspec_flux(1:nmonth, field),               & ! Eigenvalues
                Workt(1:3*nmonth),                        & ! Work array
                3*nmonth,                                 & ! Size of work array
                err)
    PRINT *, 'error flag after eigen solver for flux (field', field, ') ', err

    ! Compute the standard deviations
    PRINT *, 'Computing standard deviations'
    CVT_TemporalCors % std_tspec_flux(1:nmonth, field) = SQRT(CVT_TemporalCors % var_tspec_flux(1:nmonth, field))


  ELSE

    PRINT *, 'No temporal or spatial correlations for flux field'
    CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, 1:nmonth, field) = 0.0
    CVT_TemporalCors % var_tspec_flux(1:nmonth, field)               = 0.0
    CVT_TemporalCors % std_tspec_flux(1:nmonth, field)               = 0.0
    CVT_TemporalCors % temporal_cor_tpe(field)                       = 0
    CVT_TemporalCors % timescale_flux(field)                         = 0.0
    CVT_TemporalCors % temporalcors_flux(1:nmonth, field)            = 0.0

  END IF


  ! ----------------------------------------------------------------------------------
  ! Part III of the flux calibration: specify the spatial correlation function
  !                                   and its spectrum
  ! ----------------------------------------------------------------------------------
  ! Note for CVT_TemporalCors % temporal_covs = 1 the lengthscales are prescribed as a fn of timescale
  !      for CVT_TemporalCors % temporal_covs = 2 the lengthscales are prescribed as a fn of time
  ! ##### USER INPUT #####
  CVT_HorizCors % horizflux_cor_tpe(field) = 3
  ! 1 = Lorentzian
  ! 2 = Gaussian
  ! 3 = SOAR
  ! 4 = exponential
  ! ######################

  DO time1 = 1, nmonth
    ! time1 is either time or timescale (depending on temporal_covs(field) value)
    PRINT *, 'Dealing with flux lengthscales for field', field, ' and timescale/time', time1
    PRINT *, 'Setting horizontal correlation fn'

    ! ##### USER INPUT #####
    ! Lengthscale in metres
    CVT_HorizCors % lengthscale_flux(time1,field) = 400000.0
    ! ######################


    !PRINT *, 'Flux correlation type: ', CVT_HorizCors % horizflux_cor_tpe(field)
    !PRINT *, 'Lengthscale (m):       ', CVT_HorizCors % lengthscale_flux(time1,field)
    DO sep = 1, L+1
      CVT_HorizCors % horizcors_flux(sep,time1,field) = CorrelationFn(CVT_HorizCors % GaussianCoLats(sep-1) * deg2rad * Re, &
                                                                      CVT_HorizCors % horizflux_cor_tpe(field),             &
                                                                      CVT_HorizCors % lengthscale_flux(time1,field)) ! in metres
      !PRINT *, time1, sep, CVT_HorizCors % GaussianCoLats(sep) * deg2rad * Re, CVT_HorizCors % horizcors_flux(sep,time1,field)
    END DO


    ! Derive the horizontal variance spectrum for this time/timescale and field
    PRINT *, 'Deriving horizontal spectrum'
    CVT_HorizCors % var_hspec_flux(0:L, time1, field) = 0.0
    DO ll = 0, L
      ! Index for this total wn in CVT_assocLegPoly array
      indexll = CVT_HorizCors % Plm_index(ll, 0)
      DO lat = 1, L+1
        CVT_HorizCors % var_hspec_flux(ll,time1,field) = CVT_HorizCors % var_hspec_flux(ll, time1, field) +        &
                                                         CVT_HorizCors % horizcors_flux(lat, time1, field) *&
                                                         CVT_HorizCors % assocLegPoly(lat, indexll) *              &
                                                         CVT_HorizCors % GaussianWts(lat)
      END DO
      CVT_HorizCors % var_hspec_flux(ll,time1,field) = CVT_HorizCors % var_hspec_flux(ll,time1,field) /     &
                                                       ( 2.0 * SQRT(REAL(2*ll+1)) )
    END DO

    ! Do a ridge regression to eliminate zero and negative values
    CALL RidgeRegression (L+1, CVT_HorizCors % var_hspec_flux(0:L,time1,field), 0.01)

    ! Modify to ensure correlation
    PRINT *, 'Modifying horizontal variance spectrum to give correlation (flux)'
    lambda = 0.0
    DO ll = 0, L
      indexll = CVT_HorizCors % Plm_index(ll, 0)
      lambda  = lambda + CVT_HorizCors % var_hspec_flux(ll, time1, field) * SQRT(REAL(2*ll + 1)) *  &
                         CVT_HorizCors % assocLegPoly(1, indexll)
    END DO
    lambda = 1.0 / lambda
    CVT_HorizCors % var_hspec_flux(0:L, time1, field) = CVT_HorizCors % var_hspec_flux(0:L, time1, field) * lambda
    PRINT *, 'time = ', time1, '    Lambda = ', lambda

  END DO
  ! Compute the standard deviation of the horizontal spectrum
  PRINT *, 'Finding standard deviations'
  CVT_HorizCors % std_hspec_flux(0:L,1:nmonth,field) = SQRT(CVT_HorizCors % var_hspec_flux(0:L,1:nmonth,field))


END DO

! Compute global mean altitudes
DO lev = 1, ylev
  CVT_VertCors % alts(lev) = SUM(alt(1:ylat, 1:ylon, lev)) / (REAL(ylat) * REAL(ylon))
END DO

! Set the month indices (here for EnviFlux, these are not months, but in units of sDt
sDt = 2592000.0   !30 days
DO month = 1, nmonth
  CVT_TemporalCors % months(month) = REAL(month-1) * sDt
END DO


! ----------------------------------------------------------------------------------
! Output the CVT data
! ----------------------------------------------------------------------------------

! Output the data
CALL cvt_matrices_output ('CVT_calib.nc',                                                       &
! ***** Longitudes, latitudes, heights, and months
                          CVT_HorizCors % TOMCAT_longs(1:ylon),                                 &
                          CVT_HorizCors % TOMCAT_lats(1:ylat),                                  &
                          CVT_VertCors % alts(1:ylev),                                          &
                          CVT_HorizCors % SHtools_longs(1:2*L+1),                               &
                          CVT_HorizCors % SHtools_lats(1:L+1),                                  &
                          CVT_TemporalCors % months(1:nmonth),                                  &
                          CVT_HorizCors % GaussianCoLats(1:L+1),                                &
! ***** Spatial standard deviations
                          CVT_std % FracOfPrior_tracer,                                         &
                          CVT_std % std_tracer(1:ylat, 1:ylon, 1:ylev),                         &
                          CVT_std % FracOfPrior_flux,                                           &
                          CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),              &
! ***** Horizontal transform for tracer
                          CVT_HorizCors % horiztracer_cor_tpe,                                  &
                          CVT_HorizCors % lengthscale_tracer(1:ylev),                           &
                          CVT_HorizCors % horizcors_tracer(1:L+1, 1:ylev),                      &
                          CVT_HorizCors % var_hspec_tracer(0:L, 1:ylev),                        &
                          CVT_HorizCors % std_hspec_tracer(0:L, 1:ylev),                        &
! ***** Horizontal transform for flux
                          CVT_HorizCors % horizflux_cor_tpe(1:nfields),                         &
                          CVT_HorizCors % lengthscale_flux(1:nmonth, 1:nfields),                &
                          CVT_HorizCors % horizcors_flux(1:L+1, 1:nmonth, 1:nfields),           &
                          CVT_HorizCors % var_hspec_flux(0:L, 1:nmonth, 1:nfields),             &
                          CVT_HorizCors % std_hspec_flux(0:L, 1:nmonth, 1:nfields),             &
! ***** Vertical transform for tracer
                          CVT_VertCors % ForceVCorr,                                            &
                          CVT_VertCors % vert_covs,                                             &
                          CVT_VertCors % glob_av_vert_cov(1:ylev, 1:ylev),                      &
                          CVT_VertCors % vert_eigenvec_tracer(1:ylev, 1:ylev),                  &
                          CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1:ylon),              &
                          CVT_VertCors % std_vspec_tracer(1:ylev, 1:ylat, 1:ylon),              &
                          CVT_VertCors % vert_adjust_tracer(1:ylev, 1:ylat, 1:ylon),            &
! ***** Temporal transform for flux
                          CVT_TemporalCors % temporal_covs(1:nfields),                          &
                          CVT_TemporalCors % temporal_cor_tpe(1:nfields),                       &
                          CVT_TemporalCors % timescale_flux(1:nfields),                         &
                          CVT_TemporalCors % temporalcors_flux(1:nmonth, 1:nfields),            &
                          CVT_TemporalCors % temp_cor_matrix(1:nmonth, 1:nmonth, 1:nfields),    &
                          CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, 1:nmonth, 1:nfields), &
                          CVT_TemporalCors % var_tspec_flux(1:nmonth, 1:nfields),               &
                          CVT_TemporalCors % std_tspec_flux(1:nmonth, 1:nfields) )


! ----------------------------------------------------------------------------------
! Deallocate structures
! ----------------------------------------------------------------------------------

CALL UnInitializeHorizTrans ()
CALL DeallocateCVT (.TRUE.)


END PROGRAM Calibrate
