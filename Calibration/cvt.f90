! This file contains the routines to perform the control variable transforms
! (forward, inverse, adjoint, and inverse adjoint)
! For use with the non-balanced horizontal transform


!=================================================================================================
SUBROUTINE InitializeHorizTrans ()


! Description:
! Allocates and initialises data for the horizontal transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      11/2019  Initial version (R.N. Bannister)
!           01/2020  Use structures and read shtools stuff from file
!
! Code Description:
!   Language:           Fortran 90.


USE main_data,    ONLY: ylon, ylat

USE cparam,       ONLY: pi, rad2deg

USE onedvar_data, ONLY: L,                     &
                        ALPsize,               &
                        CVT_HorizCors,         &
                        SHtools_Normalization, &
                        SHtools_ApplyCSphase,  &
                        SHtools_CRnorm,        &
                        fft_worklen,           &
                        ALPdatafile

IMPLICIT NONE

! Local variables:
INTEGER           :: i, j, colat, mm, ll, index, err
INTEGER           :: ylat_file, ylon_file, ylev_file, L_file, ALPsize_file
CHARACTER(len=64) :: blankline


! Allocate the required arrays for associated Legendre polynomials
ALLOCATE (CVT_HorizCors % assocLegPoly(1:L+1, 1:ALPsize))  ! Associated Legendre polynomials (colat, order)
ALLOCATE (CVT_HorizCors % GaussianCosCoLats(0:L+1))        ! Cosines of the Gaussian co-latitudes
ALLOCATE (CVT_HorizCors % GaussianCoLats(0:L+1))           ! Gaussian co-latitudes
ALLOCATE (CVT_HorizCors % GaussianWts(1:L+1))              ! Gaussian weights
ALLOCATE (CVT_HorizCors % Plm_index(0:L,0:L))              ! For a given l, m this will indicate the index of CVT_assocLegPoly
ALLOCATE (CVT_HorizCors % SHtools_longs(1:2*L+1))
ALLOCATE (CVT_HorizCors % SHtools_lats(0:L+1))

! Allocate the required arrays for FFT
ALLOCATE (CVT_HorizCors % fft_wsave(1:fft_worklen))
ALLOCATE (CVT_HorizCors % fft_work(1:2*L+1))


PRINT *, 'Initialising horizonal transform'

!Read-in associated Legendre polynomials data
OPEN (13, FILE=ALPdatafile)
READ (13,'(A)')     blankline
READ (13,'(A)')     blankline
READ (13,'(A)')     blankline
READ (13,'(A24I4)') blankline, ylat_file
READ (13,'(A24I4)') blankline, L_file
READ (13,'(A)')     blankline
READ (13,'(A)')     blankline
READ (13,'(A)')     blankline
READ (13,'(A24I4)') blankline, ALPsize_file
READ (13,'(A)')     blankline
READ (13,'(A)')     blankline

! Check that this file matches the parameters in this compilation
IF ((ylat_file /= ylat) .OR. (L_file /= L) .OR. (ALPsize_file /= ALPsize)) THEN
  PRINT *, 'Error - one or more of the programs parameters does not match those in cvt_ass_legendre_poly.dat'
  PRINT *, 'Please change the following in the programs module and recompile, or use a different file'
  PRINT *, 'ylat:   ', ylat,    ' -->', ylat_file
  PRINT *, 'L:      ', L   ,    ' -->', L_file
  PRINT *, 'ALPsize:', ALPsize, ' -->', ALPsize_file
  CLOSE (13)
  STOP
END IF

DO colat = 1, L+1
  READ (13,*) CVT_HorizCors % GaussianCosCoLats(colat)
END DO
READ (13,'(A)') blankline
READ (13,'(A)') blankline
DO colat = 1, L+1
  READ (13,*) CVT_HorizCors % GaussianWts(colat)
END DO
READ (13,'(A)') blankline
READ (13,'(A)') blankline
READ (13, *)    CVT_HorizCors % Plm_index(0:L, 0:L)
READ (13,'(A)') blankline
READ (13,'(A)') blankline
DO mm = 0, L
  DO ll = mm, L
    index = CVT_HorizCors % Plm_index(ll,mm)
    READ (13,'(A)') blankline
    DO colat = 1, L+1
      READ (13,*) CVT_HorizCors % assocLegPoly(colat,index)
    END DO
  END DO
END DO
READ (13,'(A)') blankline
READ (13,'(A)') blankline
CLOSE (13)

! Compute the Gaussian co-latitudes and latitudes (SHtools latitude positions)
DO colat = 0, L+1
  CVT_HorizCors % GaussianCoLats(colat) = rad2deg * ACOS(CVT_HorizCors % GaussianCosCoLats(colat))
  CVT_HorizCors % SHtools_lats(colat)   = 90.0 - CVT_HorizCors % GaussianCoLats(colat)
END DO

! Set the SHtools longitude positions
DO i = 1, 2*L+1
  CVT_HorizCors % SHtools_longs(i) = REAL(i-1) * 360.0 / REAL(2*L+1)
END DO

! Initialize FFT structures
CALL rfft1i (2*L+1, CVT_HorizCors % fft_wsave(1:fft_worklen), fft_worklen, err)
PRINT *, 'Exit status from rfft1i: ', err


END SUBROUTINE InitializeHorizTrans





!=================================================================================================
SUBROUTINE UnInitializeHorizTrans ()


! Description:
! Un-initialises data for the horizontal transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      11/2019  Initial version (R.N. Bannister)
!           01/2020  Use structures
!
! Code Description:
!   Language:           Fortran 90.


USE onedvar_data, ONLY: CVT_HorizCors


IMPLICIT NONE

! Deallocate the required arrays for associated Legendre polynomials
DEALLOCATE (CVT_HorizCors % assocLegPoly)
DEALLOCATE (CVT_HorizCors % GaussianCosCoLats)
DEALLOCATE (CVT_HorizCors % GaussianCoLats)
DEALLOCATE (CVT_HorizCors % GaussianWts)
DEALLOCATE (CVT_HorizCors % Plm_index)
DEALLOCATE (CVT_HorizCors % SHtools_longs)
DEALLOCATE (CVT_HorizCors % SHtools_lats)

! Deallocate the required arrays for FFT
DEALLOCATE (CVT_HorizCors % fft_wsave)
DEALLOCATE (CVT_HorizCors % fft_work)


IF (ALLOCATED(CVT_HorizCors % SH_to_TOMCAT % grid_change_long_lower))    &
  DEALLOCATE(CVT_HorizCors % SH_to_TOMCAT % grid_change_long_lower)
IF (ALLOCATED(CVT_HorizCors % SH_to_TOMCAT % grid_change_long_lower_wt)) &
  DEALLOCATE(CVT_HorizCors % SH_to_TOMCAT % grid_change_long_lower_wt)
IF (ALLOCATED(CVT_HorizCors % SH_to_TOMCAT % grid_change_long_upper_wt)) &
  DEALLOCATE(CVT_HorizCors % SH_to_TOMCAT % grid_change_long_upper_wt)
IF (ALLOCATED(CVT_HorizCors % SH_to_TOMCAT % grid_change_lat_lower))     &
  DEALLOCATE(CVT_HorizCors % SH_to_TOMCAT % grid_change_lat_lower)
IF (ALLOCATED(CVT_HorizCors % SH_to_TOMCAT % grid_change_lat_lower_wt))  &
  DEALLOCATE(CVT_HorizCors % SH_to_TOMCAT % grid_change_lat_lower_wt)
IF (ALLOCATED(CVT_HorizCors % SH_to_TOMCAT % grid_change_lat_upper_wt))  &
  DEALLOCATE(CVT_HorizCors % SH_to_TOMCAT % grid_change_lat_upper_wt)

IF (ALLOCATED(CVT_HorizCors % TOMCAT_to_SH % grid_change_long_lower))    &
  DEALLOCATE(CVT_HorizCors % TOMCAT_to_SH % grid_change_long_lower)
IF (ALLOCATED(CVT_HorizCors % TOMCAT_to_SH % grid_change_long_lower_wt)) &
  DEALLOCATE(CVT_HorizCors % TOMCAT_to_SH % grid_change_long_lower_wt)
IF (ALLOCATED(CVT_HorizCors % TOMCAT_to_SH % grid_change_long_upper_wt)) &
  DEALLOCATE(CVT_HorizCors % TOMCAT_to_SH % grid_change_long_upper_wt)
IF (ALLOCATED(CVT_HorizCors % TOMCAT_to_SH % grid_change_lat_lower))     &
  DEALLOCATE(CVT_HorizCors % TOMCAT_to_SH % grid_change_lat_lower)
IF (ALLOCATED(CVT_HorizCors % TOMCAT_to_SH % grid_change_lat_lower_wt))  &
  DEALLOCATE(CVT_HorizCors % TOMCAT_to_SH % grid_change_lat_lower_wt)
IF (ALLOCATED(CVT_HorizCors % TOMCAT_to_SH % grid_change_lat_upper_wt))  &
  DEALLOCATE(CVT_HorizCors % TOMCAT_to_SH % grid_change_lat_upper_wt)

END SUBROUTINE UnInitializeHorizTrans




!=================================================================================================
SUBROUTINE AllocateCVT (ForCalibration)
! Description:
! Allocates various arrays for CVT
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      12/2019  Initial version (R.N. Bannister)
!           01/2020  Use structures
!           07/2021  Include arrays used for calibration only
!
! Code Description:
!   Language:           Fortran 90.


USE main_data,    ONLY: ylon, ylat, ylev,           &
                        nmonth,                     &
                        nfields

USE onedvar_data, ONLY: L,                          &
                        CVT_std,                    &
                        CVT_HorizCors,              &
                        CVT_VertCors,               &
                        CVT_TemporalCors

IMPLICIT NONE

! Arguments:
LOGICAL, INTENT(IN) :: ForCalibration

ALLOCATE (CVT_std % std_tracer(1:ylat, 1:ylon, 1:ylev))
ALLOCATE (CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ALLOCATE (CVT_HorizCors % TOMCAT_longs(0:ylon+1))
ALLOCATE (CVT_HorizCors % TOMCAT_lats(0:ylat+1))
ALLOCATE (CVT_HorizCors % var_hspec_tracer(0:L, 1:ylev))
ALLOCATE (CVT_HorizCors % std_hspec_tracer(0:L, 1:ylev))
ALLOCATE (CVT_HorizCors % var_hspec_flux(0:L, 1:nmonth, 1:nfields))
ALLOCATE (CVT_HorizCors % std_hspec_flux(0:L, 1:nmonth, 1:nfields))
ALLOCATE (CVT_VertCors % alts(1:ylev))
ALLOCATE (CVT_VertCors % vert_eigenvec_tracer(1:ylev,1:ylev))
ALLOCATE (CVT_VertCors % var_vspec_tracer(1:ylev, 1:ylat, 1:ylon))
ALLOCATE (CVT_VertCors % std_vspec_tracer(1:ylev, 1:ylat, 1:ylon))
ALLOCATE (CVT_VertCors % vert_adjust_tracer(1:ylev, 1:ylat, 1:ylon))
ALLOCATE (CVT_TemporalCors % months(1:nmonth))
ALLOCATE (CVT_TemporalCors % temporal_covs(1:nfields))
ALLOCATE (CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, 1:nmonth, 1:nfields))
ALLOCATE (CVT_TemporalCors % var_tspec_flux(1:nmonth, 1:nfields))
ALLOCATE (CVT_TemporalCors % std_tspec_flux(1:nmonth, 1:nfields))

IF (ForCalibration) THEN
  ! Allocate extra variables for calibration
  ALLOCATE (CVT_VertCors % glob_av_vert_cov(1:ylev,1:ylev))
  ALLOCATE (CVT_HorizCors % lengthscale_tracer(1:ylev))
  ALLOCATE (CVT_HorizCors % horizcors_tracer(1:L+1, 1:ylev))
  ALLOCATE (CVT_HorizCors % horizflux_cor_tpe(1:nfields))
  ALLOCATE (CVT_HorizCors % lengthscale_flux(1:nmonth, 1:nfields))
  ALLOCATE (CVT_HorizCors % horizcors_flux(1:L+1, 1:nmonth, 1:nfields))
  ALLOCATE (CVT_TemporalCors % timescale_flux(1:nfields))
  ALLOCATE (CVT_TemporalCors % temporal_cor_tpe(1:nfields))
  ALLOCATE (CVT_TemporalCors % temporalcors_flux(1:nmonth, 1:nfields))
  ALLOCATE (CVT_TemporalCors % temp_cor_matrix(1:nmonth, 1:nmonth, 1:nfields))
END IF

END SUBROUTINE AllocateCVT



!=================================================================================================
SUBROUTINE DeallocateCVT (ForCalibration)
! Description:
! Deallocates various arrays for CVT
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      12/2019  Initial version (R.N. Bannister)
!           01/2020  Use structures
!
! Code Description:
!   Language:           Fortran 90.


USE onedvar_data, ONLY: L,                          &
                        CVT_std,                    &
                        CVT_HorizCors,              &
                        CVT_VertCors,               &
                        CVT_TemporalCors

IMPLICIT NONE

! Arguments:
LOGICAL, INTENT(IN) :: ForCalibration

DEALLOCATE (CVT_std % std_tracer)
DEALLOCATE (CVT_std % std_flux)
DEALLOCATE (CVT_HorizCors % TOMCAT_longs)
DEALLOCATE (CVT_HorizCors % TOMCAT_lats)
DEALLOCATE (CVT_HorizCors % var_hspec_tracer)
DEALLOCATE (CVT_HorizCors % std_hspec_tracer)
DEALLOCATE (CVT_HorizCors % var_hspec_flux)
DEALLOCATE (CVT_HorizCors % std_hspec_fluX)
DEALLOCATE (CVT_VertCors % alts)
DEALLOCATE (CVT_VertCors % vert_eigenvec_tracer)
DEALLOCATE (CVT_VertCors % var_vspec_tracer)
DEALLOCATE (CVT_VertCors % std_vspec_tracer)
DEALLOCATE (CVT_VertCors % vert_adjust_tracer)
DEALLOCATE (CVT_TemporalCors % months)
DEALLOCATE (CVT_TemporalCors % temporal_covs)
DEALLOCATE (CVT_TemporalCors % temp_eigenvec_flux)
DEALLOCATE (CVT_TemporalCors % var_tspec_flux)
DEALLOCATE (CVT_TemporalCors % std_tspec_flux)
IF (ForCalibration) THEN
  DEALLOCATE (CVT_VertCors % glob_av_vert_cov)
  DEALLOCATE (CVT_HorizCors % lengthscale_tracer)
  DEALLOCATE (CVT_HorizCors % horizcors_tracer)
  DEALLOCATE (CVT_HorizCors % horizflux_cor_tpe)
  DEALLOCATE (CVT_HorizCors % lengthscale_flux)
  DEALLOCATE (CVT_HorizCors % horizcors_flux)
  DEALLOCATE (CVT_TemporalCors % timescale_flux)
  DEALLOCATE (CVT_TemporalCors % temporal_cor_tpe)
  DEALLOCATE (CVT_TemporalCors % temporalcors_flux)
  DEALLOCATE (CVT_TemporalCors % temp_cor_matrix)
END IF

END SUBROUTINE DeallocateCVT



!=================================================================================================
SUBROUTINE cvt_total (chi,           & ! IN  Control variable
                      zy,            & ! OUT Physical space increment in INVICAT format
                      do_flux,       & ! IN  true/false to deal with 2D+time flux field
                      do_tracer)       ! IN  true/false to deal with 3D intial tracer field

! Description:
! Performs the complete control variable transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Horizontal, followed by vertical (for tracer) or temporal (for flux)
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields, tidy things
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, dim_1d, nfields, gs

USE onedvar_data, ONLY: n1d

! Note : dim_1d  is the size of the physical vector comprising the tracer and flux field
!        n1d     is the size of the control vector space comprising the tracer and flux field
!        nfields is the number of types of flux field (e.g. nfield=2 for anthropogenic and natural fluxes.
!        gs      is the number of points on horizontal grid


IMPLICIT NONE

INCLUDE "cvt.interface"

! Arguments:
REAL,    INTENT(IN)  :: chi(1:n1d)
REAL,    INTENT(OUT) :: zy(1:dim_1d)
LOGICAL, INTENT(IN)  :: do_flux
LOGICAL, INTENT(IN)  :: do_tracer

! Local variables:
REAL                 :: chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: chi_vh_tracer(1:ylat, 1:ylon, 1:ylev)
REAL                 :: state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: state_tracer(1:ylat, 1:ylon, 1:ylev)


! 1. Do the horizontal transform
! ------------------------------
CALL cvt_h (chi(1:n1d),                                         &
            chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),   &
            chi_vh_tracer(1:ylat, 1:ylon, 1:ylev),              &
            do_flux,                                            &
            do_tracer)

! 2. Do the vertical transform for the tracer
! -------------------------------------------
IF (do_tracer) THEN
  CALL cvt_v (chi_vh_tracer(1:ylat, 1:ylon, 1:ylev), &
              state_tracer(1:ylat, 1:ylon, 1:ylev))
ELSE
  state_tracer(1:ylat, 1:ylon, 1:ylev) = 0.0
END IF

! 3. Do the temporal transform for the flux
! -----------------------------------------
IF (do_flux) THEN
  CALL cvt_t (chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
              state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ELSE
  state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) = 0.0
END IF

! 4. Put the fields into INVICAT format (physical space)
! ------------------------------------------------------
CALL Physical_multid2linear ( state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),   &
                              state_tracer(1:ylat, 1:ylon, 1:ylev),              &
                              zy(1:dim_1d) )

END SUBROUTINE cvt_total





!=================================================================================================
SUBROUTINE cvt_total_adj (chi_hat,           & ! OUT Control variable
                          zy_hat,            & ! IN  Physical space increment in INVICAT format
                          do_flux,           & ! IN  true/false to deal with 2D+time flux field
                          do_tracer)           ! IN  true/false to deal with 3D intial tracer field

! Description:  (adjoint of the following)
! Performs the complete control variable transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Horizontal, followed by vertical (for tracer) or temporal (for flux)
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields, tidy things
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, dim_1d, nfields, gs

USE onedvar_data, ONLY: n1d

! Note : dim_1d  is the size of the physical vector comprising the tracer and flux field
!        n1d     is the size of the control vector space comprising the tracer and flux field
!        nfields is the number of types of flux field (e.g. nfield=2 for anthropogenic and natural fluxes.
!        gs      is the number of points on horizontal grid

IMPLICIT NONE

INCLUDE "cvt.interface"

! Arguments:
REAL,    INTENT(OUT) :: chi_hat(1:n1d)
REAL,    INTENT(IN)  :: zy_hat(1:dim_1d)
LOGICAL, INTENT(IN)  :: do_flux
LOGICAL, INTENT(IN)  :: do_tracer

! Local variables:
REAL                 :: chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: chi_vh_tracer(1:ylat, 1:ylon, 1:ylev)
REAL                 :: state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: state_tracer(1:ylat, 1:ylon, 1:ylev)



! 4. Put the fields into INVICAT format
! -------------------------------------
CALL Physical_linear2multid ( state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  &
                              state_tracer(1:ylat, 1:ylon, 1:ylev),             &
                              zy_hat(1:dim_1d) )


! 3. Do the temporal transform for the flux
! -----------------------------------------
IF (do_flux) THEN
  CALL cvt_t_adj (chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                  state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ELSE
  chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) = 0.0
END IF


! 2. Do the vertical transform for the tracer
! -------------------------------------------
IF (do_tracer) THEN
  CALL cvt_v_adj (chi_vh_tracer(1:ylat, 1:ylon, 1:ylev), &
                  state_tracer(1:ylat, 1:ylon, 1:ylev))
ELSE
  chi_vh_tracer(1:ylat, 1:ylon, 1:ylev) = 0.0
END IF


! 1. Do the horizontal transform
! ------------------------------
CALL cvt_h_adj (chi_hat(1:n1d),                                     &
                chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),   &
                chi_vh_tracer(1:ylat, 1:ylon, 1:ylev),              &
                do_flux,                                            &
                do_tracer)

END SUBROUTINE cvt_total_adj





!=================================================================================================
SUBROUTINE cvt_total_inv (chi,               & ! OUT Control variable
                          zy,                & ! IN  Physical space increment in INVICAT format
                          do_flux,           & ! IN  true/false to deal with 2D+time flux field
                          do_tracer)           ! IN  true/false to deal with 3D intial tracer field

! Description:
! Performs the approximate inverse of the complete control variable transform
! (the approximateness is due to interpolation inside the horizontal transform)
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Horizontal, followed by vertical (for tracer) or temporal (for flux)
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields, tidy things
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, dim_1d, nfields, gs

USE onedvar_data, ONLY: n1d

! Note : dim_1d  is the size of the physical vector comprising the tracer and flux field
!        n1d     is the size of the control vector space comprising the tracer and flux field
!        nfields is the number of types of flux field (e.g. nfield=2 for anthropogenic and natural fluxes.
!        gs      is the number of points on horizontal grid

IMPLICIT NONE

INCLUDE "cvt.interface"

! Arguments:
REAL,    INTENT(OUT) :: chi(1:n1d)
REAL,    INTENT(IN)  :: zy(1:dim_1d)
LOGICAL, INTENT(IN)  :: do_flux
LOGICAL, INTENT(IN)  :: do_tracer

! Local variables:
REAL                 :: chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: chi_vh_tracer(1:ylat, 1:ylon, 1:ylev)
REAL                 :: state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: state_tracer(1:ylat, 1:ylon, 1:ylev)

! -------------------------------------
! 4. Put the fields into INVICAT format
! -------------------------------------

CALL Physical_linear2multid ( state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  &
                              state_tracer(1:ylat, 1:ylon, 1:ylev),             &
                              zy(1:dim_1d) )


! 3. Do the temporal transform for the flux
! -----------------------------------------
IF (do_flux) THEN
  CALL cvt_t_inv (chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                  state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ELSE
  chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) = 0.0
END IF


! 2. Do the vertical transform for the tracer
! -------------------------------------------
IF (do_tracer) THEN
  CALL cvt_v_inv (chi_vh_tracer(1:ylat, 1:ylon, 1:ylev), &
                  state_tracer(1:ylat, 1:ylon, 1:ylev))
ELSE
  chi_vh_tracer(1:ylat, 1:ylon, 1:ylev) = 0.0
END IF


! 1. Do the horizontal transform
! ------------------------------
CALL cvt_h_inv (chi(1:n1d),                                         &
                chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),   &
                chi_vh_tracer(1:ylat, 1:ylon, 1:ylev),              &
                do_flux,                                            &
                do_tracer)

END SUBROUTINE cvt_total_inv




!=================================================================================================
SUBROUTINE cvt_total_inv_adj (chi_hat,           & ! IN  Control variable
                              zy_hat,            & ! OUT Physical space increment in INVICAT format
                              do_flux,           & ! IN  true/false to deal with 2D+time flux field
                              do_tracer)           ! IN  true/false to deal with 3D intial tracer field

! Description:
! Performs the adjoint of the approximate inverse of the complete control variable transform
! (the approximateness is due to interpolation inside the horizontal transform)
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Horizontal, followed by vertical (for tracer) or temporal (for flux)
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields, tidy things (R.N. Bannister)
!           09/2019  Adapted from cvt_total_inv (R.N. Bannister)
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, dim_1d, nfields, gs

USE onedvar_data, ONLY: n1d

! Note : dim_1d  is the size of the physical vector comprising the tracer and flux field
!        n1d     is the size of the control vector space comprising the tracer and flux field
!        nfields is the number of types of flux field (e.g. nfield=2 for anthropogenic and natural fluxes.
!        gs      is the number of points on horizontal grid

IMPLICIT NONE

INCLUDE "cvt.interface"

! Arguments:
REAL,    INTENT(IN)  :: chi_hat(1:n1d)
REAL,    INTENT(OUT) :: zy_hat(1:dim_1d)
LOGICAL, INTENT(IN)  :: do_flux
LOGICAL, INTENT(IN)  :: do_tracer

! Local variables:
REAL                 :: chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev)
REAL                 :: state_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL                 :: state_tracer_hat(1:ylat, 1:ylon, 1:ylev)


! 1. Do the horizontal transform
! ------------------------------
CALL cvt_h_inv_adj (chi_hat(1:n1d),                                         &
                    chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields),   &
                    chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev),              &
                    do_flux,                                                &
                    do_tracer)

! 2. Do the vertical transform for the tracer
! -------------------------------------------
IF (do_tracer) THEN
  CALL cvt_v_inv_adj (chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev), &
                      state_tracer_hat(1:ylat, 1:ylon, 1:ylev))
ELSE
  chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev) = 0.0
END IF


! 3. Do the temporal transform for the flux
! -----------------------------------------
IF (do_flux) THEN
  CALL cvt_t_inv_adj (chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                      state_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ELSE
  chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields) = 0.0
END IF


! -------------------------------------
! 4. Put the fields into INVICAT format
! -------------------------------------
CALL Physical_multid2linear ( state_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  &
                              state_tracer_hat(1:ylat, 1:ylon, 1:ylev),             &
                              zy_hat(1:dim_1d) )

END SUBROUTINE cvt_total_inv_adj





!=================================================================================================
SUBROUTINE cvt_h (chi,           & ! IN  Control variable
                  chi_th_flux,   & ! OUT Version of chi (for flux) between temp and horiz transforms
                  chi_vh_tracer, & ! OUT Version of chi (for tracer) between vert and horiz transforms
                  do_flux,       & ! IN  true/false to deal with 2D+time flux field
                  do_tracer,     & ! IN  true/false to deal with 3D intial tracer field
                  do_nlevs,      & ! IN, OPTIONAL to do transform on first nlevs only (tracer)
                  do_nmonths)      ! IN, OPTIONAL to do transform on first nlevs only (flux)

! Description:
! Performs the horizontal control variable transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs an weighting of the input followed by an inverse spectral transform
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      04/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields
!           11/2019  Apply to non-balanced transform
!           12/2019  Add do_adjustment, do_nlevs, do_nmonths, and introduce correlation adjustments
!           12/2019  Remove do_adjustment, make transform 'balanced'
!           01/2020  Revert to unbalanced, use my own transform code
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields

USE onedvar_data, ONLY: n1d,                     &
                        L,                       &
                        ALPsize,                 &
                        CVT_HorizCors

IMPLICIT NONE

! Note: We set L=ylat
! The number of longitudes is ylon in tomcat
!                             2L+1 in shtools
! The number of latitudes is  ylat in tomcat
!                             L+1 in shtools


! Arguments:
REAL,              INTENT(IN)  :: chi(1:n1d)
REAL,              INTENT(OUT) :: chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL,              INTENT(OUT) :: chi_vh_tracer(1:ylat, 1:ylon, 1:ylev)
LOGICAL,           INTENT(IN)  :: do_flux
LOGICAL,           INTENT(IN)  :: do_tracer
INTEGER, OPTIONAL, INTENT(IN)  :: do_nlevs
INTEGER, OPTIONAL, INTENT(IN)  :: do_nmonths

! Local variables:
INTEGER                        :: ll, mm, i, lev, mon, index_1d, err, f, nlevs, nmonths, lon
REAL                           :: chi1_flux(1:2, 0:L, 0:L, 1:nmonth, 1:nfields)
REAL                           :: chi1_tracer(1:2, 0:L, 0:L, 1:ylev)
REAL                           :: chi2_flux(1:L+1, 1:2*L+1, 1:nmonth, 1:nfields)
REAL                           :: chi2_tracer(1:L+1, 1:2*L+1, 1:ylev)
REAL                           :: weight


IF (PRESENT(do_nlevs)) THEN
  nlevs = do_nlevs
ELSE
  nlevs = ylev
END IF

IF (PRESENT(do_nmonths)) THEN
  nmonths = do_nmonths
ELSE
  nmonths = nmonth
END IF


! 1. Multiply control variable by the horizontal spectrum
!    and put into a more meaningful data structure
! -------------------------------------------------------
index_1d = 0

! 1a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -

IF (do_flux) THEN
  chi1_flux(1:2, 0:L, 0:L, 1:nmonths, 1:nfields) = 0.0
  ! Loop over each flux type
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      ! Loop over total wavenumber
      DO ll = 0, L
        weight = CVT_HorizCors % std_hspec_flux(ll,mon,f)
        ! Special case for mm=0 (no Sine contribution)
        ! The cosine term
        index_1d = index_1d + 1
        chi1_flux(1,ll,0,mon,f) = chi(index_1d) * weight
        IF (ll > 0) THEN
          ! Loop over remaining meridional wavenumbers
          DO mm = 1, ll
            ! The cosine term
            index_1d = index_1d + 1
            chi1_flux(1,ll,mm,mon,f) = chi(index_1d) * weight
            ! The sine term
            index_1d = index_1d + 1
            chi1_flux(2,ll,mm,mon,f) = chi(index_1d) * weight
          END DO
        END IF
      END DO
    END DO
  END DO
END IF


! 1b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  chi1_tracer(1:2, 0:L, 0:L, 1:nlevs) = 0.0
  ! Loop over each vertical level
  DO lev = 1, nlevs
    ! Loop over total wavenumber
    DO ll = 0, L
      weight = CVT_HorizCors % std_hspec_tracer(ll,lev)
      ! Special case for mm=0 (no Sine contribution)
      ! The cosine term
      index_1d = index_1d + 1
      chi1_tracer(1,ll,0,lev) = chi(index_1d) * weight
      IF (ll > 0) THEN
        ! Loop over remaining meridional wavenumbers
        DO mm = 1, ll
          ! The cosine term
          index_1d = index_1d + 1
          chi1_tracer(1,ll,mm,lev) = chi(index_1d) * weight
          ! The sine term
          index_1d = index_1d + 1
          chi1_tracer(2,ll,mm,lev) = chi(index_1d) * weight
        END DO
      END IF
    END DO
  END DO
END IF

!PRINT *, 'Final value of index_1d =', index_1d


! 2. Perform Spherical transform operation (to real space)
! -------------------------------------------------------
! Recall that the output will be a function of **co-latitude instead of latitude**

! 2a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  ! Loop over each flux type
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      ! Do spectral transform
      CALL Spherical_forward (chi1_flux(1:2, 0:L, 0:L, mon, f),      &  ! IN
                              chi2_flux(1:L+1, 1:2*L+1, mon, f))        ! OUT
    END DO
  END DO
END IF


! 2b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  ! Loop over each vertical level
  DO lev = 1, nlevs
    ! Do spectral transform
    CALL Spherical_forward (chi1_tracer(1:2, 0:L, 0:L, lev),      &  ! IN
                            chi2_tracer(1:L+1, 1:2*L+1, lev))        ! OUT
  END DO
END IF


! 3. Interpolate to same grid as TOMCAT
! -------------------------------------------------------

! 3a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  DO f = 1, nfields
    CALL ChangeHorizGrid (nmonths,                                  & ! No of months
                          2*L+1,                                    & ! No of longitudes shtools grid
                          L+1,                                      & ! No of latitudes shtools grid
                          chi2_flux(1:L+1, 1:2*L+1, 1:nmonths, f),  & ! Output from spectral trans
                          CVT_HorizCors % SHtools_longs(1:2*L+1),   & ! Longitudes of shtools grid
                          CVT_HorizCors % SHtools_lats(1:L+1),      & ! Latitudes of shtools grid
                          ylon,                                     & ! No of longitudes tomcat grid
                          ylat,                                     & ! No of latitudes tomcat grid
                          chi_th_flux(1:ylat, 1:ylon, 1:nmonths, f),& ! On tomcat grid
                          CVT_HorizCors % TOMCAT_longs(1:ylon),     & ! Longs on tomcat grid
                          CVT_HorizCors % TOMCAT_lats(1:ylat),      & ! Lats on tomcat grid
                          CVT_HorizCors % SH_to_TOMCAT)               ! Details of the grid change
  END DO
END IF

! 3b. Tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  CALL ChangeHorizGrid (nlevs,                                 & ! No of levels
                        2*L+1,                                 & ! No of longitudes shtools grid
                        L+1,                                   & ! No of latitudes shtools grid
                        chi2_tracer(1:L+1, 1:2*L+1, 1:nlevs),  & ! Output from spectral trans
                        CVT_HorizCors % SHtools_longs(1:2*L+1),& ! Longitudes of shtools grid
                        CVT_HorizCors % SHtools_lats(1:L+1),   & ! Latitudes of shtools grid
                        ylon,                                  & ! No of longitudes tomcat grid
                        ylat,                                  & ! No of latitudes tomcat grid
                        chi_vh_tracer(1:ylat, 1:ylon, 1:nlevs),& ! On tomcat grid
                        CVT_HorizCors % TOMCAT_longs(1:ylon),  & ! Longs on tomcat grid
                        CVT_HorizCors % TOMCAT_lats(1:ylat),   & ! Lats on tomcat grid
                        CVT_HorizCors % SH_to_TOMCAT)            ! Details of the grid change
END IF

END SUBROUTINE cvt_h




!=================================================================================================
SUBROUTINE cvt_h_adj (chi_hat,           & ! OUT Control variable
                      chi_th_flux_hat,   & ! IN  Version of chi (for flux) between temp and horiz transforms
                      chi_vh_tracer_hat, & ! IN  Version of chi (for tracer) between vert and horiz transforms
                      do_flux,           & ! IN  true/false to deal with 2D+time flux field
                      do_tracer,         & ! IN  true/false to deal with 3D intial tracer field
                      do_nlevs,          & ! IN, OPTIONAL to do transform on first nlevs only (tracer)
                      do_nmonths)          ! IN, OPTIONAL to do transform on first nlevs only (flux)

! Description:  (adjoint of the following)
! Performs the horizontal control variable transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs an weighting of the input followed by an inverse spectral transform
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      04/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields
!           11/2019  Apply to non-balanced transform
!           12/2019  Add do_adjustment, do_nlevs, do_nmonths, and introduce correlation adjustments
!           12/2019  Remove do_adjustment, make transform 'balanced'
!           01/2020  Revert to unbalanced, use my own transform code
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields

USE onedvar_data, ONLY: n1d,                     &
                        L,                       &
                        ALPsize,                 &
                        CVT_HorizCors

USE cparam,       ONLY: pi


IMPLICIT NONE

! Arguments:
REAL,              INTENT(OUT) :: chi_hat(1:n1d)
REAL,              INTENT(IN)  :: chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL,              INTENT(IN)  :: chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev)
LOGICAL,           INTENT(IN)  :: do_flux
LOGICAL,           INTENT(IN)  :: do_tracer
INTEGER, OPTIONAL, INTENT(IN)  :: do_nlevs
INTEGER, OPTIONAL, INTENT(IN)  :: do_nmonths

! Local variables:
INTEGER                        :: ll, mm, i, lev, mon, index_1d, err, f, nlevs, nmonths, lon
REAL                           :: chi1_flux_hat(1:2, 0:L, 0:L, 1:nmonth, 1:nfields)
REAL                           :: chi1_tracer_hat(1:2, 0:L, 0:L, 1:ylev)
REAL                           :: chi2_flux_hat(1:L+1, 1:2*L+1, 1:nmonth, 1:nfields)
REAL                           :: chi2_tracer_hat(1:L+1, 1:2*L+1, 1:ylev)

REAL                           :: weight


IF (PRESENT(do_nlevs)) THEN
  nlevs = do_nlevs
ELSE
  nlevs = ylev
END IF

IF (PRESENT(do_nmonths)) THEN
  nmonths = do_nmonths
ELSE
  nmonths = nmonth
END IF


! 3. Interpolate to same grid as TOMCAT
! -------------------------------------------------------

! 3a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  DO f = 1, nfields
    CALL ChangeHorizGrid_adj (nmonths,                                      & ! No of months
                              2*L+1,                                        & ! No of longitudes shtools grid
                              L+1,                                          & ! No of latitudes shtools grid
                              chi2_flux_hat(1:L+1, 1:2*L+1, 1:nmonths, f),  & ! Output from spectral trans
                              CVT_HorizCors % SHtools_longs(1:2*L+1),       & ! Longitudes of shtools grid
                              CVT_HorizCors % SHtools_lats(1:L+1),          & ! Latitudes of shtools grid
                              ylon,                                         & ! No of longitudes tomcat grid
                              ylat,                                         & ! No of latitudes tomcat grid
                              chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonths, f),& ! On tomcat grid
                              CVT_HorizCors % TOMCAT_longs(1:ylon),         & ! Longs on tomcat grid
                              CVT_HorizCors % TOMCAT_lats(1:ylat),          & ! Lats on tomcat grid
                              CVT_HorizCors % SH_to_TOMCAT)                   ! Details of grid change
  END DO
END IF

! 3b. Tracer
! - - - - - - - - - - - - - - - - - - - - - - - -

IF (do_tracer) THEN
  CALL ChangeHorizGrid_adj (nlevs,                                     & ! No of levels
                            2*L+1,                                     & ! No of longitudes shtools grid
                            L+1,                                       & ! No of latitudes shtools grid
                            chi2_tracer_hat(1:L+1, 1:2*L+1, 1:nlevs),  & ! Output from spectral trans
                            CVT_HorizCors % SHtools_longs(1:2*L+1),    & ! Longitudes of shtools grid
                            CVT_HorizCors % SHtools_lats(1:L+1),       & ! Latitudes of shtools grid
                            ylon,                                      & ! No of longitudes tomcat grid
                            ylat,                                      & ! No of latitudes tomcat grid
                            chi_vh_tracer_hat(1:ylat, 1:ylon, 1:nlevs),& ! On tomcat grid
                            CVT_HorizCors % TOMCAT_longs(1:ylon),      & ! Longs on tomcat grid
                            CVT_HorizCors % TOMCAT_lats(1:ylat),       & ! Lats on tomcat grid
                            CVT_HorizCors % SH_to_TOMCAT)                ! Details of grid change
END IF


! 2. Perform Spherical transform operation (to real space)
! -------------------------------------------------------

! 2a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  ! Loop over each flux type
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      CALL Spherical_adj (chi1_flux_hat(1:2, 0:L, 0:L, mon, f), & ! OUT
                          chi2_flux_hat(1:L+1, 1:2*L+1, mon, f))  ! IN
    END DO
  END DO
END IF


! 2b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  ! Loop over each vertical level
  DO lev = 1, nlevs
    CALL Spherical_adj (chi1_tracer_hat(1:2, 0:L, 0:L, lev), & ! OUT
                        chi2_tracer_hat(1:L+1, 1:2*L+1, lev))  ! IN
  END DO
END IF


! 1. Multiply control variable by the horizontal spectrum
!    and put into a more meaningful data structure
! -------------------------------------------------------
index_1d = 0

!chi_hat(1:n1d) = 0.0

! 1a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -

IF (do_flux) THEN
  ! Loop over each month (or timescale, depending on option)
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      ! Loop over total wavenumber
      DO ll = 0, L
        weight = CVT_HorizCors % std_hspec_flux(ll,mon,f)
        ! Special case for mm=0 (no Sine contribution)
        ! The cosine term
        index_1d = index_1d + 1
        ! chi1_flux(1,ll,0,mon,f) = chi(index_1d) * weight
        chi_hat(index_1d) = chi1_flux_hat(1,ll,0,mon,f) * weight
        IF (ll > 0) THEN
          ! Loop over remaining meridional wavenumbers
          DO mm = 1, ll
            ! The cosine term
            index_1d = index_1d + 1
            ! chi1_flux(1,ll,mm,mon,f) = chi(index_1d) * weight
            chi_hat(index_1d) = chi1_flux_hat(1,ll,mm,mon,f) * weight
            ! The sine term
            index_1d = index_1d + 1
            ! chi1_flux(2,ll,mm,mon,f) = chi(index_1d) * weight
            chi_hat(index_1d) = chi1_flux_hat(2,ll,mm,mon,f) * weight
          END DO
        END IF
      END DO
    END DO
  END DO
END IF


! 1b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  ! Loop over each vertical level
  DO lev = 1, nlevs
    ! Loop over total wavenumber
    DO ll = 0, L
      weight = CVT_HorizCors % std_hspec_tracer(ll,lev)
      ! Special case for mm=0 (no Sine contribution)
      ! The cosine term
      index_1d = index_1d + 1
      ! chi1_tracer(1,ll,0,lev) = chi(index_1d) * weight
      chi_hat(index_1d) = chi1_tracer_hat(1,ll,0,lev) * weight
      IF (ll > 0) THEN
        ! Loop over remaining meridional wavenumbers
        DO mm = 1, ll
          ! The cosine term
          index_1d = index_1d + 1
          ! chi1_tracer(1,ll,mm,lev) = chi(index_1d) * weight
          chi_hat(index_1d) = chi1_tracer_hat(1,ll,mm,lev) * weight
          ! The sine term
          index_1d = index_1d + 1
          ! chi1_tracer(2,ll,mm,lev) = chi(index_1d) * weight
          chi_hat(index_1d) = chi1_tracer_hat(2,ll,mm,lev) * weight
        END DO
      END IF
    END DO
  END DO
END IF

END SUBROUTINE cvt_h_adj



!=================================================================================================
SUBROUTINE cvt_h_inv (chi,               & ! OUT Control variable
                      chi_th_flux,       & ! IN  Version of chi (for flux) between temp and horiz transforms
                      chi_vh_tracer,     & ! IN  Version of chi (for tracer) between vert and horiz transforms
                      do_flux,           & ! IN  true/false to deal with 2D+time flux field
                      do_tracer,         & ! IN  true/false to deal with 3D intial tracer field
                      do_nlevs,          & ! IN, OPTIONAL to do transform on first nlevs only (tracer)
                      do_nmonths)          ! IN, OPTIONAL to do transform on first nlevs only (flux)

! Description:
! Performs the approximate inverse horizontal control variable transform
! (the approximateness is due to the horizontal interpolation)
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs an weighting of the input followed by an inverse spectral transform
! This is an apprximate inverse (due to the change of grids)
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      04/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields
!           11/2019  Apply to non-balanced transform
!           12/2019  Add do_adjustment, do_nlevs, do_nmonths, and introduce correlation adjustments
!           12/2019  Remove do_adjustment, make transform 'balanced'
!           01/2020  Revert to unbalanced, use my own transform code
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields

USE onedvar_data, ONLY: n1d,                     &
                        L,                       &
                        ALPsize,                 &
                        CVT_HorizCors

USE cparam,       ONLY: pi

IMPLICIT NONE

! Arguments:
REAL,              INTENT(OUT) :: chi(1:n1d)
REAL,              INTENT(IN)  :: chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL,              INTENT(IN)  :: chi_vh_tracer(1:ylat, 1:ylon, 1:ylev)
LOGICAL,           INTENT(IN)  :: do_flux
LOGICAL,           INTENT(IN)  :: do_tracer
INTEGER, OPTIONAL, INTENT(IN)  :: do_nlevs
INTEGER, OPTIONAL, INTENT(IN)  :: do_nmonths

! Local variables:
INTEGER                        :: ll, mm, i, lev, mon, index_1d, err, f, nlevs, nmonths, lon
REAL                           :: chi1_flux(1:2, 0:L, 0:L, 1:nmonth, 1:nfields)
REAL                           :: chi1_tracer(1:2, 0:L, 0:L, 1:ylev)
REAL                           :: chi2_flux(1:L+1, 1:2*L+1, 1:nmonth, 1:nfields)
REAL                           :: chi2_tracer(1:L+1, 1:2*L+1, 1:ylev)
REAL                           :: weight


IF (PRESENT(do_nlevs)) THEN
  nlevs = do_nlevs
ELSE
  nlevs = ylev
END IF

IF (PRESENT(do_nmonths)) THEN
  nmonths = do_nmonths
ELSE
  nmonths = nmonth
END IF


! 3. Interpolate to same grid as TOMCAT
! -------------------------------------------------------

! 3a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  DO f = 1, nfields
    CALL ChangeHorizGrid (nmonths,                                  & ! No of months
                          ylon,                                     & ! No of longitudes tomcat grid
                          ylat,                                     & ! No of latitudes tomcat grid
                          chi_th_flux(1:ylat, 1:ylon, 1:nmonths, f),& ! Data on tomcat grid
                          CVT_HorizCors % TOMCAT_longs(1:ylon),     & ! Longs on tomcat grid
                          CVT_HorizCors % TOMCAT_lats(1:ylat),      & ! Lats on tomcat grid
                          2*L+1,                                    & ! No of longitudes shtools grid
                          L+1,                                      & ! No of latitudes shtools grid
                          chi2_flux(1:L+1, 1:2*L+1, 1:nmonths, f),  & ! Data on shtools grid
                          CVT_HorizCors % SHtools_longs(1:2*L+1),   & ! Longitudes of shtools grid
                          CVT_HorizCors % SHtools_lats(1:L+1),      & ! Latitudes of shtools grid
                          CVT_HorizCors % TOMCAT_to_SH)               ! Details of grid change
  END DO
END IF

! 3b. Tracer
! - - - - - - - - - - - - - - - - - - - - - - - -

IF (do_tracer) THEN
  CALL ChangeHorizGrid (nlevs,                                 & ! No of levels
                        ylon,                                  & ! No of longitudes tomcat grid
                        ylat,                                  & ! No of latitudes tomcat grid
                        chi_vh_tracer(1:ylat, 1:ylon, 1:nlevs),& ! Data on tomcat grid
                        CVT_HorizCors % TOMCAT_longs(1:ylon),  & ! Longs on tomcat grid
                        CVT_HorizCors % TOMCAT_lats(1:ylat),   & ! Lats on tomcat grid
                        2*L+1,                                 & ! No of longitudes shtools grid
                        L+1,                                   & ! No of latitudes shtools grid
                        chi2_tracer(1:L+1, 1:2*L+1, 1:nlevs),  & ! Data on shtools grid
                        CVT_HorizCors % SHtools_longs(1:2*L+1),& ! Longitudes of shtools grid
                        CVT_HorizCors % SHtools_lats(1:L+1),   & ! Latitudes of shtools grid
                        CVT_HorizCors % TOMCAT_to_SH)            ! Details of grid change
END IF



! 2. Perform Spherical transform operation (to spectral space)
! -------------------------------------------------------

! 2a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  ! Loop over each flux type
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      CALL Spherical_inv (chi1_flux(1:2, 0:L, 0:L, mon, f),    & ! OUT
                          chi2_flux(1:L+1, 1:2*L+1, mon, f))     ! IN
   END DO
  END DO
END IF


! 2b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  ! Loop over each vertical level
  DO lev = 1, nlevs
    CALL Spherical_inv (chi1_tracer(1:2, 0:L, 0:L, lev),    & ! OUT
                        chi2_tracer(1:L+1, 1:2*L+1, lev))     ! IN
  END DO
END IF



! 1. Divide control variable by the horizontal spectrum
!    and put into a more meaningful data structure
! -------------------------------------------------------
index_1d = 0

!chi_hat(1:n1d) = 0.0

! 1a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -

IF (do_flux) THEN
  ! Loop over each surface flux field
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      ! Loop over total wavenumber
      DO ll = 0, L
        weight = 1.0 / CVT_HorizCors % std_hspec_flux(ll,mon,f)
        ! Special case for mm=0 (no Sine contribution)
        ! The cosine term
        index_1d = index_1d + 1
        chi(index_1d) = chi1_flux(1,ll,0,mon,f) * weight
        IF (ll > 0) THEN
          ! Loop over remaining meridional wavenumbers
          DO mm = 1, ll
            ! The cosine term
            index_1d = index_1d + 1
            chi(index_1d) = chi1_flux(1,ll,mm,mon,f) * weight
            ! The sine term
            index_1d = index_1d + 1
            chi(index_1d) = chi1_flux(2,ll,mm,mon,f) * weight
          END DO
        END IF
      END DO
    END DO
  END DO
END IF


! 1b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  ! Loop over each vertical level
  DO lev = 1, nlevs
    ! Loop over total wavenumber
    DO ll = 0, L
      weight = 1.0 / CVT_HorizCors % std_hspec_tracer(ll,lev)
      ! Special case for mm=0 (no Sine contribution)
      ! The cosine term
      index_1d = index_1d + 1
      ! chi1_tracer(1,ll,0,lev) = chi(index_1d) * weight
      chi(index_1d) = chi1_tracer(1,ll,0,lev) * weight
      IF (ll > 0) THEN
        ! Loop over remaining meridional wavenumbers
        DO mm = 1, ll
          ! The cosine term
          index_1d = index_1d + 1
          chi(index_1d) = chi1_tracer(1,ll,mm,lev) * weight
          ! The sine term
          index_1d = index_1d + 1
          chi(index_1d) = chi1_tracer(2,ll,mm,lev) * weight
        END DO
      END IF
    END DO
  END DO
END IF


END SUBROUTINE cvt_h_inv





!=================================================================================================
SUBROUTINE cvt_h_inv_adj (chi_hat,               & ! IN  Control variable
                          chi_th_flux_hat,       & ! OUT Version of chi (for flux) between temp and horiz transforms
                          chi_vh_tracer_hat,     & ! OUT Version of chi (for tracer) between vert and horiz transforms
                          do_flux,               & ! IN  true/false to deal with 2D+time flux field
                          do_tracer,             & ! IN  true/false to deal with 3D intial tracer field
                          do_nlevs,              & ! IN, OPTIONAL to do transform on first nlevs only (tracer)
                          do_nmonths)              ! IN, OPTIONAL to do transform on first nlevs only (flux)

! Description:
! Performs the adjoint of the approximate inverse horizontal control variable transform
! (the approximateness is due to the horizontal interpolation)
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs an weighting of the input followed by an inverse spectral transform
! This is an apprximate inverse (due to the change of grids)
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      04/2019  Initial version (R.N. Bannister)
!           09/2019  Add nfields
!           11/2019  Apply to non-balanced transform
!           12/2019  Add do_adjustment, do_nlevs, do_nmonths, and introduce correlation adjustments
!           12/2019  Remove do_adjustment, make transform 'balanced'
!           01/2020  Revert to unbalanced, use my own transform code
!
! Code Description:
!   Language:           Fortran 90.


! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields

USE onedvar_data, ONLY: n1d,                     &
                        L,                       &
                        ALPsize,                 &
                        CVT_HorizCors

USE cparam,       ONLY: pi

IMPLICIT NONE

! Arguments:
REAL,              INTENT(IN)  :: chi_hat(1:n1d)
REAL,              INTENT(OUT) :: chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL,              INTENT(OUT) :: chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev)
LOGICAL,           INTENT(IN)  :: do_flux
LOGICAL,           INTENT(IN)  :: do_tracer
INTEGER, OPTIONAL, INTENT(IN)  :: do_nlevs
INTEGER, OPTIONAL, INTENT(IN)  :: do_nmonths

! Local variables:
INTEGER                        :: ll, mm, i, lev, mon, index_1d, err, f, nlevs, nmonths, lon
REAL                           :: chi1_flux_hat(1:2, 0:L, 0:L, 1:nmonth, 1:nfields)
REAL                           :: chi1_tracer_hat(1:2, 0:L, 0:L, 1:ylev)
REAL                           :: chi2_flux_hat(1:L+1, 1:2*L+1, 1:nmonth, 1:nfields)
REAL                           :: chi2_tracer_hat(1:L+1, 1:2*L+1, 1:ylev)
REAL                           :: weight


IF (PRESENT(do_nlevs)) THEN
  nlevs = do_nlevs
ELSE
  nlevs = ylev
END IF

IF (PRESENT(do_nmonths)) THEN
  nmonths = do_nmonths
ELSE
  nmonths = nmonth
END IF


! 1. Divide control variable by the horizontal spectrum
!    and put into a more meaningful data structure
! -------------------------------------------------------
index_1d = 0

! 1a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -

IF (do_flux) THEN
  ! Loop over each surface flux field
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      ! Loop over total wavenumber
      DO ll = 0, L
        weight = 1.0 / CVT_HorizCors % std_hspec_flux(ll,mon,f)
        ! Special case for mm=0 (no Sine contribution)
        ! The cosine term
        index_1d = index_1d + 1
        ! chi(index_1d) = chi1_flux(1,ll,0,mon,f) * weight
        chi1_flux_hat(1,ll,0,mon,f) = chi_hat(index_1d) * weight
        IF (ll > 0) THEN
          ! Loop over remaining meridional wavenumbers
          DO mm = 1, ll
            ! The cosine term
            index_1d = index_1d + 1
            ! chi(index_1d) = chi1_flux(1,ll,mm,mon,f) * weight
            chi1_flux_hat(1,ll,mm,mon,f) = chi_hat(index_1d) * weight
            ! The sine term
            index_1d = index_1d + 1
            ! chi(index_1d) = chi1_flux(2,ll,mm,mon,f) * weight
            chi1_flux_hat(2,ll,mm,mon,f) = chi_hat(index_1d) * weight
          END DO
        END IF
      END DO
    END DO
  END DO
END IF


! 1b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  ! Loop over each vertical level
  DO lev = 1, nlevs
    ! Loop over total wavenumber
    DO ll = 0, L
      weight = 1.0 / CVT_HorizCors % std_hspec_tracer(ll,lev)
      ! Special case for mm=0 (no Sine contribution)
      ! The cosine term
      index_1d = index_1d + 1
      ! chi(index_1d) = chi1_tracer(1,ll,0,lev) * weight
      chi1_tracer_hat(1,ll,0,lev) = chi_hat(index_1d) * weight
      IF (ll > 0) THEN
        ! Loop over remaining meridional wavenumbers
        DO mm = 1, ll
          ! The cosine term
          index_1d = index_1d + 1
          ! chi(index_1d) = chi1_tracer(1,ll,mm,lev) * weight
          chi1_tracer_hat(1,ll,mm,lev) = chi_hat(index_1d) * weight
          ! The sine term
          index_1d = index_1d + 1
          ! chi(index_1d) = chi1_tracer(2,ll,mm,lev) * weight
          chi1_tracer_hat(2,ll,mm,lev) = chi_hat(index_1d) * weight
        END DO
      END IF
    END DO
  END DO
END IF



! 2. Perform Spherical transform operation (to real space)
! -------------------------------------------------------

! 2a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  ! Loop over each flux type
  DO f = 1, nfields
    ! Loop over each month (or timescale, depending on option)
    DO mon = 1, nmonths
      CALL Spherical_inv_adj (chi1_flux_hat(1:2, 0:L, 0:L, mon, f), & ! IN
                              chi2_flux_hat(1:L+1, 1:2*L+1, mon, f))
    END DO
  END DO
END IF


! 2b. Initial tracer
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_tracer) THEN
  ! Loop over each vertical level
  DO lev = 1, nlevs
    CALL Spherical_inv_adj (chi1_tracer_hat(1:2, 0:L, 0:L, lev), & ! IN
                            chi2_tracer_hat(1:L+1, 1:2*L+1, lev))
  END DO
END IF



! 3. Interpolate to same grid as TOMCAT
! -------------------------------------------------------

! 3a. Surface flux
! - - - - - - - - - - - - - - - - - - - - - - - -
IF (do_flux) THEN
  DO f = 1, nfields
    CALL ChangeHorizGrid_adj (nmonths,                                      & ! No of months
                              ylon,                                         & ! No of longitudes tomcat grid
                              ylat,                                         & ! No of latitudes tomcat grid
                              chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonths, f),& ! Data on tomcat grid
                              CVT_HorizCors % TOMCAT_longs(1:ylon),         & ! Longs on tomcat grid
                              CVT_HorizCors % TOMCAT_lats(1:ylat),          & ! Lats on tomcat grid
                              2*L+1,                                        & ! No of longitudes shtools grid
                              L+1,                                          & ! No of latitudes shtools grid
                              chi2_flux_hat(1:L+1, 1:2*L+1, 1:nmonths, f),  & ! Data on shtools grid
                              CVT_HorizCors % SHtools_longs(1:2*L+1),       & ! Longitudes of shtools grid
                              CVT_HorizCors % SHtools_lats(1:L+1),          & ! Latitudes of shtools grid
                              CVT_HorizCors % TOMCAT_to_SH)                   ! Details of grid change
  END DO
END IF

! 3b. Tracer
! - - - - - - - - - - - - - - - - - - - - - - - -

IF (do_tracer) THEN
  CALL ChangeHorizGrid_adj (nlevs,                                 & ! No of levels
                        ylon,                                      & ! No of longitudes tomcat grid
                        ylat,                                      & ! No of latitudes tomcat grid
                        chi_vh_tracer_hat(1:ylat, 1:ylon, 1:nlevs),& ! Data on tomcat grid
                        CVT_HorizCors % TOMCAT_longs(1:ylon),      & ! Longs on tomcat grid
                        CVT_HorizCors % TOMCAT_lats(1:ylat),       & ! Lats on tomcat grid
                        2*L+1,                                     & ! No of longitudes shtools grid
                        L+1,                                       & ! No of latitudes shtools grid
                        chi2_tracer_hat(1:L+1, 1:2*L+1, 1:nlevs),  & ! Data on shtools grid
                        CVT_HorizCors % SHtools_longs(1:2*L+1),    & ! Longitudes of shtools grid
                        CVT_HorizCors % SHtools_lats(1:L+1),       & ! Latitudes of shtools grid
                        CVT_HorizCors % TOMCAT_to_SH)                ! Details of grid change
END IF

END SUBROUTINE cvt_h_inv_adj




!=================================================================================================
SUBROUTINE Spherical_forward (SpecCoeffs,        &
                              Field)

! Subroutine to perform forward spherical transform
! Note:
! The number of longitudes is 2L+1 in shtools
! The number of latitudes is  L+1 in shtools

USE onedvar_data, ONLY: L,                       &
                        CVT_HorizCors,           &
                        fft_worklen

IMPLICIT NONE

! Subroutine parameters
REAL, INTENT(IN)  :: SpecCoeffs(1:2, 0:L, 0:L)  ! cos, sin parts; l; m
REAL, INTENT(OUT) :: Field(1:L+1,1:2*L+1)

! Local variables
INTEGER           :: cs, mm, ll, lat, index, index_cos, index_sin, ierr
INTEGER           :: nlon
REAL              :: Inter1(1:2, 0:L, 1:L+1)    ! cos, sin parts; m; lat
REAL              :: Inter2(0:2*L)              ! inteface to FFT
REAL              :: total

nlon           = 2*L+1

! Do the Legendre transform
DO cs = 1, 2
  DO mm = 0, L
    IF ((mm > 0) .OR. (cs == 1)) THEN            ! no sin term for mm=0
      DO lat = 1, L+1
        total = 0.0
        DO ll = mm, L
          index = CVT_HorizCors % Plm_index(ll,mm)
          total = total + SpecCoeffs(cs, ll, mm) * CVT_HorizCors % assocLegPoly(lat,index)
        END DO
        Inter1(cs, mm, lat) = total
      END DO
    END IF
  END DO
END DO


! Do the Fourier transform
DO lat = 1, L+1
  ! Put the data for this latitude in the right format for the FT
  ! This involves putting cos and sin terms in a 1d array and scaling
  ! to match the definition used by fftpack
  Inter2(0) = Inter1(1,0,lat)             ! mm = 0 term (cos only)
  ! The remaining mm
  DO mm = 1, L
    index_cos         = 2*mm-1
    index_sin         = 2*mm
    Inter2(index_cos) = Inter1(1,mm,lat)  ! cos term
    Inter2(index_sin) = Inter1(2,mm,lat)  ! sin term
  END DO

  ! Do the FT back to real space
  CALL rfft1b (nlon, 1,                                           &
               Inter2(0:2*L), nlon,                               &
               CVT_HorizCors % fft_wsave,                         &
               fft_worklen,                                       &
               CVT_HorizCors % fft_work,                          &
               nlon, ierr)

  ! Put data into output array
  Field(lat,1:nlon) = Inter2(0:2*L)
END DO

END SUBROUTINE Spherical_forward



!=================================================================================================
SUBROUTINE Spherical_inv (SpecCoeffs,            &
                          Field)

! Subroutine to perform inverse spherical transform
! Note:
! The number of longitudes is 2L+1 in shtools
! The number of latitudes is  L+1 in shtools

USE onedvar_data, ONLY: L,                       &
                        CVT_HorizCors,           &
                        fft_worklen

IMPLICIT NONE

! Subroutine parameters
REAL, INTENT(OUT) :: SpecCoeffs(1:2, 0:L, 0:L)  ! cos, sin parts; l; m
REAL, INTENT(IN)  :: Field(1:L+1,1:2*L+1)

! Local variables
INTEGER           :: cs, mm, ll, lat, index, index_cos, index_sin, ierr
INTEGER           :: nlon
REAL              :: Inter1(1:2, 0:L, 1:L+1)    ! cos, sin parts; m; lat
REAL              :: Inter2(0:2*L)              ! inteface to FFT

nlon           = 2*L+1

! Do the Fourier transform
DO lat = 1, L+1
  ! Put data into temporary structure
  Inter2(0:2*L) = Field(lat,1:nlon)
  ! Call fft
  CALL rfft1f (nlon, 1,                                           &
               Inter2(0:2*L), nlon,                               &
               CVT_HorizCors % fft_wsave,                         &
               fft_worklen,                                       &
               CVT_HorizCors % fft_work,                          &
               nlon, ierr)

  ! Place in the correct place for the Legendre transform, and do scaling in prep of the Legendre step
  Inter1(1,0,lat) = Inter2(0) * 0.5              ! mm = 0 term (cos only)
  ! The remaining mm
  DO mm = 1, L
    index_cos         = 2*mm-1
    index_sin         = 2*mm
    Inter1(1,mm,lat)  = Inter2(index_cos) * 0.25 ! cos term
    Inter1(2,mm,lat)  = Inter2(index_sin) * 0.25 ! sin term
  END DO
END DO


! Do the Legendre transform
SpecCoeffs(1:2, 0:L, 0:L) = 0.0
DO cs = 1, 2
  DO mm = 0, L
    IF ((mm > 0) .OR. (cs == 1)) THEN          ! no sin term for mm=0
      DO ll = mm, L
        index                  = CVT_HorizCors % Plm_index(ll,mm)
        SpecCoeffs(cs, ll, mm) = SUM(Inter1(cs, mm, 1:L+1) *                     &
                                     CVT_HorizCors % assocLegPoly(1:L+1,index) * &
                                     CVT_HorizCors % GaussianWts(1:L+1))
      END DO
    END IF
  END DO
END DO

END SUBROUTINE Spherical_inv



!=================================================================================================
SUBROUTINE Spherical_adj (SpecCoeffs,        &
                          Field)

! Subroutine to perform adjoint spherical transform
! Note:
! The number of longitudes is 2L+1 in shtools
! The number of latitudes is  L+1 in shtools

USE onedvar_data, ONLY: L,                       &
                        CVT_HorizCors,           &
                        fft_worklen

IMPLICIT NONE

! Subroutine parameters
REAL, INTENT(OUT) :: SpecCoeffs(1:2, 0:L, 0:L)  ! cos, sin parts; l; m
REAL, INTENT(IN)  :: Field(1:L+1,1:2*L+1)

! Local variables
INTEGER           :: cs, mm, ll, lat, index, index_cos, index_sin
INTEGER           :: nlon, ierr
REAL              :: Inter1(1:2, 0:L, 1:L+1)    ! cos, sin parts; m; lat
REAL              :: Inter2(0:2*L)              ! inteface to FFT
REAL              :: rnlon, half_nlon

nlon           = 2*L+1
rnlon          = REAL(nlon)
half_nlon      = rnlon / 2.0


! Do the Fourier transform
DO lat = 1, L+1
  ! Put data into temporary structure
  Inter2(0:2*L) = Field(lat,1:nlon)
  ! Call fft
  CALL rfft1f (nlon, 1,                                           &
               Inter2(0:2*L), nlon,                               &
               CVT_HorizCors % fft_wsave,                         &
               fft_worklen,                                       &
               CVT_HorizCors % fft_work,                          &
               nlon, ierr)

  ! Place in the correct place for the Legendre transform and undo scalings done by rfft1f
  Inter1(1,0,lat) = Inter2(0) * rnlon                 ! mm = 0 term (cos only)
  ! The remaining mm
  DO mm = 1, L
    index_cos         = 2*mm-1
    index_sin         = 2*mm
    Inter1(1,mm,lat)  = Inter2(index_cos) * half_nlon ! cos term
    Inter1(2,mm,lat)  = Inter2(index_sin) * half_nlon ! sin term
  END DO
END DO


! Do the Legendre transform
DO cs = 1, 2
  DO mm = 0, L
    IF ((mm > 0) .OR. (cs == 1)) THEN          ! no sin term for mm=0
      DO ll = mm, L
        index                  = CVT_HorizCors % Plm_index(ll,mm)
        SpecCoeffs(cs, ll, mm) = SUM(Inter1(cs, mm, 1:L+1)         * &
                                     CVT_HorizCors % assocLegPoly(1:L+1,index))
      END DO
    END IF
  END DO
END DO

END SUBROUTINE Spherical_adj



!=================================================================================================
SUBROUTINE Spherical_inv_adj (SpecCoeffs,        &
                              Field)

! Subroutine to perform adjoint of inverse spherical transform
! Note:
! The number of longitudes is 2L+1 in shtools
! The number of latitudes is  L+1 in shtools

USE onedvar_data, ONLY: L,                       &
                        CVT_HorizCors,           &
                        fft_worklen

IMPLICIT NONE

! Subroutine parameters
REAL, INTENT(IN)  :: SpecCoeffs(1:2, 0:L, 0:L)  ! cos, sin parts; l; m
REAL, INTENT(OUT) :: Field(1:L+1,1:2*L+1)

! Local variables
INTEGER           :: cs, mm, ll, lat, index, index_cos, index_sin, ierr
INTEGER           :: nlon
REAL              :: Inter1(1:2, 0:L, 1:L+1)    ! cos, sin parts; m; lat
REAL              :: Inter2(0:2*L)              ! inteface to FFT
REAL              :: recip_nlon, total, fac

nlon           = 2*L+1
recip_nlon     = 1.0 / REAL(nlon)

! Do the Legendre transform
DO cs = 1, 2
  DO mm = 0, L
    IF ((mm > 0) .OR. (cs == 1)) THEN            ! no sin term for mm=0
      IF (mm == 0) THEN
        fac  = 0.5
      ELSE
        fac  = 0.25
      END IF
      DO lat = 1, L+1
        total = 0.0
        DO ll = mm, L
          index = CVT_HorizCors % Plm_index(ll,mm)
          total = total + SpecCoeffs(cs, ll, mm) * CVT_HorizCors % assocLegPoly(lat,index)
        END DO
        Inter1(cs, mm, lat) = total * CVT_HorizCors % GaussianWts(lat) * fac
      END DO
    END IF
  END DO
END DO



! Do the Fourier transform
DO lat = 1, L+1
  ! Put the data for this latitude in the right format for the FT
  ! This involves putting cos and sin terms in a 1d array and scaling
  ! to match the definition used by fftpack
  Inter2(0) = Inter1(1,0,lat)             ! mm = 0 term (cos only)
  ! The remaining mm
  DO mm = 1, L
    index_cos         = 2*mm-1
    index_sin         = 2*mm
    Inter2(index_cos) = 2.0 * Inter1(1,mm,lat)  ! cos term
    Inter2(index_sin) = 2.0 * Inter1(2,mm,lat)  ! sin term
  END DO

  ! Do the FT back to real space
  CALL rfft1b (nlon, 1,                                           &
               Inter2(0:2*L), nlon,                               &
               CVT_HorizCors % fft_wsave,                         &
               fft_worklen,                                       &
               CVT_HorizCors % fft_work,                          &
               nlon, ierr)

  ! Put data into output array and do further scaling
  Field(lat,1:nlon) = Inter2(0:2*L) * recip_nlon
END DO

END SUBROUTINE Spherical_inv_adj



!=================================================================================================
SUBROUTINE Physical_multid2linear ( state_flux,     & !IN  Multi-dimensional flux field
                                    state_tracer,   & !IN  Multi-dimensional tracer field
                                    zy )              !OUT Linear state
! Description:
! Rearranges flux and tracer information in physical space on multi-dimensional grids to
! a single 1D array as used by the rest of the INVICAT system.
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------

! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      09/2019  Initial version (R.N. Bannister)
!
! Code Description:
!   Language:           Fortran 90.


USE main_data,  ONLY  : ylon, ylat, ylev, nmonth, dim_1d, nfields


IMPLICIT NONE

! Subroutine parameters
REAL,    INTENT(IN)  :: state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL,    INTENT(IN)  :: state_tracer(1:ylat, 1:ylon, 1:ylev)
REAL,    INTENT(OUT) :: zy(1:dim_1d)


! Local variables
INTEGER              :: f, i, k, kk, m, ll, tot
INTEGER              :: nmonth_ylat, nmonth_ylat_ylon, ylon_ylat_nmonth_nfields, ylev_ylat
INTEGER              :: offset1, offset2, offset3


! Set some constants
! ------------------
nmonth_ylat              = nmonth * ylat
nmonth_ylat_ylon         = nmonth_ylat * ylon
ylon_ylat_nmonth_nfields = nmonth_ylat_ylon * nfields
ylev_ylat                = ylev * ylat


! Mirrors structure of 1d array fl_lin in main_1dvar (in main.f90) (no need to have complications due to duplicated months)
! Have query about the use of 'k' and 'kk' in the tracer and flux transfers below (respectively)

! Transfer the flux field to the correct location in the 1D array
DO f = 1, nfields
  offset3 = nmonth_ylat_ylon * (f-1)
  DO i = 1, ylon
    offset2 = nmonth_ylat * (i-1)
    DO k = 1, ylat
      kk = ylat - k + 1
      offset1 = nmonth * (kk-1)
      DO m = 1, nmonth
	tot     = m + offset1 + offset2 + offset3
	zy(tot) = state_flux(kk,i,m,f)
	!IF (fl_lin(tot).GT.ini_max) THEN
	!  ini_max = fl_lin(tot)
	!END IF
	!IF (fl_lin(tot) .LT. ini_min .AND. fl_lin(tot) .NE. 0.) THEN
	!  ini_min = fl_lin(tot)
	!END IF
      END DO
    END DO
  END DO
END DO

! Transfer the tracer field to the correct location in the 1D array
DO i = 1, ylon
  offset2 = ylev_ylat * (i-1)
  DO k = 1, ylat
    offset1 = ylev * (k-1)
    DO ll = 1, ylev
      tot     = ylon_ylat_nmonth_nfields + ll + offset1 + offset2
      zy(tot) = state_tracer(k,i,ll)
    END DO
  END DO
END DO

END SUBROUTINE Physical_multid2linear




!=================================================================================================
SUBROUTINE Physical_linear2multid ( state_flux,     & !OUT Multi-dimensional flux field
                                    state_tracer,   & !OUT Multi-dimensional tracer field
                                    zy )              !IN  Linear state
! Description:
! Rearranges flux and tracer information in physical space on a single 1D array as used by the
! rest of the INVICAT system to multi-dimensional grids
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------

! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      09/2019  Initial version (R.N. Bannister)
!
! Code Description:
!   Language:           Fortran 90.


USE main_data,  ONLY  : ylon, ylat, ylev, nmonth, dim_1d, nfields


IMPLICIT NONE

! Subroutine parameters
REAL,    INTENT(OUT) :: state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL,    INTENT(OUT) :: state_tracer(1:ylat, 1:ylon, 1:ylev)
REAL,    INTENT(IN)  :: zy(1:dim_1d)


! Local variables
INTEGER              :: f, i, k, kk, m, ll, tot
INTEGER              :: nmonth_ylat, nmonth_ylat_ylon, ylon_ylat_nmonth_nfields, ylev_ylat
INTEGER              :: offset1, offset2, offset3


! Set some constants
! ------------------
nmonth_ylat              = nmonth * ylat
nmonth_ylat_ylon         = nmonth_ylat * ylon
ylon_ylat_nmonth_nfields = nmonth_ylat_ylon * nfields
ylev_ylat                = ylev * ylat


! Mirrors structure of 1d array fl_lin in main_1dvar (in main.f90) (no need to have complications due to duplicated months)
! Have query about the use of 'k' and 'kk' in the tracer and flux transfers below (respectively)

! Transfer the flux field to the correct location in the 1D array
DO f = 1, nfields
  offset3 = nmonth_ylat_ylon * (f-1)
  DO i = 1, ylon
    offset2 = nmonth_ylat * (i-1)
    DO k = 1, ylat
      kk = ylat - k + 1
      offset1 = nmonth * (kk-1)
      DO m = 1, nmonth
	tot                  = m + offset1 + offset2 + offset3
	state_flux(kk,i,m,f) = zy(tot)
	!IF (fl_lin(tot).GT.ini_max) THEN
	!  ini_max = fl_lin(tot)
	!END IF
	!IF (fl_lin(tot) .LT. ini_min .AND. fl_lin(tot) .NE. 0.) THEN
	!  ini_min = fl_lin(tot)
	!END IF
      END DO
    END DO
  END DO
END DO

! Transfer the tracer field to the correct location in the 1D array
DO i = 1, ylon
  offset2 = ylev_ylat * (i-1)
  DO k = 1, ylat
    offset1 = ylev * (k-1)
    DO ll = 1, ylev
      tot                  = ylon_ylat_nmonth_nfields + ll + offset1 + offset2
      state_tracer(k,i,ll) = zy(tot)
    END DO
  END DO
END DO

END SUBROUTINE Physical_linear2multid



!=================================================================================================
SUBROUTINE ChangeHorizGrid (nother,      & ! size of last index (e.g. months or levels)
                            nlongs_in,   &
                            nlats_in,    &
                            data_in,     &
                            longs_in,    &
                            lats_in,     &
                            nlongs_out,  &
                            nlats_out,   &
                            data_out,    &
                            longs_out,   &
                            lats_out,    &
                            CVT_gridchange_details)

! Description:
! Performs interpolation from the input grid to the output grid (no change in vertical grid)
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs bilinear interpolation in the horizontal

! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!
! Code Description:
!   Language:           Fortran 90.


USE onedvar_data, ONLY : ChangeGrid_type

IMPLICIT NONE

! Subroutine parameters
INTEGER, INTENT(IN)  :: nother
INTEGER, INTENT(IN)  :: nlongs_in, nlats_in
REAL,    INTENT(IN)  :: data_in(1:nlats_in, 1:nlongs_in, 1:nother)
REAL,    INTENT(IN)  :: longs_in(1:nlongs_in)
REAL,    INTENT(IN)  :: lats_in(1:nlats_in)
INTEGER, INTENT(IN)  :: nlongs_out, nlats_out
REAL,    INTENT(OUT) :: data_out(1:nlats_out, 1:nlongs_out, 1:nother)
REAL,    INTENT(IN)  :: longs_out(1:nlongs_out)
REAL,    INTENT(IN)  :: lats_out(1:nlats_out)
TYPE (ChangeGrid_type), INTENT(INOUT) :: CVT_gridchange_details

! Local variables
INTEGER              :: i, j
INTEGER              :: index1, index2
REAL                 :: w11, w12, w21, w22
REAL                 :: data_in_extra(0:nlats_in+1, 0:nlongs_in+1, 1:nother) ! With Halos
REAL                 :: longs_in_extra(0:nlongs_in+1)
REAL                 :: lats_in_extra(0:nlats_in+1)

! Function used
INTEGER              :: FindLowerIndex

! Copy information from the data_in array to the new array that will contain
data_in_extra(1:nlats_in, 1:nlongs_in, 1:nother) = data_in(1:nlats_in, 1:nlongs_in, 1:nother)

! Swap information around the borders
CALL Swap_halos (nother, nlongs_in, nlats_in, &
                 data_in_extra(0:nlats_in+1, 0:nlongs_in+1, 1:nother))


! Compute the index arrays if not known already
IF (.NOT.ALLOCATED(CVT_gridchange_details % grid_change_long_lower)) THEN

  ! These specify the points on the SHtools grid are either side of each point on the TOMCAT grid
  ! and the linear interpolation weights

  ! Set the longitudes and latitudes
  longs_in_extra(1:nlongs_in) = longs_in(1:nlongs_in)
  longs_in_extra(0)           = longs_in(1) - (longs_in(2) - longs_in(1))
  longs_in_extra(nlongs_in+1) = longs_in(nlongs_in) + (longs_in(nlongs_in) - longs_in(nlongs_in-1))
  lats_in_extra(1:nlats_in)   = lats_in(1:nlats_in)
  lats_in_extra(0)            = lats_in(1) - (lats_in(2) - lats_in(1))
  lats_in_extra(nlats_in+1)   = lats_in(nlats_in) + (lats_in(nlats_in) - lats_in(nlats_in-1))
  ! PRINT *, '-----------------'
  ! PRINT *, 'CVT_SHtools_longs (extra)'
  ! PRINT *, '-----------------'
  ! PRINT *, longs_in_extra(0:nlongs_in+1)
  ! PRINT *, '-----------------'
  ! PRINT *, 'CVT_SHtools_lats (extra)'
  ! PRINT *, '-----------------'
  ! PRINT *, lats_in_extra(0:nlats_in+1)
  ! PRINT *, '-----------------'
  ! PRINT *, '-----------------'

  ! Set-up indices for longitude
  ALLOCATE(CVT_gridchange_details % grid_change_long_lower(1:nlongs_out))
  ALLOCATE(CVT_gridchange_details % grid_change_long_lower_wt(1:nlongs_out))
  ALLOCATE(CVT_gridchange_details % grid_change_long_upper_wt(1:nlongs_out))
  !PRINT *, 'Constructing reference data for longitudes'
  !PRINT *, '------------------------------------------'
  DO i = 1, nlongs_out
    index1                                          = FindLowerIndex (longs_out(i), nlongs_in, longs_in_extra(0:nlongs_in+1))
    index2                                          = index1 + 1
    CVT_gridchange_details % grid_change_long_lower(i)    = index1
    !PRINT *, i, ' longitude', longs_out(i), ' on o/p grid is between', index1, ' and', index2
    !PRINT *, '      (', longs_in_extra(index1), longs_in_extra(index2), ')'
    CVT_gridchange_details % grid_change_long_lower_wt(i) = (longs_in_extra(index2) - longs_out(i)) / &
                                                      (longs_in_extra(index2) - longs_in_extra(index1))
    CVT_gridchange_details % grid_change_long_upper_wt(i) = 1. - CVT_gridchange_details % grid_change_long_lower_wt(i)
    !PRINT *, '         with weights (', CVT_gridchange_details % grid_change_long_lower_wt(i), &
    !                                    CVT_gridchange_details % grid_change_long_upper_wt(i), ')'
  END DO


  ! Set-up indices for latitude
  ALLOCATE(CVT_gridchange_details % grid_change_lat_lower(1:nlats_out))
  ALLOCATE(CVT_gridchange_details % grid_change_lat_lower_wt(1:nlats_out))
  ALLOCATE(CVT_gridchange_details % grid_change_lat_upper_wt(1:nlats_out))
  !PRINT *, 'Constructing reference data for latitudes'
  !PRINT *, '------------------------------------------'
  DO j = 1, nlats_out
    index1                                          = FindLowerIndex (lats_out(j), nlats_in, lats_in_extra(0:nlats_in+1))
    index2                                          = index1 + 1
    CVT_gridchange_details % grid_change_lat_lower(j)     = index1
    !PRINT *, j, 'latitude', lats_out(j), ' on o/p grid is between', index1, ' and', index2
    !PRINT *, '         (', lats_in_extra(index1), lats_in_extra(index2), ')'
    CVT_gridchange_details % grid_change_lat_lower_wt(j)  = (lats_in_extra(index2) - lats_out(j)) / &
                                                      (lats_in_extra(index2) - lats_in_extra(index1))
    CVT_gridchange_details % grid_change_lat_upper_wt(j)  = 1. - CVT_gridchange_details % grid_change_lat_lower_wt(j)
    !PRINT *, '         with weights (', CVT_gridchange_details % grid_change_lat_lower_wt(j), &
    !                                    CVT_gridchange_details % grid_change_lat_upper_wt(j), ')'
  END DO
END IF


! Interpolate the data from the input grid to the output grid
DO i = 1, nlats_out
  index1 = CVT_gridchange_details % grid_change_lat_lower(i)
  DO j = 1, nlongs_out
    w11    = CVT_gridchange_details % grid_change_lat_lower_wt(i) * CVT_gridchange_details % grid_change_long_lower_wt(j)
    w12    = CVT_gridchange_details % grid_change_lat_lower_wt(i) * CVT_gridchange_details % grid_change_long_upper_wt(j)
    w21    = CVT_gridchange_details % grid_change_lat_upper_wt(i) * CVT_gridchange_details % grid_change_long_lower_wt(j)
    w22    = CVT_gridchange_details % grid_change_lat_upper_wt(i) * CVT_gridchange_details % grid_change_long_upper_wt(j)
    index2 = CVT_gridchange_details % grid_change_long_lower(j)

    data_out(i,j,1:nother) = w11 * data_in_extra(index1,  index2,  1:nother) + &
                             w12 * data_in_extra(index1,  index2+1,1:nother) + &
                             w21 * data_in_extra(index1+1,index2,  1:nother) + &
                             w22 * data_in_extra(index1+1,index2+1,1:nother)
  END DO
END DO

END SUBROUTINE ChangeHorizGrid





!=================================================================================================
SUBROUTINE ChangeHorizGrid_adj (nother,       & ! size of last index (e.g. months or levels)
                                nlongs_in,    &
                                nlats_in,     &
                                data_in_hat,  &
                                longs_in,     &
                                lats_in,      &
                                nlongs_out,   &
                                nlats_out,    &
                                data_out_hat, &
                                longs_out,    &
                                lats_out,     &
                                CVT_gridchange_details)

! Description: (adjoint of the following)
! Performs interpolation from the input grid to the output grid (no change in vertical grid)
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs bilinear interpolation in the horizontal

! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!
! Code Description:
!   Language:           Fortran 90.


USE onedvar_data, ONLY : ChangeGrid_type

IMPLICIT NONE

! Subroutine parameters
INTEGER, INTENT(IN)  :: nother
INTEGER, INTENT(IN)  :: nlongs_in, nlats_in
REAL,    INTENT(OUT) :: data_in_hat(1:nlats_in, 1:nlongs_in, 1:nother)
REAL,    INTENT(IN)  :: longs_in(1:nlongs_in)
REAL,    INTENT(IN)  :: lats_in(1:nlats_in)
INTEGER, INTENT(IN)  :: nlongs_out, nlats_out
REAL,    INTENT(IN)  :: data_out_hat(1:nlats_out, 1:nlongs_out, 1:nother)
REAL,    INTENT(IN)  :: longs_out(1:nlongs_out)
REAL,    INTENT(IN)  :: lats_out(1:nlats_out)
TYPE (ChangeGrid_type), INTENT(INOUT) :: CVT_gridchange_details

! Local variables
INTEGER              :: i, j
INTEGER              :: index1, index2
REAL                 :: w11, w12, w21, w22
REAL                 :: data_in_extra_hat(0:nlats_in+1, 0:nlongs_in+1, 1:nother) ! With Halos
REAL                 :: longs_in_extra(0:nlongs_in+1)
REAL                 :: lats_in_extra(0:nlats_in+1)

! Function used
INTEGER              :: FindLowerIndex



! Compute the index arrays if not known already
IF (.NOT.ALLOCATED(CVT_gridchange_details % grid_change_long_lower)) THEN

  ! These specify the points on the SHtools grid are either side of each point on the TOMCAT grid
  ! and the linear interpolation weights

  ! Set the longitudes and latitudes
  longs_in_extra(1:nlongs_in) = longs_in(1:nlongs_in)
  longs_in_extra(0)           = longs_in(1) - (longs_in(2) - longs_in(1))
  longs_in_extra(nlongs_in+1) = longs_in(nlongs_in) + (longs_in(nlongs_in) - longs_in(nlongs_in-1))
  lats_in_extra(1:nlats_in)   = lats_in(1:nlats_in)
  lats_in_extra(0)            = lats_in(1) - (lats_in(2) - lats_in(1))
  lats_in_extra(nlats_in+1)   = lats_in(nlats_in) + (lats_in(nlats_in) - lats_in(nlats_in-1))

  ! Set-up indices for longitude
  ALLOCATE(CVT_gridchange_details % grid_change_long_lower(1:nlongs_out))
  ALLOCATE(CVT_gridchange_details % grid_change_long_lower_wt(1:nlongs_out))
  ALLOCATE(CVT_gridchange_details % grid_change_long_upper_wt(1:nlongs_out))
  DO i = 1, nlongs_out
    index1                           = FindLowerIndex (longs_out(i), nlongs_in, longs_in_extra(0:nlongs_in+1))
    index2                           = index1 + 1
    CVT_gridchange_details % grid_change_long_lower(i)    = index1
    CVT_gridchange_details % grid_change_long_lower_wt(i) = (longs_in(index2) - longs_out(i)) / &
                                                            (longs_in(index2) - longs_in(index1))
    CVT_gridchange_details % grid_change_long_upper_wt(i) = 1. - CVT_gridchange_details % grid_change_long_lower_wt(i)
  END DO


  ! Set-up indices for latitude
  ALLOCATE(CVT_gridchange_details % grid_change_lat_lower(1:nlats_out))
  ALLOCATE(CVT_gridchange_details % grid_change_lat_lower_wt(1:nlats_out))
  ALLOCATE(CVT_gridchange_details % grid_change_lat_upper_wt(1:nlats_out))
  DO j = 1, nlats_out
    index1                           = FindLowerIndex (lats_out(j), nlats_in, lats_in_extra(0:nlats_in+1))
    index2                           = index1 + 1
    CVT_gridchange_details % grid_change_lat_lower(j)     = index1
    CVT_gridchange_details % grid_change_lat_lower_wt(j)  = (lats_in(index2) - lats_out(j)) / (lats_in(index2) - lats_in(index1))
    CVT_gridchange_details % grid_change_lat_upper_wt(j)  = 1. - CVT_gridchange_details % grid_change_lat_lower_wt(j)
  END DO
END IF

! Initialise to zero
data_in_extra_hat(0:nlats_in+1,  0:nlongs_in+1,  1:nother) = 0.0
! Interpolate the data from the input grid to the output grid
DO i = 1, nlats_out
  index1 = CVT_gridchange_details % grid_change_lat_lower(i)
  DO j = 1, nlongs_out
    w11    = CVT_gridchange_details % grid_change_lat_lower_wt(i) * CVT_gridchange_details % grid_change_long_lower_wt(j)
    w12    = CVT_gridchange_details % grid_change_lat_lower_wt(i) * CVT_gridchange_details % grid_change_long_upper_wt(j)
    w21    = CVT_gridchange_details % grid_change_lat_upper_wt(i) * CVT_gridchange_details % grid_change_long_lower_wt(j)
    w22    = CVT_gridchange_details % grid_change_lat_upper_wt(i) * CVT_gridchange_details % grid_change_long_upper_wt(j)
    index2 = CVT_gridchange_details % grid_change_long_lower(j)

!    data_out(i,j,1:nother) = w11 * data_in_extra(index1,  index2,  1:nother) + &
!                             w12 * data_in_extra(index1,  index2+1,1:nother) + &
!                             w21 * data_in_extra(index1+1,index2,  1:nother) + &
!                             w22 * data_in_extra(index1+1,index2+1,1:nother)

    data_in_extra_hat(index1,  index2,  1:nother) = data_in_extra_hat(index1,  index2,  1:nother) + &
                                                    w11 * data_out_hat(i,j,1:nother)
    data_in_extra_hat(index1,  index2+1,1:nother) = data_in_extra_hat(index1,  index2+1,1:nother) + &
                                                    w12 * data_out_hat(i,j,1:nother)
    data_in_extra_hat(index1+1,index2,  1:nother) = data_in_extra_hat(index1+1,index2,  1:nother) + &
                                                    w21 * data_out_hat(i,j,1:nother)
    data_in_extra_hat(index1+1,index2+1,1:nother) = data_in_extra_hat(index1+1,index2+1,1:nother) + &
                                                    w22 * data_out_hat(i,j,1:nother)

  END DO
END DO


! Swap information around the borders
CALL Swap_halos_adj (nother, nlongs_in, nlats_in, &
                     data_in_extra_hat(0:nlats_in+1, 0:nlongs_in+1, 1:nother))


! Copy information from the data_in array to the new array that will contain
data_in_hat(1:nlats_in, 1:nlongs_in, 1:nother) = data_in_extra_hat(1:nlats_in, 1:nlongs_in, 1:nother)

END SUBROUTINE ChangeHorizGrid_adj






!=================================================================================================
INTEGER FUNCTION FindLowerIndex (value, arraysize, array)

IMPLICIT NONE

REAL,     INTENT(IN) :: value
INTEGER,  INTENT(IN) :: arraysize
REAL,     INTENT(IN) :: array(0:arraysize+1)

INTEGER              :: testpoint, result

testpoint = -1
DO
  testpoint = testpoint + 1
  IF (testpoint > arraysize) THEN
    result = arraysize + 1
    EXIT
  ELSE
    ! Test to see if the value is between testpoint and testpoint+1
    ! (should work for ascending and descending arrays)
    IF (((value >= array(testpoint)) .AND. (value < array(testpoint+1))) .OR.  &
        ((value >= array(testpoint+1)) .AND. (value < array(testpoint)))) THEN
      result = testpoint
      EXIT
    END IF
  END IF
END DO

FindLowerIndex = result

END FUNCTION FindLowerIndex



!=================================================================================================
SUBROUTINE Swap_halos (nother, nlongs, nlats, datafield)

IMPLICIT NONE

! Subroutine parameters
INTEGER, INTENT(IN)    :: nother, nlongs, nlats
REAL,    INTENT(INOUT) :: datafield(0:nlats+1, 0:nlongs+1, 1:nother)

! Local variables
INTEGER                :: i, adjacent

! Swap halo elements in longitude
datafield(1:nlats, 0,        1:nother) = datafield(1:nlats, nlongs, 1:nother)
datafield(1:nlats, nlongs+1, 1:nother) = datafield(1:nlats, 1,      1:nother)

! Swap halo elements in latitude
DO i = 1, nlongs
  adjacent = i + nlongs/2
  IF (adjacent > nlongs) adjacent = adjacent - nlongs
  datafield(0,       i, 1:nother) = datafield(1,     adjacent, 1:nother)
  datafield(nlats+1, i, 1:nother) = datafield(nlats, adjacent, 1:nother)
END DO

! Corners
datafield(0,       0,        1:nother) = datafield(0,       nlongs, 1:nother)
datafield(0,       nlongs+1, 1:nother) = datafield(0,       1,      1:nother)
datafield(nlats+1, 0,        1:nother) = datafield(nlats+1, nlongs, 1:nother)
datafield(nlats+1, nlongs+1, 1:nother) = datafield(nlats+1, 1,      1:nother)

END SUBROUTINE Swap_halos






!=================================================================================================
SUBROUTINE Swap_halos_adj (nother, nlongs, nlats, datafield_hat)

IMPLICIT NONE

! Subroutine parameters
INTEGER, INTENT(IN)    :: nother, nlongs, nlats
REAL,    INTENT(INOUT) :: datafield_hat(0:nlats+1, 0:nlongs+1, 1:nother)

! Local variables
INTEGER                :: i, adjacent

! Corners
! datafield(0,       0,        1:nother) = datafield(0,       nlongs, 1:nother)
! datafield(0,       nlongs+1, 1:nother) = datafield(0,       1,      1:nother)
! datafield(nlats+1, 0,        1:nother) = datafield(nlats+1, nlongs, 1:nother)
! datafield(nlats+1, nlongs+1, 1:nother) = datafield(nlats+1, 1,      1:nother)

datafield_hat(0,       nlongs,  1:nother) = datafield_hat(0,       nlongs, 1:nother) + datafield_hat(0,       0,        1:nother)
datafield_hat(0,       1,       1:nother) = datafield_hat(0,       1,      1:nother) + datafield_hat(0,       nlongs+1, 1:nother)
datafield_hat(nlats+1, nlongs,  1:nother) = datafield_hat(nlats+1, nlongs, 1:nother) + datafield_hat(nlats+1, 0,        1:nother)
datafield_hat(nlats+1, 1,       1:nother) = datafield_hat(nlats+1, 1,      1:nother) + datafield_hat(nlats+1, nlongs+1, 1:nother)

! Swap halo elements in latitude
DO i = 1, nlongs
  adjacent = i + nlongs/2
  IF (adjacent > nlongs) adjacent = adjacent - nlongs
  ! datafield(0,       i, 1:nother) = datafield(1,     adjacent, 1:nother)
  ! datafield(nlats+1, i, 1:nother) = datafield(nlats, adjacent, 1:nother)
  datafield_hat(1,     adjacent, 1:nother) = datafield_hat(1,     adjacent, 1:nother) + datafield_hat(0,       i, 1:nother)
  datafield_hat(nlats, adjacent, 1:nother) = datafield_hat(nlats, adjacent, 1:nother) + datafield_hat(nlats+1, i, 1:nother)
END DO

! Swap halo elements in longitude
datafield_hat(1:nlats, nlongs, 1:nother) = datafield_hat(1:nlats, nlongs, 1:nother) + datafield_hat(1:nlats, 0,        1:nother)
datafield_hat(1:nlats, 1,      1:nother) = datafield_hat(1:nlats, 1,      1:nother) + datafield_hat(1:nlats, nlongs+1, 1:nother)

! Set halos to zero
!!datafield_hat(0,         0:nlongs+1, 1:nother) = 0.0
!!datafield_hat(nlats+1,   0:nlongs+1, 1:nother) = 0.0
!!datafield_hat(0:nlats+1, 0,          1:nother) = 0.0
!!datafield_hat(0:nlats+1, nlongs+1,   1:nother) = 0.0

END SUBROUTINE Swap_halos_adj


! =======================================================================================
REAL FUNCTION InnerProdReal (lats,                    &
                             lons,                    &
                             FieldA,                  &
                             FieldB,                  &
                             GaussWts)
! To compute inner product in real space

IMPLICIT NONE
! Declare the subroutine parameters
INTEGER,      INTENT(IN)  :: lats                     !Number of latitudes
INTEGER,      INTENT(IN)  :: lons                     !Number of longitudes
REAL,         INTENT(IN)  :: FieldA(1:lats, 1:lons)
REAL,         INTENT(IN)  :: FieldB(1:lats, 1:lons)
REAL,         INTENT(IN)  :: GaussWts(1:lats)

! Declare the local variables
REAL                      :: ip
REAL                      :: Product(1:lats, 1:lons)
REAL                      :: Integrate


Product(1:lats, 1:lons) = FieldA(1:lats, 1:lons) * FieldB(1:lats, 1:lons)
ip                      = Integrate (lats,                    &
                                     lons,                    &
                                     Product(1:lats, 1:lons), &
                                     GaussWts(1:lats))
InnerProdReal           = ip

END FUNCTION InnerProdReal



! =======================================================================================
REAL FUNCTION InnerProdSpec (L,                       &
                             FieldA,                  &
                             FieldB)
! To compute inner product in real space

IMPLICIT NONE
! Declare the subroutine parameters
INTEGER,      INTENT(IN)  :: L                     !Max wavenumber
REAL,         INTENT(IN)  :: FieldA(1:2,0:L,0:L)
REAL,         INTENT(IN)  :: FieldB(1:2,0:L,0:L)

! Declare the local variables
REAL                      :: ip
INTEGER                   :: ll, mm

ip = 0.0
DO ll = 0, L
  ip = ip + FieldA(1,ll,0) * FieldB(1,ll,0)  ! Cosine for mm=0
  IF (ll > 0) THEN
    DO mm = 1, ll
      ip = ip + FieldA(1,ll,mm) * FieldB(1,ll,mm) & ! Cosine contribution
              + FieldA(2,ll,mm) * FieldB(2,ll,mm)   ! Sine contribution
    END DO
  END IF
END DO

!ip = SUM(FieldA(1:2,0:L,0:L) * FieldB(1:2,0:L,0:L))

InnerProdSpec = ip

END FUNCTION InnerProdSpec



! =======================================================================================
REAL FUNCTION Integrate (lats,                  &
                         lons,                  &
                         Field,                 &
                         GaussWts)
! To do global integral in real space

USE cparam,       ONLY: pi

IMPLICIT NONE
! Declare the subroutine parameters
INTEGER,      INTENT(IN)  :: lats
INTEGER,      INTENT(IN)  :: lons
REAL,         INTENT(IN)  :: Field(1:lats, 1:lons)
REAL,         INTENT(IN)  :: GaussWts(1:lats)

! Declare the local variables
REAL                      :: grandtotal, total
INTEGER                   :: i, j

grandtotal = 0.0
DO j = 1, lats
  total = SUM(Field(j,1:lons))
  grandtotal = grandtotal + total * GaussWts(j)
END DO

Integrate = grandtotal * 2.0 * pi / REAL(lons)

END FUNCTION Integrate



!=================================================================================================
SUBROUTINE cvt_v (chi_vh_tracer, & ! IN  Version of chi (for tracer) between vert and horiz transforms
                  state_tracer )   ! OUT Tracer field in real space (lat, long, lev)

! Description:
! Performs the vertical control variable transform, and multiplies by the error std devs
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs the vertical transform and multiplies by the std dev of the bg error
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           11/2019  Allow for adjustment for correlation model.
!           01/2020  Use structures.
!           01/2020  Add switch to turn off adjustment
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev

USE onedvar_data, ONLY: CVT_VertCors, &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(IN)  :: chi_vh_tracer(1:ylat, 1:ylon, 1:ylev)
REAL, INTENT(OUT) :: state_tracer(1:ylat, 1:ylon, 1:ylev)

! Local variables:
INTEGER           :: i, j, k, mode
REAL              :: tracer_modes(1:ylev)
REAL              :: tracer_prescale(1:ylev)

DO j = 1, ylat
  DO i = 1, ylon
    ! Project on to vertical modes
    DO mode = 1, ylev
      tracer_modes(mode) = SUM(CVT_VertCors % vert_eigenvec_tracer(1:ylev, mode) * &
                               chi_vh_tracer(j, i, 1:ylev))
    END DO
    ! Multiply by standard deviation (a function of lat and lon)
    tracer_modes(1:ylev) = tracer_modes(1:ylev) * CVT_VertCors % std_vspec_tracer(1:ylev,j,i)
    ! Project back to real space
    DO k = 1, ylev
      tracer_prescale(k) = SUM(CVT_VertCors % vert_eigenvec_tracer(k,1:ylev) *     &
                               tracer_modes(1:ylev))
    END DO
    IF (CVT_VertCors % ForceVCorr) THEN
      ! Do the adjustment and multiply by the spatially-dependent standard deviation
      state_tracer(j,i,1:ylev) = tracer_prescale(1:ylev) *                           &
                                 CVT_std % std_tracer(j,i,1:ylev) /                  &
                                 CVT_VertCors % vert_adjust_tracer(1:ylev,j,i)
    ELSE
      ! Do not do the adjustment
      state_tracer(j,i,1:ylev) = tracer_prescale(1:ylev) *                           &
                                 CVT_std % std_tracer(j,i,1:ylev)
    END IF
  END DO
END DO

END SUBROUTINE cvt_v



!=================================================================================================
SUBROUTINE cvt_v_adj (chi_vh_tracer_hat, & ! OUT Version of chi (for tracer) between vert and horiz transforms
                      state_tracer_hat )   ! IN  Tracer field in real space (lat, long, lev)

! Description:  (adjoint of the following)
! Performs the vertical control variable transform, and multiplies by the error std devs
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs the vertical transform and multiplies by the std dev of the bg error
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           11/2019  Allow for adjustment for correlation model.
!           01/2020  Use structures.
!           01/2020  Add switch to turn off adjustment
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev

USE onedvar_data, ONLY: CVT_VertCors, &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(OUT) :: chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev)
REAL, INTENT(IN)  :: state_tracer_hat(1:ylat, 1:ylon, 1:ylev)

! Local variables:
INTEGER           :: i, j, k, mode
REAL              :: tracer_modes(1:ylev)
REAL              :: tracer_prescale(1:ylev)


DO j = 1, ylat
  DO i = 1, ylon
    IF (CVT_VertCors % ForceVCorr) THEN
      ! Do the adjustment and multiply by the spatially-dependent standard deviation
      tracer_prescale(1:ylev) = state_tracer_hat(j, i, 1:ylev) *                     &
                                CVT_std % std_tracer(j,i,1:ylev) /                   &
                                CVT_VertCors % vert_adjust_tracer(1:ylev,j,i)
    ELSE
      ! Do not do the adjustment
      tracer_prescale(1:ylev) = state_tracer_hat(j, i, 1:ylev) *                     &
                                CVT_std % std_tracer(j,i,1:ylev)
    END IF

    ! Project on to vertical modes
    DO mode = 1, ylev
      tracer_modes(mode) = SUM(CVT_VertCors % vert_eigenvec_tracer(1:ylev, mode) *   &
                               tracer_prescale(1:ylev))
    END DO
    ! Multiply by standard deviation (a function of lat)
    tracer_modes(1:ylev) = tracer_modes(1:ylev) * CVT_VertCors % std_vspec_tracer(1:ylev,j,i)
    ! Project back to real space
    DO k = 1, ylev
      chi_vh_tracer_hat(j,i,k) = SUM(CVT_VertCors % vert_eigenvec_tracer(k,1:ylev) * &
                                     tracer_modes(1:ylev))
    END DO
  END DO
END DO

END SUBROUTINE cvt_v_adj




!=================================================================================================
SUBROUTINE cvt_v_inv (chi_vh_tracer, & ! OUT Version of chi (for tracer) between vert and horiz transforms
                      state_tracer )   ! IN  Tracer field in real space (lat, long, lev)

! Description:
! Divides by the std dev of the bg error and performs the inverse vertical control variable transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs the inverse vertical transform and multiplies by the std dev of the bg error
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           11/2019  Allow for adjustment for correlation model.
!           01/2020  Use structures.
!           01/2020  Add switch to turn off adjustment
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev

USE onedvar_data, ONLY: CVT_VertCors, &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(OUT) :: chi_vh_tracer(1:ylat, 1:ylon, 1:ylev)
REAL, INTENT(IN)  :: state_tracer(1:ylat, 1:ylon, 1:ylev)

! Local variables:
INTEGER           :: i, j, k, mode
REAL              :: tracer_modes(1:ylev)
REAL              :: tracer_prescale(1:ylev)

DO j = 1, ylat
  DO i = 1, ylon
    IF (CVT_VertCors % ForceVCorr) THEN
      ! Divide by the spatially-dependent standard deviation and multiply by the adjustment
      tracer_prescale(1:ylev) = state_tracer(j, i, 1:ylev) *                       &
                                CVT_VertCors % vert_adjust_tracer(1:ylev,j,i) /    &
                                CVT_std % std_tracer(j, i, 1:ylev)
    ELSE
      ! Do not do the adjustment
      tracer_prescale(1:ylev) = state_tracer(j, i, 1:ylev)  /                      &
                                CVT_std % std_tracer(j, i, 1:ylev)
    END IF

    ! Project on to vertical modes
    DO mode = 1, ylev
      tracer_modes(mode) = SUM(CVT_VertCors % vert_eigenvec_tracer(1:ylev, mode) * &
                               tracer_prescale(1:ylev))
    END DO
    ! Divide by standard deviation (a function of lat)
    tracer_modes(1:ylev) = tracer_modes(1:ylev) / CVT_VertCors % std_vspec_tracer(1:ylev,j,i)
    ! Project back to real space
    DO k = 1, ylev
      chi_vh_tracer(j,i,k) = SUM(CVT_VertCors % vert_eigenvec_tracer(k,1:ylev) *   &
                                 tracer_modes(1:ylev))
    END DO
  END DO
END DO

END SUBROUTINE cvt_v_inv




!=================================================================================================
SUBROUTINE cvt_v_inv_adj (chi_vh_tracer_hat, & ! IN  Version of chi (for tracer) between vert and horiz transforms
                          state_tracer_hat )   ! OUT Tracer field in real space (lat, long, lev)

! Description:
! Performs the adjoint of the inverse of the vertical control variable transform, and divides by the error std devs
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs the vertical transform and multiplies by the std dev of the bg error
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      09/2019  Initial version (R.N. Bannister)
!           11/2019  Allow for adjustment for correlation model.
!           01/2020  Use structures.
!           01/2020  Add switch to turn off adjustment
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, ylev

USE onedvar_data, ONLY: CVT_VertCors, &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(IN)  :: chi_vh_tracer_hat(1:ylat, 1:ylon, 1:ylev)
REAL, INTENT(OUT) :: state_tracer_hat(1:ylat, 1:ylon, 1:ylev)

! Local variables:
INTEGER           :: i, j, k, mode
REAL              :: tracer_modes(1:ylev)
REAL              :: tracer_prescale(1:ylev)

DO j = 1, ylat
  DO i = 1, ylon
    ! Project on to vertical modes
    DO mode = 1, ylev
      tracer_modes(mode) = SUM(CVT_VertCors % vert_eigenvec_tracer(1:ylev, mode) *  &
                               chi_vh_tracer_hat(j, i, 1:ylev))
    END DO
    ! Divide by standard deviation (a function of lat)
    tracer_modes(1:ylev) = tracer_modes(1:ylev) / CVT_VertCors % std_vspec_tracer(1:ylev,j,i)
    ! Project back to real space
    DO k = 1, ylev
      tracer_prescale(k) = SUM(CVT_VertCors % vert_eigenvec_tracer(k,1:ylev) *      &
                               tracer_modes(1:ylev))
    END DO
    ! Multiply by the adjustment and divide by the spatially-dependent standard deviation
    IF (CVT_VertCors % ForceVCorr) THEN
      state_tracer_hat(j,i,1:ylev) = tracer_prescale(1:ylev) *                      &
                                     CVT_VertCors % vert_adjust_tracer(1:ylev,j,i) /&
                                     CVT_std % std_tracer(j, i, 1:ylev)
    ELSE
      state_tracer_hat(j,i,1:ylev) = tracer_prescale(1:ylev) /                      &
                                     CVT_std % std_tracer(j, i, 1:ylev)
    END IF
  END DO
END DO

END SUBROUTINE cvt_v_inv_adj





!=================================================================================================
SUBROUTINE cvt_t (chi_th_flux, & ! IN  Version of chi (for flux) between temp and horiz transforms
                  state_flux )   ! OUT Flux field in real space (lat, long, month)

! Description:
! Performs the temporal control variable transform, and multiplies by the error std devs
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs the temporal transform and multiplies by the std dev of the bg error
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Include nfields (R.N. Bannister)
!           01/2020  Use structures.
!           07/2021  Allow for different fields to have different correlation types
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, nmonth, nfields

USE onedvar_data, ONLY: CVT_TemporalCors,       &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(IN)  :: chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL, INTENT(OUT) :: state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)

! Local variables:
INTEGER           :: i, j, t, mode, f
REAL              :: flux_modes(1:nmonth)
REAL              :: flux_prescale(1:nmonth)



DO f = 1, nfields

  ! Decide on the temporal correlation option
  SELECT CASE (CVT_TemporalCors % temporal_covs(f))

  CASE (0)
    ! White noise in time (no temporal correlations)
    ! Expect input to forward routine to be a function of lat/long/time/field
    ! Output of the forward routine to be a function of lat/long/time/field
    state_flux(1:ylat, 1:ylon, 1:nmonth, f) = chi_th_flux(1:ylat, 1:ylon, 1:nmonth, f) * &
                                              CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, f)

  CASE (1)
    ! Temporal covariances, spatial scale a fn of timescale
    ! This is the non-symmetric temporal transform
    ! Expect input to forward routine to be a function of lat/long/timescale/field
    ! Output of the forward routine to be a function of lat/long/time/field
    DO i = 1, ylat
      DO j = 1, ylon
        ! Multiply by standard deviation (flux_modes a function of timescale)
        flux_modes(1:nmonth) = chi_th_flux(i, j, 1:nmonth, f) *                        &
                               CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
        ! Project back to time
        DO t = 1, nmonth
          flux_prescale(t) = SUM(CVT_TemporalCors % temp_eigenvec_flux(t, 1:nmonth, f) *&
                                 flux_modes(1:nmonth))
        END DO
        ! Multiply by the time-dependent standard deviation
        state_flux(i, j, 1:nmonth, f) = flux_prescale(1:nmonth) *                      &
                                        CVT_std % std_flux(i, j, 1:nmonth, f)
      END DO
    END DO


  CASE (2)
    ! Temporal covariances, spatial scale a fn of time
    ! This represents a symmetric square-root of the temporal correlation matrix
    ! Expect input to forward routine to be a function of lat/long/time/field
    ! Output of the forward routine to be a function of lat/long/time/field
    DO i = 1, ylat
      DO j = 1, ylon
        ! Project on to temporal modes
        DO mode = 1, nmonth
          flux_modes(mode) = SUM(CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, mode, f) * &
                                 chi_th_flux(i, j, 1:nmonth, f))
        END DO
        ! Multiply by standard deviation (flux_modes a function of timescale)
        flux_modes(1:nmonth) = flux_modes(1:nmonth) *                                       &
                               CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
        ! Project back to time
        DO t = 1, nmonth
          flux_prescale(t) = SUM(CVT_TemporalCors % temp_eigenvec_flux(t,1:nmonth, f) *     &
                                 flux_modes(1:nmonth))
        END DO
        ! Multiply by the time-dependent standard deviation
        state_flux(i,j,1:nmonth, f) = flux_prescale(1:nmonth) *                             &
                                      CVT_std % std_flux(i, j, 1:nmonth, f)
      END DO
    END DO

  END SELECT

END DO

END SUBROUTINE cvt_t




!=================================================================================================
SUBROUTINE cvt_t_adj (chi_th_flux_hat, & ! OUT Version of chi (for flux) between temp and horiz transforms
                      state_flux_hat )   ! IN  Flux field in real space (lat, long, month)

! Description:  (adjoint of the following)
! Performs the temporal control variable transform, and multiplies by the error std devs
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs the inverse vertical transform and multiplies by the std dev of the bg error
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Include nfields (R.N. Bannister)
!           01/2020  Use structures.
!           07/2021  Allow for different fields to have different correlation types
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, nmonth, nfields

USE onedvar_data, ONLY: CVT_TemporalCors, &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(OUT) :: chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL, INTENT(IN)  :: state_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)

! Local variables:
INTEGER           :: i, j, t, mode, f
REAL              :: flux_modes(1:nmonth)
REAL              :: flux_prescale(1:nmonth)

DO f = 1, nfields

  ! Decide on the temporal correlation option
  SELECT CASE (CVT_TemporalCors % temporal_covs(f))

  CASE (0)
    ! White noise in time (no temporal correlations)
    ! Expect input to forward routine to be a function of lat/long/time/field
    ! Output of the forward routine to be a function of lat/long/time/field
    chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, f) = state_flux_hat(1:ylat, 1:ylon, 1:nmonth, f) * &
                                                   CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, f)

  CASE (1)
    ! Temporal covariances, spatial scale a fn of timescale
    ! This is the non-symmetric temporal transform
    ! Expect input to forward routine to be a function of lat/long/timescale
    ! Output of the forward routine to be a function of lat/long/time
    DO i = 1, ylat
      DO j = 1, ylon
        ! Multiply by the time-dependent standard deviation
        flux_prescale(1:nmonth) = state_flux_hat(i,j,1:nmonth,f) *                          &
                                  CVT_std % std_flux(i, j, 1:nmonth, f)
        ! Project on to temporal modes
        DO mode = 1, nmonth
          flux_modes(mode) = SUM(CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, mode, f) * &
                                 flux_prescale(1:nmonth))
        END DO
        ! Multiply by standard deviation (flux_modes a function of timescale and field)
        chi_th_flux_hat(i, j, 1:nmonth, f) = flux_modes(1:nmonth) *                         &
                                             CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
      END DO
    END DO


  CASE (2)
    ! Temporal covariances, spatial scale a fn of time
    ! This represents a symmetric square-root of the temporal correlation matrix
    ! Expect input to forward routine to be a function of lat/long/time/field
    ! Output of the forward routine to be a function of lat/long/time/field
    DO i = 1, ylat
      DO j = 1, ylon
        ! Multiply by the time-dependent standard deviation
        flux_prescale(1:nmonth) = state_flux_hat(i,j,1:nmonth,f) *                                  &
                                  CVT_std % std_flux(i, j, 1:nmonth, f)
        ! Project on to temporal modes
        DO mode = 1, nmonth
          flux_modes(mode) = SUM(CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, mode, f) *         &
                                 flux_prescale(1:nmonth))
        END DO
        ! Multiply by standard deviation (flux_modes a function of timescale)
        flux_modes(1:nmonth) = flux_modes(1:nmonth) *                                               &
                               CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
        ! Project back to time
        DO t = 1, nmonth
          chi_th_flux_hat(i, j, t, f) = SUM(CVT_TemporalCors % temp_eigenvec_flux(t, 1:nmonth, f) * &
                                            flux_modes(1:nmonth))
        END DO
      END DO
    END DO

  END SELECT
END DO

END SUBROUTINE cvt_t_adj




!=================================================================================================
SUBROUTINE cvt_t_inv (chi_th_flux, & ! OUT Version of chi (for flux) between temp and horiz transforms
                      state_flux )   ! IN  Flux field in real space (lat, long, month)

! Description:
! Divides by the std of the bg err and performs the inverse temporal control variable transform
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Divides by the standard dev of the bg error and performs the inverse temporal transform
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Include nfields (R.N. Bannister)
!           01/2020  Use structures.
!           07/2021  Allow for different fields to have different correlation types
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, nmonth, nfields

USE onedvar_data, ONLY: CVT_TemporalCors, &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(OUT) :: chi_th_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL, INTENT(IN)  :: state_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)

! Local variables:
INTEGER           :: i, j, t, mode, f
REAL              :: flux_modes(1:nmonth)
REAL              :: flux_prescale(1:nmonth)


DO f = 1, nfields

  ! Decide on the temporal correlation option
  SELECT CASE (CVT_TemporalCors % temporal_covs(f))

  CASE (0)
    ! White noise in time (no temporal correlations)
    ! Expect input to forward routine to be a function of lat/long/time
    ! Output of the forward routine to be a function of lat/long/time
    chi_th_flux(1:ylat, 1:ylon, 1:nmonth, f) = state_flux(1:ylat, 1:ylon, 1:nmonth, f) / &
                                               CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, f)

  CASE (1)
    ! Temporal covariances, spatial scale a fn of timescale
    ! This is the non-symmetric temporal transform
    ! Expect input to forward routine to be a function of lat/long/timescale
    ! Output of the forward routine to be a function of lat/long/time
    DO i = 1, ylat
      DO j = 1, ylon
        ! Divide by the time-dependent standard deviation
        flux_prescale(1:nmonth) = state_flux(i,j,1:nmonth, f) /                             &
                                  CVT_std % std_flux(i, j, 1:nmonth, f)
        ! Project on to temporal modes
        DO mode = 1, nmonth
          flux_modes(mode) = SUM(CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, mode, f) * &
                                 flux_prescale(1:nmonth))
        END DO
        ! Divide by standard deviation (flux_modes a function of timescale)
        chi_th_flux(i, j, 1:nmonth, f) = flux_modes(1:nmonth) /                             &
                                         CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
      END DO
    END DO

  CASE (2)
    ! Temporal covariances, spatial scale a fn of time
    ! This represents a symmetric square-root of the temporal correlation matrix
    ! Expect input to forward routine to be a function of lat/long/time
    ! Output of the forward routine to be a function of lat/long/time
    DO i = 1, ylat
      DO j = 1, ylon
        ! Divide by the time-dependent standard deviation
        flux_prescale(1:nmonth) = state_flux(i,j,1:nmonth,f) /                                  &
                                  CVT_std % std_flux(i, j, 1:nmonth,f)
        ! Project on to temporal modes
        DO mode = 1, nmonth
          flux_modes(mode) = SUM(CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, mode, f) *     &
                                 flux_prescale(1:nmonth))
        END DO
        ! Divide by standard deviation (flux_modes a function of timescale)
        flux_modes(1:nmonth) = flux_modes(1:nmonth) /                                           &
                               CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
        ! Project back to time
        DO t = 1, nmonth
          chi_th_flux(i, j, t, f) = SUM(CVT_TemporalCors % temp_eigenvec_flux(t, 1:nmonth, f) * &
                                        flux_modes(1:nmonth))
        END DO
      END DO
    END DO
  END SELECT

END DO

END SUBROUTINE cvt_t_inv




!=================================================================================================
SUBROUTINE cvt_t_inv_adj (chi_th_flux_hat, & ! IN  Version of chi (for flux) between temp and horiz transforms
                          state_flux_hat )   ! OUT Flux field in real space (lat, long, month)

! Description:
! Performs the adjoint of the inverse temporal control variable transform, and divides by the error std devs
!
! Copyright:
!    This software was developed by Ross Bannister, NCEO
!    Copyright 2019, NCEO, All Rights Reserved.
!
! Method :
! ------
! Performs the temporal transform and multiplies by the std dev of the bg error
! References are to:
! Bannister, Modelling spatial correlations in DA
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0      07/2019  Initial version (R.N. Bannister)
!           09/2019  Include nfields (R.N. Bannister)
!           09/2019  Developed inverse of adjoint (R.N. Bannister)
!           01/2020  Use structures.
!           07/2021  Allow for different fields to have different correlation types
!
! Code Description:
!   Language:           Fortran 90.

! Declarations:

USE main_data,    ONLY: ylon, ylat, nmonth, nfields

USE onedvar_data, ONLY: CVT_TemporalCors, &
                        CVT_std

IMPLICIT NONE

! Arguments:
REAL, INTENT(IN)  :: chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL, INTENT(OUT) :: state_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)

! Local variables:
INTEGER           :: i, j, t, mode, f
REAL              :: flux_modes(1:nmonth)
REAL              :: flux_prescale(1:nmonth)


DO f = 1, nfields

  ! Decide on the temporal correlation option
  SELECT CASE (CVT_TemporalCors % temporal_covs(f))

  CASE (0)
    ! White noise in time (no temporal correlations)
    ! Expect input to forward routine to be a function of lat/long/time/field
    ! Output of the forward routine to be a function of lat/long/time/field
    state_flux_hat(1:ylat, 1:ylon, 1:nmonth, f) = chi_th_flux_hat(1:ylat, 1:ylon, 1:nmonth, f) / &
                                                  CVT_std % std_flux(1:ylat, 1:ylon, 1:nmonth, f)

  CASE (1)
    ! Temporal covariances, spatial scale a fn of timescale
    ! This is the non-symmetric temporal transform
    ! Expect input to forward routine to be a function of lat/long/timescale/field
    ! Output of the forward routine to be a function of lat/long/time/field
    DO i = 1, ylat
      DO j = 1, ylon
        ! Divide by standard deviation (flux_modes a function of timescale)
        flux_modes(1:nmonth) = chi_th_flux_hat(i, j, 1:nmonth, f) /                    &
                               CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
        ! Project back to time
        DO t = 1, nmonth
          flux_prescale(t) = SUM(CVT_TemporalCors % temp_eigenvec_flux(t,1:nmonth,f) * &
                                 flux_modes(1:nmonth))
        END DO
        ! Divide by the time-dependent standard deviation
        state_flux_hat(i,j,1:nmonth,f) = flux_prescale(1:nmonth) /                     &
                                         CVT_std % std_flux(i, j, 1:nmonth, f)
      END DO
    END DO


  CASE (2)
    ! Temporal covariances, spatial scale a fn of time
    ! This represents a symmetric square-root of the temporal correlation matrix
    ! Expect input to forward routine to be a function of lat/long/time/field
    ! Output of the forward routine to be a function of lat/long/time/field
    DO i = 1, ylat
      DO j = 1, ylon
        ! Project on to temporal modes
        DO mode = 1, nmonth
          flux_modes(mode) = SUM(CVT_TemporalCors % temp_eigenvec_flux(1:nmonth, mode, f) * &
                                 chi_th_flux_hat(i, j, 1:nmonth, f))
        END DO
        ! Divide by standard deviation (flux_modes a function of timescale)
        flux_modes(1:nmonth) = flux_modes(1:nmonth) /                                       &
                               CVT_TemporalCors % std_tspec_flux(1:nmonth, f)
        ! Project back to time
        DO t = 1, nmonth
          flux_prescale(t) = SUM(CVT_TemporalCors % temp_eigenvec_flux(t,1:nmonth, f) *     &
                                 flux_modes(1:nmonth))
        END DO
        ! Divide by the time-dependent standard deviation
        state_flux_hat(i,j,1:nmonth, f) = flux_prescale(1:nmonth) /                         &
                                          CVT_std % std_flux(i, j, 1:nmonth, f)
      END DO
    END DO

  END SELECT

END DO

END SUBROUTINE cvt_t_inv_adj




!=================================================================================================
SUBROUTINE RidgeRegression (elements,     &
                            Evals,        &
                            small)
! Do a ridge regression to remove negative (and zero) eigenvalues

IMPLICIT NONE

! Arguments:
INTEGER, INTENT(IN)    :: elements
REAL,    INTENT(INOUT) :: Evals(1:elements)
REAL,    INTENT(IN)    :: small

! Local variables
INTEGER                :: el
REAL                   :: total, addon, smallest, largest

! Find smallest and largest values
CALL ExtremeVals (elements, Evals(1:elements), smallest, largest)
PRINT *, 'Min eval =', smallest
PRINT *, 'Max eval =', largest
total = SUM(Evals(1:elements))
addon = total * small / REAL(elements)

Evals(1:elements) = Evals(1:elements) + addon

DO el = 1, elements
  IF (Evals(el) < 0.0) THEN
    PRINT *, 'Element ', el, ' zero ', Evals(el)
    Evals(el) = 0.0
  END IF
END DO

END SUBROUTINE RidgeRegression



!=================================================================================================
SUBROUTINE ExtremeVals (nels, x, smallest, largest)
! Find the minimum value
IMPLICIT NONE

! Arguments:
INTEGER, INTENT(IN)  :: nels
REAL,    INTENT(IN)  :: x(1:nels)
REAL,    INTENT(OUT) :: smallest
REAL,    INTENT(OUT) :: largest

! Local variables
INTEGER              :: i
REAL                 :: sm, la

sm = x(1)
la = x(1)
DO i = 2, nels
  IF (x(i) < sm) sm = x(i)
  IF (x(i) > la) la = x(i)
END DO

smallest = sm
largest  = la
END SUBROUTINE ExtremeVals




! The following routine is purely for testing
!=================================================================================================
SUBROUTINE ChangeHorizGrid_tomcat2shtools (Grid_in,  & ! IN  Data on tomcat grid
                                           Grid_out)   ! OUT Data on shtools grid

! Change grid from tomcat to shtools


USE main_data,    ONLY: ylon, ylat

USE onedvar_data, ONLY: L,                       &
                        CVT_HorizCors

IMPLICIT NONE

! Note: We set L=ylat
! The number of longitudes is ylon in tomcat
!                             L+1 in shtools
! The number of latitudes is  ylat in tomcat
!                             2L+1 in shtools


! Arguments:
REAL,    INTENT(IN)  :: grid_in(1:ylat, 1:ylon, 1:1)
REAL,    INTENT(OUT) :: grid_out(1:L+1, 1:2*L+1, 1:1)

CALL ChangeHorizGrid (1,                                      &
                      ylon,                                   & ! No of longitudes shtools grid
                      ylat,                                   & ! No of latitudes shtools grid
                      grid_in(1:ylat, 1:ylon, 1:1),           & ! Input data
                      CVT_HorizCors % TOMCAT_longs(1:ylon),   & ! Longs on tomcat grid
                      CVT_HorizCors % TOMCAT_lats(1:ylat),    & ! Lats on tomcat grid
                      2*L+1,                                  & ! No of longitudes shtools grid
                      L+1,                                    & ! No of latitudes shtools grid
                      grid_out(1:L+1, 1:2*L+1, 1:1),          & ! Output data
                      CVT_HorizCors % SHtools_longs(1:2*L+1), & ! Longitudes of shtools grid
                      CVT_HorizCors % SHtools_lats(1:L+1),    & ! Latitudes of shtools grid
                      CVT_HorizCors % TOMCAT_to_SH )


END SUBROUTINE ChangeHorizGrid_tomcat2shtools


! The following routine is purely for testing
!=================================================================================================
SUBROUTINE ChangeHorizGrid_shtools2tomcat (Grid_in,  & ! IN  Data on shtools grid
                                           Grid_out)   ! OUT Data on tomcat grid

! Change grid from shtools to tomcat


USE main_data,    ONLY: ylon, ylat

USE onedvar_data, ONLY: L,                       &
                        CVT_HorizCors

IMPLICIT NONE

! Note: We set L=ylat
! The number of longitudes is ylon in tomcat
!                             L+1 in shtools
! The number of latitudes is  ylat in tomcat
!                             2L+1 in shtools


! Arguments:
REAL,    INTENT(IN)  :: grid_in(1:L+1, 1:2*L+1, 1:1)
REAL,    INTENT(OUT) :: grid_out(1:ylat, 1:ylon, 1:1)


CALL ChangeHorizGrid (1,                                      &
                      2*L+1,                                  & ! No of longitudes shtools grid
                      L+1,                                    & ! No of latitudes shtools grid
                      grid_in(1:L+1, 1:2*L+1, 1:1),           & ! Output data
                      CVT_HorizCors % SHtools_longs(1:2*L+1), & ! Longitudes of shtools grid
                      CVT_HorizCors % SHtools_lats(1:L+1),    & ! Latitudes of shtools grid
                      ylon,                                   & ! No of longitudes shtools grid
                      ylat,                                   & ! No of latitudes shtools grid
                      grid_out(1:ylat, 1:ylon, 1:1),          & ! Input data
                      CVT_HorizCors % TOMCAT_longs(1:ylon),   & ! Longs on tomcat grid
                      CVT_HorizCors % TOMCAT_lats(1:ylat),    & ! Lats on tomcat grid
                      CVT_HorizCors % SH_to_TOMCAT )

END SUBROUTINE ChangeHorizGrid_shtools2tomcat




!=================================================================================================
REAL FUNCTION CorrelationFn (diff, cor_tpe, scale)
! Function to return the correlation between two points separated by diff
! where the decorrelation scale is scale.
! diff and scale should be in the same units.
! cor_type specifies the form of the correlation from the options below

IMPLICIT NONE

REAL,     INTENT(IN) :: diff
INTEGER,  INTENT(IN) :: cor_tpe
REAL,     INTENT(IN) :: scale

REAL                 :: cor
REAL                 :: scaled_diff

scaled_diff = diff / scale


SELECT CASE (cor_tpe)

CASE (1)
  ! Lorentzian
  cor = 1.0 / (1.0 + scaled_diff * scaled_diff)

CASE (2)
  ! Gaussian
  cor = EXP(-1.0 * scaled_diff * scaled_diff / 2.0)

CASE (3)
  ! SOAR
  cor = (1.0 + ABS(scaled_diff)) * EXP(-1.0 * ABS(scaled_diff))

CASE (4)
  ! Exponential
  cor = EXP(-1.0 * ABS(scaled_diff))

END SELECT

CorrelationFn = cor

END FUNCTION CorrelationFn

!=================================================================================================
