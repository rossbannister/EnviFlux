PROGRAM ImpliedCov

! This program reads in covariance file and computes implied covariances

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, nfields, dim_1d

USE onedvar_data, ONLY: n1d,                           &
                        CVT_std,                       &
                        CVT_HorizCors,                 &
                        CVT_VertCors,                  &
                        CVT_TemporalCors,              &
                        L,                             &
                        ALPsize,                       &
                        GLatEq

IMPLICIT none

INCLUDE "cvt.interface"

! Local variables

REAL              :: input_v(1:ylat, 1:ylon, 1:ylev)
REAL              :: output_v(1:ylat, 1:ylon, 1:ylev)
REAL              :: intermediate_v(1:ylat, 1:ylon, 1:ylev)

REAL              :: input_h_tracer(1:ylat, 1:ylon, 1:ylev)
REAL              :: output_h_tracer(1:ylat, 1:ylon, 1:ylev)
REAL              :: intermediate_h(1:n1d)

REAL              :: input_h_nointerp(1:L+1, 1:2*L+1, 1)
REAL              :: output_h_nointerp(1:L+1, 1:2*L+1, 1)
REAL              :: intermediate_h_nointerp(1:2, 0:L, 0:L)

REAL              :: input_h_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL              :: output_h_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)

REAL              :: input_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL              :: output_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)
REAL              :: intermediate_t(1:ylat, 1:ylon, 1:nmonth, 1:nfields)

REAL              :: input_total(1:dim_1d)
REAL              :: output_total(1:dim_1d)
REAL              :: intermediate(1:n1d)

INTEGER           :: ll, lp, m, index, j, indexll, indexlp, t, lon, err, lat, lev
REAL              :: tot, weight


! ==================================================================================

! Allocate and initialise CVT
CALL AllocateCVT (.FALSE.)
CALL InitializeHorizTrans ()



PRINT *, '===== Gaussian lat index corresponding to the equator ', GLatEq
DO lat = 1, L+1
  PRINT *, lat, CVT_HorizCors % GaussianWts(lat), CVT_HorizCors % GaussianCosCoLats(lat)
END DO

! Check an orthogonality relationship for associated Legendre polynomials
PRINT *, '===== Checking orthongonality of associated Legendre polynomials'
m = 0
DO ll = m, L
  indexll = CVT_HorizCors % Plm_index(ll,m)
  DO lp = m, L
    tot = 0.0
    indexlp = CVT_HorizCors % Plm_index(lp,m)
    DO j = 1, L+1
      tot = tot + CVT_HorizCors % assocLegPoly(j, indexll) * &
                  CVT_HorizCors % assocLegPoly(j, indexlp) * &
                  CVT_HorizCors % GaussianWts(j)
    END DO
    PRINT *, ll, lp, tot
  END DO
END DO


! ==================================================================================

! ----------------------------------------------------------------------------------
! Read-in the transform information
! ----------------------------------------------------------------------------------
PRINT *, '===== Reading in the CVT data'
CALL cvt_matrices_input ('CVT_calib.nc')


! ----------------------------------------------------------------------------------
! Test the vertical transform
! ----------------------------------------------------------------------------------
PRINT *, '===== Implied vertical transform'

CVT_std % std_tracer(1:ylat, 1:ylon, 1:ylev) = 1.0

input_v(1:ylat, 1:ylon, 1:ylev) = 0.0
input_v(15, 30, 1)              = 1.0
input_v(17, 30, 10)             = 1.0
input_v(19, 30, 50)             = 1.0
input_v(21, 30, ylev)           = 1.0

CALL cvt_v_adj ( intermediate_v(1:ylat, 1:ylon, 1:ylev),         &
                 input_v(1:ylat, 1:ylon, 1:ylev) )

CALL cvt_v     ( intermediate_v(1:ylat, 1:ylon, 1:ylev),         &
                 output_v(1:ylat, 1:ylon, 1:ylev) )

CALL Write_one_field ('Implied_v.nc',                            &
                      ylat, ylon, ylev,                          &
                      output_v(1:ylat, 1:ylon, 1:ylev),          &
                      'output_v',                                &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),       &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),      &
                      CVT_VertCors % alts(1:ylev))

OPEN (12, file='VTrans_test.dat')
DO lev = 1, ylev
  WRITE (12,*) lev, CVT_VertCors % alts(lev), output_v(15,30,lev), output_v(17,30,lev),  &
               output_v(19,30,lev), output_v(21,30,lev)
END DO
CLOSE (12)



! ----------------------------------------------------------------------------------
! Test the horizontal transform (with the interpolation step)
! ----------------------------------------------------------------------------------
PRINT *, '===== Full horizontal transform'

input_h_tracer(1:ylat, 1:ylon, 1:ylev)              = 0.0
input_h_tracer(16, 30, 1)                           = 1.0
input_h_tracer(16, 60, 2)                           = 1.0
input_h_tracer(10, 15, 3)                           = 1.0
input_h_tracer(10, 45, 4)                           = 1.0
input_h_tracer(2,  20, 5)                           = 1.0
input_h_tracer(2,  25, 6)                           = 1.0

input_h_flux  (1:ylat, 1:ylon, 1:nmonth, 1:nfields) = 0.0
input_h_flux  (16, 30, 1, 1)                        = 1.0
input_h_flux  (16, 60, 2, 1)                        = 1.0
input_h_flux  (10, 15, 3, 1)                        = 1.0
input_h_flux  (10, 45, 4, 1)                        = 1.0
input_h_flux  (2,  20, 5, 1)                        = 1.0
input_h_flux  (2,  25, 6, 1)                        = 1.0

CALL cvt_h_adj ( intermediate_h(1:n1d),                              &
                 input_h_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  &
                 input_h_tracer(1:ylat, 1:ylon, 1:ylev),             &
                 .TRUE.,                                             &
                 .TRUE. )

CALL cvt_h     ( intermediate_h(1:n1d),                              &
                 output_h_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                 output_h_tracer(1:ylat, 1:ylon, 1:ylev),            &
                 .TRUE.,                                             &
                 .TRUE. )


CALL Write_one_field ('Implied_h_tracer_interp.nc',                  &
                      ylat, ylon, ylev,                              &
                      output_h_tracer(1:ylat, 1:ylon, 1:ylev),       &
                      'output_h_tracer_interp',                      &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),           &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),          &
                      CVT_VertCors % alts(1:ylev))

CALL Write_one_field ('Implied_h_flux_interp.nc',                    &
                      ylat, ylon, nmonth,                            &
                      output_h_flux(1:ylat, 1:ylon, 1:nmonth, 1),    &
                      'output_h_flux_interp',                        &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),           &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),          &
                      CVT_TemporalCors % months(1:nmonth))



! ----------------------------------------------------------------------------------
! Test the horizontal transform (without the interpolation step -- pure horizontal transform)
! ----------------------------------------------------------------------------------
PRINT *, '===== No-interpolation horizontal transform'

input_h_nointerp(1:L+1, 1:2*L+1, 1) = 0.0


lat = 17
input_h_nointerp(lat, 18, 1) = 1.0

! ***** Adjoint step *****
CALL Spherical_adj (intermediate_h_nointerp(1:2, 0:L, 0:L), & ! OUT
                    input_h_nointerp(1:L+1, 1:2*L+1, 1))      ! IN


! ***** Multiply by the horizontal spectrum *****
! Loop over total wavenumber
DO ll = 0, L
  weight = CVT_HorizCors % var_hspec_flux(ll,1,1)
  ! Special case for m=0 (no Sine contribution)
  ! The cosine term
  intermediate_h_nointerp(1,ll,0) = intermediate_h_nointerp(1,ll,0) * weight
  intermediate_h_nointerp(2,ll,0) = intermediate_h_nointerp(2,ll,0) * weight
  IF (ll > 0) THEN
    ! Loop over remaining meridional wavenumbers
    DO m = 1, ll
      ! The cosine term
      intermediate_h_nointerp(1,ll,m) = intermediate_h_nointerp(1,ll,m) * weight
      ! The sine term
      intermediate_h_nointerp(2,ll,m) = intermediate_h_nointerp(2,ll,m) * weight
    END DO
  END IF
END DO


! ***** Forward step *****
CALL Spherical_forward (intermediate_h_nointerp(1:2, 0:L, 0:L), & ! OUT
                        output_h_nointerp(1:L+1, 1:2*L+1, 1))     ! IN


CALL Write_one_field ('Implied_h_nointerp.nc',                &
                      L+1, 2*L+1, 1,                          &
                      output_h_nointerp(1:L+1, 1:2*L+1, 1),   &
                      'output_h_nointerp',                    &
                      CVT_HorizCors % SHtools_lats(1:L+1),    &
                      CVT_HorizCors % SHtools_longs(1:2*L+1), &
                      CVT_VertCors % alts(1:1))



! ----------------------------------------------------------------------------------
! Test the temporal transform
! ----------------------------------------------------------------------------------
PRINT *, '===== Temporal transform'

input_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) = 0.0

DO t = 1, nmonth
  input_t_flux(1, t, t, 1) = 1.0
END DO

CALL cvt_t_adj ( intermediate_t(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  &
                 input_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) )

CALL cvt_t     ( intermediate_t(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  &
                 output_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) )

DO t = 1, nmonth
  PRINT *, t, output_t_flux(1, t, t, 1)
END DO

! View with y-z plot x=1 in xconv
CALL Write_one_field ('Implied_t.nc',                                &
                      ylat, ylon, nmonth,                            &
                      output_t_flux(1:ylat, 1:ylon, 1:nmonth, 1),    &
                      'output_t',                                    &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),           &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),          &
                      CVT_TemporalCors % months(1:nmonth))



! ----------------------------------------------------------------------------------
! Test the whole transform
! ----------------------------------------------------------------------------------
PRINT *, '===== The whole transform'

! Set the input field for tracer and flux

input_v(1:ylat, 1:ylon, 1:ylev)                   = 0.0
input_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) = 0.0

input_v(17, 40, 30)                               = 1.0
input_t_flux(6,  40, 12, 1)                       = 1.0

input_t_flux(27, 8,  3,  1)                       = 1.0


!input_t_flux(26, 7,  3,  1)                       = 1.0 /9.0
!input_t_flux(26, 8,  3,  1)                       = 1.0 /9.0
!input_t_flux(26, 9,  3,  1)                       = 1.0 /9.0
!input_t_flux(27, 7,  3,  1)                       = 1.0 /9.0
!input_t_flux(27, 8,  3,  1)                       = 1.0 /9.0
!input_t_flux(27, 9,  3,  1)                       = 1.0 /9.0
!input_t_flux(28, 7,  3,  1)                       = 1.0 /9.0
!input_t_flux(28, 8,  3,  1)                       = 1.0 /9.0
!input_t_flux(28, 9,  3,  1)                       = 1.0 /9.0


!input_t_flux(27, 8, 2, 1)                         = 1.0 / 3.0
!input_t_flux(27, 8, 3, 1)                         = 1.0 / 3.0
!input_t_flux(27, 8, 4, 1)                         = 1.0 / 3.0



! Print latitudes and longitudes of these impulse points
PRINT *, 'Latitude and longitudes of these impulse points'
PRINT *, 17,  40, CVT_HorizCors % TOMCAT_lats(17), CVT_HorizCors % TOMCAT_longs(40)
PRINT *, 6,   40, CVT_HorizCors % TOMCAT_lats(6),  CVT_HorizCors % TOMCAT_longs(40)
PRINT *, 27,  8,  CVT_HorizCors % TOMCAT_lats(27), CVT_HorizCors % TOMCAT_longs(8)

intermediate(1:n1d)                               = 0.0
output_total(1:dim_1d)                            = 0.0

! Place these fields into a single 1d input array
CALL Physical_multid2linear (input_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  &
                             input_v(1:ylat, 1:ylon, 1:ylev),                    &
                             input_total(1:dim_1d))

! Perform total adjoint transform
CALL cvt_total_adj (intermediate(1:n1d),   &
                    input_total(1:dim_1d), &
                    .TRUE., .TRUE.)


! Perform total forward transform
CALL cvt_total (intermediate(1:n1d),    &
                output_total(1:dim_1d), &
                .TRUE., .TRUE.)

! Extract the fields from the single 1d output array
CALL Physical_linear2multid (output_t_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                             output_v(1:ylat, 1:ylon, 1:ylev),                   &
                             output_total(1:dim_1d))

! Output the tracer field
CALL Write_one_field ('Implied_total_tracer.nc',                     &
                      ylat, ylon, ylev,                              &
                      output_v(1:ylat, 1:ylon, 1:ylev),              &
                      'implied_total_tracer',                        &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),           &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),          &
                      CVT_VertCors % alts(1:ylev))

! Output the flux field
CALL Write_one_field ('Implied_total_flux.nc',                       &
                      ylat, ylon, nmonth,                            &
                      output_t_flux(1:ylat, 1:ylon, 1:nmonth, 1),    &
                      'implied_total_flux',                          &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),           &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),          &
                      CVT_TemporalCors % months(1:nmonth))



! Deallocate CVT arrays needed
CALL UnInitializeHorizTrans ()
CALL DeallocateCVT(.FALSE.)


END PROGRAM ImpliedCov
