! For use with the non-balanced horizontal transform


PROGRAM AdjointTestsMain

USE main_data,    ONLY: ylon, ylat, ylev, nmonth, dim_1d, nfields

USE cparam,       ONLY: pi

USE onedvar_data, ONLY: L,                             &
                        n1d,                           &
                        ALPsize,                       &
                        SHtools_Normalization,         &
                        SHtools_ApplyCSphase,          &
                        SHtools_CRnorm,                &
                        CVT_std,                       &
                        CVT_HorizCors,                 &
                        CVT_VertCors,                  &
                        CVT_TemporalCors


IMPLICIT none

INCLUDE "cvt.interface"

INTEGER :: i, j, k, mon, lev, err, f
REAL    :: data_SHtools(1:L+1, 1:2*L+1, 1:ylev)
REAL    :: data_TOMCAT(1:ylat, 1:ylon, 1:ylev)
REAL    :: data_SHtools1(1:L+1, 1:2*L+1, 1:ylev)
REAL    :: data_spherical(1:2, 0:L, 0:L)
REAL    :: data_realspace(1:L+1, 1:2*L+1)
REAL    :: data_spherical1(1:2, 0:L, 0:L)
REAL    :: data_realspace1(1:L+1, 1:2*L+1)

REAL    :: intermediate(1:2, 0:L, 1:L+1)
REAL    :: intermediate1(1:2, 0:L, 1:L+1)

REAL, ALLOCATABLE    :: chi(:)
REAL, ALLOCATABLE    :: chi_hat(:)
REAL, ALLOCATABLE    :: chi_flux(:,:,:,:)
REAL, ALLOCATABLE    :: chi_tracer(:,:,:)
REAL, ALLOCATABLE    :: field_tracer(:,:,:)
REAL, ALLOCATABLE    :: field_tracer_hat(:,:,:)
REAL, ALLOCATABLE    :: chi_tracer_hat(:,:,:)
REAL, ALLOCATABLE    :: field_flux(:,:,:,:)
REAL, ALLOCATABLE    :: field_flux_hat(:,:,:,:)
REAL, ALLOCATABLE    :: chi_flux_hat(:,:,:,:)
REAL, ALLOCATABLE    :: zvar1d(:)
REAL, ALLOCATABLE    :: zvar1d_hat(:)


REAL    :: field_u(0:ylat+1, 0:ylon+1, 1:ylev)
REAL    :: field_v(0:ylat+1, 0:ylon+1, 1:ylev)
REAL    :: field_w(0:ylat+1, 0:ylon+1, 1:ylev)

REAL    :: InnerProdReal, InnerProdSpec


! Allocate CVT arrays needed
CALL AllocateCVT(.FALSE.)
CALL InitializeHorizTrans ()


! ----------------------------------------------------------------------------------
! Read-in the transform information
! ----------------------------------------------------------------------------------
CALL cvt_matrices_input ('CVT_calib.nc')


! =====================================================================================
PRINT *, '----- Testing the code that does the interpolation -----'
! Test the change of grid
DO i = 1, L+1
  DO j = 1, 2*L+1
    data_SHtools(i,j,1:ylev) = SIN(REAL(i)*6.0) + REAL(j*j)/REAL(L*L)
  END DO
END DO

PRINT *, 'Calling ChangeHorizGrid'
CALL ChangeHorizGrid (ylev,                                 &
                      2*L+1,                                &
                      L+1,                                  &
                      data_SHtools(1:L+1, 1:2*L+1, 1:ylev), &
                      CVT_HorizCors % SHtools_longs(1:2*L+1),&
                      CVT_HorizCors % SHtools_lats(1:L+1),  &
                      ylon,                                 &
                      ylat,                                 &
                      data_TOMCAT(1:ylat, 1:ylon, 1:ylev),  &
                      CVT_HorizCors % TOMCAT_longs(1:ylon), &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),  &
                      CVT_HorizCors % SH_to_TOMCAT)
PRINT *, 'Done'

PRINT *, 'Outputting data'
CALL Write_one_field ('shfield.nc',                         &
                      L+1, 2*L+1, ylev,                     &
                      data_SHtools(1:L+1, 1:2*L+1, 1:ylev), &
                      'data_SHtools',                       &
                      CVT_HorizCors % SHtools_lats(1:L+1),  &
                      CVT_HorizCors % SHtools_longs(1:2*L+1),&
                      CVT_VertCors % alts(1:ylev))

CALL Write_one_field ('tcfield.nc',                         &
                      ylat, ylon, ylev,                     &
                      data_TOMCAT(1:ylat, 1:ylon, 1:ylev),  &
                      'data_TOMCAT',                        &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),  &
                      CVT_HorizCors % TOMCAT_longs(1:ylon), &
                      CVT_VertCors % alts(1:ylev))

PRINT *, 'Done'



! =====================================================================================
! Adjoint test of swapbounds
PRINT *, '----- Adjoint test of swapbounds routine -----'
! Put random numbers in u
field_u(0:ylat+1, 0:ylon+1, 1:ylev) = 0.0
CALL RANDOM_NUMBER (field_u(1:ylat, 1:ylon, 1:ylev))
! Copy to v
field_v(0:ylat+1, 0:ylon+1, 1:ylev) = field_u(0:ylat+1, 0:ylon+1, 1:ylev)

! Put field through swap_halos
CALL Swap_halos (ylev, ylon, ylat, &
                 field_v(0:ylat+1, 0:ylon+1, 1:ylev))
! Copy to w
field_w(0:ylat+1, 0:ylon+1, 1:ylev) = field_v(0:ylat+1, 0:ylon+1, 1:ylev)

! Put field through adjoint
CALL Swap_halos_adj (ylev, ylon, ylat, &
                     field_w(0:ylat+1, 0:ylon+1, 1:ylev))

PRINT *, 'LHS = ', SUM(field_v(0:ylat+1, 0:ylon+1, 1:ylev) * &
                       field_v(0:ylat+1, 0:ylon+1, 1:ylev))
PRINT *, 'RHS = ', SUM(field_u(0:ylat+1, 0:ylon+1, 1:ylev) * &
                       field_w(0:ylat+1, 0:ylon+1, 1:ylev))



! =====================================================================================
! Adjoint test of Physical_linear2multid
PRINT *, '----- Adjoint test of Physical_multid2linear routine -----'
PRINT *, dim_1d
ALLOCATE (field_tracer(1:ylat, 1:ylon, 1:ylev))
ALLOCATE (field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ALLOCATE (zvar1d(1:dim_1d))
ALLOCATE (field_tracer_hat(1:ylat, 1:ylon, 1:ylev))
ALLOCATE (field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields))

! Put random numbers in field_tracer and field_flux
CALL RANDOM_NUMBER (field_tracer(1:ylat, 1:ylon, 1:ylev))
CALL RANDOM_NUMBER (field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))

! Forward routine
CALL Physical_multid2linear ( field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),  & !IN  Multi-dimensional flux field
                              field_tracer(1:ylat, 1:ylon, 1:ylev),             & !IN  Multi-dimensional tracer field
                              zvar1d(1:dim_1d) )                                  !OUT Linear state
! Adjoint routine
CALL Physical_linear2multid ( field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields),      & !OUT Multi-dimensional flux field
                              field_tracer_hat(1:ylat, 1:ylon, 1:ylev),                 & !OUT Multi-dimensional tracer field
                              zvar1d(1:dim_1d) )              !IN  Linear state

PRINT *, 'LHS = ', SUM(zvar1d(1:dim_1d) * zvar1d(1:dim_1d))
PRINT *, 'RHS = ', SUM(field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) *       &
                       field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)) +  &
                   SUM(field_tracer(1:ylat, 1:ylon, 1:ylev) *                  &
                       field_tracer_hat(1:ylat, 1:ylon, 1:ylev))

DEALLOCATE (field_tracer, field_flux, zvar1d, field_tracer_hat, field_flux_hat)





! =====================================================================================
! Adjoint test of ChangeHorizGrid
PRINT *, '----- Adjoint test of ChangeHorizGrid routine -----'
CALL RANDOM_NUMBER (data_SHtools(1:L+1, 1:2*L+1, 1:ylev))
! Put field through ChangeHorizGrid
CALL ChangeHorizGrid (ylev,                                 &
                      2*L+1,                                &
                      L+1,                                  &
                      data_SHtools(1:L+1, 1:2*L+1, 1:ylev), &
                      CVT_HorizCors % SHtools_longs(1:2*L+1),&
                      CVT_HorizCors % SHtools_lats(1:L+1),  &
                      ylon,                                 &
                      ylat,                                 &
                      data_TOMCAT(1:ylat, 1:ylon, 1:ylev),  &
                      CVT_HorizCors % TOMCAT_longs(1:ylon), &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),  &
                      CVT_HorizCors % SH_to_TOMCAT)
! Put field through ChangeHorizGrid_adj
CALL ChangeHorizGrid_adj (ylev,                                 &
                          2*L+1,                                &
                          L+1,                                  &
                          data_SHtools1(1:L+1, 1:2*L+1, 1:ylev),&
                          CVT_HorizCors % SHtools_longs(1:2*L+1),&
                          CVT_HorizCors % SHtools_lats(1:L+1),  &
                          ylon,                                 &
                          ylat,                                 &
                          data_TOMCAT(1:ylat, 1:ylon, 1:ylev),  &
                          CVT_HorizCors % TOMCAT_longs(1:ylon), &
                          CVT_HorizCors % TOMCAT_lats(1:ylat),  &
                          CVT_HorizCors % SH_to_TOMCAT)

PRINT *, 'LHS = ', SUM(data_TOMCAT(1:ylat, 1:ylon, 1:ylev) *    &
                       data_TOMCAT(1:ylat, 1:ylon, 1:ylev))
PRINT *, 'RHS = ', SUM(data_SHtools(1:L+1, 1:2*L+1, 1:ylev) *   &
                       data_SHtools1(1:L+1, 1:2*L+1, 1:ylev))





! =====================================================================================
! Adjoint test of spherical transform
PRINT *, '----- Adjoint test of spherical transform -----'
CALL RANDOM_NUMBER (data_spherical(1:2, 0:L, 0:L))

CALL Spherical_forward (data_spherical(1:2, 0:L, 0:L),     &  ! IN
                        data_realspace(1:L+1, 1:2*L+1))       ! OUT
PRINT *, 'LHS = ', SUM(data_realspace(1:L+1, 1:2*L+1) * data_realspace(1:L+1, 1:2*L+1))

CALL Spherical_adj (data_spherical1(1:2, 0:L, 0:L),        &  ! OUT
                    data_realspace(1:L+1, 1:2*L+1))           ! IN
PRINT *, 'RHS = ', InnerProdSpec (L,                               &
                                  data_spherical(1:2, 0:L, 0:L),   &
                                  data_spherical1(1:2, 0:L, 0:L))


! =====================================================================================
! Adjoint test of cvt_h
PRINT *, '----- Adjoint test of cvt_h -----'
PRINT *, 'Value of dim_1d =', dim_1d
PRINT *, 'Value of n1d    =', n1d
ALLOCATE (chi(1:n1d))
ALLOCATE (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ALLOCATE (chi_tracer(1:ylat, 1:ylon, 1:ylev))
ALLOCATE (chi_hat(1:n1d))

! Put random numbers in control variable
CALL RANDOM_NUMBER (chi(1:n1d))
! Input into cvt_h
CALL cvt_h (chi(1:n1d),                                        &
            chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),     &
            chi_tracer(1:ylat, 1:ylon, 1:ylev),                &
            .TRUE., .TRUE.)
! Input into cvt_h_adj
CALL cvt_h_adj (chi_hat(1:n1d),                                &
                chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                chi_tracer(1:ylat, 1:ylon, 1:ylev),            &
                .TRUE., .TRUE.)
PRINT *, 'LHS = ', SUM(chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) * &
                       chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields)) +&
                   SUM(chi_tracer(1:ylat, 1:ylon, 1:ylev) *    &
                       chi_tracer(1:ylat, 1:ylon, 1:ylev))
PRINT *, 'RHS = ', SUM(chi(1:n1d) * chi_hat(1:n1d))





! =====================================================================================
! Adjoint test of cvt_v
PRINT *, '----- Adjoint test of cvt_v -----'

DEALLOCATE (chi_tracer)
ALLOCATE (chi_tracer(1:ylat, 1:ylon, 1:ylev))
ALLOCATE (field_tracer(1:ylat, 1:ylon, 1:ylev))
ALLOCATE (chi_tracer_hat(1:ylat, 1:ylon, 1:ylev))
! Put random numbers in control variable
CALL RANDOM_NUMBER (chi_tracer(1:ylat, 1:ylon, 1:ylev))

PRINT *,'Size of chi_tracer', SUM(chi_tracer(1:ylat, 1:ylon, 1:ylev) * chi_tracer(1:ylat, 1:ylon, 1:ylev))

! Input into cvt_v
CALL cvt_v (chi_tracer(1:ylat, 1:ylon, 1:ylev),              &
            field_tracer(1:ylat, 1:ylon, 1:ylev))
! Input into cvt_v_adj
CALL cvt_v_adj (chi_tracer_hat(1:ylat, 1:ylon, 1:ylev),      &
                field_tracer(1:ylat, 1:ylon, 1:ylev))

PRINT *, 'LHS = ', SUM(field_tracer(1:ylat, 1:ylon, 1:ylev) *  &
                       field_tracer(1:ylat, 1:ylon, 1:ylev))
PRINT *, 'RHS = ', SUM(chi_tracer_hat(1:ylat, 1:ylon, 1:ylev) *  &
                       chi_tracer(1:ylat, 1:ylon, 1:ylev))




! =====================================================================================
! Run inverse tests for cvt_v
! Original field: chi_tracer
! After being put through cvt_v: field_tracer
! Use chi_tracer_hat as convenient space to hold result of inverse transform
PRINT *, '----- Inverse tests for cvt_v -----'
CALL cvt_v_inv (chi_tracer_hat(1:ylat, 1:ylon, 1:ylev),      &
                field_tracer(1:ylat, 1:ylon, 1:ylev))

CALL Write_one_field ('Vert_inverse_test_in.nc',             &
                      ylat, ylon, ylev,                      &
                      chi_tracer(1:ylat, 1:ylon, 1:ylev),    &
                      'Vert_inverse_test_in',                &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),   &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),  &
                      CVT_VertCors % alts(1:ylev))

CALL Write_one_field ('Vert_inverse_test_middle.nc',         &
                      ylat, ylon, ylev,                      &
                      field_tracer(1:ylat, 1:ylon, 1:ylev),  &
                      'Vert_inverse_test_middle',            &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),   &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),  &
                      CVT_VertCors % alts(1:ylev))

CALL Write_one_field ('Vert_inverse_test_out.nc',            &
                      ylat, ylon, ylev,                      &
                      chi_tracer_hat(1:ylat, 1:ylon, 1:ylev),&
                      'Vert_inverse_test_out',               &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),   &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),  &
                      CVT_VertCors % alts(1:ylev))

PRINT *, 'RMS difference beteen input and output for tracer'
PRINT *, SQRT(SUM((chi_tracer(1:ylat, 1:ylon, 1:ylev) - chi_tracer_hat(1:ylat, 1:ylon, 1:ylev)) * &
                  (chi_tracer(1:ylat, 1:ylon, 1:ylev) - chi_tracer_hat(1:ylat, 1:ylon, 1:ylev))))




! =====================================================================================
! Adjoint test of cvt_t
PRINT *, '----- Adjoint test of cvt_t -----'
DEALLOCATE (chi_flux)
ALLOCATE (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ALLOCATE (field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ALLOCATE (chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
! Put random numbers in control variable
CALL RANDOM_NUMBER (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
! Input into cvt_t
CALL cvt_t (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),         &
            field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
! Input into cvt_t_adj
CALL cvt_t_adj (chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))

PRINT *, 'LHS = ', SUM(field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) *  &
                       field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
PRINT *, 'RHS = ', SUM(chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields) *  &
                       chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))


! =====================================================================================
! Run inverse tests for cvt_t
PRINT *, '----- Inverse tests for cvt_t -----'
! Original field: chi_flux
! After being put through cvt_t: field_flux
! Use chi_flux_hat as convenient space to hold result of inverse transform
CALL cvt_t_inv (chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields),        &
               field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))

CALL Write_one_field ('Temporal_inverse_test_in.nc',           &
                      ylat, ylon, nmonth,                      &
                      chi_flux(1:ylat, 1:ylon, 1:nmonth, 1),   &
                      'Temporal_inverse_test_in',              &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),     &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),    &
                      CVT_TemporalCors % months(1:nmonth))

CALL Write_one_field ('Temporal_inverse_test_middle.nc',       &
                      ylat, ylon, nmonth,                      &
                      field_flux(1:ylat, 1:ylon, 1:nmonth, 1), &
                      'Temporal_inverse_test_middle',          &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),     &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),    &
                      CVT_TemporalCors % months(1:nmonth))

CALL Write_one_field ('Temporal_inverse_test_out.nc',          &
                      ylat, ylon, nmonth,                      &
                      chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1),&
                      'Temporal_inverse_test_out',             &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),     &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),    &
                      CVT_TemporalCors % months(1:nmonth))

PRINT *, 'RMS difference beteen input and output for flux'
PRINT *, SQRT(SUM((chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) - chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)) * &
                  (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) - chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields))))





! =====================================================================================
! Run adjoint test for the whole transform
PRINT *, '----- Adjoint test of cvt_total routine -----'
ALLOCATE (zvar1d(1:dim_1d))
! Put random numbers in control variable
CALL RANDOM_NUMBER (chi(1:n1d))
chi_hat(1:n1d) = 0.0
! Call the total CVT
CALL cvt_total (chi(1:n1d),             &
                zvar1d(1:dim_1d),       &
                .TRUE., .TRUE.)
! Call the adjoint of this
CALL cvt_total_adj (chi_hat(1:n1d),     &
                    zvar1d(1:dim_1d),   &
                    .TRUE., .TRUE.)
PRINT *, 'LHS = ', SUM(zvar1d(1:dim_1d) * zvar1d(1:dim_1d))
PRINT *, 'RHS = ', SUM(chi(1:n1d) * chi_hat(1:n1d))




! =====================================================================================
! Run inverse test of the ChangeHorizonalGrid
PRINT *, '----- Inverse tests for ChangeHorizonalGrid -----'
! First create functions in real space
DO i = 1, ylat
  DO j = 1, ylon
    field_tracer(i,j,1) = SIN(REAL(j)*6.0*pi/REAL(ylon)) + SIN(REAL(i)*12.0*pi/REAL(ylat))
                          !REAL(lev)/REAL(ylev)
  END DO
END DO

!ALLOCATE (field_flux_hat(1:ylat, 1:ylon, 1:nmonth))
ALLOCATE (field_flux_hat(1:L+1, 1:2*L+1, 1:nmonth, 1))
ALLOCATE (field_tracer_hat(1:ylat, 1:ylon, 1:ylev))

! Call the inverse horizontal transform
CALL ChangeHorizGrid_tomcat2shtools (field_tracer(1:ylat, 1:ylon, 1:1),   &
                                     field_flux_hat(1:L+1, 1:2*L+1, 1:1, 1))
! Go back
CALL ChangeHorizGrid_shtools2tomcat (field_flux_hat(1:L+1, 1:2*L+1, 1:1, 1), &
                                     field_tracer_hat(1:ylat, 1:ylon, 1))
CALL Write_one_field ('Horiz_inverse_test_in.nc',         &
                      ylat, ylon, 1,                       &
                      field_tracer(1:ylat, 1:ylon, 1:1),   &
                      'Horiz_inverse_test_in',             &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),&
                      CVT_HorizCors % TOMCAT_longs(1:ylon),&
                      CVT_VertCors % alts(1:1))

CALL Write_one_field ('Horiz_inverse_test_mid.nc',        &
                      L+1, 2*L+1, 1,                       &
                      field_flux_hat(1:L+1, 1:2*L+1, 1:1, 1),&
                      'Horiz_inverse_test_mid',            &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),&
                      CVT_HorizCors % TOMCAT_longs(1:ylon),&
                      CVT_VertCors % alts(1:1))

CALL Write_one_field ('Horiz_inverse_test_out.nc',        &
                     ylat, ylon, 1,                       &
                     field_tracer_hat(1:ylat, 1:ylon, 1:1),&
                     'Horiz_inverse_test_out',            &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),&
                      CVT_HorizCors % TOMCAT_longs(1:ylon),&
                      CVT_VertCors % alts(1:1))

PRINT *, 'RMS difference beteen input and output for flux'
PRINT *, SQRT(SUM((field_tracer(1:ylat, 1:ylon, 1) - field_tracer_hat(1:ylat, 1:ylon, 1)) * &
                  (field_tracer(1:ylat, 1:ylon, 1) - field_tracer_hat(1:ylat, 1:ylon, 1))))

DEALLOCATE (field_flux_hat, field_tracer_hat)




! =====================================================================================
! Inverse tests for spherical transform
PRINT *, '----- Inverse tests for spherical transform -----'
! Set-up the input function
DO i = 1, L+1
  DO j = 1, 2*L+1
    data_realspace(i,j) = SIN(REAL(j)*6.0*pi/REAL(ylon)) + SIN(REAL(i)*12.0*pi/REAL(ylat))
  END DO
END DO
! Convert to spherical harmonics
CALL Spherical_inv (data_spherical(1:2, 0:L, 0:L),     & ! OUT
                    data_realspace(1:L+1, 1:2*L+1))      ! IN

! Convert back
CALL Spherical_forward (data_spherical(1:2, 0:L, 0:L),     & ! IN
                        data_realspace1(1:L+1, 1:2*L+1))     ! OUT

! Compute the difference
PRINT *, 'rms difference: ',                          &
         SUM((data_realspace1(1:L+1, 1:2*L+1) - data_realspace(1:L+1, 1:2*L+1)) * &
             (data_realspace1(1:L+1, 1:2*L+1) - data_realspace(1:L+1, 1:2*L+1)))
! Output the input and output fields
CALL Write_one_field ('shtools_in.nc',                &
                      L+1, 2*L+1, 1,                  &
                      data_realspace(1:L+1, 1:2*L+1), &
                      'shtools_in',                   &
                      CVT_HorizCors % SHtools_lats(1:L+1),&
                      CVT_HorizCors % SHtools_longs(1:2*L+1),&
                      CVT_VertCors % alts(1:1))

CALL Write_one_field ('shtools_out.nc',               &
                      L+1, 2*L+1, 1,                  &
                      data_realspace1(1:L+1, 1:2*L+1),&
                      'shtools_out',                  &
                      CVT_HorizCors % SHtools_lats(1:L+1),&
                      CVT_HorizCors % SHtools_longs(1:2*L+1),&
                      CVT_VertCors % alts(1:1))





! =====================================================================================
! Run inverse test of the horizontal transform
PRINT *, '----- Inverse tests for cvt_h -----'
! First create functions in real space
DO i = 1, ylat
  DO j = 1, ylon
    DO lev = 1, ylev
      field_tracer(i,j,lev) = SIN(REAL(j)*6.0*pi/REAL(ylon)) + SIN(REAL(i)*12.0*pi/REAL(ylat)) + &
                              REAL(lev)/REAL(ylev)
    END DO
  END DO
END DO
DO f = 1, nfields
  DO i = 1, ylat
    DO j = 1, ylon
      DO mon = 1, nmonth
        field_flux(i,j,mon,f) = (SIN(REAL(j)*24.0*pi/REAL(ylon)) + REAL((i+8)*(i+11))/REAL(ylat*ylat)) * &
                                REAL(lev)/REAL(ylev)
      END DO
    END DO
  END DO
END DO

ALLOCATE (field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
ALLOCATE (field_tracer_hat(1:ylat, 1:ylon, 1:ylev))

! Call the inverse horizontal transform
CALL cvt_h_inv (chi(1:n1d),                                       &
               field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),   &
               field_tracer(1:ylat, 1:ylon, 1:ylev),              &
               .TRUE., .TRUE.)
! Call the forward horizontal transform (use the adj arrays as convenient structures)
CALL cvt_h (chi(1:n1d),                                           &
           field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields),   &
           field_tracer_hat(1:ylat, 1:ylon, 1:ylev),              &
           .TRUE., .TRUE.)


PRINT *, 'RMS difference beteen input and output for flux'
PRINT *, SQRT(SUM((field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) - field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)) * &
                  (field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) - field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields))))
PRINT *, 'RMS difference beteen input and output for tracer'
PRINT *, SQRT(SUM((field_tracer(1:ylat, 1:ylon, 1:ylev) - field_tracer_hat(1:ylat, 1:ylon, 1:ylev)) * &
                  (field_tracer(1:ylat, 1:ylon, 1:ylev) - field_tracer_hat(1:ylat, 1:ylon, 1:ylev))))

CALL Write_one_field ('Horiz_inverse_test_flux_in.nc',         &
                      ylat, ylon, nmonth,                      &
                      field_flux(1:ylat, 1:ylon, 1:nmonth, 1), &
                      'Horiz_inverse_test_flux_in',            &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),     &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),    &
                      CVT_TemporalCors % months(1:nmonth))

CALL Write_one_field ('Horiz_inverse_test_tracer_in.nc',       &
                      ylat, ylon, ylev,                        &
                      field_tracer(1:ylat, 1:ylon, 1:ylev),    &
                      'Horiz_inverse_test_tracer_in',          &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),     &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),    &
                      CVT_VertCors % alts(1:ylev))

CALL Write_one_field ('Horiz_inverse_test_flux_out.nc',        &
                      ylat, ylon, nmonth,                      &
                      field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1),&
                      'Horiz_inverse_test_flux_out',           &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),     &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),    &
                      CVT_TemporalCors % months(1:nmonth))

CALL Write_one_field ('Horiz_inverse_test_tracer_out.nc',      &
                      ylat, ylon, ylev,                        &
                      field_tracer_hat(1:ylat, 1:ylon, 1:ylev),&
                      'Horiz_inverse_test_tracer_out',         &
                      CVT_HorizCors % TOMCAT_lats(1:ylat),     &
                      CVT_HorizCors % TOMCAT_longs(1:ylon),    &
                      CVT_VertCors % alts(1:ylev))



! =====================================================================================
! Run adjoint test for the inverse spherical transform
PRINT *, '----- Adjoint test of inverse spherical transform -----'
CALL RANDOM_NUMBER (data_realspace(1:L+1, 1:2*L+1))
! Put field through SHExpandGLQ
CALL Spherical_inv (data_spherical(1:2, 0:L, 0:L),     & ! OUT
                    data_realspace(1:L+1, 1:2*L+1))      ! IN
PRINT *, 'LHS = ', SUM(data_spherical(1:2, 0:L, 0:L) * data_spherical(1:2, 0:L, 0:L))

! Put field through MakeGridGLQ (and scaling)
CALL Spherical_inv_adj (data_spherical(1:2, 0:L, 0:L), & ! IN
                    data_realspace1(1:L+1, 1:2*L+1))     ! OUT
PRINT *, 'RHS = ', SUM(data_realspace(1:L+1, 1:2*L+1) * data_realspace1(1:L+1, 1:2*L+1))




! =====================================================================================
! Run adjoint test for the cvt_h_inv transform
PRINT *, '----- Adjoint test of cvt_h_inv routine -----'
! Put random numbers into the input of the inverse transform
CALL RANDOM_NUMBER (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
CALL RANDOM_NUMBER (chi_tracer(1:ylat, 1:ylon, 1:ylev))
! Input into cvt_h_inv
CALL cvt_h_inv (chi(1:n1d),                                        &  ! out
                chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),     &  ! in
                chi_tracer(1:ylat, 1:ylon, 1:ylev),                &  ! in
                .TRUE., .TRUE.)
! Input chi into cvt_h_inv_adj
CALL cvt_h_inv_adj (chi(1:n1d),                                        &
                    chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields), &
                    chi_tracer_hat(1:ylat, 1:ylon, 1:ylev),            &
                    .TRUE., .TRUE.)
PRINT *, 'LHS = ', SUM(chi(1:n1d) * chi(1:n1d))
PRINT *, 'RHS = ', SUM(chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) *      &
                       chi_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields)) + &
                   SUM(chi_tracer(1:ylat, 1:ylon, 1:ylev) *                 &
                       chi_tracer_hat(1:ylat, 1:ylon, 1:ylev))





! =====================================================================================
! Run adjoint test for the cvt_v_inv transform
PRINT *, '----- Adjoint test of cvt_v_inv routine -----'
! Put random numbers into the input of the inverse transform
CALL RANDOM_NUMBER (field_tracer(1:ylat, 1:ylon, 1:ylev))
! Input into cvt_v_inv
CALL cvt_v_inv (chi_tracer(1:ylat, 1:ylon, 1:ylev),    &
                field_tracer(1:ylat, 1:ylon, 1:ylev) )
! Input chi_tracer into cvt_v_inv_adj
CALL cvt_v_inv_adj (chi_tracer(1:ylat, 1:ylon, 1:ylev),    &
                    field_tracer_hat(1:ylat, 1:ylon, 1:ylev) )
PRINT *, 'LHS = ', SUM(chi_tracer(1:ylat, 1:ylon, 1:ylev) * chi_tracer(1:ylat, 1:ylon, 1:ylev))
PRINT *, 'RHS = ', SUM(field_tracer(1:ylat, 1:ylon, 1:ylev) * field_tracer_hat(1:ylat, 1:ylon, 1:ylev))





! =====================================================================================
! Run adjoint test for the cvt_t_inv transform
PRINT *, '----- Adjoint test of cvt_t_inv routine -----'
! Put random numbers into the input of the inverse transform
CALL RANDOM_NUMBER (field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
! Input into cvt_t_inv
CALL cvt_t_inv (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),    &
                field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) )
! Input chi_flux into cvt_t_inv_adj
CALL cvt_t_inv_adj (chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields),    &
                    field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields) )
PRINT *, 'LHS = ', SUM(chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) * &
                       chi_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields))
PRINT *, 'RHS = ', SUM(field_flux(1:ylat, 1:ylon, 1:nmonth, 1:nfields) * &
                       field_flux_hat(1:ylat, 1:ylon, 1:nmonth, 1:nfields))





! =====================================================================================
! Run adjoint test for the total inverse transform
PRINT *, '----- Adjoint test of cvt_total_inv routine -----'
ALLOCATE (zvar1d_hat(1:dim_1d))
! Put random numbers into the input of the inverse transform
CALL RANDOM_NUMBER (zvar1d(1:dim_1d))
! Call the total inverse CVT
CALL cvt_total_inv (chi(1:n1d),           &
                    zvar1d(1:dim_1d),     &
                    .TRUE., .TRUE.)
! Call the adjoint of this
zvar1d_hat(1:dim_1d) = 0.0
CALL cvt_total_inv_adj (chi(1:n1d),           &
                        zvar1d_hat(1:dim_1d), &
                        .TRUE., .TRUE.)
PRINT *, 'LHS = ', SUM(chi(1:n1d) * chi(1:n1d))
PRINT *, 'RHS = ', SUM(zvar1d(1:dim_1d) * zvar1d_hat(1:dim_1d))




! =====================================================================================

CALL DeAllocateCVT(.FALSE.)
CALL UnInitializeHorizTrans ()

DEALLOCATE (field_flux_hat, field_tracer_hat)
DEALLOCATE (zvar1d, zvar1d_hat)
DEALLOCATE(chi, chi_flux, chi_tracer, chi_hat, field_tracer, chi_tracer_hat, field_flux, chi_flux_hat)

END PROGRAM AdjointTestsMain
