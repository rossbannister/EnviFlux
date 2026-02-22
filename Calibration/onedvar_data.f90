 module onedvar_data
        
! Description:
! Contains data for INVICAT (onedvar.f90 and simul.f90)
!
! Method :
! ------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0               Initial version (C. Wilson)
!           10/2020  Modified, Ross Bannister
!
! Code Description:
!   Language:           Fortran 90.
! Description:
! Contains data for 1dvar (onedvar.f90 and simul.f90)
!
!
USE main_data, ONLY : ylat,    &
                      ylev,    &
                      nmonth,  &
                      nfields

implicit none

integer                           :: nobs
real,allocatable, dimension(:)    :: v_obs, a1_obs
real,allocatable, dimension(:)    :: r_obs
real,allocatable, dimension(:)    :: var_1d, var_obs
real,allocatable, dimension(:)    :: xnorm, xjo, xjb



! R.N. BANNISTER: Structures required for non-diagonal B-matrix
! Data structures are for option 1 from table in NCEO Conference poster (2018)

integer, parameter                :: L = ylat               ! Number of associated Legendre polynomials
integer, parameter                :: GLatEq = L/2 + 1       ! Latitude index of the equator
integer, parameter                :: ALPsize=(L+1)*(L+2)/2  ! Size of space for associated LPs
real,    parameter                :: htrans_factor = SQRT(2. * REAL(2*L+1))
integer, parameter                :: fft_worklen = 2*(2*L+1) ! Much larger than needbe

character(len=320), parameter     :: ALPdatafile = './cvt_ass_legendre_poly.dat'
character(len=320), parameter     :: CVTinputfile = './cvt_input_file.nc'



! ##### USER INPUT #####
integer, parameter                :: nforecasts_calibtracer = 10
character(len=320), parameter     :: forecasts_calib_tracer_dir = '/home/ross/DataAssim/SourceSink/TOPCAT/METHANE_FORECASTS'
! ######################


! To do with SHtools - do not change these three parameters
integer, parameter                :: SHtools_Normalization = 1
integer, parameter                :: SHtools_ApplyCSphase = 1
integer, parameter                :: SHtools_CRnorm = 0


! ==============================================================================================
! ===== Structure concerning the spatial standard deviations for the tracer and flux ===========
TYPE StdDev_type
! ##### USER INPUT #####
  real                                  :: FracOfPrior_tracer = 0.0 ! 0.0 means determine consts
  real, allocatable, dimension(:,:,:)   :: std_tracer
  real                                  :: FracOfPrior_flux = 0.0   ! 0.0 means determine consts
  real, allocatable, dimension(:,:,:,:) :: std_flux
! ##### USER INPUT #####
END TYPE StdDev_type




! ==============================================================================================
! ===== Structure concerning details of converting from one horizontal grid to another =========
TYPE ChangeGrid_type
  integer, allocatable, dimension(:)     :: grid_change_long_lower
  real,    allocatable, dimension(:)     :: grid_change_long_lower_wt
  real,    allocatable, dimension(:)     :: grid_change_long_upper_wt
  integer, allocatable, dimension(:)     :: grid_change_lat_lower
  real,    allocatable, dimension(:)     :: grid_change_lat_lower_wt
  real,    allocatable, dimension(:)     :: grid_change_lat_upper_wt
END TYPE ChangeGrid_type




! ==============================================================================================
! ===== Structure concerning horizontal covariances for the flux and initial tracer fields =====
TYPE HorizCov_type
! ##### USER INPUT #####
  integer                                :: horiztracer_cor_tpe & ! Type of imposed horizontal correlations for tracer
                                              = 3                 ! 1 = Lorentzian
                                                                  ! 2 = Gaussian
                                                                  ! 3 = SOAR
                                                                  ! 4 = exponential
! ######################
  real,    allocatable, dimension(:)     :: lengthscale_tracer    ! Derived lengthscales of tracer
  real,    allocatable, dimension(:,:)   :: horizcors_tracer      ! Above in numbers




  integer, allocatable, dimension(:)     :: horizflux_cor_tpe     ! Type of imposed horizontal correlations for flux
                                                                  ! 1 = Lorentzian
                                                                  ! 2 = Gaussian
                                                                  ! 3 = SOAR
                                                                  ! 4 = exponential
  real,    allocatable, dimension(:,:)   :: lengthscale_flux      ! Imposed lengthscales of flux
  real,    allocatable, dimension(:,:,:) :: horizcors_flux        ! Above in numbers
  real,    allocatable, dimension(:,:,:) :: var_hspec_flux        ! Flux variance spectrum (fn of tot wn, month or timescale, field)
  real,    allocatable, dimension(:,:,:) :: std_hspec_flux        ! Square-root of above
  real,    allocatable, dimension(:,:)   :: var_hspec_tracer      ! Tracer variance spectrum (fn of tot wn, vert lev)
  real,    allocatable, dimension(:,:)   :: std_hspec_tracer      ! Square-root of above
  real,    allocatable, dimension(:,:)   :: assocLegPoly          ! Associated Legendre polynomials (total wn, lat)
  real,    allocatable, dimension(:)     :: GaussianCosCoLats     ! Cosine of Gaussian co-latitudes
  real,    allocatable, dimension(:)     :: GaussianCoLats        ! Actual co-latitudes (degrees)
  real,    allocatable, dimension(:)     :: GaussianWts           ! Gaussian weights
  real,    allocatable, dimension(:)     :: SHtools_longs         ! Longitude values on SSHtools grid
  real,    allocatable, dimension(:)     :: SHtools_lats          ! Latitude values on SSHtools grid
  real,    allocatable, dimension(:)     :: TOMCAT_longs          ! Longitude values on TOMCAT grid
  real,    allocatable, dimension(:)     :: TOMCAT_lats           ! Latitude values on TOMCAT grid
  integer, allocatable, dimension(:,:)   :: Plm_index             ! Mapping from l,m to CVT_assocLegPoly
  real,    allocatable, dimension(:)     :: fft_wsave             ! For fft
  real,    allocatable, dimension(:)     :: fft_work              ! For fft
  TYPE (ChangeGrid_type)                 :: SH_to_TOMCAT          ! To interpolate from SHtools grid to TOMCAT grid
  TYPE (ChangeGrid_type)                 :: TOMCAT_to_SH          ! To interpolate from TOMCAT grid to SHtools grid
END TYPE HorizCov_type




! ==============================================================================================
! ===== Structure concerning the vertical correlations for the initial tracer field ============
TYPE VertCov_type
  real, allocatable, dimension(:)        :: alts                  ! Global mean altitudes
! ##### USER INPUT #####
  integer                                :: vert_covs = 3         ! 1 = Vertical eigenvalues a fn of latitude and longitude
                                                                  ! 2 = Vertical eigenvalues a fn of latitude
                                                                  ! 3 = Global eigenvalues
  logical                                :: ForceVCorr = .FALSE.  ! Do adjustment to force correlation
! ######################
  real, allocatable, dimension(:,:)      :: glob_av_vert_cov      ! Global average vertical covariance matrix
  real, allocatable, dimension(:,:,:)    :: var_vspec_tracer      ! Tracer variance spectrum (fn of lat, lon, and vert index)
  real, allocatable, dimension(:,:,:)    :: std_vspec_tracer      ! Square-root of above
  real, allocatable, dimension(:,:)      :: vert_eigenvec_tracer  ! Vertical eigenvectors
  real, allocatable, dimension(:,:,:)    :: vert_adjust_tracer    ! To ensure vert transform implies a correlation
END TYPE VertCov_type




! ==============================================================================================
! ===== Structure concerning the temporal correlations for the flux field ======================
TYPE TemporalCov_type
  real,    allocatable, dimension(:)     :: months                ! Month values
  integer, allocatable, dimension(:)     :: temporal_covs         ! 0 = No temporal covariances for flux
                                                                  ! 1 = Temporal covariances, spatial scale a fn of timescale
                                                                  ! 2 = Temporal covariances, spatial scale a fn of time
  integer, allocatable, dimension(:)     :: temporal_cor_tpe      ! Form of temporal correlation
                                                                  ! 1 = Lorentzian
                                                                  ! 2 = Gaussian
                                                                  ! 3 = SOAR
                                                                  ! 4 = exponential
  real,    allocatable, dimension(:)     :: timescale_flux        ! Timescale
  real,    allocatable, dimension(:,:)   :: temporalcors_flux     ! The temporal correlations to be modelled
  real,    allocatable, dimension(:,:,:) :: temp_cor_matrix       ! The temporal correlation matrix itself (for each field)
  real,    allocatable, dimension(:,:)   :: var_tspec_flux        ! Flux variance spectrum (fn of timescale, field)
  real,    allocatable, dimension(:,:)   :: std_tspec_flux        ! Square-root of above
  real,    allocatable, dimension(:,:,:) :: temp_eigenvec_flux    ! Temporal eigenvectors (fn of time, timescale, field)
END TYPE TemporalCov_type




! ==============================================================================================
TYPE (StdDev_type)                :: CVT_std
TYPE (HorizCov_type)              :: CVT_HorizCors
TYPE (VertCov_type)               :: CVT_VertCors
TYPE (TemporalCov_type)           :: CVT_TemporalCors



! WARNING -- n1d = dim_1d if not using Ross's scheme
! Size of control vector
integer, parameter               :: n1d = (1+L)*(1+L) * (ylev + nmonth*nfields)
real, parameter                  :: zepsc = 1.e-12 ! security value to avoid div by zero



end module onedvar_data
