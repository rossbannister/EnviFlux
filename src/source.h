   #include <stdlib.h>
   #include <fftw3.h>
   #include <stdio.h>
   #include <math.h>
   #include <string.h>
   #include <netcdf.h>



/* -------------------------------------------------------------------------------
   Global constants
   ------------------------------------------------------------------------------- */
  const double  pi               = 3.14159265;
  const double  Re               = 6371000.0;
  const double  rad2deg          = 180.0 / pi;
  const double  deg2rad          = pi / 180.0;
  const int     L                = 32;
  const int     nlev             = 56;
  const int     ntimes_major     = 37;
  const char    ALPfilename[256] = "/home/ross/DataAssim/SourceSink/ToyModel/3D/data/cvt_ass_legendre_poly.dat";
  const int     ecmwf_nlon       = 1440;
  const int     ecmwf_nlat       = 721;
  const int     ecmwf_nlev       = 37;
  const double  H                = 7840.0;       // Scale height
  const double  g                = 9.8;          // Gravity acceleration
  const double  rho0             = 1.225;        // Density of air at sea level
  const double  m2deg            = rad2deg / Re; // metres to degrees
  const int     cubic_parts      = 3;            // Used for cubic interpolation
  const double  nearly1          = 0.99999999;   // Nearly unity
  const double  notdef           = -9999.0;      // Not defined
  const double  mu               = 0.01;         // For conjugate gradient algorithm
  const char    verbose_output   = 'v';          // (v)erbose or nearly (s)ilent
                                                 // Implemented in PenAndGrad and SemiLagrangian_1step
  const double  siderial_day     = 24.0 * 60.0 - 3.0 - 55.91/60.0;
                                                 // Earth's siderial day period in minutes


/* -------------------------------------------------------------------------------
   Data structures
   ------------------------------------------------------------------------------- */

// Convenient structure to hold the meta data
struct metadata_type
{ int    nlon = 0;                       // Number of longitudes
  int    nlat = 0;                       // Number of latitudes
  int    nlev = 0;                       // Number of levels
  int    ntimes_major = 0;               // Number of major time steps
  int    ntimes_minor = 0;               // Number of minor time steps
  int    nss = 0;                        // Number of source/sink fields
  int    L = 0;                          // Maximum total wavenumber
  double *times = NULL;                  // Times
  double *longitude = NULL;              // Longitude values
  double *latitude = NULL;               // Latitude values
  double *level = NULL;                  // Vertical level values
  double *cos_u_lat = NULL;              // Cosines of u latitudes
  double *cos_v_lat = NULL;              // Cosines of v latitudes
};

// Model field type
struct instant_tracer_type
{ int    timestep_major = 0;             // Major time step of this field
  int    timestep_minor = 0;             // Minor time step of this field
  double time = 0.0;                     // Time label for this field
  double ***tracer = NULL;               // tracer concentration
};

// Full field or increment type (initial conditions and source values)
// Can exist in physical space and/or control space (depending on which are allocated)
struct state_type
{ int    nlon = 0;                       // Number of longitudes
  int    nlat = 0;                       // Number of latitudes
  int    nlev = 0;                       // Number of levels
  int    nss = 0;                        // Number of source/sink fields
  int    L = 0;                          // Maximum total wavenumber
  char   horiz_repres;                   // 'r' (real space) or 's' (spectral space)
  char   vert_repres;                    // 'r' (real space) or 'm' (modal space)
  char   temp_repres;                    // 'r' (real space) or 'm' (modal space)
  double *times = NULL;                  // Times
  double *longitude = NULL;              // Longitude values
  double *latitude = NULL;               // Latitude values
  double *level = NULL;                  // Vertical level values
  double ***tracer0_rs = NULL;           // tracer concentration t=0 (real space: lon, lat, lev)
  double ***source_rs = NULL;            // source/sink values (real space: lon, lat, time)
  double ****tracer0_ss = NULL;          // tracer concentration t=0 (spec space: cs, l, m, z/vmode)
  double ****source_ss = NULL;           // source/sink values (spec space: cs, l, m, t/tmode)
};


// Structure holding generic data for horizontal transform
struct HorizTransData_type
{ int    L = 0;                          // Maximum total wavenumber
  int    nlon = 0;                       // Number of longitudes
  int    nlat = 0;                       // Number of latitudes
  int    ALPsize = 0;                    // For shtools
  double **assocLegPoly = NULL;          // Associated Legendre polynomials (total wn, lat)
  double *GaussianCosCoLats = NULL;      // Cosine of Gaussian co-latitudes
  double *GaussianCoLats = NULL;         // Actual co-latitudes (degrees)
  double *GaussianLats = NULL;           // Actual latitudes (degrees)
  double *GaussianWts = NULL;            // Gaussian weights
  int    **Plm_index = NULL;             // Mapping from l,m to CVT_assocLegPoly
  double *lons = NULL;                   // Longitudes
  double       *fft_input_real = NULL;   // For use with FFTw library
  double       *fft_output_real = NULL;  // For use with FFTw library
  fftw_complex *fft_input_spec = NULL;   // For use with FFTw library
  fftw_complex *fft_output_spec = NULL;  // For use with FFTw library
  fftw_plan    plan_dft_r2c_1d;          // For use with FFTw library
  fftw_plan    plan_dft_c2r_1d;          // For use with FFTw library
};

// Structure holding data for cvt
struct CVTData_type
{ int    L = 0;                          // Maximum total wavenumber
  int    nlon = 0;                       // Number of longitudes
  int    nlat = 0;                       // Number of latitudes
  int    nlev = 0;                       // Number of levels
  int    nss = 0;                        // Number of source/sink fields
  double *times = NULL;                  // Times of source/sink fields
  double *longitude = NULL;              // Longitude values
  double *latitude = NULL;               // Latitude values
  double *level = NULL;                  // Vertical level values
  double *tracer_stddev = NULL;          // Standard deviation of tracer (fn of level)
  double **tracer_hspec = NULL;          // Horizontal stddev spectrum for tracer (fn of wn, level)
  double **tracer_vertmodes = NULL;      // Vertical modes of the tracer (global)
  double *tracer_vspec = NULL;           // Vertical stddev spectrum for tracer (fn of vert mode)
  double ***source_stddev;               // Standard deviation of source (fn of lon, lat, time)
  double *source_hspec = NULL;           // Horizontal stddev spectrum for source (fn of wn)
  bool   temporal_covs;                  // temporal correlations switch
  double **source_tempmodes = NULL;      // Temporal modes for source
  double *source_tspec = NULL;           // Temporal stddev spectrum for source (fn of temporal mode)
};

// Wind-field (longitude, latitude, level)
struct Wind_type
{ int    nlon = 0;                 // Number of longitudes
  int    nlat = 0;                 // Number of latitudes
  int    nlev = 0;                 // Number of levels
  double *longitude_u = NULL;      // Longitude values
  double *longitude_v = NULL;      // Longitude values
  double *latitude_u = NULL;       // Latitude values
  double *latitude_v = NULL;       // Latitude values
  double *level_uv = NULL;         // Level values on u/v levels (heights)
  double *level_w = NULL;          // Level values on w levels (heights)
  double ***u = NULL;              // zonal wind field
  double ***v = NULL;              // meridonal wind field
  double ***w = NULL;              // vertical wind field
};

// Data needed to make user-defined fields
struct blob_type
{ double           longitude;    // longitude position of blob (deg)
  double           latitude;     // latitude position of blob (deg)
  double           level;        // height of blob (m), 0.0 for 2D blob
  double           amplitude;    // tracer units or tracer units/s
  double           size_h;       // horizontal size of blob (grid points)
  double           size_v;       // vertical size of blob (levels)
  struct blob_type *next = NULL; // pointer to next blob
};


// Trunk of observations
struct obs_trunk_type
{ double          *mass_profile = NULL;  // Mass per unit area vertical profile
  struct obs_type *first = NULL;         // Pointer to the first observation
};

// Observation network
struct obs_type
{ double longitude;
  double latitude;
  double level;
  // --- The following are to store data for spatial linear interpolation
  int    lon_index;              // Lower grid index
  int    lat_index;              // Lower grid index
  int    lev_index;              // Lower grid index
  double lon_alpha;              // Weight to lower grid index
  double lat_alpha;              // Weight to lower grid index
  double lev_alpha;              // Weight to lower grid index
  double lon_beta;               // Weight to upper grid index
  double lat_beta;               // Weight to upper grid index
  double lev_beta;               // Weight to upper grid index
  // --------------------------------------------------------
  int    day;
  int    hour;
  int    min;
  double obtime_secs;
  // --- The following are to store data for temporal linear interpolation
  int    time_index;             // set to 0 to indicate that tracer obs alpha and beta have been set
                                 // Lower grid index source timestep for flux obs
  double time_alpha;             // Weight to lower time index
  double time_beta;              // Weight to upper time index
  // --------------------------------------------------------
  char   ob_of;                  // t: ob of tracer, f: flux, x: total column
  char   obtpe;                  // i: individual ob, g: part of gridded ob, s: satellite ob
  double ob;
  double stddev;
  double variance;
  double model_ob;       
  double innov;      
  double dJo_dmodel_ob;      
  struct obs_type *next = NULL;
};


/* -------------------------------------------------------------------------------
   List of functions
   ------------------------------------------------------------------------------- */


// ========== These functions are in cvt.cpp ==========

void cvt_total (struct state_type          *specsp,     // in  spectral space (in horiz)
                struct state_type          *realsp,     // out real space (in horiz)
                struct metadata_type       *MetaData,   // in  metadata
                struct HorizTransData_type *HorizData,  // in  info about horiz transform
                struct CVTData_type        *CVTdata);   // in  info about cvt

void cvt_total_adj (struct state_type          *specsp,     // out spectral space (in horiz)
                    struct state_type          *realsp,     // in  real space (in horiz)
                    struct metadata_type       *MetaData,   // in  metadata
                    struct HorizTransData_type *HorizData,  // in  info about horiz transform
                    struct CVTData_type        *CVTdata);   // in  info about cvt

void cvt_h (struct state_type          *specsp,     // in  spectral space (in horiz)
            struct state_type          *realsp,     // out real space (in horiz)
            struct HorizTransData_type *HorizData,  // in  info about horiz transform
            struct CVTData_type        *CVTdata);   // in  info about cvt

void cvt_h_adj (struct state_type          *specsp,     // out spectral space (in horiz)
                struct state_type          *realsp,     // in  real space (in horiz)
                struct HorizTransData_type *HorizData,  // in  info about horiz transform
                struct CVTData_type        *CVTdata);   // in  info about cvt

void cvt_h_inv (struct state_type          *specsp,     // out spectral space (in horiz)
                struct state_type          *realsp,     // in  real space (in horiz)
                struct HorizTransData_type *HorizData,  // in  info about horiz transform
                struct CVTData_type        *CVTdata,    // in  info about cvt
                char                       exactpseudo);// in  'e' exact inverse
                                                        //     'p' pseudo inverse

void cvt_v (struct state_type          *output,     // out {both level space (in vert)
            struct state_type          *input,      // in  {and real space (in horiz)
            struct CVTData_type        *CVTdata);   // in  info about cvt

void cvt_v_inv (struct state_type          *output,     // out {both level space (in vert)
                struct state_type          *input,      // in  {and real space (in horiz)
                struct CVTData_type        *CVTdata);    // in  info about cvt

void cvt_t (struct state_type   *fntime,      // out source fn of time          {both real space
            struct state_type   *fntempmode,  // in  source fn of temporal mode {in horiz
            struct CVTData_type *CVTdata);    // in  info about cvt

void cvt_t_adj (struct state_type   *fntime,      // in  source fn of time          {both real space
                struct state_type   *fntempmode,  // out source fn of temporal mode {in horiz
                struct CVTData_type *CVTdata);    // in  info about cvt

void cvt_t_inv (struct state_type   *fntime,     // in  source fn of time          {both real space
                struct state_type   *fntempmode, // out source fn of temporal mode {in horiz
                struct CVTData_type *CVTdata);   // in  info about cvt

void cvt_stddev (struct state_type   *state,     // inout (real space)
                 struct CVTData_type *CVTdata);  // in  info about cvt

void cvt_stddev_inv (struct state_type   *state,     // inout (real space)
                     struct CVTData_type *CVTdata);  // in  info about cvt

void PutWhiteNoiseControlVector (struct state_type *controlsp,
                                 int               *rndseq);

// ========== These functions are in cvt_basic.cpp ==========

void InitializeHorizTrans
  ( struct HorizTransData_type *HorizData,
    const char                 ALPfilename[256] );

void DeallocateHorizTrans
  ( struct HorizTransData_type *HorizData );

void Allocate_CVT (struct CVTData_type *CVTdata,
                   int                 L,
                   int                 nlev,
                   int                 nss);

void Deallocate_CVT (struct CVTData_type *CVTdata);

void spherical_for
  ( double                     **realspace,   // out  real space (longs, lats)
    double                     ***spectral,   // in   spectral space (cs, l, m)
    struct HorizTransData_type *HorizData );  // in   Associated Legendre polynomials, etc

void spherical_adj
  ( double                     **realspace,   // in   real space (longs, lats)
    double                     ***spectral,   // out  spectral space (cs, l, m)
    struct HorizTransData_type *HorizData );  // in   Associated Legendre polynomials, etc

void spherical_inv
  ( double                     **realspace,   // in   real space (longs, lats)
    double                     ***spectral,   // out  spectral space (cs, l, m)
    struct HorizTransData_type *HorizData);   // in   Associated Legendre polynomials, etc

void spherical_inv_adj
  ( double                     **realspace,   // out  real space (longs, lats)
    double                     ***spectral,   // in   spectral space (cs, l, m)
    struct HorizTransData_type *HorizData );  // in   Associated Legendre polynomials, etc

void fft_real2spec
  ( double                     *realspace,    // in   real space (nlon)
    double                     **specspace,   // out  spectral space (2, L)
    struct HorizTransData_type *HorizData );  // in   contains data for fftw

void fft_real2spec_adj
  ( double                     *realspace,    // out  real space (nlon)
    double                     **specspace,   // in   spectral space (2, L)
    struct HorizTransData_type *HorizData );  // in   contains data for fftw

void fft_spec2real
  ( double                     *realspace,    // out  real space (nlon)
    double                     **specspace,   // in   spectral space (2, L)
    struct HorizTransData_type *HorizData );  // in   contains data for fftw

void fft_spec2real_adj
  ( double                     *realspace,    // in   real space (nlon)
    double                     **specspace,   // out  spectral space (2, L)
    struct HorizTransData_type *HorizData );  // in   contains data for fftw



// ========== These functions are in misc.cpp ==========

void Array_1d_double_create  (double **Array,
                              int    xlen);

void Array_1d_double_destroy (double **Array);

void Array_2d_double_create  (double ***Array,
                              int    xlen,
                              int    ylen);

void Array_2d_double_destroy (double ***Array,
                              int    ylen);

void Array_2d_int_create (int ***Array,
                          int    xlen,
                          int    ylen);

void Array_2d_int_destroy (int ***Array,
                           int    xlen);

void Array_3d_double_create  (double ****Array,
                              int    xlen,
                              int    ylen,
                              int    zlen);

void Array_3d_double_destroy (double ****Array,
                              int    xlen,
                              int    ylen);

void Array_3d_int_create (int ****Array,
                          int xlen,
                          int ylen,
                          int zlen);

void Array_3d_int_destroy (int ****Array,
                           int    xlen,
                           int    ylen);

void Array_4d_double_create (double *****Array,
                             int    xlen,
                             int    ylen,
                             int    zlen,
                             int    tlen);

void Array_4d_double_destroy (double *****Array,
                              int    xlen,
                              int    ylen,
                              int    zlen);

void Allocate_metadata (struct metadata_type *MetaData,
                        int                  nlon,
                        int                  nlat,
                        int                  nlev,
                        int                  nss);

void Deallocate_metadata (struct metadata_type *MetaData);

void Allocate_state (struct state_type    *state,
                     struct metadata_type *MetaData,
                     char                 horiz_repres,  // 'r' (real space) or 's' (spectral space)
                     char                 vert_repres,   // 'r' (real space) or 'm' (modal space)
                     char                 temp_repres);  // 'r' (real space) or 'm' (modal space)

void Deallocate_state (struct state_type *state);

void Allocate_instant (struct instant_tracer_type *tracer,
                       struct metadata_type       *MetaData);

void Deallocate_instant (struct instant_tracer_type *tracer,
                         struct metadata_type       *MetaData);

void Allocate_wind (struct Wind_type     *wind,
                    struct metadata_type *MetaData,
                    bool                 ecmwf);

void Deallocate_wind (struct Wind_type *wind);

void Destroy_Obs ( struct obs_type **obs );

void copy_metadata (int    nlon0,
                    int    *nlon,
                    int    nlat0,
                    int    *nlat,
                    int    nlev0,
                    int    *nlev,
                    int    nss0,
                    int    *nss,
                    int    L0,
                    int    *L,
                    double *times0,
                    double *times,
                    double *longitude0,
                    double *longitude,
                    double *latitude0,
                    double *latitude,
                    double *level0,
                    double *level,
                    double *cos_u_lat,  // If not null this array gets filled
                    double *cos_v_lat,  // If not null this array gets filled
                    bool   printdata);

void Fill_cos_lats ( double *latitude,
                     int    nlat,
                     double *cos_v_lat,
                     double *cos_u_lat,
                     bool   printdata );

double Normal ( double mean,       //in
                double stddev,     //in
                int    *rndseq );  //inout

double randomno ( int *rndseq);

double innerproduct_realsp
  ( double                     **realspace1,     // in   input state 1
    double                     **realspace2,     // in   input state 2
    struct HorizTransData_type *HorizData,       // in   contains metadata
    bool                       inc_halos );      // in   switch to include halos

double LatIntegrate
  ( double                     *latitudinal,
    struct HorizTransData_type *HorizData,
    bool                       inc_halos );

double innerproduct_spherical
  ( double                     ***specspace1,    // in   input state 1
    double                     ***specspace2,    // in   input state 2
    struct HorizTransData_type *HorizData );     // in   contains metadata

double innerproduct_general
  ( struct state_type          *state1,      // in   input state 1
    struct state_type          *state2,      // in   input state 2
    struct HorizTransData_type *HorizData,   // in   contains metadata
    bool                       inc_halos );  // in   include halos

void add_general
  ( struct state_type          *state1,      // in   input state 1
    struct state_type          *state2,      // in   input state 2
    struct state_type          *sum,         // out  sum of above
    struct HorizTransData_type *HorizData,   // in   contains metadata
    bool                       inc_halos );  // in   include halos

void add_general_fac
  ( struct state_type          *state1,      // in   input state 1
    struct state_type          *state2,      // in   input state 2
    struct state_type          *sum,         // out  state1 + factor * state2
    struct HorizTransData_type *HorizData,   // in   contains metadata
    bool                       inc_halos,    // in   include halos
    double                     factor );

void copy_general
  ( struct state_type          *state,      // in   input state
    struct state_type          *copy,       // out  copy of above
    struct HorizTransData_type *HorizData,  // in   contains metadata
    bool                       inc_halos ); // in   include halos

void multiply_general
  ( struct state_type          *state,        // inout state
    struct HorizTransData_type *HorizData,    // in    contains metadata
    bool                       inc_halos,     // in    include halos
    double                     factor_tracer, // in    multiplication factor for tracer
    double                     factor_flux);  // in    multiplication factor for flux

int Find_index_ascend ( int    Nels,       // Number of elements in Values array
                        double Values[],   // Ascending array
                        double Value );

int Find_index_descend ( int    Nels,       // Number of elements in Values array
                         double Values[],   // Descending array
                         double Value );

int Find_index_descend_mod ( int    Nels,       // Number of elements in Values array
                             double Values[],   // Descending array
                             double Value );

void construct_total_filename (char path[],
                               char name[],
                               char complete[]);

void halos_tracer (struct instant_tracer_type *tracer,
                   struct metadata_type       *MetaData);

void halos_tracer_adj (struct instant_tracer_type *tracer,
                       struct metadata_type       *MetaData);

void halos (struct state_type *field );

void halos_adj (struct state_type *field );

void maxmin_tracer (struct metadata_type       *MetaData,
                    struct instant_tracer_type *field,
                    bool                       Inc_vert,
                    double                     *max_val,
                    double                     *min_val);

int TwoToOneD (int i,
               int j,
               int Nx,
               int Ny);

void GaussEl ( double **A,  // Matrix (modified on output)
               double *x,   // Solution (output)
               double *y,   // RHS (modified on output)
               int    k );

void CalcErr (struct metadata_type *MetaData,        // Meta data
              struct state_type    *Diff,            // Difference state
              double               *dz,              // Vertical profile of layer thicknesses
              double               *density,         // Vertical profile of density
              double               *diff_tracer,     // Output: total difference for tracer (Tg)
              double               *rse_diff_tracer, // Output: root square difference for tracer (Tg)
              double               *diff_flux,       // Output: total difference for first flux field (Tg/month)
              double               *rse_diff_flux);  // Output: total root square difference for first flux field (Tg/month)

// ========== These functions are in inout.cpp ==========

void WriteSingleRealSpField
 ( int    nlon,
   int    nlat,
   double **field,
   double *lons,
   double *lats,
   char   filename[256],
   char   varname[32] );

void WriteSingleSpecSpField
 ( int    L,
   double ***field,
   char   filename[256],
   char   varname[32] );

void cvt_matrices_input (struct CVTData_type *CVTdata,
                         char                filename[256]);

void WriteStateVector
 ( struct state_type *field,
   char              filename[256] );

void ReadStateVector (struct state_type    *field,
                      struct metadata_type *MetaData,
                      bool                 Read_tracer,
                      bool                 Read_source,
                      char                 filename[256]);

bool CheckSameLengths ( char horiz_repres,      // Horizontal representation type
                        char info1, char info2, // Abbreviations of source of data
                        int  nlon1, int  nlon2, // Two versions of No. of longs
                        int  nlat1, int  nlat2, // Two versions of No. of lats
                        int  nlev1, int  nlev2, // Two versions of No. of levs
                        int  L1,    int  L2,    // Two versions of No. of wavenumbers
                        int  nss1,  int  nss2 );// Two versions of No. of source/sink fields

void ReadINVICAT (struct state_type    *field,
                  struct metadata_type *metadata,
                  bool                 Read_tracer,
                  bool                 Read_source,
                  char                 filename[256]);

void Read_ecmwf_winds (struct Wind_type *wind,
                       int              t,
                       char             filename[256]);

void Write_winds ( struct Wind_type *wind,
                   char             filename[256]);

void Read_winds (struct Wind_type     *wind,
                 struct metadata_type *MetaData,
                 double               factor_w,
                 char                 filename[256]);

void WriteTimeSeq
 ( struct metadata_type       *MetaData,
   struct instant_tracer_type *field,
   int                        ntimes_total,
   double                     deltat,
   int                        timestep,
   int                        status,       // 0=create file only,
                                            // 1=normal (write only)
   char                       filename[256]);

void Read1Time
 ( struct metadata_type       *MetaData,       // In  (to check against)
   struct instant_tracer_type *field,          // Out (assumed already allocated)
   int                        time_request,    // In  Time index requested (1 is first)
   char                       filename[256],   // In  Filename
   bool                       *ok );           // Out Successful execution flag

void WriteDeparturePoints
// Write-out semi-Lagrangian departure points for a specific time
 ( struct metadata_type *MetaData,
   double               ***dp_lon,
   double               ***dp_lat,
   double               ***dp_lev,
   int                  ***dp_index_lon,
   int                  ***dp_index_lat,
   int                  ***dp_index_lev,
   int                  nmajor,       // No of times to be output
   double               deltat,       // Separation between times (seconds, but output will be in days)
   int                  timestep,     // Timestep number
   int                  status,       // 0=create file only,
                                      // 1=normal (write only)
   char                 filename[256]);

void ReadDeparturePoints
// Read departure points for a given timestep
 ( struct metadata_type *MetaData,
   double               ***dp_lon,
   double               ***dp_lat,
   double               ***dp_lev,
   int                  ***dp_index_lon,
   int                  ***dp_index_lat,
   int                  ***dp_index_lev,
   int                  timestep,     // Timestep number
   char                 filename[256],
   bool                 *ok);

void WriteObservations ( struct metadata_type  *MetaData,
                         struct obs_trunk_type *obs_trunk,
                         char                  output_obs_file[256] );

void ReadObservations ( struct metadata_type  *MetaData,
                        struct obs_trunk_type *obs_trunk,
                        char                  input_obs_file[256] );

void Read_INVICAT_all (struct state_type *state,
                       char              filename[256]);


// ========== These functions are in program_arguments.cpp ==========

int MakeWinds_arguments (
       int           argument_count,       // in
       char          **argument_list,      // in
       char          CVTfilename[256],     // out CVT filename for meta data
       char          wind_file_list[256],  // out
       char          red_res_dir[256] );   // out

int MakeFields_arguments (
       int              argument_count,       // in
       char             **argument_list,      // in
       char             CVTfilename[256],     // out CVT filename for meta data
       char             init_cond_file[256],  // out init cond filename
       char             *init_cond_file_type, // out 'e'nvi-flux or 'i'nvicat
       double           *factor_ic,           // out Multiplication factor for initial condition
       double           *factor_flux,         // out Multiplication factor for flux
       int              *num_init_cond_blobs, // out number of init cond blobs
       struct blob_type **ic_blob1,           // out linked list of blobs for init conds
       int              *num_ss,              // out number of surface flux fields
       double           *sDt,                 // out Source timestep (seconds, converted from days)
       int              **num_ss_blobs,       // out number of surface flux blobs in each field
       struct blob_type ***ss_field1,         // out array of linked lists of blobs for init conds
       double           *min_flux_set,        // out Sets min fluxes to this value (preserve sign)
       char             output_file[256] );   // out output filename

void Destroy_MakeFields_arguments (
       struct blob_type **ic_blob1,           // out linked list of blobs for init conds
       int              num_ss,               // in  number of surface flux fields
       int              **num_ss_blobs,       // out number of surface flux blobs in each field
       struct blob_type ***ss_field1 );        // out array of linked lists of blobs for init conds

int RunTraj_arguments (
       int     argument_count,         // in
       char    **argument_list,        // in
       double  *Dt,                    // out Major timestep (between wind files, seconds)
       double  *dt,                    // out Minor timestep (for integration scheme, seconds)
       double  *w_factor,              // out Multiplication factor of vertical winds
       double  *runlength,             // out Run length (days)
       char    wind_dir[256],          // out Wind directory
       char    output_file[256]);      // out Output filename

int SemiLagrangian_arguments (
       int     argument_count,         // in
       char    **argument_list,        // in
       char    state_ic[256],          // out Filename of state (initial condition and source)
       char    wind_dir[256],          // out Directory containing driver winds
       double  *Dt,                    // out Major timestep (between wind files, seconds)
       double  *dt,                    // out Minor timestep (for integration scheme, seconds)
       double  *kappa_dt,              // out Timestep for diffusion
       double  *kappa_h,               // out Horizontal diffusion coefficient
       double  *kappa_v,               // out Vertical diffusion coefficient
       bool    *inc_vert,              // out Include vertical transport?
       double  *factor_w,              // out Multiplication factor of vertical winds
       double  *runlength,             // out Run length (seconds, converted from days)
       double  *output_freq,           // out Time separation between outputs (seconds)
       char    *interpolation_lc,      // out Linear or cubic interpolation
       char    output_file_anim[256],  // out Output filename
       char    output_file_diags[256], // out Output filename (diagnostics)
       char    file_dps[256],          // out File for departure points
       char    *dp_file_rw );          // out (r)ead or (w) departure points for forward model

int ReplaceFirstField_arguments (
       int              argument_count,            // in
       char             **argument_list,           // in
       char             State_in_filename[256],    // out Input state filename
       char             Forecast_in_filename[256], // out Input forecast filename
       int              *FcTime,                   // out Chosen time index to be extracted from fc
       char             State_out_filename[256] ); // out Output state filename

int GenerateObs_arguments (
       int             argument_count,    // in
       char            **argument_list,   // in
       char            state_ic[256],     // out Filename of state (initial condition and source)
       char            wind_dir[256],     // out Directory containing driver winds
       char            obs_file[256],     // out Output file of observations
       char            file_dps[256],     // out File for departure points
       double          *Dt,               // out Major timestep (between wind files, seconds)
       double          *dt,               // out Minor timestep (for integration scheme, seconds)
       double          *kappa_dt,         // out Timestep for diffusion
       double          *kappa_h,          // out Horizontal diffusion coefficient
       double          *kappa_v,          // out Vertical diffusion coefficient
       bool            *inc_vert,         // out Include vertical transport?
       double          *factor_w,         // out Multiplication factor of vertical winds
       double          *output_freq,      // out Time separation between outputs (seconds)
       char            *interpolation_lc, // out Linear or cubic interpolation
       char            *dp_file_rw,       // out (r)ead or (w) departure points for forward model
       double          *tracer_ob_bias,   // out bias added to observations of the tracer
       double          *flux_ob_bias,     // out bias added to observations of the flux
       double          *tc_ob_bias,       // out bias added to observation of total column tracer
       struct obs_type **obs );           // out observation network

void Secs_to_ddhhmm (double time_s,
                     int    *days,
                     int    *hrs,
                     int    *mins);

void xy_rotate( double v_in[3],
                double v_out[3],
                double angle_rad );

void xz_rotate( double v_in[3],
                double v_out[3],
                double angle_rad );

int MakeBg_arguments (
       int    argument_count,     // in
       char   **argument_list,    // in
       char   truth_in[256],      // out Filename of truth state (initial condition and source)
       char   bg_out[256],        // out Filename of background state (initial condition and source)
       char   bgpert_out[256],    // out Filename of background perturbation (initial condition and source)
       char   CVTfilename[256],   // out CVT filename
       double *factor_std_tracer, // out To multiply tracer err std
       double *factor_std_flux,   // out To multiply flux err std
       double *factor_pert_tracer,// out To multiply tracer bg pert
       double *factor_pert_flux); // out To multiply flux bg pert

int Assimilate_arguments (
       int    argument_count,           // in
       char   **argument_list,          // in
       char   bg_filename[256],         // out Filename of background state
       char   wind_dir[256],            // out Directory containing driver winds
       char   obsfile_in[256],          // out Input file for observations
       char   file_dps[256],            // out Input or output file for departure points
       char   CVTfilename[256],         // out CVT filename (specification of B)
       double *factor_std_tracer,       // out To multiply tracer err std
       double *factor_std_flux,         // out To multiply flux err std
       char   obsfile_out_bg[256],      // out Output file of observation information at the background
       char   anal_filename[256],       // out Filename of analysis state
       char   analinc_filename[256],    // out Filename of analysis increment
       char   Assim_diags_filename[256],// out Filename of assimilation diagnostics
       char   obsfile_out_anal[256],    // out Output file of observation information at the analysis
       double *Dt,                      // out Major timestep (between wind files, seconds)
       double *dt,                      // out Minor timestep (for integration scheme, seconds)
       double *kappa_dt,                // out Timestep for diffusion
       double *kappa_h,                 // out Horizontal diffusion coefficient
       double *kappa_v,                 // out Vertical diffusion coefficient
       bool   *inc_vert,                // out Include vertical transport?
       double *factor_w,                // out Multiplication factor of vertical winds
       char   *interpolation_lc,        // out Linear or cubic interpolation
       char   *descent_algorithm,       // out Descent algorithm type
       double *conv_criterion,          // out Criterion for gradient norm for convergence
       int    *max_iters );             // out Maximum number of iterations

int MakeWinds_arguments (
       int           argument_count,       // in
       char          **argument_list,      // in
       char          State1_file[256],     // out state1 filename
       char          State2_file[256] );   // out state2 filename

int CalcError_arguments (
       int    argument_count,           // in
       char   **argument_list,          // in
       char   State_file[256],          // out Filename of state
       char   Truth_file[256],          // out Filename of truth
       char   Background_file[256],     // out Filename of background
       char   CVTfilename[256],         // out CVT filename
       char   Output_file[256] );       // out Filename of data file containing error diagnostics

int Invicat2Enviflux_arguments (
       int    argument_count,           // in
       char   **argument_list,          // in
       char   CVTfilename[256],         // out CVT filename
       char   INVICAT_file[256],        // out INVICAT filename
       char   ENVIFLUX_file[256],       // out ENVIFLUX filename for output
       double *maxflux,                 // out Fluxes will be scaled to have this max value
       double *maxtracer );             // out Tracer will be scaled to have this max value



// ========== These functions are in transport_model.cpp ==========

void model_1step ( struct metadata_type       *MetaData,         // Meta data
                   struct state_type          *state,            // State vector
                   struct instant_tracer_type *state_prev,       // Tracer field at t
                   struct instant_tracer_type *state_next,       // Tracer field extrapolated to t+1
                   int                        minor_t,           // Minor timestep No.
                   double                     dt,                // Minor timestep size
                   int                        min_steps_per_maj, // No. of minor steps per major step
                   double                     kappa_h,           // Horiz diffusion coefficient
                   double                     kappa_v,           // Vert diffusion coefficient
                   double                     explicit_adv,      // Explicitness coefficient
                   double                     explicit_diff,     // Explicitness coefficient
                   bool                       upstream,          // True for upstream scheme for explicit part,
                                                                 // False for centred in space for explicit part
                   bool                       Inc_vert,          // True to include vertical transport
                                                                 // False only level 1
                   int                        diff_grid,         // How the diffusion coefficient depends on grid size
                                                                 // 1: ~1/(grid size)
                                                                 // 2: ~sqrt(1/(grid size))
                                                                 // 3: no grid size dependency
                                                                 // 4: ~grid size
                                                                 // 5: ~sqrt(grid size)
                   struct Wind_type           *windA,            // Wind fields at lower major timestep
                   struct Wind_type           *windB,            // Wind fields at next major step
                   struct Wind_type           *wind_now );       // Workspace for interpolated winds


double determine_kappa ( int    diff_grid,
                         double kappa_const,
                         double Delta,
                         double Delta_ref );

void SemiLagrangian_1step ( struct metadata_type       *MetaData,        // Meta data
                            struct state_type          *state,           // State vector
                            struct instant_tracer_type *state_prev,      // Tracer field at t-1
                            struct instant_tracer_type *state_next,      // Tracer field extrapolated to t
                            int                        major_timestep,   // Major timestep number
                            int                        sfield,           // source/sink field number
                            double                     Dt,               // Major timestep size
                            double                     dt,               // Minor timestep size
                            double                     kappa_dt,         // Timestep for diffusion
                            int                        gamma,            // Number of minor timesteps per major timestep
                            double                     kappa_h,          // Horiz diffusion coefficient
                            double                     kappa_v,          // Vert diffusion coefficient
                            bool                       adv_switch,       // Advection on/off
                            double                     fluxfactor,       // Flux factor
                            bool                       inc_vert,         // Include vertical transport
                            struct Wind_type           *windA,           // Wind fields at lower major timestep
                            struct Wind_type           *windB,           // Wind fields at next major step
                            char                       interpolate_lc,   // Linear or cubic interpolation
                            bool                       mention_dps,      // Departure points file mentioned?
                            char                       file_dps[256],    // Departure points output file
                            char                       dp_file_rw );     // (r)ead or (w)rite departure points file?


void SemiLagrangian_1step_adj ( struct metadata_type       *MetaData,         // Meta data
                                struct state_type          *state_hat,        // Out: State vector
                                struct instant_tracer_type *state_prev_hat,   // Out: Tracer field at t-1
                                struct instant_tracer_type *state_next_hat,   // In:  Tracer field extrapolated to t
                                int                        major_timestep,    // Major timestep number
                                int                        sfield,            // source/sink field number
                                double                     Dt,                // Major timestep size
                                double                     kappa_dt,          // Timestep for diffusion
                                double                     kappa_h,           // Horiz diffusion coefficient
                                double                     kappa_v,           // Vert diffusion coefficient
                                bool                       adv_switch,        // Advection on/off
                                double                     fluxfactor,        // Flux factor
                                bool                       inc_vert,          // Include vertical transport
                                char                       interpolate_lc,    // Linear or cubic interpolation
                                char                       input_file_dps[256]);// Departure points input file


// ========== These functions are in trajectory.cpp ==========

void trajectory ( int                   minor_t,           // Minor timestep No., local time (relative to windA)
                  double                dt,                // Minor timestep size
                  double                Dt,                // Major timestep size
                  struct Wind_type      *windA,            // Wind fields at start of major timestep
                  struct Wind_type      *windB,            // Wind fields at end of major timestep
                  double                lon_start,         // particle longitude at start time
                  double                lat_start,         // particle latitude at start time
                  double                lev_start,         // particle level at start time
                  double                *lon_end,          // particle longitude at end time
                  double                *lat_end,          // particle latitude at end time
                  double                *lev_end,          // particle level at end time
                  bool                  show_diags );


void RungeKutta_k ( double                dt,                // Minor timestep size
                    double                times[2],          // Time axis (relative units)
                    struct Wind_type      *windA,            // Wind fields at start of major timestep
                    struct Wind_type      *windB,            // Wind fields at end of major timestep
                    double                lon,               // particle longitude
                    double                lat,               // particle latitude
                    double                lev,               // particle level
                    double                t,                 // time
                    double                *kx,               // k value
                    double                *ky,               // k value
                    double                *kz,               // k value
                    bool                  show_diags );

void Fill_short_array ( double ***large,
                        double small[2][2][2],
                        int    index_x,
                        int    index_y,
                        int    index_z );

void Fill_short_array_adj ( double ***large,
                            double small[2][2][2],
                            int    index_x,
                            int    index_y,
                            int    index_z );

void Fill_short_dims (double *lon_u_large,
                      double *lon_v_large,
                      double *lat_u_large,
                      double *lat_v_large,
                      double *lev_uv_large,
                      double *lev_w_large,
                      double lon_u[2],
                      double lon_v[2],
                      double lat_u[2],
                      double lat_v[2],
                      double lev_uv[2],
                      double lev_w[2],
                      int    index_ux,
                      int    index_vx,
                      int    index_uy,
                      int    index_vy,
                      int    index_uvz,
                      int    index_wz );

void Fill_short_dims_cubic (int    nlon,
                            int    nlat,
                            int    nlev,
                            double *lon_large,
                            double *lat_large,
                            double *lev_large,
                            double lon_cubic[6],
                            double lat_cubic[6],
                            double lev_cubic[6],
                            int    index_lon,
                            int    index_lat,
                            int    index_lev);

void Fill_short_array_cubic ( int    nlon,
                              int    nlat,
                              int    nlev,
                              double ***large,
                              double small[6][6][6],
                              int    index_lon,
                              int    index_lat,
                              int    index_lev );

void Fill_short_array_cubic_adj ( int    nlon,
                                  int    nlat,
                                  int    nlev,
                                  double ***large,       // out
                                  double small[6][6][6], // in
                                  int    index_lon,
                                  int    index_lat,
                                  int    index_lev );

void CheckBounds ( double *lon,
                   double *lat,
                   double *lev,
                   double top );


// ========== These functions are in interpolate.cpp ==========

double Interpolate1D ( double field[2], // The 1D field to be interpolated
                       double axis[2],  // The axis values
                       double value);   // The point to interpolate to

void Interpolate1D_adj ( double field_hat[2],     // The 1D field to be interpolated (out)
                         double interpolate_hat,  // The interpolated value (in)
                         double axis[2],          // The axis values
                         double value);           // The point to interpolate to

double Interpolate2D ( double field[2][2], // The 2D field to be interpolated
                       double axisx[2],    // The axis x values
                       double axisy[2],    // The axis y values
                       double valuex,      // The x point to interpolate to
                       double valuey);     // The y point to interpolate to

void Interpolate2D_adj ( double field_hat[2][2],  // The 2D field to be interpolated (out)
                         double interpolate_hat,  // The interpolated value (in)
                         double axisx[2],         // The axis x values
                         double axisy[2],         // The axis y values
                         double valuex,           // The x point to interpolate to
                         double valuey);          // The y point to interpolate to

double Interpolate3D ( double field[2][2][2],  // The 3D field to be interpolated
                       double axisx[2],        // The axis x values
                       double axisy[2],        // The axis y values
                       double axisz[2],        // The axis z values
                       double valuex,          // The x point to interpolate to
                       double valuey,          // The y point to interpolate to
                       double valuez );        // The y point to interpolate to

void Interpolate3D_adj ( double field_hat[2][2][2],  // The 3D field to be interpolated (out)
                         double interpolate_hat,     // The interpolated value (in)
                         double axisx[2],        // The axis x values
                         double axisy[2],        // The axis y values
                         double axisz[2],        // The axis z values
                         double valuex,          // The x point to interpolate to
                         double valuey,          // The y point to interpolate to
                         double valuez );        // The y point to interpolate to

double Interpolate3Dt ( double fieldt0[2][2][2],  // The 3D field to be interpolated at time 0
                        double fieldt1[2][2][2],  // The 3D field to be interpolated at time 1
                        double axisx[2],          // The axis x values
                        double axisy[2],          // The axis y values
                        double axisz[2],          // The axis z values
                        double axist[2],          // The axis time values
                        double valuex,            // The x point to interpolate to
                        double valuey,            // The y point to interpolate to
                        double valuez,            // The z point to interpolate to
                        double valuet);           // The time point to interpolate to

void Interpolate3Dt_adj ( double fieldt0_hat[2][2][2],  // The 3D field to be interpolated at time 0 (out)
                          double fieldt1_hat[2][2][2],  // The 3D field to be interpolated at time 1 (out)
                          double interpolate_hat,       // The interpolated value (in)
                          double axisx[2],              // The axis x values
                          double axisy[2],              // The axis y values
                          double axisz[2],              // The axis z values
                          double axist[2],              // The axis time values
                          double valuex,                // The x point to interpolate to
                          double valuey,                // The y point to interpolate to
                          double valuez,                // The z point to interpolate to
                          double valuet);               // The time point to interpolate to

double interpolate_1d_cubic ( double x[6],     // grid positions
                              double vals[6],  // field values (in)
                              double pos,      // position to interpolate to
                              int    parts );  // divide domain into 2 or 3 parts

void interpolate_1d_cubic_adj ( double x[6],            // grid positions
                                double vals_hat[6],     // field values (out)
                                double interpolate_hat, // interpolated value (in)
                                double pos,             // position to interpolate to
                                int    parts );         // divide domain into 2 or 3 parts

double interpolate_2d_cubic ( double x[6],
                              double y[6],
                              double vals[6][6],
                              double posx,
                              double posy,
                              int    parts );

void interpolate_2d_cubic_adj ( double x[6],
                                double y[6],
                                double vals_hat[6][6],   // field values (out)
                                double interpolate_hat,  // interpolated value (in)
                                double posx,
                                double posy,
                                int    parts );

double interpolate_3d_cubic ( double x[6],
                              double y[6],
                              double z[6],
                              double vals[6][6][6],
                              double posx,
                              double posy,
                              double posz,
                              int    parts );

void interpolate_3d_cubic_adj ( double x[6],
                                double y[6],
                                double z[6],
                                double vals_hat[6][6][6], // field values (out)
                                double interpolate_hat,   // interpolated value (in)
                                double posx,
                                double posy,
                                double posz,
                                int    parts );

double interpolate_3d_cubic_t_linear ( double x[6],
                                       double y[6],
                                       double z[6],
                                       double t[2],
                                       double valst1[6][6][6],
                                       double valst2[6][6][6],
                                       double posx,
                                       double posy,
                                       double posz,
                                       double post,
                                       int    parts );

double interpolate_3d ( double x[2],
                        double y[2],
                        double z[2],
                        double vals[2][2][2],
                        double pos[3] );

double interpolate_2d ( double x[2],
                        double y[2],
                        double vals[2][2],
                        double pos[2] );

double interpolate_1d ( double x[2],
                        double vals[2],
                        double pos );

// ========== These functions are in ObsOperators.cpp ==========

void ObservationOperator
 ( struct obs_type            *ob,            // inout Observation and data
   struct obs_trunk_type      *ob_trunk,      // in    Containing the mass profile (for satellites)
   struct metadata_type       *MetaData,      // in    Lats, lons, etc.
   struct instant_tracer_type *field_lowert,  // in    Tracer field at lower time
   struct instant_tracer_type *field_uppert,  // in    Tracer field at upper time
   double                     lowert_tracer,  // in    Lower time of tracer field in s
   double                     uppert_tracer,  // in    Upper time of tracer field in s
   struct state_type          *state );       // in    Contains flux field

void ObservationOperator_adj
 ( struct obs_type            *ob_hat,            // in    Observation and data
   struct obs_trunk_type      *ob_trunk,          // in    Containing the mass profile (for satellites)
   struct metadata_type       *MetaData,          // in    Lats, lons, etc.
   struct instant_tracer_type *field_lowert_hat,  // inout Tracer field at lower time
   struct instant_tracer_type *field_uppert_hat,  // inout Tracer field at upper time
   struct state_type          *state_hat );       // inout Contains flux field


double ob_op_interp_tracer_3D ( struct obs_type            *ob,
                                struct instant_tracer_type *field );

void ob_op_interp_tracer_3D_adj ( struct obs_type            *ob,
                                  struct instant_tracer_type *field_hat,  // inout
                                  double                     result_hat );// in

double ob_op_interp_flux_2D ( struct obs_type *ob,
                              double          ***field,
                              int             tindex );

void ob_op_interp_flux_2D_adj ( struct obs_type *ob,
                                double          ***field_hat, // inout
                                int             tindex,
                                double          result_hat ); // in

double ob_op_interp_tracer_2D ( struct obs_type            *ob,
                                struct instant_tracer_type *field,
                                int                        lev );

void ob_op_interp_tracer_2D_adj ( struct obs_type            *ob,
                                  struct instant_tracer_type *field_hat,  // inout
                                  int                        lev,
                                  double                     result_hat );// in

// ========== These functions are in PenAndGrad.cpp ==========

void PenAndGrad
// Function to compute the penalty and gradient of the cost function
( struct state_type          *xb,             //in    background state
  struct obs_trunk_type      *obs_trunk,      //inout observation structure
  struct state_type          *chi,            //in    control variable
  bool                       zerochi,         //in    true if chi is zero
  struct metadata_type       *MetaData,       //in    metadata
  struct HorizTransData_type *HorizData,      //in    info about horiz transform
  struct CVTData_type        *CVTdata,        //in    info about cvt
  double                     *Jb,             //out   background part of cost fn
  double                     *Jo,             //out   observation part of cost fn
  double                     *J,              //out   total cost fn
  struct state_type          *grad,           //out   gradient in control space
  bool                       calc_grad,       //in    switch to compute gradient
  double                     Dt,              //in    major timestep
  double                     dt,              //in    minor timestep
  double                     kappa_dt,        //in    diffusion timestep size
  double                     kappa_h,         //in    horiz diffusion coefficient
  double                     kappa_v,         //in    vert diffusion coefficient
  bool                       inc_adv,         //in    include advection
  bool                       inc_vert,        //in    include vertical transport?
  char                       interpolation_lc,//in    interpolation type
  double                     factor_w,        //in    mult factor of vert winds
  char                       wind_dir[256],   //in    directory containing driver winds
  char                       file_dps[256],   //in    in or out file for departure points
  bool                       new_dp           //in    true if dp file needed to be made
);

