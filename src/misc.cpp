   #include <source.h>



// -------------------------------------------------------------------------------
void Array_1d_double_create (double **Array,
                             int    xlen)
{ int i;
  *Array = new double [xlen];
  for (i=0; i<xlen; i++)
  { (*Array)[i] = 0.0;
  }
}


// -------------------------------------------------------------------------------
void Array_1d_double_destroy (double **Array)
{ delete[] (*Array);
}


// -------------------------------------------------------------------------------
void Array_2d_double_create (double ***Array,
                             int    xlen,
                             int    ylen)
{ int i, j;
  *Array = new double* [xlen];
  for (i=0; i<xlen; i++)
  { (*Array)[i] = new double[ylen];
    for (j=0; j<ylen; j++)
    { (*Array)[i][j] = 0.0;
    }
  }
}


// -------------------------------------------------------------------------------
void Array_2d_double_destroy (double ***Array,
                              int    xlen)
{ int i;
  for (i=0; i<xlen; i++)
  { delete[] (*Array)[i];
  }
  delete[] (*Array);
}

// -------------------------------------------------------------------------------
void Array_2d_int_create (int ***Array,
                          int    xlen,
                          int    ylen)
{ int i, j;
  *Array = new int* [xlen];
  for (i=0; i<xlen; i++)
  { (*Array)[i] = new int[ylen];
    for (j=0; j<ylen; j++)
    { (*Array)[i][j] = 0;
    }
  }
}


// -------------------------------------------------------------------------------
void Array_2d_int_destroy (int ***Array,
                           int    xlen)
{ int i;
  for (i=0; i<xlen; i++)
  { delete[] (*Array)[i];
  }
  delete[] (*Array);
}


// -------------------------------------------------------------------------------
void Array_3d_double_create (double ****Array,
                             int    xlen,
                             int    ylen,
                             int    zlen)
{ int i, j, k;
  *Array = new double** [xlen];
  for (i=0; i<xlen; i++)
  { (*Array)[i] = new double* [ylen];
    for (j=0; j<ylen; j++)
    { (*Array)[i][j] = new double[zlen];
      for (k=0; k<zlen; k++)
      { (*Array)[i][j][k] = 0.0;
      }
    }
  }
}


// -------------------------------------------------------------------------------
void Array_3d_double_destroy (double ****Array,
                              int    xlen,
                              int    ylen)
{ int i, j;
  for (i=0; i<xlen; i++)
  { for (j=0; j<ylen; j++)
    { delete[] (*Array)[i][j];
    }
    delete[] (*Array)[i];
  }
  delete[] (*Array);
}


// -------------------------------------------------------------------------------
void Array_3d_int_create (int ****Array,
                          int xlen,
                          int ylen,
                          int zlen)
{ int i, j, k;
  *Array = new int** [xlen];
  for (i=0; i<xlen; i++)
  { (*Array)[i] = new int* [ylen];
    for (j=0; j<ylen; j++)
    { (*Array)[i][j] = new int[zlen];
      for (k=0; k<zlen; k++)
      { (*Array)[i][j][k] = 0.0;
      }
    }
  }
}


// -------------------------------------------------------------------------------
void Array_3d_int_destroy (int ****Array,
                           int    xlen,
                           int    ylen)
{ int i, j;
  for (i=0; i<xlen; i++)
  { for (j=0; j<ylen; j++)
    { delete[] (*Array)[i][j];
    }
    delete[] (*Array)[i];
  }
  delete[] (*Array);
}



// -------------------------------------------------------------------------------
void Array_4d_double_create (double *****Array,
                             int    xlen,
                             int    ylen,
                             int    zlen,
                             int    tlen)
// Also initialise to zero
{ int i, j, k, t;
  *Array = new double*** [xlen];
  for (i=0; i<xlen; i++)
  { (*Array)[i] = new double** [ylen];
    for (j=0; j<ylen; j++)
    { (*Array)[i][j] = new double* [zlen];
      for (k=0; k<zlen; k++)
      { (*Array)[i][j][k] = new double[tlen];
        for (t=0; t<tlen; t++)
        { (*Array)[i][j][k][t] = 0.0;
        }
      }
    }
  }
}


// -------------------------------------------------------------------------------
void Array_4d_double_destroy (double *****Array,
                              int    xlen,
                              int    ylen,
                              int    zlen)
{ int i, j, k;
  for (i=0; i<xlen; i++)
  { for (j=0; j<ylen; j++)
    { for (k=0; k<zlen; k++)
      { delete[] (*Array)[i][j][k];
      }
      delete[] (*Array)[i][j];
    }
    delete[] (*Array)[i];
  }
  delete[] (*Array);
}



// -------------------------------------------------------------------------------
void Allocate_metadata (struct metadata_type *MetaData,
                        int                  nlon,
                        int                  nlat,
                        int                  nlev,
                        int                  nss)
{ Array_1d_double_create (&((*MetaData).longitude), nlon+2);
  Array_1d_double_create (&((*MetaData).latitude),  nlat+2);
  Array_1d_double_create (&((*MetaData).level),     nlev+2);
  Array_1d_double_create (&((*MetaData).times),     nss);
  Array_1d_double_create (&((*MetaData).cos_u_lat), nlat+2);
  Array_1d_double_create (&((*MetaData).cos_v_lat), nlat+2);
}


// -------------------------------------------------------------------------------
void Deallocate_metadata (struct metadata_type *MetaData)
{ Array_1d_double_destroy (&((*MetaData).longitude));
  Array_1d_double_destroy (&((*MetaData).latitude));
  Array_1d_double_destroy (&((*MetaData).level));
  Array_1d_double_destroy (&((*MetaData).times));
  Array_1d_double_destroy (&((*MetaData).cos_u_lat));
  Array_1d_double_destroy (&((*MetaData).cos_v_lat));
}


// -------------------------------------------------------------------------------
void Allocate_state (struct state_type    *state,
                     struct metadata_type *MetaData,
                     char                 horiz_repres,  // 'r' (real space) or 's' (spectral space)
                     char                 vert_repres,   // 'r' (real space) or 'm' (modal space)
                     char                 temp_repres)   // 'r' (real space) or 'm' (modal space)
{ // Fill meta data and allocate field
  int x, y, z, t;

  // Meta data
  (*state).nlon         = (*MetaData).nlon;
  (*state).nlat         = (*MetaData).nlat;
  (*state).nlev         = (*MetaData).nlev;
  (*state).nss          = (*MetaData).nss;
  (*state).L            = (*MetaData).L;
  (*state).horiz_repres = horiz_repres;
  (*state).vert_repres  = vert_repres;
  (*state).temp_repres  = temp_repres;

  Array_1d_double_create (&((*state).times),
                          (*state).nss);
  for (t=0; t<(*state).nss; t++)
  { (*state).times[t] = (*MetaData).times[t];
  }

  Array_1d_double_create (&((*state).longitude),
                          (*state).nlon+2);
  for (x=0; x<(*state).nlon+2; x++)
  { (*state).longitude[x] = (*MetaData).longitude[x];
  }
  Array_1d_double_create (&((*state).latitude),
                          (*state).nlat+2);
  for (y=0; y<(*state).nlat+2; y++)
  { (*state).latitude[y] = (*MetaData).latitude[y];
    //printf ("Lat for point %u is %f\n", y, (*state).latitude[y]);
  }

  Array_1d_double_create (&((*state).level),
                          (*state).nlev+2);
  for (z=0; z<(*state).nlev+2; z++)
  { (*state).level[z] = (*MetaData).level[z];
  }

  // Allocation of main fields
  if (horiz_repres == 'r')
  { // The fields will be in real space
    Array_3d_double_create (&((*state).tracer0_rs),
                            (*state).nlon+2, (*state).nlat+2, (*state).nlev+2);
    Array_3d_double_create (&((*state).source_rs),
                            (*state).nlon+2, (*state).nlat+2, (*state).nss);
    (*state).tracer0_ss = NULL;
    (*state).source_ss  = NULL;
  }
  else
  { if (horiz_repres == 's')
    { // The fields will be in spectral space
      Array_4d_double_create (&((*state).tracer0_ss),
                              2, (*state).L+1, (*state).L+1, (*state).nlev+2);
      Array_4d_double_create (&((*state).source_ss),
                              2, (*state).L+1, (*state).L+1, (*state).nss);
    (*state).tracer0_rs = NULL;
    (*state).source_rs  = NULL;
    }
    else
    { printf ("ERROR - unknown field type given to Allocate_state\n");
      exit(0);
    }
  }
}


// -------------------------------------------------------------------------------
void Deallocate_state (struct state_type *state)
{ // Deallocate field

  Array_1d_double_destroy (&((*state).times));
  Array_1d_double_destroy (&((*state).longitude));
  Array_1d_double_destroy (&((*state).latitude));
  Array_1d_double_destroy (&((*state).level));
  if ((*state).tracer0_rs !=0)
  { Array_3d_double_destroy (&((*state).tracer0_rs),
                            (*state).nlon+2, (*state).nlat+2);
  }
  if ((*state).source_rs !=0)
  { Array_3d_double_destroy (&((*state).source_rs),
                            (*state).nlon+2, (*state).nlat+2);
  }
  if ((*state).tracer0_ss != 0)
  { Array_4d_double_destroy (&((*state).tracer0_ss),
                             2, (*state).L+1, (*state).L+1);
  }
  if ((*state).source_ss != 0)
  { Array_4d_double_destroy (&((*state).source_ss),
                             2, (*state).L+1, (*state).L+1);
  }
}


// -------------------------------------------------------------------------------
void Allocate_instant (struct instant_tracer_type *tracer,
                       struct metadata_type       *MetaData)
{ // Allocate instantaneous tracer field

  Array_3d_double_create (&((*tracer).tracer),
                          (*MetaData).nlon+2, (*MetaData).nlat+2, (*MetaData).nlev+2);

}


// -------------------------------------------------------------------------------
void Deallocate_instant (struct instant_tracer_type *tracer,
                         struct metadata_type       *MetaData)
{ // Allocate instantaneous tracer field

  Array_3d_double_destroy (&((*tracer).tracer),
                           (*MetaData).nlon+2, (*MetaData).nlat+2);

}


// -------------------------------------------------------------------------------
void Allocate_wind (struct Wind_type     *wind,
                    struct metadata_type *MetaData,
                    bool                 ecmwf)
{ // Allocate wind field
  int x, y, z;

  if (ecmwf)
  { // Allocate an ECMWF wind field
    (*wind).nlon = ::ecmwf_nlon;
    (*wind).nlat = ::ecmwf_nlat;
    (*wind).nlev = ::ecmwf_nlev;
  }
  else
  { (*wind).nlon = (*MetaData).nlon;
    (*wind).nlat = (*MetaData).nlat;
    (*wind).nlev = (*MetaData).nlev;
  }

  Array_1d_double_create (&((*wind).longitude_u),
                          (*wind).nlon+2);
  Array_1d_double_create (&((*wind).longitude_v),
                          (*wind).nlon+2);
  Array_1d_double_create (&((*wind).latitude_u),
                          (*wind).nlat+2);
  Array_1d_double_create (&((*wind).latitude_v),
                          (*wind).nlat+2);
  Array_1d_double_create (&((*wind).level_uv),
                          (*wind).nlev+2);
  Array_1d_double_create (&((*wind).level_w),
                          (*wind).nlev+2);

  Array_3d_double_create (&((*wind).u),
                          (*wind).nlon+2, (*wind).nlat+2, (*wind).nlev+2);
  Array_3d_double_create (&((*wind).v),
                          (*wind).nlon+2, (*wind).nlat+2, (*wind).nlev+2);
  Array_3d_double_create (&((*wind).w),
                          (*wind).nlon+2, (*wind).nlat+2, (*wind).nlev+2);

  // Fill-out the dimension values (if ECMWF, this is done on reading the ECMWF data)
  // Note that the MetaData structure includes halo points
  if (!ecmwf)
  { // Set the grid locations that are shared with the tracer and source
    for (x=0; x<(*wind).nlon+2; x++)
    { (*wind).longitude_v[x] = (*MetaData).longitude[x];
    }
    for (y=0; y<(*wind).nlat+2; y++)
    { (*wind).latitude_u[y] = (*MetaData).latitude[y];
    }
    for (z=0; z<(*wind).nlev+2; z++)
    { (*wind).level_uv[z] = (*MetaData).level[z];
    }

    // Set other grid locations
    for (x=0; x<(*wind).nlon+1; x++)
    { (*wind).longitude_u[x] = ((*wind).longitude_v[x] + (*wind).longitude_v[x+1]) / 2.0;
    }
    (*wind).longitude_u[(*wind).nlon+1] = (*wind).longitude_v[(*wind).nlon+1] +
                                          ((*wind).longitude_v[(*wind).nlon+1] -
                                           (*wind).longitude_v[(*wind).nlon]) / 2.0;
    for (y=1; y<(*wind).nlat+2; y++)
    { (*wind).latitude_v[y]  = ((*wind).latitude_u[y-1] + (*wind).latitude_u[y]) / 2.0;
    }
    (*wind).latitude_v[0] = (*wind).latitude_u[0] -
                            ((*wind).longitude_u[1] -
                             (*wind).longitude_u[0]) / 2.0;
    for (z=1; z<(*wind).nlev+2; z++)
    { (*wind).level_w[z] = ((*wind).level_uv[z-1] + (*wind).level_uv[z]) / 2.0;
    }
    (*wind).level_w[0] = (*wind).level_w[1] - ((*wind).level_w[1] - (*wind).level_w[0]) / 2.0;
  }
}

// -------------------------------------------------------------------------------
void Deallocate_wind (struct Wind_type *wind)
{ // Deallocate wind field

  Array_1d_double_destroy (&((*wind).longitude_u));
  Array_1d_double_destroy (&((*wind).longitude_v));
  Array_1d_double_destroy (&((*wind).latitude_u));
  Array_1d_double_destroy (&((*wind).latitude_v));
  Array_1d_double_destroy (&((*wind).level_uv));
  Array_1d_double_destroy (&((*wind).level_w));

  Array_3d_double_destroy (&((*wind).u),
                          (*wind).nlon+2, (*wind).nlat+2);
  Array_3d_double_destroy (&((*wind).v),
                          (*wind).nlon+2, (*wind).nlat+2);
  Array_3d_double_destroy (&((*wind).w),
                          (*wind).nlon+2, (*wind).nlat+2);
}


// -------------------------------------------------------------------------------
void Destroy_Obs ( struct obs_type **obs )
{ // Delete the structure created inside GenerateObs_arguments
  // Note that the obs structure is created in one of the following
  // GenerateObs_arguments
  // ReadObservations
  struct obs_type *current, *next;

  current = *obs;
  while (current)
  { next = (*current).next;
    delete current;
    current = next;
  }
}



// -------------------------------------------------------------------------------
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
                    bool   printdata)
{ int    x, y, z, t;
  // Copy meta data
  if (printdata)
  { printf ("Inside routine copy_metadata\n");
  }
  *nlon = nlon0;
  *nlat = nlat0;
  *nlev = nlev0;
  *nss  = nss0;
  *L    = L0;
  if (printdata)
  { printf ("Time dimension\n");
  }

  for (t=0; t<(*nss); t++)
  { times[t] = times0[t];
    if (printdata)
    { printf ("%f ", times[t]);
    }
  }

  for (x=0; x<(*nlon); x++)
  { longitude[x+1] = longitude0[x];
  }
  // Halo points
  longitude[0]       = longitude[1] - (longitude[2] - longitude[1]);
  longitude[nlon0+1] = longitude[nlon0] + longitude[nlon0] - longitude[nlon0-1];
  if (printdata)
  { printf ("\nLongitude dimension\n");
    for (x=0; x<(*nlon)+2; x++)
    { printf ("%f ", longitude[x]);
    }
  }

  for (y=0; y<(*nlat); y++)
  { latitude[y+1] = latitude0[y];
  }
  // Halo points
  latitude[0]       = latitude[1] - (latitude[2] - latitude[1]);
  latitude[nlat0+1] = latitude[nlat0] + latitude[nlat0] - latitude[nlat0-1];
  if (printdata)
  { printf ("\nLatitude dimension\n");
    for (y=0; y<(*nlat)+2; y++)
    { printf ("%f ", latitude[y]);
    }
  }

  for (z=0; z<(*nlev); z++)
  { level[z+1] = level0[z];
  }
  // Halo points
  level[0]       = level[1] - (level[2] - level[1]);
  level[nlev0+1] = level[nlev0] + level[nlev0] - level[nlev0-1];
  if (printdata)
  { printf ("\nLevel dimension\n");
    for (z=0; z<(*nlev)+2; z++)
    { printf ("%f ", level[z]);
    }
  }

  // Complete the cosine of latitudes
  Fill_cos_lats (latitude,
                 (*nlat),
                 cos_v_lat,
                 cos_u_lat,
                 printdata);

  if (printdata) { printf ("\n\n"); }
}

// -------------------------------------------------------------------------------
void Fill_cos_lats ( double *latitude,
                     int    nlat,
                     double *cos_v_lat,
                     double *cos_u_lat,
                     bool   printdata )
// From the latitudes, determine the cosine of the latitudes at v and u points
{ int    y;
  double lat;

  if (cos_u_lat)
  { if (printdata) { printf ("\nCompleting cos latitudes (u lats)\n"); }
    // u latitudes are the above latitudes
    for (y=0; y<nlat+2; y++)
    { cos_u_lat[y] = cos(::deg2rad * latitude[y]);
      if (printdata) { printf ("%f ", cos_u_lat[y]); }
    }
  }

  if (cos_v_lat)
  { if (printdata) { printf ("\nCompleting cos latitudes (v lats)\n"); }
    // u latitudes are the above latitudes
    for (y=1; y<nlat+2; y++)
    { // Find the in-between latitude
      lat = 0.5 * (latitude[y] + latitude[y-1]);
      cos_v_lat[y] = cos(::deg2rad * lat);
      if (printdata) { printf ("%f ", cos_v_lat[y]); }
    }
    cos_v_lat[0] = cos_v_lat[nlat+1];
  }
}


// -------------------------------------------------------------------------------
double Normal ( double mean,      //in
                double stddev,    //in
                int    *rndseq )  //inout

{ // Pull random number from Normal distribution
  // Declare local variables
  int    notyet;
  double x, y, G;

  notyet = 1;
  while (notyet)
  { // Choose a random number between -5*stddev and 5*stddev
    // Assume 'random' function gives 0 to 1
    x = 5.0 * stddev * (2.0 * randomno(rndseq) - 1.0);
    // Calculate the Gaussian probability of this random number
    G = exp(-x*x/(2.0 * stddev * stddev));
    // Calculate another random number 0 to 1
    y = randomno(rndseq);
    if (y <= G) { notyet = 0; }
  }
  x += mean;
  return x;
}


// -------------------------------------------------------------------------------
double randomno (int *rndseq)
{ // Generate a random number between 0 and 1
  // Declare local variables
  int    k;
  int    a    = 16807;
  int    m    = 2147483647;
  int    q    = 127773;
  int    r    = 2836;
  int    mask = 123459876;
  double am, result;

  am      = 1.0 / m;
  *rndseq = (*rndseq) ^ mask;
  k       = *rndseq / q;
  *rndseq = a * (*rndseq - k*q) - r*k;
  if (*rndseq < 0)
  { *rndseq += m;
  }
  result  = am * *rndseq;
  *rndseq = (*rndseq) ^ mask;
  return result;
}


// -------------------------------------------------------------------------------
double innerproduct_realsp
  ( double                     **realspace1,     // in   input state 1
    double                     **realspace2,     // in   input state 2
    struct HorizTransData_type *HorizData,       // in   contains metadata
    bool                       inc_halos )       // in   switch to include halos
{ // Computes inner product between two states in real space
  double sum;
  int    lon, lat;
  int    nlon, nlat;
  int    offsetx, offsety;
  double *latitudinal;

  if (inc_halos)
  { nlon    = (*HorizData).nlon + 2;
    nlat    = (*HorizData).nlat + 2;
    offsetx = 0;
    offsety = 0;
  }
  else
  { nlon    = (*HorizData).nlon;
    nlat    = (*HorizData).nlat;
    offsetx = 1;
    offsety = 1;
  }

/*
  // This method does a proper integral
  Array_1d_double_create (&latitudinal,
                          nlat);

  sum = 0.0;
  for (lon=0; lon<nlon; lon++)
  { // Compute product for a particular longitude, which is to be passed to the integration routine
    for (lat=0; lat<nlat; lat++)
    { latitudinal[lat] = realspace1[lon+offsetx][lat+offsety] * realspace2[lon+offsetx][lat+offsety];
    }
    sum += LatIntegrate (latitudinal,
                         HorizData,
                         inc_halos);
  }

  Array_1d_double_destroy (&latitudinal);
*/


  // This method is a straighfoward sum
  sum = 0.0;
  for (lon=offsetx; lon<(nlon+offsetx); lon++)
  { for (lat=offsety; lat<(nlat+offsety); lat++)
    { sum += realspace1[lon][lat] * realspace2[lon][lat];
    }
  }

  return sum;      //* 2.0 * ::pi / double((*HorizData).nlon);
}



// -------------------------------------------------------------------------------
double LatIntegrate
  ( double                     *latitudinal,
    struct HorizTransData_type *HorizData,
    bool                       inc_halos )
{ // Does an integral in the latitude direction
  double sum;
  int    lat;

  if (inc_halos)
  { sum = latitudinal[0] * (*HorizData).GaussianWts[1];
    for (lat=0; lat<(*HorizData).nlat; lat++)
    { sum += latitudinal[lat+1] * (*HorizData).GaussianWts[lat];
    }
    sum += latitudinal[(*HorizData).nlat+1] * (*HorizData).GaussianWts[(*HorizData).nlat-2];
  }
  else
  { sum = 0.0;
    for (lat=0; lat<(*HorizData).nlat; lat++)
    { sum += latitudinal[lat] * (*HorizData).GaussianWts[lat];
    }
  }

  return sum;
}


// -------------------------------------------------------------------------------
double innerproduct_spherical
  ( double                     ***specspace1,    // in   input state 1
    double                     ***specspace2,    // in   input state 2
    struct HorizTransData_type *HorizData )      // in   contains metadata
{ // Computes inner product between two states in spherical space
  double sum;
  int    cs, l, m;

  // Do sum over the part that is represented directly
  sum = 0.0;
  for (l=0; l<=(*HorizData).L; l++)
  { sum += specspace1[0][l][0] * specspace2[0][l][0]; // cosine for m=0
    if (l > 0)
    { for (m=1; m<=l; m++)
      { for (cs=0; cs<2; cs++)
        { sum += specspace1[cs][l][m] * specspace2[cs][l][m];
        }
      }
    }
  }
  return sum;
}


// -------------------------------------------------------------------------------
double innerproduct_general
  ( struct state_type          *state1,      // in   input state 1
    struct state_type          *state2,      // in   input state 2
    struct HorizTransData_type *HorizData,   // in   contains metadata
    bool                       inc_halos )   // in   include halos
{ // Computes inner product between two states in the space that they are configured in
  double sum;
  int    lon, lat, lev, t, l, m, cs;
  int    nlon, nlat, nlev, offsetx, offsety, offsetz;
  double *latitudinal;
  double **ms_copy1, **ms_copy2, ***ss_copy1, ***ss_copy2;


  if (inc_halos)
  { nlon    = (*state1).nlon + 2;
    nlat    = (*state1).nlat + 2;
    nlev    = (*state1).nlev + 2;
    offsetx = 0;
    offsety = 0;
    offsetz = 0;
  }
  else
  { nlon    = (*state1).nlon;
    nlat    = (*state1).nlat;
    nlev    = (*state1).nlev;
    offsetx = 1;
    offsety = 1;
    offsetz = 1;
  }

  sum = 0.0;

  if (((*state1).horiz_repres == 'r') && ((*state2).horiz_repres == 'r'))
  { // We are in real space in the horizontal
    //printf ("Inner product in model space\n");

/*
    // This way does latitude weighting
    Array_2d_double_create (&ms_copy1,
                            (*state1).nlon+2, (*state1).nlat+2);
    Array_2d_double_create (&ms_copy2,
                            (*state1).nlon+2, (*state1).nlat+2);

    // The contribution to the inner product from the initial tracer
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { ms_copy1[lon][lat] = (*state1).tracer0_rs[lon][lat][lev];
          ms_copy2[lon][lat] = (*state2).tracer0_rs[lon][lat][lev];
        }
      }
      sum += innerproduct_realsp (ms_copy1,
                                  ms_copy2,
                                  HorizData,
                                  inc_halos);
    }

    // The contribution to the inner product from the source
    for (t=0; t<(*state1).nss; t++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { ms_copy1[lon][lat] = (*state1).source_rs[lon][lat][t];
          ms_copy2[lon][lat] = (*state2).source_rs[lon][lat][t];
        }
      }
      sum += innerproduct_realsp (ms_copy1,
                                  ms_copy2,
                                  HorizData,
                                  inc_halos);
    }

    Array_2d_double_destroy (&ms_copy1,
                             nlon);
    Array_2d_double_destroy (&ms_copy2,
                             nlon);

*/

    // This way does not do latitude weighting
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { sum += (*state1).tracer0_rs[lon][lat][lev] * 
                 (*state2).tracer0_rs[lon][lat][lev];
        }
      }
    }
    for (t=0; t<(*state1).nss; t++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { sum += (*state1).source_rs[lon][lat][t] *
                 (*state2).source_rs[lon][lat][t];
        }
      }
    }
  }


  if (((*state1).horiz_repres == 's') && ((*state2).horiz_repres == 's'))
  { // We are in spectral space in the horizontal
    //printf ("Inner product in spectral space\n");
    Array_3d_double_create (&ss_copy1,
                            2, (*HorizData).L+1, (*HorizData).L+1);
    Array_3d_double_create (&ss_copy2,
                            2, (*HorizData).L+1, (*HorizData).L+1);

    // The contribution to the inner product from the initial tracer
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (l=0; l<=(*HorizData).L; l++)
      { ss_copy1[0][l][0]= (*state1).tracer0_ss[0][l][0][lev]; // cosine for m=0
        ss_copy2[0][l][0]= (*state2).tracer0_ss[0][l][0][lev]; // cosine for m=0
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { ss_copy1[cs][l][m]= (*state1).tracer0_ss[cs][l][m][lev];
              ss_copy2[cs][l][m]= (*state2).tracer0_ss[cs][l][m][lev];
            }
          }
        }
      }
      sum += innerproduct_spherical (ss_copy1,
                                     ss_copy2,
                                     HorizData);
    }

    // The contribution to the inner product from the source
    for (t=0; t<(*state1).nss; t++)
    { for (l=0; l<=(*HorizData).L; l++)
      { ss_copy1[0][l][0] = (*state1).source_ss[0][l][0][t]; // cosine for m=0
        ss_copy2[0][l][0] = (*state2).source_ss[0][l][0][t]; // cosine for m=0
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { ss_copy1[cs][l][m] = (*state1).source_ss[cs][l][m][t];
              ss_copy2[cs][l][m] = (*state2).source_ss[cs][l][m][t];
            }
          }
        }
      }
      sum += innerproduct_spherical (ss_copy1,
                                     ss_copy2,
                                     HorizData);
    }


    Array_3d_double_destroy (&ss_copy1,
                             2, (*HorizData).L+1);
    Array_3d_double_destroy (&ss_copy2,
                             2, (*HorizData).L+1);
  }

  return sum;
}


// -------------------------------------------------------------------------------
void add_general
  ( struct state_type          *state1,      // in   input state 1
    struct state_type          *state2,      // in   input state 2
    struct state_type          *sum,         // out  sum of above
    struct HorizTransData_type *HorizData,   // in   contains metadata
    bool                       inc_halos )   // in   include halos
{ // Adds two states (whatever space they are configured)
  int    lon, lat, lev, t, l, m, cs;
  int    nlon, nlat, nlev;
  int    offsetx, offsety, offsetz;

  if (inc_halos)
  { nlon    = (*state1).nlon + 2;
    nlat    = (*state1).nlat + 2;
    nlev    = (*state1).nlev + 2;
    offsetx = 0;
    offsety = 0;
    offsetz = 0;
  }
  else
  { nlon    = (*state1).nlon;
    nlat    = (*state1).nlat;
    nlev    = (*state1).nlev;
    offsetx = 1;
    offsety = 1;
    offsetz = 1;
  }

  if (((*state1).horiz_repres == 'r') && ((*state2).horiz_repres == 'r'))
  { // We are in real space in the horizontal
    //printf ("Sum in model space\n");

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*sum).tracer0_rs[lon][lat][lev] = (*state1).tracer0_rs[lon][lat][lev] +
                                             (*state2).tracer0_rs[lon][lat][lev];
        }
      }
    }
    // Deal with the flux part
    for (t=0; t<(*state1).nss; t++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*sum).source_rs[lon][lat][t] = (*state1).source_rs[lon][lat][t] +
                                          (*state2).source_rs[lon][lat][t];
        }
      }
    }
  }


  if (((*state1).horiz_repres == 's') && ((*state2).horiz_repres == 's'))
  { // We are in spectral space in the horizontal
    //printf ("Sum in spectral space\n");

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*sum).tracer0_ss[0][l][0][lev] = (*state1).tracer0_ss[0][l][0][lev] +
                                          (*state2).tracer0_ss[0][l][0][lev];
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*sum).tracer0_ss[cs][l][m][lev] = (*state1).tracer0_ss[cs][l][m][lev] +
                                                 (*state2).tracer0_ss[cs][l][m][lev];
            }
          }
        }
      }
    }

    // Deal with the flux part
    for (t=0; t<(*state1).nss; t++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*sum).source_ss[0][l][0][t] = (*state1).source_ss[0][l][0][t] +
                                       (*state2).source_ss[0][l][0][t];
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*sum).source_ss[cs][l][m][t] = (*state1).source_ss[cs][l][m][t] +
                                              (*state2).source_ss[cs][l][m][t];
            }
          }
        }
      }
    }
  }
}



// -------------------------------------------------------------------------------
void add_general_fac
  ( struct state_type          *state1,      // in   input state 1
    struct state_type          *state2,      // in   input state 2
    struct state_type          *sum,         // out  state1 + factor * state2
    struct HorizTransData_type *HorizData,   // in   contains metadata
    bool                       inc_halos,    // in   include halos
    double                     factor )
{ // Adds state1 with factor * state2 (whatever space they are configured)
  int    lon, lat, lev, t, l, m, cs;
  int    nlon, nlat, nlev;
  int    offsetx, offsety, offsetz;

  if (inc_halos)
  { nlon    = (*state1).nlon + 2;
    nlat    = (*state1).nlat + 2;
    nlev    = (*state1).nlev + 2;
    offsetx = 0;
    offsety = 0;
    offsetz = 0;
  }
  else
  { nlon    = (*state1).nlon;
    nlat    = (*state1).nlat;
    nlev    = (*state1).nlev;
    offsetx = 1;
    offsety = 1;
    offsetz = 1;
  }

  if (((*state1).horiz_repres == 'r') && ((*state2).horiz_repres == 'r'))
  { // We are in real space in the horizontal
    //printf ("Sum in model space\n");

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*sum).tracer0_rs[lon][lat][lev] = (*state1).tracer0_rs[lon][lat][lev] +
                                             factor * (*state2).tracer0_rs[lon][lat][lev];
        }
      }
    }
    // Deal with the flux part
    for (t=0; t<(*state1).nss; t++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*sum).source_rs[lon][lat][t] = (*state1).source_rs[lon][lat][t] +
                                          factor * (*state2).source_rs[lon][lat][t];
        }
      }
    }
  }


  if (((*state1).horiz_repres == 's') && ((*state2).horiz_repres == 's'))
  { // We are in spectral space in the horizontal
    //printf ("Sum in spectral space\n");

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*sum).tracer0_ss[0][l][0][lev] = (*state1).tracer0_ss[0][l][0][lev] +
                                          factor * (*state2).tracer0_ss[0][l][0][lev];
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*sum).tracer0_ss[cs][l][m][lev] = (*state1).tracer0_ss[cs][l][m][lev] +
                                                 factor * (*state2).tracer0_ss[cs][l][m][lev];
            }
          }
        }
      }
    }

    // Deal with the flux part
    for (t=0; t<(*state1).nss; t++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*sum).source_ss[0][l][0][t] = (*state1).source_ss[0][l][0][t] +
                                       factor * (*state2).source_ss[0][l][0][t];
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*sum).source_ss[cs][l][m][t] = (*state1).source_ss[cs][l][m][t] +
                                              factor * (*state2).source_ss[cs][l][m][t];
            }
          }
        }
      }
    }
  }
}



// -------------------------------------------------------------------------------
void copy_general
  ( struct state_type          *state,      // in   input state
    struct state_type          *copy,       // out  copy of above
    struct HorizTransData_type *HorizData,  // in   contains metadata
    bool                       inc_halos )  // in   include halos
{ // Takes a copy of state (whatever space they are configured)
  int    lon, lat, lev, t, l, m, cs;
  int    nlon, nlat, nlev;
  int    offsetx, offsety, offsetz;

  if (inc_halos)
  { nlon    = (*state).nlon + 2;
    nlat    = (*state).nlat + 2;
    nlev    = (*state).nlev + 2;
    offsetx = 0;
    offsety = 0;
    offsetz = 0;
  }
  else
  { nlon    = (*state).nlon;
    nlat    = (*state).nlat;
    nlev    = (*state).nlev;
    offsetx = 1;
    offsety = 1;
    offsetz = 1;
  }

  if ((*state).horiz_repres == 'r')
  { // We are in real space in the horizontal

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*copy).tracer0_rs[lon][lat][lev] = (*state).tracer0_rs[lon][lat][lev];
        }
      }
    }
    // Deal with the flux part
    for (t=0; t<(*state).nss; t++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*copy).source_rs[lon][lat][t] = (*state).source_rs[lon][lat][t];
        }
      }
    }
  }

  if ((*state).horiz_repres == 's')
  { // We are in spectral space in the horizontal

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*copy).tracer0_ss[0][l][0][lev] = (*state).tracer0_ss[0][l][0][lev];
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*copy).tracer0_ss[cs][l][m][lev] = (*state).tracer0_ss[cs][l][m][lev];
            }
          }
        }
      }
    }

    // Deal with the flux part
    for (t=0; t<(*state).nss; t++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*copy).source_ss[0][l][0][t] = (*state).source_ss[0][l][0][t];
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*copy).source_ss[cs][l][m][t] = (*state).source_ss[cs][l][m][t];
            }
          }
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------
void multiply_general
  ( struct state_type          *state,        // inout state
    struct HorizTransData_type *HorizData,    // in    contains metadata
    bool                       inc_halos,     // in    include halos
    double                     factor_tracer, // in    multiplication factor for tracer
    double                     factor_flux)   // in    multiplication factor for flux
{ // Multiplies state vector by factor
  int    lon, lat, lev, t, l, m, cs;
  int    nlon, nlat, nlev;
  int    offsetx, offsety, offsetz;

  if (inc_halos)
  { nlon    = (*state).nlon + 2;
    nlat    = (*state).nlat + 2;
    nlev    = (*state).nlev + 2;
    offsetx = 0;
    offsety = 0;
    offsetz = 0;
  }
  else
  { nlon    = (*state).nlon;
    nlat    = (*state).nlat;
    nlev    = (*state).nlev;
    offsetx = 1;
    offsety = 1;
    offsetz = 1;
  }

  if ((*state).horiz_repres == 'r')
  { // We are in real space in the horizontal

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*state).tracer0_rs[lon][lat][lev] *= factor_tracer;
        }
      }
    }
    // Deal with the flux part
    for (t=0; t<(*state).nss; t++)
    { for (lon=offsetx; lon<(nlon+offsetx); lon++)
      { for (lat=offsety; lat<(nlat+offsety); lat++)
        { (*state).source_rs[lon][lat][t] *= factor_flux;
        }
      }
    }
  }

  if ((*state).horiz_repres == 's')
  { // We are in spectral space in the horizontal

    // Deal with the tracer part
    for (lev=offsetz; lev<(nlev+offsetz); lev++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*state).tracer0_ss[0][l][0][lev] *= factor_tracer;
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*state).tracer0_ss[cs][l][m][lev] *= factor_tracer;
            }
          }
        }
      }
    }

    // Deal with the flux part
    for (t=0; t<(*state).nss; t++)
    { for (l=0; l<=(*HorizData).L; l++)
      { // cosine for m=0
        (*state).source_ss[0][l][0][t] *= factor_flux;
        if (l > 0)
        { for (m=1; m<=l; m++)
          { for (cs=0; cs<2; cs++)
            { (*state).source_ss[cs][l][m][t] *= factor_flux;
            }
          }
        }
      }
    }
  }
}


// -------------------------------------------------------------------------------
int Find_index_ascend ( int    Nels,       // Number of elements in Values array
                        double Values[],   // Ascending array
                        double Value )
{ // Find the lower index in an ascending array
  // Declare local variables
  int       result, loop;
  double    min, max, here, next;

  min    = Values[0];
  max    = Values[Nels-1];
  result = -1;
  next   = min;

  // If the value is the same as the last point in the Values array
  if (Value == max)
  { result = Nels-1;
  }
  else
  { if ((Value >= min) && (Value <= max) && (Nels > 1))
    { for (loop=0; ((result==-1) && (loop<(Nels-1))); loop++)
      { here = next;
        next = Values[loop+1];
        if (next > Value)
        { result = loop;
        }
      }
    }
  }
  return result;
}


// -------------------------------------------------------------------------------
int Find_index_descend ( int    Nels,       // Number of elements in Values array
                         double Values[],   // Descending array
                         double Value )
{ // Find the lower index in a decending array
  // Declare local variables
  int       result, loop;
  double    min, max, here, next;

  min    = Values[Nels-1];
  max    = Values[0];
  result = 0;
  next   = max;

  // If the value is the same as the last point in the Values array
  if (Value == min)
  { result = Nels-1;
  }
  else
  { if ((Value >= min) && (Value <= max) && (Nels > 1))
    { for (loop=0; ((result==0) && (loop<(Nels-1))); loop++)
      { here = next;
        next = Values[loop+1];
        //printf ("    here = %f, next = %f, value = %f\n", here, next, Value);
        if (next < Value)
        { result = loop;
        }
      }
    }
  }
  return result;
}


// -------------------------------------------------------------------------------
int Find_index_descend_mod ( int    Nels,       // Number of elements in Values array
                             double Values[],   // Descending array
                             double Value )
{ // Find the lower index in a descending array
  // This is different from Find_index_descend as it allows some points north of the northmost
  // This isn't suitable for the Semi-Lagrangian scheme, but is suitable for the obs op
  // Declare local variables
  int       result, loop;
  double    min, max, here, next;

  min    = Values[Nels-1];
  max    = Values[0];
  result = -1;
  next   = max;

  // If the value is the same as the last point in the Values array
  if (Value == min)
  { result = Nels-1;
  }
  else
  { if ((Value >= min) && (Value <= max) && (Nels > 1))
    { for (loop=0; ((result==-1) && (loop<(Nels-1))); loop++)
      { here = next;
        next = Values[loop+1];
        //printf ("Descend: here = %f, next = %f, value = %f\n", here, next, Value);
        if (next < Value)
        { result = loop;
        }
      }
    }
  }
  return result;
}




// -------------------------------------------------------------------------------
void construct_total_filename (char path[],
                               char name[],
                               char complete[])
{ // Construct the total filename path/name
  int c = 0, marker;

  while (path[c] != ' ' && path[c] != '\0' && path[c] != '\n')
  { complete[c] = path[c];
    c++;
  }
  complete[c] = '/';
  marker = c + 1;
  c = 0;
  while (name[c] != ' ' && name[c] != '\0' && name[c] != '\n')
  { complete[c+marker] = name[c];
    c++;
  }
  complete[c+marker] = '\0';
}


// -------------------------------------------------------------------------------
void halos_tracer (struct instant_tracer_type *tracer,
                   struct metadata_type       *MetaData)
{ // Sort halos in a tracer field
  int x, y, z, xp1, yp1, zp1, xpair;
  int nlon, nlat, nlev;

  nlon = (*MetaData).nlon;
  nlat = (*MetaData).nlat;
  nlev = (*MetaData).nlev;

  // Real space fields
  // (a) in the longitude direction
  for (y=0; y<nlat; y++)
  { yp1 = y + 1;
    for (z=0; z<nlev; z++)
    { zp1 = z + 1;
      (*tracer).tracer[0][yp1][zp1]      = (*tracer).tracer[nlon][yp1][zp1];
      (*tracer).tracer[nlon+1][yp1][zp1] = (*tracer).tracer[1][yp1][zp1];
    }
  }

  // (b) in the latitude direction
  for (x=0; x<nlon+2; x++)
  { xpair = x + nlon/2 + 1;
    if (xpair > nlon)
    { xpair -= nlon;
    }
    for (z=0; z<nlev; z++)
    { zp1 = z + 1;
      (*tracer).tracer[x][0][zp1]      = (*tracer).tracer[xpair][1][zp1];
      (*tracer).tracer[x][nlat+1][zp1] = (*tracer).tracer[xpair][nlat][zp1];
    }
  }

  // (c) in the height direction
  // Use Neuman boundary conditions for the tracer
  for (x=0; x<nlon+2; x++)
  { for (y=0; y<nlat+2; y++)
    { // Bottom of model
      (*tracer).tracer[x][y][0]      = (*tracer).tracer[x][y][1];
      // Top of model
      (*tracer).tracer[x][y][nlev+1] = (*tracer).tracer[x][y][nlev];
    }
  }
}


// -------------------------------------------------------------------------------
void halos_tracer_adj (struct instant_tracer_type *tracer,
                       struct metadata_type       *MetaData)
{ // Sort halos in a tracer field
  int x, y, z, xp1, yp1, zp1, xpair;
  int nlon, nlat, nlev;

  nlon = (*MetaData).nlon;
  nlat = (*MetaData).nlat;
  nlev = (*MetaData).nlev;

  // Real space fields
  // (c) in the height direction
  // Use Neuman boundary conditions for the tracer
  for (x=0; x<nlon+2; x++)
  { for (y=0; y<nlat+2; y++)
    { // Bottom of model
      (*tracer).tracer[x][y][1]    += (*tracer).tracer[x][y][0];
      // Top of model
      (*tracer).tracer[x][y][nlev] += (*tracer).tracer[x][y][nlev+1];
    }
  }

  // (b) in the latitude direction
  for (x=0; x<nlon+2; x++)
  { xpair = x + nlon/2 + 1;
    if (xpair > nlon)
    { xpair -= nlon;
    }
    for (z=0; z<nlev; z++)
    { zp1 = z + 1;
      (*tracer).tracer[xpair][1][zp1]    += (*tracer).tracer[x][0][zp1];
      (*tracer).tracer[xpair][nlat][zp1] += (*tracer).tracer[x][nlat+1][zp1];
    }
  }

  // (a) in the longitude direction
  for (y=0; y<nlat; y++)
  { yp1 = y + 1;
    for (z=0; z<nlev; z++)
    { zp1 = z + 1;
      (*tracer).tracer[nlon][yp1][zp1] += (*tracer).tracer[0][yp1][zp1];
      (*tracer).tracer[1][yp1][zp1]    += (*tracer).tracer[nlon+1][yp1][zp1];
    }
  }
}


// -------------------------------------------------------------------------------
void halos (struct state_type *field )
{ // Sort out the halos in a given field
  int x, y, z, xp1, yp1, zp1, t, xpair, l, m, cs;
  int nlon, nlat, nlev, L, nss;

  nlon = (*field).nlon;
  nlat = (*field).nlat;
  nlev = (*field).nlev;
  L    = (*field).L;
  nss  = (*field).nss;

  if ((*field).horiz_repres == 'r')
  { // Real space fields
    // (a) in the longitude direction
    for (y=0; y<nlat; y++)
    { yp1 = y + 1;
      for (z=0; z<nlev; z++)
      { zp1 = z + 1;
        (*field).tracer0_rs[0][yp1][zp1]      = (*field).tracer0_rs[nlon][yp1][zp1];
        (*field).tracer0_rs[nlon+1][yp1][zp1] = (*field).tracer0_rs[1][yp1][zp1];
      }
      for (t=0; t<nss; t++)
      { (*field).source_rs[0][yp1][t]      = (*field).source_rs[nlon][yp1][t];
        (*field).source_rs[nlon+1][yp1][t] = (*field).source_rs[1][yp1][t];
      }
    }

    // (b) in the latitude direction
    for (x=0; x<nlon+2; x++)
    { xpair = x + nlon/2 + 1;
      if (xpair > nlon)
      { xpair -= nlon;
      }
      for (z=0; z<nlev; z++)
      { zp1 = z + 1;
        (*field).tracer0_rs[x][0][zp1]      = (*field).tracer0_rs[xpair][2][zp1];
        (*field).tracer0_rs[x][nlat+1][zp1] = (*field).tracer0_rs[xpair][nlat-1][zp1];
      }
      for (t=0; t<nss; t++)
      { (*field).source_rs[x][0][t]      = (*field).source_rs[xpair][2][t];
        (*field).source_rs[x][nlat+1][t] = (*field).source_rs[xpair][nlat-1][t];
      }
    }

    // (c) in the height direction
    // Use Neuman boundary conditions for the tracer
    for (x=0; x<nlon+2; x++)
    { for (y=0; y<nlat+2; y++)
      { // Bottom of model
        (*field).tracer0_rs[x][y][0] = (*field).tracer0_rs[x][y][1];
        // Top of model
        (*field).tracer0_rs[x][y][nlev+1] = (*field).tracer0_rs[x][y][nlev];
      }
    }
  }
  else
  { // Spectral space fields
    if ((*field).vert_repres == 'r')
    { // Only halos that have to be dealt with are the vertical boundaries in the tracer0 field
      for (l=0; l<=L; l++)
      { // m=0 cos term only
        (*field).tracer0_ss[0][l][0][0]      = (*field).tracer0_ss[0][l][0][1];
        (*field).tracer0_ss[0][l][0][nlev+1] = (*field).tracer0_ss[0][l][0][nlev];
        if (l > 0)
        { // Loop over remaining zonal wavenumbers
          for (m=1; m<=l; m++)
          { // cos and sin terms
            for (cs=0; cs<2; cs++)
            { (*field).tracer0_ss[cs][l][m][0]      = (*field).tracer0_ss[cs][l][m][1];
              (*field).tracer0_ss[cs][l][m][nlev+1] = (*field).tracer0_ss[cs][l][m][nlev];
            }
          }
        }
      }
    }
  }
}



// -------------------------------------------------------------------------------
void halos_adj (struct state_type *field )
{ // Sort out the halos in a given field (adjoint operator)
  int x, y, z, xp1, yp1, zp1, t, xpair, l, m, cs;
  int nlon, nlat, nlev, L, nss;

  nlon = (*field).nlon;
  nlat = (*field).nlat;
  nlev = (*field).nlev;
  L    = (*field).L;
  nss  = (*field).nss;

  if ((*field).horiz_repres == 'r')
  { // Real space fields
    // (c) in the height direction
    // Use Neuman boundary conditions for the tracer
    for (x=0; x<nlon+2; x++)
    { for (y=0; y<nlat+2; y++)
      { // Bottom of model
        (*field).tracer0_rs[x][y][1]    += (*field).tracer0_rs[x][y][0];
        // Top of model
        (*field).tracer0_rs[x][y][nlev] += (*field).tracer0_rs[x][y][nlev+1];
      }
    }

    // (b) in the latitude direction
    for (x=0; x<nlon+2; x++)
    { xpair = x + nlon/2 + 1;
      if (xpair > nlon)
      { xpair -= nlon;
      }
      for (z=0; z<nlev; z++)
      { zp1 = z + 1;
        (*field).tracer0_rs[xpair][2][zp1]      += (*field).tracer0_rs[x][0][zp1];
        (*field).tracer0_rs[xpair][nlat-1][zp1] += (*field).tracer0_rs[x][nlat+1][zp1];
      }
      for (t=0; t<nss; t++)
      { (*field).source_rs[xpair][2][t]      += (*field).source_rs[x][0][t];
        (*field).source_rs[xpair][nlat-1][t] += (*field).source_rs[x][nlat+1][t];
      }
    }

    // (a) in the longitude direction
    for (y=0; y<nlat; y++)
    { yp1 = y + 1;
      for (z=0; z<nlev; z++)
      { zp1 = z + 1;
        (*field).tracer0_rs[nlon][yp1][zp1] += (*field).tracer0_rs[0][yp1][zp1];
        (*field).tracer0_rs[1][yp1][zp1]    += (*field).tracer0_rs[nlon+1][yp1][zp1];
      }
      for (t=0; t<nss; t++)
      { (*field).source_rs[nlon][yp1][t] += (*field).source_rs[0][yp1][t];
        (*field).source_rs[1][yp1][t]    += (*field).source_rs[nlon+1][yp1][t];
      }
    }

  }
  else
  { // Spectral space fields
    if ((*field).vert_repres == 'r')
    { // Only halos that have to be dealt with are the vertical boundaries in the tracer0 field
      for (l=0; l<=L; l++)
      { // m=0 cos term only
        (*field).tracer0_ss[0][l][0][1]    += (*field).tracer0_ss[0][l][0][0];
        (*field).tracer0_ss[0][l][0][nlev] += (*field).tracer0_ss[0][l][0][nlev+1];
        if (l > 0)
        { // Loop over remaining zonal wavenumbers
          for (m=1; m<=l; m++)
          { // cos and sin terms
            for (cs=0; cs<2; cs++)
            { (*field).tracer0_ss[cs][l][m][1]    += (*field).tracer0_ss[cs][l][m][0];
              (*field).tracer0_ss[cs][l][m][nlev] += (*field).tracer0_ss[cs][l][m][nlev+1];
            }
          }
        }
      }
    }
  }
}




// -------------------------------------------------------------------------------
void maxmin_tracer (struct metadata_type       *MetaData,
                    struct instant_tracer_type *field,
                    bool                       Inc_vert,
                    double                     *max_val,
                    double                     *min_val)
{ // Report on the maximum absolute value of the field
  int    x, y, z, nlev;
  double val, maxval, minval;

  if (Inc_vert)
  { nlev = (*MetaData).nlev;
  }
  else
  { nlev = 1;
  }

  maxval = (*field).tracer[1][1][1];
  minval = maxval;
  for (x=1; x<=(*MetaData).nlon; x++)
  { for (y=1; y<=(*MetaData).nlat; y++)
    { for (z=1; z<=nlev; z++)
      { val = fabs((*field).tracer[x][y][z]);
        if (val > maxval)
        { maxval = val;
        }
        if (val < minval)
        { minval = val;
        }
      }
    }
  }

  *max_val = maxval;
  *min_val = minval;
}



// -------------------------------------------------------------------------------
int TwoToOneD (int i,
               int j,
               int Nx,
               int Ny)
{ // 2D (i,j) to 1D index conversion.
  // Assume that 2D array indices i,j start at 1, but 1D array index starts at 0.
  // Look for boundary condition issues
  if (j < 1)
  { j = 1;
    i += Nx/2;
  }
  else
  { if (j > Ny)
    { j  = Ny;
      i += Nx/2;
    }
  }
  if (i < 1)
  { i += Nx;
  }
  else
  { if (i > Nx)
    { i -= Nx;
    }
  }

  return (i-1) * Ny + j - 1;
}



// -------------------------------------------------------------------------------
void GaussEl ( double **A,  // Matrix (modified on output)
               double *x,   // Solution (output)
               double *y,   // RHS (modified on output)
               int    k )
// Gaussian elimination for a full square matrix
// Solve A x = y
{ // Declare local variables
  int    i, j, jp;
  double recip_Aii, Aii, Aji, partialsum;
  bool   Show = false;

  if (Show)
  { // Show the matrix
    printf ("==== Gaussian Ellimination diagnostics (actual matrix) =====\n");
    for (i=0; i<k; i++)
    { printf ("%02i : ", i);
      for (j=0; j<k; j++)
      { printf ("%e ", A[i][j]);
      }
      printf (" : %e\n", y[i]);
    }
    printf ("==== ================================= =====\n");
  }

  // First go through the matrix
  for (i=0; i<(k-1); i++)
  { recip_Aii = 1.0 / A[i][i];
    // Modify later rows accordingly
    for (j=i+1; j<k; j++)
    { Aji = A[j][i];
      // Modify RHS
      y[j] -= Aji * y[i] * recip_Aii;
      // Modify matrix
      for (jp=i+1; jp<k; jp++)
      { A[j][jp] -= Aji * A[i][jp] * recip_Aii;
      }
    }
  }


  if (Show)
  { // Show the matrix (should be upper triangular)
    printf ("==== Gaussian Ellimination diagnostics (transformed matrix) =====\n");
    for (i=0; i<k; i++)
    { printf ("%02i : ", i);
      if (i > 0)
      { for (j=0; j<i; j++)
        { printf ("%e ", 0.0);
        }
      }
      for (j=i; j<k; j++)
      { printf ("%e ", A[i][j]);
      }
      printf (" : %e\n", y[i]);
    }
    printf ("==== ================================= =====\n");
  }

  // Now solve the matrix backwards
  for (i=k-1; i>=0; i--)
  { partialsum = 0.0;
    if (i<(k-1))
    { // Compute partial sum
      for (j=i+1; j<k; j++)
      { partialsum += A[i][j] * x[j];
      }
    }
    Aii       = A[i][i];
    recip_Aii = 1.0 / Aii;
    x[i]       = (y[i] - partialsum) * recip_Aii;
  }

}


// -------------------------------------------------------------------------------
void CalcErr (struct metadata_type *MetaData,        // Meta data
              struct state_type    *Diff,            // Difference state
              double               *dz,              // Vertical profile of layer thicknesses
              double               *density,         // Vertical profile of density
              double               *diff_tracer,     // Output: total difference for tracer (Tg)
              double               *rse_diff_tracer, // Output: root square difference for tracer (Tg)
              double               *diff_flux,       // Output: total difference for first flux field (Tg/month)
              double               *rse_diff_flux)   // Output: total root square difference for first flux field (Tg/month)

{ // Compute the total errors of a diff field, mass weighted
  double  dA, dV, dt, dlon_rad, dlat_rad, ten_9, ten_18, value;
  double  diff_tr, diff_fl, rse_tr, rse_fl;
  int     lon, lat, lev, t;
  int     lonp1, latp1, levp1;

  dlon_rad = ((*MetaData).longitude[1] - (*MetaData).longitude[0]) * ::deg2rad;
  ten_9    = 1.0E9;
  ten_18   = 1.0E18;
  dt       = (*MetaData).times[1] - (*MetaData).times[0];

  diff_tr = 0.0;
  rse_tr  = 0.0;
  diff_fl = 0.0;
  rse_fl  = 0.0;
  for (lon=0; lon<(*Diff).nlon; lon++)
  { lonp1 = lon + 1;
    for (lat=0; lat<(*Diff).nlat; lat++)
    { latp1    = lat + 1;
      dlat_rad = ((*MetaData).latitude[lat+2] - (*MetaData).latitude[lat]) / 2.0;
      dA       = ::Re * (*MetaData).cos_u_lat[latp1] * dlon_rad * ::Re * dlat_rad;
      for (lev=0; lev<(*Diff).nlev; lev++)
      { levp1    = lev + 1;
        dV       = dA * dz[lev];
        value    = (*Diff).tracer0_rs[lonp1][latp1][levp1] * dV * density[lev] / ten_9;
        diff_tr += value;
        rse_tr  += value * value;
      }
      //for (t=0; t<(*Diff).nss; t++)
      for (t=0; t<1; t++)
      { value    = (*Diff).source_rs[lonp1][latp1][t] * dA * dt / ten_18; // / CVTdata.source_stddev[lon][lat][t];
        diff_fl += value;
        rse_fl  += value * value;
      }
    }
  }

  *diff_tracer     = diff_tr / ten_9;
  *rse_diff_tracer = sqrt(rse_tr) / ten_9;
  *diff_flux       = diff_fl;
  *rse_diff_flux   = sqrt(rse_fl);
}
