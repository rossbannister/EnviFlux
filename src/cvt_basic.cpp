   #include <source.h>


// -------------------------------------------------------------------------------
void InitializeHorizTrans
  ( struct HorizTransData_type *HorizData,
    const char                 ALPfilename[256] )
{ // Subroutine to declare structure for horizontal data (e.g. Legendre polynomials)
  FILE*   ALPfile = NULL;
  char    blankline[] = "                                                            ";
  char    numberi[] = "      ";
  char    numberf[] = "                                ";
  int     L, ALPsize, lat, lon, l, m, index, nlon, ncmpx;
  double  inputd, dlon;
  int     inputi;

  // Start to read-in meta data from the file
  ALPfile = fopen (ALPfilename, "r");
  fgets  (blankline, 60, ALPfile);
  //printf ("1. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("2. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("3. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("4. %s", blankline);

  fgets  (blankline, 60, ALPfile);
  //printf ("5. %s", blankline);
  strncpy (numberi, blankline + 24, 6);
  //printf ("numberi = %s\n", numberi);
  L = atoi(numberi);
  printf ("L = %u\n", L);

  fgets  (blankline, 60, ALPfile);
  //printf ("6. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("7. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("8. %s", blankline);

  fgets  (blankline, 60, ALPfile);
  //printf ("9. %s", blankline);
  strncpy (numberi, blankline + 24, 6);
  //printf ("numberi = %s\n", numberi);
  ALPsize = atoi(numberi);
  printf ("ALPsize = %u\n", ALPsize);


  fgets  (blankline, 60, ALPfile);
  //printf ("%s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("%s", blankline);

  if (ALPsize != (L+1)*(L+2)/2)
  { printf ("Inconsistent settings in ALP file\n");
    printf ("%s\n", ALPfilename);
    printf ("L and ALPsize not consistent\n");
    exit(0);
  }

  // Set variables inside structure
  (*HorizData).L       = L;
  (*HorizData).nlon    = 2*L + 1;
  (*HorizData).nlat    = L + 1;
  (*HorizData).ALPsize = ALPsize;
  // Allocate arrays
  Array_2d_double_create ( &((*HorizData).assocLegPoly), L+1, ALPsize);
  Array_1d_double_create ( &((*HorizData).GaussianCosCoLats), L+2);
  Array_1d_double_create ( &((*HorizData).GaussianCoLats), L+2);
  Array_1d_double_create ( &((*HorizData).GaussianLats), L+2);
  Array_1d_double_create ( &((*HorizData).GaussianWts), L+1);
  Array_2d_int_create    ( &((*HorizData).Plm_index), L+1, L+1);
  Array_1d_double_create ( &((*HorizData).lons), 2*L+1);

  // Read in the co-latitudes and calculate the colatitudes and latitudes
  for (lat=0; lat<L+1; lat++)
  { fgets (numberf, 32, ALPfile);
    inputd = atof(numberf);
    (*HorizData).GaussianCosCoLats[lat] = inputd;
    //printf ("Gaussian cos colat %u = %22.16f\n", lat, inputd);
    (*HorizData).GaussianCoLats[lat]    = ::rad2deg * acos(inputd);
    (*HorizData).GaussianLats[lat]      = 90.0 - (*HorizData).GaussianCoLats[lat];
  }
  fgets  (blankline, 60, ALPfile);
  //printf ("10. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("11. %s", blankline);


  // Read in the Gaussian weights
  for (lat=0; lat<L+1; lat++)
  { fgets (numberf, 32, ALPfile);
    inputd = atof(numberf);
    (*HorizData).GaussianWts[lat] = inputd;
    //printf ("Gaussian weight %u = %22.16f\n", lat, inputd);
  }
  fgets  (blankline, 60, ALPfile);
  //printf ("12. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("13. %s", blankline);

  // Read in the plm indices (these say which column of the associated Legendre Polynomial
  // array corresponds to a given l and m
  for (m=0; m<L+1; m++)
  { for (l=0; l<L+1; l++)
    { fscanf (ALPfile, "%i", &inputi);
      //printf ("plm index %u %u = %i\n", l, m, inputi);
      (*HorizData).Plm_index[l][m] = inputi-1;
    }
  }
  fgets  (blankline, 60, ALPfile);
  //printf ("14. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("15. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  //printf ("16. %s", blankline);

  // Read in the associated Legendre polynomials
  for (m=0; m<=L; m++)
  { for (l=m; l<=L; l++)
    { fgets  (blankline, 60, ALPfile);
      index = (*HorizData).Plm_index[l][m];
      for (lat=0; lat<=L; lat++)
      { fgets (numberf, 32, ALPfile);
        inputd = atof(numberf);
        //printf ("alp %u %u lat %u (plmindex=%i) = %22.16f\n", l, m, lat, index, inputd);
        (*HorizData).assocLegPoly[lat][index] = inputd;
      }
    }
  }
  fgets  (blankline, 60, ALPfile);
  //printf ("17. %s", blankline);
  fgets  (blankline, 60, ALPfile);
  printf ("%s", blankline);

  fclose(ALPfile);

  // Set the longitudes
  dlon = 360.0 / double((*HorizData).nlon);
  for (lon=0; lon<(*HorizData).nlon; lon++)
  { (*HorizData).lons[lon] = double(lon) * dlon;
  }


  // Initialise some things concerning the fftw library
  nlon  = (*HorizData).nlon;
  ncmpx = nlon/2 + 1;
  // The real-space input
  (*HorizData).fft_input_real  = (double*) fftw_malloc(sizeof(double) * nlon);
  // The real-space output
  (*HorizData).fft_output_real = (double*) fftw_malloc(sizeof(double) * nlon);
  // The spectral space input
  (*HorizData).fft_input_spec  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ncmpx);
  // The spectral space output
  (*HorizData).fft_output_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ncmpx);
  // The 'plan' needed to go from real to spectral
  (*HorizData).plan_dft_r2c_1d = fftw_plan_dft_r2c_1d(nlon,
                                                      (*HorizData).fft_input_real,
                                                      (*HorizData).fft_output_spec,
                                                      FFTW_ESTIMATE);
  // The 'plan' needed to go from spectral to real
  (*HorizData).plan_dft_c2r_1d = fftw_plan_dft_c2r_1d(nlon,
                                                      (*HorizData).fft_input_spec,
                                                      (*HorizData).fft_output_real,
                                                      FFTW_ESTIMATE);

  // (last argument is either FFTW_ESTIMATE (quick to initialise, sub-optimal)
  //                       or FFTW_MEASURE  (slow to initialise, fast to run)
}


// -------------------------------------------------------------------------------
void DeallocateHorizTrans
  ( struct HorizTransData_type *HorizData )
{ // Subroutine to deallocate structure for horizontal data
  int L;

  L = (*HorizData).L;

  Array_2d_double_destroy ( &((*HorizData).assocLegPoly), L+1);
  Array_1d_double_destroy ( &((*HorizData).GaussianCosCoLats));
  Array_1d_double_destroy ( &((*HorizData).GaussianCoLats));
  Array_1d_double_destroy ( &((*HorizData).GaussianLats));
  Array_1d_double_destroy ( &((*HorizData).GaussianWts));
  Array_2d_int_destroy    ( &((*HorizData).Plm_index), L+1);
  Array_1d_double_destroy ( &((*HorizData).lons));

  fftw_free((*HorizData).fft_input_real);
  fftw_free((*HorizData).fft_output_real);
  fftw_free((*HorizData).fft_input_spec);
  fftw_free((*HorizData).fft_output_spec);
  fftw_destroy_plan((*HorizData).plan_dft_r2c_1d);
  fftw_destroy_plan((*HorizData).plan_dft_c2r_1d);
}


// -------------------------------------------------------------------------------
void Allocate_CVT (struct CVTData_type *CVTdata,
                   int                 L,
                   int                 nlev,
                   int                 nss)
{ // Allocate the CVT arrays
  (*CVTdata).L    = L;
  (*CVTdata).nlon = 2*L + 1;
  (*CVTdata).nlat = L + 1;
  (*CVTdata).nlev = nlev;
  (*CVTdata).nss  = nss;

  Array_1d_double_create (&((*CVTdata).times),
                          nss);
  Array_1d_double_create (&((*CVTdata).longitude),
                          (*CVTdata).nlon);
  Array_1d_double_create (&((*CVTdata).latitude),
                          (*CVTdata).nlat);
  Array_1d_double_create (&((*CVTdata).level),
                          nlev);

  Array_1d_double_create (&((*CVTdata).tracer_stddev),
                          nlev);
  Array_2d_double_create (&((*CVTdata).tracer_hspec),
                          L+1,
                          nlev);
  Array_2d_double_create (&((*CVTdata).tracer_vertmodes),
                          nlev,
                          nlev);
  Array_1d_double_create (&((*CVTdata).tracer_vspec),
                          nlev);
  Array_3d_double_create (&((*CVTdata).source_stddev),
                          (*CVTdata).nlon,
                          (*CVTdata).nlat,
                          nss);
  Array_1d_double_create (&((*CVTdata).source_hspec),
                          L+1);
  Array_2d_double_create (&((*CVTdata).source_tempmodes),
                          nss,
                          nss);
  Array_1d_double_create (&((*CVTdata).source_tspec),
                          nss);
}



// -------------------------------------------------------------------------------
void Deallocate_CVT (struct CVTData_type *CVTdata)
{ // Deallocate the CVT arrays
  Array_1d_double_destroy (&((*CVTdata).times));
  Array_1d_double_destroy (&((*CVTdata).longitude));
  Array_1d_double_destroy (&((*CVTdata).latitude));
  Array_1d_double_destroy (&((*CVTdata).level));
  Array_1d_double_destroy (&((*CVTdata).tracer_stddev));
  Array_2d_double_destroy (&((*CVTdata).tracer_hspec),
                           (*CVTdata).L+1);
  Array_2d_double_destroy (&((*CVTdata).tracer_vertmodes),
                           (*CVTdata).nlev);
  Array_1d_double_destroy (&((*CVTdata).tracer_vspec));
  Array_3d_double_destroy (&((*CVTdata).source_stddev),
                           (*CVTdata).nlon,
                           (*CVTdata).nlat);
  Array_1d_double_destroy (&((*CVTdata).source_hspec));
  Array_2d_double_destroy (&((*CVTdata).source_tempmodes),
                           (*CVTdata).nss);
  Array_1d_double_destroy (&((*CVTdata).source_tspec));
}


// -------------------------------------------------------------------------------
void spherical_for
  ( double                     **realspace,       // out  real space (longs, lats)
    double                     ***spectral,       // in   spectral space (cs, l, m)
    struct HorizTransData_type *HorizData )       // in   Associated Legendre polynomials, etc
{ // Subroutine to perform the forward spectral transform
  // Assume that the structures have already been allocated
  // Ross Bannister, NCEO, September 2021
  // Based on the fortran SUBROUTINE Spherical_inv (built for INVICAT)

  int          L     = (*HorizData).L;
  int          nlon  = (*HorizData).nlon;
  int          nlat  = (*HorizData).nlat;
  double       *longitudinal;
  double       **spec1d_interim;
  double       ***interim = NULL;

  int          lon, lat, l, m, cs, index;
  double       sum;


  // Allocate the interim state (after Legendre, before Fourier)
  Array_3d_double_create (&interim,
                          2, L+1, L+1);

  // Allocate other intermediate states
  Array_1d_double_create (&longitudinal, nlon);
  Array_2d_double_create (&spec1d_interim, 2, L+1);

  // Do the Legendre transform
  for (cs=0; cs<2; cs++)
  { for (m=0; m<=L; m++)
    { if ((m > 0) || (cs == 0))        // no sin term for mm=0
      { for (lat=0; lat<=L; lat++)
        { sum   = 0.0;
          for (l=m; l<=L; l++)
          { index = (*HorizData).Plm_index[l][m];
            sum  += spectral[cs][l][m] * (*HorizData).assocLegPoly[lat][index];
          }
          interim[cs][m][lat] = sum;
        }
      }
    }
  }

  // Do the Fourier transform
  for (lat=0; lat<nlat; lat++)
  { // Copy interim line into temp array
    for (m=0; m<=L; m++)
    { for (cs=0; cs<2; cs++)
      { spec1d_interim[cs][m] = interim[cs][m][lat];
      }
    }
    fft_spec2real (longitudinal,     // out
                   spec1d_interim,   // in
                   HorizData);
    // Place the output in the appropriate place
    for (lon=0; lon<nlon; lon++)
    { realspace[lon+1][lat+1] = longitudinal[lon];   // The +1 is to allow for the halos
    }
  }

  // Tidy up
  Array_3d_double_destroy(&interim, 2, L+1);
  Array_1d_double_destroy(&longitudinal);
  Array_2d_double_destroy(&spec1d_interim, 2);
}




// -------------------------------------------------------------------------------
void spherical_adj
  ( double                     **realspace,       // in   real space (longs, lats)
    double                     ***spectral,       // out  spectral space (cs, l, m)
    struct HorizTransData_type *HorizData )       // in   Associated Legendre polynomials, etc
{ // Subroutine to perform the adjoint spectral transform
  // Assume that the structures have already been allocated
  // Ross Bannister, NCEO, September 2021
  // Based on the fortran SUBROUTINE Spherical_adj (built for INVICAT)

  int          L         = (*HorizData).L;
  int          nlon      = (*HorizData).nlon;
  int          nlat      = (*HorizData).nlat;
  double       dnlon     = double(nlon);
  double       dnlat     = double(nlat);
  double       half_nlon = dnlon / 2.0;
  double       *longitudinal;
  double       **spec1d_interim;
  double       ***interim = NULL;
  int          lon, lat, l, m, cs, index;
  double       sum;


  // Allocate the interim state (after Fourier, before Legendre)
  Array_3d_double_create (&interim,
                          2, L+1, L+1);

  // Allocate other intermediate states
  Array_1d_double_create (&longitudinal, nlon);
  Array_2d_double_create (&spec1d_interim, 2, L+1);


  // Do the Fourier transform
  for (lat=0; lat<nlat; lat++)
  { // Copy longitude line into temp array
    for (lon=0; lon<nlon; lon++)
    { longitudinal[lon] = realspace[lon+1][lat+1];   // The +1 is for the halos
    }

    // Do FFT to Fourier space for this latitude
    fft_real2spec (longitudinal,     // in
                   spec1d_interim,   // out
                   HorizData);

    // Place the output in the appropriate place
    interim[0][0][lat] = spec1d_interim[0][0] * dnlon;
    for (m=1; m<=L; m++)
    { for (cs=0; cs<2; cs++)
      { interim[cs][m][lat] = spec1d_interim[cs][m] * half_nlon;
      }
    }
  }


  // Do the Legendre transform
  for (cs=0; cs<2; cs++)
  { for (m=0; m<=L; m++)
    { if ((m > 0) || (cs == 0))        // no sin term for mm=0
      { for (l=m; l<=L; l++)
        { index = (*HorizData).Plm_index[l][m];
          sum   = 0.0;
          for (lat=0; lat<=L; lat++)
          { sum += interim[cs][m][lat] * (*HorizData).assocLegPoly[lat][index];
          }
          spectral[cs][l][m] = sum;
        }
      }
    }
  }

  // Tidy up
  Array_3d_double_destroy(&interim, 2, L+1);
  Array_1d_double_destroy(&longitudinal);
  Array_2d_double_destroy(&spec1d_interim, 2);
}





// -------------------------------------------------------------------------------
void spherical_inv
  ( double                     **realspace,       // in   real space (longs, lats)
    double                     ***spectral,       // out  spectral space (cs, l, m)
    struct HorizTransData_type *HorizData )       // in   Associated Legendre polynomials, etc
{ // Subroutine to perform the inverse spectral transform
  // Assume that the structures have already been allocated
  // Ross Bannister, NCEO, September 2021
  // Based on the fortran SUBROUTINE Spherical_inv (built for INVICAT)

  int          L     = (*HorizData).L;
  int          nlon  = (*HorizData).nlon;
  int          nlat  = (*HorizData).nlat;
  double       *longitudinal;
  double       **spec1d_interim;
  double       ***interim = NULL;
  int          lon, lat, l, m, cs, index;
  double       sum;


  // Allocate the interim state (after Fourier, before Legendre)
  Array_3d_double_create (&interim,
                          2, L+1, L+1);

  // Allocate other intermediate states
  Array_1d_double_create (&longitudinal, nlon);
  Array_2d_double_create (&spec1d_interim, 2, L+1);


  // Do the Fourier transform
  for (lat=0; lat<nlat; lat++)
  { // Copy longitude line into temp array
    for (lon=0; lon<nlon; lon++)
    { longitudinal[lon] = realspace[lon+1][lat+1];   // The +1 is for the halos
    }

    // Do FFT to Fourier space for this latitude
    fft_real2spec (longitudinal,
                   spec1d_interim,
                   HorizData);

    // Place the output in the appropriate place
    interim[0][0][lat] = spec1d_interim[0][0] * 0.5;
    for (m=1; m<=L; m++)
    { for (cs=0; cs<2; cs++)
      { interim[cs][m][lat] = spec1d_interim[cs][m] * 0.25;
      }
    }
  }


  // Do the Legendre transform
  for (cs=0; cs<2; cs++)
  { for (m=0; m<=L; m++)
    { if ((m > 0) || (cs == 0))        // no sin term for mm=0
      { for (l=m; l<=L; l++)
        { index = (*HorizData).Plm_index[l][m];
          sum   = 0.0;
          for (lat=0; lat<=L; lat++)
          { sum += interim[cs][m][lat] * (*HorizData).assocLegPoly[lat][index] *
                   (*HorizData).GaussianWts[lat];
          }
          spectral[cs][l][m] = sum;
        }
      }
    }
  }

  // Tidy up
  Array_3d_double_destroy(&interim, 2, L+1);
  Array_1d_double_destroy(&longitudinal);
  Array_2d_double_destroy(&spec1d_interim, 2);
}



// -------------------------------------------------------------------------------
void spherical_inv_adj
  ( double                     **realspace,       // out  real space (longs, lats)
    double                     ***spectral,       // in   spectral space (cs, l, m)
    struct HorizTransData_type *HorizData )       // in   Associated Legendre polynomials, etc
{ // Subroutine to perform the adjoint of the inverse spectral transform
  // Assume that the structures have already been allocated
  // Ross Bannister, NCEO, September 2021
  // Based on the fortran SUBROUTINE Spherical_inv (built for INVICAT)

  int          L     = (*HorizData).L;
  int          nlon  = (*HorizData).nlon;
  int          nlat  = (*HorizData).nlat;
  double       recip_nlon = 1.0 / double(nlon);
  double       *longitudinal;
  double       **spec1d_interim;
  double       ***interim = NULL;
  int          lon, lat, l, m, cs, index;
  double       sum, fac;


  // Allocate the interim state (after Legendre, before Fourier)
  Array_3d_double_create (&interim,
                          2, L+1, L+1);

  // Allocate other intermediate states
  Array_1d_double_create (&longitudinal, nlon);
  Array_2d_double_create (&spec1d_interim, 2, L+1);


  // Do the Legendre transform
  for (cs=0; cs<2; cs++)
  { for (m=0; m<=L; m++)
    { if ((m > 0) || (cs == 0))        // no sin term for mm=0
      { if (m == 0)
        { fac = 0.5;
        }
        else
        {fac = 0.25;
        }
        for (lat=0; lat<=L; lat++)
        { sum   = 0.0;
          for (l=m; l<=L; l++)
          { index = (*HorizData).Plm_index[l][m];
            sum  += spectral[cs][l][m] * (*HorizData).assocLegPoly[lat][index];
          }
          interim[cs][m][lat] = sum * fac;
        }
      }
    }
  }


  // Do the Fourier transform
  for (lat=0; lat<nlat; lat++)
  { // Copy interim line into temp array
    spec1d_interim[0][0] = interim[0][0][lat];
    for (m=1; m<=L; m++)
    { for (cs=0; cs<2; cs++)
      { spec1d_interim[cs][m] = 2.0 * interim[cs][m][lat];
      }
    }
    fft_spec2real (longitudinal,
                   spec1d_interim,
                   HorizData);
    // Place the output in the appropriate place
    for (lon=0; lon<nlon; lon++)
    { realspace[lon+1][lat+1] = longitudinal[lon];   // The +1 is for the halos
    }
  }

  // Tidy up
  Array_3d_double_destroy(&interim, 2, L+1);
  Array_1d_double_destroy(&longitudinal);
  Array_2d_double_destroy(&spec1d_interim, 2);
}






// -------------------------------------------------------------------------------
void fft_real2spec
  ( double                     *realspace,       // in   real space (nlon)
    double                     **specspace,      // out  in spectral space (2, L+1)
    struct HorizTransData_type *HorizData )      // in   contains data for fftw

{ // Interface routine to perform the real-to-spectral FFT
  // The interface calls fftw routines, but is meant to replicate the FFT parts of
  // the transform used in INVICAT (f90 subroutine Spherical_inv)
  // Assume that the structures have already been allocated
  // factors = true is equivalent to normal inverse cvt
  // factors = false is equivalent to adjoint of forward cvt

  // Ross Bannister, NCEO, September 2021

  int          nlon  = (*HorizData).nlon;
  int          L     = (*HorizData).L;
  int          lon, k, index;
  double       recip_nlon     = 1.0 / double(nlon);
  double       two_recip_nlon = 2.0 * recip_nlon;


  // Put the input in the relevant place
  for (lon=0; lon<nlon; lon++)
  { (*HorizData).fft_input_real[lon] = realspace[lon];
  }
  // Do FFT to Fourier space for this latitude
  fftw_execute ((*HorizData).plan_dft_r2c_1d);

  // Put the output in the correct place to mimic the INVICAT code
  // k = 0
  specspace[0][0] = (*HorizData).fft_output_spec[0][0] * recip_nlon;
  specspace[1][0] = 0.0;
  // k > 0
  for (k=1; k<=L; k++)
  { specspace[0][k] =        (*HorizData).fft_output_spec[k][0] * two_recip_nlon;  // Real part
    specspace[1][k] = -1.0 * (*HorizData).fft_output_spec[k][1] * two_recip_nlon;  // Imaginary part
  }
}



// -------------------------------------------------------------------------------
void fft_real2spec_adj
  ( double                     *realspace,       // out  real space (nlon)
    double                     **specspace,      // in   spectral space (2, L+1)
    struct HorizTransData_type *HorizData )      // in   contains data for fftw
{ // Interface routine to perform the spectral-to-real FFT
  // The interface calls fftw routines, but is meant to replicate the FFT parts of
  // the transform used in INVICAT (f90 subroutine Spherical_inv)
  // Assume that the structures have already been allocated
  // Ross Bannister, NCEO, September 2021

  int          nlon  = (*HorizData).nlon;
  int          L     = (*HorizData).L;
  int          lon, k, index;
  double       recip_nlon = 1.0 / double(nlon);

  // Prepare the input to the FFT routine in the correct place to mimic the INVICAT code
  // k = 0
  (*HorizData).fft_input_spec[0][0] = specspace[0][0] * recip_nlon;
  // k > 0
  for (k=1; k<=L; k++)
  { (*HorizData).fft_input_spec[k][0] = specspace[0][k] * recip_nlon;         // Real part
    (*HorizData).fft_input_spec[k][1] = -1.0 * specspace[1][k] * recip_nlon;  // Imaginary part
  }

  // Do the Fourier transform
  // Do FFT to Fourier space for this latitude
  fftw_execute ((*HorizData).plan_dft_c2r_1d);

  // Place result in correct place
  for (lon=0; lon<nlon; lon++)
  { realspace[lon] = (*HorizData).fft_output_real[lon];
  }
}


// -------------------------------------------------------------------------------
void fft_spec2real
  ( double                     *realspace,       // out  real space (nlon)
    double                     **specspace,      // in   spectral space (2, L+1)
    struct HorizTransData_type *HorizData )      // in   contains data for fftw
{ // Interface routine to perform the spectral-to-real FFT
  // The interface calls fftw routines, but is meant to replicate the FFT parts of
  // the transform used in INVICAT (f90 subroutine Spherical_inv)
  // Assume that the structures have already been allocated
  // Ross Bannister, NCEO, September 2021

  int          nlon  = (*HorizData).nlon;
  int          L     = (*HorizData).L;
  int          lon, k, index;


  // Prepare the input to the FFT routine in the correct place to mimic the INVICAT code
  // k = 0
  (*HorizData).fft_input_spec[0][0] = specspace[0][0];
  // k > 0
  for (k=1; k<=L; k++)
  { (*HorizData).fft_input_spec[k][0] = specspace[0][k] / 2.0;         // Real part
    (*HorizData).fft_input_spec[k][1] = -1.0 * specspace[1][k] / 2.0;  // Imaginary part
  }

  // Do the Fourier transform
  // Do FFT to Fourier space for this latitude
  fftw_execute ((*HorizData).plan_dft_c2r_1d);

  // Place result in correct place
  for (lon=0; lon<nlon; lon++)
  { realspace[lon] = (*HorizData).fft_output_real[lon];
  }
}


// -------------------------------------------------------------------------------
void fft_spec2real_adj
  ( double                     *realspace,       // in   real space (nlon)
    double                     **specspace,      // out  spectral space (2, L+1)
    struct HorizTransData_type *HorizData )      // in   contains data for fftw

{ // Interface routine to perform the real-to-spectral FFT
  // The interface calls fftw routines, but is meant to replicate the FFT parts of
  // the transform used in INVICAT (f90 subroutine Spherical_inv)
  // Assume that the structures have already been allocated
  // factors = true is equivalent to normal inverse cvt
  // factors = false is equivalent to adjoint of forward cvt

  // Ross Bannister, NCEO, September 2021

  int          nlon  = (*HorizData).nlon;
  int          L     = (*HorizData).L;
  int          lon, k, index;

  // Do the Fourier transform
  for (lon=0; lon<nlon; lon++)
  { (*HorizData).fft_input_real[lon] = realspace[lon];
  }
  // Do FFT to Fourier space for this latitude
  fftw_execute ((*HorizData).plan_dft_r2c_1d);

  // Put the output in the correct place to mimic the INVICAT code
  // This is the normal operation of the code
  // k = 0
  specspace[0][0] = (*HorizData).fft_output_spec[0][0];
  specspace[1][0] = 0.0;
  // k > 0
  for (k=1; k<=L; k++)
  { specspace[0][k] =        (*HorizData).fft_output_spec[k][0];  // Real part
    specspace[1][k] = -1.0 * (*HorizData).fft_output_spec[k][1];  // Imaginary part
  }
}
