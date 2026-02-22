/* ==================================================================================
   3d source and sink code
   Convert the flux field from an INVICAT file to this system's grid

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_Invicat2Enviflux.out

   Modification history
   --------------------
   14/03/24 Adapted from Master_MakeWinds. Ross Bannister
   09/04/25 Include tracer and flux.  Ross Bannister

   Documentation
   -------------
   Run with
   ./Master_Invicat2Enviflux.out <CVT file name (for meta data)> \
                                 <input filename containing INVICAT fluxes> \
                                 <output filename for fluxes on Envi-Flux grid> \
                                 <multiplication factor for flux> \
                                 <multiplication factor for tracer>


   =============================================================================== */

   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )
{ struct HorizTransData_type HorizData;
  struct CVTData_type        CVTdata;
  struct metadata_type       MetaData;
  int                        ok, x, y, z, t;
  char                       CVTfilename[256], INVICAT_file[256], ENVIFLUX_file[256];
  double                     maxflux_set, maxtracer_set, maxval, multiplier;
  struct state_type          enviflux, invicat;
  int                        low_index_lon, *low_index_lat;
  double                     xx[2], yy[2], vals[2][2];
  int                        use_nss;


  // ===============================================================================
  // Initialisation
  // ===============================================================================

  // Find the name of the file containing the full-resolution ECMWF winds
  ok = Invicat2Enviflux_arguments (argument_count,
                                   argument_list,
                                   CVTfilename,
                                   INVICAT_file,
                                   ENVIFLUX_file,
                                   &maxflux_set,
                                   &maxtracer_set);

  if (ok)
  { // Initialise and read in generic data needed for the Gaussian latitudes
    printf ("Initialising horizontal transform for Gaussian latitudes\n");
    InitializeHorizTrans (&HorizData,
                          ::ALPfilename);

    // Read-in and allocate cvt data
    cvt_matrices_input (&CVTdata,
                        CVTfilename);
    printf ("Read-in completed\n");

    // Read-in the INVICAT data (this will also allocate the structure)
    Read_INVICAT_all (&invicat,
                      INVICAT_file);

    // Decide on the number of months in the output (EnviFlux) file
    // Use whichever is the smaller
    if (invicat.nss > CVTdata.nss)
    { use_nss = CVTdata.nss;
    }
    else
    { use_nss = invicat.nss;
    }

    // Set-up the meta data structure
    Allocate_metadata (&MetaData,
                       CVTdata.nlon,
                       CVTdata.nlat,
                       CVTdata.nlev,
                       use_nss);

    copy_metadata (CVTdata.nlon,       &(MetaData.nlon),
                   CVTdata.nlat,       &(MetaData.nlat),
                   CVTdata.nlev,       &(MetaData.nlev),
                   use_nss,            &(MetaData.nss),
                   CVTdata.L,          &(MetaData.L),
                   CVTdata.times,      MetaData.times,
                   CVTdata.longitude,  MetaData.longitude,
                   CVTdata.latitude,   MetaData.latitude,
                   CVTdata.level,      MetaData.level,
                   MetaData.cos_u_lat, MetaData.cos_v_lat,
                   false);


    // Allocate space for the state vector (only the flux part will be used)
    Allocate_state (&enviflux,
                    &MetaData,
                    'r', 'r', 'r');



    // Check validity
    if (enviflux.nlev > invicat.nlev)
    { printf ("Fatal error, there are not enough levels in the INVICAT file\n");
      printf ("Number of levels in INVICAT file     : %i\n", invicat.nlev);
      printf ("Number expected from CVT file        : %i\n", enviflux.nlev);
      exit(0);
    }

    // ===============================================================================
    // Set up some reference arrays
    // ===============================================================================

    // Make an array for efficiency
    low_index_lat = new int [MetaData.nlat+2];

    // Check INVICAT arrays
    //printf ("INVICAT longitudes\n");
    //for (x=0; x<invicat.nlon+2; x++)
    //{ printf ("%f ", invicat.longitude[x]);
    //}
    //printf ("\n");
    //printf ("INVICAT latitudes\n");
    //for (y=0; y<invicat.nlat+2; y++)
    //{ printf ("%f ", invicat.latitude[y]);
    //}
    //printf ("\n");


    // ===============================================================================
    // Interpolate the fluxes
    // ===============================================================================

    printf ("Interpolating the fluxes\n");

    for (x=1; x<=enviflux.nlon; x++)
    { low_index_lon = Find_index_ascend (invicat.nlon+2,
                                         invicat.longitude,
                                         enviflux.longitude[x]);
      //printf ("Longitude %f is index %i between %f and %f\n",
      //        enviflux.longitude[x], low_index_lon, invicat.longitude[low_index_lon], invicat.longitude[low_index_lon+1]);
      if (low_index_lon <= -1)
      { printf ("%f longitude is out of range\n", enviflux.longitude[x]);
      }
      // Set the interpolation arrays for longitude
      xx[0] = invicat.longitude[low_index_lon];
      xx[1] = invicat.longitude[low_index_lon+1];

      for (y=1; y<=enviflux.nlat; y++)
      { if (x == 1)
        { // Only need to do this once
          low_index_lat[y] = Find_index_descend (invicat.nlat+2,
                                                 invicat.latitude,
                                                 enviflux.latitude[y]);
          //printf ("Latitude %f is index %i between %f and %f\n",
          //        enviflux.latitude[y], low_index_lat[y], invicat.latitude[low_index_lat[y]], invicat.latitude[low_index_lat[y]+1]);
          if (low_index_lat[y] <= -1)
          { printf ("%f latitude is out of range\n", enviflux.latitude[y]);
          }
        }

         // Set the interpolation arrays for latitude
        yy[0] = invicat.latitude[low_index_lat[y]];
        yy[0] = invicat.latitude[low_index_lat[y]+1];


        for (t=0; t<enviflux.nss; t++)
        { // Set the interpolation arrays
          vals[0][0] = invicat.source_rs[low_index_lon]  [low_index_lat[y]]  [t];
          vals[1][0] = invicat.source_rs[low_index_lon+1][low_index_lat[y]]  [t];
          vals[0][1] = invicat.source_rs[low_index_lon]  [low_index_lat[y]+1][t];
          vals[1][1] = invicat.source_rs[low_index_lon+1][low_index_lat[y]+1][t];

          // Interpolate
          enviflux.source_rs[x][y][t] = Interpolate2D (vals, xx, yy,
                                                       enviflux.longitude[x], enviflux.latitude[y]);
        }
      }
    }
    printf ("Interpolation of fluxes done\n");

    if (maxflux_set != 999.0)
    { // Scale the fluxes to have the required maximum value
      printf ("Scaling the fluxes\n");
      // Find the current max value
      maxval = 0.0;
      for (x=1; x<=enviflux.nlon; x++)
      { for (y=1; y<=enviflux.nlat; y++)
        { for (t=0; t<enviflux.nss; t++)
          { if (fabs(enviflux.source_rs[x][y][t]) > maxval)
            { maxval = fabs(enviflux.source_rs[x][y][t]);
            }
          }
        }
      }
      multiplier = maxflux_set / maxval;

      // Do the scaling
      for (x=1; x<=enviflux.nlon; x++)
      { for (y=1; y<=enviflux.nlat; y++)
        { for (t=0; t<enviflux.nss; t++)
          { enviflux.source_rs[x][y][t] *= multiplier;
          }
        }
      }


      // Check the maximum value of EnviFlux flux
      printf ("Checking max value of flux\n");
      // Find the current max value
      maxval = 0.0;
      for (x=1; x<=enviflux.nlon; x++)
      { for (y=1; y<=enviflux.nlat; y++)
        { for (t=0; t<enviflux.nss; t++)
          { if (fabs(enviflux.source_rs[x][y][t]) > maxval)
            { maxval = fabs(enviflux.source_rs[x][y][t]);
            }
          }
        }
      }
      printf ("New maximum value %f\n", maxval);

      printf ("Scaling the fluxes done\n");
    }


    // ===============================================================================
    // Interpolate the tracer
    // ===============================================================================

    printf ("Interpolating the tracer\n");

    for (x=1; x<=enviflux.nlon; x++)
    { low_index_lon = Find_index_ascend (invicat.nlon+2,
                                         invicat.longitude,
                                         enviflux.longitude[x]);
      //printf ("Longitude %f is index %i between %f and %f\n",
      //        enviflux.longitude[x], low_index_lon, invicat.longitude[low_index_lon], invicat.longitude[low_index_lon+1]);
      if (low_index_lon <= -1)
      { printf ("%f longitude is out of range\n", enviflux.longitude[x]);
      }
      // Set the interpolation arrays for longitude
      xx[0] = invicat.longitude[low_index_lon];
      xx[1] = invicat.longitude[low_index_lon+1];

      for (y=1; y<=enviflux.nlat; y++)
      { if (x == 1)
        { // Only need to do this once
          low_index_lat[y] = Find_index_descend (invicat.nlat+2,
                                                 invicat.latitude,
                                                 enviflux.latitude[y]);
          //printf ("Latitude %f is index %i between %f and %f\n",
          //        enviflux.latitude[y], low_index_lat[y], invicat.latitude[low_index_lat[y]], invicat.latitude[low_index_lat[y]+1]);
          if (low_index_lat[y] <= -1)
          { printf ("%f latitude is out of range\n", enviflux.latitude[y]);
          }
        }

         // Set the interpolation arrays for latitude
        yy[0] = invicat.latitude[low_index_lat[y]];
        yy[0] = invicat.latitude[low_index_lat[y]+1];


        for (z=0; z<enviflux.nlev; z++)
        { // Set the interpolation arrays
          vals[0][0] = invicat.tracer0_rs[low_index_lon]  [low_index_lat[y]]  [z];
          vals[1][0] = invicat.tracer0_rs[low_index_lon+1][low_index_lat[y]]  [z];
          vals[0][1] = invicat.tracer0_rs[low_index_lon]  [low_index_lat[y]+1][z];
          vals[1][1] = invicat.tracer0_rs[low_index_lon+1][low_index_lat[y]+1][z];

          // Interpolate
          enviflux.tracer0_rs[x][y][z] = Interpolate2D (vals, xx, yy,
                                                        enviflux.longitude[x], enviflux.latitude[y]);
        }
      }
    }
    printf ("Interpolation of tracer done\n");


    if (maxtracer_set != 999.0)
    { // Scale the tracer to have the required maximum value
      printf ("Scaling the tracer\n");
      // Find the current max value
      maxval = 0.0;
      for (x=1; x<=enviflux.nlon; x++)
      { for (y=1; y<=enviflux.nlat; y++)
        { for (z=0; z<enviflux.nlev; z++)
          { if (fabs(enviflux.tracer0_rs[x][y][z]) > maxval)
            { maxval = fabs(enviflux.tracer0_rs[x][y][z]);
            }
          }
        }
      }
      multiplier = maxtracer_set / maxval;

      // Do the scaling
      for (x=1; x<=enviflux.nlon; x++)
      { for (y=1; y<=enviflux.nlat; y++)
        { for (z=0; z<enviflux.nlev; z++)
          { enviflux.tracer0_rs[x][y][z] *= multiplier;
          }
        }
      }
      printf ("Scaling the tracer done\n");

      // Check the maximum value of EnviFlux tracer
      printf ("Checking max value of tracer\n");
      // Find the current max value
      maxval = 0.0;
      for (x=1; x<=enviflux.nlon; x++)
      { for (y=1; y<=enviflux.nlat; y++)
        { for (z=0; z<enviflux.nlev; z++)
          { if (fabs(enviflux.tracer0_rs[x][y][z]) > maxval)
            { maxval = fabs(enviflux.tracer0_rs[x][y][z]);
            }
          }
        }
      }
      printf ("New maximum value %f\n", maxval);


    }

    // Sort out halos
    halos (&invicat);

    // Output the result
    WriteStateVector (&enviflux,
                      ENVIFLUX_file);

    // Tidy up
    DeallocateHorizTrans (&HorizData);
    Deallocate_CVT (&CVTdata);
    Deallocate_metadata (&MetaData);
    Deallocate_state (&invicat);
    Deallocate_state (&enviflux);
    delete[] low_index_lat;
  }
}
