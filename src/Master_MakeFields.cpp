/* ==================================================================================
   3d source and sink code
   Generate initial conditions for source sink forward model

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_MakeFields.out

   Modification history
   --------------------
   26/09/22 New Code. Ross Bannister
   15/03/24 Allow conversion from INVICAT. Ross Bannister

   Documentation
   -------------
   Run with
   ./Master_MakeFields.out <CVT filename (for metadata)>
                           <Previous state file to increment (set to "nil" to set to zero)> \
                           <Number of blobs to add on to initial condition (0 for no blobs)> \
                           For each initial condition blob:
                             <long (deg)> <lat (deg)> <lev (m)> <amplitude (tracer units)> <horiz size (deg)> <vert size (m)> \
                           <Number of surface source fields (one per sDt time step)> \
                           <source timestep (source update, days, converted to secs)>
                           For each sDt time step:
                             <Number of blobs> \
                             For each surface source blob:
                               <long (deg)> <lat (deg)> <size (tracer units per s), +/-> <horiz size (grid points)>
                           Output file name (will contain initial conditions and surface flux fields)

   =============================================================================== */

   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )
{ char                 CVTfilename[256];
  char                 init_cond_file[256];
  char                 init_cond_file_type;
  int                  num_init_cond_blobs, ss, blob;
  struct blob_type     *ic_blob1 = NULL;
  int                  nss;
  int                  *num_ss_blobs = NULL;
  struct blob_type     **ss_field1 = NULL;
  char                 output_file[256];
  double               sDt;
  double               *SStimes = NULL;
  int                  ok;
  struct state_type    state;
  struct metadata_type MetaData;
  struct CVTData_type  CVTdata;
  int                  x, y, z, t, ghostx, ghosty, yghost_low, yghost_high;
  double               blobx, bloby, blobz;
  double               distx, disty, distz, disth_sq, distv_sq, size_h2, size_v2;
  struct blob_type     *current_blob = NULL;
  double               max_field, max_blob, factor;
  double               factor_ic, factor_flux, min_flux_set, orig_sign;
  int                  num_blobs_tracer, num_blobs_flux;

  // ===============================================================================
  // Initialisation
  // ===============================================================================

  ok = MakeFields_arguments ( argument_count,
                              argument_list,
                              CVTfilename,
                              init_cond_file,
                              &init_cond_file_type,
                              &factor_ic,
                              &factor_flux,
                              &num_init_cond_blobs,
                              &ic_blob1,
                              &nss,
                              &sDt,
                              &num_ss_blobs,
                              &ss_field1,
                              &min_flux_set,
                              output_file );

  if (ok)
  { // The following are done just to get the metadata (longs, lats, levs, etc)
    // Initialise and read in generic data needed for the horizontal transform

    // Allocate and read-in cvt data
    printf ("Allocating and reading CVT data\n");
    cvt_matrices_input (&CVTdata,
                        CVTfilename);
    printf ("Read-in completed\n");

    // Set-up the source/sink times
    SStimes = new double[nss];
    for (ss=0; ss<nss; ss++)
    { SStimes[ss] = double(ss) * sDt;
    }

    // Set-up the meta data structure
    printf ("Allocating meta data\n");
    Allocate_metadata (&MetaData,
                       CVTdata.nlon,
                       CVTdata.nlat,
                       CVTdata.nlev,
                       nss);
    printf ("Done\n");

    printf ("Copying meta data\n");
    copy_metadata (CVTdata.nlon,       &(MetaData.nlon),
                   CVTdata.nlat,       &(MetaData.nlat),
                   CVTdata.nlev,       &(MetaData.nlev),
                   nss,                &(MetaData.nss),
                   CVTdata.L,          &(MetaData.L),
                   SStimes,            MetaData.times,
                   CVTdata.longitude,  MetaData.longitude,
                   CVTdata.latitude,   MetaData.latitude,
                   CVTdata.level,      MetaData.level,
                   MetaData.cos_u_lat, MetaData.cos_v_lat,
                   false);
    printf ("Done\n");

    // Free-up space not needed any more
    printf ("Freeing up space\n");
    Deallocate_CVT (&CVTdata);
    delete[] SStimes;
    printf ("Done\n");


    // Read-in the initial conditions, including allocation
    if (strcmp ("nil", init_cond_file) != 0)
    { printf ("Reading initial state from file\n");
      state.horiz_repres = 'r';
      state.vert_repres  = 'r';
      state.temp_repres  = 'r';
      switch (init_cond_file_type)
      { case 'i': // INVICAT file
          ReadINVICAT (&state,
                       &MetaData,
                       true,       // Read tracer
                       true,       // Read source
                       init_cond_file);
          printf ("Done\n");
        break;
        case 'e': // Envi-flux file
          ReadStateVector (&state,
                           &MetaData,
                           true,
                           true,
                           init_cond_file);
          printf ("Done\n");
        break;
        default:
          printf ("Error - the initial state file type must be either (i)nvicat or (e)nviflux\n");
          exit(0);
      }

      // Multiply the state vector by the specified factors
      // Deal with the tracer part
      for (z=0; z<state.nlev+2; z++)
      { for (x=0; x<state.nlon+2; x++)
        { for (y=0; y<state.nlat+2; y++)
          { state.tracer0_rs[x][y][z] *= factor_ic;
          }
        }
      }
      // Deal with the flux part
      for (ss=0; ss<state.nss; ss++)
      { for (x=0; x<state.nlon+2; x++)
        { for (y=0; y<state.nlat+2; y++)
          { state.source_rs[x][y][ss] *= factor_flux;
          }
        }
      }

      // Check min fluxes and adjust
      for (ss=0; ss<state.nss; ss++)
      { for (x=0; x<state.nlon+2; x++)
        { for (y=0; y<state.nlat+2; y++)
          { if (fabs(state.source_rs[x][y][ss]) < min_flux_set)
            { orig_sign                 = copysign(1.0, state.source_rs[x][y][ss]);
              state.source_rs[x][y][ss] = orig_sign * min_flux_set;
            }
          }
        }
      }

    }
    else
    { // Allocate space for data
      printf ("Allocating state\n");
      Allocate_state (&state,
                      &MetaData,
                      'r', 'r', 'r');
      printf ("Done\n");
    }

    // Make blobs in the initial condition field
    current_blob     = ic_blob1;
    max_blob         = 0.0;
    num_blobs_tracer = 0;
    while (current_blob)
    { printf ("Adding a new blob to the intitial field\n");
      num_blobs_tracer++;
      // Monitor the size of the blobs
      if (fabs((*current_blob).amplitude) > max_blob)
      { max_blob = fabs((*current_blob).amplitude);
      }
      // Put blob at this position (include ghost blobs to allow for periodc boundary conditions)
      size_h2 = (*current_blob).size_h * (*current_blob).size_h;
      size_v2 = (*current_blob).size_v * (*current_blob).size_v;
      // Loop over field positions
      for (x=1; x<=state.nlon; x++)
      { for (y=1; y<=state.nlat; y++)
        { for (z=1; z<=state.nlev; z++)
          { blobz = (*current_blob).level;
            distz = state.level[z] - blobz;
            for (ghostx=-1; ghostx<=1; ghostx++)
            { if (ghostx == 0)
              { yghost_low  = -1;
                yghost_high = 1;
              }
              else
              { yghost_low  = 0;
                yghost_high = 0;
              }
              for (ghosty=yghost_low; ghosty<=yghost_high; ghosty++)
              { blobx = (*current_blob).longitude + double(ghostx) * 360.0;
                if (ghosty == 0)
                { // Not a ghost
                    bloby = (*current_blob).latitude;
                }
                else
                { if (ghosty == 1)
                  { // Ghost is mirror at North Pole
                    bloby  = 180.0 - (*current_blob).latitude;
                  }
                  else
                  { // Ghost is mirror at South Pole
                    bloby = -180.0 - (*current_blob).latitude;
                  }
                  // For ghosts, need to update the x position of the blob
                 if (blobx <= 180.0)
                  { blobx += 180.0;
                  }
                  else
                  { blobx -= 180.0;
                  }
                }
                distx    = state.longitude[x] - blobx;
                disty    = state.latitude[y] - bloby;
                disth_sq = distx * distx + disty * disty;
                distv_sq = distz * distz;
                state.tracer0_rs[x][y][z] += (*current_blob).amplitude *
                                             exp(-1.0 * (disth_sq/size_h2 + distv_sq/size_v2));
              }
            }
          }
        }
      }
      current_blob = (*current_blob).next;
    }

    // Scale the fields if necessary
    if (num_blobs_tracer > 0)
    { // Find the maximum field value
      max_field = 0.0;
      for (x=1; x<=state.nlon; x++)
      { for (y=1; y<=state.nlat; y++)
        { for (z=1; z<=state.nlev; z++)
          { if (fabs(state.tracer0_rs[x][y][z]) > max_field)
            { max_field = fabs(state.tracer0_rs[x][y][z]);
            }
          }
        }
      }
      factor = max_blob / max_field;
      // Scale
      for (x=1; x<=state.nlon; x++)
      { for (y=1; y<=state.nlat; y++)
        { for (z=1; z<=state.nlev; z++)
          { state.tracer0_rs[x][y][z] *= factor;
          }
        }
      }
    }
    printf ("Initial condition set\n");


    if (nss > 0)
    { for (ss=0; ss<nss; ss++)
      { printf ("Dealing with source/sink field %i\n", ss);

        // Make blobs in this source/sink field
        current_blob   = ss_field1[ss];
        max_blob       = 0.0;
        num_blobs_flux = 0;
        while (current_blob)
        { printf ("Adding a new blob to the flux field\n");
          num_blobs_flux++;
          // Monitor the size of the blobs
          if (fabs((*current_blob).amplitude) > max_blob)
          { max_blob = fabs((*current_blob).amplitude);
          }
          // Put blob at this position (include ghost blobs to allow for periodc boundary conditions)
          size_h2 = (*current_blob).size_h * (*current_blob).size_h;

          // Loop over field positions
          for (x=1; x<=state.nlon; x++)
          { for (y=1; y<=state.nlat; y++)
            { for (ghostx=-1; ghostx<=1; ghostx++)
              { if (ghostx == 0)
                { yghost_low  = -1;
                  yghost_high = 1;
                }
                else
                { yghost_low  = 0;
                  yghost_high = 0;
                }
                for (ghosty=yghost_low; ghosty<=yghost_high; ghosty++)
                { blobx = (*current_blob).longitude + double(ghostx) * 360.0;
                  if (ghosty == 0)
                  { // Not a ghost
                      bloby = (*current_blob).latitude;
                  }
                  else
                  { if (ghosty == 1)
                    { // Ghost is mirror at North Pole
                      bloby  = 180.0 - (*current_blob).latitude;
                    }
                    else
                    { // Ghost is mirror at South Pole
                      bloby = -180.0 - (*current_blob).latitude;
                    }
                    // For ghosts, need to update the x position of the blob
                   if (blobx <= 180.0)
                    { blobx += 180.0;
                    }
                    else
                    { blobx -= 180.0;
                    }
                  }
                  distx    = state.longitude[x] - blobx;
                  disty    = state.latitude[y] - bloby;
                  disth_sq = distx * distx + disty * disty;
                  state.source_rs[x][y][ss] += (*current_blob).amplitude *
                                                 exp(-1.0 * disth_sq/size_h2);
                }
              }
            }
          }
          current_blob = (*current_blob).next;
        }
      }

      // Scale the fields if necessary
      if (num_blobs_flux > 0)
      { // Find the maximum field value
        max_field = 0.0;
        for (x=1; x<=state.nlon; x++)
        { for (y=1; y<=state.nlat; y++)
          { for (ss=0; ss<nss; ss++)
            { if (fabs(state.source_rs[x][y][ss]) > max_field)
              { max_field = fabs(state.source_rs[x][y][ss]);
              }
            }
          }
        }
        factor = max_blob / max_field;
        // Scale
        for (x=1; x<=state.nlon; x++)
        { for (y=1; y<=state.nlat; y++)
          { for (ss=0; ss<nss; ss++)
            { state.source_rs[x][y][ss] *= factor;
            }
          }
        }
      }

      printf ("Flux set\n");

    }
    printf ("Done\n");

    // Complete the halos
    printf ("Sorting the halos\n");
    halos (&state);
    printf ("Done\n");

    // Output the fields
    WriteStateVector (&state,
                      output_file);
  }


//   ===============================================================================
//   Deallocate
//   ===============================================================================

  Deallocate_state (&state);
  Deallocate_metadata (&MetaData);
  Destroy_MakeFields_arguments ( &ic_blob1,
                                 nss,
                                 &num_ss_blobs,
                                 &ss_field1 );
  printf ("Done program\n");
}
