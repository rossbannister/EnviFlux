/* ==================================================================================
   3d source and sink code
   Take a state file (ics + fluxes) and a forecast file
   Make a new state file:
     New ic: a chosen time taken from the forecast file
     Same fluxes

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_ReplaceFirstField.out

   Modification history
   --------------------
   27/11/22 New Code. Ross Bannister

   Documentation
   -------------
   Run with
   ./Master_ReplaceFirstField.out <Input state filename>
                                  <Input forecast filename>
                                  <Forecast timestep number (1 is first)>
                                  <Output state filename>

   =============================================================================== */

   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )
{ char                       State_in_filename[256];
  char                       Forecast_in_filename[256];
  int                        FcTime;
  char                       State_out_filename[256];
  struct state_type          state;
  struct metadata_type       MetaData;
  struct instant_tracer_type field;
  char                       x, y, z;
  bool                       ok;
  double                     mean, dev, stddev, npoints;

  // ===============================================================================
  // Initialisation
  // ===============================================================================

  ok = ReplaceFirstField_arguments ( argument_count,
                                     argument_list,
                                     State_in_filename,
                                     Forecast_in_filename,
                                     &FcTime,
                                     State_out_filename );

  if (ok)
  { // The following are done just to get the metadata (longs, lats, levs, etc)
    // Initialise and read in generic data needed for the horizontal transform

    // Read-in the state (ics + fluxes)
    printf ("Reading state from file\n");
    state.horiz_repres = 'r';
    state.vert_repres  = 'r';
    state.temp_repres  = 'r';

    ReadStateVector (&state,
                     &MetaData,
                     true,
                     true,
                     State_in_filename);
    printf ("Done\n");

    // Allocate a field for reading in
    Allocate_instant (&field,
                      &MetaData);

    // Read-in the field
    printf ("Reading in the replacement field\n");
    Read1Time ( &MetaData,             // In  (to check against)
                &field,                // Out (assumed already allocated)
                FcTime,                // In  Time index requested (1 is first)
                Forecast_in_filename,  // In  Filename
                &ok );                 // Out Successful execution flag
    printf ("Done\n");

    if (ok)
    { // Overwrite the ic of the state (including halos)
      printf ("Overwriting ic\n");
      printf ("nlon = %i\n", MetaData.nlon);
      printf ("nlat = %i\n", MetaData.nlat);
      printf ("nlev = %i\n", MetaData.nlev);
      printf ("nss  = %i\n", MetaData.nss);

      for (x=0; x<MetaData.nlon+2; x++)
      { for (y=0; y<MetaData.nlat+2; y++)
        { for (z=0; z<MetaData.nlev+2; z++)
          { state.tracer0_rs[x][y][z] = field.tracer[x][y][z];
          }
        }
      }
      printf ("Done\n");

      // Calculate some statistics of the state (for information)
      npoints = double(MetaData.nlon * MetaData.nlat);
      for (z=1; z<MetaData.nlev+1; z++)
      { mean = 0.0;
        for (x=1; x<MetaData.nlon+1; x++)
        { for (y=1; y<MetaData.nlat+1; y++)
          { mean += state.tracer0_rs[x][y][z];
          }
        }
        mean /= npoints;
        stddev = 0.0;
        for (x=1; x<MetaData.nlon+1; x++)
        { for (y=1; y<MetaData.nlat+1; y++)
          { dev     = state.tracer0_rs[x][y][z] - mean;
            stddev += dev * dev;
          }
        }
        stddev = sqrt(stddev / npoints);
        printf ("Level %i, mean and stddev: %f %f\n", z, mean, stddev);
      }

      // Output the new state
      printf ("Outputting the new state\n");
      WriteStateVector (&state,
                        State_out_filename);
      printf ("Done\n");
    }

    // Tidy up (deallocate field, MetaData, and state)
    printf ("Tidying up\n");
    Deallocate_instant (&field,
                        &MetaData);
    Deallocate_metadata (&MetaData);
    Deallocate_state (&state);
    printf ("Done\n");
  }

  printf ("Done program\n");
}
