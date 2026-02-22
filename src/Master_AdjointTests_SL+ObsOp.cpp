/* ==================================================================================
   3d source and sink code
   Adjoint tests for routines associated with the transport model and observation operator

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_AdjointTests_SL+ObsOp.out

   Modification history
   --------------------
   30/08/23 Adapted from Master_AdjointTests_SemiLagrangian and Master_AdjointTestsObOps.
            Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <stdio.h>
   #include <math.h>
   #include <stdlib.h>
   #include <string.h>
   #include <source.h>


int main ()
{ char                       state_ic[256] = "../data/State_56levs.nc";
  char                       wind_dir[256] = "../data/ECMWF_RedResWinds_1995_56levs";
  char                       obsfile_in[256] = "Observations.dat";
  char                       obs_file_out[256] = "Obs_AdjTest.dat";
  char                       file_dps[256] = "DeparturePoints.nc";
  char                       output_file_anim[256] = "AdjTest.nc";
  char                       output_file_anim_adj[256] = "AdjTest_adj.nc";
  char                       output_file_diags[256] = "AdjTest.dat";
  char                       wind_file[256];
  double                     Dt = 43200.0;
  double                     dt = 3600.0;
  double                     kappa_dt = 1800.0;
  double                     kappa_h = 0.0;  //4000000.0;
  double                     kappa_v = 0.0;  //2.0;
  char                       interpolate_lc = 'c';
  double                     factor_w = 1.0;
  bool                       output_dps = true;
  bool                       inc_adv = true;
  bool                       inc_vert = true;
  double                     output_freq = 86400.0;
  char                       dp_file_rw = 'w';    // Read or write departure points file?
                                                  // Write: calc departure points in forward SL scheme and output
                                                  // Read: read-in previously calc departure points
  bool                       output_anim = true;
  bool                       output_diags = true;
  int                        gamma, nminor, nmajor, epsilon, nss_needed, op_freq, n_op_times;
  int                        starttime, t, sfield, x, y, z;
  int                        nlon, nlat, nlev, windnum;
  struct metadata_type       MetaData;
  struct state_type          state, state_hat;
  struct obs_trunk_type      obs_trunk;
  double                     last_ob_time, runlength, sDt, fluxfactor, max_val, min_val;
  double                     lowert, uppert, obtime;
  double                     maxmax_val, minmin_val;
  double                     model_ob, lhs, rhs;
  struct obs_type            *currentob;
  struct Wind_type           windA, windB, *wind_lowert, *wind_uppert, *wind_temp;
  struct instant_tracer_type field1, field2, *field_lowert, *field_uppert, *field_temp, field_uppert_store;
  bool                       mention_dps;
  double                     ***dummyd;
  int                        ***dummyi;
  FILE*                      diags_file = NULL;
  int                        rndseq = 985395;   // Random number seed



  // ===============================================================================

  // Read-in the state (containing the initial conditions and the surface fluxes)
  printf ("Reading-in the state vector (initial conditions and surface flux time sequence)\n");
  state.horiz_repres = 'r';
  state.vert_repres  = 'r';
  state.temp_repres  = 'r';
  ReadStateVector (&state,
                   &MetaData,
                   true, true,
                   state_ic);

  printf ("Reading-in done\n");
  nlon = MetaData.nlon;
  nlat = MetaData.nlat;
  nlev = MetaData.nlev;
  printf ("Dimensions (nlon, nlat, nlev: %i, %i, %i)\n", nlon, nlat, nlev);

/*
  // Put random numbers into state vector
  for (x=1; x<=state.nlon; x++)
  { for (y=1; y<=state.nlat; y++)
    { for (z=1; z<=state.nlev; z++)
      { state.tracer0_rs[x][y][z] = Normal (150.0, 10.0, &rndseq);
      }
      for (t=0; t<state.nss; t++)
      { state.source_rs[x][y][t]  = Normal (0.0001, 0.00001, &rndseq);
      }
    }
  }
  halos (&state);
*/


  // Read-in observations
  printf ("Reading observations\n");
  ReadObservations (&MetaData,
                    &obs_trunk,
                    obsfile_in);

  // How many minor timesteps per major timestep?
  gamma = int(Dt/dt);
  printf ("There are %i minor timesteps per major timestep\n", gamma);

  // Determine the time step between source/sink fields
  sDt = state.times[1] - state.times[0];
  printf ("The time between surface flux fields is %f s (%f days)\n", sDt, sDt / (24.0 * 3600.0));

  // Find out the latest observation time (to determine the run length)
  last_ob_time = 0.0;
  currentob    = obs_trunk.first;
  while (currentob)
  { if ((*currentob).obtime_secs > last_ob_time)
    { last_ob_time = (*currentob).obtime_secs;
    }
    currentob = (*currentob).next;
  }
  printf ("The last obs time: %f s, = %f days\n", last_ob_time, last_ob_time/(3600.0 * 24.0));
  printf ("This may need rounding up for the run length\n");
  runlength = Dt * int(last_ob_time / Dt + ::nearly1);
  printf ("The run length: %f s, = %f days\n", runlength, runlength/(3600.0 * 24.0));

  // How many minor timesteps in the total integration?
  nminor = int(runlength / dt);
  MetaData.ntimes_minor = nminor;
  printf ("There are %i minor timesteps in the total integration\n", nminor);

  // How many major time steps are involved in this integration?
  nmajor = int(runlength / Dt + ::nearly1);
  MetaData.ntimes_major = nmajor;
  printf ("There are %i major timesteps in this integration\n", nmajor);
  printf ("There are %i wind fields required\n", nmajor+1);

  // How many major timesteps per flux timestep?
  epsilon = int(sDt/Dt);
  printf ("There are %i major timesteps per flux timestep\n", epsilon);

  // sDt must be a multiple of Dt
  if ((sDt/Dt) != double(epsilon))
  { printf ("Error: the value of sDt (%f) must be a multiple of Dt (%f)\n", sDt, Dt);
  }

  // Calculate the flux factor (the factor for the flux term)
  fluxfactor = Dt / (::rho0 * (MetaData.level[2] - MetaData.level[1]));
  printf ("The flux factor is %f\n", fluxfactor);

  // How many surface source/sink files are needed for this integration?
  nss_needed = int(runlength / sDt + ::nearly1);
  printf ("There are %i source/sink fields needed for this integration\n", nss_needed);    
  printf ("There are %i source/sink fields in the specified file\n", state.nss);

  if (nss_needed > state.nss)
  { printf ("Error: the number of source/sink fields needed is larger than those available\n");
    exit(0);
  }

  // What is the output frequency measured in timesteps?
  op_freq = int(output_freq / Dt);
  printf ("There are %i major timesteps between each output time\n", op_freq);
  if (op_freq == 0)
  { printf ("Error: Dt (%f) must be <= output_freq (%f)\n", Dt, output_freq);
    exit(0);
  }

  // How many output times will there be (1 is added because the initial conditions are reproduced in the first)?
  n_op_times = int(runlength / output_freq) + 1;
  printf ("There are %i output times\n", n_op_times);


  // Allocate space for data
  printf ("Allocating windA\n");
  Allocate_wind (&windA,
                 &MetaData,
                 false);
  printf ("Allocating windB\n");
  Allocate_wind (&windB,
                 &MetaData,
                 false);

  printf ("Allocating field1\n");
  Allocate_instant (&field1,
                    &MetaData);
  printf ("Allocating field2\n");
  Allocate_instant (&field2,
                    &MetaData);
  printf ("Allocating field_uppert_store\n");
  Allocate_instant (&field_uppert_store,
                    &MetaData);

  printf ("Done allocations\n");

  printf ("Number of levels in MetaData : %i\n", MetaData.nlev);
  printf ("Number of levels in windA    : %i\n", windA.nlev);

  mention_dps  = (strcmp ("nil", file_dps) != 0);           // Mention departure points file?

  if (output_anim)  { printf ("Will output animation file\n"); }
  if (output_diags) { printf ("Will output diagnostics file\n"); }
  if (mention_dps)  { printf ("Have mentioned departure points file\n"); }

  if ((dp_file_rw == 'r') && !mention_dps)
  { printf ("Departure points file has not been specified (is nil), but should be\n");
    exit(0);
  }

  // Set the initial conditions
  for (x=0; x<state.nlon+2; x++)
  { for (y=0; y<state.nlat+2; y++)
    { for (z=0; z<state.nlev+2; z++)
      { field1.tracer[x][y][z] = state.tracer0_rs[x][y][z];
      }
    }
  }
  field1.timestep_major = 0;
  field1.timestep_minor = 0;
  field1.time           = 0.0;

  // Set up the file for the departure points
  if (mention_dps && (dp_file_rw == 'w'))
  { printf ("Output file for departure points: %s\n", file_dps);
    WriteDeparturePoints (&MetaData,
                          dummyd,
                          dummyd,
                          dummyd,
                          dummyi,
                          dummyi,
                          dummyi,
                          nmajor,
                          Dt,
                          0,
                          0,       // 0=create file only,
                          file_dps);
  }


  // Maximum and minimum values
  maxmin_tracer (&MetaData,
                 &field1,
                 true,
                 &max_val,
                 &min_val);
  maxmax_val = max_val;
  minmin_val = min_val;

//   ===============================================================================
//   Set-up for forward operators
//   ===============================================================================

  // Set up the output file for the fields
  if (output_anim)
  { printf ("Output file for animation: %s\n", output_file_anim);
    WriteTimeSeq (&MetaData,
                  &field1,
                  n_op_times+1,
                  output_freq,
                  0,
                  0,       // 0=create file only
                  output_file_anim);

    // Write out t=0 field
    WriteTimeSeq (&MetaData,
                  &field1,
                  n_op_times+1,
                  output_freq,
                  0,
                  1,       // 1=normal (write only)
                  output_file_anim);
  }


  // Open output file for diagnostics
  if (output_diags)
  { printf ("Output file for diagnostics: %s\n", output_file_diags);
    diags_file = fopen (output_file_diags, "w");

    // Output the first maximum and mininmum values of the tracer field at this time
    fprintf (diags_file, "# Maximum and minimum values of the integration\n");
    fprintf (diags_file, "# dt         = %f\n", dt);
    fprintf (diags_file, "# kappa_h    = %f\n", kappa_h);
    fprintf (diags_file, "# kappa_v    = %f\n", kappa_v);
    fprintf (diags_file, "%f  %f\n", min_val, max_val);
  }

  // Read-in the first two windfields
  windnum = 0;
  sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
  printf ("First wind file : %s\n", wind_file);
  Read_winds (&windA, &MetaData, factor_w, wind_file);
  windnum++;
  sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
  printf ("Second wind file: %s\n", wind_file);
  Read_winds (&windB, &MetaData, factor_w, wind_file);

  // Set-up the initial configuration of the winds
  wind_lowert = &windA;
  wind_uppert = &windB;

  // Configure the pointers to the tracer fields
  field_lowert = &field1;
  field_uppert = &field2;
  starttime    = 1;


//   ===============================================================================
//   Start to run the integration of the transport model
//   ===============================================================================

  for (t=starttime; t<=nmajor; t++)
  { sfield = int(((double(t)-0.1) * Dt) / sDt);
    printf ("Major timestep %i of %i, source field %i\n", t, nmajor, sfield);
    // Run a model major timestep
    printf ("Running SemiLagrangian\n");
    SemiLagrangian_1step ( &MetaData,         // Meta data
                           &state,            // State vector
                           field_lowert,      // Tracer field at major t-1
                           field_uppert,      // Tracer field extrapolated to major t
                           t,                 // Major timestep number (new time)
                           sfield,            // Source/sink field number
                           Dt,                // Major timestep size
                           dt,                // Minor timestep size
                           kappa_dt,          // Diffusion timestep size
                           gamma,             // Number of minor timesteps per major timestep
                           kappa_h,           // Horiz diffusion coefficient
                           kappa_v,           // Vert diffusion coefficient
                           inc_adv,           // Include advection
                           fluxfactor,        // Flux factor
                           inc_vert,          // Include vertical transport?
                           wind_lowert,       // Wind fields at lower major timestep
                           wind_uppert,       // Wind fields at next major step
                           interpolate_lc,    // Interpolation type
                           mention_dps,       // Is departure point file mentioned?
                           file_dps,          // Filename of departure points file
                           dp_file_rw );      // Read or write departure points file?
    printf ("Done SemiLagrangian\n");

    // Apply the observation operator for observations between times t*Dt and (t+1)*Dt
    // Notes
    // At this stage, field_lowert corresponds to the tracer field at time t * Dt
    //                field_uppert                                     (t+1) * Dt
    // This state has source field between state.source_rs[.,.,sfield] and 
    //                                     state.source_rs[.,.,sfield+1]



    // Look for the relevant observations
    lowert    = double(t-1) * Dt;
    uppert    = double(t) * Dt;
    currentob = obs_trunk.first;
    printf ("Considering observations\n");
    while (currentob)
    { obtime = (*currentob).obtime_secs;
      printf ("Checking times %f  ** %f ** %f\n", lowert, obtime, uppert);
      if ((obtime >= lowert) && (obtime < uppert))
      { //Generate the model's version of this observation
        //printf ("Generating observation, %c %f %f %f\n",
        //          (*currentob).ob_of,
        //          (*currentob).longitude,
        //          (*currentob).latitude,
        //          (*currentob).level);
        ObservationOperator (currentob,
                             &obs_trunk,
                             &MetaData,
                             field_lowert,
                             field_uppert,
                             lowert,
                             uppert,
                             &state);
        // Pass model observation to the adjoint element (in prep for adjoint calculation later)
        (*currentob).dJo_dmodel_ob = (*currentob).model_ob;

        /*printf ("FORWARD: Found ob at t %f:  %i %f %f   %i %f %f   %i %f %f\n",
                (*currentob).obtime_secs,
                (*currentob).lon_index,
                (*currentob).lon_alpha,
                (*currentob).lon_beta,
                (*currentob).lat_index,
                (*currentob).lat_alpha,
                (*currentob).lat_beta,
                (*currentob).lev_index,
                (*currentob).lev_alpha,
                (*currentob).lev_beta );
        */
      }
      currentob = (*currentob).next;
    }


    // Output this timestep
    if (output_anim)
    { if ((t % op_freq) == 0)
      { printf ("Outputting data step %i\n", int(double(t)/double(op_freq)));
        WriteTimeSeq (&MetaData,
                      field_uppert,
                      n_op_times+1,
                      output_freq,
                      int(double(t)/double(op_freq)),
                      1,            // 1=normal (write only)
                      output_file_anim);
      }
    }

    // Maximum and minimum values
    maxmin_tracer (&MetaData,
                   field_uppert,
                   true,
                   &max_val,
                   &min_val);
    if (max_val > maxmax_val)
    { maxmax_val = max_val;
    }
    if (min_val < minmin_val)
    { minmin_val = min_val;
    }

    if (output_diags)
    { fprintf (diags_file, "%f  %f  %f\n", double(t) * dt, max_val, min_val);
    }


    if (t < nmajor)
    { // Bother with next winds fields only if there is still another step to do
      printf ("Just crossed a major timestep\n");
      windnum++;
      // Swap pointers
      wind_temp   = wind_lowert;
      wind_lowert = wind_uppert;
      wind_uppert = wind_temp;

      // Read-in the next wind field
      sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
      printf ("Wind file: %s\n", wind_file);
      Read_winds (wind_uppert, &MetaData, factor_w, wind_file);

      // Sort out the configuration of the instantateous fields for the next step
      field_temp   = field_lowert;
      field_lowert = field_uppert;
      field_uppert = field_temp;
    }
  }

  printf ("Smallest value encountered = %f\n", minmin_val);
  printf ("Largest value encountered  = %f\n", maxmax_val);


  printf ("Writing observations\n");
  WriteObservations ( &MetaData,
                      &obs_trunk,
                      obs_file_out );
  printf ("Done writing observations\n");



  // Store the last time step
  for (x=0; x<state.nlon+2; x++)
  { for (y=0; y<state.nlat+2; y++)
    { for (z=0; z<state.nlev+2; z++)
      { field_uppert_store.tracer[x][y][z] = (*field_uppert).tracer[x][y][z];
      }
    }
  }


//   ===============================================================================
//   Set-up for adjoint operators
//   ===============================================================================

  // Set-up the adjoint state
  printf ("Allocating the adj state\n");
  Allocate_state (&state_hat,
                  &MetaData,
                  'r',
                  'r',
                  'r');
  printf ("Done\n");

  printf ("Initialising the adj state\n");
  // Initialise the adjoint states to zero
  for (x=0; x<state.nlon+2; x++)
  { for (y=0; y<state.nlat+2; y++)
    { for (z=0; z<state.nlev+2; z++)
      { (*field_uppert).tracer[x][y][z] = 0.0;
        (*field_lowert).tracer[x][y][z] = 0.0;
        state_hat.tracer0_rs[x][y][z]   = 0.0;
      }
      for (t=0; t<state.nss; t++)
      { state_hat.source_rs[x][y][t]    = 0.0;
      }
    }
  }
  printf ("Done\n");

  // Set up the output file for the adjoint integration fields
  if (output_anim)
  { printf ("Output file for adjoint animation: %s\n", output_file_anim_adj);
    printf ("%i\n", n_op_times);
    printf ("%s\n", output_file_anim_adj);
    WriteTimeSeq (&MetaData,
                  field_uppert,
                  n_op_times+1,
                  Dt,
                  0,
                  0,       // 0=create file only
                  output_file_anim_adj);
    printf ("Done\n");
  }




  // ===============================================================================
  // Run the adjoint transport model
  // ===============================================================================

  for (t=nmajor; t>=starttime; t--)
  { sfield = int((fabs(double(t)-0.1) * Dt) / sDt);
    printf ("Major adjoint timestep %i of %i, source field %i\n", t, nmajor, sfield);


    if (t < nmajor)
    { printf ("Just crossed a major timestep in adjoint\n");
      // Sort out the configuration of the instantateous fields for the next step
      field_temp   = field_lowert;
      field_lowert = field_uppert;
      field_uppert = field_temp;
      // Set lower field to zero
      for (x=0; x<state.nlon+2; x++)
      { for (y=0; y<state.nlat+2; y++)
        { for (z=0; z<state.nlev+2; z++)
          { (*field_lowert).tracer[x][y][z] = 0.0;
          }
        }
      }

    }


    // Look for the relevant observations
    lowert    = double(t-1) * Dt;
    uppert    = double(t) * Dt;
    currentob = obs_trunk.first;
    printf ("Considering observations\n");
    while (currentob)
    { obtime = (*currentob).obtime_secs;
      //printf ("Checking times %f  ** %f ** %f\n", lowert, obtime, uppert);
      if ((obtime >= lowert) && (obtime < uppert))
      { ObservationOperator_adj (currentob,
                                 &obs_trunk,
                                 &MetaData,
                                 field_lowert,
                                 field_uppert,
                                 &state_hat);
        /*printf ("ADJOINT: Found ob at t %f:  %i %f %f   %i %f %f   %i %f %f\n",
                (*currentob).obtime_secs,
                (*currentob).lon_index,
                (*currentob).lon_alpha,
                (*currentob).lon_beta,
                (*currentob).lat_index,
                (*currentob).lat_alpha,
                (*currentob).lat_beta,
                (*currentob).lev_index,
                (*currentob).lev_alpha,
                (*currentob).lev_beta );
        */
      }
      currentob = (*currentob).next;
    }


    // Output this timestep
    if (output_anim)
    { if ((t % op_freq) == 0)
      { printf ("Outputting data step %i\n", int(double(t)/double(op_freq)));
        WriteTimeSeq (&MetaData,
                      field_uppert,
                      n_op_times+1,
                      output_freq,
                      int(double(t)/double(op_freq)),
                      1,            // 1=normal (write only)
                      output_file_anim_adj);
      }
    }


    // Run a model major timestep
    SemiLagrangian_1step_adj ( &MetaData,         // Meta data
                               &state_hat,        // inout: State vector
                               field_lowert,      // inout: Tracer field at major t-1
                               field_uppert,      // in:    Tracer field extrapolated to major t
                               t,                 // Major timestep number (new time)
                               sfield,            // Source/sink field number
                               Dt,                // Major timestep size
                               kappa_dt,          // Diffusion timestep size
                               kappa_h,           // Horiz diffusion coefficient
                               kappa_v,           // Vert diffusion coefficient
                               inc_adv,           // Include advection
                               fluxfactor,        // Flux factor
                               inc_vert,          // Include vertical transport?
                               interpolate_lc,    // Interpolation type
                               file_dps );        // Filename of departure points file
  }


  if (output_anim)
  { // Write out t=0 field
    WriteTimeSeq (&MetaData,
                  &field1,
                  n_op_times+1,
                  output_freq,
                  0,
                  1,       // 1=normal (write only)
                  output_file_anim_adj);
  }


  // Set the initial conditions
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { state_hat.tracer0_rs[x][y][z] = (*field_lowert).tracer[x][y][z];
      }
    }
  }


  // ===============================================================================
  // Do adjoint test calculations
  // ===============================================================================


  // Compute the LHS (this is in observation space)
  lhs = 0.0;
  currentob = obs_trunk.first;
  while (currentob)
  { model_ob = (*currentob).model_ob;
    if ((*currentob).dJo_dmodel_ob > ::notdef)
    { //printf ("lhs calculation : time of ob %f\n", (*currentob).obtime_secs);
      lhs += model_ob * model_ob;
    }
    currentob = (*currentob).next;
  }



/*
  // Compute the LHS
  lhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { lhs += field_uppert_store.tracer[x][y][z] * field_uppert_store.tracer[x][y][z];
      }
    }
  }
*/

  // Compute the RHS (this is in state space)
  rhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { rhs += state.tracer0_rs[x][y][z] * state_hat.tracer0_rs[x][y][z];
      }
    }
  }
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (t=0; t<=nss_needed; t++)
      { if ((x==0) && (y==0))
        { printf ("t = %i\n", t);
        }
        rhs += state.source_rs[x][y][t] * state_hat.source_rs[x][y][t];
      }
    }
  }

  printf ("In standard notation\n");
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);
  printf ("In exponential notation\n");
  printf ("LHS = %e\n", lhs);
  printf ("RHS = %e\n", rhs);


  if (output_diags)
  { fclose (diags_file);
  }


  // Tidy up
  Deallocate_state (&state);
  Deallocate_state (&state_hat);
  Deallocate_wind (&windA);
  Deallocate_wind (&windB);
  Deallocate_instant (&field1,
                      &MetaData);
  Deallocate_instant (&field2,
                      &MetaData);
  Deallocate_instant (&field_uppert_store,
                      &MetaData);
  Deallocate_metadata (&MetaData);
  delete[] obs_trunk.mass_profile;
  Destroy_Obs (&(obs_trunk.first));

}
