   #include <source.h>

// -------------------------------------------------------------------------------
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
)


{ struct state_type          xp;  // Equivalent of chi, but in model space
  struct state_type          xn;  // Current state
  int                        x, y, z, t;
  int                        gamma, nminor, nmajor, epsilon, nss_needed, sfield;
  int                        windnum, starttime;
  double                     last_ob_time, sDt, fluxfactor;
  double                     lowert, uppert, obtime, runlength;
  struct obs_type            *currentob;
  struct Wind_type           windA, windB;
  struct Wind_type           *wind_lowert, *wind_uppert, *wind_temp;
  struct instant_tracer_type field1, field2;
  struct instant_tracer_type *field_lowert, *field_uppert, *field_temp;
  struct state_type          state_hat;   // Real space
  double                     ***dummyd;
  int                        ***dummyi;
  char                       dp_file_rw;
  char                       wind_file[256];
  bool                       study_states = false;
  char                       ForwardFileName[256] = "PenAndGrad_ForwardIntegration.nc";
  char                       AdjointFileName[256] = "PenAndGrad_AdjointIntegration.nc";
  bool                       this_window;

  // Allocate xp and x
  Allocate_state (&xp,
                  MetaData,
                  'r',
                  'r',
                  'r');
  Allocate_state (&xn,
                  MetaData,
                  'r',
                  'r',
                  'r');


  // Compute the current state (xn = xb + xp)
  if (zerochi)
  { // chi is zero, so xn = xb + 0
    copy_general ( xb,
                   &xn,
                   HorizData,
                   true );
    // Cost function, background term
    *Jb = 0.0;
  }
  else
  { // Compute the perturbation in model space
    cvt_total ( chi,
                &xp,
                MetaData,
                HorizData,
                CVTdata );

    // Compute the current state (xn = xb + xp)
    add_general ( xb,
                  &xp,
                  &xn,
                  HorizData,
                  true );
    // Cost function, background term
    *Jb = 0.5 * innerproduct_general ( chi,
                                       chi,
                                       HorizData,
                                       true );
  }



  // Cost function, observation term (initialise)
  *Jo = 0.0;


  //   ===============================================================================
  //   Get ready to do the forward integration
  //   ===============================================================================

  // How many minor timesteps per major timestep?
  gamma = int(Dt/dt);
  if (verbose_output == 'v')
  { printf ("There are %i minor timesteps per major timestep\n", gamma);
  }

  // Find out the latest observation time (to determine the run length)
  last_ob_time = 0.0;
  currentob    = (*obs_trunk).first;
  while (currentob)
  { if ((*currentob).obtime_secs > last_ob_time)
    { last_ob_time = (*currentob).obtime_secs;
    }
    currentob = (*currentob).next;
  }

  if (verbose_output == 'v')
  { printf ("The last obs time: %f s, = %f days\n", last_ob_time, last_ob_time/(3600.0 * 24.0));
    printf ("This may need rounding up for the run length\n");
  }
  runlength = dt * int(last_ob_time / dt + ::nearly1);
  if (verbose_output == 'v')
  { printf ("The run length: %f s, = %f days\n", runlength, runlength/(3600.0 * 24.0));
  }

  // How many minor timesteps in the total integration?
  nminor = int(runlength / dt);
  (*MetaData).ntimes_minor = nminor;
  if (verbose_output == 'v')
  { printf ("There are %i minor timesteps in the total integration\n", nminor);
  }

  // How many major time steps are involved in this integration?
  nmajor = int(runlength / Dt + ::nearly1);
  (*MetaData).ntimes_major = nmajor;
  if (verbose_output == 'v')
  { printf ("There are %i major timesteps in this integration\n", nmajor);
    printf ("There are %i wind fields required\n", nmajor+1);
  }


  // Determine the time step between source/sink fields
  sDt = (*xb).times[1] - (*xb).times[0];
  if (verbose_output == 'v')
  { printf ("The time between surface flux fields is %f s (%f days)\n", sDt, sDt / (24.0 * 3600.0));
  }

  // How many major timesteps per flux timestep?
  epsilon = int(sDt/Dt);
  if (verbose_output == 'v')
  { printf ("There are %i major timesteps per flux timestep\n", epsilon);
  }

  // sDt must be a multiple of Dt
  if ((sDt/Dt) != double(epsilon))
  { printf ("Error: the value of sDt (%f) must be a multiple of Dt (%f)\n", sDt, Dt);
    exit(0);
  }

  // Calculate the flux factor (the factor for the flux term)
  fluxfactor = Dt / (::rho0 * ((*MetaData).level[2] - (*MetaData).level[1]));
  if (verbose_output == 'v')
  { printf ("The flux factor is %f\n", fluxfactor);
  }

  // How many surface source/sink files are needed for this integration?
  nss_needed = int(runlength / sDt + ::nearly1);
  if (verbose_output == 'v')
  { printf ("There are %i source/sink fields needed for this integration\n", nss_needed);    
    printf ("There are %i source/sink fields in the specified file\n", (*xb).nss);
  }

  if (nss_needed > (*xb).nss)
  { printf ("Error: the number of source/sink fields needed is larger than those available\n");
    exit(0);
  }

  if (verbose_output == 'v')
  { printf ("Number of levels in MetaData : %i\n", (*MetaData).nlev);
  }
  if (new_dp)
  { // Allocate space for wind data
    Allocate_wind (&windA,
                   MetaData,
                   false);
    Allocate_wind (&windB,
                   MetaData,
                   false);
    if (verbose_output == 'v')
    { printf ("Number of levels in windA    : %i\n", windA.nlev);
    }
  }

  // Allocate space for tracer data
  Allocate_instant (&field1,
                    MetaData);
  Allocate_instant (&field2,
                    MetaData);



  // Set the initial conditions
  for (x=0; x<(*xb).nlon+2; x++)
  { for (y=0; y<(*xb).nlat+2; y++)
    { for (z=0; z<(*xb).nlev+2; z++)
      { field1.tracer[x][y][z] = xn.tracer0_rs[x][y][z];
      }
    }
  }
  field1.timestep_major = 0;
  field1.timestep_minor = 0;
  field1.time           = 0.0;


  // Set up the file for the departure points
  if (new_dp)
  { if (verbose_output == 'v')
    { printf ("Output file for departure points: %s\n", file_dps);
    }
    WriteDeparturePoints (MetaData,
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
    dp_file_rw = 'w';
  }
  else
  { dp_file_rw = 'r';
  }


  if (new_dp)
  { // Read-in the first two windfields
    windnum = 0;
    sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
    if (verbose_output == 'v')
    { printf ("First wind file : %s\n", wind_file);
    }
    Read_winds (&windA, MetaData, factor_w, wind_file);
    windnum++;
    sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
    printf ("Second wind file: %s\n", wind_file);
    Read_winds (&windB, MetaData, factor_w, wind_file);

    // Set-up the initial configuration of the winds
    wind_lowert = &windA;
    wind_uppert = &windB;
  }
  else
  { wind_lowert = NULL;
    wind_uppert = NULL;
  }


  if (study_states)
  { printf ("Initialising output file for animation: %s\n", ForwardFileName);
      WriteTimeSeq (MetaData,
                    &field1,
                    nmajor+1,
                    1,
                    0,
                    0,       // 0=create file only
                    ForwardFileName);

      // Write out t=0 field
      WriteTimeSeq (MetaData,
                    &field1,
                    nmajor+1,
                    1,
                    0,
                    1,       // 1=normal (write only)
                    ForwardFileName);
    }



  // Configure the pointers to the tracer fields
  field_lowert = &field1;
  field_uppert = &field2;
  starttime    = 1;


  //   ===============================================================================
  //   Do the forward integration
  //   ===============================================================================

  for (t=starttime; t<=nmajor; t++)
  { // Run a model major timestep
    sfield = int(((double(t)-0.1) * Dt) / sDt);
    if (verbose_output == 'v')
    { printf ("Forward major timestep %i of %i, source field %i\n", t, nmajor, sfield);
    }
    // t-1 to t
    SemiLagrangian_1step ( MetaData,          // Meta data
                           &xn,               // State vector
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
                           interpolation_lc,  // Interpolation type
                           true,              // Is departure point file mentioned?
                           file_dps,          // Filename of departure points file
                           dp_file_rw );      // Read or write departure points file?
    //printf ("Done SemiLagrangian\n");


    if (study_states)
    { printf ("Outputting forward step\n");
      WriteTimeSeq (MetaData,
                    field_uppert,
                    nmajor+1,
                    1,
                    t,
                    1,       // 1=normal (write only)
                    ForwardFileName);
    }


    // Apply the observation operator for observations between times t*Dt and (t+1)*Dt
    // Notes
    // At this stage, field_lowert corresponds to the tracer field at time (t-1) * Dt
    //                field_uppert                                             t * Dt
    // This state has source field between state.source_rs[.,.,sfield] and 
    //                                     state.source_rs[.,.,sfield+1]

    // Look for the relevant observations to find the model obs, and the modeob adjoints
    lowert    = double(t-1) * Dt;
    uppert    = double(t) * Dt;
    currentob = (*obs_trunk).first;
    //printf ("Considering observations\n");
    while (currentob)
    { obtime = (*currentob).obtime_secs;
      //printf ("OB OP: Checking times %f  ** %f ** %f\n", lowert, obtime, uppert);
      if (t == starttime)
      { this_window = (obtime >= lowert) && (obtime <= uppert);
      }
      else
      { this_window = (obtime > lowert) && (obtime <= uppert);
      }

      if (((*currentob).time_index > -1) && this_window)
      { // Generate the model's version of this observation
        //printf ("Model ob, %c %f %f %f %f\n",
        //         (*currentob).ob_of,
        //         (*currentob).longitude,
        //         (*currentob).latitude,
        //         (*currentob).level,
        //         (*currentob).obtime_secs);
        ObservationOperator (currentob,
                             obs_trunk,
                             MetaData,
                             field_lowert,
                             field_uppert,
                             lowert,
                             uppert,
                             &xn);
        //printf ("Observation operator done\n");
        (*currentob).innov         = (*currentob).ob - (*currentob).model_ob;
        (*currentob).dJo_dmodel_ob = -1.0 * (*currentob).innov / (*currentob).variance;
        // Increment Jo ("-=" instead of "+=" to undo the "-1.0" above
        *Jo                       -= (*currentob).innov * (*currentob).dJo_dmodel_ob;

        //printf ("-------------------------------------\n");
        //printf ("Observation: %f %f %f %f %f\n",
        //        (*currentob).ob, (*currentob).model_ob, (*currentob).innov, (*currentob).dJo_dmodel_ob, *Jo);
        //printf ("-------------------------------------\n");
      }

      currentob = (*currentob).next;
    }

    if (t < nmajor)
    { // Bother with next winds fields only if there is still another step to do
      //printf ("Just crossed a major timestep\n");
      if (new_dp)
      { windnum++;
        // Swap pointers
        wind_temp   = wind_lowert;
        wind_lowert = wind_uppert;
        wind_uppert = wind_temp;

        // Read-in the next wind field
        sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
        if (verbose_output == 'v')
        { printf ("Wind file: %s\n", wind_file);
        }
        Read_winds (wind_uppert, MetaData, factor_w, wind_file);
      }

      // Sort out the configuration of the instantateous fields for the next step
      field_temp   = field_lowert;
      field_lowert = field_uppert;
      field_uppert = field_temp;
    }
  }

  *Jo *= 0.5;
  // Compute the total cost
  *J = *Jb + *Jo;


  if (calc_grad)
  { //   ===============================================================================
    //   Get ready to do the adjoint integration
    //   ===============================================================================

    Allocate_state ( &state_hat,
                     MetaData,
                     'r', 'r', 'r' );

    // Set the adjoint fields Xhat(T) to zero
    for (x=0; x<(*MetaData).nlon+2; x++)
    { for (y=0; y<(*MetaData).nlat+2; y++)
      { for (z=0; z<(*MetaData).nlev+2; z++)
        { field1.tracer[x][y][z] = 0.0;
          field2.tracer[x][y][z] = 0.0;
        }
        for (t=0; t<(*MetaData).nss; t++)
        { state_hat.source_rs[x][y][t] = 0.0;
        }
      }
    }
    // Initialise pointers
    field_lowert = &field1;
    field_uppert = &field2;


    if (study_states)
    { printf ("Initialising output file for animation: %s\n", AdjointFileName);
      WriteTimeSeq (MetaData,
                    &field1,
                    nmajor+1,
                    1,
                    0,
                    0,       // 0=create file only
                    AdjointFileName);
    }


    //   ===============================================================================
    //   Do the adjoint integration
    //   ===============================================================================

    for (t=nmajor; t>=starttime; t--)
    { sfield = int(((double(t)-0.1) * Dt) / sDt);
      if (verbose_output == 'v')
      { printf ("Adjoint major timestep %i of %i, source field %i\n", t, nmajor, sfield);
      }

      // Adjoint of observation operator
      // Look for the relevant observations to find the model obs, and the modeob adjoints
      lowert    = double(t-1) * Dt;
      uppert    = double(t) * Dt;
      currentob = (*obs_trunk).first;
      //printf ("Considering observations (adjoint)\n");
      while (currentob)
      { obtime = (*currentob).obtime_secs;
        // printf ("Checking times %f  ** %f ** %f\n", lowert, obtime, uppert);
        if (t == starttime)
        { this_window = (obtime >= lowert) && (obtime <= uppert);
        }
        else
        { this_window = (obtime > lowert) && (obtime <= uppert);
        }

        if (this_window)
        { ObservationOperator_adj (currentob,
                                   obs_trunk,
                                   MetaData,
                                   field_lowert,
                                   field_uppert,
                                   &state_hat);
        }
        currentob = (*currentob).next;
      }

      if (study_states && (t==nmajor))
      { printf ("Outputting top adjoint step\n");
        WriteTimeSeq (MetaData,
                      field_uppert,
                      nmajor+1,
                      1,
                      t,
                      1,       // 1=normal (write only)
                      AdjointFileName);
      }

      // Run a model major timestep (adjoint)
      // The lambda_t in the notes corresponds to: field_lowert (t-1)
      //                                           field_uppert (t)
      //                                           state_hat (fluxes)
      //printf ("Running SemiLagrangian (adjoint)\n");
      // Adjoint of: t-1 to t
      // Acting on lambda_t+1
      SemiLagrangian_1step_adj ( MetaData,          // Meta data
                                 &state_hat,        // inout: State vector
                                 field_lowert,      // inout: Tracer field at major t-1
                                 field_uppert,      // in:    Tracer field extrapolated to major t
                                 t,                 // Major timestep number (upper time)
                                 sfield,            // Source/sink field number
                                 Dt,                // Major timestep size
                                 kappa_dt,          // Diffusion timestep size
                                 kappa_h,           // Horiz diffusion coefficient
                                 kappa_v,           // Vert diffusion coefficient
                                 inc_adv,           // Include advection
                                 fluxfactor,        // Flux factor
                                 inc_vert,          // Include vertical transport?
                                 interpolation_lc,  // Interpolation type
                                 file_dps );        // Filename of departure points file
      //printf ("Done SemiLagrangian (adjoint)\n");

      if (study_states)
      { printf ("Outputting adjoint step\n");
        WriteTimeSeq (MetaData,
                      field_lowert,
                      nmajor+1,
                      1,
                      t-1,
                      1,       // 1=normal (write only)
                      AdjointFileName);
      }


      // Sort out the configuration of the instantateous fields for the next step
      field_temp   = field_lowert;
      field_lowert = field_uppert;
      field_uppert = field_temp;

      // Set field_lowert to zero
      for (x=0; x<(*MetaData).nlon+2; x++)
      { for (y=0; y<(*MetaData).nlat+2; y++)
        { for (z=0; z<(*MetaData).nlev+2; z++)
          { (*field_lowert).tracer[x][y][z] = 0.0;
          }
        }
      }
    }

    // Put tracer information into state gradient (real space)
    // state_hat.tracer has not been touched yet
    for (x=0; x<(*MetaData).nlon+2; x++)
    { for (y=0; y<(*MetaData).nlat+2; y++)
      { for (z=0; z<(*MetaData).nlev+2; z++)
        { state_hat.tracer0_rs[x][y][z] = (*field_uppert).tracer[x][y][z];
        }
      }
    }


    // Do adjoint of the CVT to get gradient in control space
    cvt_total_adj (grad,           // out
                   &state_hat,     // in
                   MetaData,
                   HorizData,
                   CVTdata);

    // Add on the background gradient
    add_general ( grad,
                  chi,
                  grad,
                  HorizData,
                  true );

    // Tidy up structures needed for the adjoint only
    Deallocate_state (&state_hat);

  }

  // Tidy up general structures
  Deallocate_state (&xp);
  Deallocate_state (&xn);

  if (new_dp)
  { Deallocate_wind (&windA);
    Deallocate_wind (&windB);
  }
  Deallocate_instant (&field1, MetaData);
  Deallocate_instant (&field2, MetaData);
}
