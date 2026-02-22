/* ==================================================================================
   3d source and sink code
   Adjoint tests of observation operators

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_AdjointTestsObOps.out

   Modification history
   --------------------
   14/08/23 New Code. Ross Bannister

   Documentation
   -------------

   =============================================================================== */

   #include <source.h>

int main ()
{ struct metadata_type       MetaData;
  struct state_type          state, state_hat;
  struct obs_trunk_type      obs_trunk;
  struct instant_tracer_type field1, field2, field1_hat, field2_hat;
  struct obs_type            *ob;
  int                        lon, lat, lev, t, c;
  int                        rndseq = 985395;   // Random number seed
  double                     inter, lhs, rhs, val;
  char                       state_ic[256] = "../data/State_56levs.nc";
  char                       obsfile[256] = "Obs.dat";
  char                       obsfile_copy[256] = "OneOb_copy.dat";

  // Read-in a state vector
  state.horiz_repres = 'r';
  state.vert_repres  = 'r';
  state.temp_repres  = 'r';
  ReadStateVector (&state,
                   &MetaData,
                   true, true,
                   state_ic);

  // Put random numbers into state vector
  for (lon=1; lon<=state.nlon; lon++)
  { for (lat=1; lat<=state.nlat; lat++)
    { for (lev=1; lev<=state.nlev; lev++)
      { state.tracer0_rs[lon][lat][lev] = Normal(150.0, 10.0, &rndseq);
      }
      for (t=0; t<state.nss; t++)
      { state.source_rs[lon][lat][t]  = Normal(0.0001, 0.00001, &rndseq);
      }
    }
  }
  halos (&state);




  // Create another state vector
  Allocate_state (&state_hat,
                  &MetaData,
                  'r', 'r', 'r');


  // Read-in observations
  printf ("Reading observations\n");
  ReadObservations (&MetaData,
                    &obs_trunk,
                    obsfile);

  // Write-out obervations as a check that they have been read-in correctly
  /*
  printf ("Writing observations\n");
  WriteObservations (&MetaData,
                     &obs_trunk,
                     obsfile_copy);
  */

  // Allocate fields
  Allocate_instant (&field1,
                    &MetaData);
  Allocate_instant (&field2,
                    &MetaData);
  Allocate_instant (&field1_hat,
                    &MetaData);
  Allocate_instant (&field2_hat,
                    &MetaData);




/*

  // Adjoint test of ob_op_interp_tracer_2D
  // --------------------------------------
  printf ("Adjoint test of ob_op_interp_tracer_2D\n");
  // Setup the input state, and initialise the output state
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { field1.tracer    [lon][lat][lev] = state.tracer0_rs[lon][lat][lev];
        field1_hat.tracer[lon][lat][lev] = 0.0;
      }
    }
  }
  // Call the interpolation with the first observation read in
  ob    = obs_trunk.first;
  inter = ob_op_interp_tracer_2D (ob,
                                  &field1,
                                  (*ob).lev_index);

  // Call the adjoint
  ob_op_interp_tracer_2D_adj (ob,
                              &field1_hat,
                              (*ob).lev_index,
                              inter);

  lhs = inter * inter;
  rhs = 0.0;
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { rhs += field1.tracer[lon][lat][lev] * field1_hat.tracer[lon][lat][lev];
      }
    }
  }

  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);



  // Adjoint test of ob_op_interp_flux_2D
  // ------------------------------------
  printf ("Adjoint test of ob_op_interp_flux_2D\n");
  // Setup the input state, and initialise the output state
  // Note this test uses a tracer as a flux proxy, level is interpreted as time
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { field1.tracer    [lon][lat][lev] = state.tracer0_rs[lon][lat][lev];
        field1_hat.tracer[lon][lat][lev] = 0.0;
      }
    }
  }
  // Call the interpolation with the first observation read in
  ob    = obs_trunk.first;
  inter = ob_op_interp_flux_2D (ob,
                                field1.tracer,
                                (*ob).time_index);

  // Call the adjoint
  ob_op_interp_flux_2D_adj (ob,
                            field1_hat.tracer, // inout
                            (*ob).time_index,
                            inter);        // in


  lhs = inter * inter;
  rhs = 0.0;
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { rhs += field1.tracer[lon][lat][lev] * field1_hat.tracer[lon][lat][lev];
      }
    }
  }

  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);




  // Adjoint test of ob_op_interp_tracer_3D
  // --------------------------------------
  printf ("Adjoint test of ob_op_interp_tracer_3D\n");
  // Setup the input state, and initialise the output state
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { field1.tracer   [lon][lat][lev] = state.tracer0_rs[lon][lat][lev];
        field1_hat.tracer[lon][lat][lev] = 0.0;
      }
    }
  }
  // Call the interpolation with the first observation read in
  ob    = obs_trunk.first;
  inter = ob_op_interp_tracer_3D (ob,
                                  &field1); // in

  // Call the adjoint
  ob_op_interp_tracer_3D_adj (ob,
                              &field1_hat,  // inout
                              inter );  // in

  lhs = inter * inter;
  rhs = 0.0;
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { rhs += field1.tracer[lon][lat][lev] * field1_hat.tracer[lon][lat][lev];
      }
    }
  }

  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

*/

/*
  // Version 1 of the Adjoint test: treat each observation independently
  // ===================================================================

  ob = obs_trunk.first;
  c  = 0;

  // Setup the input states, and initialise the output states
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { val                              = state.tracer0_rs[lon][lat][lev];
        field1.tracer[lon][lat][lev]     = val;
        field2.tracer[lon][lat][lev]     = Normal (val,           //in
                                                   fabs(val/5.0), //in
                                                   &rndseq);      //inout;
      }
    }
  }

  while (ob)
  { c++;
    // Adjoint test of ObservationOperator
    // ----------------------------------
    printf ("Adjoint test of ObservationOperator (separate observations)\n");
    // Setup the input states, and initialise the output states
    for (lon=0; lon<MetaData.nlon+2; lon++)
    { for (lat=0; lat<MetaData.nlat+2; lat++)
      { for (lev=0; lev<MetaData.nlev+2; lev++)
        { field1_hat.tracer[lon][lat][lev] = 0.0;
          field2_hat.tracer[lon][lat][lev] = 0.0;
        }
      }
    }

    // Set the adjoint of the state to zero
    for (lon=0; lon<MetaData.nlon+2; lon++)
    { for (lat=0; lat<MetaData.nlat+2; lat++)
      { for (t=0; t<MetaData.nss; t++)
        { state_hat.source_rs[lon][lat][t] = 0.0;
        }
      }
    }


    // Call the interpolation for this observation
    ObservationOperator (ob,            // inout Observation and data
                         &obs_trunk,    // in    Containing the mass profile (for satellites)
                         &MetaData,     // in    Lats, lons, etc.
                         &field1,       // in    Tracer field at lower time
                         &field2,       // in    Tracer field at upper time
                         0.0,           // in    Lower time of tracer field in s
                         43200.0,       // in    Upper time of tracer field in s
                         &state);       // in    Contains flux field

    // Transfer output to adjoint
    (*ob).dJo_dmodel_ob = (*ob).model_ob;

    // Call the adjoint
    ObservationOperator_adj (ob,            // in    Observation and data
                             &obs_trunk,    // in    Containing the mass profile (for satellites)
                             &MetaData,     // in    Lats, lons, etc.
                             &field1_hat,   // inout Tracer field at lower time
                             &field2_hat,   // inout Tracer field at upper time
                             &state_hat);   // inout Contains flux field

    lhs = (*ob).model_ob * (*ob).model_ob;
    rhs = 0.0;
    // The contribution from the tracer
    for (lon=0; lon<MetaData.nlon+2; lon++)
    { for (lat=0; lat<MetaData.nlat+2; lat++)
      { for (lev=0; lev<MetaData.nlev+2; lev++)
        { rhs += field1.tracer[lon][lat][lev] * field1_hat.tracer[lon][lat][lev] +
                 field2.tracer[lon][lat][lev] * field2_hat.tracer[lon][lat][lev];
        }
      }
    }
    // The contribution from the flux
    for (lon=0; lon<MetaData.nlon+2; lon++)
    { for (lat=0; lat<MetaData.nlat+2; lat++)
      { for (t=0; t<MetaData.nss; t++)
        { rhs += state.source_rs[lon][lat][t] * state_hat.source_rs[lon][lat][t];
        }
      }
    }

    printf ("Num: %i, ObOf %c, ObTpe %c\n", c, (*ob).ob_of, (*ob).obtpe);
    printf ("LHS = %e\n", lhs);
    printf ("RHS = %e\n", rhs);
    ob = (*ob).next;
  }

*/


  // Version 2 of the Adjoint test: treat all observations together
  // ===================================================================

  // Setup the input states, and initialise the output states
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { // Set the input tracer fields
        val                              = state.tracer0_rs[lon][lat][lev];
        field1.tracer[lon][lat][lev]     = val;
        field2.tracer[lon][lat][lev]     = Normal (val,           //in
                                                   fabs(val/5.0), //in
                                                   &rndseq);      //inout;
        // Set the adjoint variables
        field1_hat.tracer[lon][lat][lev] = 0.0;
        field2_hat.tracer[lon][lat][lev] = 0.0;
      }
    }
  }

  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (t=0; t<MetaData.nss; t++)
      { // The input flux fields have already been set (read-in)
        // The adjoint variables
        state_hat.source_rs[lon][lat][t] = 0.0;
      }
    }
  }

  printf ("Adjoint test of ObservationOperator\n");
  ob  = obs_trunk.first;
  c   = 0;
  lhs = 0;

  while (ob)
  { c++;

    //if ((*ob).ob_of == 't' && (*ob).obtpe != 's')
    //{
    // Call the interpolation for this observation
    ObservationOperator (ob,            // inout Observation and data
                         &obs_trunk,    // in    Containing the mass profile (for satellites)
                         &MetaData,     // in    Lats, lons, etc.
                         &field1,       // in    Tracer field at lower time
                         &field2,       // in    Tracer field at upper time
                         0.0,           // in    Lower time of tracer field in s
                         43200.0,       // in    Upper time of tracer field in s
                         &state);       // in    Contains flux field

    // Transfer output to adjoint
    (*ob).dJo_dmodel_ob = (*ob).model_ob;
    // printf ("%i %fl\n", c, (*ob).model_ob);

    // Call the adjoint
    ObservationOperator_adj (ob,            // in    Observation and data
                             &obs_trunk,    // in    Containing the mass profile (for satellites)
                             &MetaData,     // in    Lats, lons, etc.
                             &field1_hat,   // inout Tracer field at lower time
                             &field2_hat,   // inout Tracer field at upper time
                             &state_hat);   // inout Contains flux field
    lhs += (*ob).model_ob * (*ob).model_ob;
    //}
    ob   = (*ob).next;
  }

  printf ("There were %i observations\n", c);
  printf ("Computing the RHS\n");
  rhs = 0.0;
  // The contribution from the tracer
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (lev=0; lev<MetaData.nlev+2; lev++)
      { rhs += field1.tracer[lon][lat][lev] * field1_hat.tracer[lon][lat][lev] +
               field2.tracer[lon][lat][lev] * field2_hat.tracer[lon][lat][lev];
      }
    }
  }
  // The contribution from the flux
  for (lon=0; lon<MetaData.nlon+2; lon++)
  { for (lat=0; lat<MetaData.nlat+2; lat++)
    { for (t=0; t<MetaData.nss; t++)
      { rhs += state.source_rs[lon][lat][t] * state_hat.source_rs[lon][lat][t];
      }
    }
  }

  printf ("In standard notation\n");
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);
  printf ("In exponential notation\n");
  printf ("LHS = %e\n", lhs);
  printf ("RHS = %e\n", rhs);



  // Tidy up
  Deallocate_instant (&field1,
                      &MetaData);
  Deallocate_instant (&field2,
                      &MetaData);
  Deallocate_instant (&field1_hat,
                      &MetaData);
  Deallocate_instant (&field2_hat,
                      &MetaData);
  Deallocate_state (&state);
  Deallocate_state (&state_hat);
  Deallocate_metadata (&MetaData);
  delete[] obs_trunk.mass_profile;
  Destroy_Obs (&(obs_trunk.first));

}
