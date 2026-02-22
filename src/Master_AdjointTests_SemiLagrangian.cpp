/* ==================================================================================
   3d source and sink code
   Adjoint tests for routines associated with the transport model

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_AdjointTests_SemiLagrangian.out

   Modification history
   --------------------
   18/12/21 New Code. Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <stdio.h>
   #include <math.h>
   #include <stdlib.h>
   #include <string.h>
   #include <source.h>


int main ()
{ double                     field1[2], field1_hat[2], axisx[2], intpointx, value;
  double                     field2[2][2], field2_hat[2][2], axisy[2], intpointy;
  double                     field3[2][2][2], field3_hat[2][2][2], axisz[2], intpointz;
  double                     field3a[2][2][2], field3a_hat[2][2][2], axist[2], intpointt;
  double                     field1cubic[6], axisx_cubic[6], field1cubic_hat[6];
  double                     field2cubic[6][6], axisy_cubic[6], field2cubic_hat[6][6];
  double                     field3cubic[6][6][6], axisz_cubic[6], field3cubic_hat[6][6][6];
  struct metadata_type       MetaData;
  struct instant_tracer_type tracer, tracer_int, tracer_hat;
  struct state_type          state, state_hat;
  double                     spacing;
  double                     lhs, rhs;
  int                        x, y, z, t;
  int                        nlon, nlat, nlev;
  int                        rndseq = 738395;   // Random number seed
  char                       inputstatefile[256] = "../data/State_56levs.nc";
  int                        index_lon, index_lat, index_lev;
  char                       wind_dir[256] = "../data/ECMWF_RedResWinds_1995_56levs";
  char                       wind_file[256];
  double                     Dt, dt, kappa_dt, kappa_h, kappa_v, gamma, sDt, factor_w;
  double                     fluxfactor, runlength;
  char                       interpolate_lc;
  char                       dp_file_rw;
  struct Wind_type           windA, windB;
  struct Wind_type           *wind_lowert, *wind_uppert, *wind_temp;
  struct instant_tracer_type Field1, Field2, Field2_hat, Field3;
  struct instant_tracer_type *field_lowert, *field_uppert, *field_temp;

  double                     ***dummyd;
  int                        ***dummyi;
  int                        nmajor, nminor, sfield, windnum, starttime, nss_needed, epsilon;
  char                       output_file_dps[256] = "AdjointTest_dps.nc";
  bool                       output_dps = true, inc_vert, output_anim, output_diags, inc_adv;
  char                       output_file_anim[256] = "ForwardIntegration.nc";
  char                       output_file_anim_adj[256] = "AdjointIntegration.nc";


/*
  // ===============================================================================
  // Adjoint test of Interpolate1D
  // ===============================================================================
  printf ("Adjoint test of Interpolate1D\n");
  axisx[0]  = 10.0;  axisx[1]  = 20.0;
  field1[0] = 35.0;  field1[1] = -5.0;
  intpointx = 11.5;

  // Forward routine
  value = Interpolate1D (field1, axisx, intpointx);
  // Adjoint routine
  Interpolate1D_adj (field1_hat, value, axisx, intpointx);
  // Compute LHS and RHS
  lhs = value * value;
  rhs = field1[0] * field1_hat[0] + field1[1] * field1_hat[1];
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);


  // ===============================================================================
  // Adjoint test of Interpolate2D
  // ===============================================================================
  printf ("Adjoint test of Interpolate2D\n");
  axisy[0]     = 1001.0;  axisy[1]     = 1002.0;
  field2[0][0] = 35.0;    field2[0][1] = -5.0;
  field2[1][0] = 17.0;    field2[1][1] = 2.0;
  intpointy = 1001.9;

  // Forward routine
  value = Interpolate2D (field2, axisx, axisy, intpointx, intpointy);
  // Adjoint routine
  Interpolate2D_adj (field2_hat, value, axisx, axisy, intpointx, intpointy);
  // Compute LHS and RHS
  lhs = value * value;
  rhs = 0.0;
  for (x=0; x<=1; x++)
  { for (y=0; y<=1; y++)
    { rhs += field2[x][y] * field2_hat[x][y];
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);


  // ===============================================================================
  // Adjoint test of Interpolate3D
  // ===============================================================================
  printf ("Adjoint test of Interpolate3D\n");
  axisz[0]     = 10000.0;  axisz[1]        = 11000.0;
  field3[0][0][0] = 35.0;  field3[0][1][0] = -5.0;
  field3[1][0][0] = 17.0;  field3[1][1][0] = 2.0;
  field3[0][0][1] = 37.0;  field3[0][1][1] = -8.0;
  field3[1][0][1] = 12.0;  field3[1][1][1] = -2.5;
  intpointz = 10040.9;

  // Forward routine
  value = Interpolate3D (field3, axisx, axisy, axisz, intpointx, intpointy, intpointz);
  // Adjoint routine
  Interpolate3D_adj (field3_hat, value, axisx, axisy, axisz, intpointx, intpointy, intpointz);
  // Compute LHS and RHS
  lhs = value * value;
  rhs = 0.0;
  for (x=0; x<=1; x++)
  { for (y=0; y<=1; y++)
    { for (z=0; z<=1; z++)
      { rhs += field3[x][y][z] * field3_hat[x][y][z];
      }
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);


  // ===============================================================================
  // Adjoint test of Interpolate3Dt
  // ===============================================================================
  printf ("Adjoint test of Interpolate3Dt\n");
  axist[0]         = 10.0;  axist[1]         = 20.0;
  field3a[0][0][0] = 33.0;  field3a[0][1][0] = -9.0;
  field3a[1][0][0] = 18.2;  field3a[1][1][0] = 2.8;
  field3a[0][0][1] = 37.2;  field3a[0][1][1] = -11.0;
  field3a[1][0][1] = 9.4;   field3a[1][1][1] = -2.8;
  intpointt = 19.0;

  // Forward routine
  value = Interpolate3Dt (field3, field3a, axisx, axisy, axisz, axist, intpointx, intpointy, intpointz, intpointt);
  // Adjoint routine
  Interpolate3Dt_adj (field3_hat, field3a_hat, value, axisx, axisy, axisz, axist, intpointx, intpointy, intpointz, intpointt);
  // Compute LHS and RHS
  lhs = value * value;
  rhs = 0.0;
  for (x=0; x<=1; x++)
  { for (y=0; y<=1; y++)
    { for (z=0; z<=1; z++)
      { rhs += field3[x][y][z] * field3_hat[x][y][z] + field3a[x][y][z] * field3a_hat[x][y][z];
      }
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);


  // ===============================================================================
  // Adjoint test of interpolate_1d_cubic (2 parts)
  // ===============================================================================
  printf ("Adjoint test of interpolate_1d_cubic (2 parts)\n");
  field1cubic[0] = -10.0;  field1cubic[1] = -8.0;  field1cubic[2] = 0.0;
  field1cubic[3] =   3.0;  field1cubic[4] = 12.0;  field1cubic[5] = 25.0;
  axisx_cubic[0] = 5.0;    axisx_cubic[1] = 10.0;  axisx_cubic[2] = 15.0;
  axisx_cubic[3] = 20.0;   axisx_cubic[4] = 25.0;  axisx_cubic[5] = 30.0;
  spacing = (axisx_cubic[3] - axisx_cubic[2]) / 20.0;
  for (intpointx=axisx_cubic[2]; intpointx<=axisx_cubic[3]; intpointx+=spacing)
  { // Forward routine
    value = interpolate_1d_cubic (axisx_cubic,  // grid positions
                                  field1cubic,   // field values (in)
                                  intpointx,     // position to interpolate to
                                  2 );           // divide domain into 2 or 3 parts
    // Adjoint routine
    interpolate_1d_cubic_adj (axisx_cubic,      // grid positions
                              field1cubic_hat,   // field values (out)
                              value,             // interpolated value (in)
                              intpointx,         // position to interpolate to
                              2 );               // divide domain into 2 or 3 parts

    // Compute LHS and RHS
    lhs = value * value;
    rhs = 0.0;
    for (x=0; x<=5; x++)
    { rhs += field1cubic[x] * field1cubic_hat[x];
    }
    printf ("intpointx = %f\n", intpointx);
    printf ("  LHS = %f\n", lhs);
    printf ("  RHS = %f\n", rhs);
  }

  // ===============================================================================
  // Adjoint test of interpolate_1d_cubic (3 parts)
  // ===============================================================================
  printf ("Adjoint test of interpolate_1d_cubic (3 parts)\n");
  field1cubic[0] = -10.0;  field1cubic[1] = -8.0;  field1cubic[2] = 0.0;
  field1cubic[3] =   3.0;  field1cubic[4] = 12.0;  field1cubic[5] = 25.0;
  axisx_cubic[0] = 5.0;    axisx_cubic[1] = 10.0;  axisx_cubic[2] = 15.0;
  axisx_cubic[3] = 20.0;   axisx_cubic[4] = 25.0;  axisx_cubic[5] = 30.0;
  spacing = (axisx_cubic[3] - axisx_cubic[2]) / 20.0;

  for (intpointx=axisx_cubic[2]; intpointx<=axisx_cubic[3]; intpointx+=spacing)
  { // Forward routine
    value = interpolate_1d_cubic (axisx_cubic,  // grid positions
                                  field1cubic,   // field values (in)
                                  intpointx,     // position to interpolate to
                                  3 );           // divide domain into 2 or 3 parts
    // Adjoint routine
    interpolate_1d_cubic_adj (axisx_cubic,      // grid positions
                              field1cubic_hat,   // field values (out)
                              value,             // interpolated value (in)
                              intpointx,         // position to interpolate to
                              3 );               // divide domain into 2 or 3 parts

    // Compute LHS and RHS
    lhs = value * value;
    rhs = 0.0;
    for (x=0; x<=5; x++)
    { rhs += field1cubic[x] * field1cubic_hat[x];
    }
    printf ("intpointx = %f\n", intpointx);
    printf ("  LHS = %f\n", lhs);
    printf ("  RHS = %f\n", rhs);
  }


  // ===============================================================================
  // Adjoint test of interpolate_2d_cubic (3 parts)
  // ===============================================================================
  printf ("Adjoint test of interpolate_2d_cubic (3 parts)\n");

  axisy_cubic[0] = 25.0;  axisy_cubic[1] = 10.0;  axisy_cubic[2] = 5.0;
  axisy_cubic[3] = 3.0;   axisy_cubic[4] = -1.0;  axisy_cubic[5] = -10.0;
  // Fill data array
  for (x=0; x<6; x++)
  { for (y=0; y<6; y++)
    { field2cubic[x][y] = 10.0 * cos(double(2*x)) * log(double(y+3)) + randomno(&rndseq);
    }
  }

  // Forward routine
  value = interpolate_2d_cubic (axisx_cubic,  // x grid positions
                                axisy_cubic,  // y grid positions
                                field2cubic,  // field values (in)
                                16.5,         // x position to interpolate to
                                4.2,          // x position to interpolate to
                                3 );          // divide domain into 2 or 3 parts
  // Adjoint routine
  interpolate_2d_cubic_adj (axisx_cubic,      // x grid positions
                            axisy_cubic,      // y grid positions
                            field2cubic_hat,  // field values (out)
                            value,            // interpolated value (in)
                            16.5,             // x position to interpolate to
                            4.2,              // x position to interpolate to
                            3 );              // divide domain into 2 or 3 parts

  // Compute LHS and RHS
  lhs = value * value;
  rhs = 0.0;
  for (x=0; x<=5; x++)
  { for (y=0; y<=5; y++)
    { rhs += field2cubic[x][y] * field2cubic_hat[x][y];
    }
  }
  printf ("  LHS = %f\n", lhs);
  printf ("  RHS = %f\n", rhs);


  // ===============================================================================
  printf ("Adjoint test of interpolate_3d_cubic (3 parts)\n");

  axisz_cubic[0] = 1000.0;  axisz_cubic[1] = 1100.0;  axisz_cubic[2] = 1250.0;
  axisz_cubic[3] = 1450.0;  axisz_cubic[4] = 1800.0;  axisz_cubic[5] = 2400.0;
  // Fill data array
  for (x=0; x<6; x++)
  { for (y=0; y<6; y++)
    { for (z=0; z<6; z++)
      { field3cubic[x][y][z] = 100.0 * cos(double(2*x)) * log(double(y+3)) * sin(double(z)/5.0 + 0.2) + randomno(&rndseq);
      }
    }
  }

  // Forward routine
  value = interpolate_3d_cubic (axisx_cubic,  // x grid positions
                                axisy_cubic,  // y grid positions
                                axisz_cubic,  // y grid positions
                                field3cubic,  // field values (in)
                                16.5,         // x position to interpolate to
                                4.2,          // x position to interpolate to
                                1255.0,       // z position to interpolate to
                                3 );          // divide domain into 2 or 3 parts
  // Adjoint routine
  interpolate_3d_cubic_adj (axisx_cubic,      // x grid positions
                            axisy_cubic,      // y grid positions
                            axisz_cubic,      // y grid positions
                            field3cubic_hat,  // field values (out)
                            value,            // interpolated value (in)
                            16.5,             // x position to interpolate to
                            4.2,              // x position to interpolate to
                            1255.0,           // z position to interpolate to
                            3 );              // divide domain into 2 or 3 parts

  // Compute LHS and RHS
  lhs = value * value;
  rhs = 0.0;
  for (x=0; x<=5; x++)
  { for (y=0; y<=5; y++)
    { for (z=0; z<=5; z++)
      { rhs += field3cubic[x][y][z] * field3cubic_hat[x][y][z];
      }
    }
  }
  printf ("  LHS = %f\n", lhs);
  printf ("  RHS = %f\n", rhs);

*/

  // ===============================================================================
  // ROUTINES THAT REQUIRE FULL FIELDS
  state.horiz_repres = 'r';
  state.vert_repres  = 'r';
  state.temp_repres  = 'r';
  printf ("Reading a state vector in %s\n", inputstatefile);
  ReadStateVector (&state,
                   &MetaData,
                   true, true,   // Read tracer and source
                   inputstatefile);
  printf ("Done\n");
  nlon = MetaData.nlon;
  nlat = MetaData.nlat;
  nlev = MetaData.nlev;
  printf ("Dimensions (nlon, nlat, nlev: %i, %i, %i)\n", nlon, nlat, nlev);


/*

  // ===============================================================================
  printf ("Adjoint test of halos_tracer\n");

  Allocate_instant (&tracer,
                    &MetaData);
  Allocate_instant (&tracer_int,
                    &MetaData);
  Allocate_instant (&tracer_hat,
                    &MetaData);

  // Put random numbers inside the bulk
  for (x=1; x<=nlon; x++)
  { for (y=1; y<=nlat; y++)
    { for (z=1; z<=nlev; z++)
      { tracer.tracer[x][y][z] = 200.0 * (randomno(&rndseq) - 0.5);
      }
    }
  }

  // Forward routine (since this operator over-writes, take copy first)
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { tracer_int.tracer[x][y][z] = tracer.tracer[x][y][z];
      }
    }
  }

  halos_tracer (&tracer_int,
                &MetaData);

  // Adjoint routine (since this operator over-writes, take copy first)
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { tracer_hat.tracer[x][y][z] = tracer_int.tracer[x][y][z];
      }
    }
  }

  halos_tracer_adj (&tracer_hat,
                    &MetaData);

  // Compute LHS and RHS
  lhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { lhs += tracer_int.tracer[x][y][z] * tracer_int.tracer[x][y][z];
      }
    }
  }
  rhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { rhs += tracer.tracer[x][y][z] * tracer_hat.tracer[x][y][z];
      }
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);


  Deallocate_instant (&tracer,
                      &MetaData);
  Deallocate_instant (&tracer_int,
                      &MetaData);
  Deallocate_instant (&tracer_hat,
                      &MetaData);


  // ===============================================================================
  printf ("Adjoint test of Fill_short_array and Fill_short_array_cubic\n");
  // This is used with linear interpolation

  Allocate_instant (&tracer,
                    &MetaData);
  Allocate_instant (&tracer_hat,
                    &MetaData);

  // Put random numbers inside the bulk ...
  for (x=1; x<=nlon; x++)
  { for (y=1; y<=nlat; y++)
    { for (z=1; z<=nlev; z++)
      { tracer.tracer[x][y][z] = 200.0 * (randomno(&rndseq) - 0.5);
      }
    }
  }
  // ... and set halos to make field consistent
  halos_tracer (&tracer,
                &MetaData);

  // Choose a point in the domain
  index_lon = 1;
  index_lat = 32;
  index_lev = 55;


  printf ("Adjoint test of Fill_short_array\n");
  // Forward routine - prepare cube of data
  Fill_short_array (tracer.tracer,
                    field3,
                    index_lon, index_lat, index_lev);

  // Adjoint routine
  Fill_short_array_adj (tracer_hat.tracer,
                        field3,
                        index_lon, index_lat, index_lev);

  // Compute LHS and RHS
  lhs = 0.0;
  for (x=0; x<=1; x++)
  { for (y=0; y<=1; y++)
    { for (z=0; z<=1; z++)
      { lhs += field3[x][y][z] * field3[x][y][z];
      }
    }
  }
  rhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { rhs += tracer.tracer[x][y][z] * tracer_hat.tracer[x][y][z];
      }
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  printf ("Adjoint test of Fill_short_array_cubic\n");

  // Forward routine
  Fill_short_array_cubic (nlon, nlat, nlev,
                          tracer.tracer,
                          field3cubic,
                          index_lon, index_lat, index_lev);

  // Set tracer_hat to initial value of zero
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { tracer_hat.tracer[x][y][z] = 0.0;
      }
    }
  }

  // Adjoint routine
  Fill_short_array_cubic_adj (nlon, nlat, nlev,
                              tracer_hat.tracer,
                              field3cubic,
                              index_lon, index_lat, index_lev);

  // Compute LHS and RHS
  lhs = 0.0;
  for (x=0; x<=5; x++)
  { for (y=0; y<=5; y++)
    { for (z=0; z<=5; z++)
      { lhs += field3cubic[x][y][z] * field3cubic[x][y][z];
      }
    }
  }
  rhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { rhs += tracer.tracer[x][y][z] * tracer_hat.tracer[x][y][z];
      }
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  Deallocate_instant (&tracer,
                      &MetaData);
  Deallocate_instant (&tracer_hat,
                      &MetaData);

*/

/*

  // ===============================================================================
  // One step of the semi-Lagrangian model
  // ===============================================================================
  printf ("Adjoint test of SemiLagrangian_1step\n");

  // The state vector has already been read-in (see above)
  Dt             = 43200.0;
  dt             = 3600.0;
  kappa_dt       = 86.4;
  kappa_h        = 4000000.0 * 0.0;
  kappa_v        = 2.0 * 0.0;
  interpolate_lc = 'c';
  nmajor         = 1;
  factor_w       = 1.0;
  output_dps     = true;
  inc_adv        = true;
  inc_vert       = true;
  output_anim    = true;
  dp_file_rw     = 'r';    // Read or write departure points file?
                           // Write means calculate departure points in forward SL scheme and output
                           // Read means read-in previously calculated departure points

  // How many minor timesteps per major timestep?
  gamma = double(int(Dt/dt));

  // Determine the time step between source/sink fields
  sDt = state.times[1] - state.times[0];
  printf ("The time between surface flux fields in %f s (%f days)\n", sDt, sDt / (24.0 * 3600.0));

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

  // Allocate space for data
  printf ("Allocating windA\n");
  Allocate_wind (&windA,
                 &MetaData,
                 false);
  printf ("Allocating windB\n");
  Allocate_wind (&windB,
                 &MetaData,
                 false);

  printf ("Allocating Field1\n");
  Allocate_instant (&Field1,
                    &MetaData);
  printf ("Allocating Field2\n");
  Allocate_instant (&Field2,
                    &MetaData);

  printf ("Done allocations\n");

  printf ("Number of levels in MetaData : %i\n", MetaData.nlev);
  printf ("Number of levels in windA    : %i\n", windA.nlev);


  // Set the initial conditions
  for (x=0; x<state.nlon+2; x++)
  { for (y=0; y<state.nlat+2; y++)
    { for (z=0; z<state.nlev+2; z++)
      { Field1.tracer[x][y][z] = state.tracer0_rs[x][y][z];
      }
    }
  }
  Field1.timestep_major = 0;
  Field1.timestep_minor = 0;
  Field1.time           = 0.0;


  // Set up the output file for the forward integration fields
  if (output_anim)
  { printf ("Output file for forward animation: %s\n", output_file_anim);
    WriteTimeSeq (&MetaData,
                  &Field1,
                  nmajor+1,
                  Dt,
                  0,
                  0,       // 0=create file only
                  output_file_anim);

    // Write out t=0 field
    WriteTimeSeq (&MetaData,
                  &Field1,
                  nmajor+1,
                  Dt,
                  0,
                  1,       // 1=normal (write only)
                  output_file_anim);
  }



  // Set up the file for the departure points
  if (dp_file_rw == 'w')
  { printf ("Output file for departure points: %s\n", output_file_dps);
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
                          output_file_dps);
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

  starttime = 1;
  sfield    = int((double(starttime) * Dt) / sDt);

  // Run a model major timestep
  SemiLagrangian_1step ( &MetaData,         // Meta data
                         &state,            // State vector
                         &Field1,           // Tracer field at major t-1
                         &Field2,           // Tracer field extrapolated to major t
                         starttime,         // Major timestep number (new time)
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
                         &windA,            // Wind fields at lower major timestep
                         &windB,            // Wind fields at next major step
                         interpolate_lc,    // Interpolation type
                         output_dps,        // Output departure points?
                         output_file_dps,   // Filename of departure points file
                         dp_file_rw );      // Whether to (r)ead or (w)rite departure points file

  if (output_anim)
  { // Write out field
    WriteTimeSeq (&MetaData,
                  &Field2,
                  nmajor+1,
                  Dt,
                  1,
                  1,       // 1=normal (write only)
                  output_file_anim);
  }

  // Deallocate space
  Deallocate_wind (&windA);
  Deallocate_wind (&windB);

  printf ("Allocating Field3\n");
  Allocate_instant (&Field3,
                    &MetaData);

  printf ("Allocating state_hat\n");
  Allocate_state (&state_hat,
                  &MetaData,
                  'r',
                  'r',
                  'r');


  printf ("Allocating Field2_hat\n");
  Allocate_instant (&Field2_hat,
                    &MetaData);

  // Take a copy of Field2 (since it gets over-written inside SemiLagrangian_1step_adj)
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { Field2_hat.tracer[x][y][z] = Field2.tracer[x][y][z];
      }
    }
  }

  // Set up the output file for the adjoint integration fields
  if (output_anim)
  { printf ("Output file for adjoint animation: %s\n", output_file_anim_adj);
    WriteTimeSeq (&MetaData,
                  &Field2_hat,
                  nmajor+1,
                  Dt,
                  0,
                  0,       // 0=create file only
                  output_file_anim_adj);

    // Write out t=starttime field
    WriteTimeSeq (&MetaData,
                  &Field2_hat,
                  nmajor+1,
                  Dt,
                  1,
                  1,       // 1=normal (write only)
                  output_file_anim_adj);
  }


  // Run adjoint of a model timestep
  SemiLagrangian_1step_adj ( &MetaData,         // Meta data
                             &state_hat,        // In/Out: State vector
                             &Field3,           // Out: Tracer field at t-1
                             &Field2_hat,       // In: Tracer field extrapolated to t
                             starttime,         // Major timestep number
                             sfield,            // source/sink field number
                             Dt,                // Major timestep size
                             kappa_dt,          // Timestep for diffusion
                             kappa_h,           // Horiz diffusion coefficient
                             kappa_v,           // Vert diffusion coefficient
                             inc_adv,           // Include advection
                             fluxfactor,        // Flux factor
                             inc_vert,          // Include vertical transport
                             interpolate_lc,    // Linear or cubic interpolation
                             output_file_dps ); // Departure points input file

  if (output_anim)
  { // Write out field
    WriteTimeSeq (&MetaData,
                  &Field3,
                  nmajor+1,
                  Dt,
                  0,
                  1,       // 1=normal (write only)
                  output_file_anim_adj);
  }


  // Compute LHS and RHS
  lhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { lhs += Field1.tracer[x][y][z] * Field3.tracer[x][y][z];
      }
    }
  }
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { lhs += state.source_rs[x][y][sfield] * state_hat.source_rs[x][y][sfield];
    }
  }

  rhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { rhs += Field2.tracer[x][y][z] * Field2.tracer[x][y][z];
      }
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);




  // Deallocate (including metadata)
  Deallocate_instant (&Field1,
                      &MetaData);
  Deallocate_instant (&Field2,
                      &MetaData);
  Deallocate_instant (&Field3,
                      &MetaData);
  Deallocate_instant (&Field2_hat,
                      &MetaData);

  Deallocate_metadata (&MetaData);
  Deallocate_state (&state);
  Deallocate_state (&state_hat);


*/




  // ===============================================================================
  // Many steps of the semi-Lagrangian model
  // ===============================================================================
  printf ("Adjoint test of multi-step SemiLagrangian integration\n");


  // The state vector has already been read-in (see above)
  Dt             = 43200.0;
  dt             = 3600.0;
  kappa_dt       = 86.4;
  kappa_h        = 4000000.0 * 0.0;
  kappa_v        = 2.0 * 0.0;
  interpolate_lc = 'c';
  runlength      = 25.0 * 24.0 * 3600.0;
  factor_w       = 1.0;
  output_dps     = true;
  inc_adv        = true;
  inc_vert       = true;
  output_anim    = false;
  output_diags   = false;
  dp_file_rw     = 'r';    // Read or write departure points file?
                           // Write means calculate departure points in forward SL scheme and output
                           // Read means read-in previously calculated departure points

  // How many minor timesteps per major timestep?
  gamma = double(int(Dt/dt));
  printf ("There are %f minor timesteps per major timestep\n", gamma);

  // How many minor timesteps in the total integration?
  nminor = int(runlength / dt);
  MetaData.ntimes_minor = nminor;
  printf ("There are %i minor timesteps in the total integration\n", nminor);

  // How many major time steps are involved in this integration?
  nmajor = int(runlength / Dt);
  MetaData.ntimes_major = nmajor;
  printf ("There are %i major timesteps in this integration\n", nmajor);
  printf ("There are %i wind fields required\n", nmajor+1);

  // Read-in the state (containing the initial conditions and the surface fluxes)
//  printf ("Reading-in the state vector (initial conditions and surface flux time sequence)\n");
//  state.horiz_repres = 'r';
//  state.vert_repres  = 'r';
//  state.temp_repres  = 'r';
//  ReadStateVector (&state,
//                   &MetaData,
//                   true, true,
//                   inputstatefile);
//  printf ("Reading-in done\n");

  // Determine the time step between source/sink fields
  sDt = state.times[1] - state.times[0];
  printf ("The time between surface flux fields is %f s (%f days)\n", sDt, sDt / (24.0 * 3600.0));

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
  Allocate_instant (&Field1,
                    &MetaData);
  printf ("Allocating field2\n");
  Allocate_instant (&Field2,
                    &MetaData);

  printf ("Done allocations\n");

  nlon = MetaData.nlon;
  nlat = MetaData.nlat;
  nlev = MetaData.nlev;
  printf ("Dimensions (nlon, nlat, nlev: %i, %i, %i)\n", nlon, nlat, nlev);


  printf ("Number of levels in MetaData : %i\n", MetaData.nlev);
  printf ("Number of levels in windA    : %i\n", windA.nlev);

  if (output_anim)  { printf ("Will output animation file\n"); }
  if (output_diags) { printf ("Will output diagnostics file\n"); }
  if (output_dps)   { printf ("Will output departure points file\n"); }

  // Set the initial conditions
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { Field1.tracer[x][y][z] = state.tracer0_rs[x][y][z];
      }
    }
  }
  Field1.timestep_major = 0;
  Field1.timestep_minor = 0;
  Field1.time           = 0.0;



  // Set up the output file for the forward integration fields
  if (output_anim)
  { printf ("Output file for forward animation: %s\n", output_file_anim);
    WriteTimeSeq (&MetaData,
                  &Field1,
                  nmajor+1,
                  Dt,
                  0,
                  0,       // 0=create file only
                  output_file_anim);

    // Write out t=0 field
    WriteTimeSeq (&MetaData,
                  &Field1,
                  nmajor+1,
                  Dt,
                  0,
                  1,       // 1=normal (write only)
                  output_file_anim);
  }


  // Set up the file for the departure points
  if (dp_file_rw == 'w')
  { printf ("Output file for departure points: %s\n", output_file_dps);
    printf ("Setting-up file for writing departure points\n");
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
                          output_file_dps);
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
  field_lowert = &Field1;
  field_uppert = &Field2;
  starttime    = 1;

  // ===============================================================================
  // Run the transport model
  // ===============================================================================

  for (t=starttime; t<=nmajor; t++)
  { sfield = int((fabs(double(t)-0.1) * Dt) / sDt);
    printf ("Major timestep %i of %i, source field %i\n", t, nmajor, sfield);
    // Run a model major timestep
    SemiLagrangian_1step ( &MetaData,         // Meta data
                           &state,            // in:  State vector
                           field_lowert,      // in:  Tracer field at major t-1
                           field_uppert,      // out: Tracer field extrapolated to major t
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
                           output_dps,        // Output departure points?
                           output_file_dps,   // Filename of departure points file
                           dp_file_rw );      // Whether to (r)ead or (w)rite departure points file

    if (output_anim)
    { // Write out field
      WriteTimeSeq (&MetaData,
                    field_uppert,
                    nmajor+1,
                    Dt,
                    t,
                    1,       // 1=normal (write only)
                    output_file_anim);
    }



    if (t < nmajor)
    { // Bother with next winds fields only if there is still another step to do
      printf ("Just crossed a major timestep (forward)\n");
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


  // Take a copy of the state at the end (since it gets over-written inside SemiLagrangian_1step_adj)
  printf ("Allocating Field2_hat\n");
  Allocate_instant (&Field2_hat,
                    &MetaData);

  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { Field2_hat.tracer[x][y][z] = (*field_uppert).tracer[x][y][z];
      }
    }
  }

  // Set up the output file for the adjoint integration fields
  if (output_anim)
  { printf ("Output file for adjoint animation: %s\n", output_file_anim_adj);
    WriteTimeSeq (&MetaData,
                  &Field2_hat,
                  nmajor+1,
                  Dt,
                  0,
                  0,       // 0=create file only
                  output_file_anim_adj);

    // Write out t=starttime field
    WriteTimeSeq (&MetaData,
                  field_uppert,
                  nmajor+1,
                  Dt,
                  nmajor,
                  1,       // 1=normal (write only)
                  output_file_anim_adj);
  }



  // Set-up new state adjoint variable
  printf ("Allocating state_hat\n");
  Allocate_state (&state_hat,
                  &MetaData,
                  'r',
                  'r',
                  'r');

  // ===============================================================================
  // Run the adjoint transport model
  // ===============================================================================

  for (t=nmajor; t>=starttime; t--)
  { sfield = int((fabs(double(t)-0.1) * Dt) / sDt);
    printf ("Major adjoint timestep %i of %i, source field %i\n", t, nmajor, sfield);
    // Run a model major timestep

    // Set field_lowert to zero
    for (x=0; x<=nlon+1; x++)
    { for (y=0; y<=nlat+1; y++)
      { for (z=0; z<=nlev+1; z++)
        { (*field_lowert).tracer[x][y][z] = 0.0;
        }
      }
    }

    SemiLagrangian_1step_adj ( &MetaData,         // Meta data
                               &state_hat,        // inout: State vector
                               field_lowert,      // out:   Tracer field at major t-1
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
                               output_file_dps ); // Filename of departure points file

    if (output_anim)
    { // Write out field
      WriteTimeSeq (&MetaData,
                    field_lowert,
                    1,
                    Dt,
                    t-1,
                    1,       // 1=normal (write only)
                    output_file_anim_adj);
    }


    if (t > starttime)
    { // Bother with next winds fields only if there is still another step to do
      printf ("Just crossed a major timestep (adjoint)\n");
      // Sort out the configuration of the instantateous fields for the next step
      field_temp   = field_lowert;
      field_lowert = field_uppert;
      field_uppert = field_temp;
    }
  }

  // Set the initial conditions
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { state_hat.tracer0_rs[x][y][z] = (*field_lowert).tracer[x][y][z];
      }
    }
  }


  // Compute LHS and RHS
  lhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { lhs += state.tracer0_rs[x][y][z] * state_hat.tracer0_rs[x][y][z];
      }
    }
  }
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (t=0; t<nss_needed; t++)
      { if ((x==0) && (y==0))
        { printf ("t = %i\n", t);
        }
        lhs += state.source_rs[x][y][t] * state_hat.source_rs[x][y][t];
      }
    }
  }

  rhs = 0.0;
  for (x=0; x<=nlon+1; x++)
  { for (y=0; y<=nlat+1; y++)
    { for (z=0; z<=nlev+1; z++)
      { rhs += Field2_hat.tracer[x][y][z] * Field2_hat.tracer[x][y][z];
      }
    }
  }
  printf ("LHS = %f\n", lhs);
  printf ("RHS = %f\n", rhs);

  Deallocate_state (&state);
  Deallocate_state (&state_hat);
  Deallocate_wind (&windA);
  Deallocate_wind (&windB);
  Deallocate_instant (&Field1,
                      &MetaData);
  Deallocate_instant (&Field2,
                      &MetaData);
  Deallocate_instant (&Field2_hat,
                      &MetaData);
  Deallocate_metadata (&MetaData);


}
