/* ==================================================================================
   Forward trajectory code
   Run forecast model

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_ForwardTraj.out

   To run (e.g.)

./Master_ForwardTraj.out <wind directory>                             \
                         <major timestep (between winds, seconds)>    \
                         <minor timestep (integration, seconds)>      \
                         <multiplication factor of vertical winds>    \
                         <run length (days)>                          \
                         <Output file>


   Modification history
   --------------------
   23/07/22 New Code. Ross Bannister

   Documentation
   -------------

   =============================================================================== */


   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )

{ char                 wind_dir[256];
  char                 wind_file[256];
  char                 output_file[256];
  double               Dt, dt, w_factor, runlength;
  struct CVTData_type  CVTdata;
  struct metadata_type MetaData;
  int                  nminor, gamma, nmajor;
  struct Wind_type     windA, windB;
  FILE*                output_unit = NULL;
  struct Wind_type     *wind_lowert, *wind_uppert, *wind_temp;
  int                  starttime, windnum, t;
  double               particle_lon, particle_lat, particle_lev;
  double               particle_lon_new, particle_lat_new, particle_lev_new;
  bool                 ok;
  char                 CVTfilename[256] = "../Calibration/CVT_calib.nc";

//   ===============================================================================
//   Initialisation
//   ===============================================================================

  // Process command line arguments
  ok = RunTraj_arguments (argument_count,
                          argument_list,
                          &Dt,
                          &dt,
                          &w_factor,
                          &runlength,
                          wind_dir,
                          output_file);


  if (ok)
  { // How many minor timesteps per major timestep?
    gamma = int(Dt/dt);
    printf ("There are %i minor timesteps per major timestep\n", gamma);

    // How many minor timesteps in the total integration?
    nminor = int(runlength * 24.0 * 3600.0 / dt);
    MetaData.ntimes_minor = nminor;
    printf ("There are %i minor timesteps in the total integration\n", nminor);

    // How many major time steps are involved in this integration?
    nmajor = int(runlength * 24.0 * 3600.0 / Dt) + 1;
    MetaData.ntimes_major = nmajor;
    if (nmajor < 2)
    { nmajor = 2;
    }
    printf ("There are %i wind fields required\n", nmajor);
    // (This may be different from the ::ntimes_major global variable)



    // The following are done just to get the metadata (longs, lats, levs, etc)
    // Initialise and read in generic data needed for the horizontal transform

    // Allocate the cvt
    printf ("Allocating CVT\n");
    Allocate_CVT (&CVTdata,
                  ::L, ::nlev, ::ntimes_major);
    printf ("Done\n");

    // Read-in cvt data
    printf ("Reading CVT data\n");
    cvt_matrices_input (&CVTdata,
                        CVTfilename);
    printf ("Read-in completed\n");

    // Set-up the meta data structure
    printf ("Allocating meta data\n");
    Allocate_metadata (&MetaData,
                       2*::L+1, ::L+1, ::nlev, nmajor);
    printf ("Done\n");

    printf ("Copying meta data\n");
    copy_metadata (CVTdata.nlon,         &(MetaData.nlon),
                   CVTdata.nlat,         &(MetaData.nlat),
                   ::nlev,               &(MetaData.nlev),
                   nmajor,               &(MetaData.ntimes_major),
                   CVTdata.L,            &(MetaData.L),
                   CVTdata.times,        MetaData.times,
                   CVTdata.longitude,    MetaData.longitude,
                   CVTdata.latitude,     MetaData.latitude,
                   CVTdata.level,        &(MetaData.level[CVTdata.nlev-::nlev]),
                   MetaData.cos_u_lat,   MetaData.cos_v_lat,
                   true);
    printf ("Done\n");


    // Free-up space not needed any more
    printf ("Freeing up space\n");
    Deallocate_CVT (&CVTdata);
    printf ("Done\n");

    printf ("Allocating windA\n");
    Allocate_wind (&windA,
                   &MetaData,
                   false);
    printf ("Allocating windB\n");
    Allocate_wind (&windB,
                   &MetaData,
                   false);

    printf ("Done allocations\n");


    printf ("Number of levels in MetaData : %i\n", MetaData.nlev);
    printf ("Number of levels in windA    : %i\n", windA.nlev);



//  ===============================================================================
//  Start to run the particle trajectories
//  ===============================================================================

    // Set up the output file for the fields
    printf ("Output file: %s\n", output_file);
    output_unit = fopen (output_file, "w");

    // Read-in the first two windfields
    windnum = 0;
    sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
    printf ("First wind file : %s\n", wind_file);
    Read_winds (&windA, &MetaData, w_factor, wind_file);
    windnum++;
    sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
    printf ("Second wind file: %s\n", wind_file);
    Read_winds (&windB, &MetaData, w_factor, wind_file);

    // Set-up the initial configuration of the winds
    wind_lowert = &windA;
    wind_uppert = &windB;

    // Set the initial position of the particle
    particle_lon = 90.0;
    particle_lat = 20.0;
    particle_lev = 3000.0;
    fprintf (output_unit, "%f  %f  %f  %f\n", 0.0, particle_lon, particle_lat, particle_lev);

    starttime = 0;
    for (t=starttime; t<nminor; t++)
    { printf ("Minor timestep %i\n", t);

      // Are we at a stage when we have just crossed a major timestep?
      // If so we need to adjust the pair of wind fields that straddle the current state
      if (((t-starttime) % gamma == 0) && (t > starttime))
      { // Yes
        printf ("Just crossed a major timestep\n");
        windnum++;
        // Swap pointers
        wind_temp   = wind_lowert;
        wind_lowert = wind_uppert;
        wind_uppert = wind_temp;

        // Read-in the next wind field
        sprintf (wind_file, "%s/Winds%04i.nc", wind_dir, windnum);
        printf ("Wind file: %s\n", wind_file);
        Read_winds (wind_uppert, &MetaData, w_factor, wind_file);
      }


      // Increment the position of the particle using fourth order Runge-Kutta
      trajectory ( t % gamma,         // Minor timestep No., local time (relative to wind_lowert field)
                   dt,                // Minor timestep size
                   Dt,                // Major timestep size
                   wind_lowert,       // Wind fields at start of major timestep
                   wind_uppert,       // Wind fields at end of major timestep
                   particle_lon,      // particle longitude at start time
                   particle_lat,      // particle latitude at start time
                   particle_lev,      // particle level at start time
                   &particle_lon_new, // particle longitude at end time
                   &particle_lat_new, // particle latitude at end time
                   &particle_lev_new, // particle level at end time
                   false );           // Don't show diagnostics

      particle_lon = particle_lon_new;
      particle_lat = particle_lat_new;
      particle_lev = particle_lev_new;

      fprintf (output_unit, "%f  %f  %f  %f\n", dt * double(t), particle_lon, particle_lat, particle_lev);
    }

    fclose (output_unit);



//   ===============================================================================
//   Deallocate
//   ===============================================================================

    Deallocate_wind (&windA);
    Deallocate_wind (&windB);
    Deallocate_metadata (&MetaData);
    printf ("Done program\n");
  }

}
