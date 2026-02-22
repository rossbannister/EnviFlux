/* ==================================================================================
   3d source and sink code
   Make a background state from a truth with a random perturbation

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_MakeBG.out

   Modification history
   --------------------
   05/01/24 Adapted from Master_AdjointTests_SL+ObsOp and Master_AdjointTests_CVT.
            Ross Bannister
   05/08/24 Add perturbation factors

   Documentation
   -------------
   ./Master_MakeBG.out \
     <Truth state (input)> \
     <Background (output)> \
     <Background perturbation (output)> \
     <CVT file (input)> \
     <std factor ics> \
     <std factor fluxes> \
     <pert factor ics> \
     <pert factor fluxes>


   =============================================================================== */


   #include <stdio.h>
   #include <math.h>
   #include <stdlib.h>
   #include <string.h>
   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )
{ char                       truth_in[256];
  char                       bg_out[256];
  char                       bgpert_out[256];
  char                       CVTfilename[256];
  struct metadata_type       MetaData;
  struct HorizTransData_type HorizData;
  struct state_type          truth, bg, bgpert, controlpert;
  struct CVTData_type        CVTdata;
  double                     factor_std_tracer, factor_std_flux;
  double                     factor_pert_tracer, factor_pert_flux;
  int                        rndseq = 9857395;   // Random number seed
  int                        x, y, t, lev, ok;


//   ===============================================================================
//   Initialisation
//   ===============================================================================


  // Process command line arguments
  ok = MakeBg_arguments (
         argument_count,     // in
         argument_list,      // in
         truth_in,           // out Filename of truth state (initial condition and source)
         bg_out,             // out Filename of background state (initial condition and source)
         bgpert_out,         // out Filename of background perturbation (initial condition and source)
         CVTfilename,        // out CVT filename
         &factor_std_tracer, // out To multiply tracer err std
         &factor_std_flux,   // out To multiply flux err std
         &factor_pert_tracer,// out To multiply tracer bg pert
         &factor_pert_flux); // out To multiply flux bg pert

  if (ok)
  { // Initialise and read in generic data needed for the horizontal transform
    printf ("Initialising horizontal transform\n");
    InitializeHorizTrans (&HorizData,
                          ::ALPfilename);

    // Read-in cvt data (including allocation of CVT)
    printf ("Reading CVT data\n");
    cvt_matrices_input (&CVTdata,
                        CVTfilename);
    printf ("Reading CVT completed\n");

    // Adjust the standard deviations
    printf ("Adjusting the background error standard deviations\n");
    for (lev=0; lev<CVTdata.nlev; lev++)
    { CVTdata.tracer_stddev[lev] *= factor_std_tracer;
    }
    for (x=0; x<CVTdata.nlon; x++)
    { for (y=0; y<CVTdata.nlat; y++)
      { for (t=0; t<CVTdata.nss; t++)
        { CVTdata.source_stddev[x][y][t] *= factor_std_flux;
        }
      }
    }
    printf ("Adjustment completed\n");

    // Set-up the meta data structure
    printf ("Allocating MetaData\n");
    Allocate_metadata (&MetaData,
                       2*::L+1, ::L+1, ::nlev, ::ntimes_major);
    printf ("Allocating MetaData completed\n");


    printf ("Copying MetaData\n");
    copy_metadata (CVTdata.nlon,         &(MetaData.nlon),
                   CVTdata.nlat,         &(MetaData.nlat),
                   CVTdata.nlev,         &(MetaData.nlev),
                   CVTdata.nss,          &(MetaData.nss),
                   CVTdata.L,            &(MetaData.L),
                   CVTdata.times,        MetaData.times,
                   CVTdata.longitude,    MetaData.longitude,
                   CVTdata.latitude,     MetaData.latitude,
                   CVTdata.level,        MetaData.level,
                   MetaData.cos_u_lat,   MetaData.cos_v_lat,
                   true);
    printf ("Copying MetaData completed\n");

  //   ===============================================================================
  //   Input the true state
  //   ===============================================================================

    printf ("Reading-in the true state vector\n");
    truth.horiz_repres = 'r';
    truth.vert_repres  = 'r';
    truth.temp_repres  = 'r';
    ReadStateVector (&truth,
                     &MetaData,
                     true, true,
                     truth_in);
    printf ("Reading-in the true state vector completed\n");


    // Allocate the background state (real space)
    printf ("Allocating background\n");
    Allocate_state (&bg,
                    &MetaData,
                    'r',
                    'r',
                    'r');
    printf ("Allocating background done\n");

    // Allocate the background perturbation (real space)
    printf ("Allocating background pert in real space\n");
    Allocate_state (&bgpert,
                    &MetaData,
                    'r',
                    'r',
                    'r');
    printf ("Allocating background pert in real space done\n");

    // Allocate the background perturbation (control space)
    printf ("Allocating background pert in control space\n");
    Allocate_state (&controlpert,
                    &MetaData,
                    's',
                    'm',
                    'm');
    printf ("Allocating background pert in control space done\n");


  //   ===============================================================================
  //   Generate a perturbation, add to truth and output
  //   ===============================================================================

    // Generate a random state perturbation
    printf ("Generating noise\n");
    PutWhiteNoiseControlVector (&controlpert,
                                &rndseq);
    printf ("Generating noise completed\n");

    // Transform perturbation to model space
    printf ("Transforming to real space\n");
    cvt_total (&controlpert,  // in  spectral space
               &bgpert,       // out real space
               &MetaData,     // in  metadata
               &HorizData,    // in  info about horiz transform
               &CVTdata);     // in  info about cvt
    printf ("Transforming to real space completed\n");

    // Multiply the generated background perturbation by a factor
    multiply_general (&bgpert,
                      &HorizData,
                      true,
                      factor_pert_tracer,
                      factor_pert_flux);

    // Add perturbation to the truth to give the background
    printf ("Adding to the truth\n");
    add_general ( &truth,
                  &bgpert,
                  &bg,
                  &HorizData,
                  true );
    printf ("Adding to the truth completed\n");

    // Output the background perturbation
    printf ("Outputting background pert\n");
    WriteStateVector (&bgpert,
                      bgpert_out);
    printf ("Outputting background pert completed\n");

    // Output the background
    printf ("Outputting background\n");
    WriteStateVector (&bg,
                      bg_out);
    printf ("Outputting background completed\n");


  //   ===============================================================================
  //   Tidy up
  //   ===============================================================================

    printf ("Tidying up\n");
    Deallocate_state (&truth);
    Deallocate_state (&bg);
    Deallocate_state (&bgpert);
    Deallocate_state (&controlpert);
    Deallocate_metadata (&MetaData);
    DeallocateHorizTrans (&HorizData);
    Deallocate_CVT (&CVTdata);
  }
  printf ("Program completed\n");
}
