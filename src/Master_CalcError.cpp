/* ==================================================================================
   3d source and sink code
   Calculate the error between specified states
   State - Truth, State - Background
   RMS(State - Truth), RMS(State - Background)
   Calculated separately for initial tracer and surface flux
   Done for norm specified in chosen B-matrix (coded as cov file)

   Ross Bannister, Data Assimiltion Research Centre,
   University of Reading, UK

   Language: C++

   Ubuntu
   make Master_CalcError.out

   Modification history
   --------------------
   31/03/24 New Code. Ross Bannister
   01/05/25 Allow output of data to file.  Ross Bannister
   30/09/25 Normalise errors by stanard deviations.  Ross Bannister
   10/10/25 Compute difference from background as well as from truth.  Ross Bannister
   02/12/25 Formalise to calculate values with mass weighting.  Ross Bannister

   Documentation
   -------------
   Run with
   ./Master_CalcError.out <State> <Truth> <Background> <CVT file> <output file>



   =============================================================================== */

   #include <source.h>


int main ( int   argument_count,
           char  **argument_list )
{ char                       State_file[256];
  char                       Truth_file[256];
  char                       Background_file[256];
  char                       CVTfilename[256];
  char                       Output_file[256];
  struct state_type          State, Truth, Background, Diff;
  struct metadata_type       MetaData;
  struct HorizTransData_type HorizData;
  struct CVTData_type        CVTdata;
  int                        ok, lev;
  double                     *dz, *density;
  double                     diff_tracer, diff_flux, rse_diff_tracer, rse_diff_flux;
  FILE*                      Diags;

  // ===============================================================================
  // Initialisation
  // ===============================================================================

  ok = CalcError_arguments ( argument_count,
                             argument_list,
                             State_file,
                             Truth_file,
                             Background_file,
                             CVTfilename,
                             Output_file );

  if (ok)
  { // Initialise and read in generic data needed for the horizontal transform
    printf ("Initialising horizontal transform\n");
    InitializeHorizTrans (&HorizData,
                          ::ALPfilename);

    // Read-in cvt data (including allocation of CVT)
    printf ("Reading CVT data\n");
    cvt_matrices_input (&CVTdata,
                        CVTfilename);
    printf ("Read-in completed\n");

    // Set-up the meta data structure
    Allocate_metadata (&MetaData,
                       2*CVTdata.L+1, CVTdata.L+1, CVTdata.nlev, CVTdata.nss);

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
                   false);


/*
    for (lev=0; lev<MetaData.nlev; lev++)
    { printf ("%i %f\n", lev+1, CVTdata.tracer_stddev[lev]);
    }
    for (t=0; t<MetaData.nss; t++)
    { printf ("Time, %i\n", t);
      for (lat=0; lat<MetaData.nlat; lat++)
      { for (lon=0; lon<MetaData.nlon; lon++)
        { printf ("%f ", CVTdata.source_stddev[lon][lat][t]);
        }
        printf ("\n");
      }
    }
*/



    // Read-in the states
    printf ("Reading State from file\n");
    State.horiz_repres = 'r';
    State.vert_repres  = 'r';
    State.temp_repres  = 'r';
    ReadStateVector (&State,
                     &MetaData,
                     true,
                     true,
                     State_file);
    printf ("Done\n");

    printf ("Reading truth from file\n");
    Truth.horiz_repres = 'r';
    Truth.vert_repres  = 'r';
    Truth.temp_repres  = 'r';
    ReadStateVector (&Truth,
                     &MetaData,
                     true,
                     true,
                     Truth_file);
    printf ("Done\n");

    printf ("Reading background from file\n");
    Background.horiz_repres = 'r';
    Background.vert_repres  = 'r';
    Background.temp_repres  = 'r';
    ReadStateVector (&Background,
                     &MetaData,
                     true,
                     true,
                     Background_file);
    printf ("Done\n");

    // Allocate the difference state
    Allocate_state (&Diff,
                    &MetaData,
                    'r',
                    'r',
                    'r');


    // ===== Set-up vertical profiles of level thicknesses and density
    Array_1d_double_create (&dz, MetaData.nlev);
    Array_1d_double_create (&density, MetaData.nlev);
    for (lev=0; lev<Truth.nlev; lev++)
    { dz[lev]      = (MetaData.level[lev+2] - MetaData.level[lev]) / 2.0;
      density[lev] = ::rho0 * exp(-1.0 * MetaData.level[lev+1] / ::H);
      //printf ("%i  %f %f\n", lev, dz[lev], density[lev]);
    }

    // ===== Set-up diagnostics file for output
    Diags = fopen (Output_file, "w");
    fprintf (Diags, "# ==============================================\n");
    fprintf (Diags, "# TRACER (Tg)                FLUX (Tg/month)\n");
    fprintf (Diags, "# Mean diff  RMS diff        Mean diff  RMS diff\n");
    fprintf (Diags, "# ==============================================\n");

    // ===== Take the difference between the background and the truth =====
    add_general_fac (&Background,  // in   input state 1
                     &Truth,       // in   input state 2
                     &Diff,        // out  state1 - state2
                     &HorizData,   // in   not used
                     false,        // in   include halos
                     -1.0 );
    CalcErr (&MetaData,
             &Diff,
             dz,
             density,
             &diff_tracer,
             &rse_diff_tracer,
             &diff_flux,
             &rse_diff_flux);

    fprintf (Diags, "# ==========  Background - Truth     ============\n");
    fprintf (Diags, "%f %f        %f %f\n", diff_tracer, rse_diff_tracer, diff_flux, rse_diff_flux);




    // ===== Take the difference between the state and the truth =====
    add_general_fac (&State,       // in   input state 1
                     &Truth,       // in   input state 2
                     &Diff,        // out  state1 - state2
                     &HorizData,   // in   not used
                     false,        // in   include halos
                     -1.0 );

    CalcErr (&MetaData,
             &Diff,
             dz,
             density,
             &diff_tracer,
             &rse_diff_tracer,
             &diff_flux,
             &rse_diff_flux);

    fprintf (Diags, "# ==========  Analysis - Truth       ============\n");
    fprintf (Diags, "%f %f        %f %f\n", diff_tracer, rse_diff_tracer, diff_flux, rse_diff_flux);




    // ===== Take the difference between the state and the background =====
    add_general_fac (&State,       // in   input state 1
                     &Background,  // in   input state 2
                     &Diff,        // out  state1 - state2
                     &HorizData,   // in   not used
                     false,        // in   include halos
                     -1.0 );

    CalcErr (&MetaData,
             &Diff,
             dz,
             density,
             &diff_tracer,
             &rse_diff_tracer,
             &diff_flux,
             &rse_diff_flux);

    fprintf (Diags, "# ==========  Analysis - Background  ============\n");
    fprintf (Diags, "%f %f        %f %f\n", diff_tracer, rse_diff_tracer, diff_flux, rse_diff_flux);


    fclose (Diags);

//   ===============================================================================
//   Deallocate
//   ===============================================================================

    Deallocate_state (&State);
    Deallocate_state (&Truth);
    Deallocate_state (&Background);
    Deallocate_state (&Diff);
    Deallocate_metadata (&MetaData);
    DeallocateHorizTrans (&HorizData);
    Deallocate_CVT (&CVTdata);
    Array_1d_double_destroy (&dz);
    Array_1d_double_destroy (&density);

  }

  printf ("Done program\n");
}
