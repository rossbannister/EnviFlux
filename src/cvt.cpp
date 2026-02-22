   #include <source.h>


// -------------------------------------------------------------------------------
void cvt_total (struct state_type          *specsp,     // in  spectral space (in horiz)
                struct state_type          *realsp,     // out real space (in horiz)
                struct metadata_type       *MetaData,   // in  metadata
                struct HorizTransData_type *HorizData,  // in  info about horiz transform
                struct CVTData_type        *CVTdata)    // in  info about cvt
{ // Do the total transform
  // This assumes that all parameter fields have been allocated

  struct state_type post_horiz, intermediate_rs;

  Allocate_state (&post_horiz,
                  MetaData,
                  'r', (*specsp).vert_repres, (*specsp).temp_repres);

  Allocate_state (&intermediate_rs,
                  MetaData,
                  'r', 'r', 'r');

  // Do the horizontal transform -- this affects both tracer and source fields
  cvt_h ( specsp,      // in
          &post_horiz, // out
          HorizData,
          CVTdata );

  // Do the vertical transform -- this affects only the tracer field
  cvt_v (&intermediate_rs, // out
         &post_horiz,      // in
         CVTdata);

  // Do the temporal transform -- this affects only the source field
  cvt_t (&intermediate_rs, // out
         &post_horiz,      // in
         CVTdata);

  // Multiply by the standard deviations
  cvt_stddev (&intermediate_rs, // inout (real space)
              CVTdata);         // in  info about cvt

  // Deal with halos
  halos (&intermediate_rs);

  // Copy to output structure in real space
  copy_general (&intermediate_rs, // in  Before multiplying y standard deviations
                realsp,           // out copy of above
                HorizData,        // in  contains metadata
                true);            // in  include halos?


  // Tidy up
  Deallocate_state (&post_horiz);
  Deallocate_state (&intermediate_rs);
}


// -------------------------------------------------------------------------------
void cvt_total_adj (struct state_type          *specsp,     // out spectral space (in horiz)
                    struct state_type          *realsp,     // in  real space (in horiz)
                    struct metadata_type       *MetaData,   // in  metadata
                    struct HorizTransData_type *HorizData,  // in  info about horiz transform
                    struct CVTData_type        *CVTdata)    // in  info about cvt
{ // Do the total adjoint transform
  // This assumes that all parameter fields have been allocated

  struct state_type post_horiz, intermediate_rs;

  Allocate_state (&post_horiz,
                  MetaData,
                  'r', (*specsp).vert_repres, (*specsp).temp_repres);

  Allocate_state (&intermediate_rs,
                  MetaData,
                  'r', 'r', 'r');

  // Copy realsp to intermediate_rs
  copy_general (realsp,           // in   input state
                &intermediate_rs, // out  copy of above
                HorizData,        // in   contains metadata
                true);            // in   include halos?

  // Deal with halos
  halos_adj (&intermediate_rs);

  // Multiply by the standard deviations
  cvt_stddev (&intermediate_rs, // inout (real space)
              CVTdata);         // in  info about cvt

  // Do the vertical transform -- this affects only the tracer field
  cvt_v (&post_horiz,      // out
         &intermediate_rs, // in
         CVTdata);

  // Do the temporal transform -- this affects only the source field
  cvt_t_adj (&intermediate_rs, // in
             &post_horiz,      // out
             CVTdata);

  // Do the horizontal transform -- this affects both tracer and source fields
  cvt_h_adj ( specsp,      // out
              &post_horiz, // in
              HorizData,
              CVTdata );

  // Tidy up
  Deallocate_state (&post_horiz);
  Deallocate_state (&intermediate_rs);
}




// -------------------------------------------------------------------------------
void cvt_total_inv (struct state_type          *specsp,     // out spectral space (in horiz)
                    struct state_type          *realsp,     // in  real space (in horiz)
                    struct metadata_type       *MetaData,   // in  metadata
                    struct HorizTransData_type *HorizData,  // in  info about horiz transform
                    struct CVTData_type        *CVTdata)    // in  info about cvt
{ // Do the total inverse transform
  // This assumes that all parameter fields have been allocated

  struct state_type post_horiz, intermediate_rs;

  Allocate_state (&post_horiz,
                  MetaData,
                  'r', (*specsp).vert_repres, (*specsp).temp_repres);

  Allocate_state (&intermediate_rs,
                  MetaData,
                  'r', 'r', 'r');

  // Copy realsp to intermediate_rs
  copy_general (realsp,           // in   input state
                &intermediate_rs, // out  copy of above
                HorizData,        // in   contains metadata
                true);            // in   include halos?

  // Divide by the standard deviations
  cvt_stddev_inv (&intermediate_rs, // inout (real space)
                  CVTdata);         // in  info about cvt

  // Do the inverse of the vertical transform -- this affects only the tracer field
  cvt_v_inv (&post_horiz,      // out
             &intermediate_rs, // in
             CVTdata);

  // Do the inverse of the temporal transform -- this affects only the source field
  cvt_t_inv (&intermediate_rs, // in
             &post_horiz,      // out
             CVTdata);

  // Do the inverse of the horizontal transform -- this affects both tracer and source fields
  cvt_h_inv ( specsp,      // out
              &post_horiz, // in
              HorizData,
              CVTdata,
              'p' );       // pseudo inverse

  // Tidy up
  Deallocate_state (&post_horiz);
  Deallocate_state (&intermediate_rs);
}




// -------------------------------------------------------------------------------
void cvt_total_inv_adj ()
{
}




// -------------------------------------------------------------------------------
void cvt_h (struct state_type          *specsp,     // in  spectral space (in horiz)
            struct state_type          *realsp,     // out real space (in horiz)
            struct HorizTransData_type *HorizData,  // in  info about horiz transform
            struct CVTData_type        *CVTdata)    // in  info about cvt
{ // Perform the horizontal transform
  // This assumes that all parameter fields have been allocated

  int    z, zp1, t, l, m, cs, x, y, xp1, yp1;
  double weight;
  double ***interim_s, **interim_r;
  double *nil = NULL;

  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,      &((*realsp).nlon),
                 (*CVTdata).nlat,      &((*realsp).nlat),
                 (*CVTdata).nlev,      &((*realsp).nlev),
                 (*CVTdata).nss,       &((*realsp).nss),
                 (*CVTdata).L,         &((*realsp).L),
                 (*CVTdata).times,     (*realsp).times,
                 (*CVTdata).longitude, (*realsp).longitude,
                 (*CVTdata).latitude,  (*realsp).latitude,
                 (*CVTdata).level,     (*realsp).level,
                 nil,                  nil,
                 false);
  (*realsp).horiz_repres = 'r';
  (*realsp).vert_repres  = (*specsp).vert_repres;
  (*realsp).temp_repres  = (*specsp).temp_repres;

  // Allocate interim spectral field
  Array_3d_double_create (&interim_s,
                          2, (*specsp).L+1, (*specsp).L+1);
  // Allocate interim real field
  Array_2d_double_create (&interim_r,
                          (*specsp).nlon+2, (*specsp).nlat+2);

  // ----- Do initial field -----
  for (z=0; z<(*specsp).nlev; z++)
  { zp1 = z + 1;
    // Multiply by the horizontal spectrum (tracer)
    for (l=0; l<=(*specsp).L; l++)
    { weight = (*CVTdata).tracer_hspec[l][z];
      // m=0 cos term only
      interim_s[0][l][0] = (*specsp).tracer0_ss[0][l][0][zp1] * weight;
      if (l > 0)
      { // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { interim_s[cs][l][m] = (*specsp).tracer0_ss[cs][l][m][zp1] * weight;
          }
        }
      }
    }

    // Do spherical transform to real space
    spherical_for (interim_r,  // out
                   interim_s,  // in
                   HorizData);

    // Put data into correct part of the output
    for (x=0; x<(*realsp).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*realsp).nlat; y++)
      { yp1 = y + 1;
        (*realsp).tracer0_rs[xp1][yp1][zp1] = interim_r[xp1][yp1];
      }
    }
  }

  // ----- Do source field -----
  for (t=0; t<(*specsp).nss; t++)
  { // Multiply by the horizontal spectrum (flux)
    for (l=0; l<=(*specsp).L; l++)
    { weight = (*CVTdata).source_hspec[l];
      // m=0 cos term only
      interim_s[0][l][0] = (*specsp).source_ss[0][l][0][t] * weight;
      // Loop over remaining zonal wavenumbers
      for (m=1; m<=l; m++)
      { // cos and sin terms
        for (cs=0; cs<2; cs++)
        { interim_s[cs][l][m] = (*specsp).source_ss[cs][l][m][t] * weight;
        }
      }
    }

    // Do spherical transform to real space
    spherical_for (interim_r,  // out
                   interim_s,  // in
                   HorizData);

    // Put data into correct part of the output
    for (x=0; x<(*realsp).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*realsp).nlat; y++)
      { yp1 = y + 1;
        (*realsp).source_rs[xp1][yp1][t] = interim_r[xp1][yp1];
      }
    }
  }


  // Tidy up
  Array_3d_double_destroy (&interim_s,
                           2, (*specsp).L+1);
  Array_2d_double_destroy (&interim_r,
                          (*specsp).nlon+2);
}




// -------------------------------------------------------------------------------
void cvt_h_adj (struct state_type          *specsp,     // out spectral space (in horiz)
                struct state_type          *realsp,     // in  real space (in horiz)
                struct HorizTransData_type *HorizData,  // in  info about horiz transform
                struct CVTData_type        *CVTdata)    // in  info about cvt
{ // Perform the adjoint of the horizontal transform
  // This assumes that all parameter fields have been allocated

  int    z, zp1, t, l, m, cs, x, y, xp1, yp1;
  double weight;
  double ***interim_s, **interim_r;
  double *nil = NULL;

  // Copy over meta-data (in-to-out)
  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,      &((*specsp).nlon),
                 (*CVTdata).nlat,      &((*specsp).nlat),
                 (*CVTdata).nlev,      &((*specsp).nlev),
                 (*CVTdata).nss,       &((*specsp).nss),
                 (*CVTdata).L,         &((*specsp).L),
                 (*CVTdata).times,     (*specsp).times,
                 (*CVTdata).longitude, (*specsp).longitude,
                 (*CVTdata).latitude,  (*specsp).latitude,
                 (*CVTdata).level,     (*specsp).level,
                 nil,                  nil,
                 false);
  (*specsp).horiz_repres = 's';
  (*specsp).vert_repres  = (*realsp).vert_repres;
  (*specsp).temp_repres  = (*realsp).temp_repres;


  // Allocate interim spectral field
  Array_3d_double_create (&interim_s,
                          2, (*specsp).L+1, (*specsp).L+1);
  // Allocate interim real field
  Array_2d_double_create (&interim_r,
                          (*specsp).nlon+2, (*specsp).nlat+2);


  // ----- Do initial field -----
  for (z=0; z<(*specsp).nlev; z++)
  { zp1 = z + 1;
    // Put data into correct place for input into spherical adjoint operator
    for (x=0; x<(*realsp).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*realsp).nlat; y++)
      { yp1 = y + 1;
        interim_r[xp1][yp1] = (*realsp).tracer0_rs[xp1][yp1][zp1];
      }
    }

    // Do spherical transform to real space
    spherical_adj (interim_r,  // in
                   interim_s,  // out
                   HorizData);

    // Loop over total wavenumber
    for (l=0; l<=(*specsp).L; l++)
    { weight = (*CVTdata).tracer_hspec[l][z];
      // m=0 cos term only
      (*specsp).tracer0_ss[0][l][0][zp1] = interim_s[0][l][0] * weight;
      if (l > 0)
      { // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { (*specsp).tracer0_ss[cs][l][m][zp1] = interim_s[cs][l][m] * weight;
          }
        }
      }
    }
  }


  // ----- Do source field -----
  for (t=0; t<(*specsp).nss; t++)
  { // Put data into correct place for input into spherical adjoint operator
    for (x=0; x<(*realsp).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*realsp).nlat; y++)
      { yp1 = y + 1;
        interim_r[xp1][yp1] = (*realsp).source_rs[xp1][yp1][t];
      }
    }

    // Do spherical transform to real space
    spherical_adj (interim_r,  // in
                   interim_s,  // out
                   HorizData);

    // Multiply by the horizontal spectrum (tracer)
    for (l=0; l<=(*specsp).L; l++)
    { weight = (*CVTdata).source_hspec[l];
      // m=0 cos term only
      (*specsp).source_ss[0][l][0][t] = interim_s[0][l][0] * weight;
      // Loop over remaining zonal wavenumbers
      for (m=1; m<=l; m++)
      { // cos and sin terms
        for (cs=0; cs<2; cs++)
        { (*specsp).source_ss[cs][l][m][t] = interim_s[cs][l][m] * weight;
        }
      }
    }
  }

  // Tidy up
  Array_3d_double_destroy (&interim_s,
                           2, (*specsp).L+1);
  Array_2d_double_destroy (&interim_r,
                          (*specsp).nlon+2);
}



// -------------------------------------------------------------------------------
void cvt_h_inv (struct state_type          *specsp,     // out spectral space (in horiz)
                struct state_type          *realsp,     // in  real space (in horiz)
                struct HorizTransData_type *HorizData,  // in  info about horiz transform
                struct CVTData_type        *CVTdata,    // in  info about cvt
                char                       exactpseudo) // in  'e' exact inverse
                                                        //     'p' pseudo inverse
{ // Perform the inverse of the horizontal transform
  // This assumes that all parameter fields have been allocated

  int    z, zp1, t, l, m, cs, x, y, xp1, yp1;
  double weight;
  double ***interim_s, **interim_r;
  double *nil = NULL;

  // Copy over meta-data (in-to-out)
  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,      &((*specsp).nlon),
                 (*CVTdata).nlat,      &((*specsp).nlat),
                 (*CVTdata).nlev,      &((*specsp).nlev),
                 (*CVTdata).nss,       &((*specsp).nss),
                 (*CVTdata).L,         &((*specsp).L),
                 (*CVTdata).times,     (*specsp).times,
                 (*CVTdata).longitude, (*specsp).longitude,
                 (*CVTdata).latitude,  (*specsp).latitude,
                 (*CVTdata).level,     (*specsp).level,
                 nil,                  nil,
                 false);
  (*specsp).horiz_repres = 's';
  (*specsp).vert_repres  = (*realsp).vert_repres;
  (*specsp).temp_repres  = (*realsp).temp_repres;


  // Allocate interim spectral field
  Array_3d_double_create (&interim_s,
                          2, (*specsp).L+1, (*specsp).L+1);
  // Allocate interim real field
  Array_2d_double_create (&interim_r,
                          (*specsp).nlon+2, (*specsp).nlat+2);


  // ----- Do initial field -----
  for (z=0; z<(*specsp).nlev; z++)
  { zp1 = z + 1;
    // Put data into correct place for input into spherical inverse operator
    for (x=0; x<(*realsp).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*realsp).nlat; y++)
      { yp1 = y + 1;
        interim_r[xp1][yp1] = (*realsp).tracer0_rs[xp1][yp1][zp1];
      }
    }

    // Do spherical transform to real space
    spherical_inv (interim_r,  // in
                   interim_s,  // out
                   HorizData);

    // Drivide by the horizontal spectrum (tracer)
    for (l=0; l<=(*specsp).L; l++)
    { weight = (*CVTdata).tracer_hspec[l][z];
      if ((exactpseudo == 'e') || (weight != 0.0))
      { // m=0 cos term only
        (*specsp).tracer0_ss[0][l][0][zp1] = interim_s[0][l][0] / weight;
        if (l > 0)
        { // Loop over remaining zonal wavenumbers
          for (m=1; m<=l; m++)
          { // cos and sin terms
            for (cs=0; cs<2; cs++)
            { (*specsp).tracer0_ss[cs][l][m][zp1] = interim_s[cs][l][m] / weight;
            }
          }
        }
      }
    }
  }


  // ----- Do source field -----
  for (t=0; t<(*specsp).nss; t++)
  { // Put data into correct place for input into spherical adjoint operator
    for (x=0; x<(*realsp).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*realsp).nlat; y++)
      { yp1 = y + 1;
        interim_r[xp1][yp1] = (*realsp).source_rs[xp1][yp1][t];
      }
    }

    // Do spherical transform to real space
    spherical_inv (interim_r,  // in
                   interim_s,  // out
                   HorizData);

    // Divide by the horizontal spectrum (flux)
    for (l=0; l<=(*specsp).L; l++)
    { weight = (*CVTdata).source_hspec[l];
      if ((exactpseudo == 'e') || (weight != 0.0))
      { // m=0 cos term only
        (*specsp).source_ss[0][l][0][t] = interim_s[0][l][0] / weight;
        // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { (*specsp).source_ss[cs][l][m][t] = interim_s[cs][l][m] / weight;
          }
        }
      }
    }
  }

  // Tidy up
  Array_3d_double_destroy (&interim_s,
                           2, (*specsp).L+1);
  Array_2d_double_destroy (&interim_r,
                          (*specsp).nlon+2);
}



// -------------------------------------------------------------------------------
void cvt_v (struct state_type          *output,     // out {both level space (in vert)
            struct state_type          *input,      // in  {and real space (in horiz)
            struct CVTData_type        *CVTdata)    // in  info about cvt
{ // Perform the vertical transform
  // This assumes that all parameter fields have been allocated
  // This operator is self-adjoint

  int    x, y, mode, z, xp1, yp1;
  double sum;
  double tracer_modes[(*input).nlev];
  double *nil = NULL;


  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,      &((*output).nlon),
                 (*CVTdata).nlat,      &((*output).nlat),
                 (*CVTdata).nlev,      &((*output).nlev),
                 (*CVTdata).nss,       &((*output).nss),
                 (*CVTdata).L,         &((*output).L),
                 (*CVTdata).times,     (*output).times,
                 (*CVTdata).longitude, (*output).longitude,
                 (*CVTdata).latitude,  (*output).latitude,
                 (*CVTdata).level,     (*output).level,
                 nil,                  nil,
                 false);
  (*output).horiz_repres = (*input).horiz_repres;
  (*output).vert_repres  = 'r';

  for (x=0; x<(*input).nlon; x++)
  { xp1 = x + 1;
    for (y=0; y<(*input).nlat; y++)
    { yp1 = y + 1;
      // Project onto vertical modes
      for (mode=0; mode<(*input).nlev; mode++)
      { sum = 0.0;
        for (z=0; z<(*input).nlev; z++)
        { sum += (*CVTdata).tracer_vertmodes[z][mode] * (*input).tracer0_rs[xp1][yp1][z+1];
        }
        tracer_modes[mode] = sum;
      }
      // Multiply by standard deviation (a function of vertical mode)
      for (mode=0; mode<(*input).nlev; mode++)
      { tracer_modes[mode] *= (*CVTdata).tracer_vspec[mode];
      }
      // Project back to levels
      for (z=0; z<(*input).nlev; z++)
      { sum = 0.0;
        for (mode=0; mode<(*input).nlev; mode++)
        { sum += (*CVTdata).tracer_vertmodes[z][mode] * tracer_modes[mode];
        }
        (*output).tracer0_rs[xp1][yp1][z+1] = sum;
      }
    }
  }
}



// -------------------------------------------------------------------------------
void cvt_v_inv (struct state_type          *output,     // out {both level space (in vert)
                struct state_type          *input,      // in  {and real space (in horiz)
                struct CVTData_type        *CVTdata)    // in  info about cvt
{ // Perform the inverse vertical transform
  // This assumes that all parameter fields have been allocated
  // This operator is self-adjoint

  int    x, y, mode, z, xp1, yp1;
  double sum;
  double tracer_modes[(*input).nlev];
  double *nil = NULL;


  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,      &((*output).nlon),
                 (*CVTdata).nlat,      &((*output).nlat),
                 (*CVTdata).nlev,      &((*output).nlev),
                 (*CVTdata).nss,       &((*output).nss),
                 (*CVTdata).L,         &((*output).L),
                 (*CVTdata).times,     (*output).times,
                 (*CVTdata).longitude, (*output).longitude,
                 (*CVTdata).latitude,  (*output).latitude,
                 (*CVTdata).level,     (*output).level,
                 nil,                  nil,
                 false);
  (*output).horiz_repres = (*input).horiz_repres;
  (*output).vert_repres  = 'r';

  for (x=0; x<(*input).nlon; x++)
  { xp1 = x + 1;
    for (y=0; y<(*input).nlat; y++)
    { yp1 = y + 1;
      // Project onto vertical modes
      for (mode=0; mode<(*input).nlev; mode++)
      { sum = 0.0;
        for (z=0; z<(*input).nlev; z++)
        { sum += (*CVTdata).tracer_vertmodes[z][mode] * (*input).tracer0_rs[xp1][yp1][z+1];
        }
        tracer_modes[mode] = sum;
      }
      // Divide by standard deviation (a function of vertical mode)
      for (mode=0; mode<(*input).nlev; mode++)
      { tracer_modes[mode] /= (*CVTdata).tracer_vspec[mode];
      }
      // Project back to levels
      for (z=0; z<(*input).nlev; z++)
      { sum = 0.0;
        for (mode=0; mode<(*input).nlev; mode++)
        { sum += (*CVTdata).tracer_vertmodes[z][mode] * tracer_modes[mode];
        }
        (*output).tracer0_rs[xp1][yp1][z+1] = sum;
      }
    }
  }
}



// -------------------------------------------------------------------------------
void cvt_t (struct state_type   *fntime,     // out source fn of time          {both real space
            struct state_type   *fntempmode, // in  source fn of temporal mode {in horiz
            struct CVTData_type *CVTdata)    // in  info about cvt
{ // Perform the temporal transform
  // This assumes that all parameter fields have been allocated
  int    x, y, t, mode, xp1, yp1;
  double sum;
  double flux_modes[(*fntempmode).nss];
  double *nil = NULL;


  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,         &((*fntime).nlon),
                 (*CVTdata).nlat,         &((*fntime).nlat),
                 (*CVTdata).nlev,         &((*fntime).nlev),
                 (*CVTdata).nss,          &((*fntime).nss),
                 (*CVTdata).L,            &((*fntime).L),
                 (*CVTdata).times,        (*fntime).times,
                 (*CVTdata).longitude,    (*fntime).longitude,
                 (*CVTdata).latitude,     (*fntime).latitude,
                 (*CVTdata).level,        (*fntime).level,
                 nil,                     nil,
                 false);
  (*fntime).horiz_repres = (*fntempmode).horiz_repres;
  (*fntime).temp_repres  = 'r';


  if ((*CVTdata).temporal_covs)
  { // We are using temporal correlations
    for (x=0; x<(*fntempmode).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*fntempmode).nlat; y++)
      { yp1 = y + 1;
        // Multipy by temporal standard deviations (fn of timescale)
        for (mode=0; mode<(*fntempmode).nss; mode++)
        { flux_modes[mode] = (*fntempmode).source_rs[xp1][yp1][mode] * (*CVTdata).source_tspec[mode];
        }
        // Project from temporal modes to time
        for (t=0; t<(*fntempmode).nss; t++)
        { sum = 0.0;
          for (mode=0; mode<(*fntempmode).nss; mode++)
          { sum += (*CVTdata).source_tempmodes[t][mode] * flux_modes[mode];
          }
          (*fntime).source_rs[xp1][yp1][t] = sum;
        }
      }
    }
  }
  else
  { // No temporal correlations - just copy
    for (x=0; x<(*fntempmode).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*fntempmode).nlat; y++)
      { yp1 = y + 1;
        for (t=0; t<(*fntempmode).nss; t++)
        { (*fntime).source_rs[xp1][yp1][t] = (*fntempmode).source_rs[xp1][yp1][t];
        }
      }
    }
  }
}



// -------------------------------------------------------------------------------
void cvt_t_adj (struct state_type   *fntime,     // in  source fn of time          {both real space
                struct state_type   *fntempmode, // out source fn of temporal mode {in horiz
                struct CVTData_type *CVTdata)    // in  info about cvt
{ // Perform the adjoint of the temporal transform
  // This assumes that all parameter fields have been allocated
  int    x, y, t, mode, xp1, yp1;
  double sum;
  double flux_modes[(*fntime).nss];
  double *nil = NULL;

  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,      &((*fntempmode).nlon),
                 (*CVTdata).nlat,      &((*fntempmode).nlat),
                 (*CVTdata).nlev,      &((*fntempmode).nlev),
                 (*CVTdata).nss,       &((*fntempmode).nss),
                 (*CVTdata).L,         &((*fntempmode).L),
                 (*CVTdata).times,     (*fntempmode).times,
                 (*CVTdata).longitude, (*fntempmode).longitude,
                 (*CVTdata).latitude,  (*fntempmode).latitude,
                 (*CVTdata).level,     (*fntempmode).level,
                 nil,                  nil,
                 false);
  (*fntempmode).horiz_repres = (*fntime).horiz_repres;


  if ((*CVTdata).temporal_covs)
  { // We are using temporal correlations
    (*fntempmode).temp_repres  = 'm';
    for (x=0; x<(*fntime).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*fntime).nlat; y++)
      { yp1 = y + 1;
        // Project from time to temporal modes
        for (mode=0; mode<(*fntime).nss; mode++)
        { sum = 0.0;
          for (t=0; t<(*fntime).nss; t++)
          { sum += (*CVTdata).source_tempmodes[t][mode] * (*fntime).source_rs[xp1][yp1][t];
          }
          flux_modes[mode] = sum;
        }
        // Multipy by temporal standard deviations (fn of timescale)
        for (mode=0; mode<(*fntime).nss; mode++)
        { (*fntempmode).source_rs[xp1][yp1][mode] = flux_modes[mode] * (*CVTdata).source_tspec[mode];
        }
      }
    }
  }
  else
  { // No temporal correlations - just copy
    (*fntempmode).temp_repres  = 'r';
    for (x=0; x<(*fntime).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*fntime).nlat; y++)
      { yp1 = y + 1;
        for (t=0; t<(*fntime).nss; t++)
        { (*fntempmode).source_rs[xp1][yp1][t] = (*fntime).source_rs[xp1][yp1][t];
        }
      }
    }
  }
}



// -------------------------------------------------------------------------------
void cvt_t_inv (struct state_type   *fntime,     // in  source fn of time          {both real space
                struct state_type   *fntempmode, // out source fn of temporal mode {in horiz
                struct CVTData_type *CVTdata)    // in  info about cvt
{ // Perform the inverse of the temporal transform
  // This assumes that all parameter fields have been allocated
  int    x, y, t, mode, xp1, yp1;
  double sum;
  double flux_modes[(*fntime).nss];
  double *nil = NULL;

  // Copy over meta-data (in-to-out)
  copy_metadata ((*CVTdata).nlon,      &((*fntempmode).nlon),
                 (*CVTdata).nlat,      &((*fntempmode).nlat),
                 (*CVTdata).nlev,      &((*fntempmode).nlev),
                 (*CVTdata).nss,       &((*fntempmode).nss),
                 (*CVTdata).L,         &((*fntempmode).L),
                 (*CVTdata).times,     (*fntempmode).times,
                 (*CVTdata).longitude, (*fntempmode).longitude,
                 (*CVTdata).latitude,  (*fntempmode).latitude,
                 (*CVTdata).level,     (*fntempmode).level,
                 nil,                  nil,
                 false);
  (*fntempmode).horiz_repres = (*fntime).horiz_repres;


  if ((*CVTdata).temporal_covs)
  { // We are using temporal correlations
    (*fntempmode).temp_repres  = 'm';
    for (x=0; x<(*fntime).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*fntime).nlat; y++)
      { yp1 = y + 1;
        // Project from time to temporal modes
        for (mode=0; mode<(*fntime).nss; mode++)
        { sum = 0.0;
          for (t=0; t<(*fntime).nss; t++)
          { sum += (*CVTdata).source_tempmodes[t][mode] * (*fntime).source_rs[xp1][yp1][t];
          }
          flux_modes[mode] = sum;
        }
        // Divide by temporal standard deviations (fn of timescale)
        for (mode=0; mode<(*fntime).nss; mode++)
        { (*fntempmode).source_rs[xp1][yp1][mode] = flux_modes[mode] / (*CVTdata).source_tspec[mode];
        }
      }
    }
  }
  else
  { // No temporal correlations - just copy
    (*fntempmode).temp_repres  = 'r';
    for (x=0; x<(*fntime).nlon; x++)
    { xp1 = x + 1;
      for (y=0; y<(*fntime).nlat; y++)
      { yp1 = y + 1;
        for (t=0; t<(*fntime).nss; t++)
        { (*fntempmode).source_rs[xp1][yp1][t] = (*fntime).source_rs[xp1][yp1][t];
        }
      }
    }
  }
}


// -------------------------------------------------------------------------------
void cvt_stddev (struct state_type   *state,     // inout (real space)
                 struct CVTData_type *CVTdata)   // in  info about cvt
{ // Multiply by the standard deviation
  int    x, y, z, t, xp1, yp1, zp1;

  // Work on the tracer field
  for (x=0; x<(*state).nlon; x++)
  { xp1 = x + 1;
    for (y=0; y<(*state).nlat; y++)
    { yp1 = y + 1;
      for (z=0; z<(*state).nlev; z++)
      { zp1 = z + 1;
        (*state).tracer0_rs[xp1][yp1][zp1] *= (*CVTdata).tracer_stddev[z];
      }
    }
  }

  // Work on the source field
  for (x=0; x<(*state).nlon; x++)
  { xp1 = x + 1;
    for (y=0; y<(*state).nlat; y++)
    { yp1 = y + 1;
      for (t=0; t<(*state).nss; t++)
      { (*state).source_rs[xp1][yp1][t] *= (*CVTdata).source_stddev[x][y][t];
      }
    }
  }
}


// -------------------------------------------------------------------------------
void cvt_stddev_inv (struct state_type   *state,     // inout (real space)
                     struct CVTData_type *CVTdata)   // in  info about cvt
{ // Divide by the standard deviation
  int    x, y, z, t, xp1, yp1, zp1;

  // Work on the tracer field
  for (x=0; x<(*state).nlon; x++)
  { xp1 = x + 1;
    for (y=0; y<(*state).nlat; y++)
    { yp1 = y + 1;
      for (z=0; z<(*state).nlev; z++)
      { zp1 = z + 1;
        (*state).tracer0_rs[xp1][yp1][zp1] /= (*CVTdata).tracer_stddev[z];
      }
    }
  }

  // Work on the source field
  for (x=0; x<(*state).nlon; x++)
  { xp1 = x + 1;
    for (y=0; y<(*state).nlat; y++)
    { yp1 = y + 1;
      for (t=0; t<(*state).nss; t++)
      { (*state).source_rs[xp1][yp1][t] /= (*CVTdata).source_stddev[x][y][t];
      }
    }
  }
}


// -------------------------------------------------------------------------------
void PutWhiteNoiseControlVector (struct state_type *controlsp,
                                 int               *rndseq)
{ // Put N(0,1) uncorrelated white noise into control structure

  int z, zp1, l, m, cs, t;

  // ----- Add uncorrelated white noise to initial field -----
  for (z=0; z<(*controlsp).nlev; z++)
  { zp1 = z + 1;
    // Multiply by the horizontal spectrum (tracer)
    for (l=0; l<=(*controlsp).L; l++)
    { // m=0 cos term only
      (*controlsp).tracer0_ss[0][l][0][zp1] = Normal(0.0,
                                                     1.0,
                                                     rndseq);
      if (l > 0)
      { // Loop over remaining zonal wavenumbers
        for (m=1; m<=l; m++)
        { // cos and sin terms
          for (cs=0; cs<2; cs++)
          { (*controlsp).tracer0_ss[cs][l][m][zp1] = Normal(0.0,
                                                            1.0,
                                                            rndseq);
          }
        }
      }
    }
  }

  // ----- Add uncorrelated white noise to source field -----
  for (t=0; t<(*controlsp).nss; t++)
  { // Multiply by the horizontal spectrum (tracer)
    for (l=0; l<=(*controlsp).L; l++)
    { // m=0 cos term only
      (*controlsp).source_ss[0][l][0][t] = Normal(0.0,
                                                  1.0,
                                                  rndseq);
      // Loop over remaining zonal wavenumbers
      for (m=1; m<=l; m++)
      { // cos and sin terms
        for (cs=0; cs<2; cs++)
        { (*controlsp).source_ss[cs][l][m][t] = Normal(0.0,
                                                       1.0,
                                                       rndseq);
        }
      }
    }
  }

}
