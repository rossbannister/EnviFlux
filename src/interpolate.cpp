   #include <source.h>


// -------------------------------------------------------------------------------
double Interpolate1D ( double field[2], // The 1D field to be interpolated
                       double axis[2],  // The axis values
                       double value)    // The point to interpolate to
{ // Interpolate in 1D (linear interpolation)
  // Used in the semi-Lagrangian scheme
  double dx, dy, dxpoint;

  dx      = axis[1] - axis[0];
  dy      = field[1] - field[0];
  dxpoint = value - axis[0];

  return field[0] + dxpoint * dy / dx;
}

// -------------------------------------------------------------------------------
void Interpolate1D_adj ( double field_hat[2],     // The 1D field to be interpolated (out)
                         double interpolate_hat,  // The interpolated value (in)
                         double axis[2],          // The axis values
                         double value)            // The point to interpolate to
{ // Interpolate in 1D (linear interpolation, adjoint routine)
  // Used in the semi-Lagrangian scheme
  double dx, dy_hat, dxpoint;

  dx            = axis[1] - axis[0];
  dxpoint       = value - axis[0];
  field_hat[0]  = interpolate_hat;
  dy_hat        = interpolate_hat * dxpoint / dx;
  field_hat[1]  = dy_hat;
  field_hat[0] -= dy_hat;
}


// -------------------------------------------------------------------------------
double Interpolate2D ( double field[2][2], // The 2D field to be interpolated
                       double axisx[2],    // The axis x values
                       double axisy[2],    // The axis y values
                       double valuex,      // The x point to interpolate to
                       double valuey)      // The y point to interpolate to
{ // Interpolate in 2D (bi-linear interpolation)
  // Used in the semi-Lagrangian scheme
  int    xpos;
  double val[2];

  // Interpolate in the y direction
  for (xpos=0; xpos<2; xpos++)
  { val[xpos] = Interpolate1D (field[xpos],
                               axisy,
                               valuey);
  }

  // Interpolate in the x direction
  return Interpolate1D (val,
                        axisx,
                        valuex);
}

// -------------------------------------------------------------------------------
void Interpolate2D_adj ( double field_hat[2][2],  // The 2D field to be interpolated (out)
                         double interpolate_hat,  // The interpolated value (in)
                         double axisx[2],         // The axis x values
                         double axisy[2],         // The axis y values
                         double valuex,           // The x point to interpolate to
                         double valuey)           // The y point to interpolate to
{ // Interpolate in 2D (bi-linear interpolation, adjoint routine)
  // Used in the semi-Lagrangian scheme
  int    xpos;
  double val_hat[2];

  // Interpolate in the x direction
  Interpolate1D_adj (val_hat,          // out
                     interpolate_hat,  // in
                     axisx,
                     valuex);

  // Interpolate in the y direction
  for (xpos=0; xpos<2; xpos++)
  { Interpolate1D_adj (field_hat[xpos],  // out
                       val_hat[xpos],    // in
                       axisy,
                       valuey);
  }

}


// -------------------------------------------------------------------------------
double Interpolate3D ( double field[2][2][2],  // The 3D field to be interpolated
                       double axisx[2],        // The axis x values
                       double axisy[2],        // The axis y values
                       double axisz[2],        // The axis z values
                       double valuex,          // The x point to interpolate to
                       double valuey,          // The y point to interpolate to
                       double valuez )         // The y point to interpolate to
{ // Interpolate in 3D (tri-linear interpolation)
  // Used in the semi-Lagrangian scheme
  int    xpos;
  double val[2];

  // Interpolate in the y/z direction
  for (xpos=0; xpos<2; xpos++)
  { val[xpos] = Interpolate2D (field[xpos],
                               axisy,
                               axisz,
                               valuey,
                               valuez);
  }

  // Interpolate in the x direction
  return Interpolate1D (val,
                        axisx,
                        valuex);
}



// -------------------------------------------------------------------------------
void Interpolate3D_adj ( double field_hat[2][2][2],  // The 3D field to be interpolated (out)
                         double interpolate_hat,     // The interpolated value (in)
                         double axisx[2],        // The axis x values
                         double axisy[2],        // The axis y values
                         double axisz[2],        // The axis z values
                         double valuex,          // The x point to interpolate to
                         double valuey,          // The y point to interpolate to
                         double valuez )         // The y point to interpolate to
{ // Interpolate in 3D (tri-linear interpolation, adjoint routine)
  // Used in the semi-Lagrangian scheme
  int    xpos;
  double val_hat[2];

  // Interpolate in the x direction
  Interpolate1D_adj (val_hat,          // out
                     interpolate_hat,  // in
                     axisx,
                     valuex);

  // Interpolate in the y/z direction
  for (xpos=0; xpos<2; xpos++)
  { Interpolate2D_adj (field_hat[xpos],  // out
                       val_hat[xpos],    // in
                       axisy,
                       axisz,
                       valuey,
                       valuez);
  }
}


// -------------------------------------------------------------------------------
double Interpolate3Dt ( double fieldt0[2][2][2],  // The 3D field to be interpolated at time 0
                        double fieldt1[2][2][2],  // The 3D field to be interpolated at time 1
                        double axisx[2],          // The axis x values
                        double axisy[2],          // The axis y values
                        double axisz[2],          // The axis z values
                        double axist[2],          // The axis time values
                        double valuex,            // The x point to interpolate to
                        double valuey,            // The y point to interpolate to
                        double valuez,            // The z point to interpolate to
                        double valuet)            // The time point to interpolate to
{ // Interpolate in 3D and time
  // Used in the semi-Lagrangian scheme
  double val[2];

  // Interpolate to the spatial point for each field
  val[0] = Interpolate3D (fieldt0,
                          axisx,
                          axisy,
                          axisz,
                          valuex,
                          valuey,
                          valuez);
  val[1] = Interpolate3D (fieldt1,
                          axisx,
                          axisy,
                          axisz,
                          valuex,
                          valuey,
                          valuez);

  // Interpolate in time
  return Interpolate1D (val,
                        axist,
                        valuet);
}


// -------------------------------------------------------------------------------
void Interpolate3Dt_adj ( double fieldt0_hat[2][2][2],  // The 3D field to be interpolated at time 0 (out)
                          double fieldt1_hat[2][2][2],  // The 3D field to be interpolated at time 1 (out)
                          double interpolate_hat,       // The interpolated value (in)
                          double axisx[2],              // The axis x values
                          double axisy[2],              // The axis y values
                          double axisz[2],              // The axis z values
                          double axist[2],              // The axis time values
                          double valuex,                // The x point to interpolate to
                          double valuey,                // The y point to interpolate to
                          double valuez,                // The z point to interpolate to
                          double valuet)                // The time point to interpolate to
{ // Interpolate in 3D and time
  // Used in the semi-Lagrangian scheme, adjoint routine
  double val_hat[2];

  // Interpolate in time
  Interpolate1D_adj (val_hat,         // out
                     interpolate_hat, // in
                     axist,
                     valuet);

  // Interpolate to the spatial point for each field
  Interpolate3D_adj (fieldt0_hat,
                     val_hat[0],
                     axisx,
                     axisy,
                     axisz,
                     valuex,
                     valuey,
                     valuez);
  Interpolate3D_adj (fieldt1_hat,
                     val_hat[1],
                     axisx,
                     axisy,
                     axisz,
                     valuex,
                     valuey,
                     valuez);
}

// -------------------------------------------------------------------------------
double interpolate_1d_cubic ( double x[6],     // grid positions
                              double vals[6],  // field values (in)
                              double pos,      // position to interpolate to
                              int    parts )   // divide domain into 2 or 3 parts
{ // Interpolate in 1d to cubic order.  Based on Stirling's and Bessel's formulae, Lanczos, p.311
  // See also my notes CubicInterpolation.pdf
  // Used in the semi-Lagrangian scheme
  int    l, l2, l2p1, l2p2;
  char   lmr;
  double Dx, recip_Dx, Dx2d4, deltax, deltax2, dx2mDx2d4, dist_from_cent, result;
  double x_inc_halves[11];    // x at all points inc half points
  double vals_inc_halves[11]; // Values at all points inc half points
  double ddx[7];              // dval/dx at all relevant points, inc half points
  double d2dx2[5];            // d2val/dx2 at all relevant points, inc half points
  double d3dx3[3];            // d3val/dx3 at final half point

  Dx       = (x[5] - x[0]) / 5.0;
  recip_Dx = 1.0 / Dx;
  Dx2d4    = Dx * Dx / 4.0;

  // Fill the basic data at full points
  for (l=0; l<6; l++)
  { l2                  = 2*l;
    x_inc_halves[l2]    = x[l];
    vals_inc_halves[l2] = vals[l];
  }
  // Fill the basic data at half points
  for (l=0; l<5; l++)
  { l2                    = 2*l;
    l2p1                  = l2 + 1;
    l2p2                  = l2 + 2;
    x_inc_halves[l2p1]    = 0.5 * (x_inc_halves[l2] + x_inc_halves[l2p2]);
    vals_inc_halves[l2p1] = 0.5 * (vals_inc_halves[l2] + vals_inc_halves[l2p2]);
  }

  //for (l=0; l<11; l++)
  //{ printf ("x and y inc half points: %f  %f\n", x_inc_halves[l], vals_inc_halves[l]);
  //}

  // Find first derivatives
  for (l=0; l<7; l++)
  { ddx[l] = recip_Dx * (vals_inc_halves[l+3] - vals_inc_halves[l+1]);
    //printf ("dy/dx %f\n", ddx[l]);
  }

  // Find second derivatives
  for (l=0; l<5; l++)
  { d2dx2[l] = recip_Dx * (ddx[l+2] - ddx[l]);
    //printf ("d2y/dx2 %f\n", d2dx2[l]);
  }

  // Find third derivative
  for (l=0; l<3; l++)
  { d3dx3[l] = recip_Dx * (d2dx2[l+2] - d2dx2[l]);
    //printf ("d3dx3 %f\n", d3dx3[l]);
  }



  // Do cubic interpolation
  // Decide whether to use Stirling's or Bessel's formula

  dist_from_cent = (pos - x_inc_halves[5]) / Dx;
  if (parts == 3)
  { // Stirling - Bessel - Stirling (3 parts)
    if (fabs(dist_from_cent) > 0.25)
    { // One of Stirling's formulae
      if (dist_from_cent < 0.0)
      { lmr = 'l';  // Code for Stirling left
      }
      else
      { lmr = 'r';  // Code for Stirling right
      }
    }
    else
    { lmr = 'm';  // Code for Bessel
    }
  }
  else
  { // Stirling - Stirling (2 parts)
    if (dist_from_cent < 0.0)
    { lmr = 'l';  // Code for Stirling left
    }
    else
    { lmr = 'r';  // Code for Stirling right
    }
  }

  switch (lmr)
  { case 'l':
      // Left Stirling's formula
      deltax  = pos - x_inc_halves[4];
      deltax2 = deltax * deltax;
      result  = vals_inc_halves[4] + ddx[2] * deltax + d2dx2[1] * deltax2 / 2.0 +
                d3dx3[0] * deltax * (deltax2  - Dx * Dx) / 6.0;
      break;
    case 'r':
      // Right Stirling's formula
      deltax  = pos - x_inc_halves[6];
      deltax2 = deltax * deltax;
      result  = vals_inc_halves[6] + ddx[4] * deltax + d2dx2[3] * deltax2 / 2.0 +
                d3dx3[2] * deltax * (deltax2  - Dx * Dx) / 6.0;
      break;
    default:  // 'm'
      // Bessel's formula
      deltax    = pos - x_inc_halves[5];
      dx2mDx2d4 = deltax * deltax - Dx2d4;
      result    = vals_inc_halves[5] + ddx[3] * deltax + d2dx2[2] * dx2mDx2d4 / 2.0 +
                  d3dx3[1] * deltax * dx2mDx2d4 / 6.0;
  }

  return result;
}




// -------------------------------------------------------------------------------
void interpolate_1d_cubic_adj ( double x[6],            // grid positions
                                double vals_hat[6],     // field values (out)
                                double interpolate_hat, // interpolated value (in)
                                double pos,             // position to interpolate to
                                int    parts )          // divide domain into 2 or 3 parts
{ // Interpolate in 1d to cubic order (adjoint routine).  Based on Stirling's and Bessel's formulae, Lanczos, p.311
  // See also my notes CubicInterpolation.pdf
  // Used in the semi-Lagrangian scheme
  int    l, l2, l2p1, l2p2;
  char   lmr;
  double Dx, recip_Dx, Dx2d4, deltax, deltax2, dx2mDx2d4, dist_from_cent, result;
  double x_inc_halves[11];        // x at all points inc half points
  double vals_inc_halves_hat[11]; // Values at all points inc half points
  double ddx_hat[7];              // dval/dx at all relevant points, inc half points
  double d2dx2_hat[5];            // d2val/dx2 at all relevant points, inc half points
  double d3dx3_hat[3];            // d3val/dx3 at final half point

  // Initialise hat variables
  for (l=0; l<11; l++)
  { vals_inc_halves_hat[l] = 0.0;
  }
  for (l=0; l<7; l++)
  { ddx_hat[l] = 0.0;
  }
  for (l=0; l<5; l++)
  { d2dx2_hat[l] = 0.0;
  }
  for (l=0; l<3; l++)
  { d3dx3_hat[l] = 0.0;
  }

  Dx       = (x[5] - x[0]) / 5.0;
  recip_Dx = 1.0 / Dx;
  Dx2d4    = Dx * Dx / 4.0;

  // Fill the basic data at full points
  for (l=0; l<6; l++)
  { l2               = 2*l;
    x_inc_halves[l2] = x[l];
  }

  // Fill the basic data at half points
  for (l=0; l<5; l++)
  { l2                 = 2*l;
    l2p1               = l2 + 1;
    l2p2               = l2 + 2;
    x_inc_halves[l2p1] = 0.5 * (x_inc_halves[l2] + x_inc_halves[l2p2]);
  }

  // Do cubic interpolation
  // Decide whether to use Stirling's or Bessel's formula

  dist_from_cent = (pos - x_inc_halves[5]) / Dx;
  if (parts == 3)
  { // Stirling - Bessel - Stirling (3 parts)
    if (fabs(dist_from_cent) > 0.25)
    { // One of Stirling's formulae
      if (dist_from_cent < 0.0)
      { lmr = 'l';  // Code for Stirling left
      }
      else
      { lmr = 'r';  // Code for Stirling right
      }
    }
    else
    { lmr = 'm';  // Code for Bessel
    }
  }
  else
  { // Stirling - Stirling (2 parts)
    if (dist_from_cent < 0.0)
    { lmr = 'l';  // Code for Stirling left
    }
    else
    { lmr = 'r';  // Code for Stirling right
    }
  }

  switch (lmr)
  { case 'l':
      // Left Stirling's formula
      deltax  = pos - x_inc_halves[4];
      deltax2 = deltax * deltax;
      // Forward routine
      // result = vals_inc_halves[4] + ddx[2] * deltax + d2dx2[1] * deltax2 / 2.0 +
      //          d3dx3[0] * deltax * (deltax2  - Dx * Dx) / 6.0;
      vals_inc_halves_hat[4] = interpolate_hat;
      ddx_hat[2]             = interpolate_hat * deltax;
      d2dx2_hat[1]           = interpolate_hat * deltax2 / 2.0;
      d3dx3_hat[0]           = interpolate_hat * deltax * (deltax2  - Dx * Dx) / 6.0;
      break;
    case 'r':
      // Right Stirling's formula
      deltax  = pos - x_inc_halves[6];
      deltax2 = deltax * deltax;
      // Forward routine
      // result  = vals_inc_halves[6] + ddx[4] * deltax + d2dx2[3] * deltax2 / 2.0 +
      //           d3dx3[2] * deltax * (deltax2  - Dx * Dx) / 6.0;
      vals_inc_halves_hat[6] = interpolate_hat;
      ddx_hat[4]             = interpolate_hat * deltax;
      d2dx2_hat[3]           = interpolate_hat * deltax2 / 2.0;
      d3dx3_hat[2]           = interpolate_hat * deltax * (deltax2  - Dx * Dx) / 6.0;
      break;
    default:  // 'm'
      // Bessel's formula
      deltax    = pos - x_inc_halves[5];
      dx2mDx2d4 = deltax * deltax - Dx2d4;
      // Forward routine
      // result    = vals_inc_halves[5] + ddx[3] * deltax + d2dx2[2] * dx2mDx2d4 / 2.0 +
      //            d3dx3[1] * deltax * dx2mDx2d4 / 6.0;
      vals_inc_halves_hat[5] = interpolate_hat;
      ddx_hat[3]             = interpolate_hat * deltax;
      d2dx2_hat[2]           = interpolate_hat * dx2mDx2d4 / 2.0;
      d3dx3_hat[1]           = interpolate_hat * deltax * dx2mDx2d4 / 6.0;
  }

  // Find third derivative
  for (l=0; l<3; l++)
  { // Forward routine
    // d3dx3[l] = recip_Dx * (d2dx2[l+2] - d2dx2[l]);
    d2dx2_hat[l+2] += d3dx3_hat[l] * recip_Dx;
    d2dx2_hat[l]   -= d3dx3_hat[l] * recip_Dx;
  }

  // Find second derivatives
  for (l=0; l<5; l++)
  { // Forward routine
    // d2dx2[l] = recip_Dx * (ddx[l+2] - ddx[l]);
    ddx_hat[l+2] += d2dx2_hat[l] * recip_Dx;
    ddx_hat[l]   -= d2dx2_hat[l] * recip_Dx;
  }

  // Find first derivatives
  for (l=0; l<7; l++)
  { // Forward routine
    // ddx[l] = recip_Dx * (vals_inc_halves[l+3] - vals_inc_halves[l+1]);
    vals_inc_halves_hat[l+3] += ddx_hat[l] * recip_Dx;
    vals_inc_halves_hat[l+1] -= ddx_hat[l] * recip_Dx;
  }

  // Fill the basic data at half points
  for (l=0; l<5; l++)
  { l2                    = 2*l;
    l2p1                  = l2 + 1;
    l2p2                  = l2 + 2;
    // Forward routine
    // vals_inc_halves[l2p1]   = 0.5 * (vals_inc_halves[l2] + vals_inc_halves[l2p2]);
    vals_inc_halves_hat[l2]   += 0.5 * vals_inc_halves_hat[l2p1];
    vals_inc_halves_hat[l2p2] += 0.5 * vals_inc_halves_hat[l2p1];
  }

  // Fill the basic data at full points
  for (l=0; l<6; l++)
  { l2                  = 2*l;
    // Forward routine
    // vals_inc_halves[l2] = vals[l];
    vals_hat[l] = vals_inc_halves_hat[l2];
  }

}



// -------------------------------------------------------------------------------
double interpolate_2d_cubic ( double x[6],
                              double y[6],
                              double vals[6][6], // field values (in)
                              double posx,
                              double posy,
                              int    parts )
{ // Bi-interpolate in 2d to cubic order.
  // Used in the semi-Lagrangian scheme
  int    xx;
  double int_val_with_x[6];

  // Interpolate in y (1d) for each x level
  for (xx=0; xx<6; xx++)
  { int_val_with_x[xx] = interpolate_1d_cubic (y, vals[xx], posy, parts);
  }

  // Now interpolate in x (1d)
  return interpolate_1d_cubic (x, int_val_with_x, posx, parts);
}




// -------------------------------------------------------------------------------
void interpolate_2d_cubic_adj ( double x[6],
                                double y[6],
                                double vals_hat[6][6],   // field values (out)
                                double interpolate_hat,  // interpolated value (in)
                                double posx,
                                double posy,
                                int    parts )
{ // Bi-interpolate in 2d to cubic order.
  // Used in the semi-Lagrangian scheme
  int    xx;
  double int_val_with_x_hat[6];

  // Now interpolate in x (1d)
  interpolate_1d_cubic_adj (x,
                            int_val_with_x_hat, // out
                            interpolate_hat,    // in
                            posx, parts);

  // Interpolate in y (1d) for each x level
  for (xx=0; xx<6; xx++)
  { interpolate_1d_cubic_adj (y,
                              vals_hat[xx],     // out
                              int_val_with_x_hat[xx],
                              posy, parts);
  }
}


// -------------------------------------------------------------------------------
double interpolate_3d_cubic ( double x[6],
                              double y[6],
                              double z[6],
                              double vals[6][6][6], // field values (in)
                              double posx,
                              double posy,
                              double posz,
                              int    parts )
{ // Tri-interpolate in 3d to cubic order.
  // Used in the semi-Lagrangian scheme
  int    xx;
  double int_val_with_x[6];

  // Interpolate in y  and z (1d) for each x level
  for (xx=0; xx<6; xx++)
  { int_val_with_x[xx] = interpolate_2d_cubic (y, z, vals[xx], posy, posz, parts);
  }

  // Now interpolate in x (1d)
  return interpolate_1d_cubic (x, int_val_with_x, posx, parts);
}


// -------------------------------------------------------------------------------
void interpolate_3d_cubic_adj ( double x[6],
                                double y[6],
                                double z[6],
                                double vals_hat[6][6][6], // field values (out)
                                double interpolate_hat,   // interpolated value (in)
                                double posx,
                                double posy,
                                double posz,
                                int    parts )
{ // Tri-interpolate in 3d to cubic order.
  // Used in the semi-Lagrangian scheme
  int    xx;
  double int_val_with_x_hat[6];

  // Now interpolate in x (1d)
  interpolate_1d_cubic_adj (x,
                            int_val_with_x_hat,  // out
                            interpolate_hat,     // in
                            posx, parts);

  // Interpolate in y  and z (1d) for each x level
  for (xx=0; xx<6; xx++)
  { interpolate_2d_cubic_adj (y, z,
                              vals_hat[xx],           // out
                              int_val_with_x_hat[xx], // in
                              posy, posz, parts);
  }

}


// -------------------------------------------------------------------------------
double interpolate_3d_cubic_t_linear ( double x[6],
                                       double y[6],
                                       double z[6],
                                       double t[2],
                                       double valst1[6][6][6],
                                       double valst2[6][6][6],
                                       double posx,
                                       double posy,
                                       double posz,
                                       double post,
                                       int    parts )
{ // Do linear interpolation in time and cubic interpolation in space
  // Written for the semi-Lagrangian scheme, but not used
  double result;
  double int_t1, int_t2;
  bool   lower, upper, middle;

  lower  = post == t[0];
  upper  = post == t[1];
  middle = !lower && !upper;

  if (lower || middle)
  { // Do interpolation in 3d space for time level 1
    int_t1 = interpolate_3d_cubic ( x, y, z,
                                    valst1,
                                    posx, posy, posz,
                                    parts );
  }

  if (upper || middle)
  { int_t2 = interpolate_3d_cubic ( x, y, z,
                                    valst2,
                                    posx, posy, posz,
                                    parts );
  }

  if (lower)
  { // post is at lower time
    result = int_t1;
  }
  else
  { if (upper)
    { // post is at upper time
      result = int_t2;
    }
    else
    { // post is in the middle, so need to do interpolation
      result = (int_t1 * (t[1] - post) + int_t2 * (post - t[0])) /
               (t[1] - t[0]);
    }
  }
  return result;
}



// -------------------------------------------------------------------------------
double interpolate_3d ( double x[2],
                        double y[2],
                        double z[2],
                        double vals[2][2][2],
                        double pos[3] )
{ // Interpolate in 3d (linear)
  // Used in the make winds scheme
  double vals_interpolated[2];
  // Interpolate in the first y/z plane
  // At the first x position
  vals_interpolated[0] = interpolate_2d (y, z, vals[0], &(pos[1]));
  // At the second x position
  vals_interpolated[1] = interpolate_2d (y, z, vals[1], &(pos[1]));
  // Interpolate in the x direction
  return interpolate_1d(x, vals_interpolated, pos[0]);
}

// -------------------------------------------------------------------------------
double interpolate_2d ( double x[2],
                        double y[2],
                        double vals[2][2],
                        double pos[2] )
{ // Interpolate in 2d (linear)
  // Used in the make winds scheme
  double vals_interpolated[2];
  // Interpolate in the y direction
  // At the first x position
  vals_interpolated[0] = interpolate_1d (y, vals[0], pos[1]);
  // At the second x position
  vals_interpolated[1] = interpolate_1d (y, vals[1], pos[1]);
  // Interpolate in the x direction
  return interpolate_1d(x, vals_interpolated, pos[0]);
}

// -------------------------------------------------------------------------------
double interpolate_1d ( double x[2],
                        double vals[2],
                        double pos )
{ // Interpolate in 1d (linear)
  // Used in the make winds scheme
  return vals[0] + (vals[1] - vals[0]) * (pos - x[0]) / (x[1] - x[0]);
}

