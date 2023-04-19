inline void cts2qts(const double *cs,
                    const double *js,
                    const double *q0,
                          double *qs)

{
  
  double dt;

  dt = 0.0001;

  /*
  // Advance the calculation in time
  */
  qs[0] = q0[0] - dt*cs[0]/js[0];

}
