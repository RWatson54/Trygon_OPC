inline void qts2fhs(const double *qs,
                    const double *jc,
                          double *fs)

{

  /* 
  // Fit the variables for computing the fluxes in the true space
  */ 

  double fts[2];

  fts[0] = qs[0] * 1.00000;
  fts[1] = qs[0] * 0.50000;

  /* 
  // Rescale the fluxes by the |J| J^{-1}, contained in "jc"
  */ 

  fs[0] = jc[0]*fts[0] + jc[2]*fts[1];
  fs[1] = jc[1]*fts[0] + jc[3]*fts[1];

}
