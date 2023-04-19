inline void qtf2itf(const double **qf,
                    const double  *ni,
                          double **fin)

{

  double eps2, A, fiL;

  eps2 = 1.0000;

  A = abs(1.0000*ni[0]) + abs(0.5000*ni[1]);

  fiL = 0.500 * ( 1.0000 * qf[0][0] + 1.0000 * qf[1][0] ) * ni[0] +
        0.500 * ( 0.5000 * qf[0][0] + 0.5000 * qf[1][0] ) * ni[1] -
        0.500 * eps2 * A * ( qf[1][0] - qf[0][0] );

  fin[0][0] =  fiL;
  fin[1][0] = -fiL;

}
