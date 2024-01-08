/* Implementation of dcp solver for 2-dimensional scalar decomposition and
  multiscalar multiplication. Implemented using these sources:
    - [1] https://sodilinux.itd.cnr.it/sdl6x2/documentazione/pari_gp/libpari.pdf
    - [2] https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html
*/

#include "dcp.h"

int main(int argc, char *argv[]) {
  /* p, a, b, beta, k1, k2, guess, filename are arguments
  p = base-field prime
  a,b = coefficients of the curve
  beta = action of the endomorphism phi on x-coordinate
  k1,k2 = decomposed parts of the scalar
  guess = guess whether ADD(Q,P) or ADD(Q,phi^2P) is computed
  */
  long stack;
  stack = 16000000000;
  pari_init(stack, 65536);

  GEN p, a, b, beta;
  long k1, k2, guess;
  char *filename;

  p = strtoi(argv[1]);
  a = strtoi(argv[2]);
  b = strtoi(argv[3]);
  beta = strtoi(argv[4]);
  k1 = strtol(argv[5], NULL, 10);
  k2 = strtol(argv[6], NULL, 10);
  guess = strtol(argv[7], NULL, 10);
  filename = argv[8];

  dcp_semaev(p, a, b, beta, k1, k2, guess, filename);
  pari_close();
  return 0;
}
