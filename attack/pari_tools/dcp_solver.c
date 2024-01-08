/* Implementation of dcp solver for 2-dimensional scalar decomposition.
  Implemented using these sources:
    - [1] https://sodilinux.itd.cnr.it/sdl6x2/documentazione/pari_gp/libpari.pdf
    - [2] https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html
*/

#include "dcp.h"

int main(int argc, char *argv[]) {

  long stack;
  stack = 16000000000;
  pari_init(stack, 65536);

  GEN p, a, b;
  long k;
  char *filename;

  p = strtoi(argv[1]);
  k = strtol(argv[2], NULL, 10);
  a = strtoi(argv[3]);
  b = strtoi(argv[4]);
  filename = argv[5];

  dcp(p, a, b, k, filename);
  pari_close();
  return 0;
}
