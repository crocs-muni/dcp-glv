/* Implementation of Baby-step-giant-step algorithms for DCP attacks.
  Implemented using these sources:
    - [1] https://sodilinux.itd.cnr.it/sdl6x2/documentazione/pari_gp/libpari.pdf
    - [2] https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html
*/

#include "bsgs.h"
#include <array>
#include <map>
#include <pari/pari.h>
#include <time.h>

int main(int argc, char *argv[]) {

  long stack = 16000000000;
  pari_init(stack, 65536);
  GEN p, a, b, xP, yP, xQ, yQ, lam;
  GEN P, Q, E;
  long diff, l1, l2;
  char *filename;

  p = strtoi(argv[1]);
  a = strtoi(argv[2]);
  b = strtoi(argv[3]);
  xP = strtoi(argv[4]);
  yP = strtoi(argv[5]);
  xQ = strtoi(argv[6]);
  yQ = strtoi(argv[7]);
  lam = strtoi(argv[8]);
  l1 = strtol(argv[9], NULL, 10);
  l2 = strtol(argv[10], NULL, 10);
  E = ellinit(mkvec2(a, b), p, 0);
  P = ellmul(E, mkvec2(xP, yP), stoi(1));
  Q = ellmul(E, mkvec2(xQ, yQ), stoi(1));
  filename = argv[11];

  diff = l1 - l2;

  if (diff == 0 || diff == 1 || diff == -1) {
    bsgs_simple(E, P, Q, lam, l1, l2, filename);
  } else {
    bsgs_balanced(E, P, Q, lam, l1, l2, diff, filename);
  }
  pari_close();
  return 0;
}
