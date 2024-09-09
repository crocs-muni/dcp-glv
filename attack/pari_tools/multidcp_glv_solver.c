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

  GEN p, a, b, Vs, fs, V,f;
  GEN lam;
  long k1, k2, l;
  char *filename, *Vpoly_path, *fpoly_path;

  p = strtoi(argv[1]);
  a = strtoi(argv[2]);
  b = strtoi(argv[3]);
  l = strtol(argv[4], NULL, 10);
  k1 = strtol(argv[5], NULL, 10);
  lam = strtoi(argv[6]);
  k2 = strtol(argv[7], NULL, 10);
  Vpoly_path = argv[8];
  fpoly_path = argv[9];
  filename = argv[10];
  FILE *fptr;
  fptr = fopen(Vpoly_path,"r");

  if(fptr == NULL)
  {
  printf("Error!");   
  exit(1);             
  }
  Vs = gp_read_stream(fptr);
  fclose(fptr);

  fptr = fopen(fpoly_path,"r");

  if(fptr == NULL)
  {
  printf("Errr!");   
  exit(1);             
  }
  fs = gp_read_stream(fptr);
  fclose(fptr);

  GEN E, map1,map2,map0;
  E = ellinit(mkvec2(a, b), p, 0);
  map1 = ellxn(E, k1, -1);
  map2 = ellxn(E, k2, -1);
  map0 = ellxn(E, l, -1);
  pari_sp av;

  for(int i=1;i<lg(Vs);i++){
    V = gel(Vs,i);
    f = gel(fs,i);
    f = gsubst(f,varn(gp_read_str("A")),a);
    f = gsubst(f,varn(gp_read_str("B")),b);
    av = avma;
    if(multidcp_semaev(E,p,l,map0, k1, map1, lam, k2, map2, V, f, filename)) break;
    avma = av;
  }
  pari_close();
  return 0;
}
