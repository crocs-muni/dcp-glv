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

  GEN p, a, b, Vs, lam, fs, V,f;
  long k;
  char *filename, *Vpoly_path, *fpoly_path;

  p = strtoi(argv[1]);
  k = strtol(argv[2], NULL, 10);
  lam = strtoi(argv[3]);
  a = strtoi(argv[4]);
  b = strtoi(argv[5]);
  Vpoly_path = argv[6];
  fpoly_path = argv[7];
  filename = argv[8];
  
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
  printf("Error!");   
  exit(1);             
  }
  fs = gp_read_stream(fptr);
  fclose(fptr);

  GEN E, map;
  E = ellinit(mkvec2(a, b), p, 0);
  map = ellxn(E, k, -1);
  pari_sp av;
  for(int i=1;i<lg(Vs);i++){
    V = gel(Vs,i);
    f = gel(fs,i);
    f = gsubst(f,varn(gp_read_str("A")),a);
    f = gsubst(f,varn(gp_read_str("B")),b);
    av = avma;
    // pari_printf("%Ps,\n %Ps\n",V,f);
    if(dcp(E,p, k, map, lam,V,f,filename)) break;
    avma = av;
  }
  pari_close();
  return 0;
}
