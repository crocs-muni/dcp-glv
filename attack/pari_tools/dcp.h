/* Implementation of dcp solver for 2-dimensional scalar decomposition.
  Implemented using these sources:
    - [1] https://sodilinux.itd.cnr.it/sdl6x2/documentazione/pari_gp/libpari.pdf
    - [2] https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html
*/
#ifndef DCP_H_INCLUDED
#define DCP_H_INCLUDED

#include <pari/pari.h>
#include <time.h>

/* DCP solver for the polynomial x1+x2+2
  Outputs the roots (x-coordinates) into a filename.
  See dcp_pari.py for python api. */
int dcp(GEN E, GEN p, long k, GEN map, GEN lam,GEN V, GEN f, char *filename);


int multidcp(GEN E, GEN p, long k, long l, GEN mapk, GEN mapl, GEN lam,GEN V, GEN f, char *filename);


/* Sums vector elements into result*/
void sum(GEN vector, long length, GEN result);

/* Multiplies vector elements into result*/
void prod(GEN vector, long length, GEN result);


/* DCP solver for the polynomial x1+x2+2 using the Semaev polynomial
  Outputs the roots (x-coordinates) into a filename.
  See dcp_pari.py for python api.

  The polynomial is S3(-x1-2,k1(x1),k2(beta*x1)) for guess=0 and
  S3(-beta^2*x1-2,k1(x1),k2(beta*x1)) otherwise
  */
int dcp_semaev(GEN E, GEN p, long k1, GEN map1, GEN lam, long k2, GEN map2, GEN V, GEN f,
                char *filename);


int multidcp_semaev(GEN E, GEN p, long l, GEN map0, long k1, GEN map1, GEN lam, long k2, GEN map2, GEN V, GEN f,
                char *filename);

#endif