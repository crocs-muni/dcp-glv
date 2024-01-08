/* Implementation of dcp solver for 2-dimensional scalar decomposition.
  Implemented using these sources:
    - [1] https://sodilinux.itd.cnr.it/sdl6x2/documentazione/pari_gp/libpari.pdf
    - [2] https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html
*/

#include "dcp.h"
#include <pari/pari.h>
#include <time.h>

void dcp(GEN p, GEN a, GEN b, long k, char *filename) {

  /* DCP solver for the polynomial x1+x2+2
  Outputs the roots (x-coordinates) into a filename.
  See dcp_pari.py for python api. */

  GEN E, den, num, f, roots, map, x;
  FILE *file;

  /* This computes multiplication map from division polynomial. Probably
  obsolete as there is ellxn GEN E,
  den,num,phik,phikm1,phikp1,phim1,x,phikm1_0,phikp1_0,f,roots; E =
  ellinit(mkvec2(a,b),p,0); phikm1 = elldivpol(E,k-1,-1); phikp1 =
  elldivpol(E,k+1,-1); phik =elldivpol(E,k,-1);

  phim1 =
  gadd(gmul(stoi(4),gmul(x,gmul(x,x))),gadd(gmul(x,gmul(stoi(4),a)),gmul(stoi(4),b)));
  //4x^3+4ax+4b den = gsqr(phik); phikm1_0 = gdiv(phikm1,phim1); phikp1_0 =
  gdiv(phikp1,phim1); num =
  gsub(gmul(x,den),gmul(phim1,gmul(phikm1_0,phikp1_0)));
  */
  E = ellinit(mkvec2(a, b), p, 0);
  map = ellxn(E, k, -1);
  num = gel(map, 1);
  den = gel(map, 2);
  x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
  f = gadd(gadd(gmul(x, den), num), gmul(stoi(+2), den));
  roots = polrootsmod(f, p);
  file = fopen(filename, "w");
  pari_fprintf(file, "%Ps\n", roots);
  fclose(file);
}

void dcp_extension(GEN p, GEN p1, GEN p2, GEN a1, GEN a2, GEN b1, GEN b2,
                   long k, char *filename) {
  /* DCP solver for the polynomial x1+x2+2 on a curve over quadratic extension.
  Outputs the roots (x-coordinates) into a filename.
  See dcp_pari.py for python api. */

  GEN E, den, num, map, f, roots, z1, T, g, x;
  GEN a, b;
  FILE *file;

  z1 = ffinit(p, 1, -1);
  // z = gtopolyrev(mkvec2(stoi(0),z1),-1);

  T = gadd(gadd(gmul(z1, z1), gmul(p1, z1)), p2);
  g = ffgen(T, -1);
  a = gadd(gmul(a1, g), a2);
  b = gadd(gmul(b1, g), b2);
  E = ellinit(mkvec2(a, b), g, 0);

  /* This computes multiplication map from division polynomial. Probably
  obsolete as there is ellxn GEN E,
  den,num,phik,phikm1,phikp1,phim1,x,phikm1_0,phikp1_0,f,roots,z1, T,g; phikm1 =
  elldivpol(E,k-1,-1); phikp1 = elldivpol(E,k+1,-1); phik =elldivpol(E,k,-1); x
  = gtopolyrev(mkvec2(stoi(0),stoi(1)),-1); phim1 =
  gadd(gmul(stoi(4),gmul(x,gmul(x,x))),gadd(gmul(x,gmul(stoi(4),a)),gmul(stoi(4),b)));
  //4x^3+4ax+4b den = gsqr(phik); phikm1_0 = gdiv(phikm1,phim1); phikp1_0 =
  gdiv(phikp1,phim1); num =
  gsub(gmul(x,den),gmul(phim1,gmul(phikm1_0,phikp1_0)));
  */
  x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
  map = ellxn(E, k, -1);
  num = gel(map, 1);
  den = gel(map, 2);
  f = gadd(gadd(gmul(x, den), num), gmul(stoi(2), den));
  roots = polrootsmod(f, g);
  file = fopen(filename, "w");
  pari_fprintf(file, "%Ps\n", roots);
  fclose(file);
}

void sum(GEN vector, long length, GEN result) {
  /* Sums vector elements into result*/
  result = stoi(0);
  for (int i = 0; i < length; i++) {
    result = gadd(gel(vector, i), result);
  }
}

void prod(GEN vector, long length, GEN result) {
  /* Multiplies vector elements into result*/
  result = stoi(1);
  for (int i = 0; i < length; i++) {
    result = gmul(gel(vector, i), result);
  }
}

void dcp_semaev(GEN p, GEN a, GEN b, GEN beta, long k1, long k2, long guess,
                char *filename) {

  /* DCP solver for the polynomial x1+x2+2 using the Semaev polynomial
  Outputs the roots (x-coordinates) into a filename.
  See dcp_pari.py for python api.

  The polynomial is S3(-x1-2,k1(x1),k2(beta*x1)) for guess=0 and
  S3(-beta^2*x1-2,k1(x1),k2(beta*x1)) otherwise
  */

  GEN E, den, num, f, roots, map, x;
  FILE *file;

  /* This computes multiplication map from division polynomial. Probably
  obsolete as there is ellxn. Could be relevant for ellyn (which is not in
  pari).

  GEN E, den,num,phik,phikm1,phikp1,phim1,x,phikm1_0,phikp1_0,f,roots;
  E = ellinit(mkvec2(a,b),p,0);
  phikm1 = elldivpol(E,k-1,-1);
  phikp1 = elldivpol(E,k+1,-1);
  phik =elldivpol(E,k,-1);

  phim1 =
  gadd(gmul(stoi(4),gmul(x,gmul(x,x))),gadd(gmul(x,gmul(stoi(4),a)),gmul(stoi(4),b)));
  //4x^3+4ax+4b den = gsqr(phik); phikm1_0 = gdiv(phikm1,phim1); phikp1_0 =
  gdiv(phikp1,phim1); num =
  gsub(gmul(x,den),gmul(phim1,gmul(phikm1_0,phikp1_0)));
  */
  E = ellinit(mkvec2(a, b), p, 0);
  GEN beta2 = gsqr(beta);

  if (k2 == 0) {
    map = ellxn(E, k1, -1);
    num = gel(map, 1);
    den = gel(map, 2);
    x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
    if (guess == 0)
      f = gadd(gadd(gmul(x, den), num), gmul(stoi(2), den));
    else
      f = gadd(gadd(gmul(gmul(beta2, x), den), num), gmul(stoi(2), den));
    roots = polrootsmod(f, p);
    file = fopen(filename, "w");
    pari_fprintf(file, "%Ps\n", roots);
    fclose(file);
  } else {
    GEN map1, map2, g, u1, v1, u2, v2;
    map1 = ellxn(E, k1, -1);
    map2 = ellxn(E, k2, -1);

    x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
    if (guess == 0)
      g = gsub(stoi(-2), x);
    else
      g = gsub(stoi(-2), gmul(beta2, x));
    u1 = gel(map1, 1);
    v1 = gel(map1, 2);
    u2 = gel(map2, 1);
    v2 = gel(map2, 2);
    u2 = gmul(beta, u2);

    // f = u2**2*v1**2*g**2 - 2*u1*u2*v1*v2*g**2 + u1**2*v2**2*g**2 -
    // 2*u2*v1**2*v2*g*A - 2*u1*v1*v2**2*g*A + v1**2*v2**2*A**2 -
    // 4*v1**2*v2**2*g*B - 2*u1*u2**2*v1*g - 2*u1**2*u2*v2*g - 2*u1*u2*v1*v2*A -
    // 4*u2*v1**2*v2*B - 4*u1*v1*v2**2*B + u1**2*u2**2 f =
    // f0-2*f1+f2-2*f3-2*f4+f5-4*f6-2*f7-2*f8-2*f9-4*f10-4*f11+f12
    GEN f0 = stoi(0), f1 = stoi(1), f2, f3, f4, f5, f6, f7, f8, f9, f10, f11,
        f12;
    // f0: u2^2*v1^2*g^2
    f0 = gsqr(gmul(u2, gmul(v1, g)));
    // prod(mkvecn(3,u2,v1,g),3,f0);
    // f0 = gsqr(f0);

    // u1*u2*v1*v2*g^2
    f1 = gmul(gmul(u1, gmul(u2, gmul(v1, v2))), gsqr(g));
    // u1^2*v2^2*g^2
    f2 = gsqr(gmul(u1, gmul(v2, g)));
    // u2*v1^2*v2*g*A
    f3 = gmul(u2, gmul(gsqr(v1), gmul(v2, gmul(g, a))));
    // u1*v1*v2^2*g*A
    f4 = gmul(u1, gmul(gsqr(v2), gmul(v1, gmul(g, a))));
    // v1^2*v2^2*A^2
    f5 = gsqr(gmul(v1, gmul(v2, a)));
    // v1^2*v2^2*g*B
    f6 = gmul(gsqr(gmul(v1, v2)), gmul(g, b));
    // u1*u2^2*v1*g
    f7 = gmul(u1, gmul(gsqr(u2), gmul(v1, g)));
    // u1^2*u2*v2*g
    f8 = gmul(u2, gmul(gsqr(u1), gmul(v2, g)));
    // u1*u2*v1*v2*A
    f9 = gmul(u1, gmul(u2, gmul(v1, gmul(v2, a))));
    // u2*v1^2*v2*B
    f10 = gmul(u2, gmul(gsqr(v1), gmul(v2, b)));
    // u1*v1*v2^2*B
    f11 = gmul(u1, gmul(gsqr(v2), gmul(v1, b)));
    // u1^2*u2^2
    f12 = gsqr(gmul(u1, u2));
    // f = f0-2*f1+f2-2*f3-2*f4+f5-4*f6-2*f7-2*f8-2*f9-4*f10-4*f11+f12
    GEN coefs = mkvecn(13, stoi(1), stoi(-2), stoi(1), stoi(-2), stoi(-2),
                       stoi(1), stoi(-4), stoi(-2), stoi(-2), stoi(-2),
                       stoi(-4), stoi(-4), stoi(1));
    GEN monoms =
        mkvecn(13, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12);
    f = stoi(0);
    GEN mon;
    for (int i = 1; i <= 13; i++) {

      mon = gmul(gel(monoms, (long)i), gel(coefs, (long)i));
      f = gadd(f, mon);
    }
    roots = polrootsmod(f, p);
    file = fopen(filename, "w");
    pari_fprintf(file, "%Ps\n", roots);
    fclose(file);
  }
}