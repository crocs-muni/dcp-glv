/* Implementation of dcp solver for 2-dimensional scalar decomposition.
  Implemented using these sources:
    - [1] https://sodilinux.itd.cnr.it/sdl6x2/documentazione/pari_gp/libpari.pdf
    - [2] https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html
*/

#include "dcp.h"
#include <pari/pari.h>
#include <time.h>

long *create_dcp_polynomial(GEN map, GEN V){

  GEN num, den, dcp_poly;
  num = gel(map, 1);
  den = gel(map, 2);  


  GEN monom, x, xx,n,d;
  dcp_poly = stoi(0);
  x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
  for(int i=1;i<lg(V);i++){
    monom = gel(gel(V,i),1);
    xx = gpow(x,gel(gel(V,i),2),-1);
    n = gpow(num,gel(gel(V,i),3),-1);
    d = gpow(den,gel(gel(V,i),4),-1);
    monom = gmul(monom,xx);
    monom = gmul(monom,n);
    monom = gmul(monom,d);
    dcp_poly = gadd(dcp_poly,monom);

  }

  return dcp_poly;

}

int dcp(GEN E, GEN p, long k, GEN map, GEN lam,GEN V, GEN f, char *filename) {

  /* DCP solver for the polynomial x1+x2+2
  Outputs the roots (x-coordinates) into a filename.
  See dcp_pari.py for python api. */

  GEN den, num, roots, x, Q, g, dcp_poly;
  FILE *file;


  // pari_printf("%Ps\n",num);

  dcp_poly = create_dcp_polynomial(map, V);

  roots = polrootsmod(dcp_poly, p);
  file = fopen(filename, "w");
  GEN ys;
  // printf("k:%ld 0\n",k);
  // pari_printf("%Ps\n", V);
  // pari_printf("roots = %Ps\n", roots);

  for(int i=1;i<lg(roots);i++){
    ys = ellordinate(E,gel(roots,i),0);
    for(int j=1;j<lg(ys);j++){
      // pari_printf("%Ps %Ps\n", gel(roots,i), gel(ys,j));
      Q = ellmul(E, mkvec2(gel(roots,i), gel(ys,j)), stoi(k));
      Q = ellmul(E, Q, lam);
      g = gsubst(f,varn(gp_read_str("X1")),gel(roots,i));
      g = gsubst(g,varn(gp_read_str("Y1")),gel(ys,j));
      g = gsubst(g,varn(gp_read_str("X2")),gel(Q,1));
      g = gsubst(g,varn(gp_read_str("Y2")),gel(Q,2));
      if (isintzero(g[2])) {
        // printf("Yes\n");
        pari_fprintf(file, "%Ps, %Ps", gel(roots,i), gel(ys,j));
        fclose(file);
        return 1;

      }
        // printf("No\n");

    }
  }
  fclose(file);
  return 0;
}


long *create_multidcp_polynomial(GEN mapk, GEN mapl, GEN V){

  GEN den, num, den0, num0, roots, dcp_poly;

  num = gel(mapk, 1);
  den = gel(mapk, 2);
  num0 = gel(mapl, 1);
  den0 = gel(mapl, 2);

  GEN monom, x, xx,n,d, n0,d0;
  dcp_poly = stoi(0);
  x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
  for(int i=1;i<lg(V);i++){
    monom = gel(gel(V,i),1);
    n0 = gpow(num0,gel(gel(V,i),2),-1);
    d0 = gpow(den0,gel(gel(V,i),3),-1);
    n = gpow(num,gel(gel(V,i),4),-1);
    d = gpow(den,gel(gel(V,i),5),-1);
    

    monom = gmul(monom,n);
    monom = gmul(monom,d);
    monom = gmul(monom,n0);
    monom = gmul(monom,d0);
    dcp_poly = gadd(dcp_poly,monom);

  }

  return dcp_poly;

}

int multidcp(GEN E, GEN p, long k, long l, GEN mapk, GEN mapl, GEN lam,GEN V, GEN f, char *filename) {

  /* DCP solver for the polynomial x1+x2+2
  Outputs the roots (x-coordinates) into a filename.
  See dcp_pari.py for python api. */

  GEN den, num, den0, num0, roots, x, Q, g;
  FILE *file;

  num = gel(mapk, 1);
  den = gel(mapk, 2);
  num0 = gel(mapl, 1);
  den0 = gel(mapl, 2);

  GEN dcp_poly;
  dcp_poly = create_multidcp_polynomial(mapk,mapl,V);
  
  roots = polrootsmod(dcp_poly, p);
  file = fopen(filename, "w");
  GEN ys;
  for(int i=1;i<lg(roots);i++){
    ys = ellordinate(E,gel(roots,i),0);
    for(int j=1;j<lg(ys);j++){
      Q = ellmul(E, mkvec2(gel(roots,i), gel(ys,j)), stoi(l));
      g = gsubst(f,varn(gp_read_str("X1")),gel(Q,1));
      g = gsubst(g,varn(gp_read_str("Y1")),gel(Q,2));
      Q = ellmul(E, mkvec2(gel(roots,i), gel(ys,j)), stoi(k));
      Q = ellmul(E, Q, lam);
      g = gsubst(g,varn(gp_read_str("X2")),gel(Q,1));
      g = gsubst(g,varn(gp_read_str("Y2")),gel(Q,2));
      if (isintzero(g[2])) {
        pari_fprintf(file, "%Ps, %Ps", gel(roots,i), gel(ys,j));
        fclose(file);
        return 1;

      }
    }
  }
  fclose(file);
  return 0;
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



long *create_glvdcp_polynomial(GEN map1, GEN map2, GEN V){

  GEN u1, v1, u2, v2, dcp_poly;
  u1 = gel(map1, 1);
  v1 = gel(map1, 2);
  u2 = gel(map2, 1);
  v2 = gel(map2, 2);

  GEN monom, x, xx,n1,d1, n2,d2;
  dcp_poly = stoi(0);
  x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
  for(int i=1;i<lg(V);i++){
    monom = gel(gel(V,i),1);
    xx = gpow(x,gel(gel(V,i),2),-1);
    n1 = gpow(u1,gel(gel(V,i),3),-1);
    d1 = gpow(v1,gel(gel(V,i),4),-1);
    n2 = gpow(u2,gel(gel(V,i),5),-1);
    d2 = gpow(v2,gel(gel(V,i),6),-1);

    monom = gmul(monom,xx);
    monom = gmul(monom,n1);
    monom = gmul(monom,d1);
    monom = gmul(monom,n2);
    monom = gmul(monom,d2);
    dcp_poly = gadd(dcp_poly,monom);

  }

  return dcp_poly;

}


int dcp_semaev(GEN E, GEN p, long k1, GEN map1, GEN lam, long k2, GEN map2, GEN V, GEN f,
                char *filename) {


  GEN roots;
  FILE *file;


  GEN P1,P2;
  
  GEN dcp_poly,g;
  dcp_poly = create_glvdcp_polynomial(map1,map2,V);



  file = fopen(filename, "w");
  roots = polrootsmod(dcp_poly, p);
  GEN ys;


  for(int i=1;i<lg(roots);i++){
    ys = ellordinate(E,gel(roots,i),0);
    for(int j=1;j<lg(ys);j++){
      // pari_printf("%Ps %Ps\n", gel(roots,i), gel(ys,j));
      P1 = ellmul(E, mkvec2(gel(roots,i), gel(ys,j)), stoi(1));
      P2 = ellmul(E, P1, lam);
      P1 = ellmul(E,P1,stoi(k1));
      P2= ellmul(E,P2,stoi(k2));
      P2 = elladd(E,P1,P2);
      g = gsubst(f,varn(gp_read_str("X1")),gel(roots,i));
      g = gsubst(g,varn(gp_read_str("Y1")),gel(ys,j));
      g = gsubst(g,varn(gp_read_str("X2")),gel(P2,1));
      g = gsubst(g,varn(gp_read_str("Y2")),gel(P2,2));
      if (isintzero(g[2])) {
        // printf("Yes\n");
        pari_fprintf(file, "%Ps, %Ps", gel(roots,i), gel(ys,j));
        fclose(file);
        return 1;

      }
      // printf("No\n");
    }
  }
  fclose(file);
  return 0;
}


long *create_glvmultidcp_polynomial(GEN map0, GEN map1, GEN map2, GEN V){

  GEN u1, v1, u2, v2, g, u0,v0;
  u1 = gel(map1, 1);
  v1 = gel(map1, 2);
  u2 = gel(map2, 1);
  v2 = gel(map2, 2);
  u0 = gel(map0, 1);
  v0 = gel(map0, 2);


  GEN monom, x, xx,n1,d1, n2,d2, n0,d0, dcp_poly;
  dcp_poly = stoi(0);
  x = gtopolyrev(mkvec2(stoi(0), stoi(1)), -1);
  for(int i=1;i<lg(V);i++){
    monom = gel(gel(V,i),1);
    n0 = gpow(u0,gel(gel(V,i),2),-1);
    d0 = gpow(v0,gel(gel(V,i),3),-1);

    n1 = gpow(u1,gel(gel(V,i),4),-1);
    d1 = gpow(v1,gel(gel(V,i),5),-1);

    n2 = gpow(u2,gel(gel(V,i),6),-1);

    d2 = gpow(v2,gel(gel(V,i),7),-1);
  

    monom = gmul(monom,n1);
    monom = gmul(monom,d1);
    monom = gmul(monom,n2);
    monom = gmul(monom,d2);
    monom = gmul(monom,n0);
    monom = gmul(monom,d0);

    dcp_poly = gadd(dcp_poly,monom);

  }

  return dcp_poly;

}


int multidcp_semaev(GEN E, GEN p, long l, GEN map0, long k1, GEN map1, GEN lam, long k2, GEN map2, GEN V, GEN f,
                char *filename) {


  GEN roots;
  FILE *file;


  GEN P1,P2;
  
  GEN dcp_poly, g;

  dcp_poly = create_glvmultidcp_polynomial(map0,map1,map2,V);


  file = fopen(filename, "w");
  roots = polrootsmod(dcp_poly, p);
  GEN ys;
  for(int i=1;i<lg(roots);i++){
    ys = ellordinate(E,gel(roots,i),0);
    for(int j=1;j<lg(ys);j++){
      P1 = ellmul(E, mkvec2(gel(roots,i), gel(ys,j)), stoi(l));
      g = gsubst(f,varn(gp_read_str("X1")),gel(P1,1));
      g = gsubst(g,varn(gp_read_str("Y1")),gel(P1,2));
      P1 = ellmul(E, mkvec2(gel(roots,i), gel(ys,j)), stoi(1));
      P2 = ellmul(E, P1, lam);
      P1 = ellmul(E,P1,stoi(k1));
      P2= ellmul(E,P2,stoi(k2));
      P2 = elladd(E,P1,P2);
      g = gsubst(g,varn(gp_read_str("X2")),gel(P2,1));
      g = gsubst(g,varn(gp_read_str("Y2")),gel(P2,2));
      if (isintzero(g[2])) {
        pari_fprintf(file, "%Ps, %Ps", gel(roots,i), gel(ys,j));
        fclose(file);
        return 1;

      }
    }
  }
  fclose(file);
  return 0;
}
