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

typedef std::array<long, 4> longgen;

void gen_to_longs(GEN x, longgen *longs) {
  /* Converts t_INT into an array of indiviual limbs. See page 21 of [1].*/
  for (int i = 0; i <= lgefint(x) - 3; i++) {
    (*longs)[i] = *int_W(x, i);
  }
}

void longgen_print(GEN x, longgen *longs) {
  /* Prints result of @gen_to_longs.*/
  for (int i = 0; i <= lgefint(x) - 3; i++) {
    printf("%ld ", (*longs)[i]);
  }
  printf("\n");
}

long *copy_INT(GEN X, GEN Y) {

  if ((cmpii(X, stoi(0)) == 0) || LONG_MAX == (*int_MSW(Y))) {
    return gcopy(Y);
  }
  gel(X, 1) = gel(Y, 1);
  gel(X, 0) = gel(Y, 0);
  for (int i = 0; i <= lgefint(Y) - 3; i++) {
    *((long *)int_W(X, i)) = *((long *)int_W(Y, i));
  }
  return X;
}

void copy_INTMOD(GEN X, GEN Y) {

  GEN X1, Y1, X2, Y2;
  X1 = gel(X, 1);
  Y1 = gel(Y, 1);
  gel(X, 1) = copy_INT(X1, Y1);
  X2 = gel(X, 2);
  Y2 = gel(Y, 2);
  gel(X, 2) = copy_INT(X2, Y2);
}

void copy_ellpoint(GEN X, GEN Y) {
  GEN X1, Y1, X2, Y2;
  X1 = gel(X, 1);
  Y1 = gel(Y, 1);

  copy_INTMOD(X1, Y1);
  X2 = gel(X, 2);
  Y2 = gel(Y, 2);
  copy_INTMOD(X2, Y2);
}

long *add_inplace(GEN x, GEN y) {

  pari_sp ltop = avma;
  GEN z = gadd(x, y);
  pari_sp lbot = avma;
  x = copy_INT(x, z);
  return gerepile(ltop, lbot, x);
}

long *elladd_inplace(GEN e, GEN P, GEN Q) {
  if (glength(P) == 1) {
    return gcopy(Q);
  }
  pari_sp ltop = avma;
  GEN tmp = elladd(e, P, Q);
  pari_sp lbot = avma;
  copy_ellpoint(P, tmp);
  return gerepile(ltop, lbot, P);
}

long *ellsub_inplace(GEN e, GEN P, GEN Q) {
  if (glength(Q) == 1) {
    return gcopy(P);
  }
  if (glength(P) == 1) {
    return ellneg(e, Q);
  }
  pari_sp ltop = avma;
  GEN tmp = ellsub(e, P, Q);
  pari_sp lbot = avma;
  copy_ellpoint(P, tmp);
  return gerepile(ltop, lbot, P);
}

void bsgs_simple(GEN E, GEN P, GEN Q, GEN lam, long l1, long l2,
                 char *filename) {
  /* Simple Baby-step-giant-step on k = k1+@lam*k2 with |k1|<@l1 bits, |k1|<@l2
  bits.
  @P,@Q are points on @E s.t. @Q=k@P.
  Prints k1\nk2\n into filename.*/
  GEN N1, N2, lamP;
  std::map<longgen, GEN> el;
  std::map<longgen, GEN>::iterator it;
  GEN one = stoi(1), mone = stoi(-1), two = stoi(2), zero = stoi(0);

  lamP = ellmul(E, P, lam);
  N1 = gpowgs(two, l1);
  N2 = gpowgs(two, l2);

  // baby steps
  GEN i = zero, j = zero;
  GEN X, R = ellinf();
  longgen Xlongs;
  while (cmpii(N2, i) != -1) {
    Xlongs.fill(0);
    if (glength(R) == 1)
      X = mone;
    else
      X = gel(gel(R, 1), 2);
    gen_to_longs(X, &Xlongs);
    el[Xlongs] = gcopy(i);
    i = add_inplace(i, one);
    R = elladd_inplace(E, R, lamP);
  }

  // giantsteps
  FILE *file;
  GEN R1, X2, s, Q3 = gcopy(Q);
  longgen X2longs;
  X2longs.fill(0);
  while (cmpii(N1, j) != -1) {
    X2 = gel(gel(Q3, 1), 2);
    gen_to_longs(X2, &X2longs);
    it = el.find(X2longs);
    if (it != el.end()) {
      // check for sign
      R1 = ellmul(E, lamP, it->second);
      s = it->second;
      if (cmpii(gel(gel(Q3, 2), 2), gel(gel(R1, 2), 2)) != 0)
        s = gneg(it->second);
      // save the results to a file
      file = fopen(filename, "w");
      pari_fprintf(file, "\%Ps\n\%Ps\n", j, s);
      fclose(file);
      return;
    }
    j = add_inplace(j, one);
    Q3 = ellsub_inplace(E, Q3, P);
  }
}

void bsgs_balanced(GEN E, GEN P, GEN Q, GEN lam, long l1, long l2, long diff,
                   char *filename) {
  /* Balanced Baby-step-giant-step on k = k1+@lam*k2 with |k1|<@l1 bits,
  |k1|<@l2 bits. The expression k = k1+@lam*k2 is modified so that |k1|=|k2| for
  more efficiency.
  @P,@Q are points on @E s.t. @Q=k@P.
  Prints k1\nk2\n into filename.*/

  long big, shared, small;
  std::map<longgen, GEN> el;
  std::map<longgen, GEN>::iterator it;
  GEN Big, Shared, Small, N1, lamP, P2m, boundi, boundj;
  GEN i = stoi(0), j;
  GEN X, X2, T, R = ellinf(), S = ellinf();
  GEN two = stoi(2), zero = stoi(0), mone = stoi(-1), one = stoi(1),
      inf = ellinf();

  FILE *file;
  GEN Q3 = gcopy(Q);
  GEN tmp, e1, e2;
  GEN R1, R2, s, result1, result2;
  longgen X2longs;
  X2longs.fill(0);

  // Decide which scalar is larger
  if (diff < 0) {
    big = l2;
    small = l1;
    shared = -diff / 2;
  } else {
    big = l1;
    small = l2;
    shared = diff / 2;
  }
  boundi = gpowgs(two, shared);
  boundj = gpowgs(two, small);

  // pari_printf("Balancing out to complexity of bsgs: %ld,
  // %ld\n",big-shared,shared+small);

  // Shift variables
  N1 = gpowgs(two, big - shared);
  lamP = ellmul(E, P, lam);
  P2m = ellmul(E, P, N1);
  if (diff < 0) {
    Small = P;
    Shared = ellmul(E, P2m, lam);
    Big = lamP;
  } else {
    Small = lamP;
    Shared = P2m;
    Big = P;
  }
  // Baby steps
  longgen Xlongs;
  T = ellinf();
  while (cmpii(boundi, i) == 1) {
    S = inf;
    j = zero;
    // copy_INT(j,zero);
    while (cmpii(boundj, j) == 1) {
      T = ellinf();
      T = elladd_inplace(E, T, S);
      T = elladd_inplace(E, T, R);

      Xlongs.fill(0);
      if (glength(T) == 1)
        X = mone;
      else
        X = gel(gel(T, 1), 2);
      gen_to_longs(X, &Xlongs);

      el[Xlongs] = gcopy(mkvec2(i, j));

      j = add_inplace(j, one);
      S = elladd_inplace(E, S, Small);
    }
    i = add_inplace(i, one);
    R = elladd_inplace(E, R, Shared);
  }

  // Giant steps
  j = zero;
  while (cmpii(N1, j) == 1) {
    X2 = gel(gel(Q3, 1), 2);
    gen_to_longs(X2, &X2longs);
    it = el.find(X2longs);
    if (it != el.end()) {
      tmp = it->second;
      e1 = gel(tmp, 1);
      e2 = gel(tmp, 2);
      // check for sign
      R1 = elladd(E, ellmul(E, Shared, e1), ellmul(E, Small, e2));
      R2 = ellsub(E, Q, ellmul(E, Big, j));
      s = j;
      if (cmpii(gel(gel(R1, 2), 2), gel(gel(R2, 2), 2)) != 0)
        s = gneg(j);
      result1 = gadd(s, gmul(N1, e1));
      result2 = e2;
      if (diff < 0) {
        result1 = e2;
        result2 = gadd(gmul(e1, N1), s);
      }
      // save to a file
      file = fopen(filename, "w");
      pari_fprintf(file, "\%Ps\n\%Ps\n", result1, result2);
      fclose(file);
      return;
    }
    j = add_inplace(j, one);

    Q3 = ellsub_inplace(E, Q3, Big);
  }
}
