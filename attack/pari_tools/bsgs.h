/* Implementation of Baby-step-giant-step algorithms for DCP attacks.
  Implemented using these sources:
    - [1] https://sodilinux.itd.cnr.it/sdl6x2/documentazione/pari_gp/libpari.pdf
    - [2] https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html
*/

#ifndef BSGS_H_INCLUDED
#define BSGS_H_INCLUDED

#include <array>
#include <map>
#include <pari/pari.h>
#include <time.h>

typedef std::array<long, 4> longgen;

/* Converts t_INT into an array of indiviual limbs. See page 21 of [1].*/
void gen_to_longs(GEN x, longgen *longs);

/* Prints result of @gen_to_longs.*/
void longgen_print(GEN x, longgen *longs);

long *copy_INT(GEN X, GEN Y);

void copy_INTMOD(GEN X, GEN Y);

void copy_ellpoint(GEN X, GEN Y);

long *add_inplace(GEN x, GEN y);

long *elladd_inplace(GEN e, GEN P, GEN Q);

long *ellsub_inplace(GEN e, GEN P, GEN Q);

/* Simple Baby-step-giant-step on k = k1+@lam*k2 with |k1|<@l1 bits, |k1|<@l2
  bits.
  @P,@Q are points on @E s.t. @Q=k@P.
  Prints k1\nk2\n into filename.*/
void bsgs_simple(GEN E, GEN P, GEN Q, GEN lam, long l1, long l2,
                 char *filename);

/* Balanced Baby-step-giant-step on k = k1+@lam*k2 with |k1|<@l1 bits, |k1|<@l2
  bits. The expression k = k1+@lam*k2 is modified so that |k1|=|k2| for more
  efficiency.
  @P,@Q are points on @E s.t. @Q=k@P.
  Prints k1\nk2\n into filename.*/
void bsgs_balanced(GEN E, GEN P, GEN Q, GEN lam, long l1, long l2, long diff,
                   char *filename);

#endif