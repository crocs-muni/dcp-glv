#include "bsgs.h"
#include <array>
#include <assert.h>
#include <map>
#include <pari/pari.h>
#include <time.h>

void test_copy_INT() {

  GEN X = stoi(1), Y = stoi(2);
  pari_sp ltop = avma;

  X = copy_INT(X, Y);
  pari_sp lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  assert(ltop == lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  assert(equalii(X, Y));

  X = stoi(0), Y = stoi(2);
  ltop = avma;
  X = copy_INT(X, Y);
  lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  //  copy_INT does not work insitu for 0:
  assert(ltop > lbot);
  assert(equalii(X, Y));
}

void test_copy_INTMOD() {

  GEN e = ellinit(mkvec2(stoi(1), stoi(1)), stoi(101), 0);
  GEN P = ellmul(e, mkvec2(stoi(46), stoi(76)), stoi(1));
  GEN X = gel(P, 1), Y = gel(P, 2);
  pari_sp ltop = avma;

  copy_INTMOD(X, Y);
  pari_sp lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  assert(ltop == lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  assert(gequal(X, Y));

  P = ellmul(e, mkvec2(stoi(0), stoi(1)), stoi(1));
  X = gel(P, 1), Y = gel(P, 2);
  ltop = avma;
  copy_INTMOD(X, Y);
  lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  //  copy_INTMOD does not work insitu for 0:
  assert(ltop > lbot);
  assert(gequal(X, Y));
}

void test_copy_ellpoint() {

  GEN e = ellinit(mkvec2(stoi(1), stoi(1)), stoi(101), 0);
  GEN X = ellmul(e, mkvec2(stoi(46), stoi(76)), stoi(1));
  GEN Y = ellmul(e, mkvec2(stoi(11), stoi(63)), stoi(1));
  pari_sp ltop = avma;

  copy_ellpoint(X, Y);
  pari_sp lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  assert(ltop == lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  assert(gequal(X, Y));

  X = ellmul(e, mkvec2(stoi(0), stoi(1)), stoi(1));
  ltop = avma;
  copy_ellpoint(X, Y);
  lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  //  copy_ellpoint does not work insitu for points with 0:
  assert(ltop > lbot);
  assert(gequal(X, Y));
}

void test_add_inplace() {
  GEN X = stoi(1), Y = stoi(2), Z = stoi(3);
  pari_sp ltop = avma;

  X = add_inplace(X, Y);
  pari_sp lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  assert(ltop == lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  assert(gequal(X, Z));

  X = stoi(0), Y = stoi(2);
  ltop = avma;
  X = add_inplace(X, Y);
  lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  // pari_printf("%Ps %Ps\n",X,Y);
  //  add_inplace does not work insitu with 0:
  assert(ltop > lbot);
  assert(gequal(X, Y));
}

void test_elladd_inplace() {

  GEN e = ellinit(mkvec2(stoi(1), stoi(1)), stoi(101), 0);
  GEN X = ellmul(e, mkvec2(stoi(46), stoi(76)), stoi(1));
  GEN Y = ellmul(e, mkvec2(stoi(11), stoi(63)), stoi(1));
  GEN Z = ellmul(e, mkvec2(stoi(57), stoi(44)), stoi(1));
  pari_sp ltop = avma;

  X = elladd_inplace(e, X, Y);
  pari_sp lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  assert(ltop == lbot);
  // pari_printf("%Ps %Ps\n",X,Z);
  assert(gequal(X, Z));

  X = ellinf();
  Z = ellmul(e, mkvec2(stoi(60), stoi(74)), stoi(1));
  ltop = avma;
  X = elladd_inplace(e, X, Z);
  lbot = avma;
  // pari_printf("%p %p\n",ltop,lbot);
  // pari_printf("%Ps %Ps\n",X,Z);
  // elladd_inplace is not insitu for points with 0
  assert(ltop > lbot);
  assert(gequal(X, Z));
}

int main() {
  long stack = 100000;
  pari_init(stack, 65536);

  test_copy_INT();
  test_copy_INTMOD();
  test_copy_ellpoint();
  test_add_inplace();
  test_elladd_inplace();

  pari_close();

  printf("All good.\n");
  return 0;
}