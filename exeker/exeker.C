
#include <iostream>

#include "exeker.h"

class MyLoops {
 public:
  int foo_max;
  int bar_max;

  MyLoops(Executioner& ex) {
    LOOP_FROM_METHODS(ex, foo);
    LOOP_FROM_METHODS(ex, bar);
  };

  bool fooBegin(IterInfo info) {
    std::cout << "foo "<< info.curr_loop_iter << " begin\n";
    a = info.curr_loop_iter;
    return info.curr_loop_iter >= foo_max;
  };

  bool fooEnd(IterInfo info) {
    std::cout << "foo end\n";
    return false;
  };

  bool barBegin(IterInfo info) {
    std::cout << "   bar " << info.curr_loop_iter << " begin\n";
    return info.curr_loop_iter >= a;
  };

  bool barEnd(IterInfo info) {
    std::cout << "   bar end\n";
    return false;
  };

 private:
  int a;
};

class Cycle : public ExecLoop {
 public:
  int cycle_max;

  Cycle(int max) : cycle_max(max) {};

  virtual bool beginIter(IterInfo info) {
    std::cout << "      baz " << info.curr_loop_iter << " begin\n";
    return info.curr_loop_iter >= cycle_max;
  }
  virtual bool endIter(IterInfo info) {
    std::cout << "      baz end\n";
    return false;
  }
};

int main(int argc, char** argv) {
  Executioner ex;
  MyLoops loops(ex);
  loops.foo_max = 3;
  loops.bar_max = 2;
  ex.addLoop("baz", new Cycle(2));
  ex.run();

  return 0;
}

