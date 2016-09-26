
#include <iostream>

#include "exeker.h"

class MyLoops {
 public:
  int foo_max;
  int bar_max;

  bool fooBegin(IterInfo info) {
    std::cout << "foo "<< info.curr_loop_iter << " begin\n";
    return info.curr_loop_iter >= foo_max;
  };

  bool fooEnd(IterInfo info) {
    std::cout << "foo end\n";
    return false;
  };

  bool barBegin(IterInfo info) {
    std::cout << "   bar " << info.curr_loop_iter << " begin\n";
    return info.curr_loop_iter >= bar_max;
  };

  bool barEnd(IterInfo info) {
    std::cout << "   bar end\n";
    return false;
  };
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
  MyLoops loops;
  loops.foo_max = 3;
  loops.bar_max = 2;

  Executioner ex;
  EXEC_LOOP_METHOD(ex, loops, foo);
  EXEC_LOOP_METHOD(ex, loops, bar);
  ex.addLoop("baz", new Cycle(2));
  ex.run();

  return 0;
}

