
#include <iostream>

#include "exeker.h"

class MyExeker: public Exeker {
 public:
  int foo_max;
  int bar_max;

  MyExeker() : _n_foo(0), _n_bar(0){
    EXEC_LEV_BEGIN(foo);
    EXEC_LEV_END(foo);
    EXEC_LEV_BEGIN(bar);
  }

 private:
  bool fooBegin() {
    _n_bar = 0;
    _n_foo++;
    std::cout << "foo begin " << _n_foo << "\n";
    return _n_foo >= foo_max;
  };

  bool fooEnd() {
    std::cout << "foo end\n";
    return false;
  };

  bool barBegin() {
    _n_bar++;
    std::cout << "   bar begin " << _n_bar << "\n";
    return _n_bar >= bar_max;
  };

  bool barEnd() {
    std::cout << "   bar end\n";
    return false;
  };

  int _n_foo;
  int _n_bar;
};

int main(int argc, char** argv) {
  MyExeker exec;
  exec.foo_max = 3;
  exec.bar_max = 2;
  exec.run();
  
  return 0;
}

