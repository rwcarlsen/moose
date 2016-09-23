
#include <iostream>

#include "exeker.h"

class MyExeker: public Exeker {
 public:
  int foo_max;
  int bar_max;

  MyExeker() {
    EXEC_LEV(foo);
    EXEC_LEV(bar);
  }

 private:
  bool fooBegin() {
    std::cout << "foo "<< iter() << " begin\n";
    return iter() >= foo_max;
  };

  bool fooEnd() {
    std::cout << "foo end\n";
    return false;
  };

  bool barBegin() {
    std::cout << "   bar " << iter() << " begin\n";
    return iter() >= bar_max;
  };

  bool barEnd() {
    std::cout << "   bar end\n";
    return false;
  };
};

int main(int argc, char** argv) {
  MyExeker exec;
  exec.foo_max = 3;
  exec.bar_max = 2;
  exec.run();
  
  return 0;
}

