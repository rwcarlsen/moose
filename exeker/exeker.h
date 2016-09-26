#pragma once

#include <map>
#include <vector>
#include <string>
#include <functional>

struct IterInfo {
  std::string& name;
  const int& curr_loop_iter;
  const std::vector<int>& all_iters;
};

class ExecLoop {
 public:
  virtual bool beginIter(IterInfo info) = 0;
  virtual bool endIter(IterInfo info) = 0;
};

class Executioner final {
 public: 
  Executioner() : _curr_loop(0) {};

  void addLoop(std::string name, ExecLoop* loop) {
    _loop_names.push_back(name);
    _loop_counts.push_back(0);
    _loops.push_back(loop);
  }

  void run() { runLoop(0); }

 private:
  void runLoop(int loop) {
    if (loop >= _loop_names.size()) {
      return;
    }
    _curr_loop = loop;

    bool done = false;
    while (!done) {
      _loop_counts[loop]++;
      IterInfo info = {_loop_names[loop], _loop_counts[loop], _loop_counts};

      //std::cout << std::string(loop * 4, ' ') << _loop_names[loop] << " " << info.curr_loop_iter << " begin\n";
      done = _loops[loop]->beginIter(info) || done;
      runLoop(loop + 1);
      _curr_loop = loop;
      //std::cout << std::string(loop * 4, ' ') << _loop_names[loop] << " end\n";
      done = _loops[loop]->endIter(info) || done;
    }

    // reset loop counts
    for (int i = loop; i < _loop_counts.size(); i++) {
      _loop_counts[i] = 0;
    }
  }

  int _curr_loop;
  std::vector<int> _loop_counts;
  std::vector<std::string> _loop_names;
  std::vector<ExecLoop*> _loops;
};

///////////////////////////////////////////////////////////////////////////////
// all code below here is bonus for using labmdas and custom-named methods for
// loop funcs all together in a single object/class
///////////////////////////////////////////////////////////////////////////////

#define LOOP_FROM_METHODS(ex,name) ex.addLoop( #name, \
        new ExecLoopFunc([this](IterInfo info){return name##Begin(info);}, \
                         [this](IterInfo info){return name##End(info);}) \
    ); \

typedef std::function< bool(IterInfo info) > LoopFunc;

class ExecLoopFunc : public ExecLoop {
 public:
  ExecLoopFunc(LoopFunc begin, LoopFunc end) : _begin(begin), _end(end) {};
  virtual bool beginIter(IterInfo info) override { return _begin(info);}
  virtual bool endIter(IterInfo info) override { return _end(info);}

 private:
  LoopFunc _begin;
  LoopFunc _end;
};

