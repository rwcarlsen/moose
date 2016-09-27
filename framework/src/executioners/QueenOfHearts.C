
#include "QueenOfHearts.h"

QueenOfHearts::QueenOfHearts() { };

void
QueenOfHearts::addLoop(std::string name, ExecLoop* loop)
{
  _loop_names.push_back(name);
  _loop_counts.push_back(0);
  _loops.push_back(loop);
}

void
QueenOfHearts::run()
{
  runLoop(0);
}

void
QueenOfHearts::runLoop(int loop)
{
  if (loop >= _loop_names.size()) return;

  bool done = false;
  while (!done)
  {
    _loop_counts[loop]++;
    IterInfo info = {_loop_names[loop], _loop_counts[loop], _loop_counts};

    done = _loops[loop]->beginIter(info) || done;
    if (done) break;

    runLoop(loop + 1);

    done = _loops[loop]->endIter(info) || done;
  }

  // reset loop counts
  for (int i = loop; i < _loop_counts.size(); i++) _loop_counts[i] = 0;
}

