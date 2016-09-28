
#include "QueenOfHearts.h"

QueenOfHearts::QueenOfHearts() : _curr_loop(0) { };

QueenOfHearts::~QueenOfHearts() { };

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

int
QueenOfHearts::iter()
{
  return iter(_curr_loop);
}

int
QueenOfHearts::iter(int loop)
{
  if (loop >= _loop_counts.size())
    return 0;
  return _loop_counts[loop];
}

int
QueenOfHearts::iter(std::string loop)
{
  for (int i = 0; i < _loop_names.size(); i++)
  {
    if (_loop_names[i] == loop)
      return _loop_counts[i];
  }
  return 0; // TODO: throw here?
}


void
QueenOfHearts::runLoop(int loop)
{
  if (loop >= _loop_names.size()) return;
  _curr_loop = loop;

  bool done = false;
  while (!done)
  {
    _loop_counts[loop]++;
    IterInfo info = {_loop_names[loop], _loop_counts[loop], _loop_counts};

    done = _loops[loop]->beginIter(info) || done;
    if (done) break;

    runLoop(loop + 1);
    _curr_loop = loop;

    done = _loops[loop]->endIter(info) || done;
  }

  // reset loop counts
  for (int i = loop; i < _loop_counts.size(); i++) _loop_counts[i] = 0;
}

