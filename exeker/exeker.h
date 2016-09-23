#pragma once

#include <map>
#include <vector>
#include <string>

#define EXEC_LEV_BEGIN(name) addLevelBegin(#name, [this]{return name##Begin();});
#define EXEC_LEV_END(name) addLevelEnd(#name, [this]{return name##End();});

class Exeker {
 public: 
  typedef std::function< bool() > LevelFunc;

  void addLevelBegin(std::string name, LevelFunc begin) {
    createLevel(name);
    _begin_funcs[_level_index[name]] = begin;
  }

  void addLevelEnd(std::string name, LevelFunc end) {
    createLevel(name);
    _end_funcs[_level_index[name]] = end;
  }

  void run() {
    execLevel(0);
  }

 private:
  void createLevel(std::string name) {
    if (_level_index.count(name) > 0) {
      return;
    }
    int i = _level_names.size();
    _level_names.push_back(name);
    _level_index[name] = i;
    _begin_funcs.push_back([]{return false;});
    _end_funcs.push_back([]{return false;});
  }

  void execLevel(int level) {
    if (level >= _level_names.size()) {
      return;
    }

    bool done = false;
    while (!done) {
      done = _begin_funcs[level]() || done;
      execLevel(level + 1);
      done = _end_funcs[level]() || done;
    }
  }

  std::map<std::string, int> _level_index;
  std::vector<std::string> _level_names;
  std::vector<LevelFunc> _begin_funcs;
  std::vector<LevelFunc> _end_funcs;
};
