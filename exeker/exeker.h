#pragma once

#include <map>
#include <vector>
#include <string>

#define EXEC_LEV(name) addLevelBegin(#name, [this]{return name##Begin();}); \
                       addLevelEnd(#name, [this]{return name##End();});

class Exeker {
 public: 
  Exeker() : _curr_level(0) {};

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

  int iter(int level) {
    return _level_counts[level];
  }

  int iter() {
    return _level_counts[_curr_level];
  }

 private:
  void createLevel(std::string name) {
    if (_level_index.count(name) > 0) {
      return;
    }
    int i = _level_names.size();
    _level_names.push_back(name);
    _level_counts.push_back(0);
    _level_index[name] = i;
    _begin_funcs.push_back([]{return false;});
    _end_funcs.push_back([]{return false;});
  }

  void execLevel(int level) {
    if (level >= _level_names.size()) {
      return;
    }
    _curr_level = level;

    bool done = false;
    while (!done) {
      _level_counts[level]++;
      //std::cout << std::string(level * 4, ' ') << _level_names[level] << " begin\n";
      done = _begin_funcs[level]() || done;
      execLevel(level + 1);
      _curr_level = level;
      //std::cout << std::string(level * 4, ' ') << _level_names[level] << " end\n";
      done = _end_funcs[level]() || done;
    }

    // reset level counts
    for (int i = level; i < _level_counts.size(); i++) {
      _level_counts[i] = 0;
    }
  }

  int _curr_level;
  std::vector<int> _level_counts;
  std::map<std::string, int> _level_index;
  std::vector<std::string> _level_names;
  std::vector<LevelFunc> _begin_funcs;
  std::vector<LevelFunc> _end_funcs;
};
