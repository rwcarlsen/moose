
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <memory>

class Formatter
{
public:
  virtual void format(std::ostream & stream, std::map<std::string, std::string> keyvals) = 0;
};

class SimpleFormatter : public Formatter
{
public:
  SimpleFormatter(std::set<std::string> pre = {}, std::set<std::string> post = {}) : _pre(pre), _post(post) { }

  virtual void format(std::ostream & stream, std::map<std::string, std::string> keyvals)
  {
    for (auto& key : _pre)
    {
      if (keyvals.count(key) > 0)
        stream << " " << key << "=" << keyvals[key];
    }

    for (auto it = keyvals.begin(); it != keyvals.end(); ++it)
    {
      if (_pre.count(it->first) > 0 || _post.count(it->first) > 0)
        continue;
      stream << " " << it->first << "=" << it->second;
    }

    for (auto& key : _post)
    {
      if (keyvals.count(key) > 0)
        stream << " " << key << "=" << keyvals[key];
    }

    stream << "\n";
  }

private:
  std::set<std::string> _pre;
  std::set<std::string> _post;
};

class AppLevelFormatter : public Formatter
{
public:
  AppLevelFormatter(Formatter * f) : _f(f) { }

  virtual void format(std::ostream & stream, std::map<std::string, std::string> keyvals)
  {
    if (keyvals.count("appnum") > 0)
      stream << "subapp " << keyvals["appnum"] << ": ";
    keyvals.erase("appnum");
    _f->format(stream, keyvals);
  }

private:
  std::unique_ptr<Formatter> _f;
};

class Logger
{
public:
  virtual void log(std::map<std::string, std::string> keyvals) = 0;
};

class ContextLogger : public Logger
{
public:
  ContextLogger(Logger * l, std::map<std::string, std::string> keyvals) : _logger(l), _ctx(keyvals) { }

  virtual void log(std::map<std::string, std::string> keyvals)
  {
    for (auto it = _ctx.begin(); it != _ctx.end(); ++it)
      keyvals.insert(*it);
    _logger->log(keyvals);
  }
private:
  std::unique_ptr<Logger> _logger;
  std::map<std::string, std::string> _ctx;
};

class FilterLogger : public Logger
{
public:
  FilterLogger(Logger * l, std::string key, std::set<std::string> show_vals) : _l(l), _key(key), _show(show_vals) { }

  virtual void log(std::map<std::string, std::string> keyvals)
  {
    if (keyvals.count(_key) == 0 || _show.count(keyvals[_key]) == 0)
      return;
    _l->log(keyvals);
  }

private:
  std::unique_ptr<Logger> _l;
  std::string _key;
  std::set<std::string> _show;
};

void logmsg(Logger * l, std::string level, std::string msg)
{
  l->log({{"lev", level}, {"msg", msg}});
}

class StreamLogger : public Logger
{
public:
  StreamLogger(std::ostream & stream, Formatter * f) : _stream(stream), _f(f) { }

  virtual void log(std::map<std::string, std::string> keyvals)
  {
    _f->format(_stream, keyvals);
  }

private:
  std::ostream & _stream;
  std::unique_ptr<Formatter> _f;
};

int main(int argc, char** argv)
{
  StreamLogger* l = new StreamLogger(std::cout, new AppLevelFormatter(new SimpleFormatter({"lev"}, {"msg"})));
  ContextLogger* l2 = new ContextLogger(l, {{"time", "never"},{"key1", "val1"}, {"appnum", "1.4.3"}});
  FilterLogger log(l2, "lev", {"4"});

  logmsg(&log, "4", "hello from logging");

  return 0;
}

