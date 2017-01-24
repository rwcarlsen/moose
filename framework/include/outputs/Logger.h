
#include <map>
#include <set>
#include <vector>
#include <ostream>
#include <memory>

#include "Moose.h"

class Formatter
{
public:
  virtual void format(std::ostream & stream, std::map<std::string, std::string> keyvals) = 0;
};

class Logger
{
public:
  virtual void log(std::map<std::string, std::string> keyvals) = 0;
};

class SimpleFormatter : public Formatter
{
public:
  SimpleFormatter(std::vector<std::string> pre = {}, std::vector<std::string> post = {}, bool allow_others = true);
  virtual void format(std::ostream & stream, std::map<std::string, std::string> keyvals);
private:
  std::vector<std::string> _pre;
  std::vector<std::string> _post;
  std::set<std::string> _preset;
  std::set<std::string> _postset;
  bool _allow_others;
};

class PrefixFormatter : public Formatter
{
public:
  PrefixFormatter(Formatter * f, std::string key, std::string color = COLOR_DEFAULT);
  virtual void format(std::ostream & stream, std::map<std::string, std::string> keyvals);
private:
  std::unique_ptr<Formatter> _f;
  std::string _key;
  std::string _color;
};

class ContextLogger : public Logger
{
public:
  ContextLogger(Logger * l, std::map<std::string, std::string> keyvals);
  virtual void log(std::map<std::string, std::string> keyvals);
private:
  std::unique_ptr<Logger> _logger;
  std::map<std::string, std::string> _ctx;
};

class TimeLogger : public Logger
{
public:
  TimeLogger(Logger * l, const std::string & key);
  virtual void log(std::map<std::string, std::string> keyvals);
private:
  std::unique_ptr<Logger> _logger;
  std::string _key;
};

class FilterLogger : public Logger
{
public:
  FilterLogger(Logger * l, std::string key, std::set<std::string> show_vals);
  virtual void log(std::map<std::string, std::string> keyvals);
private:
  std::unique_ptr<Logger> _l;
  std::string _key;
  std::set<std::string> _show;
};

class StreamLogger : public Logger
{
public:
  StreamLogger(std::ostream & stream, Formatter * f);
  virtual void log(std::map<std::string, std::string> keyvals);
private:
  std::ostream & _stream;
  std::unique_ptr<Formatter> _f;
};

