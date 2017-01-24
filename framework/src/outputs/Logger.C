
#include "Logger.h"

// Replaces text in msg of the form "@key@" with the value stored in the
// keyvals[key] if the key exists in keyvals.
void
formatMessage(std::string & msg, std::map<std::string, std::string> & keyvals)
{
  size_t start_pos = 0;
  while((start_pos = msg.find("@", start_pos)) != std::string::npos)
  {
    size_t end_pos = msg.find("@", start_pos+1);
    auto len = end_pos - start_pos;
    if (len == 1 || end_pos == std::string::npos)
    {
      start_pos = end_pos + 1;
      continue;
    }

    std::string key = msg.substr(start_pos+1, len - 1);
    if (keyvals.count(key) == 0)
      start_pos = end_pos;
    else
    {
      msg.replace(start_pos, len+1, keyvals[key]);
      start_pos += len + 1;
    }
  }
}

void
formatMultiline(std::ostream & os, const std::string & prefix, std::string msg)
{
  msg = "\n" + msg;
  std::string from = "\n";
  std::string to = "\n" + prefix + "    ";
  size_t start_pos = 0;
  while((start_pos = msg.find(from, start_pos)) != std::string::npos)
  {
      msg.replace(start_pos, from.length(), to);
      start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
  os << msg;
}

bool
hasPrefix(std::string s, std::string prefix)
{
  auto res = std::mismatch(prefix.begin(), prefix.end(), s.begin());
  return res.first == prefix.end();
}

void
formatKey(std::ostream & os, const std::string & key, std::map<std::string, std::string> & keyvals)
{
  if (keyvals.count(key) == 0 || hasPrefix(key, "__"))
    return;

  std::string val = keyvals[key];
  if (key == "msg")
  {
    formatMessage(val, keyvals);
    os << " ";
    if (keyvals.count("__prefix") > 0 && val.find('\n') != std::string::npos)
      formatMultiline(os, keyvals["__prefix"], val);
    else
      os << val;
  }
  else if (val.find(' ') != std::string::npos)
    os << " " << key << "='" << val << "'";
  else
    os << " " << key << "=" << val;
}

SimpleFormatter::SimpleFormatter(std::vector<std::string> pre, std::vector<std::string> post, bool allow_others)
  : _pre(pre), _post(post), _allow_others(allow_others)
{
  for (unsigned int i = 0; i < pre.size(); i++)
    _preset.insert(pre[i]);
  for (unsigned int i = 0; i < post.size(); i++)
    _postset.insert(post[i]);
}

void
SimpleFormatter::format(std::ostream & stream, std::map<std::string, std::string> keyvals)
{
  stream << COLOR_DEFAULT;
  for (auto& key : _pre)
    formatKey(stream, key, keyvals);

  if (_allow_others)
  {
    for (auto & it : keyvals)
    {
      if ( _preset.count(it.first) == 0 && _postset.count(it.first) == 0)
        formatKey(stream, it.first, keyvals);
    }
  }

  for (auto& key : _post)
    formatKey(stream, key, keyvals);

  stream << "\n";
}

PrefixFormatter::PrefixFormatter(Formatter * f, std::string key, std::string color)
  : _f(f), _key(key), _color(color) { }

void
PrefixFormatter::format(std::ostream & stream, std::map<std::string, std::string> keyvals)
{
  if (keyvals.count(_key) > 0)
  {
    auto & val = keyvals[_key];
    if (_key != "appname" || val != "main")
    {
      std::string prefix = _color + val + ":" + COLOR_DEFAULT;
      stream << prefix;
      keyvals.insert({"__prefix", prefix});
    }
  }
  _f->format(stream, keyvals);
}

ContextLogger::ContextLogger(Logger * l, std::map<std::string, std::string> keyvals)
  : _logger(l), _ctx(keyvals) { }

void
ContextLogger::log(std::map<std::string, std::string> keyvals)
{
  for (auto it = _ctx.begin(); it != _ctx.end(); ++it)
    keyvals.insert(*it);
  _logger->log(keyvals);
}

TimeLogger::TimeLogger(Logger * l, const std::string & key) : _logger(l), _key(key) { }

void
TimeLogger::log(std::map<std::string, std::string> keyvals)
{
  time_t now = time(NULL);
  struct tm tstruct;
  char buf[40];
  tstruct = *localtime(&now);
  //format: HH:mm:ss
  strftime(buf, sizeof(buf), "%X", &tstruct);
  keyvals.insert({_key, std::string(buf)});
  _logger->log(keyvals);
}

FilterLogger::FilterLogger(Logger * l, std::string key, std::set<std::string> show_vals)
  : _l(l), _key(key), _show(show_vals) { }

void
FilterLogger::log(std::map<std::string, std::string> keyvals)
{
  if (keyvals.count(_key) == 0 || _show.count(keyvals[_key]) == 0)
    return;
  _l->log(keyvals);
}

StreamLogger::StreamLogger(std::ostream & stream, Formatter * f) : _stream(stream), _f(f) { }

void
StreamLogger::log(std::map<std::string, std::string> keyvals)
{
  _f->format(_stream, keyvals);
}

