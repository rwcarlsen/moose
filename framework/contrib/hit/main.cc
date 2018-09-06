
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <set>

#include "parse.h"

class Flags;
std::vector<std::string> parseOpts(int argc, char ** argv, Flags & flags);

int findParam(int argc, char ** argv);
int validate(int argc, char ** argv);
int format(int argc, char ** argv);

std::string
absPath(const std::string & path)
{
  char abspath[PATH_MAX + 1];
  realpath(path.c_str(), abspath);
  return std::string(abspath);
}


int
main(int argc, char ** argv)
{
  if (argc < 2)
  {
    std::cerr << "must specify a subcommand\n";
    return 1;
  }

  std::string subcmd(argv[1]);

  if (subcmd == "find")
    return findParam(argc - 2, argv + 2);
  else if (subcmd == "validate")
    return validate(argc - 2, argv + 2);
  else if (subcmd == "format")
    return format(argc - 2, argv + 2);

  std::cerr << "unrecognized subcommand '" << subcmd << "'\n";
  return 1;
}

struct Flag
{
  bool arg;
  bool have;
  std::string val;
  std::string help;
};

class Flags
{
public:
  Flags(const std::string & usage) : usage_msg(usage) {}
  void add(std::string name, std::string help, std::string def = "__NONE__")
  {
    if (def == "__NONE__")
      flags[name] = {false, false, "", help};
    else
      flags[name] = {true, false, def, help};
  }

  bool have(std::string flag) { return flags.count(flag) > 0 && flags[flag].have; }
  std::string val(std::string flag) { return flags[flag].val; }
  std::string usage()
  {
    std::stringstream ss;
    ss << usage_msg << "\n";
    for (auto & pair : flags)
    {
      auto flag = pair.second;
      if (flag.arg)
        ss << "-" << pair.first << " <arg>    " << flag.help << " (default='" << flag.val << "')\n";
      else
        ss << "-" << pair.first << "    " << flag.help << " (default=false)\n";
    }
    return ss.str();
  }

  std::map<std::string, Flag> flags;
  std::string usage_msg;
};

std::vector<std::string>
parseOpts(int argc, char ** argv, Flags & flags)
{
  int i = 0;
  for (; i < argc; i++)
  {
    std::string arg = argv[i];
    if (arg[0] != '-')
      break;

    std::string flagname = arg.substr(1);
    if (flagname[0] == '-')
      flagname = flagname.substr(1);

    if (flags.flags.count(flagname) == 0)
      throw std::runtime_error("unknown flag '" + arg);

    auto & flag = flags.flags[flagname];
    flag.have = true;
    if (flag.arg)
    {
      i++;
      flag.val = argv[i];
    }
  }

  std::vector<std::string> positional;
  for (; i < argc; i++)
    positional.push_back(argv[i]);
  return positional;
}

inline std::string
errormsg(std::string /*fname*/, hit::Node * /*n*/)
{
  return "";
}

template <typename T, typename... Args>
std::string
errormsg(std::string fname, hit::Node * n, T arg, Args... args)
{
  std::stringstream ss;
  if (n && fname.size() > 0)
    ss << fname << ":" << n->line() << ": ";
  else if (fname.size() > 0)
    ss << fname << ":0: ";
  ss << arg;
  ss << errormsg("", nullptr, args...);
  return ss.str();
}

class DupParamWalker : public hit::Walker
{
public:
  DupParamWalker(std::string fname) : _fname(fname) {}
  void walk(const std::string & fullpath, const std::string & /*nodepath*/, hit::Node * n) override
  {
    std::string prefix = n->type() == hit::NodeType::Field ? "parameter" : "section";

    if (_have.count(fullpath) > 0)
    {
      auto existing = _have[fullpath];
      if (_duplicates.count(fullpath) == 0)
      {
        errors.push_back(
            errormsg(_fname, existing, prefix, " '", fullpath, "' supplied multiple times"));
        _duplicates.insert(fullpath);
      }
      errors.push_back(errormsg(_fname, n, prefix, " '", fullpath, "' supplied multiple times"));
    }
    _have[n->fullpath()] = n;
  }

  std::vector<std::string> errors;

private:
  std::string _fname;
  std::set<std::string> _duplicates;
  std::map<std::string, hit::Node *> _have;
};

int
findParam(int argc, char ** argv)
{
  Flags flags("hit find [flags] <parameter-path> <file>...");
  flags.add("f", "only show file name");
  auto positional = parseOpts(argc, argv, flags);

  if (positional.size() < 2)
  {
    std::cerr << flags.usage();
    return 1;
  }

  std::string srcpath(positional[0]);

  int ret = 0;
  for (int i = 1; i < positional.size(); i++)
  {
    std::string fname(positional[i]);
    std::ifstream f(fname);
    std::string input((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    hit::Node * root = nullptr;
    try
    {
      root = hit::parse(fname, input);
    }
    catch (std::exception & err)
    {
      std::cerr << err.what() << "\n";
      ret = 1;
      continue;
    }

    auto n = root->find(srcpath);
    if (n)
    {
      if (flags.have("f"))
        std::cout << fname << "\n";
      else
        std::cout << fname << ":" << n->line() << "\n";
    }
  }

  return ret;
}

// the style file is of the format:
//
//     [format]
//       # these parameters set the correspondingly named Formatter member variables.  See them
//       # for detailed descriptions of what they do.
//       indent_string = "  "
//       line_length = 100
//       canonical_section_markers = true
//
//       # This section specifies parameter/section sorgin order as provided by the formatter's
//       # addPattern member function. The content for this section mirrors the structure of
//       # the hit file being formatted.  Each section name is a regex (limited to valid hit
//       # identifier characters). The fields and subsections within each section specify an
//       # order; any field values are ignored. See the docs for that function for more details.
//       # The sorting rules are example-based. Sorting is performed depth-first by rules in
//       # lexical order
//       [sorting]
//         [.*]
//           [.*]
//             # 'first' fields in doubly-nested section go first. This rule must go first in
//             # order to not override any later-specified ordering of higher-level sections
//             # that would occur because the higher-level section matchers for this rule
//             # would match every section.
//             first = FooBar
//           []
//         []
//         [foo]        # section 'foo' goes first (before other sections)
//           bar = bla  # field 'bar' (in the foo section) goes first
//           ** = bla   # double glob is placeholder for unordered fields/sections
//           baz = bla  # field 'baz' goes last
//         []
//       []
//     []
//
// where all fields are optional and the sorting section is also optional.  If the sorting section
// is present, you can have as many patterns as you want.
int
format(int argc, char ** argv)
{
  Flags flags("hit format [flags] <file>...");
  flags.add("h", "print help");
  flags.add("help", "print help");
  flags.add("i", "modify file(s) inplace");
  flags.add("style", "hit style file detailing format to use", ".hit-format");

  auto positional = parseOpts(argc, argv, flags);

  if (flags.have("h") || flags.have("help"))
  {
    std::cout << flags.usage();
    return 0;
  }

  if (positional.size() < 1)
  {
    std::cout << flags.usage();
    return 1;
  }

  hit::Formatter fmt;

  int ret = 0;
  for (int i = 0; i < positional.size(); i++)
  {
    std::string fname(positional[i]);
    std::ifstream f(fname);

    try
    {
      if (flags.have("i"))
      {
        std::stringstream dst;
        hit::format(fname, f, dst, flags.val("style"));
        std::ofstream output(fname);
        output << dst.str();
      }
      else
        hit::format(fname, f, std::cout, flags.val("style"));
    }
    catch (std::exception & err)
    {
      std::cerr << err.what() << "\n";
      ret = 1;
      continue;
    }
  }

  return ret;
}

int
validate(int argc, char ** argv)
{
  if (argc < 1)
  {
    std::cerr << "please pass in an input file argument\n";
    return 1;
  }

  int ret = 0;
  for (int i = 0; i < argc; i++)
  {
    std::string fname(argv[i]);
    std::ifstream f(fname);
    std::string input((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    std::unique_ptr<hit::Node> root;
    try
    {
      root.reset(hit::parse(fname, input));
    }
    catch (std::exception & err)
    {
      std::cout << err.what() << "\n";
      ret = 1;
      continue;
    }

    DupParamWalker w(fname);
    root->walk(&w, hit::NodeType::Field);
    for (auto & msg : w.errors)
      std::cout << msg << "\n";
  }
  return ret;
}
