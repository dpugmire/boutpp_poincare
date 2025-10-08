#pragma once
// cli_args.h — header-only CLI parser for --xind, --maxpunc, --apar
// C++17

#include <cmath> // std::lround
#include <stdexcept>
#include <string>
#include <type_traits> // std::is_integral
#include <vector>

namespace cli
{

template <typename T>
struct Options
{
  std::vector<T> xind;   // values from --xind
  int maxpunc = 250;     // default
  bool haveXind = false; // whether --xind was provided
  std::string apar;      // value from --apar (empty if not provided)
};

// ---------- helpers ----------
inline bool isOption(const std::string& s)
{
  // Exactly two starting dashes means "option". A single dash (e.g., "-3") is a value.
  return s.size() >= 2 && s[0] == '-' && s[1] == '-';
}

inline double parseDouble(const std::string& s, const char* name)
{
  try
  {
    size_t pos = 0;
    double v = std::stod(s, &pos);
    if (pos != s.size())
      throw std::invalid_argument("trailing");
    return v;
  }
  catch (...)
  {
    throw std::runtime_error(std::string("Invalid number for ") + name + ": '" + s + "'");
  }
}

inline long parseLong(const std::string& s, const char* name)
{
  try
  {
    size_t pos = 0;
    long v = std::stol(s, &pos, 10);
    if (pos != s.size())
      throw std::invalid_argument("trailing");
    return v;
  }
  catch (...)
  {
    throw std::runtime_error(std::string("Invalid integer for ") + name + ": '" + s + "'");
  }
}

template <typename T>
inline std::vector<T> makeIntegerRange(long x0, long x1)
{
  std::vector<T> out;
  if (x0 == x1)
  {
    out.push_back(static_cast<T>(x0));
    return out;
  }
  const long step = (x1 > x0) ? 1 : -1;
  for (long x = x0;; x += step)
  {
    out.push_back(static_cast<T>(x));
    if (x == x1)
      break;
  }
  return out;
}

template <typename T>
inline T castVal(double v)
{
  if constexpr (std::is_integral<T>::value)
  {
    return static_cast<T>(std::lround(v)); // nearest int for integral T
  }
  else
  {
    return static_cast<T>(v); // pass-through for floating T
  }
}

// Floating linspace with EXACTLY `count` samples, including endpoints.
// For count <= 0: error; for count == 1: returns {a}.
template <typename T>
inline std::vector<T> linspaceCount(double a, double b, long count)
{
  if (count < 0)
    throw std::runtime_error("linspace count must be >= 0");
  std::vector<T> out;
  if (count == 0)
    return out;
  if (count == 1)
  {
    out.push_back(castVal<T>(a));
    return out;
  }

  out.reserve(static_cast<size_t>(count));
  const double n = static_cast<double>(count - 1);
  for (long i = 0; i < count; ++i)
  {
    const double t = static_cast<double>(i) / n;
    const double v = (1.0 - t) * a + t * b;
    out.push_back(castVal<T>(v));
  }
  return out;
}

// ---------- parser ----------
template <typename T>
inline void parseArgs(int argc, char** argv, Options<T>& opts)
{
  auto needValue = [&](int& i, const char* optName) -> std::string
  {
    if (i + 1 >= argc || isOption(argv[i + 1]))
    {
      throw std::runtime_error(std::string("Missing value for ") + optName);
    }
    return argv[++i];
  };

  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];

    // --maxpunc / --maxpunc=...
    if (arg == "--maxpunc")
    {
      std::string v = needValue(i, "--maxpunc");
      opts.maxpunc = static_cast<int>(parseLong(v, "--maxpunc"));
      continue;
    }
    if (arg.rfind("--maxpunc=", 0) == 0)
    {
      opts.maxpunc = static_cast<int>(parseLong(arg.substr(10), "--maxpunc"));
      continue;
    }

    // --apar / --apar=...
    if (arg == "--apar")
    {
      opts.apar = needValue(i, "--apar");
      continue;
    }
    if (arg.rfind("--apar=", 0) == 0)
    {
      opts.apar = arg.substr(7);
      continue;
    }

    // --xind with 1/2/3 following tokens
    if (arg == "--xind")
    {
      // Collect up to 3 following NON-option tokens (negatives like -3 are OK)
      std::vector<std::string> vals;
      int j = i;
      while (j + 1 < argc)
      {
        std::string next = argv[j + 1];
        if (isOption(next))
          break; // stop at next --option
        vals.emplace_back(std::move(next));
        ++j;
        if (static_cast<int>(vals.size()) == 3)
          break;
      }
      if (vals.empty())
        throw std::runtime_error("Missing value(s) for --xind");

      if (vals.size() == 1)
      {
        // --xind x -> [round(x)]
        const double x = parseDouble(vals[0], "--xind");
        opts.xind = { castVal<T>(x) };
      }
      else if (vals.size() == 2)
      {
        // --xind x0 x1 -> integer inclusive range (rounded endpoints)
        const double x0 = parseDouble(vals[0], "--xind");
        const double x1 = parseDouble(vals[1], "--xind");
        const long i0 = std::lround(x0);
        const long i1 = std::lround(x1);
        opts.xind = makeIntegerRange<T>(i0, i1);
      }
      else
      { // 3 values
        // --xind x0 x1 n -> EXACTLY n values between x0 and x1, endpoints included
        const double x0 = parseDouble(vals[0], "--xind");
        const double x1 = parseDouble(vals[1], "--xind");
        const long n = parseLong(vals[2], "--xind n");
        if (n < 0)
          throw std::runtime_error("--xind: n must be >= 0");
        opts.xind = linspaceCount<T>(x0, x1, n);
      }

      opts.haveXind = true;
      i = j; // advance past consumed tokens
      continue;
    }

    // --xind=VAL shorthand (single value)
    if (arg.rfind("--xind=", 0) == 0)
    {
      const double x = parseDouble(arg.substr(7), "--xind");
      opts.xind = { castVal<T>(x) };
      opts.haveXind = true;
      continue;
    }

    if (arg.rfind("--", 0) == 0)
    {
      throw std::runtime_error("Unknown option: " + arg);
    }
    else
    {
      throw std::runtime_error("Unexpected positional argument: " + arg);
    }
  }
}

} // namespace cli
