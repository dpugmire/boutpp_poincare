#include <adios2.h>

#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#if defined(CODEX_USE_MPI)
#include <mpi.h>
#endif

namespace
{

struct Options
{
  std::string filePath;
  std::string divertor;
  double iline = 0.0;
};

class MpiRuntime
{
public:
  MpiRuntime(int *argc, char ***argv)
  {
#if defined(CODEX_USE_MPI)
    MPI_Init(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif
  }

  ~MpiRuntime()
  {
#if defined(CODEX_USE_MPI)
    MPI_Finalize();
#endif
  }
};

bool parseOptions(int argc, char **argv, Options &options)
{
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--file" && i + 1 < argc)
    {
      options.filePath = argv[++i];
      continue;
    }
    if (arg == "--divertor" && i + 1 < argc)
    {
      options.divertor = argv[++i];
      continue;
    }
    if (arg == "--iline" && i + 1 < argc)
    {
      options.iline = std::stod(argv[++i]);
      continue;
    }

    std::cerr << "Unknown or incomplete argument: " << arg << "\n";
    return false;
  }

  if (options.filePath.empty())
  {
    std::cerr << "Missing --file\n";
    return false;
  }
  if (options.divertor != "single" && options.divertor != "circ")
  {
    std::cerr << "Missing or invalid --divertor single|circ\n";
    return false;
  }

  return true;
}

void require(bool condition, const std::string &message,
             std::vector<std::string> &failures)
{
  if (!condition)
    failures.push_back(message);
}

template <typename T>
T readScalarAttribute(adios2::IO &io, const std::string &name,
                      std::vector<std::string> &failures)
{
  auto attribute = io.InquireAttribute<T>(name);
  if (!attribute)
  {
    failures.push_back("missing attribute " + name);
    return T{};
  }

  const std::vector<T> data = attribute.Data();
  if (data.size() != 1)
  {
    failures.push_back("attribute " + name + " is not scalar");
    return T{};
  }
  return data[0];
}

template <typename T>
std::vector<T> readVector(adios2::IO &io, adios2::Engine &engine,
                          const std::string &name,
                          std::size_t expectedCount,
                          std::vector<std::string> &failures)
{
  auto variable = io.InquireVariable<T>(name);
  if (!variable)
  {
    failures.push_back("missing variable " + name);
    return {};
  }

  const adios2::Dims shape = variable.Shape();
  if (shape.size() != 1)
  {
    failures.push_back("variable " + name + " is not a 1D array");
    return {};
  }

  if (shape[0] != expectedCount)
  {
    failures.push_back("variable " + name + " has length " +
                       std::to_string(shape[0]) + ", expected " +
                       std::to_string(expectedCount));
  }

  std::vector<T> values(shape[0]);
  if (!values.empty())
    engine.Get(variable, values.data(), adios2::Mode::Sync);
  return values;
}

template <typename T>
void requireFiniteVector(const std::string &name, const std::vector<T> &values,
                         std::vector<std::string> &failures)
{
  for (std::size_t i = 0; i < values.size(); ++i)
  {
    if (!std::isfinite(values[i]))
    {
      failures.push_back("variable " + name + " contains a non-finite value at " +
                         std::to_string(i));
      return;
    }
  }
}

int validate(const Options &options)
{
#if defined(CODEX_USE_MPI)
  adios2::ADIOS adios(MPI_COMM_WORLD);
#else
  adios2::ADIOS adios;
#endif
  adios2::IO io = adios.DeclareIO("PoincareValidation");
  adios2::Engine engine =
      io.Open(options.filePath, adios2::Mode::ReadRandomAccess);

  std::vector<std::string> failures;

  const std::uint64_t schemaVersion =
      readScalarAttribute<std::uint64_t>(io, "schema_version", failures);
  const std::uint64_t numSeeds =
      readScalarAttribute<std::uint64_t>(io, "num_seeds", failures);
  const std::uint64_t numPunctures =
      readScalarAttribute<std::uint64_t>(io, "num_punctures", failures);
  const std::string divertor =
      readScalarAttribute<std::string>(io, "divertor", failures);
  const std::string traceEngine =
      readScalarAttribute<std::string>(io, "trace_engine", failures);
  const std::string viskoresDevice =
      readScalarAttribute<std::string>(io, "viskores_device", failures);
  const std::string viskoresOutputMode =
      readScalarAttribute<std::string>(io, "viskores_output_mode", failures);
  const std::string viskoresPrecision =
      readScalarAttribute<std::string>(io, "viskores_precision", failures);

  require(schemaVersion == 1, "schema_version is not 1", failures);
  require(numSeeds == 1, "num_seeds is not 1", failures);
  require(numPunctures > 0, "num_punctures is not positive", failures);
  require(divertor == options.divertor, "divertor attribute mismatch",
          failures);
  require(!traceEngine.empty(), "trace_engine attribute is empty", failures);
  require(!viskoresDevice.empty(), "viskores_device attribute is empty",
          failures);
  require(viskoresOutputMode == "punctures",
          "viskores_output_mode is not punctures", failures);
  require(!viskoresPrecision.empty(), "viskores_precision attribute is empty",
          failures);

  const std::size_t seedCount = static_cast<std::size_t>(numSeeds);
  const std::size_t punctureCount = static_cast<std::size_t>(numPunctures);

  const std::vector<std::uint64_t> seedId =
      readVector<std::uint64_t>(io, engine, "seed_id", seedCount, failures);
  const std::vector<double> iline =
      readVector<double>(io, engine, "iline", seedCount, failures);
  const std::vector<std::int32_t> endRegion =
      readVector<std::int32_t>(io, engine, "end_region", seedCount, failures);
  const std::vector<double> length =
      readVector<double>(io, engine, "length", seedCount, failures);
  const std::vector<std::uint64_t> offset =
      readVector<std::uint64_t>(io, engine, "offset", seedCount, failures);
  const std::vector<std::uint64_t> count =
      readVector<std::uint64_t>(io, engine, "count", seedCount, failures);

  const std::vector<std::int32_t> traceStep =
      readVector<std::int32_t>(io, engine, "trace_step", punctureCount,
                               failures);
  const std::vector<double> x =
      readVector<double>(io, engine, "x", punctureCount, failures);
  const std::vector<double> y =
      readVector<double>(io, engine, "y", punctureCount, failures);
  const std::vector<double> z =
      readVector<double>(io, engine, "z", punctureCount, failures);
  const std::vector<double> theta =
      readVector<double>(io, engine, "theta", punctureCount, failures);
  const std::vector<double> psi =
      readVector<double>(io, engine, "psi", punctureCount, failures);

  if (seedId.size() == 1)
    require(seedId[0] == 0, "seed_id[0] is not 0", failures);
  if (iline.size() == 1)
  {
    require(std::abs(iline[0] - options.iline) <= 1.0e-12,
            "iline[0] does not match expected line", failures);
  }
  if (offset.size() == 1)
    require(offset[0] == 0, "offset[0] is not 0", failures);
  if (count.size() == 1)
  {
    require(count[0] == numPunctures,
            "count[0] does not match num_punctures", failures);
  }

  requireFiniteVector("iline", iline, failures);
  requireFiniteVector("length", length, failures);
  requireFiniteVector("x", x, failures);
  requireFiniteVector("y", y, failures);
  requireFiniteVector("z", z, failures);
  requireFiniteVector("theta", theta, failures);
  requireFiniteVector("psi", psi, failures);

  for (std::size_t i = 0; i < traceStep.size(); ++i)
  {
    if (traceStep[i] < 0)
    {
      failures.push_back("trace_step contains a negative value at " +
                         std::to_string(i));
      break;
    }
  }

  (void)endRegion;
  engine.Close();

  std::cout << "Validated ADIOS Poincare output: " << options.filePath << "\n";
  std::cout << "divertor=" << divertor << " num_seeds=" << numSeeds
            << " num_punctures=" << numPunctures << "\n";

  if (!failures.empty())
  {
    std::cerr << "ADIOS validation failed:\n";
    for (const std::string &failure : failures)
      std::cerr << "  " << failure << "\n";
    return 1;
  }

  return 0;
}

} // namespace

int main(int argc, char **argv)
{
  MpiRuntime mpiRuntime(&argc, &argv);

  Options options;
  if (!parseOptions(argc, argv, options))
    return 2;

  try
  {
    return validate(options);
  }
  catch (const std::exception &error)
  {
    std::cerr << "ADIOS validation error: " << error.what() << "\n";
    return 1;
  }
}
