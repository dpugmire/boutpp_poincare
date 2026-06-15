#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <netcdf.h>

#if defined(CODEX_USE_MPI)
#include <mpi.h>
#endif

#if defined(CODEX_USE_VISKORES)
#include <viskores/cont/Initialize.h>
#include <viskores/cont/RuntimeDeviceInformation.h>
#include <viskores/cont/RuntimeDeviceTracker.h>
#endif

#include "AparData.h"
#include "AparFieldModel.h"
#include "FieldLineIntegrator.h"
#include "PoincareOutput.h"
#include "ValidationSuite.h"
#if defined(CODEX_USE_ADIOS)
#include "AdiosPoincareOutput.h"
#endif
#if defined(CODEX_USE_VISKORES)
#include "ViskoresFieldLineTracer.h"
#endif

namespace
{

#if defined(CODEX_USE_MPI)
class MpiRuntime
{
public:
  MpiRuntime(int *argc, char ***argv)
  {
    MPI_Init(argc, argv);
    initialized_ = true;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  ~MpiRuntime()
  {
    if (initialized_)
      MPI_Finalize();
  }

  void barrier() const
  {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  int allreduceMaxInt(int localValue) const
  {
    int globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    return globalValue;
  }

  int allreduceMinInt(int localValue) const
  {
    int globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    return globalValue;
  }

  long long allreduceSumLongLong(long long localValue) const
  {
    long long globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    return globalValue;
  }

  double allreduceMaxDouble(double localValue) const
  {
    double globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    return globalValue;
  }

  double allreduceMinDouble(double localValue) const
  {
    double globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    return globalValue;
  }

  double allreduceSumDouble(double localValue) const
  {
    double globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    return globalValue;
  }

  void abortAll(int code) const
  {
    MPI_Abort(MPI_COMM_WORLD, code);
  }

  int rank = 0;
  int size = 1;

private:
  bool initialized_ = false;
};
#else
class MpiRuntime
{
public:
  explicit MpiRuntime(int *, char ***)
  {
  }

  void barrier() const
  {
  }

  int allreduceMaxInt(int localValue) const
  {
    return localValue;
  }

  int allreduceMinInt(int localValue) const
  {
    return localValue;
  }

  long long allreduceSumLongLong(long long localValue) const
  {
    return localValue;
  }

  double allreduceMaxDouble(double localValue) const
  {
    return localValue;
  }

  double allreduceMinDouble(double localValue) const
  {
    return localValue;
  }

  double allreduceSumDouble(double localValue) const
  {
    return localValue;
  }

  void abortAll(int code) const
  {
    std::exit(code);
  }

  int rank = 0;
  int size = 1;
};
#endif

struct TraceTask
{
  std::string divertorTag;
  Point3D seedInd;
};

struct PsiNormalization
{
  bool enabled = false;
  double psiAxis = 0.0;
  double psiBndry = 0.0;
};

void ncCheckMain(int status, const std::string &where)
{
  if (status != NC_NOERR)
    throw std::runtime_error(where + ": " + nc_strerror(status));
}

double readScalarDoubleFromNetcdf(int ncid, const std::string &name)
{
  int varid = -1;
  ncCheckMain(nc_inq_varid(ncid, name.c_str(), &varid),
              "nc_inq_varid(" + name + ")");

  double value = 0.0;
  ncCheckMain(nc_get_var_double(ncid, varid, &value),
              "nc_get_var_double(" + name + ")");
  return value;
}

PsiNormalization readPsiNormalizationFromGridFile(
    const std::string &gridFile)
{
  int ncid = -1;
  ncCheckMain(nc_open(gridFile.c_str(), NC_NOWRITE, &ncid),
              "nc_open(" + gridFile + ")");

  PsiNormalization normalization;
  try
  {
    normalization.enabled = true;
    normalization.psiAxis = readScalarDoubleFromNetcdf(ncid, "psi_axis");
    normalization.psiBndry = readScalarDoubleFromNetcdf(ncid, "psi_bndry");
  }
  catch (...)
  {
    nc_close(ncid);
    throw;
  }

  ncCheckMain(nc_close(ncid), "nc_close(" + gridFile + ")");
  return normalization;
}

const char *defaultOutputFormatChoice()
{
#if defined(CODEX_USE_ADIOS)
  return "adios";
#else
  return "text";
#endif
}

void printUsage(const char *prog)
{
  std::cout << "Usage: " << prog << " [options]\n"
            << "Options:\n"
            << "  --reference-dir DIR   MATLAB reference directory (default: "
               "/Users/dpn/proj/bout++/ben_zhu_poincare/zperiod_5)\n"
            << "  --output FILE         ADIOS BP output path (default: "
               "poincare.bp)\n"
            << "  --text-output-dir DIR Text/debug output directory "
               "(default: ./outputs)\n"
            << "  --output-dir DIR      Alias for --text-output-dir\n"
            << "  --output-format FMT   text|adios|both (default: "
            << defaultOutputFormatChoice() << ")\n"
            << "  --adios-file FILE     Alias for --output\n"
            << "  --psi-axis VALUE      Grid NetCDF psi_axis value for "
               "normalized psi output\n"
            << "  --psi-bndry VALUE     Grid NetCDF psi_bndry value for "
               "normalized psi output\n"
            << "  --grid-file FILE      Read psi_axis and psi_bndry from a "
               "BOUT++ grid NetCDF file\n"
            << "  --single-apar FILE    apar.single.nc path\n"
            << "  --circ-apar FILE      apar.circ.nc path\n"
            << "  --divertor TAG        all|single|circ (default: all)\n"
            << "  --lines LIST          comma-separated x-index values (double "
               "supported)\n"
            << "  --nlines X0 X1 N      generate N evenly-spaced lines from X0 "
               "to X1 (double supported)\n"
            << "  --direction DIR       +1 or -1 (default: 1)\n"
            << "  --np-max N            max punctures per line (default: 100)\n"
            << "  --max-steps N         max integration steps per line "
               "(default: 200 * np-max)\n"
            << "  --tol VALUE           max-abs tolerance (default: 1e-8)\n"
            << "  --trace-engine ENG    tracing backend: cpu|viskores "
               "(default: viskores)\n"
            << "  --viskores-device DEV viskores device: serial|openmp|kokkos "
               "(default: serial)\n"
            << "  --viskores-output MODE viskores output mode: "
               "punctures|rk4|states (default: punctures)\n"
            << "  --puncture-detection MODE compact Viskores puncture "
               "detection: on|off (default: on)\n"
            << "  --puncture-refine MODE compact Viskores puncture "
               "refinement: on|off (default: on)\n"
            << "  --trace-diagnostics MODE compact Viskores diagnostic "
               "counters: on|off (default: off)\n"
            << "  --per-seed-log       print one trace-result line per seed "
               "(default: off)\n"
            << "  --compare             enable MATLAB comparison mode (default "
               "is trace-only)\n"
            << "  --help                show this help\n";
}

bool hasHelpArgument(int argc, char **argv)
{
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h")
      return true;
  }
  return false;
}

std::vector<double> parseLinesCsvDoubles(const std::string &text)
{
  std::vector<double> out;
  std::stringstream ss(text);
  std::string token;
  while (std::getline(ss, token, ','))
  {
    if (token.empty())
      continue;
    out.push_back(std::stod(token));
  }
  return out;
}

std::vector<double> buildLineRange(double x0, double x1, int n)
{
  std::vector<double> out;
  if (n <= 0)
    return out;
  out.reserve(static_cast<size_t>(n));

  if (n == 1)
  {
    out.push_back(x0);
    return out;
  }

  for (int i = 0; i < n; ++i)
  {
    const double t = static_cast<double>(i) / static_cast<double>(n - 1);
    out.push_back(x0 + t * (x1 - x0));
  }
  return out;
}

int computeMaxStateCount(const TraceOptions &options)
{
  constexpr int kDefaultMaxStepsPerPuncture = 200;
  const long long derivedMaxStepsLL =
      static_cast<long long>(kDefaultMaxStepsPerPuncture) *
      static_cast<long long>(std::max(1, options.npMax));
  const int derivedMaxSteps =
      (derivedMaxStepsLL >
       static_cast<long long>(std::numeric_limits<int>::max()))
          ? std::numeric_limits<int>::max()
          : static_cast<int>(std::max(1LL, derivedMaxStepsLL));

  const int maxStepsCap =
      (options.maxSteps > 0) ? options.maxSteps : derivedMaxSteps;
  const int clampedSteps = std::max(1, maxStepsCap);
  return (clampedSteps >= std::numeric_limits<int>::max())
             ? std::numeric_limits<int>::max()
             : (clampedSteps + 1);
}

bool isIntegerVal(double v)
{
  return std::isfinite(v) && std::fabs(v - std::round(v)) < 1.0e-12;
}

enum class LineSpecMode
{
  none,
  csv,
  range
};

enum class TraceEngine
{
  Cpu,
  Viskores
};

enum class OutputFormat
{
  Text,
  Adios,
  Both
};

bool isFatalTraceStatus(TraceStatus status)
{
  return status == TraceStatus::InvalidConfiguration ||
         status == TraceStatus::OutputTooSmall;
}

int traceStatusSummaryIndex(TraceStatus status)
{
  switch (status)
  {
  case TraceStatus::Ok:
    return 0;
  case TraceStatus::InvalidSeed:
    return 1;
  case TraceStatus::InvalidConfiguration:
    return 2;
  case TraceStatus::OutputTooSmall:
    return 3;
  case TraceStatus::MaxStepLimitReached:
    return 4;
  }
  return 4;
}

TraceStatus traceStatusFromSummaryIndex(int index)
{
  switch (index)
  {
  case 0:
    return TraceStatus::Ok;
  case 1:
    return TraceStatus::InvalidSeed;
  case 2:
    return TraceStatus::InvalidConfiguration;
  case 3:
    return TraceStatus::OutputTooSmall;
  case 4:
    return TraceStatus::MaxStepLimitReached;
  }
  return TraceStatus::MaxStepLimitReached;
}

void printTraceResultSummary(
    const MpiRuntime &mpi, const std::vector<TraceStatus> &traceStatusPerSeed,
    const std::vector<int> &stateCountPerSeed,
    const std::vector<int> &endRegionPerSeed,
    const std::vector<int> &punctureCountPerSeed,
    const std::vector<double> &connectionLengthPerSeed, int maxPuncPerSeed,
    const TraceDiagnostics *diagnostics)
{
  constexpr int kStatusCount = 5;
  constexpr int kEndRegionOtherSlot = 100;
  constexpr int kEndRegionSlotCount = kEndRegionOtherSlot + 1;
  long long localStatusCounts[kStatusCount] = {};
  std::array<long long, kEndRegionSlotCount> localEndRegionCounts = {};
  long long localStepTotal = 0;
  int localStepMin = std::numeric_limits<int>::max();
  int localStepMax = 0;
  long long localPunctureTotal = 0;
  int localPunctureMin = std::numeric_limits<int>::max();
  int localPunctureMax = 0;
  long long localPunctureCapHits = 0;
  long long localLengthCount = 0;
  double localLengthTotal = 0.0;
  double localLengthMin = std::numeric_limits<double>::infinity();
  double localLengthMax = -std::numeric_limits<double>::infinity();
  long long localSignChangeCandidates = 0;
  long long localRefinementIterations = 0;
  long long localDedupRejects = 0;
  long long localYRejects = 0;
  const bool haveDiagnostics =
      diagnostics != nullptr &&
      diagnostics->signChangeCandidatesPerSeed.size() >=
          traceStatusPerSeed.size() &&
      diagnostics->refinementIterationsPerSeed.size() >=
          traceStatusPerSeed.size() &&
      diagnostics->dedupRejectsPerSeed.size() >= traceStatusPerSeed.size() &&
      diagnostics->yRejectsPerSeed.size() >= traceStatusPerSeed.size();

  for (std::size_t i = 0; i < traceStatusPerSeed.size(); ++i)
  {
    const int statusIndex = traceStatusSummaryIndex(traceStatusPerSeed[i]);
    if (statusIndex >= 0 && statusIndex < kStatusCount)
      ++localStatusCounts[statusIndex];

    const int stateCount =
        (i < stateCountPerSeed.size()) ? std::max(0, stateCountPerSeed[i]) : 0;
    const int stepCount = std::max(0, stateCount - 1);
    localStepTotal += static_cast<long long>(stepCount);
    localStepMin = std::min(localStepMin, stepCount);
    localStepMax = std::max(localStepMax, stepCount);

    const int endRegion = (i < endRegionPerSeed.size()) ? endRegionPerSeed[i] : 0;
    const int endRegionSlot =
        (endRegion >= 0 && endRegion < kEndRegionOtherSlot)
            ? endRegion
            : kEndRegionOtherSlot;
    ++localEndRegionCounts[static_cast<std::size_t>(endRegionSlot)];

    const int punctureCount =
        (i < punctureCountPerSeed.size())
            ? std::max(0, punctureCountPerSeed[i])
            : 0;
    localPunctureTotal += static_cast<long long>(punctureCount);
    localPunctureMin = std::min(localPunctureMin, punctureCount);
    localPunctureMax = std::max(localPunctureMax, punctureCount);
    if (maxPuncPerSeed > 0 && punctureCount >= maxPuncPerSeed)
      ++localPunctureCapHits;

    if (i < connectionLengthPerSeed.size() &&
        std::isfinite(connectionLengthPerSeed[i]))
    {
      ++localLengthCount;
      localLengthTotal += connectionLengthPerSeed[i];
      localLengthMin = std::min(localLengthMin, connectionLengthPerSeed[i]);
      localLengthMax = std::max(localLengthMax, connectionLengthPerSeed[i]);
    }

    if (haveDiagnostics)
    {
      localSignChangeCandidates +=
          std::max(0, diagnostics->signChangeCandidatesPerSeed[i]);
      localRefinementIterations +=
          std::max(0, diagnostics->refinementIterationsPerSeed[i]);
      localDedupRejects += std::max(0, diagnostics->dedupRejectsPerSeed[i]);
      localYRejects += std::max(0, diagnostics->yRejectsPerSeed[i]);
    }
  }

  const long long globalSeedCount =
      mpi.allreduceSumLongLong(static_cast<long long>(traceStatusPerSeed.size()));
  long long globalStatusCounts[kStatusCount] = {};
  for (int i = 0; i < kStatusCount; ++i)
    globalStatusCounts[i] = mpi.allreduceSumLongLong(localStatusCounts[i]);

  std::array<long long, kEndRegionSlotCount> globalEndRegionCounts = {};
  for (int i = 0; i < kEndRegionSlotCount; ++i)
  {
    globalEndRegionCounts[static_cast<std::size_t>(i)] =
        mpi.allreduceSumLongLong(
            localEndRegionCounts[static_cast<std::size_t>(i)]);
  }

  const long long globalStepTotal = mpi.allreduceSumLongLong(localStepTotal);
  const int globalStepMin = mpi.allreduceMinInt(localStepMin);
  const int globalStepMax = mpi.allreduceMaxInt(localStepMax);
  const long long globalPunctureTotal =
      mpi.allreduceSumLongLong(localPunctureTotal);
  const int globalPunctureMin = mpi.allreduceMinInt(localPunctureMin);
  const int globalPunctureMax = mpi.allreduceMaxInt(localPunctureMax);
  const long long globalPunctureCapHits =
      mpi.allreduceSumLongLong(localPunctureCapHits);
  const long long globalLengthCount =
      mpi.allreduceSumLongLong(localLengthCount);
  const double globalLengthTotal = mpi.allreduceSumDouble(localLengthTotal);
  const double globalLengthMin = mpi.allreduceMinDouble(localLengthMin);
  const double globalLengthMax = mpi.allreduceMaxDouble(localLengthMax);
  const long long globalSignChangeCandidates =
      mpi.allreduceSumLongLong(localSignChangeCandidates);
  const long long globalRefinementIterations =
      mpi.allreduceSumLongLong(localRefinementIterations);
  const long long globalDedupRejects =
      mpi.allreduceSumLongLong(localDedupRejects);
  const long long globalYRejects = mpi.allreduceSumLongLong(localYRejects);

  if (mpi.rank != 0)
    return;

  std::cout << "\nTrace result summary:\n";
  std::cout << "  seeds: " << globalSeedCount << "\n";
  for (int i = 0; i < kStatusCount; ++i)
  {
    const double percent =
        (globalSeedCount > 0)
            ? (100.0 * static_cast<double>(globalStatusCounts[i]) /
               static_cast<double>(globalSeedCount))
            : 0.0;
    std::cout << "  status " << traceStatusName(traceStatusFromSummaryIndex(i))
              << ": " << globalStatusCounts[i] << " (" << std::fixed
              << std::setprecision(2) << percent << "%)\n";
  }

  const double averageSteps =
      (globalSeedCount > 0)
          ? (static_cast<double>(globalStepTotal) /
             static_cast<double>(globalSeedCount))
          : 0.0;
  const int printedStepMin = (globalSeedCount > 0) ? globalStepMin : 0;
  std::cout << "  steps/seed: min=" << printedStepMin
            << ", max=" << globalStepMax << ", avg=" << std::fixed
            << std::setprecision(3) << averageSteps
            << ", total=" << globalStepTotal << "\n";

  std::cout << "  end-regions:";
  bool printedRegion = false;
  for (int i = 0; i < kEndRegionSlotCount; ++i)
  {
    const long long count =
        globalEndRegionCounts[static_cast<std::size_t>(i)];
    if (count == 0)
      continue;
    if (i == kEndRegionOtherSlot)
      std::cout << " other=" << count;
    else
      std::cout << " " << i << "=" << count;
    printedRegion = true;
  }
  if (!printedRegion)
    std::cout << " none";
  std::cout << "\n";

  const double averagePunctures =
      (globalSeedCount > 0)
          ? (static_cast<double>(globalPunctureTotal) /
             static_cast<double>(globalSeedCount))
          : 0.0;
  const int printedPunctureMin =
      (globalSeedCount > 0) ? globalPunctureMin : 0;
  std::cout << "  punctures/seed: min=" << printedPunctureMin
            << ", max=" << globalPunctureMax << ", avg=" << std::fixed
            << std::setprecision(3) << averagePunctures
            << ", total=" << globalPunctureTotal << "\n";
  const double capHitPercent =
      (globalSeedCount > 0)
          ? (100.0 * static_cast<double>(globalPunctureCapHits) /
             static_cast<double>(globalSeedCount))
          : 0.0;
  std::cout << "  puncture cap hits: " << globalPunctureCapHits << " ("
            << std::fixed << std::setprecision(2) << capHitPercent << "%)\n";

  if (globalLengthCount > 0)
  {
    const double averageLength = globalLengthTotal /
                                 static_cast<double>(globalLengthCount);
    std::cout << "  length: min=" << std::setprecision(6)
              << globalLengthMin << ", max=" << globalLengthMax
              << ", avg=" << averageLength;
    if (globalLengthCount != globalSeedCount)
      std::cout << " (" << globalLengthCount << " finite)";
    std::cout << "\n";
  }

  if (haveDiagnostics)
  {
    std::cout << "  puncture diagnostics: sign-change-candidates="
              << globalSignChangeCandidates
              << ", accepted=" << globalPunctureTotal
              << ", refinement-iterations=" << globalRefinementIterations;
    if (globalSignChangeCandidates > 0)
    {
      const double averageRefinementIterations =
          static_cast<double>(globalRefinementIterations) /
          static_cast<double>(globalSignChangeCandidates);
      std::cout << ", refine-iter/candidate=" << std::fixed
                << std::setprecision(3) << averageRefinementIterations;
    }
    std::cout << "\n";
    std::cout << "  puncture rejects: dedup=" << globalDedupRejects
              << ", y<=0=" << globalYRejects << "\n";
  }
}

std::string toLowerAscii(std::string text)
{
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c)
                 { return static_cast<char>(std::tolower(c)); });
  return text;
}

bool parseOnOffChoice(const std::string &valueLower, bool &enabled)
{
  if (valueLower == "on" || valueLower == "true" || valueLower == "1")
  {
    enabled = true;
    return true;
  }
  if (valueLower == "off" || valueLower == "false" || valueLower == "0")
  {
    enabled = false;
    return true;
  }
  return false;
}

bool isSupportedViskoresDeviceChoice(const std::string &valueLower)
{
  return valueLower == "serial" || valueLower == "openmp" ||
         valueLower == "kokkos";
}

bool isSupportedTraceEngineChoice(const std::string &valueLower)
{
  return valueLower == "cpu" || valueLower == "viskores";
}

bool isSupportedViskoresOutputModeChoice(const std::string &valueLower)
{
  return valueLower == "punctures" || valueLower == "rk4" ||
         valueLower == "states";
}

bool isSupportedOutputFormatChoice(const std::string &valueLower)
{
  return valueLower == "text" || valueLower == "adios" ||
         valueLower == "both";
}

TraceEngine parseTraceEngineChoice(const std::string &valueLower)
{
  return (valueLower == "viskores") ? TraceEngine::Viskores : TraceEngine::Cpu;
}

OutputFormat parseOutputFormatChoice(const std::string &valueLower)
{
  if (valueLower == "adios")
    return OutputFormat::Adios;
  if (valueLower == "both")
    return OutputFormat::Both;
  return OutputFormat::Text;
}

bool outputFormatWritesText(OutputFormat outputFormat)
{
  return outputFormat == OutputFormat::Text ||
         outputFormat == OutputFormat::Both;
}

bool outputFormatWritesAdios(OutputFormat outputFormat)
{
  return outputFormat == OutputFormat::Adios ||
         outputFormat == OutputFormat::Both;
}

#if defined(CODEX_USE_VISKORES)
ViskoresOutputMode parseViskoresOutputModeChoice(const std::string &valueLower)
{
  if (valueLower == "states")
    return ViskoresOutputMode::States;
  if (valueLower == "rk4")
    return ViskoresOutputMode::Rk4;
  return ViskoresOutputMode::Punctures;
}
#endif

struct ParsedCommandLineOptions
{
  ValidationConfig config;
  bool doCompare = false;
  std::string traceEngineChoice = "viskores";
  std::string viskoresDevice = "serial";
  bool viskoresDeviceSpecified = false;
  std::string viskoresOutputModeChoice = "punctures";
  bool viskoresOutputModeSpecified = false;
  std::string outputFormatChoice;
  std::string outputFile;
  bool outputFileSpecified = false;
  bool psiAxisSpecified = false;
  bool psiBndrySpecified = false;
  double psiAxis = 0.0;
  double psiBndry = 0.0;
  std::string gridFile;
  bool gridFileSpecified = false;
  bool perSeedLog = false;
  LineSpecMode lineSpecMode = LineSpecMode::none;
  std::vector<double> requestedLines;
};

ParsedCommandLineOptions makeDefaultCommandLineOptions()
{
  ParsedCommandLineOptions options;
  options.config.referenceDir =
      "/Users/dpn/proj/bout++/ben_zhu_poincare/zperiod_5";
  options.config.outputDir = "./outputs";
  options.outputFile = "poincare.bp";
  options.outputFormatChoice = defaultOutputFormatChoice();
  options.config.aparSinglePath =
      "/Users/dpn/proj/bout++/ben_zhu_poincare/apar.single.nc";
  options.config.aparCircPath =
      "/Users/dpn/proj/bout++/ben_zhu_poincare/apar.circ.nc";
  options.config.traceOptions.direction = 1;
  options.config.traceOptions.npMax = 100;
  options.config.tolerance = 1.0e-8;
  return options;
}

bool parseCommandLineArguments(int argc, char **argv, const MpiRuntime &mpi,
                               ParsedCommandLineOptions &options)
{
  options = makeDefaultCommandLineOptions();

  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--reference-dir" && i + 1 < argc)
    {
      options.config.referenceDir = argv[++i];
      continue;
    }
    if (arg == "--output-dir" && i + 1 < argc)
    {
      options.config.outputDir = argv[++i];
      continue;
    }
    if (arg == "--text-output-dir" && i + 1 < argc)
    {
      options.config.outputDir = argv[++i];
      continue;
    }
    if (arg == "--output" && i + 1 < argc)
    {
      options.outputFile = argv[++i];
      options.outputFileSpecified = true;
      continue;
    }
    if (arg == "--output-format" && i + 1 < argc)
    {
      options.outputFormatChoice = toLowerAscii(argv[++i]);
      continue;
    }
    if (arg == "--adios-file" && i + 1 < argc)
    {
      options.outputFile = argv[++i];
      options.outputFileSpecified = true;
      continue;
    }
    if (arg == "--psi-axis" && i + 1 < argc)
    {
      options.psiAxis = std::stod(argv[++i]);
      options.psiAxisSpecified = true;
      continue;
    }
    if (arg == "--psi-bndry" && i + 1 < argc)
    {
      options.psiBndry = std::stod(argv[++i]);
      options.psiBndrySpecified = true;
      continue;
    }
    if (arg == "--grid-file" && i + 1 < argc)
    {
      options.gridFile = argv[++i];
      options.gridFileSpecified = true;
      continue;
    }
    if (arg == "--single-apar" && i + 1 < argc)
    {
      options.config.aparSinglePath = argv[++i];
      continue;
    }
    if (arg == "--circ-apar" && i + 1 < argc)
    {
      options.config.aparCircPath = argv[++i];
      continue;
    }
    if (arg == "--divertor" && i + 1 < argc)
    {
      options.config.divertorFilter = argv[++i];
      continue;
    }
    if (arg == "--lines" && i + 1 < argc)
    {
      if (options.lineSpecMode == LineSpecMode::range)
      {
        if (mpi.rank == 0)
          std::cerr << "Error: cannot use both --lines and --nlines\n";
        return false;
      }
      options.lineSpecMode = LineSpecMode::csv;
      options.requestedLines = parseLinesCsvDoubles(argv[++i]);
      continue;
    }
    if (arg == "--nlines" && i + 3 < argc)
    {
      if (options.lineSpecMode == LineSpecMode::csv)
      {
        if (mpi.rank == 0)
          std::cerr << "Error: cannot use both --lines and --nlines\n";
        return false;
      }
      options.lineSpecMode = LineSpecMode::range;
      const double x0 = std::stod(argv[++i]);
      const double x1 = std::stod(argv[++i]);
      const int n = std::atoi(argv[++i]);
      if (n <= 0)
      {
        if (mpi.rank == 0)
          std::cerr << "Error: --nlines requires N > 0\n";
        return false;
      }
      options.requestedLines = buildLineRange(x0, x1, n);
      continue;
    }
    if (arg == "--direction" && i + 1 < argc)
    {
      options.config.traceOptions.direction = std::atoi(argv[++i]);
      continue;
    }
    if (arg == "--np-max" && i + 1 < argc)
    {
      options.config.traceOptions.npMax = std::atoi(argv[++i]);
      continue;
    }
    if (arg == "--max-steps" && i + 1 < argc)
    {
      options.config.traceOptions.maxSteps = std::atoi(argv[++i]);
      continue;
    }
    if (arg == "--tol" && i + 1 < argc)
    {
      options.config.tolerance = std::atof(argv[++i]);
      continue;
    }
    if (arg == "--trace-engine" && i + 1 < argc)
    {
      options.traceEngineChoice = toLowerAscii(argv[++i]);
      continue;
    }
    if (arg == "--viskores-device" && i + 1 < argc)
    {
      options.viskoresDevice = toLowerAscii(argv[++i]);
      options.viskoresDeviceSpecified = true;
      continue;
    }
    if (arg == "--viskores-output" && i + 1 < argc)
    {
      options.viskoresOutputModeChoice = toLowerAscii(argv[++i]);
      options.viskoresOutputModeSpecified = true;
      continue;
    }
    if (arg == "--puncture-detection" && i + 1 < argc)
    {
      bool enabled = true;
      const std::string valueLower = toLowerAscii(argv[++i]);
      if (!parseOnOffChoice(valueLower, enabled))
      {
        if (mpi.rank == 0)
          std::cerr << "Error: --puncture-detection must be on|off\n";
        return false;
      }
      options.config.traceOptions.punctureDetection = enabled;
      continue;
    }
    if (arg == "--puncture-refine" && i + 1 < argc)
    {
      bool enabled = true;
      const std::string valueLower = toLowerAscii(argv[++i]);
      if (!parseOnOffChoice(valueLower, enabled))
      {
        if (mpi.rank == 0)
          std::cerr << "Error: --puncture-refine must be on|off\n";
        return false;
      }
      options.config.traceOptions.punctureRefinement = enabled;
      continue;
    }
    if (arg == "--trace-diagnostics" && i + 1 < argc)
    {
      bool enabled = false;
      const std::string valueLower = toLowerAscii(argv[++i]);
      if (!parseOnOffChoice(valueLower, enabled))
      {
        if (mpi.rank == 0)
          std::cerr << "Error: --trace-diagnostics must be on|off\n";
        return false;
      }
      options.config.traceOptions.traceDiagnostics = enabled;
      continue;
    }
    if (arg == "--compare")
    {
      options.doCompare = true;
      continue;
    }
    if (arg == "--per-seed-log")
    {
      options.perSeedLog = true;
      continue;
    }

    if (mpi.rank == 0)
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
    if (mpi.rank == 0)
      printUsage(argv[0]);
    return false;
  }

  if (options.config.traceOptions.direction != 1 &&
      options.config.traceOptions.direction != -1)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --direction must be 1 or -1\n";
    return false;
  }

  if (options.config.traceOptions.npMax <= 0)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --np-max must be positive\n";
    return false;
  }

  if (options.config.traceOptions.maxSteps < 0)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --max-steps must be non-negative\n";
    return false;
  }

  if (!isSupportedViskoresDeviceChoice(options.viskoresDevice))
  {
    if (mpi.rank == 0)
      std::cerr
          << "Error: --viskores-device must be one of serial|openmp|kokkos\n";
    return false;
  }

  if (!isSupportedTraceEngineChoice(options.traceEngineChoice))
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --trace-engine must be one of cpu|viskores\n";
    return false;
  }

  if (!isSupportedViskoresOutputModeChoice(options.viskoresOutputModeChoice))
  {
    if (mpi.rank == 0)
      std::cerr
          << "Error: --viskores-output must be one of punctures|rk4|states\n";
    return false;
  }

  if (!isSupportedOutputFormatChoice(options.outputFormatChoice))
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --output-format must be one of text|adios|both\n";
    return false;
  }

  if (options.outputFileSpecified &&
      !outputFormatWritesAdios(
          parseOutputFormatChoice(options.outputFormatChoice)))
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --output/--adios-file requires "
                   "--output-format adios|both\n";
    return false;
  }

  if (options.psiAxisSpecified != options.psiBndrySpecified)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --psi-axis and --psi-bndry must be provided "
                   "together\n";
    return false;
  }

  if (options.psiAxisSpecified)
  {
    if (!std::isfinite(options.psiAxis) || !std::isfinite(options.psiBndry))
    {
      if (mpi.rank == 0)
        std::cerr << "Error: --psi-axis and --psi-bndry must be finite\n";
      return false;
    }
    if (options.psiAxis == options.psiBndry)
    {
      if (mpi.rank == 0)
        std::cerr << "Error: --psi-axis and --psi-bndry must differ\n";
      return false;
    }
  }

  return true;
}

bool validatePsiNormalization(const PsiNormalization &normalization,
                              const MpiRuntime &mpi)
{
  if (!normalization.enabled)
    return true;

  if (!std::isfinite(normalization.psiAxis) ||
      !std::isfinite(normalization.psiBndry))
  {
    if (mpi.rank == 0)
      std::cerr << "Error: psi_axis and psi_bndry must be finite\n";
    return false;
  }

  if (normalization.psiAxis == normalization.psiBndry)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: psi_axis and psi_bndry must differ\n";
    return false;
  }

  return true;
}

bool resolvePsiNormalization(const ParsedCommandLineOptions &options,
                             const MpiRuntime &mpi,
                             PsiNormalization &normalization)
{
  if (options.psiAxisSpecified)
  {
    normalization.enabled = true;
    normalization.psiAxis = options.psiAxis;
    normalization.psiBndry = options.psiBndry;
    return validatePsiNormalization(normalization, mpi);
  }

  if (!options.gridFileSpecified)
  {
    normalization = PsiNormalization{};
    return true;
  }

  try
  {
    normalization = readPsiNormalizationFromGridFile(options.gridFile);
  }
  catch (const std::exception &error)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: failed to read psi_axis/psi_bndry from grid file "
                << options.gridFile << ": " << error.what() << "\n";
    return false;
  }

  return validatePsiNormalization(normalization, mpi);
}

const char *traceEngineName(TraceEngine traceEngine)
{
  return (traceEngine == TraceEngine::Viskores) ? "viskores" : "cpu";
}

std::string viskoresPrecisionName()
{
#if defined(CODEX_USE_VISKORES)
  return CODEX_VISKORES_FLOAT_PRECISION_NAME;
#else
  return "none";
#endif
}

void printTraceBackendSelection(const MpiRuntime &mpi, TraceEngine traceEngine,
                                const std::string &viskoresDevice,
                                const std::string &viskoresOutputMode,
                                const TraceOptions &traceOptions)
{
  const bool viskoresActive = (traceEngine == TraceEngine::Viskores);

  mpi.barrier();
  for (int printer = 0; printer < mpi.size; ++printer)
  {
    if (mpi.rank == printer)
    {
      std::cout << "Rank " << mpi.rank << " runtime selection: trace-engine="
                << traceEngineName(traceEngine)
                << ", viskores-device=" << viskoresDevice
                << ", viskores-output=" << viskoresOutputMode
                << ", puncture-detection="
                << (traceOptions.punctureDetection ? "on" : "off")
                << ", puncture-refine="
                << (traceOptions.punctureRefinement ? "on" : "off")
                << ", trace-diagnostics="
                << (traceOptions.traceDiagnostics ? "on" : "off")
                << ", viskores-active=" << (viskoresActive ? "yes" : "no");
      if (!viskoresActive)
        std::cout << " (device ignored unless --trace-engine viskores)";
      std::cout << std::endl;
    }
    mpi.barrier();
  }
}

using SteadyClock = std::chrono::steady_clock;

struct TimingSummary
{
  double avgSeconds = 0.0;
  double maxSeconds = 0.0;
};

double elapsedSeconds(const SteadyClock::time_point &start,
                      const SteadyClock::time_point &end)
{
  return std::chrono::duration<double>(end - start).count();
}

TimingSummary summarizeTiming(const MpiRuntime &mpi, double localSeconds)
{
  TimingSummary summary;
  summary.avgSeconds = mpi.allreduceSumDouble(localSeconds) /
                       static_cast<double>(std::max(1, mpi.size));
  summary.maxSeconds = mpi.allreduceMaxDouble(localSeconds);
  return summary;
}

void printTimingLine(const char *label, const TimingSummary &summary)
{
  std::cout << "  " << label << ": avg=" << summary.avgSeconds
            << " s, max=" << summary.maxSeconds << " s\n";
}

void printTraceTimingSummary(const MpiRuntime &mpi, double localLoadSeconds,
                             double localTraceSeconds,
                             double localOutputSeconds,
                             double localTotalSeconds, const char *traceLabel,
                             double localHostPostprocessSeconds = -1.0)
{
  const TimingSummary load = summarizeTiming(mpi, localLoadSeconds);
  const TimingSummary trace = summarizeTiming(mpi, localTraceSeconds);
  const TimingSummary output = summarizeTiming(mpi, localOutputSeconds);
  const TimingSummary total = summarizeTiming(mpi, localTotalSeconds);
  const bool haveHostPostprocess = (localHostPostprocessSeconds >= 0.0);
  const TimingSummary hostPostprocess =
      haveHostPostprocess ? summarizeTiming(mpi, localHostPostprocessSeconds)
                          : TimingSummary{};

  if (mpi.rank != 0)
    return;

  const std::ios::fmtflags oldFlags = std::cout.flags();
  const std::streamsize oldPrecision = std::cout.precision();
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nTiming summary (trace-only, seconds)\n";
  printTimingLine("load", load);
  printTimingLine(traceLabel, trace);
  if (haveHostPostprocess)
    printTimingLine("host-postprocess", hostPostprocess);
  printTimingLine("output", output);
  printTimingLine("total", total);
  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}
void printCompareTimingSummary(const MpiRuntime &mpi,
                               double localCompareSeconds,
                               double localTotalSeconds)
{
  const TimingSummary compare = summarizeTiming(mpi, localCompareSeconds);
  const TimingSummary total = summarizeTiming(mpi, localTotalSeconds);

  if (mpi.rank != 0)
    return;

  const std::ios::fmtflags oldFlags = std::cout.flags();
  const std::streamsize oldPrecision = std::cout.precision();
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nTiming summary (compare, seconds)\n";
  printTimingLine("compare", compare);
  printTimingLine("total", total);
  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

#if defined(CODEX_USE_VISKORES)
const char *toViskoresDeviceName(const std::string &valueLower)
{
  if (valueLower == "serial")
    return "Serial";
  if (valueLower == "openmp")
    return "OpenMP";
  if (valueLower == "kokkos")
    return "Kokkos";
  return "";
}

const char *onOff(bool value)
{
  return value ? "ON" : "OFF";
}

void printViskoresDeviceDiagnostics(
    const std::string &requestedDevice,
    const viskores::cont::RuntimeDeviceInformation &info,
    viskores::cont::DeviceAdapterId forcedDeviceId)
{
  std::cout << "Viskores device requested: " << requestedDevice << "\n";
  std::cout << "Viskores device forced: " << forcedDeviceId.GetName() << "\n";

  std::cout << "Viskores runtime devices:";
  const char *deviceNames[] = {"Serial", "OpenMP", "Kokkos", "Cuda"};
  for (const char *deviceName : deviceNames)
  {
    const viskores::cont::DeviceAdapterId deviceId = info.GetId(deviceName);
    if (!deviceId.IsValueValid())
      continue;
    std::cout << " " << deviceId.GetName() << "="
              << (info.Exists(deviceId) ? "available" : "unavailable");
  }
  std::cout << "\n";

#if defined(VISKORES_ENABLE_OPENMP)
  const bool viskoresOpenMpEnabled = true;
#else
  const bool viskoresOpenMpEnabled = false;
#endif
#if defined(VISKORES_ENABLE_KOKKOS)
  const bool viskoresKokkosEnabled = true;
#else
  const bool viskoresKokkosEnabled = false;
#endif
#if defined(VISKORES_KOKKOS_CUDA)
  const bool viskoresKokkosCudaEnabled = true;
#else
  const bool viskoresKokkosCudaEnabled = false;
#endif
#if defined(VISKORES_KOKKOS_HIP)
  const bool viskoresKokkosHipEnabled = true;
#else
  const bool viskoresKokkosHipEnabled = false;
#endif

  std::cout << "Viskores build flags:"
            << " OPENMP=" << onOff(viskoresOpenMpEnabled)
            << " KOKKOS=" << onOff(viskoresKokkosEnabled)
            << " KOKKOS_CUDA=" << onOff(viskoresKokkosCudaEnabled)
            << " KOKKOS_HIP=" << onOff(viskoresKokkosHipEnabled) << "\n";
  std::cout << "Viskores trace scalar: " << CODEX_VISKORES_FLOAT_PRECISION_NAME
            << " (" << sizeof(CodeXViskoresFloat)
            << " bytes), FloatDefault=" << sizeof(viskores::FloatDefault)
            << " bytes\n";
}

bool configureViskoresDevice(const std::string &valueLower,
                             std::string &errorMsg, bool printDiagnostics)
{
  const std::string viskoresName = toViskoresDeviceName(valueLower);
  if (viskoresName.empty())
  {
    errorMsg = "Unsupported --viskores-device value '" + valueLower +
               "' (allowed: serial|openmp|kokkos)";
    return false;
  }

  viskores::cont::RuntimeDeviceInformation info;
  const viskores::cont::DeviceAdapterId deviceId = info.GetId(viskoresName);
  if (!deviceId.IsValueValid())
  {
    errorMsg = "Viskores device name '" + viskoresName + "' is not recognized";
    return false;
  }
  if (!info.Exists(deviceId))
  {
    errorMsg = "Viskores device '" + valueLower +
               "' is not available in this build/runtime";
    return false;
  }

  viskores::cont::GetRuntimeDeviceTracker().ForceDevice(deviceId);
  if (printDiagnostics)
    printViskoresDeviceDiagnostics(valueLower, info, deviceId);
  return true;
}

void configureOpenMpForViskoresSerial(const std::string &viskoresDevice)
{
  if (viskoresDevice != "serial")
    return;
  if (std::getenv("OMP_NUM_THREADS") != nullptr)
    return;
#if defined(_WIN32)
  _putenv_s("OMP_NUM_THREADS", "1");
#else
  setenv("OMP_NUM_THREADS", "1", 0);
#endif
}
#endif

std::vector<Point3D> buildSeedPoints(const std::vector<double> &requestedLines,
                                     const AparData &data)
{
  std::vector<Point3D> seeds;
  seeds.reserve(requestedLines.size());

  const double ySeed = static_cast<double>(data.jyomp + 1);
  const double zSeed = data.ziarray.empty() ? 1.0 : data.ziarray.front();
  for (double line : requestedLines)
  {
    Point3D seedInd;
    seedInd.x = line;
    seedInd.y = ySeed;
    seedInd.z = zSeed;
    seeds.push_back(seedInd);
  }
  return seeds;
}

std::vector<TraceTask> buildAllTasks(const std::vector<double> &requestedLines,
                                     bool doSingle, bool doCirc,
                                     const AparData *singleData,
                                     const AparData *circData)
{
  std::vector<TraceTask> allTasks;
  const size_t divertorCount =
      static_cast<size_t>((doSingle ? 1 : 0) + (doCirc ? 1 : 0));
  allTasks.reserve(requestedLines.size() * divertorCount);

  if (doSingle && singleData != nullptr)
  {
    const std::vector<Point3D> singleSeeds =
        buildSeedPoints(requestedLines, *singleData);
    for (const Point3D &seed : singleSeeds)
    {
      TraceTask t;
      t.divertorTag = "single";
      t.seedInd = seed;
      allTasks.push_back(t);
    }
  }

  if (doCirc && circData != nullptr)
  {
    const std::vector<Point3D> circSeeds =
        buildSeedPoints(requestedLines, *circData);
    for (const Point3D &seed : circSeeds)
    {
      TraceTask t;
      t.divertorTag = "circ";
      t.seedInd = seed;
      allTasks.push_back(t);
    }
  }

  return allTasks;
}

size_t partitionBegin(size_t totalTasks, int rank, int nranks)
{
  return (totalTasks * static_cast<size_t>(rank)) / static_cast<size_t>(nranks);
}

size_t partitionEnd(size_t totalTasks, int rank, int nranks)
{
  return (totalTasks * static_cast<size_t>(rank + 1)) /
         static_cast<size_t>(nranks);
}

void traceLocalTasksForDivertor(
    const std::string &tag, const AparData &data,
    const ValidationConfig &config, const std::vector<TraceTask> &localTasks,
    int rank, int maxStatesPerSeed, int maxTrajPerSeed, int maxPuncPerSeed,
    std::vector<double> &ilinePerSeed, std::vector<int> &endRegionPerSeed,
    std::vector<double> &connectionLengthPerSeed,
    std::vector<int> &stateCountPerSeed, std::vector<int> &trajCountPerSeed,
    std::vector<int> &punctureCountPerSeed,
    std::vector<TraceStatus> &traceStatusPerSeed,
    std::vector<TrajectoryState> &states, std::vector<Point3D> &trajectories,
    std::vector<PuncturePoint> &punctures,
    std::vector<std::uint8_t> &punctureValid, bool perSeedLog,
    bool &localFatal, std::string &localFatalMsg)
{
  AparFieldModel model(data);
  FieldLineIntegrator integrator(model, config.traceOptions);

  if (integrator.maxStatesPerSeed() != maxStatesPerSeed ||
      integrator.maxTrajPerSeed() != maxTrajPerSeed ||
      integrator.maxPuncPerSeed() != maxPuncPerSeed)
  {
    localFatal = true;
    localFatalMsg =
        "Integrator max-per-seed caps differ from preallocated array caps";
    return;
  }

  TraceOutputViews outputViews =
      makeTraceOutputViews(states, trajectories, punctures, &punctureValid);

  for (size_t localIdx = 0; localIdx < localTasks.size(); ++localIdx)
  {
    if (localTasks[localIdx].divertorTag != tag)
      continue;

    const TraceStatus traceStatus = integrator.traceLine(
        localTasks[localIdx].seedInd, localIdx, outputViews,
        stateCountPerSeed[localIdx], trajCountPerSeed[localIdx],
        punctureCountPerSeed[localIdx], endRegionPerSeed[localIdx],
        connectionLengthPerSeed[localIdx], ilinePerSeed[localIdx]);
    traceStatusPerSeed[localIdx] = traceStatus;
    if (isFatalTraceStatus(traceStatus))
    {
      localFatal = true;
      localFatalMsg = "traceLine failed on rank " + std::to_string(rank) +
                      " for local seed " + std::to_string(localIdx) +
                      " with status " + traceStatusName(traceStatus);
      return;
    }

    if (perSeedLog)
    {
      std::cout << "Rank " << rank << " traced " << tag << " line "
                << localTasks[localIdx].seedInd.x
                << ": traj=" << trajCountPerSeed[localIdx]
                << ", punctures=" << punctureCountPerSeed[localIdx];
      if (traceStatus != TraceStatus::Ok)
        std::cout << ", status=" << traceStatusName(traceStatus);
      std::cout << std::endl;
    }
  }
}

#if defined(CODEX_USE_VISKORES)
void traceLocalTasksForDivertorViskores(
    const std::string &tag, const AparData &data,
    const ValidationConfig &config, ViskoresOutputMode viskoresOutputMode,
    const std::vector<TraceTask> &localTasks, int rank, int maxStatesPerSeed,
    int maxTrajPerSeed, int maxPuncPerSeed, std::vector<double> &ilinePerSeed,
    std::vector<int> &endRegionPerSeed,
    std::vector<double> &connectionLengthPerSeed,
    std::vector<int> &stateCountPerSeed, std::vector<int> &trajCountPerSeed,
    std::vector<int> &punctureCountPerSeed,
    std::vector<TraceStatus> &traceStatusPerSeed,
    std::vector<TrajectoryState> &states, std::vector<Point3D> &trajectories,
    std::vector<PuncturePoint> &punctures,
    std::vector<std::uint8_t> &punctureValid, TraceDiagnostics *diagnostics,
    bool perSeedLog, double &localDeviceInvokeSeconds,
    double &localHostPostprocessSeconds, bool &localFatal,
    std::string &localFatalMsg)
{
  ViskoresFieldLineTracer tracer(data, config.traceOptions, viskoresOutputMode);

  if (tracer.maxStatesPerSeed() != maxStatesPerSeed ||
      tracer.maxTrajPerSeed() != maxTrajPerSeed ||
      tracer.maxPuncPerSeed() != maxPuncPerSeed)
  {
    localFatal = true;
    localFatalMsg =
        "Viskores tracer max-per-seed caps differ from preallocated array caps";
    return;
  }

  std::vector<Point3D> batchSeeds;
  std::vector<CodeXId> batchGlobalSeedIndices;
  batchSeeds.reserve(localTasks.size());
  batchGlobalSeedIndices.reserve(localTasks.size());

  for (std::size_t localIdx = 0; localIdx < localTasks.size(); ++localIdx)
  {
    if (localTasks[localIdx].divertorTag != tag)
      continue;
    batchSeeds.push_back(localTasks[localIdx].seedInd);
    batchGlobalSeedIndices.push_back(static_cast<CodeXId>(localIdx));
  }

  if (batchSeeds.empty())
    return;

  std::vector<TraceStatus> traceStatuses;
  double batchDeviceInvokeSeconds = 0.0;
  double batchHostPostprocessSeconds = 0.0;
  tracer.traceLines(
      batchSeeds, batchGlobalSeedIndices,
      makeTraceOutputViews(states, trajectories, punctures, &punctureValid),
      ilinePerSeed, endRegionPerSeed, connectionLengthPerSeed,
      stateCountPerSeed, trajCountPerSeed, punctureCountPerSeed, traceStatuses,
      &batchDeviceInvokeSeconds, &batchHostPostprocessSeconds, diagnostics);
  localDeviceInvokeSeconds += batchDeviceInvokeSeconds;
  localHostPostprocessSeconds += batchHostPostprocessSeconds;

  for (std::size_t batchIdx = 0; batchIdx < batchSeeds.size(); ++batchIdx)
  {
    const std::size_t localIdx =
        static_cast<std::size_t>(batchGlobalSeedIndices[batchIdx]);
    const TraceStatus traceStatus = traceStatuses[batchIdx];
    traceStatusPerSeed[localIdx] = traceStatus;
    if (isFatalTraceStatus(traceStatus))
    {
      localFatal = true;
      localFatalMsg = "viskores trace failed on rank " + std::to_string(rank) +
                      " for local seed " + std::to_string(localIdx) +
                      " with status " + traceStatusName(traceStatus);
      return;
    }

    if (perSeedLog)
    {
      std::cout << "Rank " << rank << " traced " << tag << " line "
                << localTasks[localIdx].seedInd.x
                << ": traj=" << trajCountPerSeed[localIdx]
                << ", punctures=" << punctureCountPerSeed[localIdx];
      if (traceStatus != TraceStatus::Ok)
        std::cout << ", status=" << traceStatusName(traceStatus);
      std::cout << std::endl;
    }
  }
}
#endif

} // namespace

int main(int argc, char **argv)
{
  if (hasHelpArgument(argc, argv))
  {
    printUsage(argv[0]);
    return 0;
  }

  MpiRuntime mpi(&argc, &argv);

#if defined(CODEX_USE_MPI)
  std::cout << "Rank " << mpi.rank << " MPI initialized (size=" << mpi.size
            << ")" << std::endl;
#endif

  ParsedCommandLineOptions options;
  if (!parseCommandLineArguments(argc, argv, mpi, options))
    return 1;

  ValidationConfig config = options.config;
  const bool doCompare = options.doCompare;
  const std::string &viskoresDevice = options.viskoresDevice;
  const bool viskoresDeviceSpecified = options.viskoresDeviceSpecified;
  const std::string &viskoresOutputModeChoice =
      options.viskoresOutputModeChoice;
  const bool viskoresOutputModeSpecified = options.viskoresOutputModeSpecified;
  const LineSpecMode lineSpecMode = options.lineSpecMode;
  const std::vector<double> &requestedLines = options.requestedLines;
  const TraceEngine traceEngine =
      parseTraceEngineChoice(options.traceEngineChoice);
  const OutputFormat requestedOutputFormat =
      parseOutputFormatChoice(options.outputFormatChoice);
  const OutputFormat outputFormat =
      doCompare ? OutputFormat::Text : requestedOutputFormat;
  const bool writeTextOutput = outputFormatWritesText(outputFormat);
  const bool writeAdiosOutput = outputFormatWritesAdios(outputFormat);
  const std::string outputFile = options.outputFile;
  const bool perSeedLog = options.perSeedLog;

#if !defined(CODEX_USE_ADIOS)
  if (writeAdiosOutput)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: ADIOS output was requested, but this binary was "
                   "built without CODEX_USE_ADIOS=ON\n";
    return 1;
  }
#endif
  if (writeAdiosOutput && config.divertorFilter != "single" &&
      config.divertorFilter != "circ")
  {
    if (mpi.rank == 0)
      std::cerr << "Error: ADIOS output requires --divertor single|circ\n";
    return 1;
  }

  PsiNormalization psiNormalization;
  if (writeAdiosOutput &&
      !resolvePsiNormalization(options, mpi, psiNormalization))
    return 1;

#if defined(CODEX_USE_VISKORES)
  const ViskoresOutputMode viskoresOutputMode =
      parseViskoresOutputModeChoice(viskoresOutputModeChoice);
  if (traceEngine == TraceEngine::Viskores)
  {
    configureOpenMpForViskoresSerial(viskoresDevice);
    viskores::cont::Initialize(argc, argv);
    std::string viskoresDeviceError;
    if (!configureViskoresDevice(viskoresDevice, viskoresDeviceError,
                                 mpi.rank == 0))
    {
      if (mpi.rank == 0)
        std::cerr << "Error: " << viskoresDeviceError << "\n";
      return 1;
    }
  }
#else
  if (traceEngine == TraceEngine::Viskores)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --trace-engine viskores was requested, but this "
                   "binary was built without Viskores support\n";
    return 1;
  }
  if (viskoresDeviceSpecified)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --viskores-device was specified, but this binary "
                   "was built without Viskores support\n";
    return 1;
  }
  if (viskoresOutputModeSpecified)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --viskores-output was specified, but this binary "
                   "was built without Viskores support\n";
    return 1;
  }
#endif

  const bool usesPunctureDiagnosticShortcut =
      !config.traceOptions.punctureDetection ||
      !config.traceOptions.punctureRefinement;
  if (usesPunctureDiagnosticShortcut && traceEngine != TraceEngine::Viskores)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --puncture-detection off / "
                   "--puncture-refine off are only supported with "
                   "--trace-engine viskores\n";
    return 1;
  }
  if (usesPunctureDiagnosticShortcut &&
      viskoresOutputModeChoice != "punctures")
  {
    if (mpi.rank == 0)
      std::cerr << "Error: puncture diagnostic shortcuts require "
                   "--viskores-output punctures\n";
    return 1;
  }
  if (config.traceOptions.traceDiagnostics && traceEngine != TraceEngine::Viskores)
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --trace-diagnostics on requires "
                   "--trace-engine viskores\n";
    return 1;
  }
  if (config.traceOptions.traceDiagnostics &&
      viskoresOutputModeChoice != "punctures")
  {
    if (mpi.rank == 0)
      std::cerr << "Error: --trace-diagnostics on requires "
                   "--viskores-output punctures\n";
    return 1;
  }

  try
  {
    const SteadyClock::time_point totalStart = SteadyClock::now();
    printTraceBackendSelection(mpi, traceEngine, viskoresDevice,
                               viskoresOutputModeChoice,
                               config.traceOptions);

    if (doCompare && mpi.size > 1)
    {
      if (mpi.rank == 0)
        std::cerr << "Error: --compare is currently only supported with a "
                     "single MPI rank\n";
      return 1;
    }

    if (doCompare)
    {
      if (traceEngine == TraceEngine::Viskores)
      {
        if (mpi.rank == 0)
          std::cerr << "Error: --compare currently uses the CPU tracer only; "
                       "run trace-only mode for viskores output comparison\n";
        return 1;
      }

      if (lineSpecMode == LineSpecMode::range)
      {
        if (mpi.rank == 0)
          std::cerr << "Error: --nlines is only supported in trace-only mode "
                       "(no --compare)\n";
        return 1;
      }

      if (lineSpecMode == LineSpecMode::csv)
      {
        config.linesFilter.clear();
        config.linesFilter.reserve(requestedLines.size());
        for (double line : requestedLines)
        {
          if (!isIntegerVal(line))
          {
            if (mpi.rank == 0)
              std::cerr
                  << "Error: --compare requires integer line indices; got "
                  << line << "\n";
            return 1;
          }
          config.linesFilter.push_back(static_cast<int>(std::llround(line)));
        }
      }

      const SteadyClock::time_point compareStart = SteadyClock::now();
      ValidationSuite suite;
      const std::vector<ValidationResult> results = suite.run(config);
      const SteadyClock::time_point compareEnd = SteadyClock::now();

      int passed = 0;
      for (const auto &r : results)
        if (r.pass)
          ++passed;

      if (mpi.rank == 0)
        std::cout << "\nValidation summary: " << passed << "/" << results.size()
                  << " cases passed\n";
      printCompareTimingSummary(mpi, elapsedSeconds(compareStart, compareEnd),
                                elapsedSeconds(totalStart, compareEnd));
      return (passed == static_cast<int>(results.size())) ? 0 : 2;
    }

    if (requestedLines.empty())
    {
      if (mpi.rank == 0)
        std::cerr << "Error: trace-only mode requires --lines or --nlines\n";
      return 1;
    }

    const bool doSingle =
        (config.divertorFilter == "all" || config.divertorFilter == "single");
    const bool doCirc =
        (config.divertorFilter == "all" || config.divertorFilter == "circ");
    const bool useCombinedOutput = (lineSpecMode == LineSpecMode::range);

    if (!doSingle && !doCirc)
    {
      if (mpi.rank == 0)
        std::cerr << "Error: --divertor must be all|single|circ\n";
      return 1;
    }

    const SteadyClock::time_point loadStart = SteadyClock::now();
    AparData singleData;
    AparData circData;
    if (doSingle)
      singleData.load(config.aparSinglePath);
    if (doCirc)
      circData.load(config.aparCircPath);
    const SteadyClock::time_point loadEnd = SteadyClock::now();
    const double localLoadSeconds = elapsedSeconds(loadStart, loadEnd);

    const SteadyClock::time_point traceStart = SteadyClock::now();
    const std::vector<TraceTask> allTasks = buildAllTasks(
        requestedLines, doSingle, doCirc, doSingle ? &singleData : nullptr,
        doCirc ? &circData : nullptr);
    const size_t localBegin =
        partitionBegin(allTasks.size(), mpi.rank, mpi.size);
    const size_t localEnd = partitionEnd(allTasks.size(), mpi.rank, mpi.size);
    std::vector<TraceTask> localTasks(
        allTasks.begin() + static_cast<std::ptrdiff_t>(localBegin),
        allTasks.begin() + static_cast<std::ptrdiff_t>(localEnd));

    if (perSeedLog)
    {
      mpi.barrier();
      for (int printer = 0; printer < mpi.size; ++printer)
      {
        if (mpi.rank == printer)
        {
          std::cout << "Rank " << mpi.rank << " seed IDs:";
          if (localBegin == localEnd)
            std::cout << " (none)";
          else
            for (size_t seedId = localBegin; seedId < localEnd; ++seedId)
              std::cout << " " << seedId;

          std::cout << std::endl;
        }
        mpi.barrier();
      }
    }

    const bool compactViskoresOutput =
        (traceEngine == TraceEngine::Viskores &&
         (viskoresOutputModeChoice == "punctures" ||
          viskoresOutputModeChoice == "rk4"));
    const int traceMaxStatesPerSeed = computeMaxStateCount(config.traceOptions);
    const int maxStatesPerSeed =
        compactViskoresOutput ? 1 : traceMaxStatesPerSeed;
    const int maxTrajPerSeed = compactViskoresOutput ? 1 : maxStatesPerSeed;
    const int maxPuncPerSeed = std::max(1, config.traceOptions.npMax);

    const size_t localTaskCount = localTasks.size();
    std::vector<double> ilinePerSeed(localTaskCount, 0.0);
    std::vector<int> endRegionPerSeed(localTaskCount, 0);
    std::vector<double> connectionLengthPerSeed(localTaskCount, 0.0);
    std::vector<int> stateCountPerSeed(localTaskCount, 0);
    std::vector<int> trajCountPerSeed(localTaskCount, 0);
    std::vector<int> punctureCountPerSeed(localTaskCount, 0);
    std::vector<TraceStatus> traceStatusPerSeed(localTaskCount,
                                                TraceStatus::Ok);
    std::vector<TrajectoryState> states(localTaskCount *
                                        static_cast<size_t>(maxStatesPerSeed));
    std::vector<Point3D> trajectories(localTaskCount *
                                      static_cast<size_t>(maxTrajPerSeed));
    std::vector<PuncturePoint> punctures(localTaskCount *
                                         static_cast<size_t>(maxPuncPerSeed));
    std::vector<std::uint8_t> punctureValid(
        localTaskCount * static_cast<size_t>(maxPuncPerSeed),
        static_cast<std::uint8_t>(0));
    TraceDiagnostics traceDiagnostics;
    traceDiagnostics.signChangeCandidatesPerSeed.assign(localTaskCount, 0);
    traceDiagnostics.refinementIterationsPerSeed.assign(localTaskCount, 0);
    traceDiagnostics.dedupRejectsPerSeed.assign(localTaskCount, 0);
    traceDiagnostics.yRejectsPerSeed.assign(localTaskCount, 0);
    TraceDiagnostics *traceDiagnosticsForSummary =
        (traceEngine == TraceEngine::Viskores &&
         viskoresOutputModeChoice == "punctures" &&
         config.traceOptions.traceDiagnostics)
            ? &traceDiagnostics
            : nullptr;

    bool localFatal = false;
    std::string localFatalMsg;
    double localDeviceInvokeSeconds = 0.0;
    double localHostPostprocessSeconds = 0.0;

    if (doSingle && !localFatal && traceEngine != TraceEngine::Viskores)
      traceLocalTasksForDivertor(
          "single", singleData, config, localTasks, mpi.rank, maxStatesPerSeed,
          maxTrajPerSeed, maxPuncPerSeed, ilinePerSeed, endRegionPerSeed,
          connectionLengthPerSeed, stateCountPerSeed, trajCountPerSeed,
          punctureCountPerSeed, traceStatusPerSeed, states, trajectories,
          punctures, punctureValid, perSeedLog, localFatal, localFatalMsg);
#if defined(CODEX_USE_VISKORES)
    if (doSingle && !localFatal && traceEngine == TraceEngine::Viskores)
      traceLocalTasksForDivertorViskores(
          "single", singleData, config, viskoresOutputMode, localTasks,
          mpi.rank, maxStatesPerSeed, maxTrajPerSeed, maxPuncPerSeed,
          ilinePerSeed, endRegionPerSeed, connectionLengthPerSeed,
          stateCountPerSeed, trajCountPerSeed, punctureCountPerSeed,
          traceStatusPerSeed, states, trajectories, punctures, punctureValid,
          traceDiagnosticsForSummary, perSeedLog, localDeviceInvokeSeconds,
          localHostPostprocessSeconds, localFatal, localFatalMsg);
#endif

    if (doCirc && !localFatal && traceEngine != TraceEngine::Viskores)
      traceLocalTasksForDivertor(
          "circ", circData, config, localTasks, mpi.rank, maxStatesPerSeed,
          maxTrajPerSeed, maxPuncPerSeed, ilinePerSeed, endRegionPerSeed,
          connectionLengthPerSeed, stateCountPerSeed, trajCountPerSeed,
          punctureCountPerSeed, traceStatusPerSeed, states, trajectories,
          punctures, punctureValid, perSeedLog, localFatal, localFatalMsg);
#if defined(CODEX_USE_VISKORES)
    if (doCirc && !localFatal && traceEngine == TraceEngine::Viskores)
      traceLocalTasksForDivertorViskores(
          "circ", circData, config, viskoresOutputMode, localTasks, mpi.rank,
          maxStatesPerSeed, maxTrajPerSeed, maxPuncPerSeed, ilinePerSeed,
          endRegionPerSeed, connectionLengthPerSeed, stateCountPerSeed,
          trajCountPerSeed, punctureCountPerSeed, traceStatusPerSeed, states,
          trajectories, punctures, punctureValid, traceDiagnosticsForSummary,
          perSeedLog, localDeviceInvokeSeconds, localHostPostprocessSeconds,
          localFatal, localFatalMsg);
#endif

    const int anyFatal = mpi.allreduceMaxInt(localFatal ? 1 : 0);
    if (anyFatal != 0)
    {
      if (localFatal)
        std::cerr << "Error: " << localFatalMsg << "\n";
      return 1;
    }

    const SteadyClock::time_point traceEnd = SteadyClock::now();
    const double localTraceSeconds = (traceEngine == TraceEngine::Viskores)
                                         ? localDeviceInvokeSeconds
                                         : elapsedSeconds(traceStart, traceEnd);
    const double localReportedHostPostprocessSeconds =
        (traceEngine == TraceEngine::Viskores) ? localHostPostprocessSeconds
                                               : -1.0;
    const char *traceTimingLabel =
        (traceEngine == TraceEngine::Viskores) ? "device-invoke" : "trace";
    printTraceResultSummary(mpi, traceStatusPerSeed, stateCountPerSeed,
                            endRegionPerSeed, punctureCountPerSeed,
                            connectionLengthPerSeed, maxPuncPerSeed,
                            traceDiagnosticsForSummary);

    const SteadyClock::time_point outputStart = SteadyClock::now();
    PoincareOutput output;
    mpi.barrier();

    if (writeAdiosOutput)
    {
#if defined(CODEX_USE_ADIOS)
      AdiosPoincareMetadata adiosMetadata;
      adiosMetadata.divertor = config.divertorFilter;
      adiosMetadata.traceEngine = traceEngineName(traceEngine);
      adiosMetadata.viskoresDevice = viskoresDevice;
      adiosMetadata.viskoresOutputMode = viskoresOutputModeChoice;
      adiosMetadata.viskoresPrecision = viskoresPrecisionName();
      adiosMetadata.hasPsiNormalization = psiNormalization.enabled;
      adiosMetadata.psiAxis = psiNormalization.psiAxis;
      adiosMetadata.psiBndry = psiNormalization.psiBndry;

      AdiosPoincareOutput adiosOutput;
      adiosOutput.writeFlatOutputs(
          ilinePerSeed, endRegionPerSeed, connectionLengthPerSeed,
          punctureCountPerSeed, localBegin, allTasks.size(), maxPuncPerSeed,
          punctures, &punctureValid, outputFile, adiosMetadata);
#endif
    }

    if (writeTextOutput)
    {
      for (int writer = 0; writer < mpi.size; ++writer)
      {
        if (mpi.rank == writer && useCombinedOutput)
          output.writeCombinedOutputsFlat(
              ilinePerSeed, trajCountPerSeed, punctureCountPerSeed,
              maxTrajPerSeed, maxPuncPerSeed, trajectories, punctures,
              &punctureValid, config.outputDir, writer != 0, writer == 0);
        if (mpi.rank == writer && !useCombinedOutput)
          for (size_t i = 0; i < localTasks.size(); ++i)
            output.writeLineOutputsFlat(
                ilinePerSeed[i], i, maxTrajPerSeed, maxPuncPerSeed,
                trajCountPerSeed[i], punctureCountPerSeed[i], trajectories,
                punctures, &punctureValid, config.outputDir,
                localTasks[i].divertorTag);

        mpi.barrier();
      }
    }
    else
    {
      mpi.barrier();
    }
    const SteadyClock::time_point outputEnd = SteadyClock::now();
    const double localOutputSeconds = elapsedSeconds(outputStart, outputEnd);
    const double localTotalSeconds = elapsedSeconds(totalStart, outputEnd);

    if (mpi.rank == 0 && writeTextOutput && useCombinedOutput)
      std::cout << "Wrote combined outputs: " << config.outputDir
                << "/ip_cxx.txt, " << config.outputDir << "/ip_cxx.TP.txt"
                << " and " << config.outputDir << "/traj_cxx.txt\n";
    if (mpi.rank == 0 && writeTextOutput && !useCombinedOutput)
      std::cout << "Wrote per-line outputs in rank order\n";
    if (mpi.rank == 0 && writeAdiosOutput)
      std::cout << "Wrote ADIOS output: " << outputFile << "\n";

    printTraceTimingSummary(mpi, localLoadSeconds, localTraceSeconds,
                            localOutputSeconds, localTotalSeconds,
                            traceTimingLabel,
                            localReportedHostPostprocessSeconds);

    if (mpi.rank == 0)
      std::cout << "\nTrace-only run complete.\n";
    return 0;
  }
  catch (const std::exception &ex)
  {
    std::cerr << "Error: " << ex.what() << "\n";
#if defined(CODEX_USE_MPI)
    if (mpi.size > 1)
      mpi.abortAll(1);
#endif
    return 1;
  }
}
