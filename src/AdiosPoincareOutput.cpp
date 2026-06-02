#include "AdiosPoincareOutput.h"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

#include <adios2.h>

#if defined(CODEX_USE_MPI)
#include <mpi.h>
#endif

namespace
{

std::uint64_t allreduceSumUInt64(std::uint64_t localValue)
{
#if defined(CODEX_USE_MPI)
  unsigned long long local = static_cast<unsigned long long>(localValue);
  unsigned long long global = local;
  MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  return static_cast<std::uint64_t>(global);
#else
  return localValue;
#endif
}

std::uint64_t exclusiveScanSumUInt64(std::uint64_t localValue)
{
#if defined(CODEX_USE_MPI)
  unsigned long long local = static_cast<unsigned long long>(localValue);
  unsigned long long prefix = 0;
  MPI_Exscan(&local, &prefix, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
             MPI_COMM_WORLD);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    prefix = 0;
  return static_cast<std::uint64_t>(prefix);
#else
  (void)localValue;
  return 0;
#endif
}

template <typename T>
adios2::Variable<T> defineVector(adios2::IO &io, const std::string &name,
                                 std::size_t globalCount,
                                 std::size_t localStart,
                                 std::size_t localCount)
{
  return io.DefineVariable<T>(name, adios2::Dims{globalCount},
                              adios2::Dims{localStart},
                              adios2::Dims{localCount});
}

template <typename T>
void putVector(adios2::Engine &engine, const adios2::Variable<T> &variable,
               const std::vector<T> &values)
{
  if (values.empty())
    return;
  engine.Put(variable, values.data(), adios2::Mode::Sync);
}

} // namespace

void AdiosPoincareOutput::writeFlatOutputs(
    const std::vector<double> &ilinePerSeed,
    const std::vector<int> &endRegionPerSeed,
    const std::vector<double> &connectionLengthPerSeed,
    const std::vector<int> &punctureCountPerSeed,
    std::size_t globalSeedBegin, std::size_t globalSeedCount,
    int maxPuncPerSeed, const std::vector<PuncturePoint> &punctures,
    const std::vector<std::uint8_t> *punctureValid,
    const std::string &outputPath,
    const AdiosPoincareMetadata &metadata) const
{
  const std::size_t localSeedCount = ilinePerSeed.size();
  if (endRegionPerSeed.size() != localSeedCount ||
      connectionLengthPerSeed.size() != localSeedCount ||
      punctureCountPerSeed.size() != localSeedCount)
  {
    throw std::runtime_error(
        "ADIOS output seed metadata arrays must all have the same size");
  }
  if (maxPuncPerSeed <= 0)
    throw std::runtime_error("ADIOS output requires maxPuncPerSeed > 0");

  const std::filesystem::path path(outputPath);
  if (!path.parent_path().empty())
    std::filesystem::create_directories(path.parent_path());

  std::vector<std::uint64_t> seedId(localSeedCount, 0);
  std::vector<std::int32_t> endRegion(localSeedCount, 0);
  std::vector<std::uint64_t> offset(localSeedCount, 0);
  std::vector<std::uint64_t> count(localSeedCount, 0);

  std::vector<std::int32_t> traceStep;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> theta;
  std::vector<double> psi;

  std::uint64_t localPunctureOffset = 0;
  for (std::size_t seed = 0; seed < localSeedCount; ++seed)
  {
    seedId[seed] = static_cast<std::uint64_t>(globalSeedBegin + seed);
    endRegion[seed] = static_cast<std::int32_t>(endRegionPerSeed[seed]);
    offset[seed] = localPunctureOffset;

    const int clampedPunctureCount =
        std::max(0, std::min(punctureCountPerSeed[seed], maxPuncPerSeed));
    const std::size_t punctureBase =
        seed * static_cast<std::size_t>(maxPuncPerSeed);
    if (punctureBase + static_cast<std::size_t>(clampedPunctureCount) >
        punctures.size())
    {
      throw std::runtime_error(
          "ADIOS output puncture array is smaller than numSeeds * maxPuncPerSeed");
    }
    if (punctureValid != nullptr &&
        punctureBase + static_cast<std::size_t>(clampedPunctureCount) >
            punctureValid->size())
    {
      throw std::runtime_error(
          "ADIOS output puncture valid-mask array is smaller than numSeeds * maxPuncPerSeed");
    }

    std::uint64_t validPuncturesForSeed = 0;
    for (int i = 0; i < clampedPunctureCount; ++i)
    {
      const std::size_t punctureArrayIndex =
          punctureBase + static_cast<std::size_t>(i);
      if (punctureValid != nullptr && (*punctureValid)[punctureArrayIndex] == 0)
        continue;

      const PuncturePoint &p = punctures[punctureArrayIndex];
      traceStep.push_back(static_cast<std::int32_t>(p.step));
      x.push_back(p.xyz.x);
      y.push_back(p.xyz.y);
      z.push_back(p.xyz.z);
      theta.push_back(p.thetaPsi.x);
      psi.push_back(p.thetaPsi.y);
      ++validPuncturesForSeed;
    }

    count[seed] = validPuncturesForSeed;
    localPunctureOffset += validPuncturesForSeed;
  }

  const std::uint64_t globalPunctureBegin =
      exclusiveScanSumUInt64(localPunctureOffset);
  const std::uint64_t globalPunctureCount =
      allreduceSumUInt64(localPunctureOffset);
  for (std::uint64_t &value : offset)
    value += globalPunctureBegin;

#if defined(CODEX_USE_MPI)
  adios2::ADIOS adios(MPI_COMM_WORLD);
#else
  adios2::ADIOS adios;
#endif
  adios2::IO io = adios.DeclareIO("PoincareOutput");

  io.DefineAttribute<std::uint64_t>("schema_version", 1);
  io.DefineAttribute<std::uint64_t>(
      "num_seeds", static_cast<std::uint64_t>(globalSeedCount));
  io.DefineAttribute<std::uint64_t>("num_punctures", globalPunctureCount);
  io.DefineAttribute<std::string>("divertor", metadata.divertor);
  io.DefineAttribute<std::string>("trace_engine", metadata.traceEngine);
  io.DefineAttribute<std::string>("viskores_device", metadata.viskoresDevice);
  io.DefineAttribute<std::string>("viskores_output_mode",
                                  metadata.viskoresOutputMode);
  io.DefineAttribute<std::string>("viskores_precision",
                                  metadata.viskoresPrecision);

  const std::size_t seedStart = globalSeedBegin;
  const auto seedIdVar =
      defineVector<std::uint64_t>(io, "seed_id", globalSeedCount, seedStart,
                                  localSeedCount);
  const auto ilineVar =
      defineVector<double>(io, "iline", globalSeedCount, seedStart,
                           localSeedCount);
  const auto endRegionVar =
      defineVector<std::int32_t>(io, "end_region", globalSeedCount, seedStart,
                                 localSeedCount);
  const auto lengthVar =
      defineVector<double>(io, "length", globalSeedCount, seedStart,
                           localSeedCount);
  const auto offsetVar = defineVector<std::uint64_t>(
      io, "offset", globalSeedCount, seedStart, localSeedCount);
  const auto countVar = defineVector<std::uint64_t>(
      io, "count", globalSeedCount, seedStart, localSeedCount);

  adios2::Variable<std::int32_t> traceStepVar;
  adios2::Variable<double> xVar;
  adios2::Variable<double> yVar;
  adios2::Variable<double> zVar;
  adios2::Variable<double> thetaVar;
  adios2::Variable<double> psiVar;

  if (globalPunctureCount > 0)
  {
    const std::size_t punctureStart =
        static_cast<std::size_t>(globalPunctureBegin);
    const std::size_t localPunctureCount = traceStep.size();
    const std::size_t totalPunctureCount =
        static_cast<std::size_t>(globalPunctureCount);
    traceStepVar =
        defineVector<std::int32_t>(io, "trace_step", totalPunctureCount,
                                   punctureStart, localPunctureCount);
    xVar = defineVector<double>(io, "x", totalPunctureCount, punctureStart,
                                localPunctureCount);
    yVar = defineVector<double>(io, "y", totalPunctureCount, punctureStart,
                                localPunctureCount);
    zVar = defineVector<double>(io, "z", totalPunctureCount, punctureStart,
                                localPunctureCount);
    thetaVar =
        defineVector<double>(io, "theta", totalPunctureCount, punctureStart,
                             localPunctureCount);
    psiVar = defineVector<double>(io, "psi", totalPunctureCount,
                                  punctureStart, localPunctureCount);
  }

  adios2::Engine engine = io.Open(outputPath, adios2::Mode::Write);
  putVector(engine, seedIdVar, seedId);
  putVector(engine, ilineVar, ilinePerSeed);
  putVector(engine, endRegionVar, endRegion);
  putVector(engine, lengthVar, connectionLengthPerSeed);
  putVector(engine, offsetVar, offset);
  putVector(engine, countVar, count);

  if (globalPunctureCount > 0)
  {
    putVector(engine, traceStepVar, traceStep);
    putVector(engine, xVar, x);
    putVector(engine, yVar, y);
    putVector(engine, zVar, z);
    putVector(engine, thetaVar, theta);
    putVector(engine, psiVar, psi);
  }
  engine.Close();
}
