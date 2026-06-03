#ifndef CODEX_CXX2_ADIOSPOINCAREOUTPUT_H
#define CODEX_CXX2_ADIOSPOINCAREOUTPUT_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "Types.h"

struct AdiosPoincareMetadata
{
  std::string divertor;
  std::string traceEngine;
  std::string viskoresDevice;
  std::string viskoresOutputMode;
  std::string viskoresPrecision;
  bool hasPsiNormalization = false;
  double psiAxis = 0.0;
  double psiBndry = 0.0;
};

class AdiosPoincareOutput
{
public:
  void writeFlatOutputs(const std::vector<double> &ilinePerSeed,
                        const std::vector<int> &endRegionPerSeed,
                        const std::vector<double> &connectionLengthPerSeed,
                        const std::vector<int> &punctureCountPerSeed,
                        std::size_t globalSeedBegin,
                        std::size_t globalSeedCount,
                        int maxPuncPerSeed,
                        const std::vector<PuncturePoint> &punctures,
                        const std::vector<std::uint8_t> *punctureValid,
                        const std::string &outputPath,
                        const AdiosPoincareMetadata &metadata) const;
};

#endif
