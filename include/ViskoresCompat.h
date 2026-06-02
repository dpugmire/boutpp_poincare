#ifndef CODEX_CXX2_VISKORESCOMPAT_H
#define CODEX_CXX2_VISKORESCOMPAT_H

#include <cstddef>

#if defined(CODEX_USE_VISKORES)
#include <viskores/Types.h>
#define CODEX_EXEC VISKORES_EXEC
#define CODEX_CONT VISKORES_CONT
using CodeXId = viskores::Id;

#if defined(CODEX_VISKORES_FLOAT_PRECISION_FLOAT32)
using CodeXViskoresFloat = viskores::Float32;
#define CODEX_VISKORES_FLOAT_PRECISION_NAME "float32"
static_assert(sizeof(CodeXViskoresFloat) == 4, "float32 trace precision must be 4 bytes");
#elif defined(CODEX_VISKORES_FLOAT_PRECISION_FLOAT64)
using CodeXViskoresFloat = viskores::Float64;
#define CODEX_VISKORES_FLOAT_PRECISION_NAME "float64"
static_assert(sizeof(CodeXViskoresFloat) == 8, "float64 trace precision must be 8 bytes");
#else
using CodeXViskoresFloat = viskores::FloatDefault;
#define CODEX_VISKORES_FLOAT_PRECISION_NAME "viskores"
#endif
#else
#define CODEX_EXEC
#define CODEX_CONT
using CodeXId = std::size_t;
#endif

#endif
