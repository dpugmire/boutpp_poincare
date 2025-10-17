// DebugIO.h - minimal compile-time switch for trace file output
#pragma once

#ifndef TRACE_IO
#define TRACE_IO 0
#endif

#if TRACE_IO
  #include <fstream>
  #include <filesystem>
  inline std::filesystem::path traceBaseDir{"."};
  // Use: TRACE( stream << "..." << '\n' );
  #define TRACE(stmt) do { stmt; } while(0)
#else
  #define TRACE(stmt) do {} while(0)
#endif
