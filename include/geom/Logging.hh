#ifndef LOGGING_HH
#define LOGGING_HH

#include "spdlog/spdlog.h"

void CreateLogger(bool to_file, int level, int flush_time);

namespace NuGeom {

// Library-scoped logger named "nugeom".  All NuGeometry library code logs
// through this instead of the spdlog default logger, so that when NuGeometry
// is embedded alongside another spdlog user (e.g. Achilles, which logs via
// the default logger), output can be told apart with the "%n" pattern field.
//
// If the host application registered a "nugeom" logger (CreateLogger or a
// driver), that one is used; otherwise the default logger's sinks, level,
// and pattern are cloned on first use.
//
// The resolved logger is cached, since Log() is called on hot paths. Any code
// that registers or replaces the "nugeom" logger (via spdlog directly) must
// call RefreshLogger() afterwards so Log() returns the new instance.
spdlog::logger &Log();

// Re-resolve and cache the "nugeom" logger. Call after registering/replacing
// it. CreateLogger() calls this for you.
void RefreshLogger();

} // namespace NuGeom

#endif
