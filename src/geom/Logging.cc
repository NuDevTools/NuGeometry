#include "geom/Logging.hh"

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

void CreateLogger(bool to_file, int level, int flush_time) {
    auto slevel = static_cast<spdlog::level::level_enum>(level);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(slevel);

    std::shared_ptr<spdlog::logger> logger;
    if(to_file) {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("geom.log", true);
        file_sink->set_level(slevel);

        spdlog::sinks_init_list sink_list{file_sink, console_sink};
        logger = std::make_shared<spdlog::logger>("nugeom", sink_list);
    } else {
        spdlog::sinks_init_list sink_list{console_sink};
        logger = std::make_shared<spdlog::logger>("nugeom", sink_list);
    }
    logger->set_level(slevel);
    logger->flush_on(spdlog::level::warn);
    spdlog::drop("nugeom"); // replace any previously registered instance
    spdlog::register_logger(logger);
    spdlog::set_default_logger(logger);
    spdlog::flush_every(std::chrono::seconds(flush_time));
    NuGeom::RefreshLogger(); // point the cached handle at the new logger
}

namespace {

// Resolve the "nugeom" logger: prefer whatever the host registered (driver,
// tests, CreateLogger); otherwise clone the default logger's sinks/level/
// pattern once and keep that fallback alive even if the registry is cleared.
std::shared_ptr<spdlog::logger> ResolveLogger() {
    if(auto registered = spdlog::get("nugeom")) return registered;
    static std::shared_ptr<spdlog::logger> fallback = []() {
        auto created = spdlog::default_logger()->clone("nugeom");
        spdlog::register_logger(created);
        return created;
    }();
    return fallback;
}

// Log() runs on hot paths (per element, per segment, per ray), so the
// resolved logger is cached to avoid the registry lookup (spdlog::get takes a
// mutex + hashes the name) on every call.  Hosts that (re)register a "nugeom"
// logger must call NuGeom::RefreshLogger() to update this handle; all in-repo
// registration sites do.
std::shared_ptr<spdlog::logger> &CachedLogger() {
    static std::shared_ptr<spdlog::logger> cached = ResolveLogger();
    return cached;
}

} // namespace

void NuGeom::RefreshLogger() {
    CachedLogger() = ResolveLogger();
}

spdlog::logger &NuGeom::Log() {
    return *CachedLogger();
}
