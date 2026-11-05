#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "runtime/logger.h"
#include "runtime/profiler.h"
#include "simulation.h"
#include "test_case.h"

namespace {
#ifdef USE_TIMERS
void disabled_and_nested_timers() {
    timer::clearFullProfile();
    timer::profiling_enabled = false;
    timer::clearTree();
    std::ostringstream restarted_tree;
    timer::printTreeProfile(restarted_tree);
    check::equal(timer::globalTimer.fieldCalls("all"), 1, "cleared root restarts while collection is disabled");
    check::require(std::isfinite(timer::globalTimer.fieldTime("all")) && timer::globalTimer.fieldTime("all") >= 0,
                   "restarted root records finite nonnegative duration");
    {
        timer::commonTimer disabled("disabled");
        timer::flatTimer not_started(timer::NoStart{});
#pragma omp parallel num_threads(2)
        { RECORD_TIMER_NAMED("disabled worker"); }
        // A timer disabled at construction stays inactive when collection resumes.
        timer::profiling_enabled = true;
    }
    timer::profiling_enabled = false;
    for (int thread = 0; thread < timer::maxThreads; ++thread)
        check::equal(timer::currEvents[thread].val, int64_t{0}, "disabled timers must not reserve events");
    check::equal(timer::globalTimer.fieldCalls("disabled"), 0, "disabled timer must not enter the tree");
    timer::profiling_enabled = true;
    timer::clearTree();
    {
        timer::commonTimer outer("outer");
        { timer::commonTimer inner("inner"); }
        timer::profiling_enabled = false;
        outer.finish();   // A timer already started must finish even if collection is disabled.
    }
    int named_events = 0;
    for (int event = 0; event < timer::currEvents[0].val; ++event) {
        const auto& recorded = timer::events[event];
        if (std::string(recorded.name) == "outer" || std::string(recorded.name) == "inner") {
            ++named_events;
            check::require(recorded.end >= recorded.start, "completed nested event timestamps");
        }
    }
    check::equal(named_events, 2, "nested timers record each named scope exactly once");
    check::equal(timer::globalTimer.fieldCalls("outer"), 1, "outer timer finishes once");
    check::equal(timer::globalTimer.fieldCalls("inner"), 1, "inner timer finishes once");
    timer::clearFullProfile();
    timer::clearTree();
    check::near(timer::globalTimer.fieldTime("all"), 0, 0, "clear resets completed root duration");
    check::equal(timer::globalTimer.fieldCalls("all"), 1, "clear restarts root with one call");
    check::equal(timer::globalTimer.fieldCalls("outer"), 0, "clear removes completed child timers");
}
#endif

class TimedSteps : public Simulation {
   public:
    explicit TimedSteps(int first_profiled_step) : first_profiled_step(first_profiled_step) {
        system_config = {{"Scheme", "timing-test"},
                         {"StartTimeStep", 0},
                         {"LastTimestep", 4},
                         {"NumCellsX", 1},
                         {"NumCellsY", 1},
                         {"NumCellsZ", 1},
                         {"Dx", 1}};
    }
    void prepare_step(int step) override {
        const bool expected = step >= first_profiled_step && step <= 3;
        check::equal(timer::profiling_enabled, expected, "collection follows inclusive profiling window");
        RECORD_TIMER_NAMED("prepare");
    }
    void make_step(int) override {
#pragma omp parallel num_threads(2)
        { RECORD_TIMER_NAMED("worker"); }
    }
    void make_diagnostic(int) override {
        RECORD_TIMER_NAMED("diagnostic");
    }

   private:
    int first_profiled_step;
};

void inclusive_window_and_step_output(bool include_zero) {
    // A bounded synthetic loop checks orchestration without a particle calculation.
    setenv("TIMER_START_STEP", include_zero ? "0" : "2", 1);
    setenv("TIMER_END_STEP", "3", 1);
    setenv("TIMERED_SIM_STEPS", "2", 1);
    for (int step : {1, 2, 3, 4}) std::filesystem::remove("profile_it_" + std::to_string(step) + ".json");
    TimedSteps simulation(include_zero ? 0 : 2);
    timer::profiling_enabled = true;
    { timer::flatTimer startup("startup before calculate"); }
    logger::init();
    simulation.calculate();
    simulation.finalize();
    for (int step : {1, 2, 3, 4}) {
#ifdef USE_TIMERS
        const bool expected = step == 2 || step == 3;
#else
        const bool expected = false;
#endif
        const std::string filename = "profile_it_" + std::to_string(step) + ".json";
        check::equal(std::filesystem::exists(filename), expected, "profile file only at period/end of window");
        if (expected) {
            std::ifstream file(filename);
            nlohmann::json events;
            file >> events;
            bool found_worker = false;
            bool found_startup = false;
            for (const auto& event : events) {
                found_worker = found_worker || event.value("name", "") == "worker";
                found_startup = found_startup || event.value("name", "") == "startup before calculate";
            }
            check::require(found_worker, "window profile contains parallel worker events");
            if (step == 2)
                check::equal(found_startup, include_zero, "step-zero window controls initialization events");
        }
    }
    std::ifstream timing("beren3d.log");
    std::string line;
    int expected_step = 1;
    while (std::getline(timing, line)) {
        const auto marker = line.find("Step wall time: ");
        if (marker == std::string::npos)
            continue;
        std::istringstream record(line);
        char open_bracket = 0, close_bracket = 0;
        int step = -1;
        std::string step_word, wall_word, time_word, unit;
        double seconds = -1;
        check::require(static_cast<bool>(record >> open_bracket >> step >> close_bracket >> step_word >> wall_word >>
                                         time_word >> seconds >> unit),
                       "wall time log record contains step and seconds");
        check::require(open_bracket == '[' && close_bracket == ']' && step_word == "Step" && wall_word == "wall" &&
                           time_word == "time:" && unit == "s",
                       "production timing log format");
        check::equal(step, expected_step++, "production timing step");
        check::require(std::isfinite(seconds) && seconds >= 0, "finite nonnegative production timing");
    }
    check::equal(expected_step, 5, "wall time exists inside and outside profiling window");
}
}   // namespace

int main(int argc, char** argv) {
    const bool include_zero = argc == 2 && std::string(argv[1]) == "--include-zero";
    return check::run([&] {
        check::require(argc == 1 || include_zero, "expected no arguments or --include-zero");
#ifdef USE_TIMERS
        if (!include_zero)
            check::scenario("disabled timers and nested timer lifetime", disabled_and_nested_timers);
#endif
        // Each window runs in its own process because envOptions caches the environment.
        check::scenario(include_zero ? "step-zero window preserves initialization events"
                                     : "inclusive profiling window preserves production timing",
                        [&] { inclusive_window_and_step_output(include_zero); });
    });
}
