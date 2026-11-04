#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "runtime/profiler.h"
#include "simulation.h"
#include "test_case.h"

namespace {
#ifdef USE_TIMERS
void disabled_and_nested_timers() {
    timer::clearFullProfile();
    timer::setProfilingEnabled(false);
    timer::clearTree();
    std::ostringstream disabled_tree;
    timer::printTreeProfile(disabled_tree);
    check::near(timer::globalTimer.fieldTime("all"), 0, 0, "disabled root print must keep zero duration");
    {
        timer::commonTimer disabled("disabled");
        timer::flatTimer not_started(timer::NoStart{});
#pragma omp parallel num_threads(2)
        { RECORD_TIMER_NAMED("disabled worker"); }
        // A timer disabled at construction stays inactive when collection resumes.
        timer::setProfilingEnabled(true);
    }
    timer::setProfilingEnabled(false);
    for (int thread = 0; thread < timer::maxThreads; ++thread)
        check::equal(timer::currEvents[thread].val, int64_t{0}, "disabled timers must not reserve events");
    check::equal(timer::globalTimer.fieldCalls("disabled"), 0, "disabled timer must not enter the tree");
    timer::setProfilingEnabled(true);
    timer::clearTree();
    {
        timer::commonTimer outer("outer");
        { timer::commonTimer inner("inner"); }
        timer::setProfilingEnabled(false);
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
}
#endif

class TimedSteps : public Simulation {
   public:
    TimedSteps() {
        system_config = {{"Scheme", "timing-test"},
                         {"StartTimeStep", 0},
                         {"LastTimestep", 4},
                         {"NumCellsX", 1},
                         {"NumCellsY", 1},
                         {"NumCellsZ", 1},
                         {"Dx", 1}};
    }
    void prepare_step(int step) override {
#ifdef USE_TIMERS
        const bool expected = step >= 2 && step <= 3;
#else
        const bool expected = false;
#endif
        check::equal(timer::profilingEnabled(), expected, "collection enabled only on inclusive steps 2..3");
        RECORD_TIMER_NAMED("prepare");
    }
    void make_step(int) override {
#pragma omp parallel num_threads(2)
        { RECORD_TIMER_NAMED("worker"); }
    }
    void make_diagnostic(int) override {
        RECORD_TIMER_NAMED("diagnostic");
    }
};

void inclusive_window_and_step_output() {
    // A bounded synthetic loop checks orchestration without a particle calculation.
    setenv("TIMER_START_STEP", "2", 1);
    setenv("TIMER_END_STEP", "3", 1);
    setenv("TIMERED_SIM_STEPS", "2", 1);
    for (int step : {1, 2, 3, 4}) std::filesystem::remove("profile_it_" + std::to_string(step) + ".json");
    TimedSteps simulation;
    simulation.calculate();
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
            for (const auto& event : events) found_worker = found_worker || event.value("name", "") == "worker";
            check::require(found_worker, "window profile contains parallel worker events");
        }
    }
    std::ifstream timing("Performance/step.time");
    std::string header;
    std::getline(timing, header);
    check::require(header == "TimeStep WallSeconds", "production timing remains enabled");
    for (int expected = 1; expected <= 4; ++expected) {
        int step = -1;
        double seconds = -1;
        check::require(static_cast<bool>(timing >> step >> seconds), "wall time exists inside and outside window");
        check::equal(step, expected, "production timing step");
        check::require(std::isfinite(seconds) && seconds >= 0, "finite nonnegative production timing");
    }
}
}   // namespace

int main() {
    return check::run([] {
#ifdef USE_TIMERS
        check::scenario("disabled timers and nested timer lifetime", disabled_and_nested_timers);
#endif
        check::scenario("inclusive profiling window preserves production timing", inclusive_window_and_step_output);
    });
}
