#pragma once
#include <sys/wait.h>
#include <unistd.h>

#include <csignal>
#include <cstdio>
#include <cstdlib>

#include "test_case.h"

namespace check {
template <typename Operation>
void assertion_fails(const Operation& operation) {
#ifndef NDEBUG
    const pid_t child = fork();
    require(child >= 0, "fork contract test");
    if (child == 0) {
        // Avoid invoking the host core-dump collector for expected failures.
        std::signal(SIGABRT, [](int) { std::_Exit(128 + SIGABRT); });
        if (std::freopen("/dev/null", "w", stderr) == nullptr)
            std::_Exit(1);
        operation();
        std::_Exit(0);
    }
    int status = 0;
    require(waitpid(child, &status, 0) == child, "wait for contract test");
    require(WIFEXITED(status) && WEXITSTATUS(status) == 128 + SIGABRT, "invalid input triggers debug assertion");
#else
    (void) operation;
    require(true, "invalid preconditions are checked only in Debug");
#endif
}
}   // namespace check
