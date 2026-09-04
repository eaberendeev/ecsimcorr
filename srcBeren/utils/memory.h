#include <memory>

#include "timer.h"

template <typename T>
struct SmartPtr : public timer::flatTimer, public std::unique_ptr<T[], decltype(&std::free)> {
    static constexpr int64_t defaultMemoryAlign = 128;

    using basePtr = std::unique_ptr<T[], decltype(&std::free)>;

    SmartPtr() : timer::flatTimer(timer::NoStart{}), basePtr(nullptr, std::free), size(0) {
    }

    SmartPtr(SmartPtr&& other) : timer::flatTimer(timer::NoStart{}), basePtr(std::move(other)), size(other.size) {
        other.size = 0;
    }

    SmartPtr(int64_t sizeIn)
        : timer::flatTimer(std::source_location::current().function_name()),
          basePtr(static_cast<T*>(std::aligned_alloc(defaultMemoryAlign, sizeof(T) * sizeIn)), std::free),
          size(sizeIn) {
        timer::flatTimer::finish();
    }

    T& operator()(int64_t ix) {
        assert(ix >= 0 && ix < size);
        return get()[ix];
    }

    T& operator[](int64_t ix) {
        assert(ix >= 0 && ix < size);
        return get()[ix];
    }

    SmartPtr& operator=(SmartPtr&& other) {
        *(static_cast<basePtr*>(this)) = std::move(other);
        size = other.size;
        other.size = 0;
        return *this;
    }

    T* get() {
        return basePtr::get();
    }

    ~SmartPtr() {
        timer::flatTimer::start(std::source_location::current().function_name(), size * sizeof(T),
                                timer::MeasureUnit::byte_no_bandwidth);
    }

    int64_t size;
};
