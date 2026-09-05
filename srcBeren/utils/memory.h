#include <memory>

#include "timer.h"

template <typename T>
struct SmartPtr : public timer::flatTimer, public std::unique_ptr<T[], decltype(&std::free)> {
    static constexpr int64_t memoryAlign = 128;

    static int64_t adjustAllocSize(int64_t minSize) {
        return (minSize + memoryAlign - 1) / memoryAlign * memoryAlign;
    }

    using basePtr = std::unique_ptr<T[], decltype(&std::free)>;

    SmartPtr() : timer::flatTimer(timer::NoStart{}), basePtr(nullptr, std::free), size(0) {
    }

    SmartPtr(SmartPtr&& other) : timer::flatTimer(timer::NoStart{}), basePtr(std::move(other)), size(other.size) {
        other.size = 0;
    }

    SmartPtr(int64_t sizeIn)
        : timer::flatTimer(std::source_location::current().function_name()),
          basePtr(static_cast<T*>(std::aligned_alloc(memoryAlign, adjustAllocSize(sizeof(T) * sizeIn))), std::free),
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

template <typename T>
struct SimpleArrayBuffer : public std::vector<T> {
    using std::vector<T>::size;
    using std::vector<T>::capacity;
    using std::vector<T>::reserve;

    T& operator()(const int64_t i) {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    const T& operator()(const int64_t i) const {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    T& operator[](const int64_t i) {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    const T& operator[](const int64_t i) const {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    // resize buffer to fit new size, does not carry about old storage
    void resizeAndReset(int64_t newSize) {
        RECORD_TIMER_PARAMS(newSize);
        if (newSize <= static_cast<int64_t>(capacity())) {
            timer::commonTimer timerResize("this resize");
            resize(newSize);
            return;
        }
        std::vector<T> other;
        other.reserve(newSize * 2);
        timer::commonTimer timerOtherResize("other resize");
        other.resize(newSize);
        timerOtherResize.finish();
        static_cast<std::vector<T>&>(*this) = std::move(other);
    }

    /// TODO: implement push_back with parallel reallocation
    /// TODO: implement emplace_back with parallel reallocation

   private:
    using std::vector<T>::resize;
};
