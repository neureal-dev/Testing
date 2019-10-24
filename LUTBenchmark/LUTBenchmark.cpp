#include "LUTBenchmark.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <random>
#include <utility>
#include <variant>
#include <vector>

#include <benchmark/benchmark.h>

struct A {

    A(uint32_t id) : id_(id), id1_{}, id2_{}, id4_{} {}

    uint32_t GetId() const { return id_; }

    uint32_t id_;
    uint64_t id1_;
    uint64_t id2_;
    std::array<uint64_t, 32> id4_;
};

namespace {

struct Comp {

    inline bool operator()(const A* s, uint32_t i) const noexcept { return s->GetId() < i; }

    inline bool operator()(uint32_t i, const A* s) const noexcept { return i < s->GetId(); }

    inline bool operator()(const A& s, uint32_t i) const noexcept { return s.GetId() < i; }

    inline bool operator()(uint32_t i, const A& s) const noexcept { return i < s.GetId(); }
};

template <class T>
constexpr void ignore(const T&) {}

} // namespace

struct VectorSearchFixture : benchmark::Fixture {

    void SetUp(const ::benchmark::State& state)
    {
        if (state.thread_index == 0) {

            searches_.resize(state.range(0), 0);
            records_.reserve(state.range(0));

            benchmark::DoNotOptimize(records_.data());

            benchmark::DoNotOptimize(searches_.data());

            std::mt19937 rng;

            rng.seed(std::random_device()());

            std::variant<std::uniform_real_distribution<double>,
                    std::normal_distribution<double>,
                    std::exponential_distribution<double>>
                    distribution;

            if (state.range(1) == 0) {
                distribution = std::uniform_real_distribution<double>(1.5, 7.5);
            } else if (state.range(1) == 1) {
                distribution = std::normal_distribution<double>(5.0, 3.0);
            } else {
                distribution = std::exponential_distribution<double>(0.5);
            }

            for (int64_t i = 0; i < state.range(0); ++i) {
                double number = std::visit([&](auto& dist) -> double { return dist(rng); }, distribution);
                if (number > 0.0 && number < 10.0) {
                    searches_[i] = static_cast<uint32_t>((number / 10.0) * state.range(0));
                }
            }

            auto tmp = searches_;

            std::partial_sort_copy(std::begin(searches_), std::end(searches_), std::begin(tmp), std::end(tmp));

            tmp.erase(std::unique(std::begin(tmp), std::end(tmp)), std::end(tmp));

            for (auto s : tmp) {
                // records_.push_back(new A(s));
                records_.emplace_back(s);
            }

            benchmark::ClobberMemory();
        }
    }

    void TearDown(const ::benchmark::State& state)
    {
        if (state.thread_index == 0) {
            // for (auto item : records_) {
            //    delete item;
            //}
            records_.clear();
            searches_.clear();
            records_.shrink_to_fit();
            searches_.shrink_to_fit();
            benchmark::ClobberMemory();
        }
    }

    // std::vector<A*> records_;
    std::vector<A> records_;
    std::vector<uint32_t> searches_;
};

template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator HybriddddInterpolationSearch(
        RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
    using difference_type = typename std::iterator_traits<RandomIterator>::difference_type;

    difference_type count = std::distance(begin, end);

#ifdef NDEBUG
    RandomIterator last = std::prev(end);
#else
    RandomIterator last = count ? std::prev(end) : end;
#endif
    while (count > 0) {
        if (!comp(*begin, key)) {
            end = !comp(key, *begin) ? begin : end;
            break;
        }

        if (!comp(key, *last)) {
            end = !comp(*last, key) ? last : end;
            break;
        }

        difference_type probe = static_cast<difference_type>((count - 1) * lerp(*begin, *last, key));

        if (comp(key, begin[probe])) {
            // probe -= count;
            // std::advance(last, probe - count);
            // count = probe;

            /*/
            probe = count >> 1;
            if (!comp(key, begin[probe])) {
                    count -= probe;
                    std::advance(begin, probe);
            } else {
                    probe -= count;
                    std::advance(last, probe);
                    count += probe;
            }
            /*/
            probe = count >> 1;
            if (!comp(key, begin[probe])) {
                std::advance(begin, probe);
                count -= probe;
            } else {
                std::advance(last, probe - count);
                count = probe;
            }
            //*/
        } else if (comp(begin[probe], key)) {
            // count -= ++probe;
            // std::advance(begin, probe);

            /*/
            probe = count >> 1;
            if (!comp(key, begin[probe])) {
            count -= probe;
            std::advance(begin, probe);
            } else {
            probe -= count;
            std::advance(last, probe);
            count += probe;
            }
            /*/
            probe = count >> 1;
            if (!comp(key, begin[probe])) {
                std::advance(begin, probe);
                count -= probe;
            } else {
                std::advance(last, probe - count);
                count = probe;
            }
            //*/
        } else {
            end = begin;
            std::advance(end, probe);
            break;
        }
    }
    return end;
}
/**
 * Implementation of Hybrid Interpolation search to optimize access path for
 * number ordered arrarys
 */
template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator HybridInterpolationSearch(
        RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
    using difference_type = typename std::iterator_traits<RandomIterator>::difference_type;

    difference_type count = std::distance(begin, end);

#ifdef NDEBUG
    RandomIterator last = std::prev(end);
#else
    RandomIterator last = count ? std::prev(end) : end;
#endif

    while (count > 0) {
        if (!comp(*begin, key)) {
            end = !comp(key, *begin) ? begin : end;
            break;
        }

        if (!comp(key, *last)) {
            end = !comp(*last, key) ? last : end;
            break;
        }

        auto probe = static_cast<difference_type>((count - 1) * lerp(*begin, *last, key));

        if (comp(key, begin[probe])) {
            auto mid = probe >> 1;
            if (!comp(key, begin[mid])) {
                std::advance(last, probe - count);
                count = probe - mid;
                std::advance(begin, mid);
            } else {
                std::advance(last, mid - count);
                count = mid;
            }
        } else if (comp(begin[probe], key)) {
            auto mid = (probe + count) >> 1;
            if (!comp(key, begin[mid])) {
                std::advance(begin, mid);
                count -= mid;
            } else {
                std::advance(last, mid - count);
                count = mid - probe;
                std::advance(begin, probe);
            }
        } else {
            end = begin;
            std::advance(end, probe);
            break;
        }
    }
    return end;
}

/**
 * Implementation of Interpolation search to optimize access path for number
 * ordered arrarys
 */
template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator InterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
    using difference_type = typename std::iterator_traits<RandomIterator>::difference_type;

    difference_type count = std::distance(begin, end);

#ifdef NDEBUG
    RandomIterator last = std::prev(end);
#else
    RandomIterator last = count ? std::prev(end) : end;
#endif

    while (count > 0) {
        if (!comp(*begin, key)) {
            end = !comp(key, *begin) ? begin : end;
            break;
        }

        if (!comp(key, *last)) {
            end = !comp(*last, key) ? last : end;
            break;
        }

        difference_type probe = static_cast<difference_type>((count - 1) * lerp(*begin, *last, key));

        if (comp(key, begin[probe])) {
            std::advance(last, probe - count);
            count = probe;
        } else if (comp(begin[probe], key)) {
            std::advance(begin, ++probe);
            count -= probe;
        } else {
            end = std::next(begin, probe);
            break;
        }
    }
    return end;
}

template <typename Integer>
class FibonacciIterator {
public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = Integer;
    using difference_type = Integer;
    using pointer = Integer*;
    using reference = Integer&;

    explicit FibonacciIterator() : curr(0), next(1) {}

    explicit FibonacciIterator(Integer target) : curr(0), next(1)
    {
        while (next <= target) {
            next = curr + next;
            curr = next - curr;
        }
    }

    const Integer& operator*() const { return curr; }

    const Integer* operator->() const { return &**this; }

    FibonacciIterator& operator++()
    {
        next = curr + next;
        curr = next - curr;
        return *this;
    }

    const FibonacciIterator operator++(int)
    {
        FibonacciIterator result = *this;
        ++*this;
        return result;
    }

    FibonacciIterator& operator--()
    {
        curr = next - curr;
        next = next - curr;
        return *this;
    }

    const FibonacciIterator operator--(int)
    {
        FibonacciIterator result = *this;
        --*this;
        return result;
    }

private:
    Integer curr;
    Integer next;
};

/* Comparison functions for FibonacciIterator. */
template <typename Integer>
constexpr bool operator==(const FibonacciIterator<Integer>& lhs, const FibonacciIterator<Integer>& rhs)
{
    return lhs.curr == rhs.curr;
}

/* Disequality implemented in terms of equality. */
template <typename Integer>
bool operator!=(const FibonacciIterator<Integer>& lhs, const FibonacciIterator<Integer>& rhs)
{
    return !(lhs == rhs);
}

/**
 * Implementation of Fibonacci search uses a Fibonacci iterator to access
 * the consecutive Fibonacci numbers
 */
template <typename RandomIterator, typename Value, typename Comparator>
RandomIterator FibonacciSearch(RandomIterator begin, RandomIterator end, const Value& key, Comparator comp)
{
    using difference_type = typename std::iterator_traits<RandomIterator>::difference_type;

    FibonacciIterator<difference_type> itr(std::distance(begin, end));

    while (*itr > 0) {
        if (comp(key, begin[*itr - 1])) {
            --itr;
        } else if (comp(begin[*itr - 1], key)) {
            std::advance(begin, *itr);
            do {
                --itr;
            } while (*itr > std::distance(begin, end));
        } else {
            std::advance(begin, *itr - 1);
            return begin;
        }
    }
    return (begin != end && !(comp(*begin, key)) && !(comp(key, *begin))) ? begin : end;
}

template <class ForwardIt, typename Value, typename Comparator>
ForwardIt BranchLessBinarySearch(ForwardIt begin, ForwardIt end, const Value& key, Comparator comp)
{
    using difference_type = typename std::iterator_traits<ForwardIt>::difference_type;

    difference_type count = std::distance(begin, end);

    while (difference_type half = count >> 1) {
        if (!comp(key, begin[half])) {
            std::advance(begin, half);
        }
        count -= half;
    }
    return (begin != end && !comp(*begin, key) && !comp(key, *begin)) ? begin : end;
}

template <class ForwardIt, class Value, typename Comparator>
ForwardIt BranchFullBinarySearch(ForwardIt begin, ForwardIt end, const Value& key, Comparator comp)
{
    using difference_type = typename std::iterator_traits<ForwardIt>::difference_type;

    difference_type count = std::distance(begin, end);

    while (difference_type half = count >> 1) {
        if (comp(begin[half], key)) {
            std::advance(begin, ++half);
        } else {
            if (!comp(key, begin[half])) {
                std::advance(begin, half);
                break;
            }
        }
        count -= half;
    }
    return (begin != end && !comp(*begin, key) && !comp(key, *begin)) ? begin : end;
}

BENCHMARK_DEFINE_F(VectorSearchFixture, HybridInterpolationIt)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = InterpolationSearch(std::cbegin(records_),
                std::cend(records_),
                rid,
                Comp{},
                //            [](const A* first, const A* last, uint32_t key) -> double
                //            { return (double(key) - first->GetId()) / (last->GetId() -
                //            first->GetId()); });
                [](const A& first, const A& last, uint32_t key) -> double {
                    return (double(key) - first.id_) / (last.id_ - first.id_);
                });

        if (low != std::cend(records_)) {
            if ((*low).id_ != rid) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*low).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, HybridInterpolationIt)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_DEFINE_F(VectorSearchFixture, InterpolationIt)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = InterpolationSearch(std::cbegin(records_),
                std::cend(records_),
                rid,
                Comp{},
                //            [](const A* first, const A* last, uint32_t key) -> double
                //            { return (double(key) - first->GetId()) / (last->GetId() -
                //            first->GetId()); });
                [](const A& first, const A& last, uint32_t key) -> double {
                    return (double(key) - first.id_) / (last.id_ - first.id_);
                });

        if (low != std::cend(records_)) {
            if ((*low).id_ != rid) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*low).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, InterpolationIt)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_DEFINE_F(VectorSearchFixture, FibonacciIt)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = FibonacciSearch(std::cbegin(records_), std::cend(records_), rid, Comp{});

        if (low != std::cend(records_)) {
            if ((*low).id_ != rid) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*low).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, FibonacciIt)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_DEFINE_F(VectorSearchFixture, BranchLessIt)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = BranchLessBinarySearch(std::cbegin(records_), std::cend(records_), rid, Comp{});

        if (low != std::cend(records_)) {
            if ((*low).id_ != rid) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*low).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, BranchLessIt)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_DEFINE_F(VectorSearchFixture, BranchFullIt)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = BranchFullBinarySearch(std::cbegin(records_), std::cend(records_), rid, Comp{});

        if (low != std::cend(records_)) {
            if ((*low).id_ != rid) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*low).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, BranchFullIt)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_DEFINE_F(VectorSearchFixture, std_equal_range)
(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto _ : state) {
        ignore(_);
        auto id = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto first = std::equal_range(std::cbegin(records_), std::cend(records_), id, Comp{});
        if (first.first != first.second) {
            if ((*first.first).id_ != id) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*first.first).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_equal_range)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_DEFINE_F(VectorSearchFixture, std_lower_bound)
(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto _ : state) {
        ignore(_);
        auto id = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto first =
                std::lower_bound(std::cbegin(records_), std::cend(records_), id, [](const auto& e, auto id) -> bool {
                    return e.id_ < id;
                });
        if (first != std::cend(records_)) {
            if ((*first).id_ != id) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*first).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_lower_bound)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_DEFINE_F(VectorSearchFixture, std_partition_point)
(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto _ : state) {
        ignore(_);
        auto id = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto first = std::partition_point(
                std::cbegin(records_), std::cend(records_), [=](const auto& e) -> bool { return e.id_ < id; });
        if (first != std::cend(records_)) {
            if ((*first).id_ != id) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*first).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_partition_point)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->Threads(6);

BENCHMARK_MAIN();
