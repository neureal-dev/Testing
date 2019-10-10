// Test.cpp : Defines the entry point for the application.
//

#include "LUTBenchmark.h"
#include <algorithm>
#include <array>
#include <benchmark/benchmark.h>
#include <chrono>
#include <mutex>
#include <random>
#include <shared_mutex>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

class A {
public:
    A()
        : id_(0)
    {
    }
    A(uint32_t id)
        : id_(id)
    {
    }
    uint32_t GetId() const { return id_; }

    uint32_t id_;
    uint64_t id1_;
    uint64_t id2_;
    //uint64_t id3_;
    std::array<uint64_t, 20> id4_;
};

class VectorSearchFixture : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State& state)
    {
        //std::cout << searches_.size() << std::endl;
        //records_ = std::vector<A*>(state.range_x());
        //searches_ = std::vector<uint32_t>(state.range_x(), 0);
        //benchmark::DoNotOptimize(records_.data());

        //for (int64_t i = state.range_x() - 1; i >= 0; --i) {
        //    records_[i] = new A(i + i + 65535);
        //}

        if (searches_.size() < state.range_x()) {
            uint32_t strt = std::max<uint32_t>(searches_.size(), 1);
            searches_.resize(state.range_x(), 0);
            benchmark::DoNotOptimize(searches_.data());
            std::mt19937 rng;
            rng.seed(std::random_device()());
            
			/*
			//std::uniform_real_distribution<double> distribution(0.0, 10.0);
			//std::normal_distribution<double> distribution(5.0, 3.0);
            //std::exponential_distribution<double> distribution(3.5);

            for (int i = 0; i < state.range_x(); ++i) {
                double number = distribution(rng);
				if (number > 0.0 && number < 10.0) {
					number /= 10.0;
					number *= state.range_x();
					searches_[i] = static_cast<uint64_t>(number);
				}
            }
			*/
			std::uniform_int_distribution<uint64_t> distribution(0, state.range_x());
			for (int i = 0; i < state.range_x(); ++i) {
				searches_[i] = static_cast<uint64_t>(distribution(rng));
			}

			auto tmp = searches_;
			std::partial_sort_copy(std::begin(searches_), std::end(searches_), std::begin(tmp), std::end(tmp));
			tmp.erase(std::unique(std::begin(tmp), std::end(tmp)), std::end(tmp));
			std::sort(std::begin(tmp), std::end(tmp));
			for (auto s : tmp) {
				//std::cout << state.range_x() << " "<< s << std::endl;
				records_.push_back(new A(s));
			}
            //sorting searches increasing cache hits
			//std::sort(std::begin(searches_), std::end(searches_));

            //for (int64_t i = searches_.size(); i > strt; ) {
//            for (int64_t i = searches_.size(), j = strt; i > j && i > strt; --i) {
//                searches_[i - 1] = uniform_dist(rng) % (j * 16);
            //}
            //}
        }
		benchmark::DoNotOptimize(records_.data());
		benchmark::DoNotOptimize(searches_.data());
        //		for (auto e : searches_) {
        //		std::cout << e << std::endl;
        //}
    }

    void TearDown(const ::benchmark::State& state)
    {
        for (auto item : records_) {
            delete item;
        }
        records_.clear();
		searches_.clear();
    }

    //    private:
    //std::vector<std::vector<A*>, >
    std::vector<A*> records_;
    std::vector<uint32_t> searches_;
};

enum locate_t { EQUAL,
    LEFT,
    RIGHT };

int nary_search(const std::vector<A*>& a, int const key, int const N)
{
    std::vector<int> mid(N + 1);
    std::vector<locate_t> locate(N + 2);

    locate[0] = RIGHT;
    locate[N + 1] = LEFT;

    int lo = 0;
    int hi = a.size() - 1;
    int pos = -1;

    while (lo <= hi && pos == -1) {
        mid[0] = lo - 1;

        double const step = (hi - lo + 1) / (N + 1);

        for (int i = 1; i <= N; i++) {
            int const offset = step * i + (i - 1);
            int const lmid = mid[i] = lo + static_cast<int>(offset);

            if (lmid <= hi) {
                if (a[lmid]->id_ > key) {
                    locate[i] = LEFT;
                } else if (a[lmid]->id_ < key) {
                    locate[i] = RIGHT;
                } else {
                    locate[i] = EQUAL;
                    pos = lmid;
                }
            } else {
                mid[i] = hi + 1;
                locate[i] = LEFT;
            }
        }
        for (int i = 1; i <= N; i++) {
            if (locate[i] != locate[i - 1]) {
                lo = mid[i - 1] + 1;
                hi = mid[i] - 1;
            }
        }
        if (locate[N] != locate[N + 1]) {
            lo = mid[N] + 1;
        }
    }
    return pos;
}

int binary_cmov(const std::vector<A*>& records, uint32_t rid)
{
    size_t min = 0, max = records.size();
    while (min < max) {
        size_t middle = (min + max) >> 1;
        //middle /= 2;
        auto val = records[middle]->id_;
        min = rid > val ? middle + 1 : min;
        max = rid > val ? max : middle;
    }
    return min;
}

int binary_cmov2(const std::vector<A*>& records, uint32_t const rid)
{
    size_t min = 0, max = records.size();
    while (min < max) {
        size_t middle = (min + max) >> 1;
        size_t middle1 = middle + 1;
        min = rid > records[middle]->id_ ? middle1 : min;
        max = rid <= records[middle]->id_ ? middle : max;
    }
    return min;
}

size_t b2search(const std::vector<A*>& records, uint32_t const rid)
{
    size_t size = records.size();
    size_t low = 0;
    while (size_t half = size >> 1) {
        size_t middle = low + half;
        //auto val = records[middle]->id_;
        low = (records[middle]->GetId() <= rid) ? middle : low;
        size = size - half;
    }
    return (size && records[low]->id_ == rid) ? low : -1;
}

int lsearch(const std::vector<A*>& records, uint32_t key)
{
    for (size_t i = 0; i < records.size(); i++) {
        if (records[i]->id_ == key) {
            return i;
        }
    }
    return -1;
}
struct Comp {
    inline bool operator()(const A* s, uint32_t i) const noexcept { return s->GetId() < i; }
    inline bool operator()(uint32_t i, const A* s) const noexcept { return i < s->GetId(); }
};

namespace qb {

template <class ForwardIt, class T>
inline ForwardIt branch_less_binary_find_n(ForwardIt first, size_t size, const T& key)
{
    using difference_type = std::iterator_traits<ForwardIt>::difference_type;

    while (difference_type half = size >> 1) {
        ForwardIt middle = first;
        std::advance(middle, half);
        if (!((*middle)->GetId() > key)) {
            std::advance(first, half);
        }
        size = size - half;
    }
    return first;
}


/**
* Implementation of Hybrid Interpolation search to optimize access path for number ordered arrarys
*/
template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator HybridInterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

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
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, probe - count);
				count = probe;
			}
		} else if (comp(begin[probe], key)) {
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, probe - count);
				count = probe;
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
* Implementation of Interpolation search to optimize access path for number ordered arrarys
*/
template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator InterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

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

/**
* Implementation of Fibonacci search uses a Fibonacci iterator to access
* the consecutive Fibonacci numbers
*/
template <typename RandomIterator, typename Value, typename Comparator>
RandomIterator FibonacciSearch(RandomIterator begin, RandomIterator end,
	const Value& key, Comparator comp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

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


/**
* Implementation of Binary search uses branch less logic to access
* small arrays
*/
template <class ForwardIt, typename Value, typename Comparator>
inline ForwardIt BranchLessBinarySearch(ForwardIt begin, ForwardIt end, const Value& key, Comparator comp)
{
	using difference_type = std::iterator_traits<ForwardIt>::difference_type;
	difference_type size = std::distance(begin, end);

	while (difference_type half = size >> 1) {
		if (!comp(key, begin[half])) {
			std::advance(begin, half);
		}
		size = size - half;
	}
	return (begin != end && !comp(*begin, key) && !comp(key, *begin)) ? begin : end;
}

template <class ForwardIt, typename Value, typename Comparator>
inline ForwardIt BranchLessBinarySearch(ForwardIt begin, ForwardIt end, const Value& key, Comparator comp)
{
	using difference_type = std::iterator_traits<ForwardIt>::difference_type;
	difference_type size = std::distance(begin, end);

	while (difference_type half = size >> 1) {
		if (!comp(key, begin[half])) {
			std::advance(begin, half);
		}
		size = size - half;
	}
	return (begin != end && !comp(*begin, key) && !comp(key, *begin)) ? begin : end;
}


template <class ForwardIt, class T>
inline ForwardIt branch_full_binary_find_n(ForwardIt first, ForwardIt last, size_t size, const T& key)
{
    using difference_type = std::iterator_traits<ForwardIt>::difference_type;
    difference_type half = std::distance(first, last) >> 1;
    ForwardIt middle = std::next(first, half);

    while (first < last) {
        auto val = (*middle)->GetId();

        if (key > val) {
            first = middle;
            std::advance(first, 1);
        } else {
            last = middle;
            if (key == val) {
                //ind = middle;
                break;
            }
        }
        half = std::distance(first, last) >> 1;
        middle = first;
        std::advance(middle, half);
    }
    return last;
}

template <class ForwardIt, class T>
inline ForwardIt equal_range(ForwardIt first, ForwardIt last, const T& key)
{
    using difference_type = std::iterator_traits<ForwardIt>::difference_type;
    difference_type size = std::distance(first, last);
    return (size == 0) ? last : size < 0x8FFFF ? branch_less_binary_find_n(first, size, key) : branch_full_binary_find_n(first, last, size, key);
}

};

size_t combine_search(std::vector<A*> const& data, uint32_t const key)
{
	auto ind = qb::BranchLessBinarySearch(std::cbegin(data), std::cend(data), key, Comp{});
	return ind != std::cend(data)? std::distance(std::cbegin(data), ind) : -1;
		//auto ind = qb::equal_range(std::cbegin(data), std::cend(data), key);
    //auto ind = std::lower_bound(std::cbegin(data), std::cend(data), key,
    //	[](const auto& lhs, const auto& i) -> bool { return lhs->GetId() < i; });
    //: btsearch(std::cbegin(data), std::cend(data), rid);
		//return (ind != std::cend(data) && !((*ind)->GetId() < key) && !((*ind)->GetId() > key)) ? std::distance(std::cbegin(data), ind) : -1;
    //auto ind = std::partition_point(std::cbegin(data), std::cend(data),
    //		[=](const auto& lhs) -> bool { return lhs->GetId() < rid; });
    //return (ind != std::cend(data) && (*ind)->GetId() == rid )? std::distance(std::cbegin(data), ind) : -1;
}

size_t bsearch(const std::vector<A*>& data, uint32_t const key)
{
    size_t low = 0;

    size_t high = data.size();
    size_t ind = -1;
    size_t mid = high >> 1;

    while (low < high) {
        const auto val = data[mid]->GetId();
        if (key > val) {
            low = mid + 1;
        } else {
            if (key == val) {
                ind = mid;
                break;
            }
            high = mid;
        }
        mid = (low + high) >> 1;
    }
    return ind;
}

// size_t bsearch(const std::vector<A*>& data, uint32_t const key) { return bbsearch(data, key, 0, data.size()); }

size_t fallback_bsearch(const std::vector<A*>& data, uint32_t key)
{
    size_t low = 0;
    size_t high = data.size();
    size_t mid = high >> 1;
    // size_t ind = -1;
    size_t size = high;
    while (size > 0xFF) {
        // auto mid = (low + high) >> 1;
        auto val = data[mid]->id_;
        if (key > val) {
            low = mid + 1;
        } else {
            if (key == val) {
                return mid;
                // low = high = mid;
                // break;
            }
            high = mid;
        }
        mid = (low + high) >> 1;
        size = high - low;
    }
    return -1;
    //return bsearch(data, key, low, high);
}

size_t b2search_original(const std::vector<A*>& records, uint32_t rid)
{
    int64_t numrecs = records.size();
    int64_t foundindex = -1;
    int64_t a = 0;
    int64_t z = numrecs - 1;
    int64_t k;
    uint32_t id;

    if (numrecs == 0)
        return records.size();
    while (a < z) {
        k = a + ((z - a) >> 1);

        id = records[k]->id_;

        if (id > rid) {
            z = k;
        } else if (id < rid) {
            a = k + 1;
        } else {
            foundindex = k;
            break;
        }

        if (z - a <= 1) {
            if (a == 0 && records[a]->id_ == rid)
                foundindex = a;
            if (z == numrecs - 1 && records[z]->id_ == rid)
                foundindex = z;
            break;
        }
    }

    return foundindex != -1 ? foundindex : records.size();
}

/*
BENCHMARK_DEFINE_F(VectorSearchFixture, lookup_table)(benchmark::State& state)
{
	size_t itr{}, sum{};
	for (auto _ : state) {
		auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
		benchmark::DoNotOptimize(sum += records_[rid]->id_);
	}
}
BENCHMARK_REGISTER_F(VectorSearchFixture, lookup_table)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();
*/

size_t fast_upper_bound4(const std::vector<A*>& data, uint32_t value)
{
    size_t size = data.size();
    size_t low = 0;

    while (size >= 16) {
        size_t half = size / 2;
        size_t other_half = size - half;
        size_t probe = low + half;
        size_t other_low = low + other_half;
        auto v = data[probe]->GetId();
        size = half;
        low = v >= value ? low : other_low;

        half = size / 2;
        other_half = size - half;
        probe = low + half;
        other_low = low + other_half;
        v = data[probe]->GetId();
        size = half;
        low = v >= value ? low : other_low;

        half = size / 2;
        other_half = size - half;
        probe = low + half;
        other_low = low + other_half;
        v = data[probe]->GetId();
        size = half;
        low = v >= value ? low : other_low;
    }

    while (size > 0) {
        size_t half = size / 2;
        size_t other_half = size - half;
        size_t probe = low + half;
        size_t other_low = low + other_half;
        auto v = data[probe]->GetId();
        size = half;
        low = v >= value ? low : other_low;
    }

    return low;
}

size_t fibMonaccianSearch(const std::vector<A*>& arr, int x)
{
    int64_t n = arr.size() - 1;
    /* Initialize fibonacci numbers */
    int64_t fibMMm2 = 0; // (m-2)'th Fibonacci No.
    int64_t fibMMm1 = 1; // (m-1)'th Fibonacci No.
    int64_t fibM = fibMMm2 + fibMMm1; // m'th Fibonacci

    /* fibM is going to store the smallest Fibonacci 
								   Number greater than or equal to n */
    while (fibM <= n) {
        fibMMm2 = fibMMm1;
        fibMMm1 = fibM;
        fibM = fibMMm2 + fibMMm1;
    }

    // Marks the eliminated range from front
    int64_t offset = -1;

    /* while there are elements to be inspected. Note that 
	we compare arr[fibMm2] with x. When fibM becomes 1, 
	fibMm2 becomes 0 */
    while (fibM > 1) {
        // Check if fibMm2 is a valid location
        int64_t i = std::min(offset + fibMMm2, n - 1);

        /* If x is greater than the value at index fibMm2, 
		cut the subarray array from offset to i */
        if (arr[i]->GetId() < x) {
            fibM = fibMMm1;
            fibMMm1 = fibMMm2;
            fibMMm2 = fibM - fibMMm1;
            offset = i;
        }

        /* If x is greater than the value at index fibMm2, 
		cut the subarray after i+1  */
        else if (arr[i]->GetId() > x) {
            fibM = fibMMm2;
            fibMMm1 = fibMMm1 - fibMMm2;
            fibMMm2 = fibM - fibMMm1;
        }

        /* element found. return index */
        else
            return i;
    }

    /* comparing the last element with x */
    if (fibMMm1 && arr[offset + 1]->GetId() == x)
        return offset + 1;

    /*element not found. return -1 */
    return -1;
}

#include <iterator>

/**
* An iterator class capable of navigating across the Fibonacci sequence using
* a user-specified integer type.
*/
template <typename Integer>
class FibonacciIterator : public std::iterator<std::bidirectional_iterator_tag, const Integer> {
public:

    explicit FibonacciIterator() : curr(0), next(1) { }

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

//#include "FibonacciIterator.hh"
#include <functional> // For std::less
#include <iterator> // For std::distance, std::iterator_traits

/* Main implementation of Fibonacci search uses a Fibonacci iterator to access
* the consecutive Fibonacci numbers as efficiently as possible.
*/
template <typename RandomIterator, typename Value, typename Comparator>
RandomIterator FibonacciSearch(RandomIterator begin, RandomIterator end,
    const Value& value, Comparator comp)
{
    /* See what type we use to keep track of iterator distances, then use that
	* to build our Fibonacci iterator.
	*/
    typedef typename std::iterator_traits<RandomIterator>::difference_type Integer;

    /* Create a Fibonacci iterator so that we can quickly access Fibonacci values
	* in sequence.
	*/
    FibonacciIterator<Integer> itr(std::distance(begin, end));

    /* Find the smallest Fibonacci number greater than the number of elements in
	* the range, then back it up one step.  This loop always executes at least
	* once, since when the iterator starts off it has value zero, which can't be
	* bigger than the number of elements.  For this reason, we write this as a
	* do ... while loop to make it clearer what's going on.
	*/
    //    const auto s = std::distance(begin, end);
    //do
    //  ++itr;
    //while (*itr <= s);

    /* Back up a step, which is well-defined because we took at least one
	* step.
	*/
    //--itr;

    /* Until the iterator has reached zero (meaning that there's nothing left to
	* check), apply the Fibonacci search to the range.
	*/
    while (*itr > 0) {
        /* See if this element is in range.  If it isn't, then we need to back the
		* iterator up a step.  Note that we allow *itr to be equal to end - begin
		* because we always read right before the given index (see below).
		*/
        //if (*itr > std::distance(begin, end)) {
        //--itr;
        //continue;
        //}
        /* Otherwise, compare the element at the given index to the value we're
		* looking for.  If the current element is larger, then we look in the
		* first half of the array.
		*
		* Notice that we compare the value at begin[*itr - 1] rather than at
		* begin[*itr].  This is because the array is zero-indexed, but *itr is
		* giving us back one-indexed locations.
		*/
        //auto mid = std::next(begin, *itr);

        //        const auto& c =  begin[*itr - 1];

        if (comp(value, begin[*itr - 1])) {
            --itr;
        }
        /*/
		else {
		if (comp(begin[*itr], value)) {
		std::advance(begin, *itr);
		--itr;
		--itr;
		auto n = std::distance(begin, end);
		while (*itr >= n) {
		--itr;
		}
		} else {
		return std::next(begin, *itr);
		}
		}

		/*/
        else if (comp(begin[*itr - 1], value)) {
            std::advance(begin, *itr);
            std::advance(itr, -2);
            while (*itr > std::distance(begin, end)) {
                --itr;
            }
            //--itr;
            //--itr;
            //const auto n = std::distance(begin, end);
            //while (*itr > n) {
            //	--itr;
            //}
        } else {
            std::advance(begin, *itr - 1);
            return begin;
            //return std::next(begin, *itr);
        }
        //*/
    }

    /* If we made it here, we didn't find the value in question. */
    return (begin != end && !(comp(*begin, value)) && !(comp(value, *begin))) ? begin : end;
}

/* Non-comparator version implemented in terms of the comparator version. */
template <typename RandomIterator, typename Value>
bool FibonacciSearch(RandomIterator begin, RandomIterator end,
    const Value& value)
{
    return FibonacciSearch(begin, end, value,
        std::less<typename std::iterator_traits<RandomIterator>::value_type>());
}

template <class F1, class F2>
struct overload_set : F1, F2 {
    overload_set(F1 f1, F2 f2)
        : F1(f1)
        , F2(f2)
    {
    }

    using F1::operator();
    using F2::operator();
};

template <class F1, class F2>
overload_set<F1, F2> overload(F1 f1, F2 f2)
{
    return overload_set<F1, F2>(f1, f2);
}

template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator FallbackInterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

	difference_type count = std::distance(begin, end);

#ifdef NDEBUG
	RandomIterator last = std::prev(end);
#else
	RandomIterator last = count ? std::prev(end) : end;
#endif
	int cnt = 6;
	while (count > 0 && --cnt) {

		if (!comp(*begin, key)) {
			end = !comp(key, *begin) ? begin : end;
			break;
		}

		if (!comp(key, *last)) {
			end = !comp(*last, key) ? last : end;
			break;
		}

		difference_type probe = static_cast<difference_type>(lerp(*begin, *last, key) * (count - 1));

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
	if (!cnt) {
		return FibonacciSearch(begin, end, key, std::move(comp));
	}
	return end;

}


template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator InterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

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

			//probe -= count;
			//std::advance(last, probe - count);
			//count = probe;

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
			//count -= ++probe;
			//std::advance(begin, probe);

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

template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator InterpolationSearchB(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

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

		difference_type probe = static_cast<difference_type>(double(lerp(*begin, *last, key)) * double(count - 1));

		if (comp(key, begin[probe])) {

			probe -= count;
			std::advance(last, probe);
			count += probe;

			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, probe - count);
				count = probe;
			}
		} else if (comp(begin[probe], key)) {
			count -= ++probe;
			std::advance(begin, probe);

			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, probe - count);
				count = probe;
			}
		} else {
			end = begin;
			std::advance(end, probe);
			break;
		}

	}
	return end;
}


template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator InterpolationSearch2(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter conv)
{
    using difference_type = std::iterator_traits<RandomIterator>::difference_type;

	difference_type count = std::distance(begin, end);

	//RandomIterator last = std::prev(end);

    while (count > 0) {
		
        if (!comp(*begin, key)) {
            begin = !comp(key, *begin) ? begin : end;
			break;
        }
		/*
        if (!comp(key, *(end - 1))) {
            begin = !comp(*(end - 1), key) ? (end -1) : end;
			break;

        difference_type probe = (double(key) - s) * (count - 1) / (e - s);
        //difference_type probe = std::distance(begin, end) >> 1;

        auto p = conv(begin[probe]);

        if (p < key) {
            std::advance(begin, ++probe);
            count -= probe;
        } else if (key < p) {
            std::advance(last, probe - count);
            count = probe;
        } else {
            std::advance(begin, probe);
            return begin;
        }
		*/
		difference_type probe = conv(*begin, *(end -1), key) * (count - 1);
		/*
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
		*/
		/*
		if (comp(key, begin[probe])) {
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				last = begin + --probe;
				//count = std::distance(begin, last);
				count -= count - probe;
			}
		} else if (comp(begin[probe], key)) {
			probe = count >> 1;
			if (!comp(begin[probe], key)) {
				last = begin + probe;
				//count = std::distance(begin, last);
				count -= count - probe;
			} else {
				std::advance(begin, ++probe);
				count -= probe;
			}
			*/
		/*
		if (comp(key, begin[probe])) {
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, probe - count);
				count = probe;
			}
		} else if (comp(begin[probe], key)) {
			probe = count >> 1;
			if (!comp(begin[probe], key)) {
				std::advance(last, probe - count);
				count = probe;
			} else {
				std::advance(begin, ++probe);
				count -= probe;
			}
		} else {
			end = std::next(begin, probe);
			break;
		}*/
		if (comp(key, begin[probe])) {
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(end, probe - count);
				count = probe;
			}
		} else if (comp(begin[probe], key)) {
			probe = count >> 1;
			if (!comp(begin[probe], key)) {
//				end = begin + probe;
//				count = std::distance(begin, end);;

				std::advance(end, probe - count);
				count = probe;
			} else {
				std::advance(begin, ++probe);
				count -= probe;
			}
		} else {
			std::advance(begin, probe);
			break;
		}
    }
    return begin;
}

template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator HybridInterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

	difference_type count = std::distance(begin, end), mid;

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

		difference_type probe = static_cast<difference_type>(lerp(*begin, *last, key) * (count - 1));

		if (comp(key, begin[probe])) {
			probe -= count;
			mid = count >> 1;
			count += probe;
			if (!comp(key, begin[mid])) {
				std::advance(begin, mid);
				std::advance(last, probe);
				count -= mid;
			} else {
				mid -= count;
				std::advance(last, probe + mid);
				count += mid;
			}
		} else if (comp(begin[probe], key)) {
			mid = (probe + count) >> 1;
			if (!comp(key, begin[mid])) {
				std::advance(begin, mid);
				count -= mid;
			} else {
				std::advance(begin, probe);
				std::advance(last, mid - count);
				count = mid - probe;
			}
		} else {
			end = std::next(begin, probe);
			break;
		}

	}
	return end;
}

int HybridSearch(const std::vector<A*>& a, uint32_t x)
{
    int left = 0, right = a.size() - 1;
    int Inter, Mid;
    while (left < right) {
        //double s = double((x - a[left]->id_)) * (right - left) / (a[right]->id_ - a[left]->id_);
        Inter = left + double((x - a[left]->id_)) * (right - left) / (a[right]->id_ - a[left]->id_);
        if (Inter > right || Inter < left)
            break;
        if (x > a[Inter]->id_) {
            Mid = (Inter + right) / 2;
            if (x <= a[Mid]->id_) {
                left = Inter + 1;
                right = Mid;

            } else {
                left = Mid + 1;
            }
        } else if (x < a[Inter]->id_) {
            Mid = (Inter + left) / 2;
            if (x >= a[Mid]->id_) {
                left = Mid;
                right = Inter - 1;
            } else {
                right = Mid - 1;
            }
        } else {
            return Inter;
        }
    }
    return -1;
}


BENCHMARK_DEFINE_F(VectorSearchFixture, combine_search)
(benchmark::State& state)
{
	uint64_t sum {}, itr {};
	for (auto _ : state) {
		auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
		auto low = combine_search(records_, rid);
		if (low < records_.size()) {
			benchmark::DoNotOptimize(sum += records_[low]->id_);
		}
	}
}
BENCHMARK_REGISTER_F(VectorSearchFixture, combine_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();


BENCHMARK_DEFINE_F(VectorSearchFixture, hib_search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = HybridSearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, hib_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();


size_t interpolationSearch(const std::vector<A*>& arr, int x)
{
    // Find indexes of two corners
    int64_t lo = 0, hi = (arr.size() - 1);

    // Since array is sorted, an element present
    // in array must be in range defined by corner
    while (lo <= hi && x >= arr[lo]->GetId() && x <= arr[hi]->GetId()) {
        if (lo == hi) {
            if (arr[lo]->GetId() == x)
                return lo;
            return -1;
        }
        // Probing the position with keeping
        // uniform distribution in mind.
        int pos = lo + (((double)(hi - lo) / (double(arr[hi]->GetId()) - arr[lo]->GetId())) * (x - arr[lo]->GetId()));
        if (pos < lo || pos > hi)
            break;
        // Condition of target found
        if (arr[pos]->GetId() == x)
            return pos;

        // If x is larger, x is in upper part
        if (arr[pos]->GetId() < x)
            lo = pos + 1;

        // If x is smaller, x is in the lower part
        else
            hi = pos - 1;
    }
    return -1;
}

BENCHMARK_DEFINE_F(VectorSearchFixture, fint_search)
(benchmark::State& state)
{
	uint64_t sum {}, itr {};
	//records_ // (double(key) - conv(*begin)) * (count - 1) / (conv(*last) - conv(*begin));
	for (auto _ : state) {
		auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
		auto low = FallbackInterpolationSearch(std::cbegin(records_), std::cend(records_), rid,
			Comp{},
			[](const A* first, const A* last, uint32_t key) -> double { return (double(key) - first->GetId()) / (last->GetId() - first->GetId()); });
		//[](const A* e) { return e->GetId(); });

		if (low != std::cend(records_)) {
			if ((*low)->id_ != rid) {
				std::cout << "error" << std::endl;
			}
			benchmark::DoNotOptimize(sum += (*low)->id_);
		} else {
			std::cout << "not found error" << std::endl;
		}
	}
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, fint_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, intb_search)
(benchmark::State& state)
{
	uint64_t sum {}, itr {};
	//records_ // (double(key) - conv(*begin)) * (count - 1) / (conv(*last) - conv(*begin));
	for (auto _ : state) {
		auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
		auto low = InterpolationSearch(std::cbegin(records_), std::cend(records_), rid,
			Comp{},
			[](const A* first, const A* last, uint32_t key) -> double { return (double(key) - first->GetId()) / (last->GetId() - first->GetId()); });
		//[](const A* e) { return e->GetId(); });

		if (low != std::cend(records_)) {
			if ((*low)->id_ != rid) {
				std::cout << "error" << std::endl;
			}
			benchmark::DoNotOptimize(sum += (*low)->id_);
		} else {
			std::cout << "not found error" << std::endl;
		}
	}
}
BENCHMARK_REGISTER_F(VectorSearchFixture, intb_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();



BENCHMARK_DEFINE_F(VectorSearchFixture, int_search)
(benchmark::State& state)
{
	uint64_t sum {}, itr {};
	//records_ // (double(key) - conv(*begin)) * (count - 1) / (conv(*last) - conv(*begin));
	for (auto _ : state) {
		auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
		auto low = InterpolationSearch(std::cbegin(records_), std::cend(records_), rid,
			Comp{},
			[](const A* first, const A* last, uint32_t key) -> double { return (double(key) - first->GetId()) / (last->GetId() - first->GetId()); });
		//[](const A* e) { return e->GetId(); });

		if (low != std::cend(records_)) {
			if ((*low)->id_ != rid) {
				std::cout << "error" << std::endl;
			}
			benchmark::DoNotOptimize(sum += (*low)->id_);
		} else {
			std::cout << "not found error" << std::endl;
		}
	}
}
BENCHMARK_REGISTER_F(VectorSearchFixture, int_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();


BENCHMARK_DEFINE_F(VectorSearchFixture, hint_search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
	//records_ // (double(key) - conv(*begin)) * (count - 1) / (conv(*last) - conv(*begin));
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
		auto low = HybridInterpolationSearch(std::cbegin(records_), std::cend(records_), rid,
			Comp{},
			[](const A* first, const A* last, uint32_t key) -> double { return (double(key) - first->GetId()) / (last->GetId() - first->GetId()); });
			//[](const A* e) { return e->GetId(); });

        if (low != std::cend(records_)) {
			if ((*low)->id_ != rid) {
				std::cout << "error" << std::endl;
			}
            benchmark::DoNotOptimize(sum += (*low)->id_);
        } else {
			std::cout << "not found error" << std::endl;
		}

    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, hint_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, fib3_search)
(benchmark::State& state)
{
    size_t itr {}, sum {};
    for (auto&& _ : state) {
        auto id = searches_[itr++ % state.range_x()];
        /*
		auto first = std::lower_bound(std::cbegin(records_), std::cend(records_), id,
		overload(
		[](const auto& lhs, const auto& i) -> bool { return lhs->id_ < i; },
		[](const auto& i, const auto& lhs) -> bool { return lhs->id_ < i; }
		));
		*/
        auto first = FibonacciSearch(std::cbegin(records_), std::cend(records_), id,
            overload(
                [](const A* lhs, const uint32_t i) -> bool const { return lhs->id_ < i; },
                [](const uint32_t i, const A* lhs) -> bool const { return i < lhs->id_; }));
        if (first != std::cend(records_)) {
            benchmark::DoNotOptimize(sum += (*first)->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, fib3_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();


BENCHMARK_DEFINE_F(VectorSearchFixture, int2_search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = interpolationSearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, int2_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, fib_search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = fibMonaccianSearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, fib_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, binary_cmov)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = binary_cmov(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, binary_cmov)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, nary_search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = nary_search(records_, rid, 3);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, nary_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, lsearch)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = lsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, lsearch)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, binary_fallback)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = fallback_bsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, binary_fallback)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, binary_cmov2)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = binary_cmov2(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, binary_cmov2)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, bsearch)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = bsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, bsearch)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, b2search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = b2search(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, b2search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, b2search_original)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = b2search_original(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, b2search_original)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

/*
BENCHMARK_F(VectorSearchFixture, test_vector_bsearch)(benchmark::State& state)
{
	uint64_t sum{};
	size_t k, a = 0, z = records_.size(), foundindex, numrecs = records_.size(), itr = 0;
	for (auto _ : state) {
		auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
		// auto first = std::lower_bound(std::begin(vect), std::end(vect), id, [](const A lhs, auto i) . bool { return
		// lhs->id_ < i; });
		a = 0;
		z = records_.size();
		while (z - a > 1) {
			k = (a + z) >> 1;

			uint32_t id = records_[k]->id_;
			if (id < rid) {
				a = k;
			} else if (id > rid) {
				z = k;
			} else {
				foundindex = k;
				break;
			}
		}

		if (z - a <= 1) {
			if (a == 0 && records_[a]->id_ == rid)
				foundindex = a;
			if (z == numrecs - 1 && records_[z]->id_ == rid)
				foundindex = z;
			//    break;
		}

		sum += records_[foundindex]->id_;
		//  }
		// std::cout << sum << std::endl;
	}
}
*/
BENCHMARK_DEFINE_F(VectorSearchFixture, fast_upperbound)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = fast_upper_bound4(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, fast_upperbound)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, std_equal_range)
(benchmark::State& state)
{
    size_t itr {}, sum {};
    for (auto _ : state) {
        auto id = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto first = std::equal_range(std::cbegin(records_), std::cend(records_), id, Comp {});
        if (first.first != first.second) {
            benchmark::DoNotOptimize(sum += (*first.first)->GetId());
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_equal_range)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, std_lower_bound)
(benchmark::State& state)
{
    size_t itr {}, sum {};
    for (auto _ : state) {
        auto id = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto first = std::lower_bound(std::cbegin(records_), std::cend(records_), id, [](const auto& e, auto id) -> bool { return e->GetId() < id; });
        if (first != std::cend(records_) && (*first)->GetId() == id) {
            benchmark::DoNotOptimize(sum += (*first)->GetId());
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_lower_bound)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

BENCHMARK_DEFINE_F(VectorSearchFixture, std_partition_point)
(benchmark::State& state)
{
    size_t itr {}, sum {};
    for (auto _ : state) {
        auto id = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto first = std::partition_point(std::cbegin(records_), std::cend(records_), [=](const auto& e) -> bool { return e->GetId() < id; });
        if (first != std::cend(records_) && (*first)->GetId() == id) {
            benchmark::DoNotOptimize(sum += (*first)->GetId());
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_partition_point)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity();

/*
void test_umap()
{
	const size_t test_size = 1000000;
	std::unordered_map<uint32_t, A> umap;
	umap.reserve(test_size);
	for (size_t i = 0; i < test_size; i++) {
		umap[static_cast<uint32_t>(i)] = new A(static_cast<uint32_t>(i));
	}
	std::random_device r;
	std::default_random_engine e1(r());
	std::uniform_int_distribution<int> uniform_dist(1, test_size - 1);
	uint64_t sum = 0;
	for (size_t i = 0; i < test_size; i++) {
		auto item = umap[uniform_dist(e1)];
		sum += item->id_;
	}
}
*/
BENCHMARK_MAIN();
