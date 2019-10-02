// BinarySearch.cpp : Defines the entry point for the application.
//

#include "BinarySearch.h"
#include <algorithm>
#include <vector>

using namespace std;

namespace qb {

template <class ForwardIt, class T>
inline ForwardIt branch_less_binary_find_n(ForwardIt first, size_t size, const T& key)
{
    using difference_type = std::iterator_traits<ForwardIt>::difference_type;

    while (difference_type half = size >> 1) {
        ForwardIt middle = first;
        std::advance(middle, half);
        if (!((*middle) > key)) {
            std::advance(first, half);
        }
        size = size - half;
    }
    return first;
}

template <class ForwardIt, class T>
inline ForwardIt branch_full_binary_find_n(ForwardIt first, ForwardIt last, size_t size, const T& key)
{
    using difference_type = std::iterator_traits<ForwardIt>::difference_type;
    /*
	while (first < last) {
	auto half = size / 2;
	ForwardIt middle = first;
	std::advance(middle, half);
	if ((*middle)->GetId() < key) {
	std::advance(first, half + 1);
	//first = ++middle;
	size -= half + 1;
	} else {
	if (!(key < (*middle)->GetId())) {
	return middle;// std::advance(first, half);
	//break;
	}
	std::advance(last, -half);
	//last = middle;
	size -= half;
	}
	//size -= half;
	}
	return first;
	*/
    //typename std::iterator_traits<ForwardIt>::difference_type half;
    typename std::iterator_traits<ForwardIt>::difference_type half;
    half = std::distance(first, last) >> 1;
    ForwardIt middle = first, ind = last;
    std::advance(middle, half);
    while (first < last) {
        auto val = (*middle);

        if (key > val) {
            first = middle;
            std::advance(first, 1);
        } else {
            if (key == val) {
                ind = middle;
                break;
            }
            last = middle;
        }
        half = std::distance(first, last) >> 1;
        middle = first;
        std::advance(middle, half);
    }
    return ind;
}

template <class ForwardIt, class T>
inline ForwardIt equal_range(ForwardIt first, ForwardIt last, const T& key)
{
    using difference_type = std::iterator_traits<ForwardIt>::difference_type;
    difference_type size = std::distance(first, last);
    return (size == 0) ? last : size < 0x0 ? branch_less_binary_find_n(first, size, key) : branch_full_binary_find_n(first, last, size, key);
}

};

size_t bfsearch(const std::vector<uint32_t>& arr, uint32_t x)
{
    //	auto ind = qb::equal_range(std::cbegin(vec), std::cend(vec), key);
    //	return ind != std::cend(vec) ? std::distance(std::cbegin(vec), ind) : -1;
    //	uint64_t fibMonaccianSearch(const std::vector<const SBrecord*>& arr, uint32_t x)
    //	{
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
        if (arr[i] < x) {
            fibM = fibMMm1;
            fibMMm1 = fibMMm2;
            fibMMm2 = fibM - fibMMm1;
            offset = i;
        }

        /* If x is greater than the value at index fibMm2, 
			cut the subarray after i+1  */
        else if (arr[i] > x) {
            fibM = fibMMm2;
            fibMMm1 = fibMMm1 - fibMMm2;
            fibMMm2 = fibM - fibMMm1;
        }

        /* element found. return index */
        else
            return i;
    }

    /* comparing the last element with x */
    //std::cout << "o" << offset << std::endl;
    if (fibMMm1 && arr[offset + 1] == x)
        return offset + 1;

    /*element not found. return -1 */
    return -1;
    //	}
}

#include <functional> // For std::plus, std::minus
#include <iterator> // For std::bidirectional_iterator_tag, std::iterator

/**
* An iterator class capable of navigating across the Fibonacci sequence using
* a user-specified integer type.
*/
template <typename Integer, typename Plus = std::plus<Integer>,
    typename Minus = std::minus<Integer>>
class FibonacciIterator : public std::iterator<std::bidirectional_iterator_tag,
                              const Integer> {
public:
    /**
		* Constructor: FibonacciIterator(Integer zero = Integer(0),
		*                                Integer one = Integer(1),
		*                                Plus p = Plus(), Minus m = Minus())
		* Usage: FibonacciIterator<int> itr;
		* --------------------------------------------------------------------------
		* Constructs a new Fibonacci iterator traversing the Fibonacci sequence
		* whose first two terms are zero and one and that uses the specified plus
		* and minus function objects to navigate the sequence.
		*/
    // explicit FibonacciIterator(Integer zero = Integer(0),
    //   Integer one = Integer(1));
    //Plus p = Plus(), Minus m = Minus());
    explicit FibonacciIterator(Integer tasrget);

    /**
		* operator*  () const;
		* operator-> () const;
		* Usage: cout << *itr << endl;
		* --------------------------------------------------------------------------
		* Dereferences and returns the current integer in the sequence.  You should
		* not modify the values returned as they are not guaranteed to be valid
		* after the iterator advances.  Moreover, you should not hold pointers or
		* references to these values, as the memory will be recycled after the
		* iterator is incremented or decremented.
		*/
    const Integer& operator*() const;
    const Integer* operator->() const;

    /**
		* operator++ ();
		* operator++ (int);
		* operator-- ();
		* operator-- (int);
		* Usage: ++itr; --itr; itr++; itr--;
		* --------------------------------------------------------------------------
		* Moves the iterator one step forward or backward in the Fibonacci sequence.
		* If integer overflow occurs, the results depend on the type of the integer
		* being used as a counter.  If the iterator is backed up while at 0, the
		* results are mathematically well-defined but depend on the underlying type
		* of the integer for correctness.
		*/
    FibonacciIterator& operator++();
    const FibonacciIterator operator++(int);

    FibonacciIterator& operator--();
    const FibonacciIterator operator--(int);

private:
    /* The current and next Fibonacci values in the sequence. */
    Integer curr, next;

    /* The plus and minus operators. */
    //Plus plus;
    //Minus minus;
};

/* Comparison functions for FibonacciIterator. */
template <typename Integer, typename Plus, typename Minus>
constexpr bool operator==(const FibonacciIterator<Integer, Plus, Minus>& lhs,
    const FibonacciIterator<Integer, Plus, Minus>& rhs);
template <typename Integer, typename Plus, typename Minus>
constexpr bool operator!=(const FibonacciIterator<Integer, Plus, Minus>& lhs,
    const FibonacciIterator<Integer, Plus, Minus>& rhs);

/* * * * * Implementation Below This Point * * * * */

/* Constructor sets up the internal fields based on the parameters. */
//template <typename Integer, typename Plus, typename Minus>
//FibonacciIterator<Integer, Plus, Minus>::FibonacciIterator(Integer zero,
//Integer one
//Plus plus,
//    Minus minus
//)
//  : curr(zero)
//, next(one)
//, plus(plus)
//, minus(minus)
//{
// Handled in initializer list.
///}

template <typename Integer, typename Plus, typename Minus>
FibonacciIterator<Integer, Plus, Minus>::FibonacciIterator(Integer target)
    : curr(0)
    , next(1)
{
    //Integer newNext = curr + next;
    //Integer newNext;
    while (next <= target) {
        //newNext = curr + next;
        //curr = next;
        //next = newNext;

        next = curr + next;
        curr = next - curr;
    }
}

/* Dereferencing to a value just returns the current value in the sequence. */
template <typename Integer, typename Plus, typename Minus>
const Integer& FibonacciIterator<Integer, Plus, Minus>::operator*() const
{
    return curr;
}
template <typename Integer, typename Plus, typename Minus>
const Integer* FibonacciIterator<Integer, Plus, Minus>::operator->() const
{
    return &**this;
}

/* Incrementing the Fibonacci iterator walks forward one step in the Fibonacci
* series.
*/
template <typename Integer, typename Plus, typename Minus>
FibonacciIterator<Integer, Plus, Minus>&
FibonacciIterator<Integer, Plus, Minus>::operator++()
{
    //Integer newNext = curr + next;
    //curr = next;
    //next = newNext;
    next = curr + next;
    curr = next - curr;

    return *this;
}
template <typename Integer, typename Plus, typename Minus>
const FibonacciIterator<Integer, Plus, Minus>
FibonacciIterator<Integer, Plus, Minus>::operator++(int)
{
    FibonacciIterator result = *this;
    ++*this;
    return result;
}

/* Decrementing the Fibonacci iterator backs it up one step in the sequence. */
template <typename Integer, typename Plus, typename Minus>
FibonacciIterator<Integer, Plus, Minus>&
FibonacciIterator<Integer, Plus, Minus>::operator--()
{
    //Integer prev = next  - curr;
    //next = curr;
    //curr = prev;
    curr = next - curr;
    next = next - curr;

    return *this;
}
template <typename Integer, typename Plus, typename Minus>
const FibonacciIterator<Integer, Plus, Minus>
FibonacciIterator<Integer, Plus, Minus>::operator--(int)
{
    FibonacciIterator result = *this;
    --*this;
    return result;
}

/* Equality comparisons just check if the two values are equal. */
template <typename Integer, typename Plus, typename Minus>
constexpr bool operator==(const FibonacciIterator<Integer, Plus, Minus>& lhs,
    const FibonacciIterator<Integer, Plus, Minus>& rhs)
{
    return lhs.curr == rhs.curr; //    &&lhs.next == rhs.next;
}

/* Disequality implemented in terms of equality. */
template <typename Integer, typename Plus, typename Minus>
constexpr bool operator!=(const FibonacciIterator<Integer, Plus, Minus>& lhs,
    const FibonacciIterator<Integer, Plus, Minus>& rhs)
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
    typedef typename std::iterator_traits<RandomIterator>::difference_type Integer;

    FibonacciIterator<Integer> itr(std::distance(begin, end));

    while (*itr > 0) {
        if (comp(value, begin[*itr - 1])) {
            --itr;
        } else if (comp(begin[*itr - 1], value)) {
            std::advance(begin, *itr);
            std::advance(itr, -2);
            while (*itr > std::distance(begin, end)) {
                --itr;
            }
        } else {
            std::advance(begin, *itr - 1);
            return begin;
        }
    }
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

size_t bfibsearch(const std::vector<uint32_t>& records_, uint32_t id)
{
    auto first = FibonacciSearch(std::cbegin(records_), std::cend(records_), id,
        std::less<uint32_t> {});
    if (first != std::cend(records_)) {
        return std::distance(std::cbegin(records_), first);
    }
    return -1;
}

size_t bintsearch(const std::vector<uint32_t>& arr, uint32_t x)
{
    // Find indexes of two corners
    int64_t lo = 0, hi = (arr.size() - 1);

    // Since array is sorted, an element present
    // in array must be in range defined by corner
    while (lo <= hi && x >= arr[lo] && x <= arr[hi]) {
        if (lo == hi) {
            if (arr[lo] == x)
                return lo;
            return -1;
        }
        // Probing the position with keeping
        // uniform distribution in mind.
        int pos = lo + (((double)(hi - lo) / (arr[hi] - arr[lo])) * (x - arr[lo]));
        if (pos < lo || pos > hi)
            break;
        // Condition of target found
        if (arr[pos] == x)
            return pos;

        // If x is larger, x is in upper part
        if (arr[pos] < x)
            lo = pos + 1;

        // If x is smaller, x is in the lower part
        else
            hi = pos - 1;
    }
    return -1;
}

template <typename RandomIterator, typename Value, typename Converter>
RandomIterator InterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Converter conv)
{
    using difference_type = std::iterator_traits<RandomIterator>::difference_type;
    RandomIterator last = end;
    difference_type count = std::distance(begin, end);

    while (count > 0) {
        auto e = conv(*(end - 1));
        auto s = conv(*begin);
        if (key < s || e < key) {
            break;
        } else if (e - s == 0) {
            return begin;
        }

        difference_type probe = (double(key) - s) * (count - 1) / (e - s);

        auto p = conv(begin[probe]);

        if (p < key) {
            std::advance(begin, ++probe);
            count -= probe;
        } else if (key < p) {
            std::advance(end, probe - count);
            count = probe;
        } else {
            std::advance(begin, probe);
            return begin;
        }
    }
    return last;
}

size_t bsearch(const std::vector<uint32_t>& arr, uint32_t x)
{
    auto itr = InterpolationSearch(std::cbegin(arr), std::cend(arr), x, [](const auto& e) { return e; });
    //auto itr = FibonacciSearch(std::cbegin(arr), std::cend(arr), x, std::less<uint32_t>{});
    if (itr != std::cend(arr)) {
        return std::distance(std::cbegin(arr), itr);
    }
    return -1;
}

#include <random>

int main()
{
    //std::vector<uint32_t> vec = {10, 11, 11, 12, 18, 110, 111};
    //std::vector<uint32_t> vec = {10, 11, 12, 14, 16, 18, 110};
    std::vector<uint32_t> vec = { 10, 11, 12, 14, 16, 18, 120 };
    std::cout << bsearch(vec, 1) << std::endl;
    std::cout << bsearch(vec, 10) << std::endl;
    std::cout << bsearch(vec, 11) << std::endl;
    std::cout << bsearch(vec, 12) << std::endl;
    std::cout << bsearch(vec, 14) << std::endl;
    std::cout << bsearch(vec, 16) << std::endl;
    std::cout << bsearch(vec, 18) << std::endl;
    std::cout << bsearch(vec, 110) << std::endl;
    std::cout << bsearch(vec, 120) << std::endl;
    std::cout << bsearch(vec, 111) << std::endl;
    std::cout << bsearch(vec, 112) << std::endl;
    while (1) {
        for (FibonacciIterator<uint32_t> itr(0); *itr < 500000000; ++itr) {
            std::vector<uint32_t> vec(*itr, 0);
            std::mt19937 rng;
            rng.seed(std::random_device()());
            std::uniform_int_distribution<uint32_t> dist(*itr / 2, *itr * 2);
            for (auto& e : vec) {
                e = dist(rng);
            }
            std::sort(std::begin(vec), std::end(vec));
            vec.erase(std::unique(std::begin(vec), std::end(vec)), std::end(vec));
            std::sort(std::begin(vec), std::end(vec));
            std::cout << "test: " << *itr << " " << vec.size() << std::endl;
            for (auto& e : vec) {
                if (vec[bsearch(vec, e)] != e) {
                    std::cout << "error: " << e << " :";
                    for (auto& x : vec) {
                        std::cout << x << ",";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
    cout << "Hello CMake." << endl;
    return 0;
}
