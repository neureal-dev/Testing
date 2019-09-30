// Test.cpp : Defines the entry point for the application.
//

#include "LUTBenchmark.h"
#include <benchmark/benchmark.h>
#include <algorithm>
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
    A() : id_(0) {}
    A(uint32_t id) : id_(id) {}
    uint32_t getId() const { return id_; }

    uint32_t id_;
    uint64_t id1_;
    uint64_t id2_;
    uint64_t id3_;
    uint64_t id4_;
};

class VectorSearchFixture : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State& state)
    {
        records_ = std::vector<A*>(state.range_x());
        searches_ = std::vector<uint32_t>(state.range_x(), 0);
        benchmark::DoNotOptimize(records_.data());
        benchmark::DoNotOptimize(searches_.data());
        std::mt19937_64 rng;
        rng.seed(std::random_device()());
        std::uniform_int_distribution<uint32_t> uniform_dist(65536, static_cast<uint32_t>(state.range_x()) * 2 + 65536);

        for (int32_t i = state.range_x() - 1; i >= 0; --i) {
            records_[i] = new A(i + i + 65536);
            searches_[i] = uniform_dist(rng);
        }
        for (int64_t i = state.range_x() - 1; i >= 0; i-=3) {
            auto tmp = records_[i];
            records_[i] = new A(i + i + 65536);
            delete tmp;
            //searches_[i] = uniform_dist(rng);
        }

    }

    void TearDown(const ::benchmark::State& state)
    {
        for (auto item : records_) {
            delete item;
        }
        records_.clear();
    }

//    private:
    std::vector<A*> records_;
    std::vector<uint32_t> searches_;
};

enum locate_t { EQUAL, LEFT, RIGHT };

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


size_t binary_cmov(const std::vector<A*>& records, uint32_t rid)
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

size_t binary_cmov2(std::vector<A*>& records, uint32_t const rid)
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

size_t b2search(std::vector<A*>& records, uint32_t const rid)
{
    size_t size = records.size() - 1;
    size_t low = 0;
    //size_t middle;
    /*
    while (size_t half = size >> 1) {
        //size_t middle = low + half;
        //low = (records[middle]->id_ <= rid) ? middle : low;
        //size = size - half;

        low += half;
        low = records[low]->id_ <= rid ? low : low - half;
        size = size - half;

    }
    */
    while (low != size) {
        auto n = low + size;
        auto middle = n / 2 + (n & 1);
        if (records[middle]->id_ > rid) {
            size = middle - 1;
        } else {
            low = middle;
        }
    }
    return (size && records[low]->id_ == rid) ? low : -1;
}

inline size_t bb2search(std::vector<A*> const& records, uint32_t const rid, size_t low, size_t high)
{
    //high -= low;
    while (size_t half = high >> 1) {
        size_t middle = low + half;
        auto val = records[middle]->id_;
        low = (val <= rid) ? middle : low;
        high = high - half;
    }
    return (high && records[low]->id_ == rid) ? low : -1;
}


int lsearch(const std::vector<A*>& records, uint32_t key)
{
    auto itr = std::find_if(std::cbegin(records), std::cend(records), [=](const auto& e) -> bool { return key == e->id_; });
    return itr == std::cend(records) ? -1 : std::distance(std::cbegin(records), itr);
}

inline size_t bbsearch(std::vector<A*> const& resource, uint32_t key, size_t low, size_t high)
{
    size_t ind = -1;
    
    //high >> 1;
    
 
    return ind;
}

/* Returns index of x if present,  else returns -1 */
int fibMonaccianSearch(const std::vector<A*>& arr, uint32_t x)
{
    size_t n = arr.size()-1;
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
        int64_t i = std::min<int64_t>(offset + fibMMm2, n - 1);

        /* If x is greater than the value at index fibMm2, 
           cut the subarray array from offset to i */
        if (arr[i]->id_ < x) {
            fibM = fibMMm1;
            fibMMm1 = fibMMm2;
            fibMMm2 = fibM - fibMMm1;
            offset = i;
        }

        /* If x is greater than the value at index fibMm2, 
           cut the subarray after i+1  */
        else if (arr[i]->id_ > x) {
            fibM = fibMMm2;
            fibMMm1 = fibMMm1 - fibMMm2;
            fibMMm2 = fibM - fibMMm1;
        }

        /* element found. return index */
        else
            return i;
    }

    /* comparing the last element with x */
    if (fibMMm1 && arr[offset + 1]->id_ == x)
        return offset + 1;

    /*element not found. return -1 */
    return -1;
} 

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
/*
 * If val is found in arr, return the index of its location in arr.
 * Otherwise, return the index of the smallest element greater than val
 */
static int binsrch_geq(const int* arr, int n, int val)
{
    int low, high, mid;
    int geq;

    low = 0;
    high = n - 1;
    geq = -1;

    /* binary search for finding the element with value val */
    while (low <= high) {
        mid = (low + high) >> 1; //(low+high)/2;
        if (val < arr[mid]) {
            high = mid - 1;
            geq = mid;
        } else if (val > arr[mid])
            low = mid + 1;
        else
            return mid; /* found */
    }

    return geq; /* not found */
}

/*
  Fibonaccian search for locating the index of "val" in an array "arr" of size "n"
  that is sorted in ascending order.
 
  Algorithm description
  -----------------------------------------------------------------------------
  Let Fk represent the k-th Fibonacci number where Fk+2=Fk+1 + Fk for k>=0 and
  F0 = 0, F1 = 1. To test whether an item is in a list of n = Fm ordered numbers,
  proceed as follows:
  a) Set k = m.
  b) If k = 0, finish - no match.
  c) Test item against entry in position Fk-1.
  d) If match, finish.
  e) If item is less than entry Fk-1, discard entries from positions Fk-1 + 1 to n.
     Set k = k - 1 and go to b).
  f) If item is greater than entry Fk-1, discard entries from positions 1 to Fk-1.
     Renumber remaining entries from 1 to Fk-2, set k = k - 2 and go to b)
 
  If Fm>n then the original array is augmented with Fm-n numbers larger
  than val and the above algorithm is applied.
 */

int fibsrch(const std::vector<A*>& arr, int val)
{
    int n = arr.size() - 1;
    int k, idx, offs;
    int prevn = -1, prevk = -1;

    /*  Precomputed Fibonacci numbers F0 up to F46. This implementation assumes that the size n
     *  of the input array fits in 4 bytes. Note that F46=1836311903 is the largest Fibonacci
     *  number that is less or equal to the 4-byte INT_MAX (=2147483647). The next Fibonacci
     *  number, i.e. F47, is 2971215073 and is larger than INT_MAX, implying that it does not
     *  fit in a 4 byte integer. Note also that the last array element is INT_MAX rather than
     *  F47. This ensures correct operation for n>F46.
     */
    const  int Fib[47 + 1] = { 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765,
        10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578,
        5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141, 267914296,
        433494437, 701408733, 1134903170, 1836311903, INT_MAX };

    /* find the smallest fibonacci number that is greater or equal to n. Store this
     * number to avoid recomputing it in the case of repetitive searches with identical n.
     */
    if (n != prevn) {
#if 0
        k = (n > 1) ? binsrch_geq(Fib, sizeof(Fib) / sizeof(int), n) : 1;
#else /* the above binary search call can be substituted by the following code fragment: */
        {
            int f0, f1, t;

            for (f0 = 0, f1 = 1, k = 1; f1 < n; t = f1, f1 += f0, f0 = t, ++k)
                ;
        }
#endif
        prevk = k;
        prevn = n;
    } else
        k = prevk;

    /* If the sought value is larger than the largest Fibonacci number less than n,
     * care must be taken top ensure that we do not attempt to read beyond the end
     * of the array. If we do need to do this, we pretend that the array is padded
     * with elements larger than the sought value.
     */
    for (offs = 0; k > 0;) {
        idx = offs + Fib[--k];

        /* note that at this point k  has already been decremented once */
        if (idx >= n || val < arr[idx]->getId()) // index out of bounds or val in 1st part
            continue;
        else if (val > arr[idx]->getId()) // val in 2nd part
        {
            offs = idx;
            --k;
        } else // val==arr[idx], found
            return idx;
    }

    return -1; // not found
}

//template <typename RandomIterator, typename Value>
//bool FibonacciSearch(RandomIterator begin, RandomIterator end,
    //const Value& value);

/**
 * bool FibonacciSearch(RandomIterator begin, RandomIterator end, 
 *                      const Value& value, Comparator comp);
 * Usage: if (FibonacciSearch(v.begin(), v.end(), value, myComparator)) ...
 * ----------------------------------------------------------------------------
 * Uses the Fibonacci search technique to scan the range [begin, end) for the
 * given value, returning whether or not it was found.  The elements are
 * assumed to be sorted in ascending order according to comp.
 */
//template <typename RandomIterator, typename Value, typename Comparator>
//bool FibonacciSearch(RandomIterator begin, RandomIterator end,
    //const Value& value, Comparator comp);

/* * * * * Implementation Below This Point * * * * */
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
FibonacciIterator<Integer, Plus, Minus>::FibonacciIterator(Integer target) : curr(0), next(1)
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
    return lhs.curr == rhs.curr;//    &&lhs.next == rhs.next;
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
        if (comp(value, begin[*itr])) {
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
        else if (comp(begin[*itr], value)) {
            std::advance(begin, *itr);
            --itr;
            --itr;
            const auto n = std::distance(begin, end);
            while (*itr >= n) {
                --itr;
            }
        } else {
            std::advance(begin, *itr);
            return begin;
            //return std::next(begin, *itr);
        }
        //*/
    }

    /* If we made it here, we didn't find the value in question. */
    return end;
}

/* Non-comparator version implemented in terms of the comparator version. */
template <typename RandomIterator, typename Value>
bool FibonacciSearch(RandomIterator begin, RandomIterator end,
    const Value& value)
{
    return FibonacciSearch(begin, end, value,
        std::less<typename std::iterator_traits<RandomIterator>::value_type>());
}



size_t bsearch(const std::vector<A*>& data, uint32_t const key)
{
    size_t low = 0;
    
	size_t high = data.size();
    size_t ind = -1;
    size_t mid = high >> 1;

    while (low < high) {
        const auto val = data[mid]->id_;
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
    while (size > 0xFFFFF) {
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
    return bb2search(data, key, low, size);
}

size_t combine_search(std::vector<A*> const& data, uint32_t const rid)
{
    return (data.size() > 0xFFFFFull) ? bbsearch(data, rid, 0ull, data.size()) : bb2search(data, rid, 0ull, data.size());
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


   /* log_2 ceiling */
#ifdef _MSC_VER

#define FORCEINLINE __forceinline

#define NOINLINE __declspec(noinline)

#define ALIGN(n) __declspec(align(n))

FORCEINLINE uint32_t bsr(uint32_t x)
{

    unsigned long res;

    _BitScanReverse(&res, x);

    return res;
}

FORCEINLINE uint32_t bsf(uint32_t x)
{

    unsigned long res;

    _BitScanForward(&res, x);

    return res;
}

#else

#define FORCEINLINE __attribute__((always_inline)) inline

#define NOINLINE __attribute__((noinline))

#define ALIGN(n) __attribute__((aligned(n)))

FORCEINLINE uint32_t bsr(uint32_t x)
{

    return 31 - __builtin_clz(x);
}

FORCEINLINE uint32_t bsf(uint32_t x)
{

    return __builtin_ctz(x);
}

#endif

size_t bl3search(const std::vector<A*>& records, uint32_t key)
{
    //std::cout << bsr(10) << std::endl;
    //unsigned* low = vector;
    size_t low = 0, size = records.size();
    for (unsigned i = bsr(size); i != 0; i--) {
        size /= 2;
        unsigned mid = low + size;
        if (records[mid]->id_ <= key)
            low += size;
    }

    return records[low]->id_ == key ? -1 : low;
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

//#include "sse-binsearch-block.h"
#include <smmintrin.h>

#ifdef _MSC_VER
#include <intrin.h>

static uint32_t __inline ctz(uint32_t x)
{
    unsigned long r = 0;
    _BitScanReverse(&r, x);
    return r;
}

static uint32_t __inline clz(uint32_t x)
{
    unsigned long r = 0;
    _BitScanForward(&r, x);
    return r;
}
#endif

int SSEBinSearchBlock(const std::vector<A*>& data, uint32_t key)
{

    const __m128i keys = _mm_set1_epi32(key);
    __m128i v;

    int limit = data.size() - 1;
    int a = 0;
    int b = limit;

    while (a <= b) {
        const int c = (a + b) / 2;

        if (data[c]->getId() == key) {
            return c;
        }

        if (key < data[c]->getId()) {
            b = c - 1;

            if (b >= 4) {
                v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[b - 4]));
                v = _mm_cmpeq_epi32(v, keys);
                const uint16_t mask = _mm_movemask_epi8(v);
                if (mask) {
                    return b - 4 + ctz(mask) / 4;
                }
            }
        } else {
            a = c + 1;

            if (a + 4 < limit) {
                v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[a]));
                v = _mm_cmpeq_epi32(v, keys);
                const uint16_t mask = _mm_movemask_epi8(v);
                if (mask) {
                    return a + ctz(mask) / 4;
                }
            }
        }
    }

    return -1;
}

/*
BENCHMARK_DEFINE_F(VectorSearchFixture, lookup_table)(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        benchmark::DoNotOptimize(sum += records_[rid]->id_);
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, lookup_table)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);
*/
/*
BENCHMARK_DEFINE_F(VectorSearchFixture, binary_cmov)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = binary_cmov(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, binary_cmov)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, nary_search)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = nary_search(records_, rid, 3);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, nary_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, lsearch)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = lsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, lsearch)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, binary_cmov2)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = binary_cmov2(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, binary_cmov2)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


*/
/*
BENCHMARK_F(VectorSearchFixture, test_vector_bsearch)(benchmark::State& state)
{
    uint64_t sum{};
    size_t k, a = 0, z = records_.size(), foundindex, numrecs = records_.size(), itr = 0;
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
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
BENCHMARK_REGISTER_F(VectorSearchFixture, fib3_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);

BENCHMARK_DEFINE_F(VectorSearchFixture, bb_search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = bbsearch(records_, rid, 0, records_.size() - 1);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, bb_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, sse_search)
(benchmark::State& state)
{

    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto first = SSEBinSearchBlock(records_, rid);

            if (first != -1) {
                benchmark::DoNotOptimize(sum += records_[first]->id_);
            }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, sse_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, fib2_search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = fibsrch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
//BENCHMARK_REGISTER_F(VectorSearchFixture, fib2_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, fib_search)(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = fibMonaccianSearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, fib_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, b2search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = b2search(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, b2search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);

BENCHMARK_DEFINE_F(VectorSearchFixture, bsearch)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = bsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, bsearch)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);



BENCHMARK_DEFINE_F(VectorSearchFixture, std_lower_bound)(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto&& _ : state) {
        auto id = searches_[itr++ % state.range_x()];
        /*
        auto first = std::lower_bound(std::cbegin(records_), std::cend(records_), id,
                                      overload(
                                          [](const auto& lhs, const auto& i) -> bool { return lhs->id_ < i; },
                                          [](const auto& i, const auto& lhs) -> bool { return lhs->id_ < i; }
                                          ));
        */
        auto first = std::equal_range(std::cbegin(records_), std::cend(records_), id,
            overload(
                [](const A* lhs, const uint32_t i) -> bool const { return lhs->id_ < i; },
                [](const uint32_t i, const A* lhs) -> bool const { return i < lhs->id_; }));

        if (first.first != first.second) {
            benchmark::DoNotOptimize(sum += (*first.first)->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_lower_bound)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);


BENCHMARK_DEFINE_F(VectorSearchFixture, binary_fallback)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = fallback_bsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, binary_fallback)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);

BENCHMARK_DEFINE_F(VectorSearchFixture, bl3search)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = bl3search(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, bl3search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);

BENCHMARK_DEFINE_F(VectorSearchFixture, b2search_original)
(benchmark::State& state)
{
    uint64_t sum {}, itr {};
    for (auto _ : state) {
        auto rid = searches_[itr++ % state.range_x()];
        auto low = b2search_original(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low]->id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, b2search_original)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1)->Complexity()->MinTime(10);
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
