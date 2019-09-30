// LUTSearch.cpp : Defines the entry point for the application.
//

#include "LUTSearch.h"
#include <thread>
#include <vector>
#include <algorithm>

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
FibonacciIterator<Integer, Plus, Minus>::FibonacciIterator(Integer target)
    : curr(0)
    , next(1)
{
    //Integer newNext = curr + next;
    while (next <= target) {
        //curr = next;
        //next = newNext;
        //newNext = curr + next;
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
        if (value > begin[*itr]) {
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
        else if (begin[*itr] > value) {
            std::advance(begin, *itr);
            --itr;
            --itr;
            const auto n = std::distance(begin, end);
            while (*itr >= n) {
                --itr;
            }
        } else
            return std::next(begin, *itr);
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



size_t fallback_bsearch(const std::vector<uint32_t>& data, uint32_t key, size_t low, size_t high)
{
    auto first = FibonacciSearch(std::crbegin(data), std::crend(data), key,
        overload(
            [](const uint32_t &lhs, const uint32_t i) -> bool const { return lhs < i; },
            [](const uint32_t i, const uint32_t &lhs) -> bool const { return i < lhs; }));
    if (first != std::crend(data)) {
        return std::distance(std::crbegin(data), first);
    }
    return -1;
}



int main()
{
    std::vector<uint32_t> rec = {0, 3, 5, 6, 10, 100};
    /*
	std::vector<uint32_t> rec(0xFFFF, 0);
    long i = 0;
    for (auto& a : rec) {
        a = ++i;
    }
    i += 100;
	for (; i > 0; i--) {
        auto j = fallback_bsearch(rec, i, 0, rec.size() - 1);
        //auto j = bbsearch(rec, i, 0, rec.size());
        if (i != j + 1) {
            std::cout << i << " " << j << std::endl;
        }
    }
      */

    std::cout << fallback_bsearch(rec, 0, 0, rec.size() - 1) << std::endl;
    std::cout << fallback_bsearch(rec, 3, 0, rec.size() -1 ) << std::endl;
    std::cout << fallback_bsearch(rec, 5, 0, rec.size() - 1) << std::endl;
    std::cout << fallback_bsearch(rec, 6, 0, rec.size() - 1) << std::endl;
    std::cout << fallback_bsearch(rec, 10, 0, rec.size() - 1) << std::endl;
    std::cout << fallback_bsearch(rec, 11, 0, rec.size() - 1) << std::endl;
    std::cout << fallback_bsearch(rec, 16, 0, rec.size() - 1) << std::endl;
    return 0;
}
