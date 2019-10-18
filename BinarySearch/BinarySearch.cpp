// BinarySearch.cpp : Defines the entry point for the application.
//

#include "BinarySearch.h"
#include <algorithm>
#include <iterator>
#include <vector>

//using namespace std;

namespace qb {


template <typename Integer>
class FibonacciIterator : public std::iterator<std::bidirectional_iterator_tag, const Integer> {
public:
	explicit FibonacciIterator()
		: curr(0)
		, next(1)
	{
	}

	explicit FibonacciIterator(Integer target)
		: curr(0)
		, next(1)
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


template <class ForwardIt, class Value, typename Comparator>
ForwardIt Search(ForwardIt begin, ForwardIt end, const Value& key, Comparator comp)
{
	using difference_type = std::iterator_traits<ForwardIt>::difference_type;

	difference_type count = std::distance(begin, end);

	while (difference_type half = count >> 1) {
		
		if (comp(begin[half], key)) {
			std::advance(begin, ++half);
		}
		else {
			if (!comp(key, begin[half])) {
				std::advance(begin, half);
				break;
			}
			//std::advance(end, half - count);
		}
		count -= half;
	}
	return (begin != end && !comp(*begin, key) && !comp(key, *begin)) ? begin : end;
}

}

size_t bsearch(const std::vector<uint32_t>& arr, uint32_t x)
{
    auto itr = qb::Search(std::cbegin(arr), std::cend(arr), x,
		std::less<uint32_t>{}//, [](uint32_t first, uint32_t last, uint32_t key) -> double { return (double(key) - first) / (last - first); }
		);
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
    std::vector<uint32_t> vec = { 10, 12, 14,16, 18, 110, 111};
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
	std::cout << "finish " << std::endl;
    while (1) {
        for (qb::FibonacciIterator<uint32_t> itr(0); *itr < 500000000; ++itr) {
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
    std::cout << "Hello CMake." << std::endl;
    return 0;
}
