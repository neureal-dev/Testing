// BinarySearch.cpp : Defines the entry point for the application.
//

#include "BinarySearch.h"
#include <algorithm>
#include <array>
#include <iterator>
#include <string>
#include <vector>
#include <functional>

using charstr = std::basic_string<unsigned char>;

const std::array<char, 5> ga{ "test" };

//using namespace std;
/*
#include <experimental/coroutine>
#include <experimental/generator>

std::experimental::generator<long long> fibonacciGenerator(long long end)
{
    const char a[] = "test";
    long long curr = 0, next = 1;
    //*    {
        for (; next <= end;) {
            next = curr + next;
            curr = next - curr;
        }
    }
///*
    while (curr > 0) {
        curr = next - curr;
        next = next - curr;
        co_yield curr;
    }
}
*/
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
        std::less<uint32_t> {} //, [](uint32_t first, uint32_t last, uint32_t key) -> double { return (double(key) - first) / (last - first); }
    );
    //auto itr = FibonacciSearch(std::cbegin(arr), std::cend(arr), x, std::less<uint32_t>{});
    if (itr != std::cend(arr)) {
        return std::distance(std::cbegin(arr), itr);
    }
    return -1;
}

#include <random>

struct strr {

    strr()
    {
        std::cout << "strr const" << std::endl;
    }

    strr(const char *)
    {
        std::cout << "strr char * const" << std::endl;
    }

    strr(strr&&)
    {
        std::cout << "strr move const" << std::endl;
    }

    strr(const strr&)
    {
        std::cout << "strr copy const" << std::endl;
    }

    ~strr()
    {
        std::cout << "strr dest" << std::endl;
    }
};


struct A {
/*
    A(strr a) : str(std::move(a))
    {
        std::cout << "str copy const" << std::endl;
    }
/*/
    
    A(strr&& a)
        : str(std::move(a))
    {
        std::cout << "str move const" << std::endl;
    }
    
    A(const A& a) : str(a.str)
    {
        std::cout << "copy const" << std::endl;
    }
  //*/  
    A(A&& a)
        : str(std::move(a.str))
    {
        std::cout << "move const" << std::endl;
    }


    ~A()
    {
        std::cout << "dest" << std::endl;
    }


    strr str;
};


struct tttt {
    int32_t a;
    uint32_t b;
    void* p;

    uint64_t z;
};

int main()
{
    std::cout << sizeof(tttt) << std::endl;
        strr sss("asdfasdf");
        //A aaa("asfasdf");
        A aaa(std::move(sss));
    
    //strr rrrr = aaa.str;
        std::string aaaa;
        decltype(aaaa)::size_type s;

    return EXIT_SUCCESS;

    uint8_t a = 245;

    unsigned char z = a;

    char b = z;

    unsigned char y = b;

    if (auto i = std::max<unsigned char>(y, z); i != 10) {

    }

    std::cout << std::hex << int(a) << int(b) << int(z) << int(y) << std::endl;
    return EXIT_SUCCESS;

    //std::vector<uint32_t> vec = {10, 11, 11, 12, 18, 110, 111};
    //std::vector<uint32_t> vec = {10, 11, 12, 14, 16, 18, 110};
    std::vector<uint32_t> vec = { 10, 12, 14, 16, 18, 110, 111 };
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
