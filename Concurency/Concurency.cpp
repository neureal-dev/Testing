#include "Concurency.h"
#include <benchmark/benchmark.h>
#include <future>
#include <vector>
#include <omp.h>

size_t X(int a)
{
    size_t result = 0;

    #pragma omp simd
    for (size_t  i = 1; i < 1000000000; i++) {
        result += i;
        result /= i;
        result *= i;
    }
    return a&result;
}

int main()
{
    size_t result = 0;

    std::vector<std::future<size_t>> v(100);

    for (auto&& f : v) {
        f = std::async(std::launch::async, X, result++);
    }

    result = result * (result - 1) / 2;

    std::cout << result << std::endl;

    for (auto&& f : v) {
        result -= f.get();
        std::cout << result << std::endl;
    }

    std::cout << result << std::endl;

    return static_cast<int>(result);
}