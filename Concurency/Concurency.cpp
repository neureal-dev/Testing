// Performance.cpp : Defines the entry point for the application.
//

#include "Concurency.h"
#include <benchmark/benchmark.h>
#include <future>
#include <vector>

size_t X(int a)
{
    //std::cout << "thr: " << std::this_thread::get_id() << std::endl;
    size_t result = 0;
    for (size_t  i = 1; i < 1000000000; i++) {
        result += i;
        result /= i;
        result *= i;
    }
    return a&result;
}


//static void BM_SomeFunction(benchmark::State& state)
//{
//    std::vector<std::future<size_t>> v(100);
//    for (auto&& f : v) {
//        f = std::async(X, 43ul);
//    }
//}
/*
BENCHMARK(BM_SomeFunction)->Threads(1);

BENCHMARK_MAIN();
*/

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