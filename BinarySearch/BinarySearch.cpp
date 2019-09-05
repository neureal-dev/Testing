// BinarySearch.cpp : Defines the entry point for the application.
//

#include "BinarySearch.h"
#include <vector>
#include <algorithm>

using namespace std;

size_t bsearch(const std::vector<uint32_t>& resource, uint32_t key)
{
    size_t low, high, mid;

    low = 0;
    high = std::min<size_t>(key, resource.size());
    mid = low + (high - low) / 2;

    while (low < high) {
        auto val = resource[mid];
        if (key < val) {
            high = mid;
        } else {
            if (key == val) {
                break;
            }
            low = mid + 1;
        }
        mid = low + (high - low) / 2;
    }
    return resource[mid] == key ? mid : resource.size();
}


int main()
{
    std::vector<uint32_t> vec = {};
    std::cout << bsearch(vec, 0) << std::endl;
    std::cout << bsearch(vec, 1) << std::endl;
    std::cout << bsearch(vec, 2) << std::endl;
    std::cout << bsearch(vec, 4) << std::endl;
    std::cout << bsearch(vec, 6) << std::endl;
    std::cout << bsearch(vec, 8) << std::endl;
    std::cout << bsearch(vec, 10) << std::endl;
	cout << "Hello CMake." << endl;
	return 0;
}
