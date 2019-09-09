// BinarySearch.cpp : Defines the entry point for the application.
//

#include "BinarySearch.h"
#include <vector>
#include <algorithm>

using namespace std;

size_t bsearch(const std::vector<uint32_t>& resource, uint32_t key)
{
    size_t low = 0;
    size_t high = resource.size();
    if (high > 1) {
        high = std ::min<size_t>(1ull + key - resource[low], high - 1);
        //high = std ::min<size_t>(key - resource[low], high);
    }
    //size_t high = resource.size() - 1;
    size_t ind = -1;
    size_t mid = high >> 1;

    while (low < high) {
        //auto mid = low + ((high - low) >> 1);
        auto val = resource[mid];

        if (key < val) {
            high = mid;
        } else {
            if (key == val) {
                ind = mid;
                break;
            }
            low = mid + 1;
        }

        mid = low + ((high - low) >> 1);
    }
    return ind;
    //return resource.size() && resource[mid]->id_ == key ? mid : resource.size();
}


int main()
{
    std::vector<uint32_t> vec = {10, 11, 12, 18, 110};
    std::cout << bsearch(vec, 1) << std::endl;
    std::cout << bsearch(vec, 10) << std::endl;
    std::cout << bsearch(vec, 11) << std::endl;
    std::cout << bsearch(vec, 12) << std::endl;
    std::cout << bsearch(vec, 14) << std::endl;
    std::cout << bsearch(vec, 16) << std::endl;
    std::cout << bsearch(vec, 18) << std::endl;
    std::cout << bsearch(vec, 110) << std::endl;
    std::cout << bsearch(vec, 120) << std::endl;
	cout << "Hello CMake." << endl;
	return 0;
}
