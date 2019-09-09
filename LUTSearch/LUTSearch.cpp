// LUTSearch.cpp : Defines the entry point for the application.
//

#include "LUTSearch.h"
#include <thread>
#include <vector>
#include <algorithm>


inline size_t bbsearch(const std::vector<uint32_t>& resource, uint32_t const key, size_t low, size_t high)
{
    size_t ind = -1;
    size_t mid = high >> 1;

    while (low < high) {
        auto val = resource[mid];
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

//size_t bsearch(const std::vector<uint32_t>& data, uint32_t const key) { return bbsearch(data, key, 0, data.size()); }

size_t fallback_bsearch(const std::vector<uint32_t>& data, uint32_t key)
{
    size_t low = 0;
    size_t high = data.size();
    size_t mid = high >> 1;
    // size_t ind = -1;
    size_t size = high;
    while (size > 0x8F) {
        // auto mid = (low + high) >> 1;
        auto val = data[mid];
        if (key > val) {
            low = mid + 1;
        } else {
            if (key == val) {
                return mid;
                //low = high = mid;
                //break;
            }
            high = mid;
        }
        mid = (low + high) >> 1;
        size = high - low;
    }
    return bbsearch(data, key, low, high);
}



size_t b2search(const std::vector<uint32_t>& records, uint32_t rid)
{
    size_t size = records.size();
    size_t low = 0;
    size_t middle;
    while (size_t half = size >> 1) {
        middle = low + half;
        auto val = records[middle];
        low = (val <= rid) ? middle : low;
        size = size - half;
    }
    return (size && records[low] == rid) ? low : -1;
}
size_t b2search2(const std::vector<uint32_t>& resource, uint32_t key)
{
    size_t low = 0;
    // size_t high = std::min<size_t>(key, (resource.size() > 1 ? resource.size() - 1 : resource.size()));
    size_t high = resource.size();
    //if (high > 2) {
        //high -= 1;
    //} else {
        //return -1;
    //}
    size_t ind = -1;
    //size_t mid = high >> 1;

    while (low < high) {
        auto mid = (low + high) >> 1;
        auto val = resource[mid];
        if (key > val) {
            low = mid + 1;
        } else {
            if (key == val) {
                ind = mid;
                break;
            }
            high = mid - 1;
        }
        mid = (low + high) >> 1;
    }
    return ind;
}



int main()
{
    //std::vector<uint32_t> rec = {5};
    
	std::vector<uint32_t> rec(0xFFFF, 0);
    long i = 0;
    for (auto& a : rec) {
        a = ++i;
    }
    i += 100;
	for (; i > 0; i--) {
        auto j = fallback_bsearch(rec, i);
        //auto j = bbsearch(rec, i, 0, rec.size());
        if (i != j + 1) {
            std::cout << i << " " << j << std::endl;
        }
    }
      

    std::cout << fallback_bsearch(rec, 0) << std::endl;
    std::cout << fallback_bsearch(rec, 3) << std::endl;
    std::cout << fallback_bsearch(rec, 5) << std::endl;
    std::cout << fallback_bsearch(rec, 6) << std::endl;
    std::cout << fallback_bsearch(rec, 10) << std::endl;
    std::cout << fallback_bsearch(rec, 11) << std::endl;
    std::cout << fallback_bsearch(rec, 16) << std::endl;
    std::cin.ignore();
    return 0;
}
