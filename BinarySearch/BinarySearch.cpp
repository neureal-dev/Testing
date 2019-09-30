// BinarySearch.cpp : Defines the entry point for the application.
//

#include "BinarySearch.h"
#include <vector>
#include <algorithm>

using namespace std;

namespace qb {

template <class ForwardIt, class T>
inline ForwardIt branch_less_binary_find_n(ForwardIt first, size_t size, const T& key)
{
	using difference_type = std::iterator_traits<ForwardIt>::difference_type;

	while (difference_type half = size >> 1)
	{
		ForwardIt middle = first;
		std::advance(middle, half);
		if (!((*middle) > key)) {
			std::advance(first, half);
		}
		size = size - half;
	}
	return first;
}


template< class ForwardIt, class T>
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
	half             = std::distance(first, last) >> 1;
	ForwardIt middle = first, ind = last;
	std::advance(middle, half);
	while (first < last)
	{
		auto val = (*middle);

		if (key > val)
		{
			first = middle;
			std::advance(first, 1);
		}
		else
		{
			if (key == val)
			{
				ind = middle;
				break;
			}
			last = middle;
		}
		half   = std::distance(first, last) >> 1;
		middle = first;
		std::advance(middle, half);
	}
	return ind;
}

template<class ForwardIt, class T>
inline ForwardIt equal_range(ForwardIt first, ForwardIt last, const T& key)
{
	using difference_type = std::iterator_traits<ForwardIt>::difference_type;
	difference_type size = std::distance(first, last);
	return (size == 0) ? last : size < 0x0
		? branch_less_binary_find_n(first, size, key)
		: branch_full_binary_find_n(first, last, size, key);
}

};


size_t bsearch(const std::vector<uint32_t>& arr, uint32_t x)
{
//	auto ind = qb::equal_range(std::cbegin(vec), std::cend(vec), key);
//	return ind != std::cend(vec) ? std::distance(std::cbegin(vec), ind) : -1;
//	uint64_t fibMonaccianSearch(const std::vector<const SBrecord*>& arr, uint32_t x)
//	{ 
		int64_t n = arr.size()-1;
		/* Initialize fibonacci numbers */
		int64_t fibMMm2 = 0;   // (m-2)'th Fibonacci No. 
		int64_t fibMMm1 = 1;   // (m-1)'th Fibonacci No. 
		int64_t fibM = fibMMm2 + fibMMm1;  // m'th Fibonacci 

										   /* fibM is going to store the smallest Fibonacci 
										   Number greater than or equal to n */
		while (fibM <= n) 
		{ 
			fibMMm2 = fibMMm1; 
			fibMMm1 = fibM; 
			fibM  = fibMMm2 + fibMMm1; 
		} 

		// Marks the eliminated range from front 
		int64_t offset = -1; 

		/* while there are elements to be inspected. Note that 
		we compare arr[fibMm2] with x. When fibM becomes 1, 
		fibMm2 becomes 0 */
		while (fibM > 1) 
		{ 
			// Check if fibMm2 is a valid location 
			int64_t i = std::min(offset+fibMMm2, n-1); 

			/* If x is greater than the value at index fibMm2, 
			cut the subarray array from offset to i */
			if (arr[i] < x) 
			{ 
				fibM  = fibMMm1; 
				fibMMm1 = fibMMm2; 
				fibMMm2 = fibM - fibMMm1; 
				offset = i; 
			} 

			/* If x is greater than the value at index fibMm2, 
			cut the subarray after i+1  */
			else if (arr[i] > x) 
			{ 
				fibM  = fibMMm2; 
				fibMMm1 = fibMMm1 - fibMMm2; 
				fibMMm2 = fibM - fibMMm1; 
			} 

			/* element found. return index */
			else return i; 
		} 

		/* comparing the last element with x */
		//std::cout << "o" << offset << std::endl;
		if(fibMMm1 && arr[offset+1]==x)return offset+1; 

		/*element not found. return -1 */
		return -1; 
//	} 
}

int main()
{
    std::vector<uint32_t> vec = {10, 11, 11, 12, 18, 110};
    std::cout << bsearch(vec, 1) << std::endl;
    std::cout << bsearch(vec, 10) << std::endl;
    std::cout << bsearch(vec, 11) << std::endl;
    std::cout << bsearch(vec, 12) << std::endl;
    std::cout << bsearch(vec, 14) << std::endl;
    std::cout << bsearch(vec, 16) << std::endl;
    std::cout << bsearch(vec, 18) << std::endl;
    std::cout << bsearch(vec, 110) << std::endl;
    std::cout << bsearch(vec, 120) << std::endl;
/*
	std::cout << (btsearch(vec.begin(), vec.end(), 1) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 10) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 11) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 12) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 13) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 14) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 16) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 18) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 110) == vec.end()) << std::endl;
	std::cout << (btsearch(vec.begin(), vec.end(), 120) == vec.end()) << std::endl;
*/

	cout << "Hello CMake." << endl;
	return 0;
}
