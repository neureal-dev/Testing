// Performance.cpp : Defines the entry point for the application.
//

#include "Performance.h"
#include <benchmark/benchmark.h>
#include <random>
#include <vector>
#include <array>

struct A {
    A(uint32_t id, uint32_t ind)
        : id_(id), ind_(ind)
    {
    }
    uint32_t id_;
	uint32_t ind_;
	std::array<uint64_t, 8> ids;
};

/**
* Implementation of Interpolation search to optimize access path for number ordered arrarys
*/
template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator InterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

	difference_type count = std::distance(begin, end);

#ifdef NDEBUG
	RandomIterator last = std::prev(end);
#else
	RandomIterator last = count ? std::prev(end) : end;
#endif

	while (count > 0) {
		if (!comp(*begin, key)) {
			end = !comp(key, *begin) ? begin : end;
			break;
		}

		if (!comp(key, *last)) {
			end = !comp(*last, key) ? last : end;
			break;
		}

		difference_type probe = lerp(*begin, *last, key) * double(count - 1);

		if (comp(key, begin[probe])) {
			std::advance(last, probe - count);
			count = probe;
		} else if (comp(begin[probe], key)) {
			std::advance(begin, ++probe);
			count -= probe;
		} else {
			end = std::next(begin, probe);
			break;
		}
	}
	return end;
}

template <class ForwardIt, typename Value, typename Comparator>
inline ForwardIt BranchLessBinarySearch(ForwardIt begin, ForwardIt end, const Value& key, Comparator comp)
{
	using difference_type = std::iterator_traits<ForwardIt>::difference_type;
	difference_type size = std::distance(begin, end);

	while (difference_type half = size >> 1) {
		if (!comp(key, begin[half])) {
			std::advance(begin, half);
		}
		size = size - half;
	}
	return (begin != end && !comp(*begin, key) && !comp(key, *begin)) ? begin : end;
}

template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator HybridInterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

	difference_type count = std::distance(begin, end);

#ifdef NDEBUG
	RandomIterator last = std::prev(end);
#else
	RandomIterator last = count ? std::prev(end) : end;
#endif
	while (count > 0) {

		if (!comp(*begin, key)) {
			end = !comp(key, *begin) ? begin : end;
			break;
		}

		if (!comp(key, *last)) {
			end = !comp(*last, key) ? last : end;
			break;
		}

		difference_type probe = static_cast<difference_type>((count - 1) * lerp(*begin, *last, key));

		if (comp(key, begin[probe])) {
			//probe -= count;
			//std::advance(last, probe);
			//count += probe;

			//if (comp(key, begin[probe])) {
			//	probe -= count;
			//	std::advance(last, probe);
			//	count += probe;

			//	probe = count >> 1;
			//	if (!comp(key, begin[probe])) {
			//		count -= probe;
			//		std::advance(begin, probe);
			//	} else {
			//		probe -= count;
			//		std::advance(last, probe);
			//		count += probe;
			//	}
			/*} else if (comp(begin[probe], key)) {
				count -= ++probe;
				std::advance(begin, probe);

				probe = count >> 1;
				if (!comp(key, begin[probe])) {
					count -= probe;
					std::advance(begin, probe);
				} else {
					probe -= count;
					std::advance(last, probe);
					count += probe;
				}
			} else {
				end = begin;
				std::advance(end, probe);
				break;
			}
			*/

			/*/
			difference_type mid = (count - probe) >> 1;
			if (!comp(key, begin[mid])) {
				std::advance(begin, mid);
				std::advance(last, probe - count);
				count = (probe - mid);
			} else {
				std::advance(last, mid - count);
				count = mid;
			}
			/*/
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, probe - count);
				count = probe;
			}
			//*/
		} else if (comp(begin[probe], key)) {
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, probe - count);
				count = probe;
			}
			/*
			difference_type mid = probe + (count - probe) >> 1;
			if (!comp(key, begin[mid])) {
				std::advance(begin, mid);
				//std::advance(last, mid - count);
				count -= mid;
				//count = (mid - probe);
			} else {
				std::advance(begin, probe);
				std::advance(last, mid - count);
				count = mid - probe;
			}
			*/
			/*
			difference_type mid = probe + (count - probe) >> 1;
			if (!comp(begin[mid], key)) {
				std::advance(begin, probe);
				std::advance(last, mid - count);
				count = mid - probe;
			} else {

				std::advance(begin, mid);
				count -= mid;

			}
			*/
			/*
			probe = count >> 1;
			if (!comp(begin[probe], key)) {
				std::advance(last, probe - count);
				count = probe;
			} else {
				std::advance(begin, ++probe);
				count -= probe;
			}
			*/
			//probe = count >> 1;

/*
			if (!comp(key, begin[mid])) {

				std::advance(begin, probe);
				count -= probe ? probe : 1;
			} else {
				std::advance(last, probe - count);
				count = probe;
			}
			*/
		} else {
			end = std::next(begin, probe);
			break;
		}

	}
	return end;
}

template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator HHybridInterpolationSearch(RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;
	//std::cout << key << std::endl;
	difference_type count = std::distance(begin, end);

//#ifdef NDEBUG
	//RandomIterator last = std::prev(end);
//#else
	//RandomIterator last = count ? std::prev(end) : end;
//#endif
	/*
	if (!comp(*begin, key)) {
		return  !comp(key, *begin) ? begin : end;
	//	break;
	}

	if (!comp(key, *last)) {
		return  !comp(*last, key) ? last : end;
//		break;
	}
	*/
	//int s = 0;
	int cnd = 0;
	while (count > 0) {
		//s++;
		
		if (!comp(*begin, key)) {
			end = !comp(key, *begin) ? begin : end;
			break;
		}

		if (!comp(key, *(end - 1))) {
			end = !comp(*(end -1), key) ? (end - 1) : end;
			break;
		}
		
		difference_type probe = lerp(*begin, *(end - 1), key) * (count - 1);
		/*
		if (s > 18) {
			std::cout << key << "  " << s << " " << count << " " << probe << std::endl;
			std::cout << "    " << (*begin)->id_ << "-" << (*last)->id_ << std::endl;
			//std::cout << "!!!" << s << " " << key << (*end)->id_  <<std::endl;
		}
		*/
		/*
		if (comp(key, begin[probe])) {
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(last, --probe - count);
				count = probe;
			}
		} else if (comp(begin[probe], key)) {
			probe = count >> 1;
			if (!comp(begin[probe], key)) {
				std::advance(last, probe - count);
				count = probe;
			} else {
				std::advance(begin, ++probe);
				count -= probe;
			}
		} else {
			//end = std::next(begin, probe);
			std::advance(begin, probe);
			break;
		}
		*/
		if (comp(key, begin[probe])) {
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count -= probe;
			} else {
				std::advance(end, probe - count);
				count = probe;
			}
		} else if (comp(begin[probe], key)) {
			//*/
			probe = count >> 1;
			if (!comp(begin[probe], key)) {
				std::advance(end, ++probe - count);
				count = probe;
			} else {
				std::advance(begin, ++probe);
				count -= probe;
			}

			/*/
			difference_type mid = count >> 1;
			if (!comp(begin[mid], key)) {
				std::advance(begin, probe);
				std::advance(end, mid - count);
				count = mid - probe;
			} else {
				std::advance(begin, ++mid);
				count -= mid;
			}
			//*/
		} else {
			end = std::next(begin, probe);
			break;
		}
		
/*
		if (comp(key, begin[probe])) {
			std::advance(last, probe - count);
			count = probe;
			probe = count >> 1;
			if (!comp(key, begin[probe])) {
				std::advance(begin, probe);
				count = probe;
				cnd = 0;
			} else {
				last = begin + probe;
				//std::advance(last, -probe);
				count = probe + 1;
				//count -= count & 0x1;
				cnd = 1;
			}
			
		} else if (comp(begin[probe], key)) {
			std::advance(begin, ++probe);
			count -= probe;
			probe = count >> 1;
			if (!comp(begin[probe], key)) {
				last = begin + probe;
				count = probe + 1;
				//std::advance(last, probe - count);
				cnd = 2;
			} else {
				std::advance(begin, ++probe);
				count -= probe;
				cnd = 3;
			}
		} else {
			end = std::next(begin, probe);
			break;
		}
		*/
		if (count != ((*(end - 1))->ind_ - (*begin)->ind_) + 1) {
			std::cout << "error condition " << cnd <<  " " << count << " " << probe << " "<< (((*(end - 1))->ind_ - (*begin)->ind_) + 1 ) << std::endl;
		}
	}
	//std::cout << count << std::endl;
	//std::cout << (*end)->id_ << std::endl;
	return end;
}

/**
* An iterator class capable of navigating across the Fibonacci sequence using
* a user-specified integer type.
*/
template <typename Integer>
class FibonacciIterator : public std::iterator<std::bidirectional_iterator_tag, const Integer> {
public:

	explicit FibonacciIterator() : curr(0), next(1) { }

	explicit FibonacciIterator(Integer target) : curr(0), next(1)
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
bool operator==(const FibonacciIterator<Integer>& lhs, const FibonacciIterator<Integer>& rhs)
{
	return lhs.curr == rhs.curr;
}

/* Disequality implemented in terms of equality. */
template <typename Integer>
bool operator!=(const FibonacciIterator<Integer>& lhs, const FibonacciIterator<Integer>& rhs)
{
	return !(lhs == rhs);
}


/**
* Implementation of Fibonacci search uses a Fibonacci iterator to access
* the consecutive Fibonacci numbers
*/
template <typename RandomIterator, typename Value, typename Comparator>
RandomIterator FibonacciSearch(RandomIterator begin, RandomIterator end,
	const Value& value, Comparator comp)
{
	using difference_type = std::iterator_traits<RandomIterator>::difference_type;

	FibonacciIterator<difference_type> itr(std::distance(begin, end));
	uint32_t s = 0;
	while (*itr > 0) {
		//if (++s > 27) {
			//std::cout << value << "  " << s << " " << *itr << " " << std::endl;
			//std::cout << "    " << (*begin)->id_ << "-" << (*last)->id_ << std::endl;
			//std::cout << "!!!" << s << " " << key << (*end)->id_  <<std::endl;
		//}

		if (comp(value, begin[*itr - 1])) {
			--itr;
		} else if (comp(begin[*itr - 1], value)) {
			std::advance(begin, *itr);
			do {
				--itr;
			} while (*itr > std::distance(begin, end));
		} else {
			std::advance(begin, *itr - 1);
			return begin;
		}
	}
	return (begin != end && !(comp(*begin, value)) && !(comp(value, *begin))) ? begin : end;
}

struct Comp
{
	inline bool operator() ( const A* s, uint32_t i ) const { return s->id_ < i; }
	inline bool operator() ( uint32_t i, const A* s ) const { return i < s->id_; }
};

size_t HybridSearch(const std::vector<A*>& a, uint32_t x)
{
	int left = 0, right = a.size() - 1;
	int Inter, Mid;
	while (left < right) {
		//double s = double((x - a[left]->id_)) * (right - left) / (a[right]->id_ - a[left]->id_);
		Inter = left + double((x - a[left]->id_)) * (right - left) / (a[right]->id_ - a[left]->id_);
		if (Inter > right || Inter < left) {

		break;
	}
		if (x > a[Inter]->id_) {
			Mid = (Inter + right) / 2;
			if (x <= a[Mid]->id_) {
				left = Inter + 1;
				right = Mid;

			} else {
				left = Mid + 1;
			}
		} else if (x < a[Inter]->id_) {
			Mid = (Inter + left) / 2;
			if (x >= a[Mid]->id_) {
				left = Mid;
				right = Inter - 1;
			} else {
				right = Mid - 1;
			}
		} else {
			return Inter;
		}
	}
	std::cout << Inter << right << left << std::endl;
	return -1;
}

static void BM_SomeFunction(benchmark::State& state)
{
	std::cout << " saetase " << std::endl;
	/*/
	static std::vector<uint64_t> vecn ={1642, 1652, 1653, 1663, 1664, 1668, 1669, 1670, 1671, 1672, 1673, 1674, 1675, 1677,
		1681, 1689, 1694, 1695, 1706, 1935, 2012, 2401, 2525, 2644, 2711, 2930, 3368, 3376,
		3604, 4257, 8046, 87451, 88474, 88475, 88476, 88477, 88478, 88479, 88480, 88481, 88482,
		88483, 88484, 88485, 88486, 88487, 88488, 88489, 88490, 88491, 88492, 88493, 88494, 88495};
/*/
	/*
		[55]	88496	unsigned __int64
		[56]	88497	unsigned __int64
		[57]	88498	unsigned __int64
		[58]	88499	unsigned __int64
		[59]	88500	unsigned __int64
		[60]	88501	unsigned __int64
		[61]	88502	unsigned __int64
		[62]	88503	unsigned __int64
		[63]	88504	unsigned __int64
		[64]	91272	unsigned __int64
		[65]	92908	unsigned __int64
		[66]	92909	unsigned __int64
		[67]	92910	unsigned __int64
		[68]	92911	unsigned __int64
		[69]	92912	unsigned __int64
		[70]	92913	unsigned __int64
		[71]	92914	unsigned __int64
		[72]	92915	unsigned __int64
		[73]	92916	unsigned __int64
		[74]	92917	unsigned __int64
		[75]	92918	unsigned __int64
		[76]	92919	unsigned __int64
		[77]	92920	unsigned __int64
		[78]	92921	unsigned __int64
		[79]	92922	unsigned __int64
		[80]	92923	unsigned __int64
		[81]	92924	unsigned __int64
		[82]	92925	unsigned __int64
		[83]	92926	unsigned __int64
		[84]	92927	unsigned __int64
		[85]	92928	unsigned __int64
		[86]	92929	unsigned __int64
		[87]	92930	unsigned __int64
		[88]	96660	unsigned __int64
		[89]	98390	unsigned __int64
		[90]	98698	unsigned __int64
		[91]	98730	unsigned __int64
		[92]	98992	unsigned __int64
		[93]	99006	unsigned __int64
		[94]	99036	unsigned __int64
		[95]	99052	unsigned __int64
		[96]	99413	unsigned __int64
		[97]	99414	unsigned __int64
		[98]	99415	unsigned __int64
		[99]	99416	unsigned __int64
		[100]	99417	unsigned __int64
		[101]	99418	unsigned __int64
		[102]	99419	unsigned __int64
		[103]	99420	unsigned __int64
		[104]	99421	unsigned __int64
		[105]	99569	unsigned __int64
		[106]	99574	unsigned __int64
		[107]	99578	unsigned __int64
		[108]	99582	unsigned __int64
		[109]	99603	unsigned __int64
		[110]	99619	unsigned __int64
*/
	static std::vector<uint64_t> vecn;
	static std::vector<A *> vect;
	uint32_t i;
	if (vect.empty()) {
//		vect.reserve(1'000'000);
		
		while (std::cin >> i) {
			vecn.emplace_back(i);
		}
		std::sort(std::begin(vecn), std::end(vecn));
		std::cout << vecn.size() << std::endl;
		i = 0;
		for (auto e : vecn)
			vect.emplace_back(new A(e, i++));
		benchmark::DoNotOptimize(vect.data());
		benchmark::DoNotOptimize(vecn.data());
	}
	
	std::sort(std::begin(vect), std::end(vect), [](A* first, A* last) {return first->id_ < last->id_; });
	i = 0;
	for (auto _ : state) {
		//benchmark::DoNotOptimize(
		//auto a = FibonacciSearch(std::cbegin(vect), std::cend(vect), vecn[++i % vecn.size()], Comp{});
		//auto r = HybridSearch(vect, vecn[++i % vecn.size()]);
		//auto a = BranchLessBinarySearch(std::cbegin(vect), std::cend(vect), vecn[++i % vecn.size()],Comp{});
		auto a = HybridInterpolationSearch(std::cbegin(vect), std::cend(vect), vecn[++i % vecn.size()], Comp{}, [](A* first, A* last, uint32_t key) { return double(key - first->id_) / double(last->id_ - first->id_); });
		//if (r < vecn.size()) {
		if (a != std::end(vect)) {
			if ((*a)->id_ != vecn[i % vecn.size()])
			std::cout << "!!!err" << (*a)->id_ << " " << vecn[i % vecn.size()] <<  std::endl;
			//if (vect[r]->id_ != vecn[i % vecn.size()]) {
			//	std::cout << "error " << std::endl;
			//}
		}
		else {
			std::cout << vecn[i % vecn.size()] << "errror " << (a != std::end(vect)) << std::endl;
			//std::cout << "search error " << r << std::endl;
		}
		//);
    }
	std::cout << ".... saetase " << i % vect.size() << std::endl;
}
// Register the function as a benchmark
BENCHMARK(BM_SomeFunction);
// Run the benchmark
BENCHMARK_MAIN();