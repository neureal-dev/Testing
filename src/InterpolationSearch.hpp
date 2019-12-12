#pragma once

namespace mtl {

template <typename RandomIterator, typename Value, typename Comparator, typename Converter>
RandomIterator HybridInterpolationSearch(
        RandomIterator begin, RandomIterator end, Value key, Comparator comp, Converter lerp)
{
    using difference_type = typename std::iterator_traits<RandomIterator>::difference_type;

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

        auto probe = static_cast<difference_type>((count - 1) * lerp(*begin, *last, key));

        if (comp(key, begin[probe])) {
            auto mid = probe >> 1;
            if (!comp(key, begin[mid])) {
                std::advance(last, probe - count);
                count = probe - mid;
                std::advance(begin, mid);
            } else {
                std::advance(last, mid - count);
                count = mid;
            }
        } else if (comp(begin[probe], key)) {
            auto mid = (probe + count) >> 1;
            if (!comp(key, begin[mid])) {
                std::advance(begin, mid);
                count -= mid;
            } else {
                std::advance(last, mid - count);
                count = mid - probe;
                std::advance(begin, probe);
            }
        } else {
            end = begin;
            std::advance(end, probe);
            break;
        }
    }
    return end;
}

} // namespace mtl
