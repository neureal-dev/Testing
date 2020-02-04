#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

namespace mtl {

template <typename T, unsigned W, unsigned N>
class PersistentVector {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>::type;
    using leaf_type = PersistentVector<T, W, N - 1>;
    using storage_type = std::array<std::unique_ptr<leaf_type>, 0x1 << W>;

    class Iterator {
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = value_type;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type*;
        using reference = value_type&;

        explicit Iterator(PersistentVector& vect, size_t level) : root_(&vect), level_(level) {}

        explicit Iterator(PersistentVector& vect)
                : root_(&vect), level_(N - 1)
        {
            p[level_].first = root_->leafs_.begin();
            p[level_].second = root_->leafs_.end();
            for (; level_ > 0;) {
                for (; p[level_].first != p[level_].second; ++p[level_].first) {
                    if (*p[level_].first) {
                        leaf_type* l = p[level_].first->get();

                        //std::cout << "intr:" << level_ << ":" << l << std::endl;
                        --level_;
                        //                        p[level_].first = l->leafs_.begin();
                        //p[level_].second = l->leafs_.end();
                        p[level_].first = *reinterpret_cast<typename storage_type::iterator*>(&(l->leafs_.begin()));
                        p[level_].second = *reinterpret_cast<typename storage_type::iterator*>(&(l->leafs_.end()));
                        if (level_ == 0) {
                            for (; p[level_].first != p[level_].second; ++p[level_].first) {
                                if (*p[level_].first) {
                                    //std::cout << "intr:" << level_ << ":" << p[level_].first->get() << std::endl;
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }
                if (p[level_].first == p[level_].second) {
                    level_++;
                }

                if (level_ == N) {
                    break;
                }
            }
        }

        const reference operator*() const
        {
            //std::cout << p[0].first->get() << std::endl;
            return *((value_type*)(p[0].first->get()));
        }

        const pointer operator->() const { return &**this; }

        Iterator& operator++()
        {

            do {
                ++p[0].first;
            } while (p[0].first != p[0].second && !*p[0].first);

            if (p[0].first == p[0].second) {
                do {
                    level_++;
                    do {
                        ++p[level_].first;
                    } while (p[level_].first != p[level_].second && !*p[level_].first);
                } while (level_ < N && p[level_].first == p[level_].second);

                for (; level_ > 0;) {

                    for (; level_ < N && p[level_].first != p[level_].second; ++p[level_].first) {
                        if (*p[level_].first) {
                            leaf_type* l = p[level_].first->get();

                            //std::cout << "++:" << level_ << ":" << l << std::endl;
                            --level_;
                            //                        p[level_].first = l->leafs_.begin();
                            //p[level_].second = l->leafs_.end();
                            p[level_].first = *reinterpret_cast<typename storage_type::iterator*>(&(l->leafs_.begin()));
                            p[level_].second = *reinterpret_cast<typename storage_type::iterator*>(&(l->leafs_.end()));
                            if (level_ == 0) {
                                for (; p[level_].first != p[level_].second; ++p[level_].first) {
                                    if (*p[level_].first) {
                              //          std::cout << "intr:" << level_ << ":" << p[level_].first->get() << std::endl;
                                        break;
                                    }
                                }
                            }
                            break;
                        }
                    }
                    if (p[level_].first == p[level_].second) {
                        level_++;
                    }

                    if (level_ == N) {
                        break;
                    }
                }
            }
            return *this;
        }

        const Iterator operator++(int)
        {
            Iterator result = *this;
            ++*this;
            return result;
        }

        Iterator& operator--()
        {
            return *this;
        }

        const Iterator operator--(int)
        {
            Iterator result = *this;
            --*this;
            return result;
        }

        bool operator==(Iterator& rit)
        {
            //std::cout << root_ << rit.root_ << level_ << rit.level_ << std::endl;
            return root_ == rit.root_ && level_ == rit.level_;
        }

        bool operator!=(Iterator& rit)
        {
            return !operator==(rit);
        }

    private:
        PersistentVector* root_;
        size_t level_;
        std::array<std::pair<typename storage_type::iterator, typename storage_type::iterator>, N> p;
    };

    using iterator = Iterator;

    void addNode(uint32_t id, pointer_type element)
    {
        size_t ref_id = getReferenceId(id);
        if (!leafs_[ref_id]) {
            leafs_[ref_id] = std::make_unique<value_type>();
        }
        leafs_[ref_id]->addNode(id, element);
    }

    template <class... Args>
    reference_type emplace_back(uint32_t id, Args&&... args)
    {
        size_t ref_id = getReferenceId(id);
        if (!leafs_[ref_id]) {
            leafs_[ref_id] = std::make_unique<leaf_type>();
            //std::cout << "crt:" << N - 1 << ":" << ref_id << ":" << leafs_[ref_id].get() << std::endl;
        }
        return leafs_[ref_id]->emplace_back(id, std::forward<Args>(args)...);
    }

    void push_back(value_type&& value)
    {
        emplace_back(value.id_, std::move(value));
    }

    void push_back(const value_type& value)
    {
        emplace_back(value.id_, value);
    }

    [[nodiscard]] pointer_type getNode(uint32_t id) const
    {
        auto& node = leafs_[getReferenceId(id)];
        return node ? node->getNode(id) : nullptr;
    }

    iterator begin()
    {
        return Iterator(*this);
    }

    iterator end()
    {
        return Iterator(*this, N);
    }

    void clear() noexcept
    {
        for (auto& a : leafs_) {
            if (a) {
                a.reset();
            }
        }
    }

    void shrink_to_fit() {}

public:
    [[nodiscard]] static constexpr size_t getReferenceId(uint32_t id)
    {
        constexpr uint64_t xoffset = W * (N - 1);
        constexpr uint64_t xbits = (0x1 << W) - 1;

        return {(id >> xoffset) & xbits};
    }

    storage_type leafs_;
};

template <typename T, unsigned W>
class PersistentVector<T, W, 1> {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<value_type>::type;
    using reference_type = typename std::reference_wrapper<value_type>::type;
    using storage_type = std::array<std::unique_ptr<value_type>, 0x1 << W>;
    //using storage_type = std::array<pointer_type, 0x1 << W>;
    /*
    class Iterator {
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = typename storage_type::iterator::value_type;
        using difference_type = void;
        using pointer = typename storage_type::iterator::pointer;
        using reference = typename storage_type::iterator::reference;

        explicit Iterator(typename storage_type::iterator begin, typename storage_type::iterator end)
                : root_itr_(begin), root_end_(end)
        {
            for (; root_itr_ != root_end_; ++root_itr_) {
                if (*root_itr_) {
                    break;
                }
            }
        }

        const reference operator*() const { return *root_itr_; }

        const pointer operator->() const { return &**this; }

        Iterator& operator++()
        {
            do {
                ++root_itr_;
            } while (root_itr_ != root_end_ && !*root_itr_);
            return *this;
        }

        const Iterator
        operator++(int)
        {
            Iterator result = *this;
            ++*this;
            return result;
        }

        Iterator& operator--()
        {
            do {
                --root_itr_;
            } while (!*root_itr_);
            return *this;
        }

        const Iterator operator--(int)
        {
            Iterator result = *this;
            --*this;
            return result;
        }

        bool operator==(Iterator& rit)
        {
            return root_itr_ == rit.root_itr_;
        }

        bool operator!=(Iterator& rit)
        {
            return root_itr_ != rit.root_itr_;
        }

    private:
        typename storage_type::iterator root_itr_;
        typename storage_type::iterator root_end_;
    };
    */
    //using iterator = Iterator;

    void addNode(uint32_t id, pointer_type element)
    {
        leafs_[getReferenceId(id)] = element;
    }

    [[nodiscard]] pointer_type getNode(uint32_t id) const
    {
        return leafs_[getReferenceId(id)].get();
    }

    template <class... Args>
    reference_type emplace_back(uint32_t id, Args&&... args)
    {
        leafs_[getReferenceId(id)] = std::make_unique<T>(std::forward<Args>(args)...); //new T(std::forward<Args>(args)...); //
        //std::cout << "create:" << id << ":" << leafs_[getReferenceId(id)].get() << std::endl;
        return *leafs_[getReferenceId(id)].get();
    }

    template <uint32_t>
    reference_type emplace_back(uint32_t id)
    {
        leafs_[getReferenceId(id)] = std::make_unique<T>(id);
        return *leafs_[getReferenceId(id)].get();
    }

    void push_back(value_type&& value)
    {
        emplace_back(value.id_, std::move(value));
    }

    void push_back(const value_type& value)
    {
        emplace_back(value.id_, value);
    }
    /*
    iterator begin()
    {
        //        for (auto itr = leafs_.begin(); itr != leafs_.end(); ++itr) {
        //            if (itr->get()) {
        //                return Iterator(itr, leafs_.end());
        //            }
        //        }
        return Iterator(leafs_.begin(), leafs_.end());
    }

    iterator end()
    {
        return Iterator(leafs_.end(), leafs_.end());
    }
    */
    void clear() noexcept
    {
        for (auto& a : leafs_) {
            if (a) {
                // a.reset();
            }
        }
    }

    void shrink_to_fit() {}

    //private:
    [[nodiscard]] static inline size_t getReferenceId(uint32_t id)
    {
        constexpr size_t xbits = ((0x1u << W) - 1);
        //return {(id >> W) & xbits};
        return id & xbits;
    }

    storage_type leafs_;
}; // namespace mtl

} // namespace mtl
