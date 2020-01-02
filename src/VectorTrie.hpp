#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <vector>

namespace mtl {

template <typename T, unsigned W, unsigned N>
class PersistentVector {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>;
    using leaf_type = PersistentVector<T, W, N - 1>;
    using storage_type = std::array<std::unique_ptr<leaf_type>, 0x1 << W>;

    class Iterator {
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = value_type;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type*;
        using reference = value_type&;

        explicit Iterator() : root_(0)
        {
        }

        explicit Iterator(PersistentVector& vect, bool)
                : root_(&vect), root_itr_(root_->leafs_.end()), root_ind_(0x01 << W) {}

        explicit Iterator(PersistentVector& vect)
                : root_(&vect), root_itr_(root_->leafs_.begin()), root_ind_(0)
        {
            for (; root_ind_ != 0x1 << W; ++root_itr_, ++root_ind_) {
                if (*root_itr_) {
                    leaf_itr_ = typename leaf_type::iterator(*root_itr_->get());
                    break;
                }
            }
        }

        const reference operator*() const
        {
            return *leaf_itr_;
        }

        const pointer operator->() const { return &**this; }

        Iterator& operator++()
        {
            if (!leaf_itr_.adjust()) {
                for (++root_ind_, ++root_itr_; root_ind_ != 0x1 << W; ++root_ind_, ++root_itr_) {
                    if (*root_itr_) {
                        leaf_itr_ = typename leaf_type::iterator(*root_itr_->get());
                    }
                }
            }
            return *this;
        }

        bool adjust()
        {
            if (!leaf_itr_.adjust()) {
                for (++root_ind_, ++root_itr_; root_ind_ != 0x1 << W; ++root_ind_, ++root_itr_) {
                    if (*root_itr_) {
                        leaf_itr_ = typename leaf_type::iterator(*root_itr_->get());
                        return true;
                    }
                }
                return false;
            }
            return true;
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
            return root_ == rit.root_ && root_itr_ == rit.root_itr_;
        }

        bool operator!=(Iterator& rit)
        {
            return !operator==(rit);
        }

    private:
        PersistentVector* root_;
        typename storage_type::iterator root_itr_;
        typename leaf_type::iterator leaf_itr_;
        size_t root_ind_;
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
        return Iterator(*this, false);
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
    [[nodiscard]] static inline size_t getReferenceId(uint32_t id)
    {
        constexpr uint64_t xoffset = W * (N - 1);
        constexpr uint64_t xbits = (0x1 << W) - 1;

        return {(id >> xoffset) & xbits};
    }

    storage_type leafs_;
};

template <typename T, unsigned W>
class PersistentVector<T, W, 0> {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<value_type>::type;
    using reference_type = typename std::reference_wrapper<value_type>::type;
    using storage_type = std::array<std::unique_ptr<value_type>, 0x1 << W>;
    //using storage_type = std::array<pointer_type, 0x1 << W>;

    class Iterator {
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = typename storage_type::iterator::value_type;
        using difference_type = void;
        using pointer = typename storage_type::iterator::pointer;
        using reference = typename storage_type::iterator::reference;

        explicit Iterator()
                : root_(0)
        {
        }

        explicit Iterator(typename storage_type::iterator begin, typename storage_type::iterator end)
                : root_itr_(begin), root_emd_(end)
        {
 
        }

        const reference operator*() const { return *root_itr_; }

        const pointer operator->() const { return &**this; }
        
        Iterator& operator++()
        {
            ++root_itr_;
            
            //for (++root_itr_; root_itr_ != root_emd_; ++root_itr_) {
            //    if (*root_itr_) {
            //        break;
            //    }
            //}
            
            return *this;
        }
        
        bool adjust()
        {
            for (; root_ind_ < 0x01 << W;) {
                if (root_->leafs_[++root_ind_]) {
                    return true;
                }
            }
            return false;
        }
        //}
        
    const Iterator
    operator++(int)
    {
        Iterator result = *this;
        ++*this;
        return result;
    }

    Iterator& operator--()
    {
        --root_itr_;
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
        //PersistentVector* root_;
        typename storage_type::iterator root_itr_;
        typename storage_type::iterator root_emd_;
    };

    using iterator = Iterator;

    void addNode(uint32_t id, pointer_type element)
    {
        leafs_[getReferenceId(id)] = element;
    }

    [[nodiscard]] pointer_type getNode(uint32_t id) const
    {
        return leafs_[getReferenceId(id)];
        //.get();
    }

    template <class... Args>
    reference_type emplace_back(uint32_t id, Args&&... args)
    {
        leafs_[getReferenceId(id)] = std::make_unique<T>(std::forward<Args>(args)...);//new T(std::forward<Args>(args)...); //
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
