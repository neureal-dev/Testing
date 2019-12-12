#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

namespace mtl {

template <typename T, int W, int N>
class PersistentVector {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>::type;
    using leaf_type = PersistentVector<T, W, N - 1>;
    using storage_type = std::array<std::unique_ptr<leaf_type>, 0x1 << W>;

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

    void clear() noexcept
    {
        for (auto& a : leafs_) {
            if (a) {
                a.reset();
            }
        }
    }

public:
    [[nodiscard]] static inline size_t getReferenceId(uint32_t id)
    {
        constexpr uint64_t xoffset = W * (N - 1);
        constexpr uint64_t xbits = (0x1 << W) - 1;

        return {(id >> xoffset) & xbits};
    }

    storage_type leafs_;
};

template <typename T, int W>
class PersistentVector<T, W, 0> {
public:
    using type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>::type;
    using storage_type = std::array<std::unique_ptr<type>, 0x1 << W>;

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
        leafs_[getReferenceId(id)] = std::make_unique<T>(std::forward<Args>(args)...);
        return *leafs_[getReferenceId(id)].get();
    }

    template <uint32_t>
    reference_type emplace_back(uint32_t id)
    {
        leafs_[getReferenceId(id)] = std::make_unique<T>(id);
        return *leafs_[getReferenceId(id)].get();
    }

private:
    [[nodiscard]] static inline size_t getReferenceId(uint32_t id)
    {
        constexpr size_t xbits = ((0x1 << W) - 1);

        return {(id >> W) & xbits};
    }

    storage_type leafs_;
};

} // namespace mtl
