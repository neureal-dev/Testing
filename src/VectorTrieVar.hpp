#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace mtl {

template <typename T, size_t W, size_t N>
class PersistentVectorVar {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>::type;

    //class Leafs {
    //public:
    // /   std::array<, 0x1 << W> nodes;
    //}/;

    //using Node = leaf_type;

    class Nodes {
    public:
        std::array<std::variant<std::unique_ptr<Nodes>, std::unique_ptr<value_type>>, 0x1 << W> nodes;
    };

    using leaf_type = Nodes;
    using storage_type = leaf_type;

    template <class... Args>
    reference_type emplace_back(uint32_t id, Args&&... args)
    {
        leaf_type* leaf = &leafs_;
        for (size_t l = N - 1; l > 0; --l) {
            size_t ref_id = getReferenceId(id, l);
            //std::cout << "for" << l << " " << ref_id << " " << id << std::endl;
            //if (!std::holds_alternative<std::unique_ptr<Nodes>>(leaf->nodes[ref_id])) {
            auto& n = std::get<std::unique_ptr<Nodes>>(leaf->nodes[ref_id]);
            if (!n) {
                //std::cout  << "create"<< l << " " << ref_id << " " << id << std::endl;
                leaf->nodes[ref_id] = std::make_unique<Nodes>();
            }
            leaf = std::get<std::unique_ptr<Nodes>>(leaf->nodes[ref_id]).get();
        }
        constexpr size_t xbits = ((0x1u << W) - 1);

        leaf->nodes[id & xbits] = std::make_unique<T>(std::forward<Args>(args)...);
        return *std::get<std::unique_ptr<value_type>>(leaf->nodes[id & xbits]).get();
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
        const leaf_type* node = &leafs_;
        for (size_t l = N - 1; node && l > 0; --l) {
            node = std::get<std::unique_ptr<Nodes>>(node->nodes[getReferenceId(id, l)]).get();
            //node = node->nodes[getReferenceId(id, l)].get();
        }
        constexpr size_t xbits = ((0x1u << W) - 1);

        return node ? std::get<std::unique_ptr<value_type>>(node->nodes[id & xbits]).get() : nullptr;
    }

    void clear() noexcept
    {
        //for (auto& a : leafs_) {
        //    if (a) {
        //        a.reset();
        //    }
        //}
    }

    void shrink_to_fit() {}

public:
    [[nodiscard]] static inline size_t getReferenceId(uint32_t id, size_t n)
    {
        uint64_t xoffset = W * n;
        constexpr uint64_t xbits = (0x1 << W) - 1;

        return (id >> xoffset) & xbits;
    }

    storage_type leafs_;
};

} // namespace mtl
