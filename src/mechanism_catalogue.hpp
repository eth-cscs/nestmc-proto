#pragma once

#include <map>
#include <stdexcept>
#include <string>

#include <mechanism.hpp>
#include <mechanisms/hh.hpp>
#include <mechanisms/pas.hpp>
#include <mechanisms/expsyn.hpp>
#include <mechanisms/exp2syn.hpp>

#ifdef WITH_CUDA
    #include <mechanisms/gpu/hh.hpp>
    #include <mechanisms/gpu/pas.hpp>
    #include <mechanisms/gpu/expsyn.hpp>
    #include <mechanisms/gpu/exp2syn.hpp>
#endif

namespace nest {
namespace mc {
namespace mechanisms {

enum class targetKind {host, gpu};

template <typename MemoryTraits>
class catalogue {
public:

    using memory_traits = MemoryTraits;
    using view        = typename memory_traits::view;
    using const_iview = typename memory_traits::const_iview;
    using mechanism_type = mechanism<memory_traits>;
    using mechanism_ptr_type = mechanism_ptr<memory_traits>;

    template <typename Indices>
    static mechanism_ptr_type make(
        const std::string& name, view vec_v, view vec_i,
        const Indices& node_indices)
    {
        auto entry = mech_map_.find(name);
        if (entry==mech_map_.end()) {
            throw std::out_of_range("no mechanism in database : " + name);
        }

        return entry->second(vec_v, vec_i,  memory::make_const_view(node_indices));
    }

    static bool has(const std::string& name) {
        return mech_map_.count(name)>0;
    }

private:
    using maker_type = mechanism_ptr_type (*)(view, view, const_iview);
    static const std::map<std::string, maker_type> mech_map_;

    template <template <typename> class Mech>
    static mechanism_ptr_type maker(view vec_v, view vec_i, const_iview node_indices) {
        return make_mechanism<Mech<memory_traits>>(vec_v, vec_i, node_indices);
    }
};

template <typename MemoryTraits>
const std::map<std::string, typename catalogue<MemoryTraits>::maker_type> catalogue<MemoryTraits>::mech_map_ = {
    { "pas",     maker<pas::mechanism_pas> },
    { "hh",      maker<hh::mechanism_hh> },
    { "expsyn",  maker<expsyn::mechanism_expsyn> },
    { "exp2syn", maker<exp2syn::mechanism_exp2syn> }
};

#ifdef WITH_CUDA
namespace gpu {
    template <typename T, typename I>
    struct catalogue {
        using view_type  = typename gpu::mechanism<T, I>::view_type;
        using index_view = typename gpu::mechanism<T, I>::index_view;
        using index_type = typename gpu::mechanism<T, I>::index_type;
        using const_index_view = typename memory::HostVector<I>::const_view_type;

        template <typename Indices>
        static mechanism_ptr<T, I> make(
            const std::string& name,
            view_type vec_v,
            view_type vec_i,
            const Indices& node_indices)
        {
            auto entry = mech_map.find(name);
            if (entry==mech_map.end()) {
                throw std::out_of_range("no such mechanism");
            }

            return entry->second(vec_v, vec_i, memory::make_const_view(node_indices));
        }

        static bool has(const std::string& name) {
            return mech_map.count(name)>0;
        }

    private:
        using maker_type = gpu::mechanism_ptr<T, I> (*)(view_type, view_type, const_index_view);
        static const std::map<std::string, maker_type> mech_map;

        template <template <typename, typename> class Mech>
        static gpu::mechanism_ptr<T, I> maker(
            view_type vec_v,
            view_type vec_i,
            const_index_view node_indices)
        {
            return gpu::make_mechanism<Mech<T, I>>(vec_v, vec_i, node_indices);
        }
    };

    template <typename T, typename I>
    const std::map<std::string, typename catalogue<T, I>::maker_type> catalogue<T, I>::mech_map = {
        { "pas",     maker<pas::mechanism_pas> },
        { "hh",      maker<hh::mechanism_hh> },
        { "expsyn",  maker<expsyn::mechanism_expsyn> },
        { "exp2syn", maker<exp2syn::mechanism_exp2syn> }
    };
}
#endif

} // namespace mechanisms
} // namespace mc
} // namespace nest
