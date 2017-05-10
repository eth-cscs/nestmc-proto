#pragma once

#pragma once
#pragma once

#include <common_types.hpp>
#include <util/debug.hpp>
#include <vector>
#include <iostream>

namespace nest {
namespace mc {


class cell_interface {
public:
    virtual ~cell_interface() = default;

    /// Return the kind of cell, used for grouping into cell_groups
    virtual cell_kind get_cell_kind() const = 0;

    ///// Collect all spikes until tfinal.
    //// updates the internal time state to tfinal as a side effect
    //virtual std::vector<time_type> spikes_until(time_type tfinal) = 0;

    ///// reset internal state;
    //virtual void reset() = 0;

};


struct cell_description {
    template <typename T, typename = typename std::enable_if<!std::is_reference<T>::value>::type>
    cell_description(T&& concrete_cell) :
        cellptr(new T(std::move(concrete_cell))) {
            kind = concrete_cell.get_cell_kind();
    }

    template <typename T>
    T& as() {
        auto ptr = dynamic_cast<T*>(cellptr.get());

        if (!ptr) {
            std::cout << "Trying to cast cell_description a the wrong type!!";
            throw 1;
        }

        return *ptr;
    }

    cell_kind kind;

private:
    std::unique_ptr<cell_interface> cellptr;
};

} // namespace mc
} // namespace nest
