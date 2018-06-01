#pragma once

#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

#include <cassert>

#include <mpi.h>

#include <algorithms.hpp>
#include <communication/gathered_vector.hpp>
#include <util/debug.hpp>
#include <profiling/profiler.hpp>


namespace arb {
namespace mpi {

// prototypes
void init(int *argc, char ***argv);
void finalize();
int rank(MPI_Comm);
int size(MPI_Comm);
void barrier(MPI_Comm);

void handle_mpi_error(const char* msg, int code);

// Exception class to be thrown when MPI API calls return a error code other
// than MPI_SUCCESS.
class mpi_error: public std::exception {
public:
    mpi_error(const char* msg, int code);
    const char* what() const throw() override;
    int error_code() const;

private:
    std::string message_;
    int error_code_;
};

struct scoped_guard {
    scoped_guard(int *argc, char ***argv);
    ~scoped_guard();
};

// Type traits for automatically setting MPI_Datatype information for C++ types.
template <typename T>
struct mpi_traits {
    constexpr static size_t count() {
        return sizeof(T);
    }
    constexpr static MPI_Datatype mpi_type() {
        return MPI_CHAR;
    }
    constexpr static bool is_mpi_native_type() {
        return false;
    }
};

#define MAKE_TRAITS(T,M)     \
template <>                 \
struct mpi_traits<T> {  \
    constexpr static size_t count()            { return 1; } \
    /* constexpr */ static MPI_Datatype mpi_type()   { return M; } \
    constexpr static bool is_mpi_native_type() { return true; } \
};

MAKE_TRAITS(double, MPI_DOUBLE)
MAKE_TRAITS(float,  MPI_FLOAT)
MAKE_TRAITS(int,    MPI_INT)
MAKE_TRAITS(long int, MPI_LONG)
MAKE_TRAITS(char,   MPI_CHAR)
MAKE_TRAITS(unsigned int, MPI_UNSIGNED)
MAKE_TRAITS(size_t, MPI_UNSIGNED_LONG)
static_assert(sizeof(size_t)==sizeof(unsigned long),
              "size_t and unsigned long are not equivalent");

// Gather individual values of type T from each rank into a std::vector on
// the root rank.
// T must be trivially copyable.
template<typename T>
std::vector<T> gather(T value, int root, MPI_Comm comm) {
    using traits = mpi_traits<T>;
    auto buffer_size = (rank(comm)==root) ? size(comm) : 0;
    std::vector<T> buffer(buffer_size);

    handle_mpi_error("MPI_Gather",
    MPI_Gather( &value,        traits::count(), traits::mpi_type(), // send buffer
                buffer.data(), traits::count(), traits::mpi_type(), // receive buffer
                root, comm));

    return buffer;
}

// Gather individual values of type T from each rank into a std::vector on
// the every rank.
// T must be trivially copyable
template <typename T>
std::vector<T> gather_all(T value, MPI_Comm comm) {
    using traits = mpi_traits<T>;
    std::vector<T> buffer(size(comm));

    handle_mpi_error("MPI_Allgather",
        MPI_Allgather(
            &value,        traits::count(), traits::mpi_type(), // send buffer
            buffer.data(), traits::count(), traits::mpi_type(), // receive buffer
            comm));

    return buffer;
}

// Specialize gather for std::string.
inline std::vector<std::string> gather(std::string str, int root, MPI_Comm comm) {
    using traits = mpi_traits<char>;

    auto counts = gather_all(int(str.size()), comm);
    auto displs = algorithms::make_index(counts);

    std::vector<char> buffer(displs.back());

    // const_cast required for MPI implementations that don't use const* in
    // their interfaces.
    std::string::value_type* ptr = const_cast<std::string::value_type*>(str.data());
    handle_mpi_error("MPI_Gatherv",
        MPI_Gatherv(
            ptr, counts[rank(comm)], traits::mpi_type(),                       // send
            buffer.data(), counts.data(), displs.data(), traits::mpi_type(),   // receive
            root, comm));

    // Unpack the raw string data into a vector of strings.
    std::vector<std::string> result;
    auto nranks = size(comm);
    result.reserve(nranks);
    for (auto i=0; i<nranks; ++i) {
        result.push_back(std::string(buffer.data()+displs[i], counts[i]));
    }
    return result;
}

template <typename T>
std::vector<T> gather_all(const std::vector<T>& values, MPI_Comm comm) {

    using traits = mpi_traits<T>;
    auto counts = gather_all(int(values.size()), comm);
    for (auto& c : counts) {
        c *= traits::count();
    }
    auto displs = algorithms::make_index(counts);

    std::vector<T> buffer(displs.back()/traits::count());
    handle_mpi_error("MPI_Allgatherv",
        MPI_Allgatherv(
            // const_cast required for MPI implementations that don't use const* in their interfaces
            const_cast<T*>(values.data()), counts[rank(comm)], traits::mpi_type(),  // send buffer
            buffer.data(), counts.data(), displs.data(), traits::mpi_type(), // receive buffer
            comm));

    return buffer;
}

/// Gather all of a distributed vector
/// Retains the meta data (i.e. vector partition)
template <typename T>
gathered_vector<T> gather_all_with_partition(const std::vector<T>& values, MPI_Comm comm) {
    using gathered_type = gathered_vector<T>;
    using count_type = typename gathered_vector<T>::count_type;
    using traits = mpi_traits<T>;

    // We have to use int for the count and displs vectors instead
    // of count_type because these are used as arguments to MPI_Allgatherv
    // which expects int arguments.
    auto counts = gather_all(int(values.size()), comm);
    for (auto& c : counts) {
        c *= traits::count();
    }
    auto displs = algorithms::make_index(counts);

    std::vector<T> buffer(displs.back()/traits::count());

    handle_mpi_error("MPI_Allgatherv",
        MPI_Allgatherv(
            // const_cast required for MPI implementations that don't use const* in their interfaces
            const_cast<T*>(values.data()), counts[rank(comm)], traits::mpi_type(), // send buffer
            buffer.data(), counts.data(), displs.data(), traits::mpi_type(), // receive buffer
            comm));

    for (auto& d : displs) {
        d /= traits::count();
    }

    return gathered_type(
        std::move(buffer),
        std::vector<count_type>(displs.begin(), displs.end())
    );
}

template <typename T>
T reduce(T value, MPI_Op op, int root, MPI_Comm comm) {
    using traits = mpi_traits<T>;
    static_assert(traits::is_mpi_native_type(),
                  "can only perform reductions on MPI native types");

    T result;

    handle_mpi_error("MPI_Reduce",
        MPI_Reduce(&value, &result, 1, traits::mpi_type(), op, root, comm));

    return result;
}

template <typename T>
T reduce(T value, MPI_Op op, MPI_Comm comm) {
    using traits = mpi_traits<T>;
    static_assert(traits::is_mpi_native_type(),
                  "can only perform reductions on MPI native types");

    T result;

    MPI_Allreduce(&value, &result, 1, traits::mpi_type(), op, comm);

    return result;
}

template <typename T>
std::pair<T,T> minmax(T value) {
    return {reduce<T>(value, MPI_MIN), reduce<T>(value, MPI_MAX)};
}

template <typename T>
std::pair<T,T> minmax(T value, int root) {
    return {reduce<T>(value, MPI_MIN, root), reduce<T>(value, MPI_MAX, root)};
}

template <typename T>
T broadcast(T value, int root, MPI_Comm comm) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "broadcast can only be performed on trivally copyable types");

    using traits = mpi_traits<T>;

    handle_mpi_error("MPI_Bcast",
        MPI_Bcast(&value, traits::count(), traits::mpi_type(), root, comm));

    return value;
}

template <typename T>
T broadcast(int root, MPI_Comm comm) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "broadcast can only be performed on trivally copyable types");

    using traits = mpi_traits<T>;
    T value;

    handle_mpi_error("MPI_Bcast",
        MPI_Bcast(&value, traits::count(), traits::mpi_type(), root, comm));

    return value;
}

} // namespace mpi
} // namespace arb
