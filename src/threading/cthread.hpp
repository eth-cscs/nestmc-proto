#pragma once

#if !defined(NMC_HAVE_CTHREAD)
    #error "this header can only be loaded if NMC_HAVE_CTHREAD is set"
#endif

#include <thread>
#include <mutex>
#include <algorithm>
#include <array>
#include <chrono>
#include <string>
#include <vector>
#include <type_traits>
#include <functional>
#include <condition_variable>
#include <utility>
#include <unordered_map>
#include <deque>

#include <cstdlib>

namespace nest {
namespace mc {
namespace threading {

// Forward declare task_group at bottom of this header
class task_group;

namespace impl {

using nest::mc::threading::task_group;
using std::mutex;
using lock = std::unique_lock<mutex>;
using std::condition_variable;

using task = std::pair<std::function<void()>, task_group*>;
using task_queue = std::deque<task>;

using thread_list = std::vector<std::thread>;
using thread_map = std::unordered_map<std::thread::id, std::size_t>;

class task_pool {
private:
    // lock and signal on task availability change
    // this is the crucial bit
    mutex _tasks_mutex;
    condition_variable _tasks_available; 

    // fifo of pending tasks
    task_queue _tasks;

    // thread resource
    thread_list _threads;
    // threads -> index
    thread_map _thread_ids;
    // flag to handle exit from all threads
    bool _quit = false;
    
    // internals for taking tasks as a resource
    // and running them (updating above)
    // They get run by a thread in order to consume
    // tasks
    struct run_task;
    // run tasks until a task_group tasks are done
    // for wait
    void run_tasks_while(task_group*);
    // loop forever for secondary threads
    // until quit is set
    void run_tasks_loop();

    // Create nthreads-1 new c std threads
    // must be > 0
    // singled only created in static get_global_task_pool()
    task_pool(std::size_t nthreads);
    
    // task_pool is a singleton 
    task_pool(const task_pool&) = delete;
    task_pool& operator=(const task_pool&) = delete;

    // set quit and wait for secondary threads to end
    ~task_pool();

    // reference task_group specifics (in_flight)
    friend task_group;
  
public:
    // Like tbb calls: run queues a task,
    // wait waits for all tasks to be done
    // cancel dumps queued tasks
    // and is_canceling return true iff
    //   cancel() has been called
    //   but a task is still in flight
    void run(const task&);
    void run(task&&);
    void wait(class task_group*);
  
    // includes master thread
    int get_num_threads() {
        return _threads.size() + 1;
    }

    // get a stable integer for the current thread that
    // is 0..nthreads
    std::size_t get_current_thread() {
        return _thread_ids[std::this_thread::get_id()];
    }
    
    static task_pool& get_global_task_pool();
};
} //impl

///////////////////////////////////////////////////////////////////////
// types
///////////////////////////////////////////////////////////////////////
template <typename T>
class enumerable_thread_specific {
    impl::task_pool& global_task_pool;

    using storage_class = std::vector<T>;
    storage_class data;
  
public :
    using iterator = typename storage_class::iterator;
    using const_iterator = typename storage_class::const_iterator;

    enumerable_thread_specific():
        global_task_pool{impl::task_pool::get_global_task_pool()},
        data{std::vector<T>(global_task_pool.get_num_threads())}
    {}

    enumerable_thread_specific(const T& init):
        global_task_pool{impl::task_pool::get_global_task_pool()},
        data{std::vector<T>(global_task_pool.get_num_threads(), init)}
    {}

    T& local() {
      return data[global_task_pool.get_current_thread()];
    }
    const T& local() const {
      return data[global_task_pool.get_current_thread()];
    }

    auto size() -> decltype(data.size()) const { return data.size(); }

    iterator begin() { return data.begin(); }
    iterator end()   { return data.end(); }

    const_iterator begin() const { return data.begin(); }
    const_iterator end()   const { return data.end(); }

    const_iterator cbegin() const { return data.cbegin(); }
    const_iterator cend()   const { return data.cend(); }
};

template <typename T>
class parallel_vector {
    using value_type = T;
    std::vector<value_type> data_;

private:
    // lock the parallel_vector to update
    impl::mutex mutex;

    // call a function of type X f() in a lock
    template<typename F>
    auto critical(F f) -> decltype(f()) {
        impl::lock lock{mutex};
        return f();
    }

public:
    parallel_vector() = default;
    using iterator = typename std::vector<value_type>::iterator;
    using const_iterator = typename std::vector<value_type>::const_iterator;

    iterator begin() { return data_.begin(); }
    iterator end()   { return data_.end(); }

    const_iterator begin() const { return data_.begin(); }
    const_iterator end()   const { return data_.end(); }

    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator cend()   const { return data_.cend(); }

    void push_back (value_type&& val) {
        critical([&] {
            data_.push_back(std::move(val));
        });
    }
};

inline std::string description() {
    return "CThread Pool";
}

struct timer {
    using time_point = std::chrono::time_point<std::chrono::system_clock>;

    static inline time_point tic() {
        return std::chrono::system_clock::now();
    }

    static inline double toc(time_point t) {
        return std::chrono::duration<double>{tic() - t}.count();
    }

    static inline double difference(time_point b, time_point e) {
        return std::chrono::duration<double>{e-b}.count();
    }
};

constexpr bool multithreaded() { return true; }

class task_group {
private:
    std::size_t in_flight = 0;
    impl::task_pool& global_task_pool;
    // task pool manipulates in_flight
    friend impl::task_pool;
  
public:
    task_group():
        global_task_pool{impl::task_pool::get_global_task_pool()}
    {}
    
    task_group(const task_group&) = delete;
    task_group& operator=(const task_group&) = delete;

    // send function void f() to threads
    template<typename F>
    void run(const F& f) {
        global_task_pool.run(impl::task{f, this});
    }

    template<typename F>
    void run(F&& f) {
        global_task_pool.run(impl::task{std::move(f), this});
    }

    // run function void f() and then wait on all threads in group
    template<typename F>
    void run_and_wait(const F& f) {
        f();
        global_task_pool.wait(this);
    }

    template<typename F>
    void run_and_wait(F&& f) {
        f();
        global_task_pool.wait(this);
    }

    // wait till all tasks in this group are done
    void wait() {
        global_task_pool.wait(this);
    }

    // Make sure that all tasks are done before clean up
    ~task_group() {
        wait();
    }
};

///////////////////////////////////////////////////////////////////////
// algorithms
///////////////////////////////////////////////////////////////////////
struct parallel_for {
    template <typename F>
    static void apply(int left, int right, F f) {
        task_group g;
        for(int i = left; i < right; ++i) {
          g.run([=] {f(i);});
        }
        g.wait();
    }
};

} // threading
} // mc
} // nest

// paralele stable sort uses threading
#include "cthread_parallel_stable_sort.h"

// and sorts use chtread_parallel_stable_sort
namespace nest {
namespace mc {
namespace threading {

template <typename RandomIt>
void sort(RandomIt begin, RandomIt end) {
    pss::parallel_stable_sort(begin, end);
}

template <typename RandomIt, typename Compare>
void sort(RandomIt begin, RandomIt end, Compare comp) {
    pss::parallel_stable_sort(begin, end ,comp);
}

template <typename Container>
void sort(Container& c) {
    pss::parallel_stable_sort(c.begin(), c.end());
}

}
}
}
