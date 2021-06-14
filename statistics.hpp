#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include "compat.hpp"
#include <atomic>
#include <string>
#include <sstream>
#include <ostream>
#include <type_traits>


// The file is taken from the g6k library https://github.com/fplll/g6k
// and has been appropriately modified to meet the needs of the
// Polard-type mutithreaded walk

//
// To enable statistics run make with make FLAGS+=-DENABLE_STATS=1
//

// If you want to add more stattics:
//  -- Add the members variables to the class (private and public)
//  -- Add / Modify the getters / incrementers
//  -- Add it to clear_statistics
//  -- Optionlly: add it to print_statistics

#if defined ENABLE_STATS
    #define COLLECT_STATISTICS 1
#else
    #define COLLECT_STATISTICS 0
#endif

// prints macro values, TODO: remove
//#define XSTR(x) STR(x)
//#define STR(x) #x
//#pragma message "The value of COLLECT_STATISTICS: " XSTR(COLLECT_STATISTICS)


// collect overjumps that can occure in JumpOneStep function
#ifndef COLLECT_STATISTICS_OVERJUMPS
  #define COLLECT_STATISTICS_OVERJUMPS COLLECT_STATISTICS
#endif



class MCTStatistics
{

private:
#if COLLECT_STATISTICS_OVERJUMPS
    std::atomic_ullong stats_overjumps_B1min;
    std::atomic_ullong stats_overjumps_B1max;
    std::atomic_ullong stats_overjumps_B2min;
    std::atomic_ullong stats_overjumps_B2max;
#else
    static constexpr unsigned long long stats_overjumps_B1min = 0;
    static constexpr unsigned long long stats_overjumps_B1max = 0;
    static constexpr unsigned long long stats_overjumps_B2min = 0;
    static constexpr unsigned long long stats_overjumps_B2max = 0;
#endif

#define MAKE_ATOMIC_INCREMENTER(INCNAME, STAT) \
void inc_stats_ ## INCNAME( mystd::decay_t<decltype(stats_##STAT.load())> how_much = 1) noexcept { stats_##STAT.fetch_add(how_much, std::memory_order_relaxed); }

#define MAKE_ATOMIC_DECREMENTER(INCNAME, STAT) \
void dec_stats_ ## INCNAME( mystd::decay_t<decltype(stats_##STAT.load())> how_much = 1) noexcept { stats_##STAT.fetch_sub(how_much, std::memory_order_relaxed); }

#define MAKE_ATOMIC_GETTER(GETTERNAME, STAT) \
FORCE_INLINE mystd::decay_t<decltype(stats_##STAT.load())> get_stats_##GETTERNAME() const noexcept { return stats_##STAT.load(); }

#define MAKE_ATOMIC_SETTER(SETTERNAME, STAT) \
void set_stats_##SETTERNAME(mystd::decay_t<decltype(stats_##STAT.load())> const new_val) noexcept {stats_##STAT.store(new_val); }

#define MAKE_NOOP_INCREMENTER(INCNAME) \
template<class Arg> FORCE_INLINE static void inc_stats_##INCNAME(Arg) noexcept {} \
FORCE_INLINE static void inc_stats_##INCNAME() noexcept {}

#define MAKE_NOOP_DECREMENTER(INCNAME) \
template<class Arg> FORCE_INLINE static void dec_stats_##INCNAME(Arg) noexcept {} \
FORCE_INLINE static void dec_stats_##INCNAME() noexcept {}

#define MAKE_NOOP_SETTER(SETTERNAME) \
template<class Arg> FORCE_INLINE static void set_stats_##SETTERNAME(Arg) noexcept {} \
FORCE_INLINE static void set_stats_##SETTERNAME() noexcept {}


/** Totally evil hackery to work around lack of C++ constexpr if (or #if's inside macro definitions...)
    Some gcc version might actually allow #if's inside macros, but we prefer portability. **/

#define MAKE_INCREMENTER_FOR(INCNAME, STAT, NONTRIVIAL) MAKE_INCREMENTER_AUX(INCNAME, STAT, NONTRIVIAL) // to macro-expand the name "NONTRIVIAL", such that token-pasting in the following macro is done AFTER macro expansion.
#define MAKE_INCREMENTER(INCNAME, NONTRIVIAL) MAKE_INCREMENTER_AUX(INCNAME, INCNAME, NONTRIVIAL)
#define MAKE_INCREMENTER_AUX(INCNAME, STAT, NONTRIVIAL) MAKE_INCREMENTER_##NONTRIVIAL(INCNAME, STAT) // This is evil
#define MAKE_INCREMENTER_0(INCNAME, STAT) MAKE_NOOP_INCREMENTER(INCNAME)
#define MAKE_INCREMENTER_1(INCNAME, STAT) MAKE_ATOMIC_INCREMENTER(INCNAME, STAT)
#define MAKE_INCREMENTER_2(INCNAME, STAT) MAKE_ATOMIC_INCREMENTER(INCNAME, STAT)
#define MAKE_INCREMENTER_3(INCNAME, STAT) MAKE_ATOMIC_INCREMENTER(INCNAME, STAT)
#define MAKE_INCREMENTER_4(INCNAME, STAT) MAKE_ATOMIC_INCREMENTER(INCNAME, STAT)

#define MAKE_DECREMENTER(NAME, NONTRIVIAL) MAKE_DECREMENTER_AUX(NAME, NAME, NONTRIVIAL)
#define MAKE_DECREMENTER_AUX(NAME, STAT, NONTRIVIAL) MAKE_DECREMENTER_##NONTRIVIAL(NAME, STAT)
#define MAKE_DECREMENTER_0(NAME, STAT) MAKE_NOOP_DECREMENTER(NAME)
#define MAKE_DECREMENTER_1(NAME, STAT) MAKE_ATOMIC_DECREMENTER(NAME, STAT)
#define MAKE_DECREMENTER_2(NAME, STAT) MAKE_ATOMIC_DECREMENTER(NAME, STAT)
#define MAKE_DECREMENTER_3(NAME, STAT) MAKE_ATOMIC_DECREMENTER(NAME, STAT)
#define MAKE_DECREMENTER_4(NAME, STAT) MAKE_ATOMIC_DECREMENTER(NAME, STAT)

#define MAKE_GETTER_FOR(GETTERNAME, STAT, NONTRIVIAL) MAKE_GETTER_AUX(GETTERNAME, STAT, NONTRIVIAL)
#define MAKE_GETTER(GETTERNAME, NONTRIVIAL) MAKE_GETTER_AUX(GETTERNAME, GETTERNAME, NONTRIVIAL)
#define MAKE_GETTER_AUX(GETTERNAME, STAT, NONTRIVIAL) MAKE_GETTER_##NONTRIVIAL(GETTERNAME, STAT)
#define MAKE_GETTER_0(GETTERNAME, STAT) \
FORCE_INLINE static constexpr auto get_stats_##GETTERNAME() noexcept -> mystd::remove_cv_t<decltype(stats_##STAT)> { return stats_##STAT; }
#define MAKE_GETTER_1(GETTERNAME, STAT) MAKE_ATOMIC_GETTER(GETTERNAME, STAT)
#define MAKE_GETTER_2(GETTERNAME, STAT) MAKE_ATOMIC_GETTER(GETTERNAME, STAT)
#define MAKE_GETTER_3(GETTERNAME, STAT) MAKE_ATOMIC_GETTER(GETTERNAME, STAT)
#define MAKE_GETTER_4(GETTERNAME, STAT) MAKE_ATOMIC_GETTER(GETTERNAME, STAT)

#define MAKE_SETTER_FOR(SETTERNAME, STAT, NONTRIVIAL) MAKE_SETTER_AUX(SETTERNAME, STAT, NONTRIVIAL)
#define MAKE_SETTER(SETTERNAME, NONTRIVIAL) MAKE_SETTER_AUX(SETTERNAME, SETTERNAME, NONTRIVIAL)
#define MAKE_SETTER_AUX(SETTERNAME, STAT, NONTRIVIAL) MAKE_SETTER_##NONTRIVIAL(SETTERNAME, STAT)
#define MAKE_SETTER_0(SETTERNAME, STAT) MAKE_NOOP_SETTER(SETTERNAME)
#define MAKE_SETTER_1(SETTERNAME, STAT) MAKE_ATOMIC_SETTER(SETTERNAME, STAT)
#define MAKE_SETTER_2(SETTERNAME, STAT) MAKE_ATOMIC_SETTER(SETTERNAME, STAT)
#define MAKE_SETTER_3(SETTERNAME, STAT) MAKE_ATOMIC_SETTER(SETTERNAME, STAT)
#define MAKE_SETTER_4(SETTERNAME, STAT) MAKE_ATOMIC_SETTER(SETTERNAME, STAT)

#define MAKE_GETTER_AND_INCREMENTER(NAME, NONTRIVIAL) \
MAKE_INCREMENTER(NAME, NONTRIVIAL) \
MAKE_GETTER(NAME, NONTRIVIAL)

public:
  static constexpr int collect_statistics_level = COLLECT_STATISTICS;

  static constexpr bool collect_statistics_overjumps_total  = (COLLECT_STATISTICS_OVERJUMPS >= 1);
  static constexpr bool collect_statistics_overjumps_B1min  = (COLLECT_STATISTICS_OVERJUMPS >= 1);
  static constexpr bool collect_statistics_overjumps_B1max  = (COLLECT_STATISTICS_OVERJUMPS >= 1);
  static constexpr bool collect_statistics_overjumps_B2min  = (COLLECT_STATISTICS_OVERJUMPS >= 1);
  static constexpr bool collect_statistics_overjumps_B2max  = (COLLECT_STATISTICS_OVERJUMPS >= 1);

  MAKE_GETTER_AND_INCREMENTER(overjumps_B1min, COLLECT_STATISTICS_OVERJUMPS)
  MAKE_GETTER_AND_INCREMENTER(overjumps_B1max, COLLECT_STATISTICS_OVERJUMPS)
  MAKE_GETTER_AND_INCREMENTER(overjumps_B2min, COLLECT_STATISTICS_OVERJUMPS)
  MAKE_GETTER_AND_INCREMENTER(overjumps_B2max, COLLECT_STATISTICS_OVERJUMPS)
  unsigned long long get_stats_overjumps_total() const
  {
    return get_stats_overjumps_B1min() + get_stats_overjumps_B1max()+
    get_stats_overjumps_B2min() + get_stats_overjumps_B2max();
  }

  inline void clear_statistics() noexcept
  {
    #if COLLECT_STATISTICS_OVERJUMPS
        stats_overjumps_B1min = 0;
        stats_overjumps_B1max = 0;
        stats_overjumps_B2min = 0;
        stats_overjumps_B2max = 0;
    #endif

  }

  void call_print_statistics(std::ostream &os = std::cout)
  {
        return print_statistics(os);
  }

#define STATS_PRINT_IF(NAME, STRING)\
if(collect_statistics_##NAME) { os << STRING << get_stats_##NAME(); }

    void print_statistics(std::ostream &os = std::cout)
    {
      os << "Statistics:" << std::endl;
      if(collect_statistics_overjumps_total)
      {
          os << "Overjumps: " << get_stats_overjumps_total();
          if(collect_statistics_overjumps_B1min)
              os << ", overjumps_B1min: " << get_stats_overjumps_B1min();
          if(collect_statistics_overjumps_B1max)
              os << ", overjumps_B1max: " << get_stats_overjumps_B1max();
          if(collect_statistics_overjumps_B2min)
              os << ", overjumps_B1min: " << get_stats_overjumps_B2min();
          if(collect_statistics_overjumps_B2max)
              os << ", overjumps_B1max: " << get_stats_overjumps_B2max();
          os << "\n";
      }


    }

    //std::string get_statistics_string(int alg) { return get_statistics_string(static_cast<StatisticsOutputForAlg>(alg)); }
    std::string get_statistics_string()
    {
        std::stringstream out;
        print_statistics(out);
        return out.str();
    }


}; //class MCTStatistics

/**
    MergeOnExit<Int,Functor> wraps an Int x (with very limited functionality as of now) and Functor fun
    and call fun(x) upon destruction.
    The intented use case is to create a thread-local counter local_counter_foo
    with a lambda fun = [](Int x){ global_counter_foo += x; }

    To pass parameters, use the merge_on_exit helper function.
    E.g.
    auto &&local_counter_foo = merge_on_exit<int>( [](int x){global_counter_foo += x;} );
    auto &&local_statistic_foo = merge_on_exit<unsigned long> ( [this](unsigned long x) { this->statistics.inc_stats_foo(x); } );
    (where this->statistics is of type SieveStatistics )
    Note that the && is MANDATORY pre C++17 (since MergeOnExit is non-movable, non-copyable, this
    uses lifetime-extension of temporaries).
*/
template<class Integer, class Functor>
class MergeOnExit
{
    static_assert(std::is_integral<Integer>::value, "Wrong template argument");
    // Alas, is_invocable<Functor> is C++17...
    public:
    Functor const fun;
    Integer val;
    constexpr MergeOnExit(Functor fun_) : fun(fun_), val(0) {}
    constexpr MergeOnExit(Functor fun_, Integer val_) : fun(fun_), val(val_) {}
    MergeOnExit &operator=(Integer new_val) { val = new_val; return *this; }
    MergeOnExit &operator=(MergeOnExit const &) = delete;
    MergeOnExit &operator=(MergeOnExit &&) = delete;
    MergeOnExit(MergeOnExit const &) = delete;
    MergeOnExit(MergeOnExit &&) = delete;
    template<class Arg> MergeOnExit &operator+=(Arg &&arg) { val+=(std::forward<Arg>(arg)); return *this; }
    template<class Arg> MergeOnExit &operator-=(Arg &&arg) { val+=(std::forward<Arg>(arg)); return *this; }
    template<class Arg> MergeOnExit &operator*=(Arg &&arg) { val+=(std::forward<Arg>(arg)); return *this; }
    Integer operator++(int) { return val++; }
    MergeOnExit& operator++() { ++val; return *this; }
    Integer operator--(int) { return val--; }
    MergeOnExit& operator--() { --val; return *this; }
    ~MergeOnExit()
    {
        fun(val);
    }
    constexpr operator Integer() const { return val; }
};

// Lacking deduction guides (pre C++17), to deduce template parameters, we use a helper function:
template<class Integer, class Functor> MergeOnExit<Integer, Functor> merge_on_exit(Functor const &fun_) { return {fun_}; }
template<class Integer, class Functor> MergeOnExit<Integer, Functor> merge_on_exit(Functor const &fun_, Integer val_) { return {fun_, val_}; }




#endif
