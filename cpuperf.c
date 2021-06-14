#include <iostream>
#include <vector>

#ifdef PERFORMANCE_COUNTING
std::vector<uint64_t> perfcounters(1000,0);
std::vector<std::atomic<uint64_t>> atomic_perfcounters(1000);
cpu::performance_counter_manager perfmanager;

struct cpuperf_init_t {
    cpuperf_init_t()
    {
        /* LMPCT.h */
        perfmanager.add_performance_counter(atomic_perfcounters[101], "SelectRandomStartingPoint");
        perfmanager.add_performance_counter(atomic_perfcounters[102], "hashvalue");
        perfmanager.add_performance_counter(atomic_perfcounters[103], "hashvalue2");
    }
} cpuperf_init;

void show_cpu_stats()
{
    uint64_t now = cpu::cpu_timestamp();
    uint64_t total = *perfmanager._counters[0];
    if (total > uint64_t(1)<<57)
        total += now;
    for (unsigned i = 0; i < perfmanager._counters.size(); ++i)
    {
        uint64_t cnt = *perfmanager._counters[i];
        if (0 == cnt)
            continue;
        if (cnt > uint64_t(1)<<57)
            cnt += now;
        std::cout << i << "\t: " << perfmanager._descriptions[i] << std::endl;
        std::cout << i << "\t: " << double(cnt)*100.0/double(total) << " %, \t cycles=" << cnt << std::endl;
    }
    for (unsigned i = 0; i < perfmanager._atomic_counters.size(); ++i)
    {
        uint64_t cnt = *perfmanager._atomic_counters[i];
        if (0 == cnt)
            continue;
        if (cnt > uint64_t(1)<<57)
            cnt += now;
        std::cout << i + perfmanager._counters.size() << "\t: " << perfmanager._atomic_descriptions[i] << std::endl;
        std::cout << i + perfmanager._counters.size() << "\t: " << double(cnt)*100.0/double(total) << " %, \t cycles=" << cnt << std::endl;
    }
}

#else
void show_cpu_stats()
{
}
#endif
