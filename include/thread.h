#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

// 线程池类
class ThreadPool {
public:
    explicit ThreadPool(size_t num_threads);
    ~ThreadPool();

    // 将任务加入线程池队列并返回结果
    template <class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<typename std::invoke_result<F, Args...>::type>;

private:
    std::vector<std::thread> workers;               // 工作线程
    std::queue<std::function<void()>> tasks;        // 任务队列
    std::mutex queue_mutex;                         // 任务队列的互斥锁
    std::condition_variable condition;              // 条件变量
    bool stop;                                      // 停止标志
};

// 模板方法的实现必须放在头文件中
template <class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) -> std::future<typename std::invoke_result<F, Args...>::type> {
    using return_type = typename std::invoke_result<F, Args...>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    std::future<return_type> res = task->get_future();

    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if (stop) {
            throw std::runtime_error("enqueue on stopped ThreadPool");
        }
        tasks.emplace([task]() { (*task)(); });
    }
    condition.notify_one();
    return res;
}

#endif // THREAD_POOL_H