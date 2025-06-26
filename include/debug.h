// debug.h
#pragma once

#include <iostream>
#include <mutex>
#include <chrono>
#include <iomanip>
#include <sstream>

extern std::mutex debug_mutex;

#ifdef DEBUG
    #define DEBUG_PRINT(x) do { \
        std::lock_guard<std::mutex> lock(debug_mutex); \
        auto now = std::chrono::system_clock::now(); \
        std::time_t now_time = std::chrono::system_clock::to_time_t(now); \
        std::ostringstream oss; \
        oss << std::put_time(std::localtime(&now_time), "%Y-%m-%d %H:%M:%S"); \
        std::cerr << oss.str() << " - " << x << std::endl; \
    } while (0)
#else
    #define DEBUG_PRINT(x)
#endif
