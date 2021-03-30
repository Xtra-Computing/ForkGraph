/* Copyright (C) DataInz Technologies Pte. Ltd. - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited.
 * Proprietary and confidential.
 */

#ifndef MARLINDB_TIMEMEASURER_HPP
#define MARLINDB_TIMEMEASURER_HPP
#include <iostream>
#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::microseconds;
using std::chrono::nanoseconds;

class Timer {
 public:
  Timer() {
    ns_count_ = 0;
  }
  ~Timer() {}

  void tic() {
    start_time_ = high_resolution_clock::now();
  }

  void toc() {
    end_time_ = high_resolution_clock::now();
  }

  void resume() {
    tic();
  }
  void pause() {
    toc();
    ns_count_ += time_ns();
  }

  long long time_ms() const {
    return std::chrono::duration_cast<milliseconds>(end_time_ - start_time_).count();
  }

  long long time_us() const {
    return std::chrono::duration_cast<microseconds>(end_time_ - start_time_).count();
  }

  long long time_ns() const {
    return std::chrono::duration_cast<nanoseconds>(end_time_ - start_time_).count();
  }

  long long get_total_ns() const{
    return ns_count_;
  }

  void print_s(std::string msg) const {
    std::cout << msg << ": ";
    print_s();
  }

  void print_s() const {
    std::cout << std::chrono::duration_cast<milliseconds>(end_time_ - start_time_).count() / 1000 << " s" << std::endl;
  }

  void print_ms(std::string msg) const {
    std::cout << msg << ": ";
    print_ms();
  }

  void print_ms() const {
    std::cout << std::chrono::duration_cast<milliseconds>(end_time_ - start_time_).count() << " ms" << std::endl;
  }

  void print_us(std::string msg) const {
    std::cout << msg << ": ";
    print_us();
  }
  void print_us() const {
    std::cout << std::chrono::duration_cast<microseconds>(end_time_ - start_time_).count() << " us" << std::endl;
  }

  void print_ns(std::string msg) const {
    std::cout << msg << ": ";
    print_ns();
  }
  void print_ns() const {
    std::cout << std::chrono::duration_cast<nanoseconds>(end_time_ - start_time_).count() << " ns" << std::endl;
  }

 private:
  Timer(const Timer &);
  Timer &operator=(const Timer &);

 private:
  high_resolution_clock::time_point start_time_;
  high_resolution_clock::time_point end_time_;
  long long ns_count_;
};
#endif //MARLINDB_TIMEMEASURER_HPP
