#pragma once
#include <string>
#include <random>
#include <vector>
#include "mm_reliability.h"

namespace rrap
{
  enum struct TEST_TYPE {
    MIXED       = 0,
    EXAMPLE_ONE = 1,
    EXAMPLE_TWO = 2,
    TEST_ONE    = 3,
    TEST_TWO    = 4,
    ERROR_TYPE  = 5
  };

  TEST_TYPE get_test_type(int test_type);

  // function to calculate the mean of a set of data
  inline double instance_calculate_mean(const std::vector<double>& data) {
    double sum = 0.0;
    for (const double value : data) {
      sum += value;
    }
    return sum / static_cast<double>(data.size());
  }

  // function to calculate the variance
  inline double instance_calculate_variance(const std::vector<double>& data, const double mean) {
    double variance = 0.0;
    for (const double value : data) {
      variance += (value - mean) * (value - mean);
    }
    return variance /  static_cast<double>(data.size());  // Population variance
  }

  // function to calculate the standard deviation
  inline double instance_calculate_standard_deviation(const double variance) {
    return std::sqrt(variance);
  }

  // function to calculate the Mean Absolute Deviation (MAD)
  inline double instance_calculate_MAD(const std::vector<double>& data, const double mean) {
    double mad = 0.0;
    for (const double value : data) {
      mad += fabs(value - mean);  // absolute deviation
    }
    return mad /  static_cast<double>(data.size());
  }

  inline double instance_compute_median(std::vector<double> data) {
    // sort the data
    std::sort(data.begin(), data.end());

    // if the size is odd, return the middle element
    if (const auto size = data.size(); size % 2 != 0) {
      return data[size / 2];
    } else { // if the size is even, return the average of the two middle elements
      return (data[size / 2 - 1] + data[size / 2]) / 2.0;
    }
  }

  template<typename T>
  void compute_tree_lengths(const std::vector<T>&  systems, unsigned& min, unsigned& max, double& avg,
          unsigned& total, unsigned& branches)
  {
    auto min_length = std::numeric_limits<unsigned>::max();
    auto max_length = std::numeric_limits<unsigned>::min();
    auto total_length = 0u;
    auto num_of_branches = 0u;
    for (auto l = 1u; l < systems.size(); ++l) {
      if (systems[l]->m_children.empty()) {
        // leaf node
        auto length = 1u;
        auto f = systems[l]->m_father->m_id;
        while (f != 0) {
          f = systems[f]->m_father->m_id;
          ++length;
        }
        min_length = std::min(min_length, length);
        max_length = std::max(max_length, length);
        total_length += length;
        ++num_of_branches;
      }
    }
    const auto avg_length = static_cast<double>(total_length) / static_cast<double>(num_of_branches);

    min = min_length;
    max = max_length;
    total = total_length;
    avg = avg_length;
    branches = num_of_branches;
  }

  struct ProbInstance
  {
    std::string inst_file;

    size_t num_resource_types{};
    size_t num_subsystems{};
    size_t num_component_types{};

    TEST_TYPE test_type;
    int system_identifier;

    std::vector<double> amount_of_each_resources;
    std::vector<std::vector<double>> component_reliability_at_subsystem;
    std::vector<std::vector<std::vector<double>>> resource_for_component_at_subsystem;

    double min_component_resource;
    std::vector<size_t> num_component_types_at_subsystem;
    std::vector<std::vector<int>> component_lb_at_subsystem;
    std::vector<std::vector<int>> component_ub_at_subsystem;

    // Q-form reliability function
    QSystemReliability system_reliability;
    QSystemTree logic_system_tree;
    QSortedSystemTree sorted_logic_system_tree;

    // general-form reliability function
    QRSystemReliability qr_system_reliability;
    QRSystemTree qr_system_tree;
    QRSortedSystemTree qr_sorted_system_tree;

    std::set<std::pair<int,int>> dominated_components;

    ProbInstance(int test_type_, int system_index_, const std::string &inst_file);


  private:
    void read_data(const std::string &inst_file);
    void calculate_component_lb_at_subsystem();
    void calculate_component_ub_at_subsystem();
    void check_component_domination ();
    auto calculate_min_resource_usage(int j) -> std::vector<double>;
    auto compute_system_reliability() -> QSystemReliability;
    auto compute_qr_system_reliability() -> QRSystemReliability;
  };

  std::ostream &operator<<(std::ostream &out, const ProbInstance &inst);
} // namespace rrap
