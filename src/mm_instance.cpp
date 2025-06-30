#include <fstream>
#include <iostream>
#include <cmath>
#include "mm_instance.h"

namespace rrap {
    TEST_TYPE get_test_type(const int test_type) {
        if (test_type == 0) {
            return TEST_TYPE::MIXED;
        }
        if (test_type == 1) {
            return TEST_TYPE::EXAMPLE_ONE;
        }
        if (test_type == 2) {
            return TEST_TYPE::EXAMPLE_TWO;
        }
        if (test_type == 3) {
            return TEST_TYPE::TEST_ONE;
        }
        if (test_type == 4) {
            return TEST_TYPE::TEST_TWO;
        }
        return TEST_TYPE::ERROR_TYPE;
    }

    ProbInstance::ProbInstance(const int test_type_, const int system_index_, const std::string &inst_file) : inst_file(inst_file), system_identifier(system_index_) {
        test_type = get_test_type(test_type_);

        if (test_type == TEST_TYPE::ERROR_TYPE) {
            std::cerr << "ERROR TYPE" << std::endl;
        }

        // read data from file
        read_data(inst_file);

        // compute the minimum amount of resource required for a single component
        min_component_resource = std::numeric_limits<double>::max();
        for (int i = 0; i < num_resource_types; ++i)
        {
            for (int j = 0; j < num_subsystems; ++j)
            {
                for (int h = 0; h < num_component_types_at_subsystem[j]; h++)
                {
                    min_component_resource = std::min(min_component_resource, resource_for_component_at_subsystem[i][j][h]);
                }
            }
        }

        // collect the dominated component types for each subsystem
        check_component_domination();

        // initialize the reliability function
        system_reliability = compute_system_reliability();
        qr_system_reliability = compute_qr_system_reliability();

        logic_system_tree = QSystemTree(system_reliability);
        sorted_logic_system_tree = QSortedSystemTree(system_reliability);
        qr_system_tree = QRSystemTree(qr_system_reliability);
        qr_sorted_system_tree = QRSortedSystemTree(qr_system_reliability);
    }


    void ProbInstance::read_data(const std::string &inst_file)
    {
        std::ifstream ifs{inst_file};
        if (ifs.fail())
        {
            std::cerr << "ERROR: cannot open instance file " << inst_file << "\n";
            std::exit(EXIT_FAILURE);
        }

        // read 1st line:
        ifs >> num_resource_types;
        ifs >> num_subsystems;
        ifs >> num_component_types;

        // assign the number of component types at each subsystem
        num_component_types_at_subsystem.resize(num_subsystems, 0);
        for (auto j = 0u; j < num_subsystems; ++j) {
            num_component_types_at_subsystem[j] = num_component_types;
        }
        
        // read 2nd line: the amount of each type of resource
        amount_of_each_resources.resize(num_resource_types);
        for (auto i = 0u; i < num_resource_types; ++i)
        {
            ifs >> amount_of_each_resources[i];
        }

        // - the examples E1 and E2 are hard-coded in the implementation
        // - the data for tests T1 and T2 are the same with MIXED, except that the lower band upper bounds of
        // -- the number of components for T1 and T2 are given as inputs while those for MIXED are computed
        if (test_type == TEST_TYPE::TEST_ONE || test_type == TEST_TYPE::TEST_TWO || test_type == TEST_TYPE::MIXED) {
            // for each subsystem, read the reliability of each component-type: r_jh
            component_reliability_at_subsystem.resize(num_subsystems);
            for (auto j = 0u; j < num_subsystems; ++j)
            {
                component_reliability_at_subsystem[j].resize(num_component_types_at_subsystem[j]);
            }
            for (auto j = 0u; j < num_subsystems; ++j)
            {
                for (auto h = 0u; h < num_component_types_at_subsystem[j]; ++h)
                {
                    ifs >> component_reliability_at_subsystem[j][h];
                }
            }

            // for each type of resource, for each subsystem, read the resource of each component-type: r_jh
            resource_for_component_at_subsystem.resize(num_resource_types);
            for (auto i = 0u; i < num_resource_types; ++i)
            {
                resource_for_component_at_subsystem[i].resize(num_subsystems);
                for (auto j = 0u; j < num_subsystems; ++j)
                {
                    resource_for_component_at_subsystem[i][j].resize(num_component_types_at_subsystem[j]);
                }
            }
            for (int i = 0; i < num_resource_types; ++i)
            {
                for (int j = 0; j < num_subsystems; ++j)
                {
                    for (int h = 0; h < num_component_types_at_subsystem[j]; h++)
                    {
                        ifs >> resource_for_component_at_subsystem[i][j][h];
                    }
                }
            }
        }

        // for Ha and Kuo (2006)'s examples and tests, the lower band upper bounds are given as inputs
        if (test_type == TEST_TYPE::EXAMPLE_ONE || test_type == TEST_TYPE::EXAMPLE_TWO || test_type == TEST_TYPE::TEST_ONE || test_type == TEST_TYPE::TEST_TWO) {
            component_lb_at_subsystem.resize(num_subsystems);
            for (int j = 0; j < num_subsystems; j++)
            {
                component_lb_at_subsystem[j].resize(num_component_types_at_subsystem[j]);
                for (int h = 0; h < num_component_types_at_subsystem[j]; ++h) {
                    ifs >> component_lb_at_subsystem[j][h];
                }
            }

            component_ub_at_subsystem.resize(num_subsystems);
            for (int j = 0; j < num_subsystems; j++)
            {
                component_ub_at_subsystem[j].resize(num_component_types_at_subsystem[j]);
                for (int h = 0; h < num_component_types_at_subsystem[j]; ++h) {
                    ifs >> component_ub_at_subsystem[j][h];
                }
            }
        }
        // for MIXED instances, the lower band upper bounds are computed
        else if (test_type == TEST_TYPE::MIXED) {
            calculate_component_lb_at_subsystem();
            calculate_component_ub_at_subsystem();
        }
    }

    void ProbInstance::calculate_component_lb_at_subsystem() {
        component_lb_at_subsystem.resize(num_subsystems);

        for (int j = 0; j < num_subsystems; j++)
        {
            component_lb_at_subsystem[j].resize(num_component_types_at_subsystem[j], 0);
        }
    }

    void ProbInstance::calculate_component_ub_at_subsystem() {
        component_ub_at_subsystem.resize(num_subsystems);

        for (int j = 0; j < num_subsystems; j++)
        {
            const auto minUsage = calculate_min_resource_usage(j);
            component_ub_at_subsystem[j].resize(num_component_types_at_subsystem[j],std::numeric_limits<int>::max());

            for (int h = 0; h < num_component_types_at_subsystem[j]; h++)
            {
                for (int i = 0; i < num_resource_types; i++)
                {
                    component_ub_at_subsystem[j][h] = std::min(component_ub_at_subsystem[j][h],
                                                               static_cast<int>(std::floor((amount_of_each_resources[i] - minUsage[i]) / resource_for_component_at_subsystem[i][j][h])));
                }
            }
        }
    }

    std::vector<double> ProbInstance::calculate_min_resource_usage(const int j)
    {
        std::vector<double> min_res_usage(num_resource_types, 0.0);

        for (int i = 0; i < num_resource_types; ++i)
        {
            for (int ja = 0; ja < num_subsystems; ++ja)
            {
                if (ja == j) continue;

                min_res_usage[i] += *std::min_element(resource_for_component_at_subsystem[i][ja].begin(),
                    resource_for_component_at_subsystem[i][ja].end());
            }

            /* original: min_res_usage[i] = Math.Round(min_res_usage[i], 2); */
            min_res_usage[i] = std::round(min_res_usage[i] * 100) / 100.0;
        }

        return min_res_usage;
    }

    void ProbInstance::check_component_domination() {
        for (auto j = 0; j < num_subsystems; j++) {
            for (auto h1 = 0u; h1 < num_component_types_at_subsystem[j]; ++h1) {
                for (auto h2 = h1 + 1u; h2 < num_component_types_at_subsystem[j]; ++h2) {

                    if (component_reliability_at_subsystem[j][h1] <= component_reliability_at_subsystem[j][h2]) {
                        auto is_dominated = true;
                        for (auto i = 0u; i < num_resource_types; ++i) {
                            if (resource_for_component_at_subsystem[i][j][h1] <
                                resource_for_component_at_subsystem[i][j][h2]) {
                                is_dominated = false;
                                break;
                                }
                        }
                        if (is_dominated) {
                            dominated_components.emplace(j, h1);
#ifndef NDEBUG
                            std::cout << j << ": " << h1 << " << " << h2 << std::endl;
#endif
                            continue;
                        }
                    }


                    if (component_reliability_at_subsystem[j][h2] <= component_reliability_at_subsystem[j][h1]) {
                        auto is_dominated = true;
                        for (auto i = 0u; i < num_resource_types; ++i) {
                            if (resource_for_component_at_subsystem[i][j][h2] <
                                resource_for_component_at_subsystem[i][j][h1]) {
                                is_dominated = false;
                                break;
                                }
                        }
                        if (is_dominated) {
                            dominated_components.emplace(j, h2);
#ifndef NDEBUG
                            std::cout << j << ": " << h2 << " << " << h1 << std::endl;
#endif
                        }
                    }
                }
            }
        }
#ifndef NDEBUG
        std::cout << dominated_components.size() << " components are dominated\n";
        for (auto p : dominated_components) {
            std::cout << p.first << " " << p.second << std::endl;
        }
#endif
    }


    auto ProbInstance::compute_system_reliability() -> QSystemReliability {
        std::vector<QSystemReliability> Q;
        std::vector<QSystemReliability> R;
        Q.reserve(num_subsystems);
        R.reserve(num_subsystems);

        for (auto i = 0u; i < num_subsystems; ++i) {
            Q.emplace_back(num_subsystems, i, 1);
            R.push_back(1 - Q[i]);
        }

        switch (system_identifier) {
            case 1: // 5 subsystems
            {
                auto S1 = R[4] * (1 - Q[0] * Q[2]) * (1 - Q[1] * Q[3]) +
                          Q[4] * (1 - (1 - R[0] * R[1]) * (1 - R[2] * R[3]));
                std::cout << "S1: " << S1 << std::endl;
                return S1;
            }
            case 2: // 5 subsystems
            {
                auto S2 = R[4] * (1 - Q[1] * Q[3]) + Q[4] * (1 - (1 - R[0] * R[1]) * (1 - R[2] * R[3]));
                std::cout << "S2: " << S2 << std::endl;
                return S2;
            }
            case 3: // 6 subsystems
            {
                auto S3 = R[5] * (1 - Q[0] * Q[2] * Q[3]) * (1 - Q[1] * Q[4]) +
                          Q[5] * (1 - (1 - R[0] * R[1]) * (1 - (1 - Q[2] * Q[3]) * R[4]));
                std::cout << "S3: " << S3 << std::endl;
                return S3;
            }
            case 4: // 7 subsystems
            {
                auto Q234 = Q[2] + R[2] * Q[3] * Q[4];
                auto R234 = 1 - Q234;
                auto S4 = Q[6] * (1 - (1 - R[0] * R[1]) * (1 - R234 * R[5])) +
                          R[6] * (1 - Q[0] * Q234) * (1 - Q[1] * Q[5]);
                std::cout << "S4: " << S4 << std::endl;
                return S4;
            }
            case 5: { // 7 subsystems
                auto S5 = Q[2] * Q[4] * (1 - (1 - R[0] * R[3] * R[5]) * (1 - R[1] * R[6]))
                          + R[2] * R[4] * (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6])
                          + R[2] * Q[4] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[5]) * Q[6])
                          + Q[2] * R[4] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[5] * Q[6]);
                std::cout << "S5: " << S5 << std::endl;
                return S5;
            }
            case 6: { // 8 subsystems
                auto S6 = Q[3] * Q[5] * (1 - (1 - R[0] * R[2] * R[6]) * (1 - R[1] * R[4] * R[7]))
                          + R[3] * R[5] * (1 - Q[0] * (1 - R[1] * R[4])) * (1 - Q[6] * Q[7])
                          + R[3] * Q[5] * (1 - Q[0] * (1 - R[1] * R[4])) * (1 - (1 - R[2] * R[6]) * Q[7])
                          + Q[3] * R[5] * (1 - (1 - R[0] * R[2]) * (1 - R[1] * R[4])) * (1 - Q[6] * Q[7]);
                std::cout << "S6: " << S6 << std::endl;
                return S6;
            }
            case 7: { // 8 subsystems
                auto S7 = Q[3] * Q[4] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[6]) * (1 - R[2] * R[7]))
                          + R[3] * R[4] * (1 - Q[0] * Q[1] * Q[2]) * (1 - Q[5] * Q[6] * Q[7])
                          + R[3] * Q[4] * (1 - (1 - R[2] * R[7]) * (1 - (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6])))
                          + Q[3] * R[4] * (1 - (1 - R[0] * R[5]) * (1 - (1 - Q[1] * Q[2]) * (1 - Q[6] * Q[7])));
                std::cout << "S7: " << S7 << std::endl;
                return S7;
            }
            case 8: { // 9 subsystems
                auto S8 = Q[4] * Q[2] * Q[6] * (1 - (1 - R[0] * R[3] * R[7]) * (1 - R[1] * R[5] * R[8]))
                          + Q[4] * R[2] * Q[6] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[7]) * (1 - R[5] * R[8]))
                          + Q[4] * Q[2] * R[6] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) * (1 - Q[7] * Q[8])
                          + Q[4] * R[2] * R[6] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) * (1 - Q[7] * Q[8])
                          + R[4] * Q[2] * Q[6] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - (1 - R[5] * R[8]) * Q[7])
                          + R[4] * R[2] * Q[6] * (1 - Q[0] * Q[1]) * (1 - (1 - R[5] * R[8]) * Q[7])
                          + R[4] * Q[2] * R[6] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[7] * Q[8])
                          + R[4] * R[2] * R[6] * (1 - Q[0] * Q[1]) * (1 - Q[7] * Q[8]);
                std::cout << "S8: " << S8 << std::endl;
                return S8;
            }
            case 9: // 10 subsystems
            {
                auto S9 = (1 - (1 - (1 - Q[2] * (1 - R[0] * R[1])) * R[3]) * (1 - R[4] * R[5])) * (1 - Q[6] * Q[7] * Q[8]) *
                          R[9];
                std::cout << "S9: " << S9 << std::endl;
                return S9;
            }
            case 10: { // 10 subsystems
                auto num_sub_per_hsp = 5;
                auto num_hsp = num_subsystems / num_sub_per_hsp;

                auto x1 = (1 -
                               (1 -
                                Q[0u] * Q[1u] * R[2u]) *
                               (1 - R[3u] * R[4u]));
                auto x2 = x1 * (1 -
                               (1 -
                                Q[num_sub_per_hsp + 0u] * Q[num_sub_per_hsp + 1u] * R[num_sub_per_hsp + 2u]) *
                               (1 - R[num_sub_per_hsp + 3u] * R[num_sub_per_hsp + 4u]));
                return x2;
                //
                //
                //
                // QSystemReliability S10(system_identifier, 1);
                // for (auto i = 0u; i < num_hsp; ++i) {
                //     S10 = S10 * (1 -
                //                (1 -
                //                 Q[i * num_sub_per_hsp + 0u] * Q[i * num_sub_per_hsp + 1u] * R[i * num_sub_per_hsp + 2u]) *
                //                (1 - R[i * num_sub_per_hsp + 3u] * R[i * num_sub_per_hsp + 4u]));
                // }
                // std::cout << "S10: " << S10 << std::endl;
                // return S10;
            }
            case 11: { // 12 subsystems
                auto S11 = Q[2] * Q[4] * Q[9] * (1 - (1 - R[0] * R[3] * (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[1] * R[5] * R[11]))
                           + Q[2] * Q[4] * R[9] * Q[7] * Q[8] *
                             (1 - (1 - R[0] * R[3] * R[6] * R[10]) * (1 - R[1] * R[5] * R[11]))
                           + Q[2] * Q[4] * R[9] * Q[7] * R[8] * (1 - (1 - R[0] * R[3] * R[6]) * (1 - R[1] * R[5])) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * Q[4] * R[9] * R[7] * Q[8] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + Q[2] * Q[4] * R[9] * R[7] * R[8] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * R[4] * Q[9] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + Q[2] * R[4] * R[9] * Q[7] * Q[8] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * (1 - R[5] * R[11]))
                           + Q[2] * R[4] * R[9] * Q[7] * R[8] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[5] * Q[6]) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * R[4] * R[9] * R[7] * Q[8] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + Q[2] * R[4] * R[9] * R[7] * R[8] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[10] * Q[11])
                           + R[2] * Q[4] * Q[9] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * Q[4] * R[9] * Q[7] * Q[8] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - R[3] * R[6] * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * Q[4] * R[9] * Q[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[6]) * Q[5]) *
                             (1 - Q[10] * Q[11])
                           + R[2] * Q[4] * R[9] * R[7] * Q[8] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + R[2] * Q[4] * R[9] * R[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) * (1 - Q[10] * Q[11])
                           + R[2] * R[4] * Q[9] * (1 - Q[0] * Q[1]) * (1 - (1 - (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * R[4] * R[9] * Q[7] * Q[8] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * R[4] * R[9] * Q[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6]) * (1 - Q[10] * Q[11])
                           + R[2] * R[4] * R[9] * R[7] * Q[8] * (1 - Q[0] * Q[1]) * (1 - (1 - R[6] * R[10]) * Q[11])
                           + R[2] * R[4] * R[9] * R[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[10] * Q[11]);
                std::cout << "S11: " << S11 << std::endl;
                return S11;
            }
            case 12: { // 12 subsystems
                auto S12 = (1 - R[6] * R[7]) * (1 - R[3] * R[4]) * Q[9] *
                           (1 - (1 - R[0] * R[5] * R[10]) * (1 - R[1] * R[2] * R[8] * R[11]))
                           + (1 - R[6] * R[7]) * R[3] * R[4] * Q[9] * (1 - Q[0] * (1 - R[1] * R[2])) *
                             (1 - (1 - R[5] * R[10]) * (1 - R[8] * R[11]))
                           + (1 - R[6] * R[7]) * (1 - R[3] * R[4]) * R[9] *
                             (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2] * R[8])) * (1 - Q[10] * Q[11])
                           + (1 - R[6] * R[7]) * R[3] * R[4] * R[9] * (1 - Q[0] * (1 - R[1] * R[2])) * (1 - Q[5] * Q[8]) *
                             (1 - Q[10] * Q[11])
                           + R[6] * R[7] * (1 - R[3] * R[4]) * Q[9] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2])) *
                             (1 - (1 - R[8] * R[11]) * Q[10])
                           + R[6] * R[7] * R[3] * R[4] * Q[9] * (1 - Q[0] * (1 - R[1] * R[2])) *
                             (1 - (1 - R[8] * R[11]) * Q[10])
                           + R[6] * R[7] * (1 - R[3] * R[4]) * R[9] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2])) *
                             (1 - Q[10] * Q[11])
                           + R[6] * R[7] * R[3] * R[4] * R[9] * (1 - Q[0] * (1 - R[1] * R[2])) * (1 - Q[10] * Q[11]);
                std::cout << "S12: " << S12 << std::endl;
                return S12;
            }
            case 13: {
                auto num_blocks = 5;
                auto num_sub_per_block = num_subsystems / num_blocks;

                std::vector<QSystemReliability> QQ;
                std::vector<QSystemReliability> RR;
                QQ.reserve(num_blocks);
                RR.reserve(num_blocks);

                for (auto j = 0u; j < num_blocks; ++j) {
                    RR[j] = R[j];
                }
                for (auto j = 0u; j < num_blocks; ++j) {
                    for (auto i = 1u; i < num_sub_per_block; ++i) {
                        RR[j] =  RR[j] * R[num_blocks * i + j];
                    }
                }

                for (auto j = 0u; j < num_blocks; ++j) {
                    QQ[j] = 1 - RR[j];
                }

                auto S13 = RR[4] * (1 - QQ[0] * QQ[2]) * (1 - QQ[1] * QQ[3]) + QQ[4] * (1 - (1 - RR[0] * RR[1]) * (1 - RR[2] * RR[3]));

                std::cout << "S13: " << S13 << std::endl;
                return S13;
            }
            default:
                std::cout << "ERROR: invalid number of subsystems: " << system_identifier << std::endl;
            return QSystemReliability{};
        }
    }

    auto ProbInstance::compute_qr_system_reliability() -> QRSystemReliability {
        std::vector<QRSystemReliability> Q;
        std::vector<QRSystemReliability> R;
        Q.reserve(num_subsystems);
        R.reserve(num_subsystems);

        for (auto i = 0u; i < num_subsystems; ++i) {
            Q.emplace_back(num_subsystems, i, false, 1);
            R.emplace_back(num_subsystems, i, true, 1);
        }

        switch (system_identifier) {
            case 1: {
                auto S1 = R[4] * (1 - Q[0] * Q[2]) * (1 - Q[1] * Q[3]) +
                          Q[4] * (1 - (1 - R[0] * R[1]) * (1 - R[2] * R[3]));
                std::cout << "S1: " << S1 << std::endl;
                return S1;
            }
            case 2: {
                auto S2 = R[4] * (1 - Q[1] * Q[3]) + Q[4] * (1 - (1 - R[0] * R[1]) * (1 - R[2] * R[3]));
                std::cout << "S2: " << S2 << std::endl;
                return S2;
            }
            case 3: {
                auto S3 = R[5] * (1 - Q[0] * Q[2] * Q[3]) * (1 - Q[1] * Q[4]) +
                          Q[5] * (1 - (1 - R[0] * R[1]) * (1 - (1 - Q[2] * Q[3]) * R[4]));
                std::cout << "S3: " << S3 << std::endl;
                return S3;
            }
            case 4: {
                // auto Q234 = Q3 + R3 * (Q4 + Q5 - Q4 * Q5);
                auto Q234 = Q[2] + R[2] * Q[3] * Q[4];
                auto R234 = 1 - Q234;
                auto S4 = Q[6] * (1 - (1 - R[0] * R[1]) * (1 - R234 * R[5])) +
                          R[6] * (1 - Q[0] * Q234) * (1 - Q[1] * Q[5]);
                std::cout << "S4: " << S4 << std::endl;
                return S4;
            }
            case 5: { // 7 subsystems
                auto S5 = Q[2] * Q[4] * (1 - (1 - R[0] * R[3] * R[5]) * (1 - R[1] * R[6]))
                          + R[2] * R[4] * (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6])
                          + R[2] * Q[4] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[5]) * Q[6])
                          + Q[2] * R[4] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[5] * Q[6]);
                std::cout << "S5: " << S5 << std::endl;
                return S5;
            }
            case 6: { // 8 subsystems
                auto S6 = Q[3] * Q[5] * (1 - (1 - R[0] * R[2] * R[6]) * (1 - R[1] * R[4] * R[7]))
                          + R[3] * R[5] * (1 - Q[0] * (1 - R[1] * R[4])) * (1 - Q[6] * Q[7])
                          + R[3] * Q[5] * (1 - Q[0] * (1 - R[1] * R[4])) * (1 - (1 - R[2] * R[6]) * Q[7])
                          + Q[3] * R[5] * (1 - (1 - R[0] * R[2]) * (1 - R[1] * R[4])) * (1 - Q[6] * Q[7]);
                std::cout << "S6: " << S6 << std::endl;
                return S6;
            }
            case 7: { // 8 subsystems
                auto S7 = Q[3] * Q[4] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[6]) * (1 - R[2] * R[7]))
                          + R[3] * R[4] * (1 - Q[0] * Q[1] * Q[2]) * (1 - Q[5] * Q[6] * Q[7])
                          + R[3] * Q[4] * (1 - (1 - R[2] * R[7]) * (1 - (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6])))
                          + Q[3] * R[4] * (1 - (1 - R[0] * R[5]) * (1 - (1 - Q[1] * Q[2]) * (1 - Q[6] * Q[7])));
                std::cout << "S7: " << S7 << std::endl;
                return S7;
            }
            case 8: { // 9 subsystems
                auto S8 = Q[4] * Q[2] * Q[6] * (1 - (1 - R[0] * R[3] * R[7]) * (1 - R[1] * R[5] * R[8]))
                          + Q[4] * R[2] * Q[6] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[7]) * (1 - R[5] * R[8]))
                          + Q[4] * Q[2] * R[6] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) * (1 - Q[7] * Q[8])
                          + Q[4] * R[2] * R[6] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) * (1 - Q[7] * Q[8])
                          + R[4] * Q[2] * Q[6] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - (1 - R[5] * R[8]) * Q[7])
                          + R[4] * R[2] * Q[6] * (1 - Q[0] * Q[1]) * (1 - (1 - R[5] * R[8]) * Q[7])
                          + R[4] * Q[2] * R[6] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[7] * Q[8])
                          + R[4] * R[2] * R[6] * (1 - Q[0] * Q[1]) * (1 - Q[7] * Q[8]);
                std::cout << "S8: " << S8 << std::endl;
                return S8;
            }
            case 9: // 10 subsystems
            {
                auto S9 = (1 - (1 - (1 - Q[2] * (1 - R[0] * R[1])) * R[3]) * (1 - R[4] * R[5])) * (1 - Q[6] * Q[7] * Q[8]) *
                          R[9];
                std::cout << "S9: " << S9 << std::endl;
                return S9;
            }
            case 10: { // 10 subsystems
                auto num_sub_per_hsp = 5;
                auto num_hsp = num_subsystems / num_sub_per_hsp;
                QRSystemReliability S10(system_identifier, 1);
                for (auto i = 0u; i < num_hsp; ++i) {
                    S10 = S10 * (1 -
                               (1 -
                                Q[i * num_sub_per_hsp + 0u] * Q[i * num_sub_per_hsp + 1u] * R[i * num_sub_per_hsp + 2u]) *
                               (1 - R[i * num_sub_per_hsp + 3u] * R[i * num_sub_per_hsp + 4u]));
                }
                std::cout << "S10: " << S10 << std::endl;
                return S10;
            }
            case 11: { // 12 subsystems
                auto S11 = Q[2] * Q[4] * Q[9] * (1 - (1 - R[0] * R[3] * (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[1] * R[5] * R[11]))
                           + Q[2] * Q[4] * R[9] * Q[7] * Q[8] *
                             (1 - (1 - R[0] * R[3] * R[6] * R[10]) * (1 - R[1] * R[5] * R[11]))
                           + Q[2] * Q[4] * R[9] * Q[7] * R[8] * (1 - (1 - R[0] * R[3] * R[6]) * (1 - R[1] * R[5])) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * Q[4] * R[9] * R[7] * Q[8] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + Q[2] * Q[4] * R[9] * R[7] * R[8] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * R[4] * Q[9] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + Q[2] * R[4] * R[9] * Q[7] * Q[8] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * (1 - R[5] * R[11]))
                           + Q[2] * R[4] * R[9] * Q[7] * R[8] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[5] * Q[6]) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * R[4] * R[9] * R[7] * Q[8] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + Q[2] * R[4] * R[9] * R[7] * R[8] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[10] * Q[11])
                           + R[2] * Q[4] * Q[9] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * Q[4] * R[9] * Q[7] * Q[8] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - R[3] * R[6] * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * Q[4] * R[9] * Q[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[6]) * Q[5]) *
                             (1 - Q[10] * Q[11])
                           + R[2] * Q[4] * R[9] * R[7] * Q[8] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + R[2] * Q[4] * R[9] * R[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) * (1 - Q[10] * Q[11])
                           + R[2] * R[4] * Q[9] * (1 - Q[0] * Q[1]) * (1 - (1 - (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * R[4] * R[9] * Q[7] * Q[8] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * R[4] * R[9] * Q[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6]) * (1 - Q[10] * Q[11])
                           + R[2] * R[4] * R[9] * R[7] * Q[8] * (1 - Q[0] * Q[1]) * (1 - (1 - R[6] * R[10]) * Q[11])
                           + R[2] * R[4] * R[9] * R[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[10] * Q[11]);
                std::cout << "S11: " << S11 << std::endl;
                return S11;
            }
            case 12: { // 12 subsystems
                auto S12 = (1 - R[6] * R[7]) * (1 - R[3] * R[4]) * Q[9] *
                           (1 - (1 - R[0] * R[5] * R[10]) * (1 - R[1] * R[2] * R[8] * R[11]))
                           + (1 - R[6] * R[7]) * R[3] * R[4] * Q[9] * (1 - Q[0] * (1 - R[1] * R[2])) *
                             (1 - (1 - R[5] * R[10]) * (1 - R[8] * R[11]))
                           + (1 - R[6] * R[7]) * (1 - R[3] * R[4]) * R[9] *
                             (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2] * R[8])) * (1 - Q[10] * Q[11])
                           + (1 - R[6] * R[7]) * R[3] * R[4] * R[9] * (1 - Q[0] * (1 - R[1] * R[2])) * (1 - Q[5] * Q[8]) *
                             (1 - Q[10] * Q[11])
                           + R[6] * R[7] * (1 - R[3] * R[4]) * Q[9] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2])) *
                             (1 - (1 - R[8] * R[11]) * Q[10])
                           + R[6] * R[7] * R[3] * R[4] * Q[9] * (1 - Q[0] * (1 - R[1] * R[2])) *
                             (1 - (1 - R[8] * R[11]) * Q[10])
                           + R[6] * R[7] * (1 - R[3] * R[4]) * R[9] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2])) *
                             (1 - Q[10] * Q[11])
                           + R[6] * R[7] * R[3] * R[4] * R[9] * (1 - Q[0] * (1 - R[1] * R[2])) * (1 - Q[10] * Q[11]);
                std::cout << "S12: " << S12 << std::endl;
                return S12;
            }
            case 13: {
                auto num_blocks = 5;
                auto num_sub_per_block = num_subsystems / num_blocks;

                std::vector<QRSystemReliability> QQ;
                std::vector<QRSystemReliability> RR;
                QQ.reserve(num_blocks);
                RR.reserve(num_blocks);

                for (auto j = 0u; j < num_blocks; ++j) {
                    RR[j] = R[j];
                }
                for (auto j = 0u; j < num_blocks; ++j) {
                    for (auto i = 1u; i < num_sub_per_block; ++i) {
                        RR[j] =  RR[j] * R[num_blocks * i + j];
                    }
                }

                for (auto j = 0u; j < num_blocks; ++j) {
                    QQ[j] = 1 - RR[j];
                }

                auto S13 = RR[4] * (1 - QQ[0] * QQ[2]) * (1 - QQ[1] * QQ[3]) +
                       QQ[4] * (1 - (1 - RR[0] * RR[1]) * (1 - RR[2] * RR[3]));

                std::cout << "S13: " << S13 << std::endl;
                return S13;
            }
            default:
                std::cout << "ERROR: invalid number of subsystems: " << system_identifier << std::endl;
            return QRSystemReliability{};
        }
    }

    std::ostream &operator<<(std::ostream &out, const ProbInstance &inst)
    {
        out << "num_resource_types: " << inst.num_resource_types << "\n";
        out << "num_subsystems: " << inst.num_subsystems << "\n";
        out << "num_component_types: " << inst.num_component_types << "\n";

        out << "amount_of_each_resources: " << "\n";
        for (auto i = 0u; i < inst.num_resource_types; ++i)
        {
            out << "\t resource " << i + 1u << ": " << inst.amount_of_each_resources[i] << "\n";
        }

        out << "component_reliability_at_subsystem: " << "\n";
        for (auto j = 0u; j < inst.num_subsystems; ++j)
        {
            for (auto h = 0u; h < inst.num_component_types_at_subsystem[j]; ++h)
            {
                out << "\t component " << h + 1u << ": " << inst.component_reliability_at_subsystem[j][h] << "\t";
            }
            out << "\n";
        }

        out << "resource_for_component_at_subsystem: " << "\n";
        for (int i = 0; i < inst.num_resource_types; ++i)
        {
            out << "\t resource " << i + 1u << "\n";
            for (int j = 0; j < inst.num_subsystems; ++j)
            {
                for (int h = 0; h < inst.num_component_types_at_subsystem[j]; h++)
                {
                    out << "\t[s" << j + 1u << ", c" << h + 1u << "] = " << inst.resource_for_component_at_subsystem[i][j][h] << "\t";
                }
                out << "\n";
            }
        }

        out << "calculate_component_lb_at_subsystem: " << "\n";
        for (int j = 0; j < inst.num_subsystems; ++j)
        {
            for (int h = 0; h < inst.num_component_types_at_subsystem[j]; h++)
            {
                out << "\tL[s" << j << ", c" << h << "] = " << inst.component_lb_at_subsystem[j][h] << "\t";
            }
            out << "\n";
        }

        out << "calculate_component_ub_at_subsystem: " << "\n";
        for (int j = 0; j < inst.num_subsystems; ++j)
        {
            for (int h = 0; h < inst.num_component_types_at_subsystem[j]; h++)
            {
                out << "\tU[s" << j << ", c" << h << "] = " << inst.component_ub_at_subsystem[j][h] << "\t";
            }
            out << "\n";
        }

        return out;
    }
}