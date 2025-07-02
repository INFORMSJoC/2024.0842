#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <cassert>
#include "cxxopts.hpp"
#include "bb_instance.h"
#include "bb_solver.h"

/*
 * This project implements a Branch-and-Bound algorithm for solving the Reliability Allocation
 * Problem (RAP), adapted from the work of Ha and Kuo (2006).
 * 1. Transform a subsystem with H types of components in a heterogeneous system into a series of
 * H parallel subsystems, each configured with exactly one type of components.
 * 2. Set the lower bounds to zero and enumerate all possible combinations: for each (original)
 * subsystem, we can select any of the available component types and include at least one component
 * of that type. Then, for each combination, we run the branch-and-bound algorithm. The number of
 * possible combinations can be large.
 */

int main(int argc, char *argv[]) {
    cxxopts::Options options("ra_bb", "Branch-and-bound algorithm for RAP, adapted from Ha and Kuo (2006).");
    options.add_options()("h,help", "display usage information")
            ("s,system", "system identifier associated with the instance: 1, 2,..., 12", cxxopts::value<int>()->default_value("0"))
            ("f,file", "path to the instance data file", cxxopts::value<std::string>())
            ("x,time-limit", "time limit in seconds, 0 for unlimited", cxxopts::value<double>()->default_value("0"))
            ("t,test", "test mode: 0, 1, 2, 3, 4", cxxopts::value<int>()->default_value("0"))
            ("v,version", "print the version", cxxopts::value<std::string>());

    /* sparse the parameters */
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return EXIT_SUCCESS;
    }

    if (result.count("version")) {
        std::cout << "version: 1.0" << std::endl;
        return EXIT_SUCCESS;
    }

    // if (!result.count("test")) {
    //     std::cout << "ERROR: instance type must be provided." << std::endl;
    //     std::cout << options.help() << std::endl;
    //     return EXIT_FAILURE;
    // }
    auto test_type = result["test"].as<int>();

    auto system_identifier = result["system"].as<int>();

    if (!result.count("file")) {
        std::cout << "ERROR: data file of instance must be provided." << std::endl;
        std::cout << options.help() << std::endl;
        return EXIT_FAILURE;
    }
    auto data_file = result["file"].as<std::string>();
    auto time_limit = result["time-limit"].as<double>();
    // std::cout << "time limit is: " << time_limit << " seconds" << std::endl;

    // read instance
    ProbInstance problem(test_type, system_identifier, data_file);
#ifndef NDEBUG
    std::cout << problem << std::endl;
#endif

    if (problem.test_type == TEST_TYPE::TEST_ONE || problem.test_type == TEST_TYPE::TEST_TWO) {
        assert(problem.num_component_types == 1);
        assert(problem.num_resource_types == 3);
    }

    // solve
    BabSolver solver(problem);
    const auto start = std::chrono::high_resolution_clock::now();
    solver.branch_and_bound(start, time_limit);
    const auto end = std::chrono::high_resolution_clock::now();

    const auto used_time = std::chrono::duration<double>(end - start).count();

    // output result
    std::cout << "=================\n";
    std::cout << "Optimal Solution:\t";
    const auto solution = solver.get_best_solution();
    for (const auto m: solution) {
        std::cout << m << "\t";
    }
    std::cout << "\nOptimal Solution Value: " << solver.get_best_solution_value() << std::endl;

    std::cout << "Explored nodes: " << solver.get_node_count() << std::endl;
    std::cout << "used_time: " << used_time << std::endl;


    // output the solution
    const std::string method_name("BB");
    std::filesystem::path out_dir(data_file);
    auto instance_name = out_dir.stem().string();
    auto out_file_name = "_" + out_dir.stem().string() + "_" + method_name + "_sol.txt";

    if (problem.test_type == TEST_TYPE::MIXED) {
        out_file_name.insert(0, std::string("s") + std::to_string(system_identifier));
    }
    else if (problem.test_type == TEST_TYPE::EXAMPLE_ONE) {
        out_file_name.insert(0, std::string("e1"));
    }
    else if (problem.test_type == TEST_TYPE::EXAMPLE_TWO) {
        out_file_name.insert(0, std::string("e2"));
    }
    else if (problem.test_type == TEST_TYPE::TEST_ONE) {
        out_file_name.insert(0, std::string("t1"));
    }
    else if (problem.test_type == TEST_TYPE::TEST_TWO) {
        out_file_name.insert(0, std::string("t2"));
    }
    else {
        std::cerr << "error test type..." << std::endl;
    }

    out_dir.replace_filename(out_file_name);
    std::cout << "--> output directory: " << out_dir.string() << std::endl;
    std::ofstream ofs(out_dir.string(), std::ios::out);

    ofs << "instance_name" << ","
        << "method_name" << ","
        << "used_time" << ","
        << "initial_node_count" << ","
        << "node_count" << ","
        << "remained_node" << ","
        << "time_limit" << ","
        << "best_solution_value" << ","
        << "best_solution" << "\n";
    ofs << instance_name << ","
        << method_name << ","
        << used_time << ","
        << solver.get_initial_node_count() << ","
        << solver.get_node_count() << ","
        << solver.get_remained_node_count() << ","
        << (solver.is_time_limit()? 1 : 0) << ","
        << solver.get_best_solution_value() << ",";

    for (const auto m: solution) {
        ofs << m << "\t";
    }
    ofs << "\n";
    ofs.close();

    return EXIT_SUCCESS;
}
