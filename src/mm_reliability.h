#pragma once
#include <cassert>
#include <algorithm>
#include <list>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace rrap {

    // Q-based system reliability

    struct QSystemReliability {
        size_t num_subsystems{};

        int const_term{};
        std::map<std::set<int>, int> m_other_terms;  // 2 * Q1Q2Q5 + ... --> { <{1,2,5}, 2>, ... }

        QSystemReliability() = default;

        explicit QSystemReliability(size_t num_subsystems, int constant = 0);

        QSystemReliability(size_t num_subsystems, int s, int coeff, int constant = 0) :
                    num_subsystems {num_subsystems}, const_term {constant} {

            std::set<int> term {s};
            m_other_terms[term] = coeff;
        }

        QSystemReliability(size_t num_subsystems, const std::vector<int>& ss, const std::vector<int>& coeffs, int constant = 0) :
                num_subsystems {num_subsystems}, const_term {constant} {

            assert(ss.size() == coeffs.size());

            for (auto s = 0; s < ss.size(); ++s) {
                std::set<int> term{ss[s]};
                m_other_terms[term] = coeffs[s];
            }
        }


        QSystemReliability(std::initializer_list<int> subsystems, int constant = 0);

        QSystemReliability(const QSystemReliability &from);

        QSystemReliability &operator=(const QSystemReliability &from);

        QSystemReliability &operator+=(const QSystemReliability &rhs);

        QSystemReliability &operator-=(const QSystemReliability &rhs);

        QSystemReliability &operator*=(int multi);

        friend QSystemReliability operator+(const QSystemReliability &left, const QSystemReliability &right);

        friend QSystemReliability operator-(const QSystemReliability &left, const QSystemReliability &right);

        friend QSystemReliability operator*(const QSystemReliability &left, const QSystemReliability &right);

        friend QSystemReliability operator+(int num, const QSystemReliability &other);

        friend QSystemReliability operator+(const QSystemReliability &other, int num);

        friend QSystemReliability operator-(int num, const QSystemReliability &other);

        friend QSystemReliability operator-(const QSystemReliability &other, int num);

        friend QSystemReliability operator*(int num, const QSystemReliability &other);

        friend QSystemReliability operator*(const QSystemReliability &other, int num);

        friend std::ostream &operator<<(std::ostream &stream, const QSystemReliability &expression);
    };



    // Q-Tree-based system reliability

    struct QTreeNode {
        size_t m_id;                                         // the order this node is added
        std::optional<size_t> m_subsystem_id;                // the (original) subsystem associated with this node

        std::shared_ptr<QTreeNode> m_father;                  // father of this node
        std::vector<std::shared_ptr<QTreeNode>> m_children;   // children of the node
        std::set<int> m_subsystems_on_route;                          // the subsystems on route from the root to current node

        bool m_logic;                                        // is this node represents a real logic subsystem ?
        int m_coeff;                                         // the coefficient of a logic subsystem

        QTreeNode(const size_t id, const std::optional<size_t> subsystem_id, std::shared_ptr<QTreeNode> father,
                  std::set<int> subsystems, bool logic_subsystem, int coeff) :
                m_id(id), m_subsystem_id(subsystem_id), m_father(std::move(father)),
                m_subsystems_on_route(std::move(subsystems)), m_logic(logic_subsystem), m_coeff(coeff) {};
    };

    struct QSystemTree {
        std::vector<std::shared_ptr<QTreeNode>> systems;

        QSystemTree() = default;
        explicit QSystemTree(const QSystemReliability &s) {
            /* create the root node */
            auto root = std::make_shared<QTreeNode>(0u, std::nullopt,
                                                    nullptr, std::set<int>(), true, s.const_term);
            systems.emplace_back(root);

            for (auto &term: s.m_other_terms) {

                std::vector<int> subsystems;
                std::copy(term.first.begin(), term.first.end(), std::back_inserter(subsystems));
                auto coeff = term.second;

                auto curr_node = root;
                auto idx = 0L;
                for (const auto &sub: subsystems) {
                    /* check if the <sub> is already exist */
                    auto iter = std::find_if(curr_node->m_children.begin(), curr_node->m_children.end(),
                                             [&sub](auto &child) {
                                                 return child->m_subsystem_id && child->m_subsystem_id.value() == sub;
                                             });

                    /* if so, we process the next sub */
                    if (iter != curr_node->m_children.end() && idx + 1u != subsystems.size()) {
                        curr_node = *iter;
                        ++idx;
                        continue;
                    }

                    /* otherwise, we create all the nodes corresponding to the subsequent subsystems */
                    for (auto i = idx; i < subsystems.size(); ++i) {
                        const auto node_id = systems.size();
                        std::set<int> current_subsystems;
                        std::copy(subsystems.begin() + i, subsystems.end(),
                                  std::inserter(current_subsystems, current_subsystems.end()));
                        auto child_node = std::make_shared<QTreeNode>(node_id, subsystems[i], curr_node,
                                                                      current_subsystems, false, 0);
                        systems.push_back(child_node);
                        curr_node->m_children.push_back(child_node);
                        curr_node = child_node;
                    }

                    // current node represents a logic subsystem
                    curr_node->m_logic = true;
                    curr_node->m_coeff = coeff;
                    break;
                }
            }
        }
    };


    // Q Sorted Tree based system reliability

    struct QSortedTreeNode {
        size_t m_id;                                         // the order this node is added
        std::optional<size_t> m_subsystem_id;                // the (original) subsystem associated with this node

        std::shared_ptr<QSortedTreeNode> m_father;                  // father of this node
        std::vector<std::shared_ptr<QSortedTreeNode>> m_children;   // children of the node
        std::set<int> m_subsystems_on_route;                          // the subsystems on route from the root to current node

        std::vector<std::pair<std::vector<int>, int>> m_allocated_logic_subsystems;

        bool m_logic;                                        // is this node represents a real logic subsystem ?
        int m_coeff;                                         // the coefficient of a logic subsystem

        QSortedTreeNode(const size_t                        id,
                        const std::optional<size_t>         subsystem_id,
                        std::shared_ptr<QSortedTreeNode>    father,
                        std::set<int>                       subsystems,
                        bool                                logic_subsystem,
                        int                                 coeff) :
                m_id(id), m_subsystem_id(subsystem_id), m_father(std::move(father)),
                m_subsystems_on_route(std::move(subsystems)), m_logic(logic_subsystem), m_coeff(coeff)
                {}
    };

    struct QSortedSystemTree {
        std::vector<std::shared_ptr<QSortedTreeNode>> systems;

        QSortedSystemTree() = default;

        explicit QSortedSystemTree(const QSystemReliability &s) {
            std::vector<std::shared_ptr<QSortedTreeNode>> remaining_nodes;

            /* create the root node */
            auto root = std::make_shared<QSortedTreeNode>(0u, std::nullopt,
                                                    nullptr, std::set<int>(), true, s.const_term);

            // copy the subsystems into a vector
            std::vector<std::pair<std::vector<int>, int>> logical_subsystems;
            for (auto &term: s.m_other_terms) {

                std::vector<int> subsystem;
                std::copy(term.first.begin(), term.first.end(), std::back_inserter(subsystem));

                logical_subsystems.emplace_back(subsystem, term.second);
            }

            root->m_allocated_logic_subsystems = logical_subsystems;

            systems.emplace_back(root);
            remaining_nodes.push_back(root);


            while (!remaining_nodes.empty()) {
                auto current_node = remaining_nodes.back();
                remaining_nodes.pop_back();

                auto current_logic_subsystem = current_node->m_allocated_logic_subsystems;

                while (!current_logic_subsystem.empty()) {
                    /* sort the subsystems according to the number of their occurrence in logic subsystems */
                    std::map<int, int> num_occurrence;
                    for (auto &subsystems: current_logic_subsystem) {
                        for (auto &term: subsystems.first) {
                            num_occurrence[term] += 1;
                        }
                    }
#if 0
                    for (int x = 0; x < s.num_subsystems; ++x) {
                        std::cout << "num_occurrence[" << x << "]=" << num_occurrence[x] << std::endl;
                    }
#endif
                    auto pos = std::max_element(num_occurrence.begin(), num_occurrence.end(),
                                                [](const auto &a, const auto &b) { return a.second < b.second; });

                    auto max = pos->first;

#ifndef NDEBUG
                    std::cout << "max = " << max << "\n";
#endif
                    // split the logical subsystems into two groups
                    // - the first one contains the subsystems containing <max>
                    // - the other one contains the remaining subsystems

                    std::vector<std::pair<std::vector<int>, int>> allocated_logical_subsystems;
                    std::vector<std::pair<std::vector<int>, int>> remaining_logical_subsystems;
                    auto is_logic = false;
                    auto logic_coeff = 0;
                    for (auto &[subsystems, coeff]: current_logic_subsystem) {
                        if (std::find(subsystems.begin(), subsystems.end(), max) != subsystems.end()) {
                            auto ss = subsystems;
                            ss.erase(std::remove(ss.begin(), ss.end(), max), ss.end());

                            if (ss.empty()) {
                                // this represents a logic subsystem
                                is_logic = true;
                                logic_coeff = coeff;
                            } else {
                                allocated_logical_subsystems.emplace_back(ss, coeff);
                            }
                        } else {
                            remaining_logical_subsystems.emplace_back(subsystems, coeff);
                        }
                    }
#if 0
                    std::cout << "allocated_logical_subsystems: \n";
                    for (auto& [c, coeff]: allocated_logical_subsystems) {
                        for (auto cc: c) {
                            std::cout << cc << " ";
                        }
                        std::cout << "\n";
                    }

                    std::cout << "remaining_logical_subsystems: \n";
                    for (auto& [c, coeff]: remaining_logical_subsystems) {
                        for (auto cc: c) {
                            std::cout << cc << " ";
                        }
                        std::cout << "\n";
                    }
#endif
                    // create a new child node
                    const auto node_id = systems.size();
                    auto child_node = std::make_shared<QSortedTreeNode>(node_id, max, current_node,
                                                                        std::set<int>(), is_logic, logic_coeff);

                    child_node->m_allocated_logic_subsystems = allocated_logical_subsystems;

                    current_node->m_children.push_back(child_node);

                    systems.push_back(child_node);

                    if (!allocated_logical_subsystems.empty()) {
                        // the child node need to be further processed
                        remaining_nodes.push_back(child_node);
                    }

                    current_logic_subsystem = remaining_logical_subsystems;
                }
            }
#if 0
            for (auto & system : systems) {
                std::cout << "m_id: " << system->m_id << "\n";
                std::cout << "m_father: ";
                if (system->m_father) std::cout << system->m_father->m_id;
                std::cout << "\n";
                std::cout << "m_logic: " << system->m_logic << "\n";
                std::cout << "m_subsystem_id: ";
                if(system->m_subsystem_id.has_value()) { std::cout << system->m_subsystem_id.value(); }
                std::cout << "\n";
                std::cout << "m_coeff: " << system->m_coeff << "\n";
                std::cout << "m_subsystems_on_route: \n";
                //  for (auto c: system->m_subsystems_on_route) {
                //      std::cout << c << " ";
                //  }
                std::cout << "m_subsystems_on_route: \n";
                for (auto& [c, coeff]: system->m_allocated_logic_subsystems) {
                    for (auto cc: c) {
                        std::cout << cc << " ";
                    }
                    std::cout << "\n";
                }
                std::cout << "m_children: ";
                for (const auto &c: system->m_children) {
                    std::cout << c->m_id << "[" << c->m_subsystem_id.value_or(-1) << "] ";
                }
                std::cout << "\n" << std::endl;
            }
#endif
        }
    };



    // QR-based system reliability

    struct QRSystemReliability {
        size_t num_subsystems{};

        int const_term{};

        std::map<std::set<std::pair<int, bool>>, int> m_other_terms;  // 2 * Q1Q2Q5 + ... --> { <{1,2,5}, 2>, ... }

        QRSystemReliability() = default;

        explicit QRSystemReliability(size_t num_subsystems, int constant = 0);

        QRSystemReliability(size_t num_subsystems, int s, bool qr, int coeff, int constant = 0);

        QRSystemReliability(size_t num_subsystems, const std::vector<int>& ss, const std::vector<bool>& qrs, const std::vector<int>& coeffs, int constant = 0);

        QRSystemReliability(std::initializer_list<std::pair<int,bool>> subsystems, int constant = 0);

        QRSystemReliability(const QRSystemReliability &from);

        QRSystemReliability &operator=(const QRSystemReliability &from);

        QRSystemReliability &operator+=(const QRSystemReliability &rhs);

        QRSystemReliability &operator-=(const QRSystemReliability &rhs);

        QRSystemReliability &operator*=(int multi);


        friend QRSystemReliability operator+(const QRSystemReliability &left, const QRSystemReliability &right);

        friend QRSystemReliability operator-(const QRSystemReliability &left, const QRSystemReliability &right);

        friend QRSystemReliability operator*(const QRSystemReliability &left, const QRSystemReliability &right);

        friend QRSystemReliability operator+(int num, const QRSystemReliability &other);

        friend QRSystemReliability operator+(const QRSystemReliability &other, int num);

        friend QRSystemReliability operator-(int num, const QRSystemReliability &other);

        friend QRSystemReliability operator-(const QRSystemReliability &other, int num);

        friend QRSystemReliability operator*(int num, const QRSystemReliability &other);

        friend QRSystemReliability operator*(const QRSystemReliability &other, int num);

        friend std::ostream &operator<<(std::ostream &stream, const QRSystemReliability &expression);
    };



    // QR Tree Structure

    struct QRTreeNode {
        size_t m_id;                                           // the order this node is added
        std::optional<size_t> m_subsystem_id;                  // the (original) subsystem associated with this node
                                                               // - root node has no associated subsystem
        bool m_qr;                                             // is this an R node? false for Q node

        std::shared_ptr<QRTreeNode> m_father;                  // father of this node
        std::vector<std::shared_ptr<QRTreeNode>> m_children;   // children of the node
        std::set<std::pair<int,bool>> m_subsystems_on_route;   // the subsystems on route from the root to current node

        bool m_logic;                                          // is this node represents a real logic subsystem ?
        int m_coeff;                                           // the coefficient of a logic subsystem

        // bool m_root;

        [[nodiscard]] bool is_root() const { return m_father == nullptr; }
        [[nodiscard]] bool is_logic() const { return m_logic; }
        [[nodiscard]] bool is_q() const { return m_qr; }
        [[nodiscard]] bool is_r() const { return !is_q(); }
        [[nodiscard]] std::shared_ptr<QRTreeNode> get_father() const { return m_father; }
        [[nodiscard]] std::vector<std::shared_ptr<QRTreeNode>> get_children() const { return m_children; }
        [[nodiscard]] int get_coeff() const {return m_coeff; }
        [[nodiscard]] size_t get_subsystem() const { return m_subsystem_id.value(); }
        [[nodiscard]] bool has_child() const { return !m_children.empty(); }

        QRTreeNode(const size_t id,
                   const std::optional<size_t> subsystem_id,
                   bool qr,
                   std::shared_ptr<QRTreeNode> father,
                   std::set<std::pair<int,bool>> subsystems,
                   bool logic_subsystem,
                   int coeff) :
                m_id(id), m_subsystem_id(subsystem_id), m_qr(qr), m_father(std::move(father)),
                m_subsystems_on_route(std::move(subsystems)), m_logic(logic_subsystem), m_coeff(coeff) {};
    };


    struct QRSystemTree {
        std::vector<std::shared_ptr<QRTreeNode>> systems;

        QRSystemTree() = default;
        explicit QRSystemTree(const QRSystemReliability &s) {
            /* create the root node */
            auto root = std::make_shared<QRTreeNode>(0u, std::nullopt, false,
                                                    nullptr, std::set<std::pair<int,bool>>(), true, s.const_term);
            systems.emplace_back(root);

            for (auto &term: s.m_other_terms) {

                std::vector<std::pair<int,bool>> subsystems;
                std::copy(term.first.begin(), term.first.end(), std::back_inserter(subsystems));
                auto coeff = term.second;

                auto curr_node = root;
                auto idx = 0L;
                for (const auto &sub: subsystems) {
                    /* check if the <sub> is already exist */
                    auto iter = std::find_if(curr_node->m_children.begin(), curr_node->m_children.end(),
                                             [&sub](auto &child) {
                                                 return child->m_subsystem_id && child->m_qr == sub.second && child->m_subsystem_id.value() == sub.first;
                                             });

                    /* if so, we process the next sub */
                    if (iter != curr_node->m_children.end() && idx + 1u != subsystems.size()) {
                        curr_node = *iter;
                        ++idx;
                        continue;
                    }

                    /* otherwise, we create all the nodes corresponding to the subsequent subsystems */
                    for (auto i = idx; i < subsystems.size(); ++i) {
                        const auto node_id = systems.size();
                        std::set<std::pair<int, bool>> current_subsystems;
                        std::copy(subsystems.begin(), subsystems.begin() + i,
                                  std::inserter(current_subsystems, current_subsystems.end()));


                        auto child_node = std::make_shared<QRTreeNode>(node_id, subsystems[i].first, subsystems[i].second, curr_node,
                                                                       current_subsystems, false, 0);
                        systems.push_back(child_node);
                        curr_node->m_children.push_back(child_node);
                        curr_node = child_node;
                    }

                    // current node represents a logic subsystem
                    curr_node->m_logic = true;
                    curr_node->m_coeff = coeff;
                    break;
                }
            }

#if 0
            for (auto & system : systems) {
                std::cout << "m_id: " << system->m_id << "\n";
                std::cout << "m_father: ";
                if (system->m_father) std::cout << system->m_father->m_id;
                std::cout << "\n";
                std::cout << "m_logic: " << system->m_logic << "\n";
                std::cout << "m_subsystem_id: ";
                if(system->m_subsystem_id.has_value()) {
                    if (system->m_qr) { std::cout << "R_"; }
                    else { std::cout << "Q_"; }
                    std::cout << system->m_subsystem_id.value();
                }
                std::cout << "\n";
                std::cout << "m_coeff: " << system->m_coeff << "\n";
                std::cout << "m_subsystems_on_route: \n";
                //  for (auto c: system->m_subsystems_on_route) {
                //      std::cout << c << " ";
                //  }
                std::cout << "m_children: ";
                for (const auto &c: system->m_children) {
                    std::cout << c->m_id << "[" << (c->m_qr ? "R" : "Q") << c->m_subsystem_id.value_or(-1) << "] ";
                }
                std::cout << "\n" << std::endl;
            }
#endif
        }
    };


    // QR Sorted Tree Structure

    struct QRSortedTreeNode {
        size_t m_id;                                         // the order this node is added
        std::optional<size_t> m_subsystem_id;                // the (original) subsystem associated with this node
        bool m_qr;                                           // is this an R node? false for Q node

        std::shared_ptr<QRSortedTreeNode> m_father;                  // father of this node
        std::vector<std::shared_ptr<QRSortedTreeNode>> m_children;   // children of the node
        std::set<int> m_subsystems_on_route;                          // the subsystems on route from the root to current node

        std::vector<std::pair<std::vector<std::pair<int, bool>>, int>> m_allocated_logic_subsystems;

        bool m_logic;                                        // is this node represents a real logic subsystem ?
        int m_coeff;                                         // the coefficient of a logic subsystem

        QRSortedTreeNode(const size_t                        id,
                         const std::optional<size_t>         subsystem_id,
                         bool                                qr,
                         std::shared_ptr<QRSortedTreeNode>   father,
                         std::set<int>                       subsystems,
                         bool                                logic_subsystem,
                         int                                 coeff) :
                m_id(id), m_subsystem_id(subsystem_id), m_qr(qr), m_father(std::move(father)),
                m_subsystems_on_route(std::move(subsystems)), m_logic(logic_subsystem), m_coeff(coeff)
        {}
    };



    struct QRSortedSystemTree {
        std::vector<std::shared_ptr<QRSortedTreeNode>> systems;

        QRSortedSystemTree() = default;

        [[nodiscard]] size_t get_num_nodes() const { return systems.size(); }

        explicit QRSortedSystemTree(const QRSystemReliability &s) {
            std::vector<std::shared_ptr<QRSortedTreeNode>> remaining_nodes;

            /* create the root node */
            auto root = std::make_shared<QRSortedTreeNode>(0u, std::nullopt, false,
                                                          nullptr, std::set<int>(), true, s.const_term);

            // copy the subsystems into a vector
            std::vector<std::pair<std::vector<std::pair<int, bool>>, int>> logical_subsystems;
            for (auto &term: s.m_other_terms) {

                std::vector<std::pair<int, bool>> subsystem;
                std::copy(term.first.begin(), term.first.end(), std::back_inserter(subsystem));

                logical_subsystems.emplace_back(subsystem, term.second);
            }

            root->m_allocated_logic_subsystems = logical_subsystems;

            systems.push_back(root);
            remaining_nodes.push_back(root);


            while (!remaining_nodes.empty()) {
                auto current_node = remaining_nodes.back();
                remaining_nodes.pop_back();

                auto current_logic_subsystem = current_node->m_allocated_logic_subsystems;

                while (!current_logic_subsystem.empty()) {
                    /* sort the subsystems according to the number of their occurrence in logic subsystems */
                    std::map<std::pair<int, bool>, int> num_occurrence;
                    for (auto &subsystems: current_logic_subsystem) {
                        for (auto &term: subsystems.first) {
                            num_occurrence[term] += 1;
                        }
                    }
#if 0
                    for (int x = 0; x < s.num_subsystems; ++x) {
                        std::cout << "num_occurrence[" << x << "]=" << num_occurrence[x] << std::endl;
                    }
#endif
                    auto pos = std::max_element(num_occurrence.begin(), num_occurrence.end(),
                                                [](const auto &a, const auto &b) { return a.second < b.second; });

                    auto max = pos->first;
#ifndef NDEBUG
                    std::cout << "max = " << (max.second ? "R" : "Q") << max.first << "\n";
#endif
                    // split the logical subsystems into two groups
                    // - the first one contains the subsystems containing <max>
                    // - the other one contains the remaining subsystems

                    std::vector<std::pair<std::vector<std::pair<int,bool>>, int>> allocated_logical_subsystems;
                    std::vector<std::pair<std::vector<std::pair<int,bool>>, int>> remaining_logical_subsystems;
                    auto is_logic = false;
                    auto logic_coeff = 0;
                    for (auto &[subsystems, coeff]: current_logic_subsystem) {
                        if (std::find(subsystems.begin(), subsystems.end(), max) != subsystems.end()) {
                            auto ss = subsystems;
                            ss.erase(std::remove(ss.begin(), ss.end(), max), ss.end());

                            if (ss.empty()) {
                                // this represents a logic subsystem
                                is_logic = true;
                                logic_coeff = coeff;
                            } else {
                                allocated_logical_subsystems.emplace_back(ss, coeff);
                            }
                        } else {
                            remaining_logical_subsystems.emplace_back(subsystems, coeff);
                        }
                    }
#if 0
                    std::cout << "allocated_logical_subsystems: \n";
                    for (auto& [c, coeff]: allocated_logical_subsystems) {
                        for (auto cc: c) {
                            std::cout << cc << " ";
                        }
                        std::cout << "\n";
                    }

                    std::cout << "remaining_logical_subsystems: \n";
                    for (auto& [c, coeff]: remaining_logical_subsystems) {
                        for (auto cc: c) {
                            std::cout << cc << " ";
                        }
                        std::cout << "\n";
                    }
#endif
                    // create a new child node
                    const auto node_id = systems.size();
                    auto child_node = std::make_shared<QRSortedTreeNode>(node_id, max.first, max.second, current_node,
                                                                        std::set<int>(), is_logic, logic_coeff);

                    child_node->m_allocated_logic_subsystems = allocated_logical_subsystems;

                    current_node->m_children.push_back(child_node);

                    systems.push_back(child_node);

                    if (!allocated_logical_subsystems.empty()) {
                        // the child node need to be further processed
                        remaining_nodes.push_back(child_node);
                    }

                    current_logic_subsystem = remaining_logical_subsystems;
                }
            }
#if 0
            for (auto & system : systems) {
                std::cout << "m_id: " << system->m_id << "\n";
                std::cout << "m_father: ";
                if (system->m_father) std::cout << system->m_father->m_id;
                std::cout << "\n";
                std::cout << "m_logic: " << system->m_logic << "\n";
                std::cout << "m_subsystem_id: ";
                if(system->m_subsystem_id.has_value()) { std::cout << system->m_subsystem_id.value(); }
                std::cout << "\n";
                std::cout << "m_coeff: " << system->m_coeff << "\n";
                std::cout << "m_subsystems_on_route: \n";
                //  for (auto c: system->m_subsystems_on_route) {
                //      std::cout << c << " ";
                //  }
                std::cout << "m_subsystems_on_route: \n";
                for (auto& [c, coeff]: system->m_allocated_logic_subsystems) {
                    for (auto cc: c) {
                        std::cout << (cc.second ? "R" : "Q") << cc.first << " ";
                    }
                    std::cout << "\n";
                }
                std::cout << "m_children: ";
                for (const auto &c: system->m_children) {
                    std::cout << c->m_id << "[" << (c->m_qr ? "R" : "Q") << c->m_subsystem_id.value_or(-1) << "] ";
                }
                std::cout << "\n" << std::endl;
            }
#endif
        }
    };

}