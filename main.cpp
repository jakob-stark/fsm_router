#include <map>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string_view>

#include <argparse/argparse.hpp>
#include <fmt/core.h>

#include <fsm.hpp>

using namespace std::literals;

int main(int argc, char* argv[]) {
    
    argparse::ArgumentParser parser("url_test");
    parser.add_argument("url")
        .help("url to parse");

    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    //fmt::print("hallo {}\n", parser.get<std::string>("url"));

    eNFA enfa {
        {
            { {0, {'a'}}, {1} },
            { {1, {'b'}}, {2,3} },
            { {2, {   }}, {0} }
        },
        {0},
        {3}
    };

    NFA nfa = enfa.to_nfa();
    nfa.print_dot();

    /*
    for ( auto [k, v] : enfa.T ) {
        std::cout << k.first << "(" << k.second.value_or('#') << ") -> {";
        for ( auto e : v ) {
            std::cout << e << ", ";
        }
        std::cout << "}\n";
    }

    std::cout << '\n';

    std::cout << "E({0}) = {";
    for ( auto e : enfa.E({0}) ) {
        std::cout << e << ", ";
    }
    std::cout << "}\n";

    std::cout << "E({1}) = {";
    for ( auto e : enfa.E({1}) ) {
        std::cout << e << ", ";
    }
    std::cout << "}\n";

    std::cout << "E({}) = {";
    for ( auto e : enfa.E({}) ) {
        std::cout << e << ", ";
    }
    std::cout << "}\n";

    std::cout << "E({0,1}) = {";
    for ( auto e : enfa.E({0,1}) ) {
        std::cout << e << ", ";
    }
    std::cout << "}\n";

    std::cout << '\n'; 
    for ( auto [k, v] : nfa.T ) {
        std::cout << k.first << "(" << k.second << ") -> {";
        for ( auto e : v ) {
            std::cout << e << ", ";
        }
        std::cout << "}\n";
    }*/

    return 0;

}
