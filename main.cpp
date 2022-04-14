#include <map>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string_view>

#include <argparse/argparse.hpp>

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

    eNFA enfa2 {
        {
            { {1, {'0'}}, {2} },
            { {1, {   }}, {3} },
            { {2, {'1'}}, {2,4} },
            { {3, {'0'}}, {4} },
            { {3, {   }}, {2} },
            { {4, {'0'}}, {3} }
        },
        {1},
        {3,4}
    };

    NFA nfa = enfa2.to_nfa();
    //nfa.print_dot();

    fmt::print("\n");
    DFA dfa = nfa.powerset();
    fmt::print("\n");

    dfa.print_dot();

    return 0;

}
