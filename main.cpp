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

    eNFA enfa3 {
        {
            { {1, {'0'}}, {2} },
            { {1, {'1'}}, {3} },
            { {2, {'0'}}, {1} },
            { {2, {'1'}}, {4} },
            { {3, {'0'}}, {5} },
            { {3, {'1'}}, {6} },
            { {4, {'0'}}, {5} },
            { {4, {'1'}}, {6} },
            { {5, {'0'}}, {5} },
            { {5, {'1'}}, {6} },
            { {6, {'0'}}, {6} },
            { {6, {'1'}}, {6} },
        },
        {1},
        {3,4,5}
    };

    eNFA enfa4 {
        std::string("hallo welt")
    };

    eNFA enfa5 {
        std::string("xyz")
    };

    eNFA enfa6 = enfa_number | enfa4 | enfa5 | eNFA("xzc"sv);

    enfa6.powerset().print_dot();

    //DFA dfa = nfa.powerset().reverse().powerset().reverse().powerset();
    //fmt::print("\n");

    //dfa.print_dot();

    return 0;

}
