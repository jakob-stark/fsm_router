#include <map>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string_view>

#include <fmt/core.h>
#include <argparse/argparse.hpp>

#include <fsm2.hpp>

using namespace std::literals;
using fmt::print;

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

    fsm2::FSM fsm = fsm2::from_string("abcde");

    std::cout << fsm.to_graphviz();


    return 0;

}
