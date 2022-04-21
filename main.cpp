#include <map>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string_view>

#include <argparse/argparse.hpp>

#include <router.hpp>

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

    return 0;

}
