#include <algorithm>

#include <fsm2.hpp>

namespace fsm2 {

using namespace std::literals;
namespace stdr = std::ranges;


std::string FSM::to_graphviz() const {
    std::string result;
    result += "digraph fsm {\n"s                                ;
    result += "rankdir=LR;\n"s                                  ;
    for ( index_t i : stdr::iota_view{0ul,states.size()} ) {
        result += std::to_string(i) + " ["s                     ;
        if ( states[i].transitions.contains(leaving) ) {
            result += "shape = doublecircle, "s                 ;
        } else {
            result += "shape = circle "s                        ;
        }
        if ( i == 0ul ) {
            result += "style = filled"s                         ;
        } else {
            result += "style = \"\""s                           ;
        }
        result += "];\n"s                                       ;
    }
            
    index_t i {0};
    for ( const auto& state : states ) {
        for ( const auto& [symbol, target] : state.transitions) {
            if ( symbol != leaving ) {
                result += std::to_string(i) + " -> "s;
                result += std::to_string(target.state) + " [label=\""s;
                if ( symbol == epsilon ) {
                    result += "\u03b5"s;
                } else {
                    result += symbol;
                }
                result += "\"];\n"s;
            }
        }
        i++;
    }       
    result += "}\n"s                                          ;
    return result;

}

FSM powerset(const FSM& fsm) {

}

FSM from_string(const std::string& str) {
    FSM fsm {};

    FSM::index_t i {0};
    for ( FSM::symbol_t c : str ) {
        fsm.states.push_back({{{c, {++i, {}}}}});
    }
    fsm.states.push_back({{{FSM::leaving, {0ul, {}}}}});
    return fsm;
}

FSM concatenate(const FSM& lhs, const FSM& rhs) {

}

}
