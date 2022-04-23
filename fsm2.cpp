#include <ranges>
#include <algorithm>
#include <type_traits>
#include <iostream>

#include <fsm2.hpp>


namespace fsm2 {

using namespace std::literals;
namespace stdr = std::ranges;
namespace stdv = std::views;

constexpr inline auto operator"" _z ( unsigned long long int literal ) {
    return std::size_t(literal);
}

void FSM::remap(FSM& fsm, FSM::index_t state_offset) {
    for ( auto&& [symbol, target] : stdv::join(fsm.states) ) {
        if ( symbol != FSM::leaving ) {
            target.state += state_offset;
        }
    }
}

void FSM::relink(FSM& fsm, FSM::index_t leave_to) {
    for ( auto&& state : fsm.states ) {
        while ( auto node_handle = state.extract(FSM::leaving) ) {
            node_handle.key() = FSM::epsilon;
            node_handle.mapped().state = leave_to;
            state.insert(std::move(node_handle));
        }
    }
}

std::string FSM::to_graphviz() const {
    std::string result;
    result += "digraph fsm {\n"s                                ;
    result += "rankdir=LR;\n"s                                  ;
    for ( index_t i : stdv::iota(0_z,states.size()) ) {
        result += std::to_string(i) + " ["s                     ;
        if ( states[i].contains(leaving) ) {
            result += "shape = doublecircle, "s                 ;
        } else {
            result += "shape = circle "s                        ;
        }
        if ( i == 0_z ) {
            result += "style = filled"s                         ;
        } else {
            result += "style = \"\""s                           ;
        }
        result += "];\n"s                                       ;
    }
            
    index_t i {0_z};
    for ( const auto& state : states ) {
        for ( const auto& [symbol, target] : state) {
            if ( symbol != leaving ) {
                result += std::to_string(i) + " -> "s;
                result += std::to_string(target.state) + " [label=\""s;
                if ( symbol == epsilon ) {
                    result += "\u03f5"s;
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

FSM::state_set_t FSM::epsilon_closure(const FSM::state_set_t& input_set) const {
    state_set_t closure {input_set};
    std::queue<index_t> Q {{closure.begin(), closure.end()}};

    while ( !Q.empty() ) {
        auto&& state = states[Q.front()]; Q.pop();

        auto&& [ebegin, eend] = state.equal_range(epsilon);
        for ( auto&& [_, target] : stdr::subrange(ebegin, eend) ) {
            if ( closure.insert(target.state).second ) {
                Q.push(target.state);
            }
        }
    }

    return closure;
}

FSM FSM::powerset() const {
    FSM fsm {};

    if ( states.empty() ) {
        return fsm;
    }

    fsm.states.reserve(states.size() * 2);

    std::map<state_set_t, index_t> smap {{epsilon_closure({0_z}), 0_z}};
    std::queue<stdr::iterator_t<decltype(smap)>> Q;
    Q.push(smap.begin());
    fsm.states.push_back({});

    while ( !Q.empty() ) {
        auto&& [pq, pi] = *Q.front(); Q.pop();
        std::map<symbol_t, std::pair<state_set_t,std::vector<action_t>>> pt {};

        for ( index_t q : pq ) {
            for ( auto&& [symbol, target] : states[q] ) {
                if ( symbol != epsilon ) {
                    auto&& [state_set, action_list] {pt[symbol]};
                    state_set.insert(target.state);
                    stdr::copy(target.actions, std::back_inserter(action_list));
                }
            }
        }

        for ( auto&& [s, pp] : pt ) {
            auto&& [map_entry, inserted] = smap.insert(
                    {epsilon_closure(pp-), fsm.states.size()});
            if ( inserted ) {
                Q.push(map_entry);
                fsm.states.push_back({});
            }

            fsm.states[pi].insert({s, {map_entry->second,{}}});
        }
    }
    
    fsm.states.shrink_to_fit(); 
    return fsm;
}


FSM FSM::reverse() const {
    FSM fsm {};

    if ( states.empty() ) {
        return fsm;
    }
    
    fsm.states = std::vector<state_t> {states.size() + 1_z, state_t{}};

    index_t si {1_z};
    for ( auto&& state : states ) {
        for ( auto&& [symbol, target] : state ) {
            if ( symbol != leaving ) {
                fsm.states[target.state+1_z].insert(
                        {symbol, {si, target.actions}});
            } else {
                fsm.states[0_z].insert(
                        {epsilon, {si, target.actions}});
            }
        }
        si++;
    }

    fsm.states[1_z].insert({leaving, {0_z, {}}});

    fsm.states.shrink_to_fit();
    return fsm;
}

FSM FSM::brzozowski() const {
    return this->reverse().powerset().reverse().powerset();
}

bool FSM::match(const std::string& str) const {
    if ( states.empty() ) {
        return false;
    }

    index_t cs {0_z};
    std::size_t pos {0_z};
    for ( FSM::symbol_t c : str ) {
        if ( auto&& t { states[cs].find(c) }; t == states[cs].end() ) {
            return false;
        } else {
            cs = t->second.state;
        }
    }

    return states[cs].contains(leaving);
}

FSM from_string(const std::string& str) {
    FSM fsm {};

    fsm.states.reserve(str.size() + 1_z);
    
    fsm.states.push_back({});
    for ( FSM::symbol_t c : str ) {
        fsm.states.back().insert({c, {fsm.states.size(), {}}});
        fsm.states.push_back({});
    }
    fsm.states.back().insert({FSM::leaving, {0_z, {}}});
    return fsm;
}

FSM concatenate(FSM lhs, FSM rhs) {
    FSM fsm {};

    const auto lss {lhs.states.size()};
    const auto rss {rhs.states.size()};
    
    fsm.states.reserve( lss + rss );

    FSM::relink(lhs, lss);
    FSM::remap(rhs, lss);

    stdr::move(lhs.states, std::back_inserter(fsm.states));
    stdr::move(rhs.states, std::back_inserter(fsm.states));

    return fsm;
}

FSM alternative(FSM lhs, FSM rhs) {
    FSM fsm {};

    const auto lss {lhs.states.size()};
    const auto rss {rhs.states.size()};

    fsm.states.reserve( lss + rss + 2_z );

    FSM::remap(lhs, 1_z);
    FSM::remap(rhs, lss + 1_z);
    FSM::relink(lhs, lss + rss + 1_z);
    FSM::relink(rhs, lss + rss + 1_z);

    fsm.states.push_back({{
            {FSM::epsilon, {1_z      , {}}},
            {FSM::epsilon, {lss + 1_z, {}}}
    }});

    stdr::move(lhs.states, std::back_inserter(fsm.states));
    stdr::move(rhs.states, std::back_inserter(fsm.states));

    fsm.states.push_back({{{FSM::leaving, {}}}});
    return fsm;
}

FSM kleene_star(FSM lhs) {
    FSM fsm {};

    const auto lss {lhs.states.size()};

    FSM::remap(lhs, 1_z);
    FSM::relink(lhs, lss + 1_z);

    fsm.states.push_back({{
            {FSM::epsilon, {1_z      , {}}},
            {FSM::epsilon, {lss + 1_z, {}}}
    }});

    stdr::move(lhs.states, std::back_inserter(fsm.states));

    fsm.states.push_back({{
            {FSM::epsilon, {1_z, {}}},
            {FSM::leaving, {{},  {}}}
    }});

    return fsm;
}

}
