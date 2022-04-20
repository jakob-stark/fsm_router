#include <map>
#include <set>
#include <queue>
#include <string>

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <fmt/format.h>

using namespace std::literals;
namespace stdr = std::ranges;

template<typename state_t = unsigned int, typename symbol_t = char>
class DFA;

template<typename state_t = unsigned int, typename symbol_t = char>
class NFA;

template<typename state_t = unsigned int, typename symbol_t = char>
class eNFA;

/**
 *  A nondeterministic finite automaton (NFA) with epsilon connections
 *
 *  Like a normal NFA a eNFA is defined by
 *
 *   - a set of states @f$ Q @f$
 *     (implicitly given by the template parameter `state_t`)
 *   - an alphabet or set of input symbols @f$ \Sigma @f$
 *     (implicitly given by the template parameter `symbol_t`)
 *   - a set of initial states @f$ S @f$
 *     (stored in the member #S)
 *   - a set of finial states @f$ F @f$
 *     (stored in the member #F)
 *   - a transition function
 *     @f$ T : Q \times ( \Sigma \cup \{ \epsilon \} )
 *       \rightarrow \mathcal{P}(Q) @f$
 *     (stored in the member #T)
 *
 *  compared to a normal NFA the transition function is defined on a broader
 *  domain to include epsilon transitions.
 *
 *  The transition function is modeled by a std::map, which stores all mappings
 *  to non-empty sets of states. The domain of the transition function is @f$
 *  Q \times ( \Sigma \cup \{ \epsilon \} ) @f$. These values (that are used as
 *  keys in the map) are modeled by the #key_t type, which is a pair of a state
 *  and an optional symbol. The absence of a symbol corresponds to the epsilon
 *  transitions.
 */
template<typename state_t, typename symbol_t>
class eNFA {
    public:
        /** Set of states */
        using set_t = std::set<state_t>;
        /** Domain of the transition function, epsilon is modeled by the empty
         * std::optional */
        using key_t = std::pair<state_t, std::optional<symbol_t>>;

    private:
        /** empty optional symbol is used to represent the epsilon transitions
         * */
        static constexpr std::optional<symbol_t> e {};
        
        /**
         *  Transition function
         *
         *  pairs that are not in the map are mapped implicitely to the empty
         *  set
         */
        std::map<key_t, set_t> T;
       
        /** Set of start states */
        set_t S;
        /** Set of final states */
        set_t F;

    public:
        eNFA(std::map<key_t, set_t> T, set_t S, set_t F);

        /**
         *  Print dot code that represents the eNFA
         */
        void print_dot() const;

        /** 
         *  Compute the epsilon closure E(P) off a set of states P.
         *
         *  The epsilon closure is defined as the set of all states @f$ q \in Q
         *  @f$ that are 'reachable' by epsilon connections from any state @f$
         *  p \in P @f$.
         *
         *  Does a breadth first search along epsilon transitions starting from
         *  all the states in P.
         */
        set_t E(const set_t& P) const;

        /**
         *  Convert the eNFA to a normal NFA.
         *
         *  This works by computing the epsilon closure of the set of starting
         *  states and of each of the transition target sets.
         */
        NFA<state_t, symbol_t> to_nfa() const;

        /**
         *  Convert the eNFA directly into a DFA
         *
         *  This uses the powerset or subset construction.
         */
        DFA<state_t, symbol_t> powerset() const;

    public:
        eNFA(const std::string& s);

        friend eNFA<state_t, symbol_t> operator | <state_t, symbol_t>(
                const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs);
        friend eNFA<state_t, symbol_t> operator & <state_t, symbol_t>(
                const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs);

};

/**
 *  A nondeterministic finite automaton (NFA)
 *
 *  A NFA is defined by
 *
 *   - a set of states @f$ Q @f$
 *     (implicitly given by the template parameter `state_t`)
 *   - an alphabet or set of input symbols @f$ \Sigma @f$
 *     (implicitly given by the template parameter `symbol_t`)
 *   - a set of initial states @f$ S @f$
 *     (stored in the member #S)
 *   - a set of finial states @f$ F @f$
 *     (stored in the member #F)
 *   - a transition function
 *     @f$ T : Q \times \Sigma \rightarrow \mathcal{P}(Q) @f$
 *     (stored in the member #T)
 *
 *  The transition function is modeled by a std::map, which stores all mappings
 *  to non-empty sets of states. The domain values are modeled by a std::pair
 *  of `state_t` and `symbol_t`.
 */
template<typename state_t, typename symbol_t>
class NFA {
    public:
        /** Set of states */
        using set_t = std::set<state_t>;
        /** Domain values of the transition function */
        using key_t = std::pair<state_t, symbol_t>;

    private:
        /** Transition function
         *
         *  pairs that are not in the map are mapped implicitely to the empty
         *  set
         */
        std::map<key_t, set_t> T;
       
        /** Set of start states */
        set_t S;
        /** Set of final states */
        set_t F;

    public:
        NFA(std::map<key_t, set_t> T, set_t S, set_t F);
                
        /**
         *  Print dot code that represents the eNFA
         */
        void print_dot() const;

        DFA<state_t, symbol_t> powerset() const;
};


/**
 *  A deterministic finite automaton (DFA)
 */
template<typename state_t, typename symbol_t>
class DFA {
    public:
        using set_t = std::set<state_t>;
        using key_t = std::pair<state_t, symbol_t>;

    private:
        std::map<key_t, state_t> T;
        state_t S;
        set_t F;

    public:
        DFA(std::map<key_t, state_t> T, state_t S, set_t F);
        void print_dot() const;
        DFA<state_t, symbol_t> minimize() const;
        NFA<state_t, symbol_t> reverse() const;
};

    



/* *********************** *
 * implementation          *
 * *********************** */

template<typename state_t, typename symbol_t>
DFA<state_t, symbol_t>::DFA(std::map<key_t, state_t> T, state_t S, set_t F):
    T(std::move(T)),
    S(std::move(S)),
    F(std::move(F))
{};

template<typename state_t, typename symbol_t>
void DFA<state_t, symbol_t>::print_dot() const {
    using fmt::print;
    print("digraph nfa {{\n");
    print("  rankdir=LR;\n");
    if ( F.size() ) {
        print("  node [shape = doublecircle]; ");
        for ( const auto& y : F ) {
            print("{} ",y);
        }
        print(";\n");
    }
    print("  node [shape = circle, style = filled]; {};\n", S);
    print("  node [shape = circle, style = \"\"];\n");
    for ( const auto& [key, y] : T ) {
        const auto& [x, a] = key;
        print("  {} -> {} [label=\"'{}'\"];\n", x, y, a);
    }
    print("}}\n");
}

template<typename state_t, typename symbol_t>
DFA<state_t, symbol_t> DFA<state_t, symbol_t>::minimize() const {
    using dfa_t = DFA<state_t, symbol_t>;


    return dfa_t(T, S, F);
}

template<typename state_t, typename symbol_t>
NFA<state_t, symbol_t> DFA<state_t, symbol_t>::reverse() const {
    using nfa_t = NFA<state_t, symbol_t>;
    std::map<key_t, set_t> RT;
    for ( const auto& [key, target]: T ) {
        RT[{target, key.second}].insert(key.first);
    }
    return nfa_t(RT, F, {S});
}

template<typename state_t, typename symbol_t>
NFA<state_t, symbol_t>::NFA(std::map<key_t, set_t> T, set_t S, set_t F):
            T(std::move(T)),
            S(std::move(S)),
            F(std::move(F))
{};

template<typename state_t, typename symbol_t>
void NFA<state_t, symbol_t>::print_dot() const {
    using fmt::print;
    print("digraph nfa {{\n");
    print("  rankdir=LR;\n");
    if ( F.size() ) {
        print("  node [shape = doublecircle]; ");
        for ( const auto& y : F ) {
            print("{} ",y);
        }
        print(";\n");
    }
    print("  node [shape = circle, style = filled]; ");
    for ( const auto& y : S ) {
        print("{} ",y);
    }
    print(";\n");
    print("  node [shape = circle, style = \"\"];\n");
    for ( const auto& [key, set] : T ) {
        const auto& [x, a] = key;
        for ( const auto& y : set ) {
            print("  {} -> {} [label=\"'{}'\"];\n", x, y, a);
        }
    }
    print("}}\n");
}

template<typename state_t, typename symbol_t>
DFA<state_t, symbol_t> NFA<state_t, symbol_t>::powerset() const {
    using dfa_t = DFA<state_t, symbol_t>;
    constexpr symbol_t symbol_min = std::numeric_limits<symbol_t>::min();
    constexpr symbol_t symbol_max = std::numeric_limits<symbol_t>::max();

    std::map<std::pair<set_t, symbol_t>, set_t> PT;
    std::map<set_t, state_t> PQ {{S,0}};
    std::queue<set_t> Q; Q.push(S);
    state_t i = std::numeric_limits<state_t>::min() + 1;

    while ( !Q.empty() ) {
        auto pq = std::move(Q.front()); Q.pop();
        std::map<symbol_t, set_t> pt;
        for ( const auto& q : pq ) {
            for ( const auto& [key, p] : stdr::subrange(
                    T.lower_bound({q, symbol_min}), T.upper_bound({q, symbol_max})) ) {
                auto& s = pt[key.second];
                for ( const auto& x : p ) {
                    s.insert(x);
                }
            }
        }

        for ( auto [ sym, set] : pt ) {
            PT.insert({{pq,sym}, set});
            if ( PQ.insert({set, i}).second ) {
                i++;
                Q.push(set);
            }
        }
    }

    std::map<std::pair<state_t, symbol_t>, state_t> DT;
    std::set<state_t> DF;
    state_t DS {PQ.at(S)};

    for ( const auto& [pq, dq] : PQ ) {
        for ( const auto& q: pq ) {
            if ( F.contains(q) ) {
                DF.insert(dq);
                break;
            }
        }

        for ( const auto& [key , pp] : stdr::subrange(
                PT.lower_bound({pq, symbol_min}), PT.upper_bound({pq, symbol_max})) ) {
            DT.insert({{dq, key.second}, PQ.at(pp)});
        }
    }

    return dfa_t(DT, DS, DF);
}

template<typename state_t, typename symbol_t>
eNFA<state_t,symbol_t>::eNFA(std::map<key_t, set_t> T, set_t S, set_t F):
    T(std::move(T)),
    S(std::move(S)),
    F(std::move(F)) 
{};

template<typename state_t, typename symbol_t>
eNFA<state_t, symbol_t>::set_t eNFA<state_t, symbol_t>::E(const set_t& P) const {
    using queue_t = std::queue<state_t>;
    set_t   EP {P};
    queue_t Q  {{EP.begin(), EP.end()}};

    while ( !Q.empty() ) {
        auto q = std::move(Q.front()); Q.pop();
        // loop over all states, that can be reached via an epsilon
        // transition
        if ( auto it = T.find({q, e}); it != T.end() ) {
            for ( auto& p : it->second ) {
                if ( EP.insert(p).second ) {
                    Q.push(p);
                }
            }
        }
    }

    return EP;
}

template<typename state_t, typename symbol_t>
DFA<state_t, symbol_t> eNFA<state_t, symbol_t>::powerset() const {
    using dfa_t = DFA<state_t, symbol_t>;
    constexpr symbol_t symbol_min = std::numeric_limits<symbol_t>::min();
    constexpr symbol_t symbol_max = std::numeric_limits<symbol_t>::max();

    constexpr state_t s { std::numeric_limits<state_t>::min() };

    // initialize the new transition function map PT and the
    // new final state set PF
    std::map<std::pair<state_t, symbol_t>, state_t> PT {};
    set_t PF {};

    // initialize the new state counter i, the bfs queue Q, and
    // the powerset to new state map PQ
    state_t i {s};
    std::queue<set_t> Q; Q.push(E(S));
    std::map<set_t, state_t> PQ {{Q.front(),i++}};

    while ( !Q.empty() ) {
        auto pq = std::move(Q.front()); Q.pop();
        std::map<symbol_t, set_t> pt;

        // loop over all states in this powerset
        for ( state_t q : pq ) {
            // loop over all transitions from this state
            // (excluding epsilon transitions)
            for ( const auto& [key, pp] : stdr::subrange(
                    T.lower_bound({q, symbol_min}), T.upper_bound({q, symbol_max})) ) {
                auto& ps = pt[*key.second];
                for ( state_t p : pp ) {
                    ps.insert(p);
                }
            }
        }

        for ( auto& [s, pp] : pt ) {
            // try to insert the new state E(pp) into to mapping
            // if it was not included in PQ this means, that we
            // visited this state for the first time and must add
            // it to the queue and increase the state counter.
            auto [map, ins] = PQ.insert({E(pp), i});
            if ( ins ) {
                i++;
                Q.push(map->first);
            }

            // add the transition
            PT.insert({{PQ[pq], s}, map->second});
            
            // look if new state is a final state and add it to
            // the set of final states if so.
            for ( auto f : F ) {
                if ( map->first.contains(f) ) {
                    PF.insert(map->second);
                    break;
                }
            }
        }
    }

    return dfa_t(PT, s, PF);
}

template<typename state_t, typename symbol_t>
eNFA<state_t, symbol_t>::eNFA(const std::string& s):
    T{}, S{}, F{}
{
    state_t i = std::numeric_limits<state_t>::min();

    S.insert(i);
    for ( auto c : s ) {
        T.insert({{i,c},{i+1}});
        i++;
    }
    F.insert(i);
}


template<typename state_t, typename symbol_t>
eNFA<state_t, symbol_t> operator | (
        const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs) {
    using enfa_t = eNFA<state_t, symbol_t>;
    using key_t = typename enfa_t::key_t;
    using set_t = typename enfa_t::set_t;

    state_t i {std::numeric_limits<state_t>::min()};
    std::map<state_t, state_t> vmap {};
    auto get_mapped = [&vmap, &i](state_t p) -> state_t {
        auto [map, ins] = vmap.insert({p,i});
        if ( ins ) {
            i++;
        }
        return map->second;
    };

    std::map<key_t, set_t> NT {};
    state_t s {i++};
    state_t f {i++};

    auto& start { NT[{s, enfa_t::e}] };
    for ( auto& xhs : {lhs, rhs} ) {
        vmap.clear();
        for ( auto p : xhs.S ) {
            start.insert(get_mapped(p));
        }
        for ( auto& [key, target] : xhs.T ) {
            auto& current = NT[{get_mapped(key.first), key.second}];
            for ( auto p : target ) {
                current.insert(get_mapped(p));
            }
        }
        for ( auto p : xhs.F ) {
            NT[{get_mapped(p), enfa_t::e}].insert(f);
        }
    }

    return enfa_t {NT, {s}, {f}};
}

template<typename state_t, typename symbol_t>
eNFA<state_t, symbol_t> operator & (
        const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs) {

    using enfa_t = eNFA<state_t, symbol_t>;
    using key_t = typename enfa_t::key_t;
    using set_t = typename enfa_t::set_t;

    state_t i {std::numeric_limits<state_t>::min()};
    std::map<state_t, state_t> vmap {};
    auto get_mapped = [&vmap, &i](state_t p) -> state_t {
        auto [map, ins] = vmap.insert({p,i});
        if ( ins ) {
            i++;
        }
        return map->second;
    };

    std::map<key_t, set_t> NT {};
    set_t NS {};
    set_t NF {};

    // handle lhs
    for ( auto p : lhs.S ) {
        NS.insert(get_mapped(p));
    }

    for ( auto& [key, target] : lhs.T ) {
        auto& current = NT.insert({{get_mapped(key.first), key.second}, {}}).first->second;
        for ( auto p : target ) {
            current.insert(get_mapped(p));
        }
    }

    // handle connection
    state_t y {i++};
    for ( auto p : lhs.F ) {
        NT[{get_mapped(p), enfa_t::e}].insert(y);
    }

    // handle rhs
    vmap.clear();
    auto& conn = NT.insert({{y, enfa_t::e}, {}}).first->second;
    for ( auto p : rhs.S ) {
        conn.insert(get_mapped(p));
    }
    for ( auto& [key, target] : rhs.T ) {
        auto [current, _] = NT.insert({{get_mapped(key.first), key.second}, {}});
        for ( auto p : target ) {
            current->second.insert(get_mapped(p));
        }
    }
    for ( auto p : rhs.F ) {
        NF.insert(get_mapped(p));
    }

    return enfa_t {NT, NS, NF};
}


template<typename state_t, typename symbol_t>
void eNFA<state_t,symbol_t>::print_dot() const {
    using fmt::print;
    print("digraph nfa {{\n");
    print("  rankdir=LR;\n");

    if ( F.size() ) {
        print("  node [shape = doublecircle]; ");
        for ( const auto& y : F ) {
            print("{} ",y);
        }
        print(";\n");
    }

    print("  node [shape = circle, style = filled]; ");
    for ( const auto& y : S ) {
        print("{} ",y);
    }
    print(";\n");
    print("  node [shape = circle, style = \"\"];\n");
    for ( const auto& [key, set] : T ) {
        const auto& [x, a] = key;
        for ( const auto& y : set ) {
            std::string s = a.has_value() ? fmt::format("'{}'",*a) : "\u03b5"s;
            print("  {} -> {} [label=\"{}\"];\n", x, y, s);
        }
    }
    print("}}\n");
}

