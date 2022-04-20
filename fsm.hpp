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
class eNFA;

template<typename R, typename symbol_t>
concept symbol_range =
    std::ranges::common_range<R> &&    
    std::convertible_to<std::ranges::range_value_t<R>,symbol_t>;

enum class action_type {
    _FIRST, ENTER, LEAVE, STAY, _LAST
};

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

        /* types for action map */
        using akey_t = std::pair<state_t, action_type>;
        using alist_t = std::list<std::function<void(std::size_t)>>;

    public:
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

        /**
         *  Action map
         */
        std::map<akey_t, alist_t> A;

    public:
        /**
         *  Construct using explicit listing
         */
        eNFA(std::map<key_t, set_t> T, set_t S, set_t F, std::map<akey_t, alist_t> A);

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
         *  Convert the eNFA directly into a DFA
         *
         *  This uses the powerset or subset construction.
         */
        DFA<state_t, symbol_t> powerset() const;

    public:
        /**
         *  Construct using a range of symbol_t
         *
         *  This constructs a machine that matches the exact sequence of
         *  symbols specified in the range
         */
        template<symbol_range<symbol_t> R>
        explicit eNFA(R&& r);

        /**
         *  Construct the union of two enfa using Thompson's construction
         */
        friend eNFA<state_t, symbol_t> operator | <state_t, symbol_t>(
                const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs);
        /**
         *  Construct the concatenation of two enfa using Thompson's construction
         */
        friend eNFA<state_t, symbol_t> operator & <state_t, symbol_t>(
                const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs);

};

/**
 *  A deterministic finite automaton (DFA)
 */
template<typename state_t, typename symbol_t>
class DFA {
    public:
        using set_t = std::set<state_t>;
        using key_t = std::pair<state_t, symbol_t>;

        using akey_t = std::pair<state_t, action_type>;
        using alist_t = std::list<std::function<void(std::size_t)>>;
    public:
        std::map<key_t, state_t> T;
        state_t S;
        set_t F;

        /**
         *  Action map
         */
        std::map<akey_t, alist_t> A;

    public:
        DFA(std::map<key_t, state_t> T, state_t S, set_t F, std::map<akey_t, alist_t> A);
        void print_dot() const;
        DFA<state_t, symbol_t> minimize() const;
        eNFA<state_t, symbol_t> reverse() const;

        template<symbol_range<symbol_t> R>
        bool match(R&& r);
};

/* *********************** *
 * implementation          *
 * *********************** */

template<typename state_t, typename symbol_t>
DFA<state_t, symbol_t>::DFA(std::map<key_t, state_t> T, state_t S, set_t F,
        std::map<akey_t, alist_t> A):
    T(std::move(T)),
    S(std::move(S)),
    F(std::move(F)),
    A(std::move(A))
{};

template<typename state_t, typename symbol_t>
void DFA<state_t, symbol_t>::print_dot() const {
    using fmt::print;
    print("digraph nfa {{\n");
    print("  rankdir=LR;\n");

    set_t FuS {};
    set_t FwS {};
    set_t SwF {};
    std::ranges::set_intersection(F,set_t{S},std::inserter(FuS, FuS.end()));
    std::ranges::set_difference(F,set_t{S},std::inserter(FwS, FwS.end()));
    std::ranges::set_difference(set_t{S},F,std::inserter(SwF, SwF.end()));

    if ( FuS.size() ) {
        print("  node [shape = doublecircle, style = filled]; ");
        for ( auto y : FuS ) {
            print("{} ",y);
        }
        print(";\n");
    }

    if ( FwS.size() ) {
        print("  node [shape = doublecircle, style = \"\"]; ");
        for ( auto y : FwS ) {
            print("{} ",y);
        }
        print(";\n");
    }

    if ( SwF.size() ) {
        print("  node [shape = circle, style = filled]; ");
        for ( auto y : SwF ) {
            print("{} ",y);
        }
        print(";\n");
    }

    print("  node [shape = circle, style = \"\"];\n");
    for ( const auto& [key, y] : T ) {
        const auto& [x, a] = key;
        print("  {} -> {} [label=\"'{}'\"];\n", x, y, a);
    }
    for ( const auto& [key, _] : A ) {
        print("  {} -> a;\n", key.first);
    }
    print("}}\n");
}

template<typename state_t, typename symbol_t>
DFA<state_t, symbol_t> DFA<state_t, symbol_t>::minimize() const {
    using dfa_t = DFA<state_t, symbol_t>;


    return dfa_t(T, S, F);
}

template<typename state_t, typename symbol_t>
eNFA<state_t, symbol_t> DFA<state_t, symbol_t>::reverse() const {
    using enfa_t = eNFA<state_t, symbol_t>;
    std::map<typename enfa_t::key_t, set_t> RT;
    for ( const auto& [key, target]: T ) {
        RT[{target, key.second}].insert(key.first);
    }
    return enfa_t(RT, F, {S});
}

template<typename state_t, typename symbol_t>
template<symbol_range<symbol_t> R>
bool DFA<state_t, symbol_t>::match(R&& r) {
    state_t q {S};
    std::size_t pos {0};
    for ( symbol_t c : r ) {
        if ( auto t = T.find({q, c}); t == T.end() ) {
            return false;
        } else {
            if ( q == t->second ) {
                if ( auto a = A.find({q, action_type::STAY}); a != A.end() ) {
                    for ( const auto& f : a->second ) {
                        f(pos);
                    }
                }
            } else {
                if ( auto a = A.find({q, action_type::LEAVE}); a != A.end() ) {
                    for ( const auto& f : a->second ) {
                        f(pos);
                    }
                }
                if ( auto a = A.find({t->second, action_type::ENTER}); a != A.end() ) {
                    for ( const auto& f : a->second ) {
                        f(pos);
                    }
                }
            }

            q = t->second;
        }
        pos++;
    }
    return F.contains(q);
}

template<typename state_t, typename symbol_t>
eNFA<state_t,symbol_t>::eNFA(std::map<key_t, set_t> T, set_t S, set_t F,
        std::map<akey_t, alist_t> A):
    T(std::move(T)),
    S(std::move(S)),
    F(std::move(F)),
    A(std::move(A))
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
    decltype(dfa_t::T) PT {};
    decltype(dfa_t::A) PA {};
    decltype(dfa_t::F) PF {};

    // initialize the new state counter i, the bfs queue Q, and
    // the powerset to new state map PQ
    state_t i {s};
    std::map<set_t, state_t> PQ {{E(S),i++}};
    std::queue<typename decltype(PQ)::iterator> Q; Q.push(PQ.begin());

    while ( !Q.empty() ) {
        const auto& [pq, pi] = *Q.front(); Q.pop();
        std::map<symbol_t, set_t> pt;

        for ( state_t q : pq ) {
            // look if the state is a final state
            if ( F.contains(q) ) {
                PF.insert(pi);
            }

            // look if the state contains an action
            static_assert(false, "move this down, so that we can tranfer actions more precisely");
            for ( const auto& [key, alist] : stdr::subrange(
                        A.lower_bound({q, action_type::_FIRST}),
                        A.upper_bound({q, action_type::_LAST})) ) {
                auto& palist = PA[{pi, key.second}];
                stdr::copy(alist, std::back_inserter(palist));
            }

            // loop over all transitions from this state
            // (excluding epsilon transitions)
            for ( const auto& [key, pp] : stdr::subrange(
                    T.lower_bound({q, symbol_min}), T.upper_bound({q, symbol_max})) ) {
                auto& ps = pt[*key.second];
                stdr::copy(pp, std::inserter(ps, ps.end()));
            }
        }

        for ( auto& [s, pp] : pt ) {
            // try to insert the new state E(pp) into to mapping
            // if it was not included in PQ this means, that we
            // visited this state for the first time and must add
            // it to the queue and increase the state counter.
            auto [map, ins] = PQ.insert({E(pp), i});
            if ( ins ) {
                Q.push(map);
                i++;
            }

            // add the transition
            PT.insert({{pi, s}, map->second});
        }
    }

    return dfa_t(PT, s, PF, PA);
}

template<typename state_t, typename symbol_t>
template<symbol_range<symbol_t> R>
eNFA<state_t, symbol_t>::eNFA(R&& r):
    T{}, S{}, F{}, A{}
{
    state_t i = std::numeric_limits<state_t>::min();

    S.insert(i);
    for ( auto c : r ) {
        T.insert({{i,c},{i+1}});
        i++;
    }
    F.insert(i);
}


template<typename state_t, typename symbol_t>
eNFA<state_t, symbol_t> operator | (
        const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs) {
    using enfa_t = eNFA<state_t, symbol_t>;

    state_t i {std::numeric_limits<state_t>::min()};
    std::map<state_t, state_t> vmap {};
    auto get_mapped = [&vmap, &i](state_t p) -> state_t {
        auto [map, ins] = vmap.insert({p,i});
        if ( ins ) {
            i++;
        }
        return map->second;
    };

    decltype(enfa_t::T) NT {};
    decltype(enfa_t::A) NA {};
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
        for ( auto& [key, action] : xhs.A ) {
            NA.insert({{get_mapped(key.first), key.second}, action});
        }
    }

    return enfa_t {NT, {s}, {f}, NA};
}

template<typename state_t, typename symbol_t>
eNFA<state_t, symbol_t> operator & (
        const eNFA<state_t, symbol_t>& lhs, const eNFA<state_t, symbol_t>& rhs) {

    using enfa_t = eNFA<state_t, symbol_t>;

    state_t i {std::numeric_limits<state_t>::min()};
    std::map<state_t, state_t> vmap {};
    auto get_mapped = [&vmap, &i](state_t p) -> state_t {
        auto [map, ins] = vmap.insert({p,i});
        if ( ins ) {
            i++;
        }
        return map->second;
    };

    decltype(enfa_t::T) NT {};
    decltype(enfa_t::A) NA {};
    decltype(enfa_t::S) NS {};
    decltype(enfa_t::F) NF {};

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
    for ( auto& [key, action] : lhs.A ) {
        NA.insert({{get_mapped(key.first), key.second}, action});
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
    for ( auto& [key, action] : rhs.A ) {
        NA.insert({{get_mapped(key.first), key.second}, action});
    }

    return enfa_t {NT, NS, NF, NA};
}


template<typename state_t, typename symbol_t>
void eNFA<state_t,symbol_t>::print_dot() const {
    using fmt::print;
    print("digraph nfa {{\n");
    print("  rankdir=LR;\n");

    set_t FuS {};
    set_t FwS {};
    set_t SwF {};
    std::ranges::set_intersection(F,S,std::inserter(FuS, FuS.end()));
    std::ranges::set_difference(F,S,std::inserter(FwS, FwS.end()));
    std::ranges::set_difference(S,F,std::inserter(SwF, SwF.end()));

    if ( FuS.size() ) {
        print("  node [shape = doublecircle, style = filled]; ");
        for ( auto y : FuS ) {
            print("{} ",y);
        }
        print(";\n");
    }

    if ( FwS.size() ) {
        print("  node [shape = doublecircle, style = \"\"]; ");
        for ( auto y : FwS ) {
            print("{} ",y);
        }
        print(";\n");
    }

    if ( SwF.size() ) {
        print("  node [shape = circle, style = filled]; ");
        for ( auto y : SwF ) {
            print("{} ",y);
        }
        print(";\n");
    }

    print("  node [shape = circle, style = \"\"];\n");
    for ( const auto& [key, set] : T ) {
        const auto& [x, a] = key;
        for ( const auto& y : set ) {
            std::string s = a.has_value() ? fmt::format("'{}'",*a) : "\u03b5"s;
            print("  {} -> {} [label=\"{}\"];\n", x, y, s);
        }
    }

    for ( const auto& [key, _] : A ) {
        print("  {} -> a;\n", key.first);
    }
    print("}}\n");
}

/* *************
 * Constants
 * ************/

eNFA<> enfa_number(std::size_t& begin, std::size_t& end) {
    return {
        {
            { {0, {'0'}}, {0} },
            { {0, {'1'}}, {0} },
            { {0, {'2'}}, {0} },
            { {0, {'3'}}, {0} },
            { {0, {'4'}}, {0} },
            { {0, {'5'}}, {0} },
            { {0, {'6'}}, {0} },
            { {0, {'7'}}, {0} },
            { {0, {'8'}}, {0} },
            { {0, {'9'}}, {0} },
        },
        {0},
        {0},
        {
            { {0,action_type::STAY}, { [&end  ](std::size_t pos){end   = pos;} }},
            { {0,action_type::ENTER}, { [&begin](std::size_t pos){begin = pos;} }}
        }
    };
}


