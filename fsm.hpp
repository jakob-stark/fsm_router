#include <map>
#include <set>
#include <queue>
#include <fmt/core.h>

class DFA;
    
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
template<typename state_t = unsigned int, typename symbol_t = char>
class NFA {
    public:
        /** Set of states */
        using set_t = std::set<state_t>;
        /** Domain values of the transition function */
        using key_t = std::pair<state_t, symbol_t>;

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
        NFA(decltype(T) T, decltype(S) S, decltype(F) F):
            T(std::move(T)),
            S(std::move(S)),
            F(std::move(F)) {};

        void print_dot() const {
            using fmt::print;
            
            print("digraph nfa {{\n");
            print("  rankdir=LR;\n");

            print("  node [shape = doublecircle]; ");
            for ( const auto& y : F ) {
                print("{} ",y);
            }
            print(";\n");

            print("  node [shape = circle, style = filled]; ");
            for ( const auto& y : S ) {
                print("{} ",y);
            }
            print(";\n");

            print("  node [shape = circle, style = \"\"];\n");
            for ( const auto& [key, set] : T ) {
                const auto& [x, a] = key;
                for ( const auto& y : set ) {
                    print("  {} -> {} [label={}];\n", x, y, a);
                }
            }
            print("}}\n");
        }
        //DFA powerset;
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
template<typename state_t = unsigned int, typename symbol_t = char>
class eNFA {
    public:

        /** Set of states */
        using set_t = std::set<state_t>;
        /** Domain of the transition function, epsilon is modeled by the empty
         * std::optional */
        using key_t = std::pair<state_t, std::optional<symbol_t>>;

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
        eNFA(decltype(T) T, decltype(S) S, decltype(F) F):
            T(std::move(T)),
            S(std::move(S)),
            F(std::move(F)) {};

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
        set_t E(const set_t& P) const {
            using queue_t = std::queue<state_t>;
            set_t   EP {P};
            queue_t Q  {{EP.begin(), EP.end()}};

            while ( !Q.empty() ) {
                auto q = std::move(Q.front()); Q.pop();
                // loop over all states, that can be reached via an epsilon
                // transition
                if ( auto it = T.find(std::make_pair(q, e)); it != T.end() ) {
                    for ( auto& p : it->second ) {
                        if ( EP.insert(p).second ) {
                            Q.push(p);
                        }
                    }
                }
            }

            return EP;
        }

        /**
         *  Convert the eNFA to a normal NFA.
         *
         *  This works by computing the epsilon closure of the set of starting
         *  states and of each of the transition target sets.
         */
        auto to_nfa() const {
            using nfa_t = NFA<state_t, symbol_t>;
            
            // construct T'(q,a) = E(T(q,a))
            std::map<typename nfa_t::key_t, set_t> ET {};

            for ( const auto& [key, P] : T ) {
                const auto& [s, c] = key;
                if ( c.has_value() ) {
                    ET.insert({{s, *c},E(P)});
                }
            }

            return nfa_t{ET, E(S), F};
        };
};


class DFA {


};

