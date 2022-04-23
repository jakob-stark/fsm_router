#ifndef FSM2_HPP
#define FSM2_HPP

#include <vector>
#include <set>
#include <functional>
#include <map>
#include <queue>

namespace fsm2 {

struct FSM {

    using index_t  = std::size_t;
    using symbol_t = char;
    using sflag_t = unsigned char;
    //using action_t = std::function<void(std::size_t pos)>;
    using action_t = char;

    static_assert(std::is_signed_v<symbol_t>,
            "implementation requires char to be signed");

    static constexpr symbol_t epsilon = -1;
    static constexpr symbol_t leaving = -2;

    struct target_t {
        index_t state;
        std::vector<action_t> actions;
    };

    using state_t = std::multimap<symbol_t, target_t>;

    // states list contains indices into transitions list. First state is the
    // start state, last states (behind final_state_cut) are the final states.
    std::vector<state_t> states;



    static inline void remap(FSM& fsm, index_t state_offset);
    static inline void relink(FSM& fsm, index_t leave_to);

    bool contains_epsilon_transitions() const;
    bool is_deterministic() const;

    std::string to_graphviz() const;

    using state_set_t = std::set<index_t>;
    state_set_t epsilon_closure(const state_set_t& input_set) const;

    FSM powerset() const;
    FSM reverse() const;
    FSM brzozowski() const;

    bool match(const std::string& str) const;

    friend FSM from_string(const std::string& str, char x);
    friend FSM concatenate(FSM lhs, FSM rhs);
    friend FSM alternative(FSM lhs, FSM rhs);
    friend FSM kleene_star(FSM lhs);
};

FSM from_string(const std::string& str, char x);

}

#endif//FSM2_HPP
