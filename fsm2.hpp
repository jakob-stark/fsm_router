#ifndef FSM2_HPP
#define FSM2_HPP

#include <vector>
#include <set>
#include <functional>
#include <map>

namespace fsm2 {

struct FSM {

    using index_t  = std::size_t;
    using symbol_t = char;
    using sflag_t = unsigned char;
    using action_t = std::function<void(std::size_t pos)>;

    static_assert(std::is_signed_v<symbol_t>, "implementation requires char to be signed");

    static constexpr symbol_t epsilon = -1;
    static constexpr symbol_t leaving = -2;

    struct target_t {
        index_t state;
        std::vector<index_t> actions;
    };

    struct state_t {
        std::multimap<symbol_t, target_t> transitions;
    };

    // states list contains indices into transitions list. First state is the
    // start state, last states (behind final_state_cut) are the final states.
    std::vector<state_t> states;
    std::vector<action_t> actions;

    bool contains_epsilon_transitions() const;
    bool is_deterministic() const;

    std::string to_graphviz() const;

    friend FSM powerset(const FSM& fsm);

    friend FSM from_string(const std::string& str);
    friend FSM concatenate(const FSM& lhs, const FSM& rhs);
    friend FSM alternative(const FSM& lhs, const FSM& rhs);

};

FSM from_string(const std::string& str);


}

#endif//FSM2_HPP
