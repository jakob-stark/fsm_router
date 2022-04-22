#ifndef ROUTER_HPP
#define ROUTER_HPP

#include <functional>
#include <string_view>

class router {
    private:
        std::vector<std::function<void()>> handlers;
        unsigned int route_id;
    public:
        router() {}
        
        template<typename handler_t>
        void add_route(std::string_view route, handler_t&& handler) {
        }

        template<typename request_t>
        void match(std::string_view path, request_t&& request);

};

#endif//ROUTER_HPP
