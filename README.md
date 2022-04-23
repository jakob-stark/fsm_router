# fsm_router
An url/path router with an optimized deterministic finite automaton for fast route matching.

In many web servers and especially REST toolkits, one of the central parts is an
url router. If a new request must be processed, a handler function must be called
based on the url path that is requested. There are different approaches to make
this decision.

A naive approach would be to store each route and the corresponding handler function
in a table and search for the route each time a new request is made. There are two
problems with this approach however:

1. Let's take an example. Say we have the following routes
   ```
   route_a = "abcdefghijklmnopqrstuvwxyz/a";
   route_b = "abcdefghijklmnopqrstuvwxyz/b";
   
   request = "abcdefghijklmnopqrstuvwxyz/z";
   ```
   A search algorithm will have to compare the request with both routes to ultimately
   decide, that it does not match either of them. Although the large common prefix did
   already match, it is checked again while the second entry in the routing table is
   scanned.
   
2. Often it is desirable that the routes include wildcards, that can match e.g. any
   number or an arbitrary string.
   ```
   route_c = "foo/<int>/bar/<string>";
   ```
   It is impossible to store each possible value of the wildcard in a table. Thus some
   kind of dynamic pattern matching must be applied.
   
The goal of this project is to construct the theoretical optimal matcher for a given
set of routes. Each route can be specified with a regular expression and a deterministic
finite state machine is contructed that matches those routes as much simultaneously
as possible.

# Build

### Build requirements

 - [meson](https://mesonbuild.com/)
 - [conan](https://conan.io/)
 - ninja
 - cmake
 - A C++20 compiler and standard library

### How to build

```
$ meson setup build
$ cd build
$ ninja all
```
