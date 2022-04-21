project('beast_server', 'cpp', default_options: 'cpp_std=c++20')

conan_pkgs = {
  'boost'    : '1.78.0',
  'argparse' : '2.2',
  'fmt'      : '8.1.1'
}

conan_deps = []
foreach pkg_name, conan_version : conan_pkgs
  module_path = meson.current_build_dir() / 'conan_cmake'
  run_command('conan', 'install', pkg_name + '/' + conan_version + '@',
    '-if', module_path,
    '-g', 'cmake_find_package', check:true)
  conan_deps += dependency(pkg_name, method: 'cmake', cmake_module_path: module_path)
endforeach

cmake = import('cmake')
opt_var = cmake.subproject_options()
opt_var.add_cmake_defines({
  'CMAKE_MODULE_PATH': meson.current_build_dir() / 'conan_cmake',
  'BOOST_URL_BUILD_TESTS': 'OFF',
})
sub = cmake.subproject('url', options : opt_var)

executable('server', ['main.cpp', 'router.cpp'],
  dependencies: conan_deps + [sub.dependency('boost_url')]
)
