project('beast_server', 'c++', default_options: 'cpp_std=c++20')

conan_pkgs = [
  'boost'    : '1.78.0',
  'argparse' : '2.2',
  'fmt'      : '8.1.1'
]

conan_deps = []
foreach pkg_name, conan_version : conan_pkgs
  module_path = meson.current_build_dir() / 'conan_cmake'
  run_command('conan', 'install', pkg_name + conan_version + '@',
    '-if', module_path,
    '-g', 'cmake_find_package', check:true)
  deps += dependency(pkgname, method: 'cmake', cmake_module_path: module_path)
endforeach

executable('server', ['main.cpp'],
  dependencies: deps
)
