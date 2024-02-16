project('tskit', ['c', 'cpp'],
    version: files('VERSION.txt'),    
    default_options: ['c_std=c99', 'cpp_std=c++11']
)

kastore_proj = subproject('kastore')
kastore_dep = kastore_proj.get_variable('kastore_dep')
kastore_inc = kastore_proj.get_variable('kastore_inc')

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)
lib_deps = [m_dep, kastore_dep]

extra_c_args = [
    '-Wall', '-Wextra', '-Werror', '-Wpedantic', '-W',
    '-Wmissing-prototypes',  '-Wstrict-prototypes',
    '-Wconversion', '-Wshadow', '-Wpointer-arith', '-Wcast-align',
    '-Wcast-qual', '-Wwrite-strings', '-Wnested-externs',
    '-fshort-enums', '-fno-common']

lib_sources = [
    'tskit/core.c', 'tskit/tables.c', 'tskit/trees.c',
    'tskit/genotypes.c', 'tskit/stats.c', 'tskit/convert.c', 'tskit/haplotype_matching.c']
lib_headers = [
    'tskit/core.h', 'tskit/tables.h', 'tskit/trees.h',
    'tskit/genotypes.h', 'tskit/stats.h', 'tskit/convert.h', 'tskit/haplotype_matching.h']

# Subprojects use the static library for simplicity.
tskit_inc = [kastore_inc, include_directories(['.'])]
tskit_lib = static_library('tskit',
    sources: lib_sources, dependencies: lib_deps)
tskit_dep = declare_dependency(include_directories:tskit_inc, link_with: tskit_lib)

if not meson.is_subproject()

    # Shared library install target.
    shared_library('tskit',
        sources: lib_sources, dependencies: lib_deps, c_args: extra_c_args, install: true)
    install_headers('tskit.h')
    install_headers(lib_headers, subdir: 'tskit')

    cunit_dep = dependency('cunit')
    # We don't specify extra C args here as CUnit won't pass the checks.
    test_lib = static_library('testlib',
        sources: ['tests/testlib.c'], dependencies: [cunit_dep, kastore_dep, tskit_dep])

    test_core = executable('test_core',
        sources: ['tests/test_core.c'],
        link_with: [tskit_lib, test_lib],
        c_args: extra_c_args+['-DMESON_PROJECT_VERSION="@0@"'.format(meson.project_version())],
        dependencies: kastore_dep,
        )
    test('core', test_core)

    test_tables = executable('test_tables',
        sources: ['tests/test_tables.c'],
        link_with: [tskit_lib, test_lib], c_args: extra_c_args, dependencies: kastore_dep)
    test('tables', test_tables)

    test_trees = executable('test_trees',
        sources: ['tests/test_trees.c'],
        link_with: [tskit_lib, test_lib], c_args: extra_c_args, dependencies: kastore_dep)
    test('trees', test_trees)

    test_genotypes = executable('test_genotypes',
        sources: ['tests/test_genotypes.c'],
        link_with: [tskit_lib, test_lib], c_args: extra_c_args, dependencies: kastore_dep)
    test('genotypes', test_genotypes)

    test_convert = executable('test_convert',
        sources: ['tests/test_convert.c'],
        link_with: [tskit_lib, test_lib], c_args: extra_c_args, dependencies: kastore_dep)
    test('convert', test_convert)

    test_stats = executable('test_stats',
        sources: ['tests/test_stats.c'],
        link_with: [tskit_lib, test_lib], c_args: extra_c_args, dependencies: kastore_dep)
    test('stats', test_stats)

    test_haplotype_matching = executable('test_haplotype_matching',
        sources: ['tests/test_haplotype_matching.c'],
        link_with: [tskit_lib, test_lib], c_args: extra_c_args, dependencies: kastore_dep)
    test('haplotype_matching', test_haplotype_matching)

    test_file_format = executable('test_file_format',
        sources: ['tests/test_file_format.c'],
        link_with: [tskit_lib, test_lib], c_args: extra_c_args, dependencies: kastore_dep)
    test('file_format', test_file_format)

    test_minimal_cpp = executable('test_minimal_cpp',
        sources: ['tests/test_minimal_cpp.cpp'], link_with: [tskit_lib],
        dependencies: kastore_dep)
    test('minimal_cpp', test_minimal_cpp)

    if get_option('build_examples')
      # These example programs use less portable features,
      # and we don't want to always compile them. Use, e.g.,
      # meson build -Dbuild_examples=false
      executable('api_structure',
          sources: ['examples/api_structure.c'], 
          link_with: [tskit_lib], dependencies: lib_deps)
      executable('error_handling',
          sources: ['examples/error_handling.c'], 
          link_with: [tskit_lib], dependencies: lib_deps)
      executable('tree_iteration',
          sources: ['examples/tree_iteration.c'], 
          link_with: [tskit_lib], dependencies: lib_deps)
      executable('tree_traversal',
          sources: ['examples/tree_traversal.c'], 
          link_with: [tskit_lib], dependencies: lib_deps)
      executable('streaming',
          sources: ['examples/streaming.c'], 
          link_with: [tskit_lib], dependencies: lib_deps)
      executable('cpp_sorting_example',
          sources: ['examples/cpp_sorting_example.cpp'], 
          link_with: [tskit_lib], dependencies: lib_deps)
      executable('haploid_wright_fisher',
          sources: ['examples/haploid_wright_fisher.c'], 
          link_with: [tskit_lib], dependencies: lib_deps)
    endif
endif
