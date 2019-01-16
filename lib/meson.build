project('msprime', 'c')

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
gsl_dep = dependency('gsl')
cunit_dep = dependency('cunit')
config_dep = dependency('libconfig')

extra_c_args = [
    '-std=c99', '-Wall', '-Wextra', '-Werror', '-Wpedantic', '-W',
    '-Wmissing-prototypes',  '-Wstrict-prototypes',
    '-Wconversion', '-Wshadow', '-Wpointer-arith', '-Wcast-align',
    '-Wcast-qual', '-Wwrite-strings', '-Wnested-externs',
    '-fshort-enums', '-fno-common']

tskit_dir=join_paths(meson.source_root(), 'tskit/c')
kastore_dir=join_paths(tskit_dir, 'kastore/c')

add_global_arguments(['-I' + kastore_dir, '-I' + tskit_dir], language: 'c')

# We don't use all of tskit, so no point in building it all.
tskit_sources = [
    join_paths(kastore_dir, 'kastore.c'), 
    join_paths(tskit_dir, 'tsk_core.c'), 
    join_paths(tskit_dir, 'tsk_tables.c'), 
    join_paths(tskit_dir, 'tsk_trees.c')]
    
msprime_sources =[
    'msprime.c', 'fenwick.c', 'util.c', 'mutgen.c', 'object_heap.c',
    'recomb_map.c']

avl_lib = static_library('avl', sources: ['avl.c'])
tskit_lib = static_library('tskit', 
    sources: tskit_sources, dependencies: m_dep, c_args: extra_c_args)
msprime_lib = static_library('msprime', 
    sources: msprime_sources, dependencies: [m_dep, gsl_dep], 
    c_args: extra_c_args, link_with:[avl_lib, tskit_lib])

unit_tests = executable('tests', 
    sources: ['tests/tests.c'], 
    link_with: [msprime_lib], dependencies:[cunit_dep])
test('Unit tests', unit_tests)

# The development CLI. Don't use extra C args because argtable code won't pass
executable('dev-cli', 
    sources: ['dev-tools/dev-cli.c', 'dev-tools/argtable3.c'], 
    link_with: [msprime_lib], dependencies: [config_dep],
    c_args:['-Dlint'])
