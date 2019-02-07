project('msprime', 'c')

kastore_proj = subproject('kastore')
kastore_dep = kastore_proj.get_variable('kastore_dep')
tskit_proj = subproject('tskit')
tskit_dep = tskit_proj.get_variable('tskit_dep')

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


# add_global_arguments(['-I' + kastore_dir, '-I' + tskit_dir], language: 'c')
    
msprime_sources =[
    'msprime.c', 'fenwick.c', 'util.c', 'mutgen.c', 'object_heap.c',
    'recomb_map.c']

avl_lib = static_library('avl', sources: ['avl.c'])
msprime_lib = static_library('msprime', 
    sources: msprime_sources, dependencies: [m_dep, gsl_dep, kastore_dep, tskit_dep], 
    c_args: extra_c_args, link_with:[avl_lib])

unit_tests = executable('tests', 
    sources: ['tests/tests.c'], 
    link_with: [msprime_lib], dependencies:[cunit_dep, kastore_dep, tskit_dep])
test('Unit tests', unit_tests)

# The development CLI. Don't use extra C args because argtable code won't pass
executable('dev-cli', 
    sources: ['dev-tools/dev-cli.c', 'dev-tools/argtable3.c'], 
    link_with: [msprime_lib], dependencies: [config_dep, kastore_dep, tskit_dep],
    c_args:['-Dlint'])
