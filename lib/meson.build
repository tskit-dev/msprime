project('msprime', 'c')

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
    
msprime_sources =[
    'msprime.c', 'fenwick.c', 'util.c', 'mutgen.c', 'object_heap.c',
    'likelihood.c', 'rate_map.c']

avl_lib = static_library('avl', sources: ['avl.c'])
msprime_lib = static_library('msprime', 
    sources: msprime_sources, dependencies: [m_dep, gsl_dep, tskit_dep], 
    c_args: extra_c_args, link_with:[avl_lib])

# Unit tests
test_lib = static_library('testlib',
    sources: ['tests/testlib.c'], dependencies:[cunit_dep, tskit_dep])

test_core = executable('test_core',
    sources: ['tests/test_core.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('core', test_core)

test_ancestry = executable('test_ancestry',
    sources: ['tests/test_ancestry.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('ancestry', test_ancestry)

test_pedigrees = executable('test_pedigrees',
    sources: ['tests/test_pedigrees.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('pedigrees', test_pedigrees)

test_mutations = executable('test_mutations',
    sources: ['tests/test_mutations.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('mutations', test_mutations)

test_likelihood = executable('test_likelihood',
    sources: ['tests/test_likelihood.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('likelihood', test_likelihood)

test_fenwick = executable('test_fenwick',
    sources: ['tests/test_fenwick.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('fenwick', test_fenwick)

test_rate_map = executable('test_rate_map',
    sources: ['tests/test_rate_map.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('rate_map', test_rate_map)

test_sweeps = executable('test_sweeps',
    sources: ['tests/test_sweeps.c'], 
    link_with: [msprime_lib, test_lib], dependencies: [cunit_dep, tskit_dep])
test('sweeps', test_sweeps)

# The development CLI. Don't use extra C args because argtable code won't pass
executable('dev-cli', 
    sources: ['dev-tools/dev-cli.c', 'dev-tools/argtable3.c'], 
    link_with: [msprime_lib], dependencies: [config_dep, tskit_dep],
    c_args:['-Dlint'])
