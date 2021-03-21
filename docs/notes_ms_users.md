(sec_notes_ms_users)=
# Notes for ms users

Msprime began as an efficient reimplementation of the classical ``ms``
program, and therefore largely follows the same underlying models
as ``ms``. If you wish to use ``msprime`` as a direct replacement
for ``ms``, please use the {ref}`mspms<sec_mspms>` program. This
aims to be 100% compatible with ``ms``, and should be substantially
faster when simulating large regions.

However, simulations using ``mspms`` will still be hampered by the
text output of ``ms``, which is very inefficient for large
simulations. For larger simulations the [tskit](https://tskit.dev/tskit)
output produced by ``msprime``'s Python API will be many times
smaller and faster to process than text based formats.

The Python APIs for {ref}`describing demographic models<sec_demography>`,
running {ref}`ancestry<sec_ancestry>` and
{ref}`mutation<sec_mutations>` simulations are extensively documented.
There are a few basic differences to ``ms`` that are worth keeping
in mind:

- Time is measured in generations rather than coalescent units.
- Population sizes are defined as absolute values, not scaled relative to Ne.
- All rates are absolute rates per generation (i.e., not scaled)
- Recombination rates are per unit of sequence length (not total rates
  along the genome)
- Gene conversion rates are absolute, and do not depend on the recombination
  rate.
- Mutations are at finite discrete sites by default, following the Jukes-Cantor
  nucleotide model (but infinite sites 0/1 mutations at continuous locations
  can be produced, if required).
