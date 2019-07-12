import unittest
import msprime as msp


class TestZeroWaitingTime(unittest.TestCase):
    """
    Test of a rare occurrence, which is a zero
    time to the next common ancestor using the
    Hudson algorithm.  The simulation here
    is based on input from David Lawrie.
    """
    def _ancient_sample_test(
            num_modern=1000, anc_pop=0, anc_num=1, anc_time=200,
            split_time_anc=400, Ne0=10000, Ne1=10000, length=1000):
        samples = [msp.Sample(population=0, time=0)]*num_modern
        samples.extend(
            [msp.Sample(population=anc_pop, time=anc_time)]*(2*anc_num))
        pop_config = [msp.PopulationConfiguration(
            initial_size=Ne0), msp.PopulationConfiguration(initial_size=Ne1)]
        divergence = [msp.MassMigration(
            time=split_time_anc, source=1, destination=0, proportion=1.0)]
        seed = 94320219
        sims = msp.simulate(
            samples=samples, Ne=Ne0, population_configurations=pop_config,
            demographic_events=divergence, length=length,
            random_seed=seed)
        return sims

    def test_zero_waitin_time(self):
        self._ancient_sample_test(num_modern=100, anc_pop=0,
                                  anc_num=10, Ne0=3000,
                                  Ne1=3000, anc_time=5419,
                                  split_time_anc=5919, length=500)


if __name__ == "__main__":
    unittest.main()
