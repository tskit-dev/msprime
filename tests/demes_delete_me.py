import collections
import itertools
import math

import demes
import hypothesis as hyp
import hypothesis.strategies as st
import numpy as np


@st.composite
def deme_names(draw, max_length=20):
    """
    A hypothesis strategy for creating a valid deme name.
    """
    name = draw(st.text(min_size=1, max_size=max_length))
    # Names must be valid Python identifiers.
    hyp.assume(name.isidentifier())
    return name


@st.composite
def yaml_strings(draw, min_size=1, max_size=100):
    """
    A hypothesis strategy for creating a valid YAML string.

    From https://yaml.org/spec/1.2/spec.html#id2770814

        To ensure readability, YAML streams use only the printable subset of
        the Unicode character set. The allowed character range explicitly
        excludes the C0 control block #x0-#x1F (except for TAB #x9, LF #xA,
        and CR #xD which are allowed), DEL #x7F, the C1 control block #x80-#x9F
        (except for NEL #x85 which is allowed), the surrogate block #xD800-#xDFFF,
        #xFFFE, and #xFFFF.

        On input, a YAML processor must accept all Unicode characters except
        those explicitly excluded above.

        On output, a YAML processor must only produce acceptable characters.
        Any excluded characters must be presented using escape sequences.
        In addition, any allowed characters known to be non-printable should
        also be escaped. This isn’t mandatory since a full implementation would
        require extensive character property tables.
    """
    return draw(
        st.text(
            alphabet=st.characters(
                blacklist_categories=(
                    "Cc",  # control block (C0 and C1)
                    "Cs",  # surrogate block
                ),
                blacklist_characters=("\ufffe", "\uffff"),
                whitelist_characters=("\x09", "\x0a", "\x0d", "\x85"),
            ),
            min_size=min_size,
            max_size=max_size,
        )
    )


@st.composite
def epochs_lists(draw, start_time=math.inf, max_epochs=5):
    """
    A hypothesis strategy for creating lists of Epochs for a deme.

    :param float start_time: The start time of the deme.
    :param int max_epochs: The maximum number of epochs in the list.
    """
    assert max_epochs >= 2
    times = draw(
        st.lists(
            st.floats(min_value=0, max_value=start_time, exclude_max=True),
            unique=True,
            min_size=1,
            max_size=max_epochs,
        )
    )
    times.sort(reverse=True)
    epochs = []

    for i, end_time in enumerate(times):
        start_size = draw(st.floats(min_value=1, allow_infinity=False))
        if i == 0 and math.isinf(start_time):
            end_size = start_size
        else:
            end_size = draw(st.floats(min_value=1, allow_infinity=False))
        cloning_rate = draw(st.floats(min_value=0, max_value=1))
        selfing_rate = draw(st.floats(min_value=0, max_value=1 - cloning_rate))

        epochs.append(
            dict(
                end_time=end_time,
                start_size=start_size,
                end_size=end_size,
                cloning_rate=cloning_rate,
                selfing_rate=selfing_rate,
            )
        )

    return epochs


@st.composite
def migration_matrices(
    draw, graph, max_migrations=10, max_additional_migration_intervals=5
):
    """
    A hypothesis strategy for creating migration matrices for a graph.
    """
    n = len(graph.demes)
    assert n > 0

    uniq_deme_times = {deme.start_time for deme in graph.demes}
    uniq_deme_times.update(deme.end_time for deme in graph.demes)
    start_time, *end_times = sorted(uniq_deme_times, reverse=True)

    # Identify the first time at which 2 or more demes exist simultaneously.
    for end_time in end_times:
        if sum(1 for deme in graph.demes if deme.start_time <= start_time) > 1:
            break
        start_time = end_time

    if start_time == end_times[-1]:
        # No two demes exist simultaneously.
        return [[[0] * n for _ in range(n)]], math.inf, [0]

    saved_start_time = start_time

    # Partition time intervals even further.
    additional_times = draw(
        st.lists(
            st.floats(min_value=end_times[-1], max_value=start_time, exclude_max=True),
            unique=True,
            min_size=0,
            max_size=max_additional_migration_intervals,
        )
    )
    end_times = sorted(set(end_times + additional_times), reverse=True)

    mm_list = [[[0] * n for _ in range(n)] for _ in range(len(end_times))]
    n_migrations = draw(st.integers(min_value=0, max_value=max_migrations))

    for migration_matrix, end_time in zip(mm_list, end_times):
        # Find demes alive in this interval.
        deme_indices = [
            j
            for j, deme in enumerate(graph.demes)
            if (
                deme.start_time >= start_time > deme.end_time
                and deme.start_time > end_time >= deme.end_time
            )
        ]
        if len(deme_indices) < 2:
            continue

        # Select pairs of demes for migration.
        pairs = list(itertools.permutations(deme_indices, 2))
        pair_indices = draw(
            st.lists(
                st.integers(min_value=0, max_value=len(pairs) - 1),
                unique=True,
                min_size=0,
                max_size=min(len(pairs), n_migrations),
            )
        )

        for k in pair_indices:
            a, b = pairs[k]
            assert migration_matrix[a][b] == 0
            max_rate = 1 - sum(migration_matrix[a])
            if math.isclose(max_rate, 0):
                continue
            n_migrations -= 1
            rate = draw(st.floats(min_value=0, max_value=max_rate, exclude_min=True))
            migration_matrix[a][b] = rate

        if n_migrations == 0:
            break
        start_time = end_time

    return mm_list, saved_start_time, end_times


@st.composite
def migrations_lists(draw, graph, max_migrations=10):
    """
    A hypothesis strategy for creating a migration list for a graph.
    """
    mm_list, start_time, end_times = draw(
        migration_matrices(graph, max_migrations=max_migrations)
    )
    assert len(mm_list) == len(end_times)
    migrations = []
    for migration_matrix, end_time in zip(mm_list, end_times):
        for j, row in enumerate(migration_matrix):
            for k, rate in enumerate(row):
                if rate > 0:
                    migration = demes.AsymmetricMigration(
                        source=graph.demes[k].name,
                        dest=graph.demes[j].name,
                        start_time=start_time,
                        end_time=end_time,
                        rate=rate,
                    )
                    migrations.append(migration)
        start_time = end_time
    return migrations


@st.composite
def pulses_lists(draw, graph, max_pulses=10):
    """
    A hypothesis strategy for creating a pulses list for a graph.
    """
    n_pulses = draw(st.integers(min_value=0, max_value=max_pulses))
    pulses = []
    ingress_proportions = collections.defaultdict(lambda: 0)
    for j, deme_j in enumerate(graph.demes[:-1]):
        for deme_k in graph.demes[j + 1 :]:
            time_lo = max(deme_j.end_time, deme_k.end_time)
            time_hi = min(deme_j.start_time, deme_k.start_time)

            # We wish to draw times for the pulses. They must be in the open
            # interval (time_lo, time_hi) to ensure the pulse doesn't happen
            # at any deme's start_time or end_time, which could be invalid.
            # So we check there is at least one floating point number between
            # time_lo and time_hi.
            if time_hi <= np.nextafter(time_lo, np.inf, dtype=float):
                continue
            n = draw(st.integers(min_value=0, max_value=n_pulses))
            for _ in range(n):
                source, dest = deme_j.name, deme_k.name
                if draw(st.booleans()):
                    source, dest = dest, source
                time = draw(
                    st.floats(
                        min_value=time_lo,
                        max_value=time_hi,
                        exclude_min=True,
                        exclude_max=True,
                    )
                )
                max_proportion = 1 - ingress_proportions[(dest, time)]
                if math.isclose(max_proportion, 0):
                    continue
                proportion = draw(
                    st.floats(
                        min_value=0,
                        max_value=max_proportion,
                        exclude_min=True,
                        exclude_max=True,
                    )
                )
                ingress_proportions[(dest, time)] += proportion
                pulse = dict(
                    source=source,
                    dest=dest,
                    time=time,
                    proportion=proportion,
                )
                pulses.append(pulse)
                n_pulses -= 1
            if n_pulses == 0:
                break
        if n_pulses == 0:
            break
    return pulses


@st.composite
def graphs(draw, max_demes=5, max_epochs=10, max_migrations=10, max_pulses=10):
    """
    A hypothesis strategy for creating a Graph.

    .. code-block::

        @hypothesis.given(graphs())
        def test_something(self, graph: demes.Graph):
            # Do something with the ``graph``.
            ...

    :param int max_demes: The maximum number of demes in the graph.
    :param int max_epochs: The maximum number of epochs per deme.
    :param int max_migrations: The maximum number of migrations in the graph.
    :param int max_pulses: The maximum number of pulses in the graph.
    """
    generation_time = draw(st.none() | st.floats(min_value=1e-9, max_value=1e6))
    # JK: setting time_units to always be generations to avoid rounding problems
    # with Graph.in_generations() at the very extreme floating point values.
    # if generation_time is None:
    #     time_units = "generations"
    # else:
    #     time_units = draw(yaml_strings())
    time_units = "generations"
    b = demes.Builder(
        description=draw(yaml_strings()),
        generation_time=generation_time,
        time_units=time_units,
        doi=draw(st.lists(yaml_strings(), max_size=3)),
    )

    for deme_name in draw(st.sets(deme_names(), min_size=1, max_size=max_demes)):
        ancestors = []
        proportions = []
        start_time = math.inf
        n_demes = 0 if "demes" not in b.data else len(b.data["demes"])
        if n_demes > 0:
            # draw indices into demes list to use as ancestors
            anc_idx = draw(
                st.lists(
                    st.integers(min_value=0, max_value=n_demes - 1),
                    unique=True,
                    max_size=n_demes,
                )
            )
            if len(anc_idx) > 0:
                time_hi = min(b.data["demes"][j]["start_time"] for j in anc_idx)
                time_lo = max(
                    b.data["demes"][j]["epochs"][-1]["end_time"] for j in anc_idx
                )
                # If time_hi > time_lo, the proposed ancestors exist
                # at the same time. So we draw a number for the deme's
                # start_time, which must be in the half-open interval
                # [time_lo, time_hi), with the further constraint that the
                # start_time cannot be 0.
                # However, there may not be any floating point numbers between
                # 0 and time_hi even if time_hi > 0, so we check that time_hi
                # is greater than the smallest positive number.
                if (time_lo > 0 and time_hi > time_lo) or (
                    time_lo == 0 and time_hi > np.finfo(float).tiny
                ):
                    # Draw a start time and the ancestry proportions.
                    start_time = draw(
                        st.floats(
                            min_value=time_lo,
                            max_value=time_hi,
                            exclude_max=True,
                            # Can't have start_time=0.
                            exclude_min=time_lo == 0,
                        )
                    )
                    ancestors = [b.data["demes"][j]["name"] for j in anc_idx]
                    if len(ancestors) == 1:
                        proportions = [1.0]
                    else:
                        proportions = draw(
                            st.lists(
                                st.integers(min_value=1, max_value=10 ** 9),
                                min_size=len(ancestors),
                                max_size=len(ancestors),
                            )
                        )
                        psum = sum(proportions)
                        proportions = [p / psum for p in proportions]
        b.add_deme(
            name=deme_name,
            description=draw(st.none() | yaml_strings()),
            ancestors=ancestors,
            proportions=proportions,
            epochs=draw(epochs_lists(start_time=start_time, max_epochs=max_epochs)),
            start_time=start_time,
        )

    graph = b.resolve()
    graph.migrations = draw(migrations_lists(graph, max_migrations=max_migrations))
    graph.pulses = draw(pulses_lists(graph, max_pulses=max_pulses))
    # Resolve the graph again. This is not strictly necessary, but has only
    # a small computational overhead and serves to catch simple errors in
    # the migrations_lists()/pulses_lists() implementations.
    graph = demes.Builder.fromdict(graph.asdict()).resolve()
    return graph
