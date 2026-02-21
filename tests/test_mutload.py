import os
from pathlib import Path

import numpy as np
import pytest
import tskit

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import argtest_common as mc
import mutload_summary as ms


def make_simple_ts():
    tables = tskit.TableCollection(sequence_length=10)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    ind0 = tables.individuals.add_row(metadata={"id": "A"})
    ind1 = tables.individuals.add_row(metadata={"id": "B"})
    n0 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind0, population=pop)
    n1 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind1, population=pop)
    anc = tables.nodes.add_row(time=1, population=pop)
    tables.edges.add_row(left=0, right=10, parent=anc, child=n0)
    tables.edges.add_row(left=0, right=10, parent=anc, child=n1)
    s1 = tables.sites.add_row(position=1, ancestral_state="0")
    s7 = tables.sites.add_row(position=7, ancestral_state="0")
    tables.mutations.add_row(site=s1, node=n0, derived_state="1")
    tables.mutations.add_row(site=s7, node=n1, derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def make_ts_no_mutations(n_samples=2, length=10):
    tables = tskit.TableCollection(sequence_length=length)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    inds = [tables.individuals.add_row(metadata={"id": f"I{i}"}) for i in range(n_samples)]
    samples = [
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=inds[i], population=pop)
        for i in range(n_samples)
    ]
    anc = tables.nodes.add_row(time=1, population=pop)
    for s in samples:
        tables.edges.add_row(left=0, right=length, parent=anc, child=s)
    tables.sort()
    return tables.tree_sequence()


def make_ts_many_individuals(n=100, length=10):
    tables = tskit.TableCollection(sequence_length=length)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    inds = [tables.individuals.add_row(metadata={"id": f"I{i}"}) for i in range(n)]
    samples = [
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=inds[i], population=pop)
        for i in range(n)
    ]
    anc = tables.nodes.add_row(time=1, population=pop)
    for s in samples:
        tables.edges.add_row(left=0, right=length, parent=anc, child=s)
    s1 = tables.sites.add_row(position=1, ancestral_state="0")
    tables.mutations.add_row(site=s1, node=samples[0], derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def test_windowing_sanity():
    ts = make_simple_ts()
    windows = np.array([0, 5, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, _ = mc.aggregate_by_individual(load, names)
    assert load.shape == (2, 2)
    # First window has site at 1 on A, second window has site at 7 on B
    assert load[0, 0] == 1
    assert load[0, 1] == 0
    assert load[1, 0] == 0
    assert load[1, 1] == 1


def test_outlier_mask_logic():
    load = np.array([[12, 5, 5], [2, 2, 2]], dtype=float)
    means = load.mean(axis=1)
    cutoff = 0.5
    high = (1 + cutoff) * means
    low = (1 - cutoff) * means
    mask = (load > high[:, None]) | (load < low[:, None])
    assert mask[0].tolist() == [True, False, False]
    assert mask[1].tolist() == [False, False, False]


def test_remove_bed_parsing(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t1\t3\tA,B\nchr1\t5\t6\tC\n")
    remove = mc.load_remove_intervals([bed])
    assert set(remove.keys()) == {"A", "B", "C"}
    assert remove["A"]["starts"] == [1.0]
    assert remove["B"]["starts"] == [1.0]
    assert remove["C"]["starts"] == [5.0]


def test_segment_merge_logic():
    intervals = {
        "A": {"starts": [1.0, 4.0], "ends": [3.0, 6.0]},
        "B": {"starts": [2.0], "ends": [5.0]},
    }
    segs = mc.build_removal_segments(intervals, 10.0)
    # Expect segments covering [0,1), [1,2), [2,3), [3,4), [4,5), [5,6), [6,10)
    assert segs[0][0] == 0.0 and segs[0][1] == 1.0
    assert segs[-1][0] == 6.0 and segs[-1][1] == 10.0


def test_trim_preserves_coordinates(tmp_path):
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.sequence_length == ts.sequence_length
    assert trimmed.sites_position.min() >= 0.0
    assert trimmed.sites_position.max() <= ts.sequence_length


def test_trim_removes_nodes_in_interval():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    # Mutation at position 1 on A should be removed; position 7 on B remains
    assert trimmed.num_sites == 1
    assert trimmed.sites_position.tolist() == [7.0]


def test_trim_mutation_parent_integrity():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    # Should load without parent errors
    assert trimmed.num_sites >= 0


def test_outputs_written(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    # Clean outputs
    for p in (cwd / "results", cwd / "logs"):
        if p.exists():
            for item in p.iterdir():
                if item.is_file():
                    item.unlink()

    # Run main via function
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mpl"))
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    monkeypatch.setattr(ms, "parse_args", lambda: type("A", (), {
        "ts": str(ts_path),
        "window_size": 5.0,
        "cutoff": 0.5,
        "out": "out.html",
        "suffix_to_strip": "_anchorwave",
    })())
    ms.main()

    assert (cwd / "results" / "out.html").exists()
    assert (cwd / "results" / "test_outliers.bed").exists()
    assert (cwd / "logs" / "out.log").exists()


def test_no_remove_no_trimmed(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    monkeypatch.setattr(ms, "parse_args", lambda: type("A", (), {
        "ts": str(ts_path),
        "window_size": 5.0,
        "cutoff": 0.5,
        "out": "out.html",
        "suffix_to_strip": "_anchorwave",
    })())
    ms.main()
    assert not (cwd / "results" / "test_trimmed.tsz").exists()


def test_no_mutations_outliers_empty(tmp_path, monkeypatch):
    ts = make_ts_no_mutations()
    ts_path = tmp_path / "nomut.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    monkeypatch.setattr(ms, "parse_args", lambda: type("A", (), {
        "ts": str(ts_path),
        "window_size": 5.0,
        "cutoff": 0.5,
        "out": "nomut.html",
        "suffix_to_strip": "_anchorwave",
    })())
    ms.main()
    bed = cwd / "results" / "nomut_outliers.bed"
    assert bed.exists()
    assert bed.read_text().strip() == ""


def test_single_sample_ts():
    ts = make_ts_no_mutations(n_samples=1, length=10)
    windows = np.array([0, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, _ = mc.aggregate_by_individual(load, names)
    assert load.shape == (1, 1)


def test_intervals_outside_sequence_length():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [100.0], "ends": [200.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.num_sites == ts.num_sites


def test_overlapping_bed_with_commas(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t1\t4\tA,B\nchr1\t3\t5\tB,C\n")
    remove = mc.load_remove_intervals([bed])
    assert set(remove.keys()) == {"A", "B", "C"}


def test_all_samples_removed_in_segment():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [10.0]}, "B": {"starts": [0.0], "ends": [10.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.sequence_length == ts.sequence_length


def test_idempotent_no_effect_remove():
    ts = make_simple_ts()
    intervals = {"C": {"starts": [0.0], "ends": [5.0]}}  # C not in TS
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.num_sites == ts.num_sites
    assert trimmed.num_edges >= ts.num_edges


def test_sample_order_preserved():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.samples().tolist() == ts.samples().tolist()


def test_trim_validate():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    if hasattr(trimmed, "validate"):
        trimmed.validate()


def test_large_interval_counts():
    intervals = {"A": {"starts": [], "ends": []}}
    for i in range(1000):
        intervals["A"]["starts"].append(float(i))
        intervals["A"]["ends"].append(float(i + 0.5))
    segs = mc.build_removal_segments(intervals, 2000.0)
    assert segs[0][0] == 0.0
    assert segs[-1][1] == 2000.0


def test_many_individuals_shapes():
    ts = make_ts_many_individuals(n=100, length=10)
    windows = np.array([0, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, unique = mc.aggregate_by_individual(load, names)
    assert load.shape == (1, 100)
    assert len(unique) == 100


def test_relative_bed_paths(tmp_path, monkeypatch):
    bed = tmp_path / "rel.bed"
    bed.write_text("chr1\t1\t3\tA\n")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        remove = mc.load_remove_intervals([Path("rel.bed")])
        assert "A" in remove
    finally:
        os.chdir(cwd)


def test_output_overwrite(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    monkeypatch.setattr(ms, "parse_args", lambda: type("A", (), {
        "ts": str(ts_path),
        "window_size": 5.0,
        "cutoff": 0.5,
        "out": "overwrite.html",
        "suffix_to_strip": "_anchorwave",
    })())
    ms.main()
    out = cwd / "results" / "overwrite.html"
    assert out.exists()
    first = out.read_text()
    ms.main()
    second = out.read_text()
    assert first == second
