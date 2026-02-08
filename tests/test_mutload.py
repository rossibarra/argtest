import os
from pathlib import Path

import numpy as np
import pytest
import tskit

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mutload as ml


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


def test_windowing_sanity():
    ts = make_simple_ts()
    windows = np.array([0, 5, 10], dtype=float)
    load = ml.mutational_load(ts, windows=windows)
    names = ml.sample_names(ts)
    load, _ = ml.aggregate_by_individual(load, names)
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
    remove = ml.load_remove_intervals([bed])
    assert set(remove.keys()) == {"A", "B", "C"}
    assert remove["A"]["starts"] == [1.0]
    assert remove["B"]["starts"] == [1.0]
    assert remove["C"]["starts"] == [5.0]


def test_segment_merge_logic():
    intervals = {
        "A": {"starts": [1.0, 4.0], "ends": [3.0, 6.0]},
        "B": {"starts": [2.0], "ends": [5.0]},
    }
    segs = ml.build_removal_segments(intervals, 10.0)
    # Expect segments covering [0,1), [1,2), [2,3), [3,4), [4,5), [5,6), [6,10)
    assert segs[0][0] == 0.0 and segs[0][1] == 1.0
    assert segs[-1][0] == 6.0 and segs[-1][1] == 10.0


def test_trim_preserves_coordinates(tmp_path):
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = ml.name_to_nodes_map(ts)
    trimmed = ml.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.sequence_length == ts.sequence_length
    assert trimmed.sites_position.min() >= 0.0
    assert trimmed.sites_position.max() <= ts.sequence_length


def test_trim_removes_nodes_in_interval():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = ml.name_to_nodes_map(ts)
    trimmed = ml.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    # Mutation at position 1 on A should be removed; position 7 on B remains
    assert trimmed.num_sites == 1
    assert trimmed.sites_position.tolist() == [7.0]


def test_trim_mutation_parent_integrity():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = ml.name_to_nodes_map(ts)
    trimmed = ml.trim_ts_by_intervals(ts, intervals, name_to_nodes)
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
    monkeypatch.setattr(ml, "load_ts", lambda _: ts)
    monkeypatch.setattr(ml, "parse_args", lambda: type("A", (), {
        "ts": str(ts_path),
        "window_size": 5.0,
        "remove": None,
        "cutoff": 0.5,
        "out": "out.html",
        "suffix_to_strip": "_anchorwave",
    })())
    ml.main()

    assert (cwd / "results" / "out.html").exists()
    assert (cwd / "results" / "test_outliers.bed").exists()
    assert (cwd / "logs" / "out.log").exists()


def test_no_remove_no_trimmed(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ml, "load_ts", lambda _: ts)
    monkeypatch.setattr(ml, "parse_args", lambda: type("A", (), {
        "ts": str(ts_path),
        "window_size": 5.0,
        "remove": None,
        "cutoff": 0.5,
        "out": "out.html",
        "suffix_to_strip": "_anchorwave",
    })())
    ml.main()
    assert not (cwd / "results" / "test_trimmed.tsz").exists()
