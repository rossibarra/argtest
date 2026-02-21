#!/usr/bin/env python
from __future__ import annotations

import sys
from bisect import bisect_right
from pathlib import Path

import numpy as np
import tskit

try:
    import tszip
except Exception:  # pragma: no cover - optional dependency
    tszip = None

try:
    import msprime
except Exception:  # pragma: no cover - optional dependency
    msprime = None


def _require_msprime():
    if msprime is None:
        raise RuntimeError("msprime is required for accessibility mask and collapse helpers")


def load_ts(path: Path) -> tskit.TreeSequence:
    # Handle compressed tree sequences if requested.
    if path.suffix == ".tsz":
        if tszip is None:
            raise RuntimeError("tszip is required to load .tsz files")
        return tszip.load(str(path))
    return tskit.load(str(path))


def dump_ts(ts: tskit.TreeSequence, out_path: Path) -> None:
    # Mirror input format when writing if extension is .tsz.
    if out_path.suffix == ".tsz":
        if tszip is None:
            raise RuntimeError("tszip is required to write .tsz files")
        tszip.compress(ts, out_path)
        return
    ts.dump(str(out_path))


def get_individual_name(ind, suffix_to_strip="_anchorwave") -> str:
    # Prefer metadata id when present; otherwise fall back to a stable synthetic name.
    nm = None
    try:
        if isinstance(ind.metadata, dict):
            nm = ind.metadata.get("id")
    except Exception:
        nm = None
    if nm is None:
        nm = f"ind{ind.id}"
    if isinstance(nm, bytes):
        nm = nm.decode()
    return nm.replace(suffix_to_strip, "")


def sample_names(ts: tskit.TreeSequence, suffix_to_strip="_anchorwave"):
    # Map each sample node to its individual name (or node id if missing).
    names = []
    for u in ts.samples():
        ind_id = ts.node(u).individual
        if ind_id != tskit.NULL:
            nm = get_individual_name(ts.individual(ind_id), suffix_to_strip=suffix_to_strip)
        else:
            nm = f"node{u}"
        names.append(nm)
    return names


def aggregate_by_individual(load, names):
    # Collapse per-sample loads into per-individual loads by name.
    unique = []
    idx_map = {}
    for i, nm in enumerate(names):
        if nm not in idx_map:
            idx_map[nm] = len(unique)
            unique.append(nm)

    if load.ndim == 1:
        agg = np.zeros(len(unique), dtype=float)
        for i, nm in enumerate(names):
            agg[idx_map[nm]] += load[i]
        return agg, unique

    agg = np.zeros((load.shape[0], len(unique)), dtype=float)
    for i, nm in enumerate(names):
        agg[:, idx_map[nm]] += load[:, i]
    return agg, unique


def name_to_nodes_map(ts: tskit.TreeSequence, suffix_to_strip="_anchorwave"):
    # Build lookup from individual name to all node ids for that individual.
    mapping = {}
    for ind_id, ind in enumerate(ts.individuals()):
        nm = get_individual_name(ind, suffix_to_strip=suffix_to_strip)
        nodes = list(ts.individual(ind_id).nodes)
        mapping[nm] = nodes
    return mapping


def mutational_load(
    ts: tskit.TreeSequence,
    windows: np.ndarray | None = None,
    remove_intervals=None,
    name_to_nodes=None,
) -> np.ndarray:
    # Compute derived mutation counts per sample, optionally per window and with removals.
    genome_windows = np.array([0, ts.sequence_length]) if windows is None else windows
    assert genome_windows[0] == 0 and genome_windows[-1] == ts.sequence_length
    load = np.zeros((genome_windows.size - 1, ts.num_samples))

    segments = [(0.0, ts.sequence_length, frozenset())]
    if remove_intervals and name_to_nodes:
        # Convert BED intervals into contiguous genome segments with active drop sets.
        segments = build_removal_segments(remove_intervals, ts.sequence_length)

    for left, right, drop_names in segments:
        if right <= left:
            continue
        drop_nodes = set()
        for nm in drop_names:
            drop_nodes.update(name_to_nodes.get(nm, []))

        sub = ts.keep_intervals([(left, right)], simplify=False)
        if drop_nodes:
            # Simplify to retained samples but preserve original ids for load accumulation.
            keep = [u for u in ts.samples() if u not in drop_nodes]
            sub, node_map = sub.simplify(samples=keep, keep_unary=True, map_nodes=True)
            rev_map = np.full(sub.num_nodes, tskit.NULL, dtype=int)
            for u in keep:
                new_id = node_map[u]
                if new_id != tskit.NULL:
                    rev_map[new_id] = u
        else:
            rev_map = None

        site_windows = None
        for tree in sub.trees(sample_lists=True):
            for s in tree.sites():
                if site_windows is None:
                    # Cache site->window indices once per sub-TS.
                    site_windows = np.digitize(sub.sites_position, genome_windows) - 1
                window = int(site_windows[s.id])
                if window < 0 or window >= genome_windows.size - 1:
                    continue
                for m in s.mutations:
                    if m.edge != tskit.NULL:
                        samples = list(tree.samples(m.node))
                        if rev_map is not None:
                            samples = [rev_map[u] for u in samples if rev_map[u] != tskit.NULL]
                        load[window, samples] += 1.0
    return load.squeeze(0) if windows is None else load


def build_removal_segments(remove_intervals, sequence_length):
    # Sweep line to build segments where a set of names is "dropped".
    events = {}
    for name, intervals in remove_intervals.items():
        for start, end in zip(intervals["starts"], intervals["ends"]):
            events.setdefault(start, []).append((name, 1))
            events.setdefault(end, []).append((name, -1))

    positions = sorted(set([0.0, float(sequence_length)] + list(events.keys())))
    active = {}
    segments = []
    for i, pos in enumerate(positions[:-1]):
        for name, delta in events.get(pos, []):
            active[name] = active.get(name, 0) + delta
            if active[name] <= 0:
                active.pop(name, None)
        next_pos = positions[i + 1]
        if next_pos > pos:
            segments.append((pos, next_pos, frozenset(active.keys())))
    return segments


def build_segments_with_drop_nodes(remove_intervals, name_to_nodes, sequence_length):
    # Convert name-based segments into node-id drop sets for trimming.
    segments = build_removal_segments(remove_intervals, sequence_length)
    out = []
    for left, right, drop_names in segments:
        drop_nodes = set()
        for nm in drop_names:
            drop_nodes.update(name_to_nodes.get(nm, []))
        out.append((left, right, drop_nodes))
    return out


def trim_ts_by_intervals(ts: tskit.TreeSequence, remove_intervals, name_to_nodes):
    # Remove edges and mutations that overlap dropped individuals in each segment.
    segments = build_segments_with_drop_nodes(remove_intervals, name_to_nodes, ts.sequence_length)
    seg_lefts = [s[0] for s in segments]
    seg_rights = [s[1] for s in segments]

    tables = ts.dump_tables()
    tables.edges.clear()
    tables.sites.clear()
    tables.mutations.clear()

    # Edges: split by segments, remove if parent/child dropped in that segment
    for edge in ts.edges():
        i = bisect_right(seg_rights, edge.left)
        while i < len(segments) and seg_lefts[i] < edge.right:
            left, right, drop_nodes = segments[i]
            seg_l = max(edge.left, left)
            seg_r = min(edge.right, right)
            if seg_r > seg_l:
                if edge.parent not in drop_nodes and edge.child not in drop_nodes:
                    tables.edges.add_row(
                        left=seg_l,
                        right=seg_r,
                        parent=edge.parent,
                        child=edge.child,
                        metadata=edge.metadata,
                    )
            i += 1

    # Sites and mutations: drop mutations whose node is removed in the segment
    for site in ts.sites():
        i = bisect_right(seg_rights, site.position)
        if i >= len(segments) or seg_lefts[i] > site.position:
            i = bisect_right(seg_lefts, site.position) - 1
        if i < 0 or i >= len(segments):
            continue
        _, _, drop_nodes = segments[i]
        kept = []
        for mut in site.mutations:
            if mut.node in drop_nodes:
                continue
            kept.append(mut)
        if not kept:
            continue
        new_site_id = tables.sites.add_row(
            position=site.position,
            ancestral_state=site.ancestral_state,
            metadata=site.metadata,
        )
        for mut in kept:
            tables.mutations.add_row(
                site=new_site_id,
                node=mut.node,
                derived_state=mut.derived_state,
                parent=tskit.NULL,
                time=mut.time,
                metadata=mut.metadata,
            )

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def assert_sample_ids_preserved(orig_ts: tskit.TreeSequence, trimmed_ts: tskit.TreeSequence):
    # Ensure we only removed edges/mutations, not samples or ids.
    orig_samples = orig_ts.samples()
    trimmed_samples = trimmed_ts.samples()
    if len(orig_samples) != len(trimmed_samples):
        msg = "Trimmed tree sequence changed the number of samples"
        print(f"ERROR: {msg}", file=sys.stderr)
        raise RuntimeError(msg)
    if not np.array_equal(orig_samples, trimmed_samples):
        msg = "Trimmed tree sequence changed sample IDs/order"
        print(f"ERROR: {msg}", file=sys.stderr)
        raise RuntimeError(msg)


def validate_trimmed_ts(ts: tskit.TreeSequence):
    # Ensure internal indexes and topology are consistent
    if hasattr(ts, "check_index"):
        try:
            ts.check_index()
        except TypeError:
            # Older/newer tskit versions require explicit index/length args.
            index = None
            if hasattr(ts, "index"):
                index = ts.index
            elif hasattr(ts, "tables") and hasattr(ts.tables, "index"):
                index = ts.tables.index
            if index is None:
                print("WARNING: trimmed ts skipped check_index (no index available)", file=sys.stderr)
                return
            try:
                ts.check_index(index, ts.sequence_length)
            except Exception as e:
                print(f"ERROR: trimmed ts failed check_index: {e}", file=sys.stderr)
                raise
        except Exception as e:
            print(f"ERROR: trimmed ts failed check_index: {e}", file=sys.stderr)
            raise
    if hasattr(ts, "validate"):
        try:
            ts.validate()
        except Exception as e:
            print(f"ERROR: trimmed ts failed validate: {e}", file=sys.stderr)
            raise


def load_remove_intervals(paths):
    # Load BED intervals into per-name start/end arrays.
    remove = {}
    for path in paths:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Remove BED not found: {p}")
        with open(p, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    raise ValueError(f"Invalid BED line in {p}: {line}")
                start = float(parts[1])
                end = float(parts[2])
                if end <= start:
                    continue
                if len(parts) >= 4:
                    raw = parts[3]
                    names = [n.strip() for n in raw.split(",") if n.strip()]
                else:
                    names = [p.stem]
                for name in names:
                    remove.setdefault(name, []).append((start, end))

    intervals = {}
    for name, spans in remove.items():
        spans.sort()
        starts = [s for s, _ in spans]
        ends = [e for _, e in spans]
        intervals[name] = {"starts": starts, "ends": ends}
    return intervals


def merge_intervals(intervals):
    # Merge overlapping or adjacent half-open intervals [left, right).
    if len(intervals) == 0:
        return []
    intervals = np.asarray(intervals, dtype=float)
    intervals = intervals[np.argsort(intervals[:, 0])]
    merged = [intervals[0].tolist()]
    for left, right in intervals[1:]:
        if left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return merged


def accessible_intervals_from_mu(mu):
    # Convert mutation-rate map into explicit accessible intervals.
    pos = np.asarray(mu.position, dtype=float)
    rate = np.asarray(mu.rate)
    keep = rate > 0
    lefts = pos[:-1][keep]
    rights = pos[1:][keep]
    return np.column_stack([lefts, rights])


def overlap_lengths(intervals, windows):
    # Return total overlap length per window for sorted half-open intervals.
    intervals = np.asarray(intervals, dtype=float)
    windows = np.asarray(windows, dtype=float)
    k = len(windows) - 1
    out = np.zeros(k, dtype=float)
    i = 0
    m = len(intervals)
    for w in range(k):
        wl, wr = windows[w], windows[w + 1]
        while i < m and intervals[i, 1] <= wl:
            i += 1
        j = i
        while j < m and intervals[j, 0] < wr:
            ol = max(wl, intervals[j, 0])
            or_ = min(wr, intervals[j, 1])
            if or_ > ol:
                out[w] += or_ - ol
            if intervals[j, 1] >= wr:
                break
            j += 1
    return out


def ratemap_from_keep_intervals(keep_intervals, sequence_length):
    # Build a binary RateMap over [0, sequence_length] with 1 in keep intervals.
    _require_msprime()
    keep = sorted(
        [
            [max(0.0, float(left)), min(float(sequence_length), float(right))]
            for left, right in keep_intervals
            if right > left
        ],
        key=lambda x: x[0],
    )
    breaks = [0.0]
    for left, right in keep:
        breaks.extend([left, right])
    breaks.append(float(sequence_length))
    pos = np.unique(np.asarray(breaks, dtype=float))
    pos.sort()
    rate = np.zeros(len(pos) - 1, dtype=float)
    k = 0
    for i in range(len(rate)):
        seg_l, seg_r = pos[i], pos[i + 1]
        while k < len(keep) and keep[k][1] <= seg_l:
            k += 1
        if k < len(keep) and keep[k][0] < seg_r and keep[k][1] > seg_l:
            rate[i] = 1.0
    return msprime.RateMap(position=pos, rate=rate)


def and_ratemaps_binary(a: msprime.RateMap, b: msprime.RateMap):
    # Return binary RateMap that is 1 where both inputs are 1.
    _require_msprime()
    assert a.sequence_length == b.sequence_length
    pos = np.unique(np.concatenate([np.asarray(a.position), np.asarray(b.position)]))
    pos.sort()
    rate = (a.get_rate(pos[:-1]).astype(bool) & b.get_rate(pos[:-1]).astype(bool)).astype(float)
    return msprime.RateMap(position=pos, rate=rate)


def collapse_masked_intervals(ts: tskit.TreeSequence, accessible: msprime.RateMap):
    # Collapse masked intervals and compact coordinates to unmasked sequence length.
    assert np.all(np.logical_or(accessible.rate == 0.0, accessible.rate == 1.0))
    assert accessible.sequence_length == ts.sequence_length
    tables = ts.dump_tables()
    tables.sequence_length = accessible.get_cumulative_mass(ts.sequence_length)
    tables.edges.left = accessible.get_cumulative_mass(tables.edges.left)
    tables.edges.right = accessible.get_cumulative_mass(tables.edges.right)
    tables.edges.keep_rows(tables.edges.right > tables.edges.left)
    is_connected = np.full(tables.nodes.num_rows, False)
    is_connected[tables.edges.parent] = True
    is_connected[tables.edges.child] = True
    node_map = tables.nodes.keep_rows(is_connected)
    tables.edges.parent = node_map[tables.edges.parent]
    tables.edges.child = node_map[tables.edges.child]
    site_map = tables.sites.keep_rows(accessible.get_rate(tables.sites.position).astype(bool))
    tables.sites.position = accessible.get_cumulative_mass(tables.sites.position)
    tables.mutations.node = node_map[tables.mutations.node]
    tables.mutations.site = site_map[tables.mutations.site]
    tables.mutations.keep_rows(
        np.logical_and(
            tables.mutations.site != tskit.NULL,
            tables.mutations.node != tskit.NULL,
        )
    )
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def collapse_masked_and_low_access_windows(ts, mu, window_size, cutoff_bp):
    # Build combined accessibility mask and collapse TS in one step.
    mask, good_intervals, windows, acc_bp = build_shared_mask(
        ts=ts,
        mu=mu,
        window_size=window_size,
        cutoff_bp=cutoff_bp,
    )
    ts2 = collapse_masked_intervals(ts, mask)
    return ts2, mask, good_intervals, windows, acc_bp


def build_shared_mask(ts, mu, window_size, cutoff_bp):
    # Compute the binary mask from mutation-map accessibility and window accessibility.
    _require_msprime()
    sequence_length = float(ts.sequence_length)
    acc_mu = msprime.RateMap(
        position=np.asarray(mu.position, dtype=float),
        rate=(np.asarray(mu.rate) > 0).astype(float),
    )
    assert acc_mu.sequence_length == sequence_length
    acc_intervals = accessible_intervals_from_mu(mu)
    windows = np.arange(0, sequence_length + window_size, window_size, dtype=float)
    if windows[-1] > sequence_length:
        windows[-1] = sequence_length
    acc_bp = overlap_lengths(acc_intervals, windows)
    assert len(acc_bp) == len(windows) - 1
    assert np.all(acc_bp >= 0)
    assert np.all(acc_bp[:-1] <= window_size + 1e-6)
    good = [[windows[i], windows[i + 1]] for i in range(len(windows) - 1) if acc_bp[i] >= cutoff_bp]
    good_intervals = merge_intervals(good)
    if len(good_intervals) == 0:
        raise ValueError("No windows pass the accessibility cutoff; nothing to keep.")
    acc_win = ratemap_from_keep_intervals(good_intervals, sequence_length)
    mask = and_ratemaps_binary(acc_mu, acc_win)
    return mask, good_intervals, windows, acc_bp
