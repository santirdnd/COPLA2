#!/usr/bin/env python3

import os
import sys
import json
import pickle
import pathlib
import warnings
import argparse
import multiprocessing
import numpy as np
import graph_tool.all as gt

from psani import get_databases_from_json

prog = os.path.basename(sys.argv[0])
version = '2.0.0'
# Format version of the manifest.json of the reference databases (see copla2_setup_db.py)
MANIFEST_FORMAT_VERSION = 1

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Predict a plasmid PTU from precomputed ANI values')
    parser.add_argument('query_db', type=str,
        help='JSON file of the query database created by psani.py makedb. Every query plasmid in it is reported, including those with no ANI value')
    parser.add_argument('ref_db', type=str,
        help='reference database directory prepared by copla2_setup_db.py (e.g. databases/Copla_RS224_v1.0). Its manifest.json gives the reference network')
    parser.add_argument('ani', type=str,
        help='TSV file with the ANI values computed by psani.py dist between the query database and the reference database')
    parser.add_argument('outdir', type=str,
        help='output directory')
    parser.add_argument('--size-ratio', type=float,
        default=2,
        help='maximum ratio of plasmid sizes for an ANI link between a query and a reference plasmid. The reference networks were built with 2. Set to 0 to deactivate this constraint. Default: 2')
    parser.add_argument('--seed', type=int,
        default=41,
        help='random seed for the MCMC refinement of the partition, set for each query so its result does not depend on the other queries or on -n. Default: 41')
    parser.add_argument('-n', '--n-jobs', type=int,
        default=0,
        help='number of queries processed in parallel. Default: 0 (use all available CPUs)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {version}')

    return parser.parse_args()


def get_available_cpus() -> int:
    """
    Return the number of CPUs available to this process.

    On Linux, respect CPU affinity (cgroups, sched_setaffinity, etc.).
    On other systems, fallback to os.cpu_count().
    Always return at least 1.
    """
    # Linux with sched_getaffinity
    if hasattr(os, "sched_getaffinity"):
        try:
            return len(os.sched_getaffinity(0))
        except (OSError, AttributeError):
            pass

    # Portable fallback
    count = os.cpu_count()
    return count if (count is not None and count > 0) else 1


def validate_arguments(args):
    available_cpus = get_available_cpus()

    load_reference_db(args)

    if not os.path.isfile(args.ani):
        raise SystemExit(f"ANI TSV file not found: {args.ani}")

    if not os.path.isfile(args.query_db):
        raise SystemExit(f"Query database file not found: {args.query_db}")

    if args.n_jobs < 0:
        print(f"{prog}: error: number of jobs ({args.n_jobs}) must be >= 0", file=sys.stderr)
        sys.exit(2)

    if args.n_jobs == 0:
        args.n_jobs = available_cpus

    if (args.size_ratio != 0) and (args.size_ratio < 1):
        print(f"{prog}: error: ratio of plasmid sizes ({args.size_ratio}) must be 0 or >= 1", file=sys.stderr)
        sys.exit(2)


def load_ani_by_query(fname):
    ani_by_query = {}

    with open(fname, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')

            # Skip comments and empty lines
            if (line == '') or line.startswith('#'):
                continue

            items = line.split('\t')
            qry = items[0]

            if qry not in ani_by_query:
                ani_by_query[qry] = []
            ani_by_query[qry].append(items)

    return ani_by_query


def load_reference_db(args):
    # The manifest.json of the reference database directory gives the path of the network
    # pickle (relative to the directory) and the database identifier and version
    manifest_path = os.path.join(args.ref_db, 'manifest.json')
    if not os.path.isfile(manifest_path):
        raise SystemExit(f"{args.ref_db} is not a COPLA2 reference database directory (no manifest.json)")
    try:
        with open(manifest_path, 'rt', encoding='utf-8') as fh:
            manifest = json.load(fh)
        format_version = manifest.get('format_version')
        args.ref_graph = os.path.join(args.ref_db, manifest['files']['graph']['path'])
        args.ref_label = f"{manifest['id']} version {manifest['version']}"
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise SystemExit(f"Invalid manifest {manifest_path}: {exc}")
    if format_version != MANIFEST_FORMAT_VERSION:
        found = 'has no format version' if format_version is None else f'has format version {format_version}'
        raise SystemExit(f"{manifest_path} {found}, but {prog} reads version {MANIFEST_FORMAT_VERSION}")

    if not os.path.isfile(args.ref_graph):
        raise SystemExit(f"Reference graph file not found: {args.ref_graph}. Prepare the database with copla2_setup_db.py")


def sanitize_query_id(query_id):
    return os.path.basename(query_id).replace('/', '_')


def insert_ani_edges(graph, ani_rows, v_qry, query_id, size_ratio):
    idx = {}
    for v in graph.vertices():
        idx[graph.vp.AccessionVersion[v]] = v

    query_size = None

    for items in ani_rows:
        src_len = int(items[2])
        trg_len = int(items[3])

        if query_size is None:
            query_size = src_len
        elif query_size != src_len:
            raise SystemExit(
                f"{prog}: error: inconsistent query size for query "
                f"{query_id} ({query_size} vs {src_len})"
            )

        # Only link plasmids of similar size (size_ratio 0 links them all)
        if (size_ratio != 0) and ((src_len > size_ratio * trg_len) or (trg_len > size_ratio * src_len)):
            continue

        ref = items[1]

        if ref not in idx:
            continue

        # The links are not weighted: the ANI value is not stored in the network
        graph.add_edge(v_qry, graph.vertex(idx[ref]))

    return query_size


def block_annotation(graph, state):
    levels = state.get_levels()

    # Find the informative hierarchical levels (i.e. the non-redundant levels, those with non-equivalent block assignment)
    def check_level_redundancy(l):
        x = state.project_partition(l, 0).a
        y = state.project_partition(l+1, 0).a
        return gt.partition_overlap(x, y, norm=True) == 1
    L = len(levels) - 1
    redundant_levels = [False]+list(map(check_level_redundancy, reversed(range(L))))
    nr_levels = [L-i for i, x in enumerate(redundant_levels) if not x]

    b = levels[0].get_blocks()
    bcc = graph.new_vertex_property('string')
    bcc4 = graph.new_vertex_property('string')
    for i in np.unique(b.a):
        b_filter = (b.a == i)
        u = gt.GraphView(graph, vfilt=b_filter)
        tmp = []
        r = u.get_vertices()[0]
        for l in range(len(levels)):
            r = levels[l].get_blocks()[r]
            if l in nr_levels:
                tmp.append(str(r))
        tmp.reverse()
        comp, hist = gt.label_components(u)
        for v in u.vertices():
            tag = '_'.join(tmp + [str(comp[v])])
            bcc[v] = tag
            bcc4[v] = tag if (hist[comp[int(v)]] >= 4) else '-'

    return((b, bcc, bcc4))


def ptu_annotation(graph):
    # ------------------------------------------------------------------
    # Initial cluster setup from BlockCC4 labels
    # ------------------------------------------------------------------
    bcc4_arr = np.array(list(graph.vp.BlockCC4))
    labels = np.unique(bcc4_arr)
    labels = labels[labels != '-']
    k = len(labels)

    # Vertex -> initial cluster id (-1 for unassigned '-' vertices)
    label_to_cid = {label: i for i, label in enumerate(labels)}
    num_vertices_total = graph.num_vertices()
    cluster_of = np.full(num_vertices_total, -1, dtype=np.int64)
    for v_idx, label in enumerate(bcc4_arr):
        cid = label_to_cid.get(label, -1)
        if cid >= 0:
            cluster_of[v_idx] = cid

    # Precompute intra- and inter-cluster edge counts in ONE pass over the graph.
    # This avoids rebuilding a GraphView + num_edges() for every candidate pair,
    # which was the main bottleneck (O(k^2 * |E|) -> O(|E| + k^2)).
    if k > 0:
        sizes = np.bincount(cluster_of[cluster_of >= 0], minlength=k).astype(np.int64)
    else:
        sizes = np.zeros(0, dtype=np.int64)
    intra_edges_arr = np.zeros(k, dtype=np.float64)
    # E[i,j] for i!=j = number of edges between clusters i and j (symmetric)
    E = np.zeros((k, k), dtype=np.float64)

    for e in graph.edges():
        s = int(e.source())
        t = int(e.target())
        if s == t:
            continue  # skip self-loops (ANI=1.0 autoloops)
        cs = cluster_of[s]
        ct = cluster_of[t]
        if cs < 0 or ct < 0:
            continue
        if cs == ct:
            intra_edges_arr[cs] += 1
        else:
            E[cs, ct] += 1
            E[ct, cs] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        pair_count = sizes * (sizes - 1) / 2.0
        intra_density = np.where(pair_count > 0, intra_edges_arr / pair_count, 0.0)

    # Per-cluster metadata (mutated across merges)
    clusters_of_cid = [[labels[i]] for i in range(k)]
    base_hier = [labels[i] for i in range(k)]
    members = {i: np.where(cluster_of == i)[0] for i in range(k)}
    active = np.ones(k, dtype=bool)
    checked = np.zeros(k, dtype=bool)

    # ------------------------------------------------------------------
    # Main merge loop
    #   - No DataFrame rebuilt inside the loop (plain NumPy sorts).
    #   - Vectorized evaluation of all candidates against `current`.
    #   - Merges update the cluster-level matrix E in O(k) instead of O(|E|).
    # ------------------------------------------------------------------
    while True:
        unchecked_mask = active & ~checked
        if not unchecked_mask.any():
            break

        unchecked = np.where(unchecked_mask)[0]
        # Sort unchecked by (intra_density desc, sizes desc) via stable lexsort
        if len(unchecked) > 1:
            order = unchecked[np.lexsort((-sizes[unchecked], -intra_density[unchecked]))]
        else:
            order = unchecked
        current = int(order[0])
        rest = order[1:]

        if len(rest) == 0:
            checked[current] = True
            continue

        # Vectorized evaluation of all candidates against `current`
        inter = E[current, rest]
        cand_sizes = sizes[rest]
        cand_intra_density = intra_density[rest]

        with np.errstate(divide='ignore', invalid='ignore'):
            inter_density = inter / (sizes[current] * cand_sizes)
            denom = intra_density[current] * cand_intra_density
            inter_density_adjusted = np.where(denom > 0, inter_density / denom, 0.0)

        accepted = inter_density_adjusted > 0.4
        if not accepted.any():
            checked[current] = True
            continue

        # Pick candidate maximizing (inter_density desc, size desc) — matches original
        acc_idx = np.where(accepted)[0]
        acc_order = acc_idx[np.lexsort((
            -cand_sizes[acc_idx],
            -inter_density[acc_idx],
        ))]
        cand = int(rest[int(acc_order[0])])

        # ----- Merge `cand` into `current` -----
        merged_inter = E[current, cand]
        new_intra_edges = intra_edges_arr[current] + intra_edges_arr[cand] + merged_inter
        new_size = sizes[current] + sizes[cand]

        intra_edges_arr[current] = new_intra_edges
        sizes[current] = new_size
        pc = new_size * (new_size - 1) / 2.0
        intra_density[current] = new_intra_edges / pc if pc > 0 else 0.0

        # Update inter-edges to all other clusters in O(k):
        #   1. zero the (current, cand) link (it's now intra)
        #   2. accumulate cand's row/column into current's
        #   3. retire cand
        E[current, cand] = 0
        E[cand, current] = 0
        E[current, :] += E[cand, :]
        E[:, current] += E[:, cand]
        E[cand, :] = 0
        E[:, cand] = 0
        active[cand] = False

        # Cluster metadata
        clusters_of_cid[current].extend(clusters_of_cid[cand])
        base_hier[current] = os.path.commonprefix(
            [base_hier[current], base_hier[cand]]
        ).rstrip('_')
        members[current] = np.concatenate([members[current], members[cand]])
        del members[cand]
        # `current` stays unchecked so it can merge again next iteration.

    # ------------------------------------------------------------------
    # Build v_sHSBM vertex property (label = first original cluster of each group)
    # ------------------------------------------------------------------
    v_sHSBM = graph.new_vertex_property('string')
    for cid, vs in members.items():
        label = clusters_of_cid[cid][0]
        for v_idx in vs:
            v_sHSBM[int(v_idx)] = label
    for v in graph.vertices():
        if v_sHSBM[v] == '':
            v_sHSBM[v] = '-'

    # ------------------------------------------------------------------
    # Assign orphan plasmids (v_sHSBM == '-') to a unique neighbor cluster
    # if the connection is strong enough.
    #
    # Optimized: reuses per-cluster size/intra_density from the merge loop
    # above, so no GraphView / num_edges() calls are needed here.
    # ------------------------------------------------------------------
    label_info = {}
    for cid in members:
        lab = clusters_of_cid[cid][0]
        label_info[lab] = (sizes[cid], intra_density[cid])

    for v in graph.vertices():
        if v_sHSBM[v] != '-':
            continue
        v_idx = int(v)

        # Count edges from v to each surviving neighbor cluster in one pass
        cluster_edge_count = {}
        for w in v.out_neighbors():
            w_idx = int(w)
            if w_idx == v_idx:
                continue  # skip self-loop
            lab = v_sHSBM[w]
            if lab == '-':
                continue
            cluster_edge_count[lab] = cluster_edge_count.get(lab, 0) + 1

        # Evaluate each candidate cluster; require exactly one to pass the threshold
        qualifying_count = 0
        qualifying_label = ''
        for lab, e_count in cluster_edge_count.items():
            info = label_info.get(lab)
            if info is None:
                continue
            cluster_size, cluster_intra_density = info
            if cluster_size <= 0 or cluster_intra_density <= 0:
                continue
            inter_density = e_count / cluster_size
            inter_density_adjusted = inter_density / cluster_intra_density
            if inter_density_adjusted > 0.4:
                qualifying_count += 1
                qualifying_label = lab

        if qualifying_count == 1:
            v_sHSBM[v] = qualifying_label

    # Rename as many PTUs as possible aligning the partition labels
    v_PtuRef = graph.vp.PtuRef.copy()
    ptuRef_index={}
    for i, p in enumerate(np.unique(list(v_PtuRef))):
        ptuRef_index[p] = i
    sHSBM_index={}
    for i, p in enumerate(np.unique(list(v_sHSBM))):
        sHSBM_index[p] = i
    ptuRef_num = []
    sHSBM_num = []
    for v in graph.vertices():
        ptuRef_num.append(ptuRef_index[v_PtuRef[v]])
        sHSBM_num.append(sHSBM_index[v_sHSBM[v]])
    sHSBM_aligned = gt.align_partition_labels(sHSBM_num, ptuRef_num)
    v_Ptu = graph.new_vertex_property('string')
    goodPtuLabels = {}
    for v in graph.vertices():
        if v_sHSBM[v] == '-':
            v_Ptu[v] = '-'
        elif sHSBM_aligned[int(v)] == ptuRef_num[int(v)]:
            v_Ptu[v] = v_PtuRef[v]
            if v_sHSBM[v] not in goodPtuLabels:
                goodPtuLabels[v_sHSBM[v]] = v_PtuRef[v]
        else:
            v_Ptu[v] = '?' + v_sHSBM[v]
    for v in graph.vertices():
        if v_Ptu[v].startswith('?'):
            if v_sHSBM[v] in goodPtuLabels:
                v_Ptu[int(v)] = goodPtuLabels[v_sHSBM[v]]

    return(v_sHSBM, v_Ptu)


def write_classes(graph, filename):
    with open(filename, 'wt') as fh:
        header = "AccessionVersion\tPTU_Ref\tCComp\tBlock\tBlockCC\tBlockCC4"
        if 'sHSBM' in graph.vp.keys():
            header += "\tsHSBM\tPTU"

        fh.write(header + "\n")

        for v in graph.vertices():
            AccessionVersion = graph.vp.AccessionVersion[v]
            PtuRef = str(graph.vp.PtuRef[v])
            CComp = str(graph.vp.CComp[v])
            Block = str(graph.vp.Block[v])
            BlockCC = graph.vp.BlockCC[v]
            BlockCC4 = graph.vp.BlockCC4[v]

            if 'sHSBM' in graph.vp.keys():
                sHSBM = graph.vp.sHSBM[v]
                Ptu = graph.vp.Ptu[v]
                fh.write("\t".join(
                    (AccessionVersion, PtuRef, CComp, Block, BlockCC, BlockCC4, sHSBM, Ptu)
                ) + "\n")
            else:
                fh.write("\t".join(
                    (AccessionVersion, PtuRef, CComp, Block, BlockCC, BlockCC4)
                ) + "\n")


def get_related_plasmids(graph, vfilt, qry_acc):
    ref = []
    new = []
    ref_index = {}
    new_index = {}

    u = gt.GraphView(graph, vfilt=vfilt)
    for v in u.vertices():
        if u.vp.AccessionVersion[v] == qry_acc:
            # As the query is not present in the reference network, we get rid of it to calculate the fit of both partitions
            continue
        if u.vp.sHSBMRef2[v] not in ref_index:
            if len(list(ref_index.values())) > 0:
                ref_index[u.vp.sHSBMRef2[v]] = np.max(list(ref_index.values())) + 1
            else:
                ref_index[u.vp.sHSBMRef2[v]] = 0
        ref.append(ref_index[u.vp.sHSBMRef2[v]])
        if u.vp.sHSBM2[v] not in new_index:
            if len(list(new_index.values())) > 0:
                new_index[u.vp.sHSBM2[v]] = np.max(list(new_index.values())) + 1
            else:
                new_index[u.vp.sHSBM2[v]] = 0
        new.append(new_index[u.vp.sHSBM2[v]])

    return(ref, new)

def run_single_query(task):
    query_id, ani_rows, args = task

    summary_lines = []

    def log(message):
        summary_lines.append(message)

    # One process per query; keep graph-tool single-threaded inside each worker, which
    # is also needed for the MCMC to be reproducible with a fixed seed
    gt.openmp_set_num_threads(1)
    gt.seed_rng(args.seed)
    np.random.seed(args.seed)

    output_prefix = os.path.join(args.outdir, sanitize_query_id(query_id))

    # Load the reference plasmid network. The block model state of the pickle keeps a
    # 'sampling' argument that graph-tool no longer uses and warns about (ignoring it) when
    # the state is rebuilt: silence only that warning while loading, and drop the argument
    # so that later copies of the state do not repeat it
    with open(args.ref_graph, 'rb') as fh, warnings.catch_warnings():
        warnings.filterwarnings('ignore', message=r"unrecognized keyword arguments: \['sampling'\]",
                                category=UserWarning, module=r'graph_tool\.inference')
        [g, state] = pickle.load(fh)
    state.state_args.pop('sampling', None)

    # Insert query into the network
    qry_acc = query_id
    v_qry = g.add_vertex()
    g.vp.AccessionVersion[v_qry] = qry_acc

    # Add query autoloop edge
    g.add_edge(v_qry, v_qry)

    # Populate query edges with database plasmids and recover query size from ANI TSV
    query_size = insert_ani_edges(g, ani_rows, v_qry, query_id, args.size_ratio)

    if 'Size' in g.vp:
        g.vp.Size[v_qry] = query_size if query_size is not None else 0

    has_external_neighbors = any(w != v_qry for w in v_qry.out_neighbors())

    # If query has no external neighbors, it will be isolated (only self-loop)
    if not has_external_neighbors:
        log("PTU could not be assigned")
        if query_id in args.no_fragments:
            log(f"Query has no fragments of {args.window} bp, so it could not be compared with the reference plasmids.")
            notes = f"No fragments of {args.window} bp to compare with reference plasmids; query is a singleton"
        else:
            log("Query has no ANI edges to reference plasmids.")
            notes = "No ANI links to reference plasmids; query is a singleton"
        with open(output_prefix + '.ptu_prediction.tsv', 'wt') as fh:
            fh.write("#QueryID\tPredicted_PTU\tScore\tNotes\n")
            fh.write(f"{query_id}\t-\t1.0000\t{notes}\n")
        with open(output_prefix + '.related_plasmids.tsv', 'wt') as fh:
            fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
            fh.write(f"{qry_acc}\t-\t-\n")
        with open(output_prefix + '.summary.txt', 'wt') as fh:
            fh.write("\n".join(summary_lines + [args.run_info]) + "\n")
        return

    # Use the most frequent block among the neighbors for query block initialitation
    neigh_blocks = list()
    for w in v_qry.out_neighbors():
        if g.vp.AccessionVersion[w] == qry_acc:
            continue
        neigh_blocks.append(g.vp.Block[w])

    def most_frequent_element(lst):
        unique, counts = np.unique(lst, return_counts=True)
        index = np.argmax(counts)
        return unique[index]

    # Save as our reference the existing clustering and PTU assignment from the pickle file
    qry_null = '<NULL>' # Using a placeholder until PTU is assigned
    g.vp.Block[v_qry] = most_frequent_element(neigh_blocks) # Block is a numeric value
    g.vertex_properties['BlockRef'] = g.vp.Block.copy()
    g.vp.BlockCC[v_qry] = qry_null
    g.vertex_properties['BlockCCRef'] = g.vp.BlockCC.copy()
    g.vp.BlockCC4[v_qry] = qry_null
    g.vertex_properties['BlockCC4Ref'] = g.vp.BlockCC4.copy()
    g.vp.sHSBM[v_qry] = qry_null
    g.vertex_properties['sHSBMRef'] = g.vp.sHSBM.copy()
    g.vp.Ptu[v_qry] = qry_null
    g.vertex_properties['PtuRef'] = g.vp.Ptu.copy()

    # Check if the query belongs to a graph component of 4 or less members, in which case we can quit
    comp, hist = gt.label_components(g)
    g.vertex_properties['CComp'] = comp
    if hist[comp[v_qry]] <= 4:
        i_filter = (comp.a == comp[v_qry])
        u = gt.GraphView(g, vfilt=i_filter)
        cl_size = u.num_vertices()
        if cl_size < 4:
            ptu_pred = '-'
            log('PTU could not be assigned')
            log(f"Query is part of a graph component of size {cl_size}")
            log('However, at least four members are required for PTU assignation')
            log('This plasmid could form part of a new, still unnamed, PTU')
        else:
            ptu_pred = 'PTU-?'
            log('New (putative) PTU')
            log(f"Query is part of a graph component of size {cl_size}")
            log('This plasmid could form part of a new, still unnamed, PTU')
        with open(output_prefix+'.ptu_prediction.tsv', 'wt') as fh:
            fh.write("#QueryID\tPredicted_PTU\tScore\tNotes\n")
            score = 1.0
            if cl_size < 4:
                notes = f"PTU could not be assigned. Query is part of a graph component of size {cl_size}"
            else:
                notes = f"New (putative) PTU. Query is part of a graph component of size {cl_size}"
            fh.write(f"{query_id}\t{ptu_pred}\t{score:.4f}\t{notes}\n")
        with open(output_prefix+'.related_plasmids.tsv', 'wt') as fh:
            fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
            for v in u.vertices():
                note = '*'
                fh.write(f"{u.vp.AccessionVersion[v]}\t{u.vp.PtuRef[v]}\t{note}\n")
        with open(output_prefix + '.summary.txt', 'wt') as fh:
            fh.write("\n".join(summary_lines + [args.run_info]) + "\n")
        return

    # The SBM state keeps its own vertex and edge weights, and the query vertex and edges
    # get 0: the model would ignore the query edges and the MCMC would never move the
    # query. Give them the weight of the reference ones (1), including the self-loop
    state.get_levels()[0].vweight[v_qry] = 1
    eweight = state.get_levels()[0].eweight
    for e in v_qry.out_edges():
        eweight[e] = 1

    # Apply SBM to cluster the query. The query starts in the most frequent block of its
    # neighbors: the state does not read g.vp.Block, and would leave it in block 0
    bs = state.get_bs()
    bs[0][int(v_qry)] = g.vp.Block[v_qry]
    new_state = state.copy(bs=bs)
    new_state.multiflip_mcmc_sweep(d=0, psplit=0, pmerge=0, pmergesplit=0, niter=100)
    new_state.multiflip_mcmc_sweep(d=0, psplit=0, pmerge=0, pmergesplit=0, niter=100)
    new_state.multiflip_mcmc_sweep(d=0, psplit=0, pmerge=0, pmergesplit=0, niter=100)
    new_state.multiflip_mcmc_sweep(d=0, psplit=0, pmerge=0, pmergesplit=0, niter=100)
    new_state.multiflip_mcmc_sweep(d=0, psplit=0, pmerge=0, pmergesplit=0, niter=100)

    # Apply sHSBM algorithm to the new graph with the query
    (v_Block, v_BlockCC, v_BlockCC4) = block_annotation(g, new_state)
    g.vertex_properties['Block'] = v_Block
    g.vertex_properties['BlockCC'] = v_BlockCC
    g.vertex_properties['BlockCC4'] = v_BlockCC4
    (v_sHSBM, v_Ptu) = ptu_annotation(g)
    g.vertex_properties['sHSBM'] = v_sHSBM
    g.vertex_properties['Ptu'] = v_Ptu
    write_classes(g, output_prefix+'.sHSBM.tsv')

    # TODO: Code is easier if sHSBM column has the data of all clusters and not only those with 4 or more members
    v_sHSBM2 = g.vp.sHSBM.copy()
    v_sHSBMRef2 = g.vp.sHSBMRef.copy()
    for v in g.vertices():
        if v_sHSBM2[v] == '-':
            v_sHSBM2[v] = g.vp.BlockCC[v]
        if v_sHSBMRef2[v] == '-':
            v_sHSBMRef2[v] = g.vp.BlockCCRef[v]
    g.vertex_properties['sHSBM2'] = v_sHSBM2
    g.vertex_properties['sHSBMRef2'] = v_sHSBMRef2

    # Get all plasmids of the PTU assigned to the query
    strict_filter = (np.array(list(g.vp.sHSBM2)) == g.vp.sHSBM2[v_qry])

    # Use the reference clustering of these plasmids to expand the selection to all plasmids of all PTUs with a member in the new PTU assigned to the query
    sHSBMRef2_list = set()
    u = gt.GraphView(g, vfilt=strict_filter)
    for v in u.vertices():
        sHSBMRef2_list.add(u.vp.sHSBMRef2[v])
    expand_filter = np.full_like(strict_filter, False)
    for i in sHSBMRef2_list:
        i_filter = (np.array(list(g.vp.sHSBMRef2)) == i)
        expand_filter = np.logical_or(expand_filter, i_filter)
    (ref_expand, new_expand) = get_related_plasmids(g, expand_filter, qry_acc)
    if len(ref_expand) == 0:
        # Query was assigned to a singleton
        overlap_expand = 1.0
    else:
        overlap_expand = gt.partition_overlap(ref_expand, new_expand, norm=True)

    # Predicted PTU is the most frequent reference label among the new cluster
    w = gt.GraphView(g, vfilt=expand_filter)
    cl_size = u.num_vertices()
    tmp_ptu = np.delete(np.array(list(u.vp.PtuRef)), np.argwhere(np.array(list(u.vp.PtuRef)) == qry_null))
    if cl_size < 4:
        ptu_pred = '-'
        ptu_related = '-'
        # Make sure ptu_related is not '-' if there is another option
        # This works because cl_size <= 3 and one of them is the previously deleted qry_null
        for i in tmp_ptu:
            if i != '-':
                ptu_related = i
        log('PTU could not be assigned')
        log('Query is part of a sHSBM cluster of size {}'.format(cl_size))
        log('However, at least four members are required for PTU assignation')
        log('This plasmid could form part of a new, still unnamed, PTU')
        if ptu_related != '-':
            log('Query is related to {} plasmids'.format(ptu_related))
        with open(output_prefix+'.ptu_prediction.tsv', 'w') as fh:
            fh.write("#QueryID\tPredicted_PTU\tScore\tNotes\n")
            score = overlap_expand
            notes = 'PTU could not be assigned. Query is part of a sHSBM cluster of size {}'.format(cl_size)
            fh.write("{}\t{}\t{:.4f}\t{}\n".format(query_id, ptu_pred, score, notes))
        with open(output_prefix+'.related_plasmids.tsv', 'w') as fh:
            fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
            for v in w.vertices():
                note = '*' if g.vp.sHSBM[v] == g.vp.sHSBM[v_qry] else ''
                fh.write("{}\t{}\t{}\n".format(w.vp.AccessionVersion[v], w.vp.PtuRef[v], note))
    else:
        unique, counts = np.unique(list(tmp_ptu), return_counts=True)
        ptu_pred = unique[counts.argmax()]
        ptu_related = '-'
        if (ptu_pred == '-') or (ptu_pred != u.vp.Ptu[v_qry]):
            if ptu_pred != u.vp.Ptu[v_qry]:
                ptu_related = ptu_pred
            ptu_pred = 'PTU-?'
            log('New (putative) PTU')
            log('Query is part of a sHSBM cluster of size {}'.format(cl_size))
            log('This plasmid could form part of a new, still unnamed, PTU')
            if ptu_related != '-':
                log('Query is related to {} plasmids'.format(ptu_related))
            with open(output_prefix+'.ptu_prediction.tsv', 'w') as fh:
                fh.write("#QueryID\tPredicted_PTU\tScore\tNotes\n")
                score = overlap_expand
                notes = 'New (putative) PTU. Query is part of a sHSBM cluster of size {}'.format(cl_size)
                fh.write("{}\t{}\t{:.4f}\t{}\n".format(query_id, ptu_pred, score, notes))
            with open(output_prefix+'.related_plasmids.tsv', 'w') as fh:
                fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
                for v in w.vertices():
                    note = '*' if g.vp.sHSBM[v] == g.vp.sHSBM[v_qry] else ''
                    fh.write("{}\t{}\t{}\n".format(w.vp.AccessionVersion[v], w.vp.PtuRef[v], note))
        else:
            log('Query is a {} plasmid'.format(ptu_pred))
            log('Query is part of a sHSBM cluster of size {}'.format(cl_size))
            with open(output_prefix+'.ptu_prediction.tsv', 'w') as fh:
                fh.write("#QueryID\tPredicted_PTU\tScore\tNotes\n")
                score = overlap_expand
                notes = 'Query is a {} plasmid'.format(ptu_pred)
                fh.write("{}\t{}\t{:.4f}\t{}\n".format(query_id, ptu_pred, score, notes))
            with open(output_prefix+'.related_plasmids.tsv', 'w') as fh:
                fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
                for v in w.vertices():
                    note = '*' if g.vp.sHSBM[v] == g.vp.sHSBM[v_qry] else ''
                    fh.write("{}\t{}\t{}\n".format(w.vp.AccessionVersion[v], w.vp.PtuRef[v], note))

    with open(output_prefix + '.summary.txt', 'wt') as fh:
        fh.write("\n".join(summary_lines + [args.run_info]) + "\n")


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    args.run_info = f"COPLA2 {version}; reference database {args.ref_label}"

    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)

    ani_by_query = load_ani_by_query(args.ani)
    query_db = get_databases_from_json(args.query_db)
    query_ids = list(query_db['seqs'])
    args.window = query_db['window']
    args.no_fragments = {query_id for query_id, desc in query_db['seqs'].items() if desc['fragments'] == 0}

    tasks = [(query_id, ani_by_query.get(query_id, []), args) for query_id in query_ids]

    if len(query_ids) == 0:
        raise SystemExit(f"{prog}: error: query database contains no plasmids")

    n_jobs = min(args.n_jobs, len(tasks))

    with multiprocessing.Pool(n_jobs) as pool:
        pool.map(run_single_query, tasks)

    # Combine per-query PTU prediction files into a single file
    combined_file = os.path.join(args.outdir, "ptu_predictions.tsv")

    first = True
    with open(combined_file, 'wt') as out_fh:
        for query_id, _, _ in tasks:
            prefix = os.path.join(args.outdir, sanitize_query_id(query_id))
            fname = prefix + ".ptu_prediction.tsv"

            if not os.path.exists(fname):
                continue

            with open(fname, 'rt') as fh:
                for line in fh:
                    # Write header only once
                    if line.startswith("#"):
                        if first:
                            out_fh.write(line)
                            first = False
                        continue
                    out_fh.write(line)

    if first:
        raise SystemExit("No PTU prediction files were generated")


if __name__ == '__main__':
    main()
