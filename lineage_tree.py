#!/usr/bin/env python3

import argparse
import collections
import glob
import itertools
import os
import re
import shutil
import subprocess

import numpy as np

MIN_NODE_READS = 2       # reads a subclone needs to exist
MIN_MUT_READS = 2        # reads a mutation needs before it can define one
MIN_STEP_MUTATIONS = 1   # new mutations a node must add over its ancestor
BCR_V = ("IGHV", "IGKV", "IGLV")   # only BCR clonotypes are analysed, not TCR
CORE_FRAC = 0.9          # share of reads that must cover a position to compare on it



MATCH, MISMATCH, GAP_OPEN, GAP_EXT = 2, -4, -5, -1
NEG = -10 ** 8


def align(ref, read):
    """Affine-gap alignment with free end gaps. Returns span, per-ref-column bases,
    insertions keyed by the ref position they precede."""
    n, m = len(ref), len(read)
    ra = np.frombuffer(ref.encode(), dtype=np.uint8)
    qa = np.frombuffer(read.encode(), dtype=np.uint8)

    M = np.full((n + 1, m + 1), NEG, dtype=np.int32)
    E = np.full((n + 1, m + 1), NEG, dtype=np.int32)
    F = np.full((n + 1, m + 1), NEG, dtype=np.int32)
    M[0, 0] = 0
    E[0, 1:] = 0          # leading read bases consumed for free
    F[1:, 0] = 0          # leading ref bases consumed for free

    jj = np.arange(m + 1, dtype=np.int32)
    for i in range(1, n + 1):
        prev = np.maximum(np.maximum(M[i - 1], E[i - 1]), F[i - 1])
        M[i, 1:] = prev[:-1] + np.where(qa == ra[i - 1], MATCH, MISMATCH)
        F[i, :] = np.maximum(M[i - 1] + GAP_OPEN + GAP_EXT, F[i - 1] + GAP_EXT)
        # E depends on M within this row -> running-max formulation
        b = M[i] + GAP_OPEN - GAP_EXT * jj
        run = np.maximum.accumulate(b)
        E[i, 1:] = GAP_EXT * jj[1:] + run[:-1]
        E[i, 0] = NEG

    # free trailing gaps: best cell on last row or last column
    best, bi, bj, bs = NEG, n, m, "M"
    for mat, name in ((M, "M"), (E, "E"), (F, "F")):
        j = int(np.argmax(mat[n]))
        if mat[n, j] > best:
            best, bi, bj, bs = int(mat[n, j]), n, j, name
        i = int(np.argmax(mat[:, m]))
        if mat[i, m] > best:
            best, bi, bj, bs = int(mat[i, m]), i, m, name

    cols, ins = {}, {}
    i, j, state = bi, bj, bs
    pend = []
    while i > 0 and j > 0:
        if state == "M":
            cols[i - 1] = read[j - 1]
            if pend:
                ins[i] = "".join(reversed(pend))
                pend = []
            score = M[i, j] - (MATCH if read[j - 1] == ref[i - 1] else MISMATCH)
            i, j = i - 1, j - 1
            for cand in ("M", "E", "F"):
                if {"M": M, "E": E, "F": F}[cand][i, j] == score:
                    state = cand
                    break
        elif state == "F":
            cols[i - 1] = ""
            if M[i - 1, j] + GAP_OPEN + GAP_EXT == F[i, j]:
                state = "M"
            i -= 1
        else:
            pend.append(read[j - 1])
            if M[i, j - 1] + GAP_OPEN + GAP_EXT == E[i, j]:
                state = "M"
            j -= 1
    if pend and i >= 0:
        ins[i] = "".join(reversed(pend))

    covered = sorted(cols)
    if not covered:
        return None
    return {"span": (covered[0], covered[-1] + 1), "cols": cols, "ins": ins}


def signature(ref, al):
    """Variant signature relative to ref, restricted to the aligned span."""
    s, e = al["span"]
    sig = []
    run = []
    for i in range(s, e):
        b = al["cols"].get(i)
        if b == "":
            run.append(i)
            continue
        if run:
            sig.append(("del", run[0], len(run)))
            run = []
        if b is not None and b != ref[i]:
            sig.append(("sub", i, ref[i], b))
    if run:
        sig.append(("del", run[0], len(run)))
    for i, seq in al["ins"].items():
        if s <= i <= e:
            sig.append(("ins", i, seq))
    return frozenset(sig)




_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[i] for i, (a, b, c) in enumerate(itertools.product(_B, repeat=3))}

# IMGT unique numbering: codon ranges are fixed across all V genes.
REGIONS = [("FR1", 1, 26), ("CDR1", 27, 38), ("FR2", 39, 55),
           ("CDR2", 56, 65), ("FR3", 66, 104)]
REGION_NAMES = [r[0] for r in REGIONS]


def region_of(codon):
    for name, lo, hi in REGIONS:
        if lo <= codon <= hi:
            return name
    return "CDR3"


def load_germline(path, patterns=("IG*V*.fasta*", "IG*V*.fa", "IG*V*.fasta")):
    """Load a germline V reference. Accepts a FASTA file or a directory of them.

    Returns (alleles, is_gapped). Allele keys use the pipeline's naming
    (IGHV1-69D_01), so IMGT '*' separators are converted.
    """
    if os.path.isdir(path):
        files = []
        for pat in patterns:
            files += glob.glob(os.path.join(path, pat))
        if not files:
            files = [p for p in glob.glob(os.path.join(path, "*.fa*"))]
        if not files:
            raise SystemExit(f"no FASTA files found in {path}")
    else:
        if not os.path.exists(path):
            raise SystemExit(f"germline reference not found: {path}")
        files = [path]

    ref, name = {}, None
    for fp in sorted(files):
        for line in open(fp, errors="replace"):
            line = line.rstrip("\n")
            if line.startswith(">"):
                h = line[1:]
                # IMGT: accession|allele|species|...   plain: allele [desc]
                f = h.split("|")
                name = (f[1] if len(f) > 1 and f[1] else h.split()[0]).replace("*", "_")
                ref.setdefault(name, "")
            elif name is not None:
                ref[name] += line.strip()
    is_gapped = any("." in s for s in ref.values())
    return ref, is_gapped



_PARALOG_D = re.compile(r"(\d)D(?=_\d+$)")


def germline_for(ref, allele):
 
    if allele in ref:
        return allele, ref[allele]
    key = _PARALOG_D.sub(r"\1", allele)
    for name in ref:
        if _PARALOG_D.sub(r"\1", name) == key:
            return name, ref[name]
    return None, None


def imgt_codon_map(gapped):
 
    return [gi // 3 + 1 for gi, ch in enumerate(gapped) if ch != "."]


def classify_sub(germ_ungapped, pos, newbase):
 
    ci = (pos // 3) * 3
    codon = germ_ungapped[ci:ci + 3].upper()
    if len(codon) < 3:
        return None
    mut = list(codon)
    mut[pos - ci] = newbase
    a, b = CODON.get(codon), CODON.get("".join(mut))
    if a is None or b is None:
        return None
    return "S" if a == b else "R"


def load_rows(path):
    rows = []
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == len(h):
                rows.append(dict(zip(h, p)))
    return rows


def pick_v_column(rows):
 
    if rows and "v_seq_corrected" in rows[0]:
        return "v_seq_corrected"
    return "v_seq"


# ------------------------------------------------------------------ subclone calls

def describe_relation(sig, root):
    """How a subclone's mutation set relates to sub1's (from subclone_v1.py)."""
    if not sig:
        return "germline_unmutated"
    if sig == root:
        return "root"
    if root and root < sig:
        return f"derived_from_sub1(+{len(sig - root)})"
    if not root:
        return f"derived_from_sub1(+{len(sig)})"
    shared = len(sig & root)
    if shared == 0:
        return "parallel_disjoint"
    return f"partial_overlap(shared={shared})"


def call_subclones(reads, germ, vcol, min_node_reads=MIN_NODE_READS,
                   min_mut_reads=MIN_MUT_READS, min_step=MIN_STEP_MUTATIONS,
                   core_frac=CORE_FRAC):
  
    cache, per = {}, []
    for r in reads:
        s = r.get(vcol, ".")
        if s in (".", ""):
            per.append((r, None, None)); continue
        if s not in cache:
            al = align(germ, s)
            cache[s] = ((al["span"], frozenset(v for v in signature(germ, al)
                                               if v[0] == "sub")) if al else None)
        got = cache[s]
        per.append((r, got[0], got[1]) if got else (r, None, None))

    spans = [sp for _, sp, _ in per if sp]
    if not spans:
        return {}, None
    cov = collections.Counter()
    for s, e in spans:
        for i in range(s, e):
            cov[i] += 1
    keep = [i for i, c in cov.items() if c >= core_frac * len(spans)]
    if not keep:
        return {}, None
    win = (min(keep), max(keep) + 1)

    sup = collections.Counter()
    trimmed = {}
    for r, sp, sig in per:
        if sp is None or sp[0] > win[0] or sp[1] < win[1]:
            continue
        s2 = frozenset(v for v in sig if win[0] <= v[1] < win[1])
        trimmed[id(r)] = s2
        for v in s2:
            sup[v] += 1
    conf = {v for v, c in sup.items() if c >= min_mut_reads}

    groups = collections.defaultdict(list)
    outside = []
    for r, sp, sig in per:
        if id(r) not in trimmed:
            outside.append(r); continue
        groups[frozenset(trimmed[id(r)] & conf)].append(r)

    accepted = []
    small = []
    for sig, rs in sorted(groups.items(), key=lambda kv: len(kv[0])):
        if len(rs) < min_node_reads:
            small.extend(rs); continue
        parent = None
        for a in accepted:
            if a[0] < sig and (parent is None or len(a[0]) > len(parent[0])):
                parent = a
        if not accepted or parent is None or len(sig - parent[0]) >= min_step:
            accepted.append([sig, list(rs)])
        else:
            parent[1].extend(rs)
    if not accepted:
        big = max(groups.items(), key=lambda kv: len(kv[1]))
        accepted = [[big[0], list(big[1])]]
        small = [r for s, rs in groups.items() if s != big[0] for r in rs]

    for r in small + outside:
        sig = trimmed.get(id(r), frozenset())
        best = min(accepted, key=lambda a: len(sig ^ a[0]))
        best[1].append(r)

    accepted.sort(key=lambda a: (len(a[0]), -len(a[1])))
    return {i + 1: (a[0], a[1]) for i, a in enumerate(accepted)}, win


def fmt_mutations(sig):
    """('sub', pos, from, to) -> 'A12G', sorted by position, 1-based."""
    return ",".join(f"{g}{p + 1}{n}" for _, p, g, n in sorted(sig, key=lambda v: v[1]))


def annotate(sig, germ, cmap):
 
    reg = collections.Counter()
    rs = collections.Counter()
    for _, pos, _, new in sig:
        reg[region_of(cmap[pos]) if pos < len(cmap) else "CDR3"] += 1
        k = classify_sub(germ, pos, new)
        if k:
            rs[k] += 1
    return reg, rs


def build_edges(subclones, root_name="germline"):
    """Parent of each node is the largest proper subset of its mutation set."""
    items = [(lab, sig) for lab, (sig, _) in subclones.items()]
    edges = []
    for lab, sig in sorted(items, key=lambda kv: len(kv[1])):
 
        best_lab, best_sig = root_name, None
        for lab2, sig2 in items:
            if lab2 == lab or not sig2 < sig:
                continue
            if best_sig is None or len(sig2) > len(best_sig):
                best_lab, best_sig = lab2, sig2
        edges.append((best_lab, lab, sig - (best_sig or frozenset())))
    return edges


def to_newick(edges, root="germline"):
    children = collections.defaultdict(list)
    blen = {}
    for parent, child, gained in edges:
        children[parent].append(child)
        blen[child] = len(gained)

    def render(node):
        kids = children.get(node)
        if not kids:
            return node
        inner = ",".join(f"{render(k)}:{blen[k]}" for k in kids)
        return f"({inner}){node}"

    return render(root) + ";"


def apply_mutations(germ, sig):
    """Reconstruct a node's V sequence by writing its substitutions onto the germline,
    so every record shares germline coordinates and comes out already aligned."""
    seq = list(germ)
    for _, pos, _, new in sig:
        if 0 <= pos < len(seq):
            seq[pos] = new
    return "".join(seq)


# ------------------------------------------------------------------------ plotting

# plot_lineage.R ships next to this module; it draws from the nodes/edges tables
PLOT_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_lineage.R")


def plot_lineage(trees, prefix, out_pdf, top=20):

    if shutil.which("Rscript") is None:
        return None, "Rscript is not on PATH"
    if not os.path.exists(PLOT_SCRIPT):
        return None, f"{PLOT_SCRIPT} not found"
    pick = [t["clone"] for t in sorted(trees, key=lambda t: -t["reads"])[:top]]
    if not pick:
        return 0, None
    res = subprocess.run(["Rscript", PLOT_SCRIPT, prefix, "--clonotypes", ",".join(pick),
                          "--top", str(len(pick)), "--out", out_pdf],
                         capture_output=True, text=True)
    if res.returncode != 0:
        err = [l for l in res.stderr.strip().splitlines() if l.strip()]
        return None, err[-1] if err else f"Rscript exited with {res.returncode}"
    done = re.search(r"done: (\d+) pages", res.stderr)
    return (int(done.group(1)) if done else len(pick)), None


def run_subclone_tree(corrected_tsv, germline, out_prefix=None,
                      min_node_reads=MIN_NODE_READS, min_mut_reads=MIN_MUT_READS,
                      min_step=MIN_STEP_MUTATIONS, core_frac=CORE_FRAC,
                      multi_only=False, plot_top=20, verbose=True):
    """Call subclones and emit lineage-tree inputs, plus a PDF of the `plot_top`
    largest multi-subclone trees (0 skips the plot). Returns a summary dict."""
    prefix = out_prefix or corrected_tsv.rsplit(".", 1)[0]
    out_tsv = f"{prefix}.subclones_tree.tsv"
    out_nodes = f"{prefix}.lineage_nodes.tsv"
    out_edges = f"{prefix}.lineage_edges.tsv"
    out_newick = f"{prefix}.lineage.newick"
    out_fasta = f"{prefix}.lineage_sequences.fasta"
    out_pdf = f"{prefix}.lineage_trees.pdf"
    say = print if verbose else (lambda *a, **k: None)

    ref, is_gapped = load_germline(germline)
    rows = load_rows(corrected_tsv)
    vcol = pick_v_column(rows)
    header = list(rows[0].keys()) if rows else []
    # error_correct_new names the columns clonotype/clonotype_size; older tables clone
    ckey, skey = (("clonotype", "clonotype_size") if "clonotype" in header
                  else ("clone", "clone_size"))
    # BCR only (heavy and light chain): TCR reads are dropped from every output
    n_in = len(rows)
    rows = [r for r in rows if r.get("v_gene", "").startswith(BCR_V)]
    n_other = n_in - len(rows)

    byc = collections.defaultdict(list)
    for r in rows:
        for c in ("subclone", "subclone_size", "subclone_mut",
                  "subclone_relation", "subclone_mutations"):
            r[c] = "." if c != "subclone_size" else "0"
        if r.get(ckey, ".") != ".":
            byc[r[ckey]].append(r)

    nodes, edges_all, newicks, fasta, trees = [], [], [], [], []
    skipped = collections.Counter()
    via = collections.Counter()   # alleles found only under their D-paralog name
    n_multi = 0
    for clone in sorted(byc, key=lambda c: -len(byc[c])):
        reads = byc[clone]
        allele = collections.Counter(r["v_gene"] for r in reads).most_common(1)[0][0]
        germ_name, g = germline_for(ref, allele)
        if g is None:
            skipped["no germline for allele"] += 1
            continue
        if germ_name != allele:
            via[f"{allele} -> {germ_name}"] += 1
        germ = g.replace(".", "").upper()
        cmap = imgt_codon_map(g) if is_gapped else []
        subs, win = call_subclones(reads, germ, vcol, min_node_reads, min_mut_reads,
                                   min_step, core_frac)
        if not subs:
            skipped["no alignable V sequence"] += 1
            continue
        # reads identical to the germline are not a subclonotype: they are counted on
        # the germline node, and only mutated groups are numbered sub1, sub2, ...
        unmut = []
        if not subs[1][0]:
            unmut = subs.pop(1)[1]
            subs = {i: v for i, v in enumerate(subs.values(), 1)}
        for r in unmut:
            r["subclone"] = f"{clone}_germline"
            r["subclone_size"] = str(len(unmut))
            r["subclone_mut"] = "0"
            r["subclone_relation"] = "germline_unmutated"
            r["subclone_mutations"] = "none"

        labelled = {f"{clone}_sub{i}": (sig, rs) for i, (sig, rs) in subs.items()}
        root_sig = subs[1][0] if subs else frozenset()
        is_multi = len(subs) > 1
        n_multi += is_multi
        for i, (sig, rs) in subs.items():
            label = f"{clone}_sub{i}"
            rel = "root" if i == 1 else describe_relation(sig, root_sig)
            for r in rs:
                r["subclone"] = label
                r["subclone_size"] = str(len(rs))
                r["subclone_mut"] = str(len(sig))
                r["subclone_relation"] = rel
                r["subclone_mutations"] = fmt_mutations(sig) if sig else "none"

        if multi_only and not is_multi:
            continue

        edges = build_edges({k: (v[0], len(v[1])) for k, v in labelled.items()},
                            root_name=f"{clone}_germline")
        parent_of = {c: p for p, c, _ in edges}

        # the germline root gets its own row so a plot has something to draw from
        # its n_reads are the reads identical to germline (0 when none were seen)
        nodes.append([f"{clone}_germline", clone, allele, ".", "yes", len(unmut), 0,
                      "germline",
                      "."] + [0] * len(REGION_NAMES) + [0, 0])
        for label, (sig, rs) in labelled.items():
            reg, rs_counts = annotate(sig, germ, cmap)
            iso = collections.Counter(r["c_gene"].split("_")[0] for r in rs
                                      if r.get("c_gene", ".") not in (".", ""))
            rel = rs[0]["subclone_relation"] if rs else "."
            nodes.append([label, clone, allele, parent_of.get(label, "."), "no",
                          len(rs), len(sig), rel,
                          iso.most_common(1)[0][0] if iso else "."]
                         + [reg[x] for x in REGION_NAMES]
                         + [rs_counts["R"], rs_counts["S"]])

        for parent, child, gained in edges:
            reg, rs_counts = annotate(gained, germ, cmap)
            cdr = sum(reg[x] for x in ("CDR1", "CDR2"))
            fr = sum(reg[x] for x in ("FR1", "FR2", "FR3"))
            edges_all.append([clone, allele, parent, child, len(gained),
                              fmt_mutations(gained) if gained else "none",
                              len(labelled[child][1]), cdr, fr,
                              rs_counts["R"], rs_counts["S"],
                              ",".join(f"{x}:{reg[x]}" for x in REGION_NAMES if reg[x])
                              or "none"])

        if is_multi:
            trees.append({"clone": clone, "reads": len(unmut)
                          + sum(len(rs) for _, rs in subs.values())})

        newicks.append((clone, allele, len(subs), to_newick(edges, f"{clone}_germline")))
        fasta.append((f"{clone}_germline",
                      f"clonotype={clone} allele={allele} role=root mut=0 "
                      f"reads={len(unmut)}", germ))
        for label, (sig, rs) in labelled.items():
            fasta.append((label, f"clonotype={clone} allele={allele} role=node "
                                 f"mut={len(sig)} reads={len(rs)}",
                          apply_mutations(germ, sig)))

    new_cols = ["subclone", "subclone_size", "subclone_mut", "subclone_relation",
                "subclone_mutations"]
    out_cols = [c for c in header if c not in new_cols]
    pos = out_cols.index(skey) + 1 if skey in out_cols else 1
    out_cols = out_cols[:pos] + new_cols + out_cols[pos:]
    with open(out_tsv, "w") as f:
        f.write("\t".join(out_cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, ".")) for c in out_cols) + "\n")

    with open(out_nodes, "w") as f:
        f.write("\t".join(["node", "clonotype", "v_allele", "parent", "is_root", "n_reads",
                           "mut_from_germline", "relation", "isotype"]
                          + [f"mut_{x}" for x in REGION_NAMES] + ["R", "S"]) + "\n")
        for n in nodes:
            f.write("\t".join(str(x) for x in n) + "\n")

    with open(out_edges, "w") as f:
        f.write("\t".join(["clonotype", "v_allele", "parent", "child", "n_mutations_gained",
                           "mutations_gained", "child_reads", "cdr_gained", "fr_gained",
                           "R_gained", "S_gained", "regions_gained"]) + "\n")
        for e in edges_all:
            f.write("\t".join(str(x) for x in e) + "\n")

    with open(out_newick, "w") as f:
        f.write("# one tree per clonotype; branch lengths are mutations gained\n")
        f.write("# clonotype\tv_allele\tn_subclones\tnewick\n")
        for clone, allele, n, nw in newicks:
            f.write(f"{clone}\t{allele}\t{n}\t{nw}\n")

    with open(out_fasta, "w") as f:
        for name, desc, seq in fasta:
            f.write(f">{name} {desc}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    n_pages = 0
    if plot_top:
        n_pages, why = plot_lineage(trees, prefix, out_pdf, plot_top)
        if n_pages is None:
            # printed even when quiet: the user asked for a figure and is not getting one
            print(f"lineage plot skipped: {why}")
            n_pages, out_pdf = 0, None
        elif n_pages == 0:
            out_pdf = None

    say(f"TCR / unassigned reads   : {n_other} (excluded; BCR only)")
    say(f"clonotypes examined      : {len(byc)}")
    say(f"clonotypes with >1 subclone: {n_multi}")
    say(f"nodes written            : {len(nodes)}   edges: {len(edges_all)}   trees: {len(newicks)}")
    if not is_gapped:
        say("germline is NOT IMGT-gapped: region columns are all zero")
    if via:
        say("germline via D-paralog   : " + ", ".join(f"{k} ({v})" for k, v in via.most_common()))
    if skipped:
        say("skipped                  : " + ", ".join(f"{v} ({k})" for k, v in skipped.most_common()))
    say(f"\n  {out_nodes}\n  {out_edges}\n  {out_newick}\n  {out_fasta}\n  {out_tsv}")
    if out_pdf:
        say(f"  {out_pdf}  ({n_pages} trees)")
    return {"clones": len(byc), "clones_multi_subclone": n_multi, "n_nodes": len(nodes),
            "n_edges": len(edges_all), "n_trees": len(newicks), "n_plotted": n_pages,
            "out_nodes": out_nodes, "out_edges": out_edges, "out_newick": out_newick,
            "out_fasta": out_fasta, "out_tsv": out_tsv, "out_pdf": out_pdf,
            "rows": rows}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("corrected_tsv", help="corrected table from clone_correct_v2")
    p.add_argument("-g", "--germline", required=True,
                   help="germline V reference: a FASTA file or a directory of them; "
                        "IMGT-gapped enables the FR/CDR annotation")
    p.add_argument("--out-prefix")
    p.add_argument("--min-node-reads", type=int, default=MIN_NODE_READS,
                   help=f"reads a subclone needs to exist (default {MIN_NODE_READS})")
    p.add_argument("--min-mut-reads", type=int, default=MIN_MUT_READS,
                   help=f"reads a mutation needs before it can define a subclone "
                        f"(default {MIN_MUT_READS})")
    p.add_argument("--min-step-mutations", type=int, default=MIN_STEP_MUTATIONS,
                   help=f"new mutations a node must add over its ancestor; below this "
                        f"it merges into that ancestor (default {MIN_STEP_MUTATIONS})")
    p.add_argument("--core-frac", type=float, default=CORE_FRAC)
    p.add_argument("--multi-only", action="store_true",
                   help="write tree files only for clonotypes with >1 subclone")
    p.add_argument("--plot-top", type=int, default=20,
                   help="draw the N largest multi-subclone trees to "
                        "<prefix>.lineage_trees.pdf; 0 skips the plot (default 20)")
    p.add_argument("-q", "--quiet", action="store_true")
    a = p.parse_args()
    run_subclone_tree(a.corrected_tsv, a.germline, out_prefix=a.out_prefix,
                      min_node_reads=a.min_node_reads, min_mut_reads=a.min_mut_reads,
                      min_step=a.min_step_mutations, core_frac=a.core_frac,
                      multi_only=a.multi_only, plot_top=a.plot_top,
                      verbose=not a.quiet)


if __name__ == "__main__":
    main()
