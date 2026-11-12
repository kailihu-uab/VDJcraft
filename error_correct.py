#!/usr/bin/env python3

import argparse
import collections
import re

import numpy as np

SEQ_COLS = ["v_seq", "d_seq", "j_seq", "cdr_sequence"]
CDR3_MAX_DIVERGENCE = 0.15
CDR3_MIN_OVERLAP = 24
STRICT_CDR3_LENGTH = False

# --- alignment, orientation, consensus primitives ---------------------------
MIN_HAP_READS = 3

MIN_HAP_FRAC = 0.05

INDEL_MIN_READS = 5

INDEL_MIN_FRAC = 0.10

SHM_MIN_NESTING = 0.25

MAX_ASSIGN_DIST = 5

COV_FRAC = 0.5

KMER = 12

MATCH, MISMATCH, GAP_OPEN, GAP_EXT = 2, -4, -5, -1

NEG = -10 ** 8

COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")

def rc(s):
    return s.translate(COMP)[::-1]

def gene_of(allele):
    return re.sub(r"_\d+$", "", allele)

def align(ref, read):
 
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
            best, bi, bj, bs = int(mat[i, m], ), i, m, name

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

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

def is_hotspot(ref, pos):
    W, R, Y = set("AT"), set("AG"), set("CT")
    b = ref[pos]
    if b == "C" and pos >= 2 and ref[pos - 2] in W and ref[pos - 1] in R:
        return True
    if b == "G" and pos + 2 < len(ref) and ref[pos + 1] in Y and ref[pos + 2] in W:
        return True
    return False

def shm_score(ref, sigs):
  
    var = collections.Counter()
    for sig, w in sigs:
        for v in sig:
            var[v] += w
    rec = {v: c for v, c in var.items() if c >= 3}
    subs = [v for v in rec if v[0] == "sub"]
    nz = [s for s, _ in sigs if s]
    nested = comp = 0
    for i in range(len(nz)):
        for j in range(i + 1, len(nz)):
            comp += 1
            if nz[i] <= nz[j] or nz[j] <= nz[i]:
                nested += 1
    nest = nested / comp if comp else 0.0
    stats = {"n_rec_sub": len(subs), "titv": 0.0, "hotspot": 0.0, "nesting": nest}
    if len(subs) < 5:
        return False, stats
    ti = sum(1 for v in subs if (v[2], v[3]) in TRANSITIONS)
    tv = len(subs) - ti
    stats["titv"] = ti / tv if tv else float(ti)
    stats["hotspot"] = sum(1 for v in subs if is_hotspot(ref, v[1])) / len(subs)
    flag = (stats["titv"] > 1.4 and stats["hotspot"] > 0.25
            and nest > SHM_MIN_NESTING)
    return flag, stats

def build_consensus(ref, members):
    """members: list of (alignment, weight). Returns per-ref-position pieces and span."""
    n = len(ref)
    total = sum(w for _, w in members)
    cov = collections.Counter()
    base_votes = collections.defaultdict(collections.Counter)
    ins_votes = collections.defaultdict(collections.Counter)
    for al, w in members:
        s, e = al["span"]
        for i in range(s, e):
            cov[i] += w
            base_votes[i][al["cols"].get(i, "")] += w
        for i in range(s, e + 1):
            ins_votes[i][al["ins"].get(i, "")] += w

    keep = [i for i in range(n) if cov[i] >= COV_FRAC * total]
    if not keep:
        return None, None
    span = (keep[0], keep[-1] + 1)
    pieces = {}
    for i in range(span[0], span[1]):
        ins = ins_votes[i].most_common(1)[0][0] if ins_votes[i] else ""
        base = base_votes[i].most_common(1)[0][0] if base_votes[i] else ref[i]
        pieces[i] = ins + base
    return pieces, span

def phase(ref, uniq):
    """uniq: list of (seq, count, alignment, signature). Returns list of haplotypes."""
    total = sum(c for _, c, _, _ in uniq)
    by_sig = collections.defaultdict(lambda: {"count": 0, "members": []})
    for seq, cnt, al, sig in uniq:
        g = by_sig[sig]
        g["count"] += cnt
        g["members"].append((al, cnt))

    var = collections.Counter()
    for sig, g in by_sig.items():
        for v in sig:
            var[v] += g["count"]

    called = []
    for sig, g in by_sig.items():
        n = g["count"]
        if n < MIN_HAP_READS:
            continue
        if not sig:
            called.append(sig)                      # the seed haplotype itself
        elif len(sig) >= 2:
            called.append(sig)                      # co-occurring block: coincidence implausible
        elif sig and next(iter(sig))[0] == "sub":
            if n >= MIN_HAP_FRAC * total:           # systematic substitution error is rare in CCS
                called.append(sig)
        elif n >= INDEL_MIN_READS and n >= INDEL_MIN_FRAC * total:
            called.append(sig)                      # lone indel: guard against homopolymer artefact
    if not called:
        called = [max(by_sig, key=lambda s: by_sig[s]["count"])]

    # Only absorb a read into a haplotype it plausibly belongs to. CCS error runs 0-3
    # variants per read; anything further is a genuinely different sequence (rare allele,
    # paralog, chimera) that lacks the depth to stand alone. Rewriting it would fabricate
    # data, so it is left untouched and flagged instead.
    assign = {}
    for sig in by_sig:
        best = min(called, key=lambda c: len(sig ^ c))
        assign[sig] = best if len(sig ^ best) <= MAX_ASSIGN_DIST else None

    haps = collections.defaultdict(list)
    for seq, cnt, al, sig in uniq:
        if assign[sig] is not None:
            haps[assign[sig]].append((al, cnt, seq, sig))
    return haps, assign

STOPS = {"TAA", "TAG", "TGA"}

def longest_orf(s):
    best = 0
    for f in range(3):
        run = 0
        for i in range(f, len(s) - 2, 3):
            if s[i:i + 3] in STOPS:
                run = 0
            else:
                run += 1
                best = max(best, run)
    return best

def is_coding(s):
    """A V exon translates without stops on the coding strand; its complement does not."""
    return longest_orf(s) >= longest_orf(rc(s))

def orient(seqs_with_counts, prefer_coding=False):
    """Pick modal sequence as seed, flip reads whose reverse complement matches better."""
    mode = max(seqs_with_counts, key=lambda kv: kv[1])[0]
    if prefer_coding and not is_coding(mode):
        mode = rc(mode)
    kset = {mode[i:i + KMER] for i in range(len(mode) - KMER + 1)}

    def score(s):
        if len(s) < KMER:
            return 0
        return sum(1 for i in range(len(s) - KMER + 1) if s[i:i + KMER] in kset)

    out = {}
    for s, _ in seqs_with_counts:
        r = rc(s)
        out[s] = (s, "+") if score(s) >= score(r) else (r, "-")
    return out

def correct_segment(rows, seq_col, group_key, label, qc, fasta, tag, anchor="orf"):
    """Correct one segment column in place, writing per-read results onto each row."""
    groups = collections.defaultdict(list)
    for r in rows:
        s = r[seq_col]
        if s in (".", ""):
            continue
        k = group_key(r)
        if k is None:
            continue
        groups[k].append(r)

    for gname, members in sorted(groups.items()):
        counts = collections.Counter(r[seq_col] for r in members)
        flips = orient(list(counts.items()), prefer_coding=(anchor == "orf"))
        # A J or junction segment is too short to read a frame from, but it came off the
        # same molecule as its V, so anchor it to that read's already-resolved V strand.
        if anchor == "read_strand":
            agree = dis = 0
            for r in members:
                if r["v_strand"] in ("+", "-"):
                    if flips[r[seq_col]][1] == r["v_strand"]:
                        agree += 1
                    else:
                        dis += 1
            if dis > agree:
                flips = {k: (rc(v), "-" if s == "+" else "+") for k, (v, s) in flips.items()}
        for r in members:
            r[tag + "_strand"] = flips[r[seq_col]][1]
            r["_" + seq_col] = flips[r[seq_col]][0]

        oriented = collections.Counter(r["_" + seq_col] for r in members)
        if len(members) < MIN_HAP_READS:
            for r in members:
                r[seq_col] = r["_" + seq_col]
                r["_flags"].append(f"{tag}:low_depth({len(members)})")
            qc.append(f"{label} {gname}: n={len(members)} SKIPPED (depth < {MIN_HAP_READS})")
            continue

        seed = oriented.most_common(1)[0][0]
        uniq = []
        for seq, cnt in oriented.items():
            al = align(seed, seq) if seq != seed else align(seed, seed)
            if al is None:
                continue
            uniq.append((seq, cnt, al, signature(seed, al)))
        if not uniq:
            continue

        shm, st = shm_score(seed, [(sig, c) for _, c, _, sig in uniq])
        if shm:
            for r in members:
                r[seq_col] = r["_" + seq_col]
                r["_flags"].append(f"{tag}:shm_suspected")
            qc.append(f"{label} {gname}: n={len(members)} NOT COLLAPSED (SHM-like: "
                      f"Ti/Tv={st['titv']:.2f} hotspot={st['hotspot']:.0%} "
                      f"nesting={st['nesting']:.0%})")
            continue

        haps, assign = phase(seed, uniq)
        seq2hap = {}
        hap_ids = {}
        for k, (sig, mem) in enumerate(sorted(haps.items(), key=lambda kv: -sum(c for _, c, _, _ in kv[1])), 1):
            hid = f"{gname}|hap{k}"
            hap_ids[sig] = hid
            pieces, span = build_consensus(seed, [(al, c) for al, c, _, _ in mem])
            if pieces is None:
                continue
            support = sum(c for _, c, _, _ in mem)
            cons = "".join(pieces[i] for i in range(span[0], span[1]))
            fasta.append((hid, support, cons))
            for al, c, seq, sig2 in mem:
                s, e = al["span"]
                a, b = max(s, span[0]), min(e, span[1])
                seq2hap[seq] = (hid, support, "".join(pieces[i] for i in range(a, b)) if a < b else cons)
            qc.append(f"{label} {gname}|hap{k}: support={support} len={len(cons)} "
                      f"variants={len(sig)} distinct_input={len(mem)}")

        for r in members:
            o = r["_" + seq_col]
            if o not in seq2hap:
                r[seq_col] = o
                r["_flags"].append(f"{tag}:too_divergent_left_unchanged")
                continue
            hid, support, corrected = seq2hap[o]
            r[tag + "_haplotype"] = hid
            r[tag + "_support"] = str(support)
            r[tag + "_edits"] = str(0 if corrected == o else 1)
            if corrected != o:
                r["_" + tag + "_changed"] = True
            r[seq_col] = corrected

# --- loading, strand normalisation, read splicing ----------------------------
def load(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) != len(header):
                continue
            r = dict(zip(header, p))
            r["_flags"] = []
            for c in SEQ_COLS:
                r[c] = r.get(c, ".").upper().replace("-", "")
            # exactly as the input carried it, in the read's own orientation --
            # needed to locate each segment inside the original read FASTA
            r["_raw"] = {c: r[c] for c in SEQ_COLS}
            rows.append(r)
    return header, rows

def orient_column(rows, col, keyfn, prefer_coding=False, anchor_to_strand=False):
 
    groups = collections.defaultdict(list)
    for r in rows:
        if r[col] in (".", ""):
            continue
        k = keyfn(r)
        if k is not None:
            groups[k].append(r)
    flipped = 0
    for mem in groups.values():
        counts = collections.Counter(r[col] for r in mem)
        flips = orient(list(counts.items()), prefer_coding=prefer_coding)
   
        if anchor_to_strand:
            agree = tot = 0
            for r in mem:
                if r.get("strand") in ("+", "-"):
                    tot += 1
                    agree += flips[r[col]][1] == r["strand"]
            if tot and agree / tot < 0.15:
                flips = {k: (rc(v), "-" if s == "+" else "+") for k, (v, s) in flips.items()}
        for r in mem:
            new, st = flips[r[col]]
            if new != r[col]:
                flipped += 1
            r[col] = new
            if col == "v_seq":
                r["strand"] = st
    return flipped

def assign_strands(rows):
    vkey = lambda r: gene_of(r["v_gene"]) if r["v_gene"] not in (".", "") else None
    vjkey = lambda r: ((gene_of(r["v_gene"]), gene_of(r["j_gene"]))
                       if r["v_gene"] not in (".", "") and r["j_gene"] not in (".", "") else None)
    stats = {"v_seq": orient_column(rows, "v_seq", vkey, prefer_coding=True)}
    for c in ("d_seq", "j_seq"):
        stats[c] = orient_column(rows, c, vjkey, anchor_to_strand=True)
    stats["cdr_sequence"] = orient_column(rows, "cdr_sequence", vjkey)
    for r in rows:
        r.setdefault("strand", ".")
    return stats

def cdr3_similar(a, b, max_divergence=CDR3_MAX_DIVERGENCE,
                 min_overlap=CDR3_MIN_OVERLAP, strict_length=STRICT_CDR3_LENGTH):
    """Identity over the aligned overlap, tolerant of shifted CDR3 windows."""
    if strict_length and len(a) != len(b):
        return False
    if a == b:
        return True
    al = align(a, b)
    if al is None:
        return False
    s, e = al["span"]
    ov = e - s
    if ov < min_overlap or ov < 0.5 * min(len(a), len(b)):
        return False
    diffs = 0
    for v in signature(a, al):
        diffs += 1 if v[0] == "sub" else len(v[2]) if v[0] == "ins" else v[2]
    return diffs / ov <= max_divergence

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs

def _match_orientation(original, corrected):
    """The pipeline strand-normalises segments, so a corrected segment may be the
    reverse complement of the original it came from. Return it facing the original."""
    if original == corrected:
        return corrected
    r = rc(corrected)
    a, b = align(original, corrected), align(original, r)
    na = len(signature(original, a)) if a else 10 ** 6
    nb = len(signature(original, b)) if b else 10 ** 6
    return corrected if na <= nb else r

def splice_read(read, segments):

    up = read.upper()
    repl, taken = {}, set()
    placed = missing = clipped = 0
    for orig, corr in segments:
        if not orig or orig in (".", ""):
            continue
        corr = _match_orientation(orig, corr)
        i = up.find(orig)
        if i >= 0:
            ro, rcorr = orig, corr
        else:
            ro = rc(orig)
            i = up.find(ro)
            if i < 0:
                missing += 1
                continue
            rcorr = rc(corr)
        placed += 1
        if ro == rcorr:                       # unchanged: still claim the territory
            for off in range(len(ro)):
                taken.add(i + off)
            continue
        al = align(ro, rcorr)
        if al is None:
            continue
        a0, a1 = al["span"]
        if len(rcorr) > (a1 - a0):
            clipped += 1
        last_claimed = None
        for off in range(len(ro)):
            p = i + off
            if p in taken:
                continue
            taken.add(p)
            last_claimed = p
            base = al["cols"].get(off, ro[off])
            repl[p] = al["ins"].get(off, "") + base
        # an insertion restored at the segment's 3' end is keyed one past the last
        # base, so it would otherwise be dropped
        tail = al["ins"].get(len(ro), "")
        if tail and last_claimed is not None:
            repl[last_claimed] = repl[last_claimed] + tail
    if not repl:
        return read, 0, placed, missing, clipped
    out = [repl.get(k, ch) for k, ch in enumerate(read)]
    changed = sum(1 for k, v in repl.items() if v != read[k])
    return "".join(out), changed, placed, missing, clipped

# --- clonotype clustering and correction -------------------------------------
CDR_ALIASES = ["cdr_sequence", "cdr_ntseq", "cdr3_ntseq", "cdr3_nt", "cdr3_sequence",
               "junction", "junction_nt"]

MIN_CLONOTYPE_READS = 10       # support<=2 only means "rare" once depth allows it

MAX_ERROR_SUPPORT = 2          # variants carried by <= this many reads are error

MAX_ERROR_FREQ = 0.005         # ...or by under this fraction of the reads covering them

MAX_INDEL_LEN = 5              # only SMALL indels are ever treated as error

def detect_cdr_column(header):
    for c in CDR_ALIASES:
        if c in header:
            return c
    raise SystemExit("no CDR3 nucleotide column found; expected one of: "
                     + ", ".join(CDR_ALIASES))

def cluster_clonotypes(rows, max_divergence=CDR3_MAX_DIVERGENCE,
                       min_overlap=CDR3_MIN_OVERLAP, strict_length=False):
 
    groups = collections.defaultdict(list)
    unassigned = []
    for r in rows:
        if (r["v_gene"] in (".", "") or r["j_gene"] in (".", "")
                or r["cdr_sequence"] in (".", "")):
            unassigned.append(r)
            continue
        groups[(gene_of(r["v_gene"]), gene_of(r["j_gene"]))].append(r)

    clonotypes = []
    for key, mem in groups.items():
        counts = collections.Counter(r["cdr_sequence"] for r in mem)
        centroids = []
        for cdr, _ in counts.most_common():
            hit = next((c for c in centroids
                        if cdr3_similar(c[0], cdr, max_divergence, min_overlap,
                                        strict_length)), None)
            if hit is None:
                centroids.append((cdr, []))
                hit = centroids[-1]
            hit[1].append(cdr)
        by_cdr = collections.defaultdict(list)
        for r in mem:
            by_cdr[r["cdr_sequence"]].append(r)
        for cdr, members in centroids:
            clonotypes.append((key, cdr, [r for c in members for r in by_cdr[c]]))

    clonotypes.sort(key=lambda c: -len(c[2]))
    for i, (key, cdr, reads) in enumerate(clonotypes, 1):
        for r in reads:
            r["clonotype"] = f"clonotype{i}"
            r["clonotype_size"] = str(len(reads))
    for r in unassigned:
        r["clonotype"], r["clonotype_size"] = ".", "0"
        r["_flags"].append("no_clonotype_key")
    return clonotypes

def splice_read_all(read, segments):
  
    up = read.upper()
    repl, owner = {}, {}
    placed = missing = conflicts = 0
    for name, orig, corr in segments:
        if not orig or orig in (".", ""):
            continue
        corr = _match_orientation(orig, corr)
        i = up.find(orig)
        if i >= 0:
            ro, rcorr = orig, corr
        else:
            ro = rc(orig)
            i = up.find(ro)
            if i < 0:
                missing += 1
                continue
            rcorr = rc(corr)
        placed += 1
        if ro == rcorr:
            continue                      # nothing to write, and nothing reserved
        al = align(ro, rcorr)
        if al is None:
            continue
        last = None
        for off in range(len(ro)):
            p = i + off
            new = al["ins"].get(off, "") + al["cols"].get(off, ro[off])
            if new == ro[off]:
                continue                  # unchanged base: leave the position free
            if p in repl:
                if repl[p] != new:
                    conflicts += 1
                continue                  # first (higher priority) writer keeps it
            repl[p] = new
            owner[p] = name
            last = p
        tail = al["ins"].get(len(ro), "")
        if tail and last is not None:
            repl[last] = repl[last] + tail
    if not repl:
        return read, 0, placed, missing, conflicts
    out = [repl.get(k, ch) for k, ch in enumerate(read)]
    changed = sum(1 for k, v in repl.items() if v != read[k])
    return "".join(out), changed, placed, missing, conflicts

def correct_segment(reads, col, max_error_support=MAX_ERROR_SUPPORT,
                    min_reads=MIN_CLONOTYPE_READS, max_indel_len=MAX_INDEL_LEN,
                    max_error_freq=MAX_ERROR_FREQ):

    present = [r for r in reads if r[col] not in (".", "")]
    if len(present) < min_reads:
        for r in present:
            r["_flags"].append(f"{col}:low_depth({len(present)})")
        return {}, collections.Counter()

    counts = collections.Counter(r[col] for r in present)
    seed = counts.most_common(1)[0][0]
    aligns = {}
    for s in counts:
        al = align(seed, s)
        if al:
            aligns[s] = (al, signature(seed, al))

    support = collections.Counter()
    votes = collections.defaultdict(collections.Counter)
    for s, (al, sig) in aligns.items():
        w = counts[s]
        for v in sig:
            support[v] += w
        for i, b in al["cols"].items():
            votes[i][b] += w
    consensus = {i: c.most_common(1)[0][0] for i, c in votes.items()}
    depth = {i: sum(c.values()) for i, c in votes.items()}

    def coverage(v):
        p = v[1]
        return depth.get(p) or depth.get(p - 1) or sum(counts.values())

    fixed = {}
    skipped_large = {}          # keyed by sequence, so only that read is flagged
    for s, (al, sig) in aligns.items():
        # A long singleton indel is not a plausible CCS miscall -- more likely a
        # chimera or a misplaced boundary -- so it is left alone rather than rewritten.
        bad = []
        for v in sig:
            if not (support[v] <= max_error_support
                    or support[v] / coverage(v) < max_error_freq):
                continue
            size = len(v[2]) if v[0] == "ins" else v[2] if v[0] == "del" else 0
            if size > max_indel_len:
                skipped_large[s] = skipped_large.get(s, 0) + 1
                continue
            bad.append(v)
        if not bad:
            continue
        cols, ins = dict(al["cols"]), dict(al["ins"])
        for v in bad:
            if v[0] == "sub":
                cols[v[1]] = consensus.get(v[1], seed[v[1]])
            elif v[0] == "del":
                for k in range(v[1], v[1] + v[2]):
                    cols[k] = consensus.get(k, seed[k] if k < len(seed) else "")
            else:
                ins.pop(v[1], None)
        lo, hi = al["span"]

        def render(colmap, insmap):
            out = []
            for i in range(lo, hi):
                out.append(insmap.get(i, ""))
                out.append(colmap.get(i, ""))
            out.append(insmap.get(hi, ""))
            return "".join(out)

        # The alignment only covers the part of the read that maps to the seed. Splice
        # the rebuilt region back into the read so bases beyond the seed's extent are
        # kept -- rendering the span alone would silently truncate the read.
        new_span = render(cols, ins)
        old_span = render(al["cols"], al["ins"])
        at = s.find(old_span)
        new = (s[:at] + new_span + s[at + len(old_span):]) if at >= 0 else new_span
        if new != s:
            fixed[s] = (new, len(bad))

    for r in present:
        got = fixed.get(r[col])
        if got:
            r[col + "_nfix"] = str(got[1])
            r[col] = got[0]
            r["_changed"] = True
        if skipped_large.get(r[col]):
            r["_flags"].append(f"{col}:large_indel_left_unchanged")
    return fixed, support

def run_correction(infile, out_prefix=None, query_fasta=None,
                   out_tsv=None, out_fix=None, out_qc=None, out_reads=None,
                   out_consensus=None,
                   max_error_support=MAX_ERROR_SUPPORT,
                   min_reads=MIN_CLONOTYPE_READS, max_indel_len=MAX_INDEL_LEN,
                   max_error_freq=MAX_ERROR_FREQ,
                   cdr3_max_divergence=CDR3_MAX_DIVERGENCE,
                   cdr3_min_overlap=CDR3_MIN_OVERLAP, strict_cdr3_length=False,
                   include_unprocessed=True, verbose=True):
    """Cluster clonotypes, correct singleton-support variants, write all outputs."""
    if out_prefix is None:
        out_prefix = infile.rsplit(".", 1)[0]
    out_tsv = out_tsv or f"{out_prefix}.corrected.tsv"
    out_fix = out_fix or f"{out_prefix}.corrections.tsv"
    out_qc = out_qc or f"{out_prefix}.qc_report.txt"
    out_reads = out_reads or f"{out_prefix}.corrected_reads.fasta"
    say = print if verbose else (lambda *a, **k: None)

    header, rows = load(infile)
    cdr_col = detect_cdr_column(header)
    say(f"loaded {len(rows)} reads from {infile}   (CDR3 column: {cdr_col})")
    if cdr_col != "cdr_sequence":
        for r in rows:
            r["cdr_sequence"] = r.get(cdr_col, ".").upper().replace("-", "")
            r["_raw"]["cdr_sequence"] = r["cdr_sequence"]
    # the input may already carry a strand column of its own; keep it and report the
    # orientation this run derived under a separate name
    had_strand = {id(r): r.get("strand") for r in rows} if "strand" in header else {}
    ostats = assign_strands(rows)
    for r in rows:
        r["v_strand"] = r.get("strand", ".")
        if id(r) in had_strand:
            r["strand"] = had_strand[id(r)]
    say("  reoriented per column: " + "  ".join(f"{k}={v}" for k, v in ostats.items()))
    for r in rows:
        r["_orig"] = {c: r[c] for c in SEQ_COLS}

    clonotypes = cluster_clonotypes(rows, cdr3_max_divergence, cdr3_min_overlap,
                                    strict_cdr3_length)
    sizes = [len(c[2]) for c in clonotypes]
    say(f"clonotypes: {len(clonotypes)}  largest={sizes[0] if sizes else 0}")

    strict = {(gene_of(r["v_gene"]), gene_of(r["j_gene"]), len(r["cdr_sequence"]))
              for r in rows if r.get("clonotype", ".") != "."}
    # reads without a V gene, J gene or CDR3 cannot be put in a clonotype, so nothing
    # was corrected in them; they are left out of every output rather than written
    # with an empty clonotype
    out_rows = [r for r in rows if r["clonotype"] != "."]
    n_dropped = len(rows) - len(out_rows)
    say(f"reads without V, J or CDR3 (not written): {n_dropped}")

    qc, kept = [], collections.Counter()
    for i, (key, cdr, reads) in enumerate(clonotypes, 1):
        name = f"clonotype{i}"
        if len(reads) < min_reads:
            for r in reads:
                r["_flags"].append(f"clonotype_too_small({len(reads)})")
            continue
        line = [f"{name}: {key[0]} / {key[1]}  n={len(reads)}  CDR3 len={len(cdr)}"]
        for col in SEQ_COLS:
            fixed, support = correct_segment(reads, col, max_error_support, min_reads,
                                             max_indel_len, max_error_freq)
            if support:
                single = sum(1 for v, c in support.items() if c <= max_error_support)
                multi = len(support) - single
                kept[col] += multi
                if single or multi:
                    line.append(f"    {col}: {len(support)} variants "
                                f"-> {single} corrected, {multi} kept (>= 2 reads)")
        qc.append("\n".join(line))

    # ---- corrected table -------------------------------------------------------
    lead = [c for c in ["clonotype", "clonotype_size", "Query",
                        "v_gene", "d_gene", "j_gene", "c_gene"]
            if c in ("clonotype", "clonotype_size") or c in header]
    name_of = {c: (cdr_col if c == "cdr_sequence" else c) for c in SEQ_COLS}

    plain = [name_of[c] for c in SEQ_COLS]
    cor = [name_of[c] + "_corrected" for c in SEQ_COLS]
    org = [name_of[c] + "_original" for c in SEQ_COLS]
    meta = ["v_strand"] + [name_of[c] + "_nfix" for c in SEQ_COLS]
    back = {name_of[c]: c for c in SEQ_COLS}
    # an input that is itself a corrected table carries the old clone columns; drop them
    extras = [c for c in header if c not in lead and c not in SEQ_COLS and c != cdr_col
              and c not in ("clone", "clone_size", "clonotype", "clonotype_size")]
    out_cols = lead + plain + cor + org + meta + extras + ["flags"]

    def value(r, c):
        for suf, src in (("_corrected", None), ("_original", "_orig"), ("_nfix", None)):
            if c.endswith(suf):
                key = back.get(c[:-len(suf)], c[:-len(suf)])
                if src:
                    return r["_orig"].get(key, ".")
                return r.get(key + "_nfix", ".") if suf == "_nfix" else r.get(key, ".")
        if c in back:

            return r.get(back[c], ".")
        return r.get(c, ".")

    with open(out_tsv, "w") as f:
        f.write("\t".join(out_cols) + "\n")
        for r in out_rows:
            r["flags"] = ";".join(r["_flags"]) if r["_flags"] else "."
            f.write("\t".join(str(value(r, c)) for c in out_cols) + "\n")

    # ---- the specific corrections ---------------------------------------------
    n_fix = 0
    with open(out_fix, "w") as f:
        f.write("Query\tclonotype\tclonotype_size\tsegment\tn_variants_reverted\t"
                "len_original\tlen_corrected\toriginal_seq\tcorrected_seq\n")
        for r in out_rows:
            for c in SEQ_COLS:
                if c + "_nfix" not in r:
                    continue
                o, n = r["_orig"][c], r[c]
                f.write("\t".join(str(x) for x in [
                    r.get("Query", "."), r.get("clonotype", "."), r.get("clonotype_size", "."),
                    name_of.get(c, c), r[c + "_nfix"], len(o), len(n), o, n]) + "\n")
                n_fix += 1

    # ---- updated read FASTA ----------------------------------------------------
    fa_stats = None
    if query_fasta:
        fa = read_fasta(query_fasta)
        fa_stats = collections.Counter()
        seen = set()
        with open(out_reads, "w") as out:
            for r in out_rows:
                q = r.get("Query")
                if q in seen:
                    fa_stats["duplicate read id (written once)"] += 1
                    continue
                seq = fa.get(q)
                if seq is None:
                    fa_stats["read id not in query fasta"] += 1
                    continue
                seen.add(q)
                segs = [(c, r["_raw"][c], r[c]) for c in
                        ("v_seq", "j_seq", "d_seq", "cdr_sequence")
                        if r["_raw"].get(c) not in (".", "", None)]
                new, nch, _, missing, conf = splice_read_all(seq, segs)
                fa_stats["reads written"] += 1
                fa_stats["reads changed" if nch else "reads unchanged"] += 1
                fa_stats["segments not locatable in read"] += missing
                fa_stats["positions two segments both wanted to change"] += conf
                out.write(f">{q} status=corrected clonotype={r.get('clonotype','.')} "
                          f"corrected_bases={nch}\n")
                for i in range(0, len(new), 60):
                    out.write(new[i:i + 60] + "\n")
            # every other read -- no V, J or CDR3 to place it in a clonotype, or not in
            # the table at all -- is copied through unchanged, so the FASTA keeps the
            # full read set while the corrected table keeps only clonotyped reads
            if include_unprocessed:
                for q, seq in fa.items():
                    if q in seen:
                        continue
                    fa_stats["unprocessed reads copied through"] += 1
                    out.write(f">{q} status=unprocessed\n")
                    for i in range(0, len(seq), 60):
                        out.write(seq[i:i + 60] + "\n")

    if out_consensus:
        # one record per clonotype per segment: the modal corrected sequence, which
        # after correction is the clonotype consensus
        with open(out_consensus, "w") as f:
            for i, (key, cdr, reads) in enumerate(clonotypes, 1):
                if len(reads) < min_reads:
                    continue
                for col in SEQ_COLS:
                    vals = collections.Counter(r[col] for r in reads
                                               if r[col] not in (".", ""))
                    if not vals:
                        continue
                    seq, n = vals.most_common(1)[0]
                    f.write(f">clonotype{i}|{col} v={key[0]} j={key[1]} reads={len(reads)} "
                            f"support={n} len={len(seq)}\n")
                    for k in range(0, len(seq), 60):
                        f.write(seq[k:k + 60] + "\n")

    changed = sum(1 for r in out_rows if r.get("_changed"))
    flags = collections.Counter(x for r in out_rows for x in r["flags"].split(";") if x != ".")
    with open(out_qc, "w") as f:
        f.write(f"input: {infile}\n")
        f.write(f"reads in: {len(rows)}   reads out: {len(out_rows)}\n")
        f.write(f"  not written, no V gene, J gene or CDR3 to place them in a "
                f"clonotype: {n_dropped}\n")
        f.write(f"clonotypes: {len(clonotypes)}   largest: {sizes[0] if sizes else 0}\n")
        f.write(f"  reads in clonotypes of >={min_reads}: "
                f"{sum(s for s in sizes if s >= min_reads)}\n\n")
        f.write(f"indel size cap: indels longer than {max_indel_len}nt are never corrected\n")
        f.write(f"clonotype floor: only clonotypes of >= {min_reads} reads are corrected\n")
        f.write("RULE: a variant is treated as sequencing error and reverted to the\n"
                f"consensus when it is carried by <= {max_error_support} read(s) OR by under\n"
                f"{max_error_freq:.2%} of the reads covering its position. Anything above both\n"
                "thresholds is left untouched, on the basis that agreement between\n"
                "independent reads is evidence of real sequence.\n\n")
        f.write(f"reads with >=1 segment corrected: {changed}\n")
        f.write(f"segment corrections written     : {n_fix}\n")
        f.write("variants KEPT because >=2 reads carried them:\n")
        for k, v in kept.most_common():
            f.write(f"  {k}: {v}\n")
        f.write(f"\nCDR3 clustering: {'strict equal length' if strict_cdr3_length else 'similarity over aligned overlap'}"
                f" (max divergence {cdr3_max_divergence:.0%})\n")
        f.write(f"  clonotypes called this way        : {len(clonotypes)}\n")
        f.write(f"  clonotypes an equal-length rule gives: {len(strict)}\n")
        if fa_stats:
            f.write(f"\ncorrected-read FASTA: {out_reads}\n")
            for k, v in fa_stats.most_common():
                if v:
                    f.write(f"  {k}: {v}\n")
        f.write("\nflags:\n")
        for k, v in flags.most_common():
            f.write(f"  {k}: {v}\n")
        f.write("\n" + "\n\n".join(qc) + "\n")

    say(f"reads corrected: {changed}   segment corrections: {n_fix}")
    say(f"variants kept (>=2 reads): {dict(kept)}")
    say(f"wrote {out_tsv}, {out_fix}, {out_qc}" + (f", {out_reads}" if query_fasta else ""))
    return {"n_reads": len(out_rows), "n_dropped": n_dropped, "n_clonotypes": len(clonotypes),
            "largest": sizes[0] if sizes else 0, "reads_corrected": changed,
            "segment_corrections": n_fix, "variants_kept": dict(kept),
            "out_tsv": out_tsv, "out_corrections": out_fix, "out_qc": out_qc,
            "out_reads_fasta": out_reads if query_fasta else None, "rows": out_rows}


def run_clone_correction(infile, out_tsv=None, out_fasta=None, out_qc=None,
                         out_corrections=None, query_fasta=None,
                         out_reads_fasta=None, verbose=True, **kwargs):

    prefix = (out_tsv or infile).rsplit(".", 1)[0]
    return run_correction(infile, out_prefix=prefix, out_tsv=out_tsv,
                          out_fix=out_corrections, out_qc=out_qc,
                          out_consensus=out_fasta, out_reads=out_reads_fasta,
                          query_fasta=query_fasta, verbose=verbose, **kwargs)
