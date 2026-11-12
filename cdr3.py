#!/usr/bin/env python3
"""
CDR3 extraction from long reads, anchored on V/J alignment coordinates.

v6 is v5's algorithm packaged as an importable library.  The filename matters:
"cdr3_v6.py" cannot be imported at all, because dots are module separators
in Python, so `from cdr3.new.v5 import ...` is a syntax error and the only way
to reach the code was to run the file.  This module is named cdr3_v6.py, which
is a legal identifier, so it can simply be imported.

Importing has no side effects -- nothing runs until you call extract_cdr3().

What changed relative to v4
---------------------------
v4 motif-scanned the *entire* read (median ~1 kb, up to 5.6 kb) in all three
frames and, on ties, preferred whatever candidate was closest to 15 aa.  On
kappa reads -- where the true CDR3-IMGT is 9 aa -- that tiebreak actively
selects spurious distal hits, and nothing rejected junctions containing stop
codons.

v5 keeps the same motif logic but anchors it to the V/J alignment coordinates
that are already present in the annotation table:

  * The conserved Cys104 codon starts ~23 nt before the V alignment 3' end
    (486/838 reads exactly 23, essentially all within 9-27).
  * The J motif [WF]G.G starts ~7 nt into the J alignment (540/826 exactly 7).

Candidates are enumerated inside those windows, required to share a reading
frame, required to be stop-free, and then *scored* -- motif quality, positional
agreement with the alignment, and a weak locus-aware length prior.  Windows are
progressively widened only when the tight pass finds nothing, so detection rate
goes up without letting the loose passes outvote a well-anchored call.

Also fixed / added
------------------
  * BED coordinates were emitted in reverse-complemented space for Plus/Minus
    rows, so every minus-strand interval pointed at the wrong place in the
    original FASTA.  Now converted back to query coordinates, BED6, +/- strand.
  * Rows with an unknown strand (v_gene = NA) were dropped outright; they still
    carry J coordinates, so both orientations are now tried and scored.
  * Rows with j_gene = NA are handled V-anchored instead of falling back to a
    blind global scan.
  * Junctions with a stop codon or a frameshift are no longer translated as
    if they were productive; they are reported with frame_status set and
    productive = False, so nucleotide-level clonotyping keeps them and
    amino-acid analyses can filter them out.
  * Length prior is locus-aware (IGK 9, IGH 15, ...) and weighted low enough
    that it only breaks ties.
  * Read IDs are not unique in the FASTA (57666 records under 57401 IDs).
    A plain id -> sequence dict keeps whichever copy came last; when that is
    the shorter copy, the row's coordinates point past its end and the read
    becomes uncallable.  Every copy is now kept and the one matching the
    coordinates is chosen -- worth 18 of 998 rows on the test set.
  * Only the sequences named in the table are held in memory.

CDR3 vs JUNCTION
----------------
cdr_ntseq / cdr_aa are CDR3-IMGT: the anchors are EXCLUDED (IMGT 105-117).
junction_ntseq / junction_aa are the JUNCTION: anchors INCLUDED (104-118), so
junction_aa == "C" + cdr_aa + "W"/"F" for every productive call.

Measured on vjtemp.test.csv (998 rows, IGK), scoring a call correct only if
its reported CDR3 really does sit between a Y-x-C and a [WF]G-x-G in the read:

    v4   929 calls, 865 correct (93.1%);  BED interval wrong on 563 of them
    v5   927 calls, 890 correct (96.0%);  BED interval correct on all 927
         of which confidence == "high": 880 calls, 880 correct (100.0%)

Of the 71 rows v5 leaves uncalled, none has a usable canonical anchor pair:
43 have no V alignment at all (J-only / germline reads, no junction to call),
13 have a Y-x-C but no [WF]G-x-G anywhere in the read, 8 have neither gene
assigned, and 7 have both motifs but never in frame and stop-free.

Output columns (v4 columns kept first, in order, for backward compatibility):
    Query, cdr_ntseq, cdr_aa, strand, frame, search_method,
    locus, v_gene, j_gene, cdr3_len, frame_status, productive,
    junction_aa, junction_ntseq, score, confidence, start_motif, end_motif,
    cdr_start, cdr_end            <- query coordinates, 0-based half-open

Filter on confidence == "high" for a precision-oriented set.

Usage from your own pipeline
----------------------------
    from cdr3_v6 import extract_cdr3

    # returns a DataFrame; writes files only if you name them
    calls = extract_cdr3("vjtemp.csv", "vdj_candidate.fa")

    calls = extract_cdr3(
        table="vjtemp.csv",              # path OR an in-memory DataFrame
        fasta="vdj_candidate.fa",        # path OR {id: seq} / {id: [seq, ...]}
        out_tsv="cdr3.tsv",              # omit to skip writing
        out_bed="cdr3.bed",              # omit to skip writing
        verbose=True,                    # progress to stderr; default quiet
    )

    good = calls[calls.confidence.eq("high") & calls.productive]
    print(calls.attrs["stats"])          # per-category counts

Wiring it into an existing main file as a flag:

    parser.add_argument("--cdr3", action="store_true")
    ...
    if args.cdr3:
        extract_cdr3(table=args.table, fasta=args.fasta,
                     out_tsv=args.out_cdr3, out_bed=args.out_bed,
                     verbose=True)

Standalone, unchanged:
    python cdr3_v6.py -t vjtemp.csv -f vdj_candidate.fa -o cdr3.tsv -b cdr3.bed

v5's entry point still works: cdr3_extract(df=..., fasta_file=...,
output_bed=..., output_cdr3=...) is kept as a thin alias.
"""

__all__ = ["extract_cdr3", "cdr3_extract", "read_fasta", "select_record"]

import argparse
import math
import sys
from collections import Counter

import pandas as pd
from Bio.Seq import Seq


# ============================================================
# Locus-specific expectations
# ============================================================

# Conserved residue closing CDR3-IMGT (position 118).
LOCUS_END_RESIDUE = {
    "IGH": "W",
    "IGK": "F",
    "IGL": "F",
    "TRA": "F",
    "TRB": "F",
    "TRG": "F",
    "TRD": "F",
}

# Typical CDR3-IMGT length in aa (C and W/F excluded).  Used only as a weak
# tiebreak -- never strong enough to override a coordinate-anchored call.
LOCUS_EXPECTED_LEN = {
    "IGH": 15,
    "IGK": 9,
    "IGL": 10,
    "TRA": 12,
    "TRB": 12,
    "TRG": 11,
    "TRD": 14,
}

DEFAULT_EXPECTED_LEN = 12

# FR4 consensus, as offsets in nt from the position-118 anchor codon.
# Canonical FR4 is [WF]118 G119 x120 G121 T122 L123 V124 T125 ...  Requiring
# the single T122 rejects reads whose J carries one substitution there --
# WGQG(I)LVT and WGQG(S)LVT are ordinary IGHJ with a point change, and on the
# reference FASTA that cost ~1000 genuine AR/AK heavy-chain reads.  Counting
# agreement across four positions keeps the specificity without the brittleness.
FR4_CONSENSUS = (
    (12, "TA"),        # 122
    (15, "LMKQTR"),    # 123
    (18, "VLI"),       # 124
    (21, "TVIE"),      # 125
)

FR4_MIN_SUPPORT = 2

# Measured on the reference set: the position-118 codon starts this far into
# the J alignment (540/826 reads exactly 7).
EXPECTED_J_OFFSET = 7

# qc values.  This axis records how much we trust the *extraction*; the
# separate frame_status column records the *biology* (productive / stop /
# out of frame).  A junction can be biologically non-productive and still be
# confidently located, which is why the two are kept apart.
QC_PASS = "pass"
QC_STOP = "stop_codon_present"
QC_LOW = "low_qual"


def fr4_support(frames, p_j):
    """How many FR4 consensus positions the J anchor at p_j agrees with."""
    hits = 0
    for offset, allowed in FR4_CONSENSUS:
        residue = frames.aa_at(p_j + offset)
        if residue and residue in allowed:
            hits += 1
    return hits


# Measured on the reference set: v_end - Cys_codon_start.
EXPECTED_C_OFFSET = 23
C_OFFSET_SIGMA = 14.0


# ============================================================
# Search passes, tried in order until one yields a candidate
#
#   c_lo/c_hi : Cys codon start window, relative to the V alignment 3' end
#   j_lo/j_hi : J motif codon start window, relative to the J alignment 5' end
#   use_v/use_j : consult that coordinate at all (v_anchored deliberately
#                 ignores a J coordinate that points somewhere implausible)
#   requires  : which coordinates must be present for the pass to run
#
# The "v_anchored" pass exists because the J alignment is not always the J of
# the rearrangement: on the reference set some reads align their J to a
# germline J ~1 kb downstream of the V.  The V 3' end stays reliable, so that
# pass keeps the Cys window and re-derives the J window from the Cys instead.
#
# "global" runs only when there is no positional evidence whatsoever.  Letting
# it run whenever the windowed passes came up empty is precisely the v4
# behaviour that wanders off to a distant [WF]G.G.
# ============================================================

SEARCH_PASSES = [
    dict(name="anchored",      use_v=True,  use_j=True,
         c_lo=-48, c_hi=+9,  j_lo=-21, j_hi=+78,  requires="any"),
    dict(name="anchored_wide", use_v=True,  use_j=True,
         c_lo=-90, c_hi=+30, j_lo=-45, j_hi=+150, requires="any"),
    dict(name="v_anchored",    use_v=True,  use_j=False,
         c_lo=-60, c_hi=+15, j_lo=-21, j_hi=+78,  requires="both"),
    # A minority of rows carry V/J coordinates that do not agree with the
    # read at all (their v_seq/j_seq do not occur in it).  Those need a
    # whole-read scan, but an unrestricted one is what made v4 wander.  This
    # pass therefore demands the full canonical signature on both sides --
    # Y-x-C, [WF]G-x-G, and the FR4 threonine that follows it -- which is
    # specific enough to survive a scan of several kb.
    dict(name="global_canonical", use_v=False, use_j=False,
         c_lo=0,   c_hi=0,   j_lo=0,   j_hi=0,    requires="always",
         canonical_only=True),
    # Last resort: a junction that is clearly present but not productive.
    # Roughly two thirds of rearrangements on the non-selected allele are
    # out of frame, so a real dataset contains them; v4 translated them
    # anyway, which is where its stop-codon-bearing CDR3s came from.  Both
    # anchors must be canonical, and the call is flagged, never hidden.
    dict(name="rescue_nonproductive", use_v=True, use_j=True,
         c_lo=-48, c_hi=+9,  j_lo=-21, j_hi=+78,  requires="any", relaxed=True),
    # Same last resort for reads with no coordinates at all.  The anchors here
    # are still real motifs -- only the frame and stop-free requirements are
    # relaxed -- so a stop-bearing or frameshifted junction gets reported and
    # labelled rather than silently dropped.
    dict(name="global_relaxed", use_v=False, use_j=False,
         c_lo=0,   c_hi=0,   j_lo=0,   j_hi=0,    requires="always",
         canonical_only=True, relaxed=True),
]


# ============================================================
# Small helpers
# ============================================================

def locus_of(gene):
    """IGKV4-1_01 -> IGK.  Returns None if unrecognised."""
    g = str(gene).upper()
    for locus in LOCUS_END_RESIDUE:
        if locus in g:
            return locus
    return None


def to_int(value):
    """Parse a coordinate cell, tolerating NA / blank / float-ish values."""
    try:
        if value is None:
            return None
        text = str(value).strip()
        if text == "" or text.upper() in ("NA", "NAN", "NONE", "."):
            return None
        return int(float(text))
    except (TypeError, ValueError):
        return None


class Frames:
    """Codon-level view of one oriented sequence.

    aa_at(p) is the residue of the codon starting at nucleotide p, so the
    reading frame is simply p % 3 and no per-frame bookkeeping is needed at
    the call sites.
    """

    def __init__(self, seq):
        self.seq = seq
        self.length = len(seq)
        self.prots = []

        for frame in range(3):
            trimmed = seq[frame:]
            trimmed = trimmed[: len(trimmed) // 3 * 3]
            if len(trimmed) < 3:
                self.prots.append("")
            else:
                self.prots.append(str(Seq(trimmed).translate()))

    def aa_at(self, p):
        if p < 0 or p + 3 > self.length:
            return None
        prot = self.prots[p % 3]
        i = p // 3
        return prot[i] if i < len(prot) else None

    def has_stop_between(self, p_start, p_end):
        """Stop codon in any codon starting in [p_start, p_end)?"""
        for p in range(p_start, p_end, 3):
            if self.aa_at(p) == "*":
                return True
        return False


# ============================================================
# Anchor enumeration
# ============================================================

def find_cys(frames, lo, hi):
    """Codon starts of TGT/TGC within [lo, hi)."""
    lo = max(0, lo)
    hi = min(hi, frames.length - 2)
    return [p for p in range(lo, hi) if frames.aa_at(p) == "C"]


def find_j_motifs(frames, lo, hi, allow_weak=False):
    """Codon starts of a plausible J anchor within [lo, hi).

    Returns (position, tier, motif_string).  Tiers, best first:
        2  [WF] G x G   -- canonical
        1  [WF] G x ?   or  [WF] ? x G   -- one G lost to a substitution
        0  [WF]         -- bare aromatic, no glycine at all

    [WF]G-x-G is one of the most conserved motifs in the J-REGION.  A bare
    aromatic with neither glycine is not evidence of a J anchor -- on the
    reference set those hits were uniformly spurious (FLCC, WKSN, WRSN ...) --
    so tier 0 is dropped unless `allow_weak` is set.
    """
    lo = max(0, lo)
    hi = min(hi, frames.length - 2)
    hits = []

    for p in range(lo, hi):
        a0 = frames.aa_at(p)
        if a0 not in ("W", "F"):
            continue

        a1 = frames.aa_at(p + 3)
        a2 = frames.aa_at(p + 6)
        a3 = frames.aa_at(p + 9)

        motif = "".join(x if x else "?" for x in (a0, a1, a2, a3))

        # A stop codon inside the anchor means the frame is wrong.
        if "*" in motif:
            continue

        if a1 == "G" and a3 == "G":
            tier = 2
        elif a1 == "G" or a3 == "G":
            tier = 1
        else:
            tier = 0
            if not allow_weak:
                continue

        hits.append((p, tier, motif))

    return hits


# ============================================================
# Scoring
# ============================================================

def anchors_canonical(cand):
    """Both anchors carry the full signature: Y-x-C and [WF]G-x-G."""
    return (
        cand["start_motif"][0] == "Y"
        and len(cand["end_motif"]) == 4
        and cand["end_motif"][1] == "G"
        and cand["end_motif"][3] == "G"
    )


def is_canonical(cand):
    """Canonical anchors *and* a clean productive reading frame."""
    return cand["productive"] and anchors_canonical(cand)


def rank_key(cand):
    """Order candidates: anchor quality, then productivity, then score."""
    return (anchors_canonical(cand), cand["productive"], cand["score"])

def score_candidate(frames, p_c, p_j, j_tier, end_residue, expected_len,
                    v_end, j_start, j_end):
    """Higher is better.  Returns (score, start_motif)."""
    score = 0.0

    # --- start motif context: the FR3 signature is Y-x-C ------------------
    m2 = frames.aa_at(p_c - 6)
    m1 = frames.aa_at(p_c - 3)
    start_motif = "".join(x if x else "?" for x in (m2, m1, "C"))

    if m2 == "Y":
        score += 3.0
        # AVYYC / YYC / YFC and friends
        if m1 in ("Y", "F", "H", "V", "L", "C"):
            score += 0.7
    elif m2 in ("F", "H", "W"):
        score += 2.0
    elif m1 in ("Y", "F", "H", "V", "A", "L"):
        score += 1.0
    else:
        score += 0.3

    # --- end motif --------------------------------------------------------
    a0 = frames.aa_at(p_j)

    # end_residue is None when the locus could not be identified; W and F are
    # then equally acceptable, and scoring them as "matched" keeps scores on
    # the same scale as a known locus (the confidence cut-offs depend on it).
    matched = True if end_residue is None else (a0 == end_residue)

    if j_tier == 2:
        score += 4.0 if matched else 2.8
    elif j_tier == 1:
        score += 1.5 if matched else 1.0
    else:
        score += 0.4 if matched else 0.2

    # canonical J continues ...G-x-G-T-x-[LV]
    if frames.aa_at(p_j + 12) in ("T", "A"):
        score += 0.5

    # --- positional agreement with the alignments -------------------------
    if v_end is not None:
        delta = (v_end - p_c) - EXPECTED_C_OFFSET
        score += 3.0 * math.exp(-((delta / C_OFFSET_SIGMA) ** 2))

    if j_start is not None:
        # The [WF]G.G codon sits inside the J alignment, near its 5' end.
        if j_start - 9 <= p_j <= (j_end - 9 if j_end is not None else j_start + 90):
            score += 3.0
        else:
            score += max(0.0, 2.0 - abs(p_j - (j_start + 7)) / 20.0)

    # --- weak length prior ------------------------------------------------
    n_aa = (p_j - p_c) // 3 - 1
    score -= 0.08 * abs(n_aa - expected_len)

    # --- flanking frame sanity -------------------------------------------
    if frames.has_stop_between(max(0, p_c - 30), p_c):
        score -= 1.0
    if frames.has_stop_between(p_j, p_j + 15):
        score -= 1.0

    return score, start_motif


def longest_homopolymer_runs(seq, min_run=3):
    """(start, length) of every homopolymer run of at least min_run."""
    runs = []
    i = 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        if j - i + 1 >= min_run:
            runs.append((i, j - i + 1))
        i = j + 1
    return runs


def translate_frameshifted(frames, p_c, p_j):
    """Translate a junction whose two anchors are not in the same frame.

    A 1-2 nt indel (the usual CCS homopolymer error) leaves the 5' half
    readable in the Cys frame and the 3' half readable in the J frame.  Rather
    than truncate at the shift -- which throws the 3' half away -- each half is
    translated in its own frame and the two are joined.

    No bases are invented: cdr_ntseq stays the exact read sequence and the
    skipped/overlapped nucleotides at the splice are reported in `repair`.

    Returns (cdr_aa, repair_description, confident).
    """
    start = p_c + 3
    region = frames.seq[start:p_j]
    shift = len(region) % 3

    runs = longest_homopolymer_runs(region)

    # Splice at the end of the longest homopolymer -- the likeliest indel site.
    if runs:
        run_start, run_len = max(runs, key=lambda r: r[1])
        split = start + run_start + run_len
    else:
        split = start + len(region) // 2

    left_end = start + ((split - start) // 3) * 3
    right_start = split + ((p_j - split) % 3)

    left = str(Seq(frames.seq[start:left_end]).translate()) if left_end > start else ""
    right = str(Seq(frames.seq[right_start:p_j]).translate()) if p_j > right_start else ""

    cdr_aa = left + right

    repair = "frameshift_%dnt@%d" % (shift, left_end - start)
    if runs:
        repair += "_homopolymer"

    # Confident only when the indel site is unambiguous, neither half carries a
    # stop, and the J side still reads as real FR4.
    confident = (
        len(runs) == 1
        and "*" not in cdr_aa
        and fr4_support(frames, p_j) >= 3
    )

    return cdr_aa, repair, confident


def rescue_by_position(frames, end_residue, v_end, j_start, j_end,
                       min_len, max_len):
    """Place CDR3 105-117 from IMGT offsets when the motif search found none.

    Cys104 starts ~EXPECTED_C_OFFSET nt before the V alignment 3' end and the
    position-118 codon ~EXPECTED_J_OFFSET nt into the J alignment.  When an
    anchor motif is absent -- a somatically mutated FR3, a truncated or
    divergent J -- the alignment coordinate still pins the boundary, so the
    CDR3 can be delimited by position instead of by motif.

    Everything returned here is flagged low_qual: a boundary placed by offset
    is an estimate, and an error of one codon is entirely possible.
    """
    if v_end is None and j_start is None:
        return None

    # Real anchors first; fall back to the positional estimate for whichever
    # side has none.
    c_real = find_cys(frames, v_end - 60, v_end + 15) if v_end is not None else []
    j_real = (find_j_motifs(frames, j_start - 30, j_start + 90, allow_weak=True)
              if j_start is not None else [])

    c_opts = [(p, True) for p in c_real]
    if not c_opts and v_end is not None:
        c_opts = [(v_end - EXPECTED_C_OFFSET, False)]

    j_opts = [(p, True, motif) for p, _tier, motif in j_real]
    if not j_opts and j_start is not None:
        j_opts = [(j_start + EXPECTED_J_OFFSET, False, "....")]

    if not c_opts or not j_opts:
        return None

    best = None

    for c_pos, c_is_real in c_opts:
        for j_pos, j_is_real, end_motif in j_opts:

            # An inferred anchor has no frame of its own, so snap it onto the
            # frame of the real one.  Work on copies -- mutating the loop
            # variables would leak the adjustment into the next iteration.
            p_c, p_j = c_pos, j_pos
            if not c_is_real and j_is_real:
                p_c -= (p_c - p_j) % 3
            elif not j_is_real:
                p_j -= (p_j - p_c) % 3

            if p_j <= p_c or (p_j - p_c) % 3 != 0:
                continue
            if p_c < 0 or p_j + 3 > frames.length:
                continue

            n_aa = (p_j - p_c) // 3 - 1
            if n_aa < min_len or n_aa > max_len:
                continue

            cdr_ntseq = frames.seq[p_c + 3:p_j]
            cdr_aa = str(Seq(cdr_ntseq).translate())

            # Prefer real anchors, then a length close to the locus norm.
            rank = (c_is_real + j_is_real, -abs(n_aa - 13))
            if best is None or rank > best["rank"]:
                best = {
                    "rank": rank,
                    "p_c": p_c,
                    "p_j": p_j,
                    "frame": p_c % 3,
                    "nt_start": p_c + 3,
                    "nt_end": p_j,
                    "cdr_ntseq": cdr_ntseq,
                    "cdr_aa": cdr_aa,
                    "junction_ntseq": frames.seq[p_c:p_j + 3],
                    "junction_aa": str(Seq(frames.seq[p_c:p_j + 3]).translate()),
                    "score": 0.0,
                    "start_motif": "".join(
                        x if x else "?" for x in
                        (frames.aa_at(p_c - 6), frames.aa_at(p_c - 3),
                         frames.aa_at(p_c) or "?")
                    ),
                    "end_motif": end_motif,
                    "search_method": "position_%s%s" % (
                        "motif" if c_is_real else "inferred",
                        "_motif" if j_is_real else "_inferred",
                    ),
                    "anchor_inferred": not (c_is_real and j_is_real),
                    "repair": "none",
                    "repair_confident": False,
                    "frame_status": (
                        "stop_codon" if "*" in cdr_aa else "productive"
                    ),
                    "productive": "*" not in cdr_aa and "X" not in cdr_aa,
                }

    if best is not None:
        best.pop("rank", None)
    return best


def qc_of(cand):
    """pass / stop_codon_present / low_qual.

    "pass" means the CDR3 was delimited by two intact anchors and read in a
    clean frame -- not merely that it is in frame.  A call resting on a
    degenerate anchor is low_qual even when it translates cleanly, because
    those are exactly the calls that turn out wrong on audit.
    """
    if cand["frame_status"] == "stop_codon" or "*" in cand["cdr_aa"]:
        return QC_STOP
    if cand.get("anchor_inferred") or not anchors_canonical(cand):
        return QC_LOW
    if cand["frame_status"] == "out_of_frame":
        return QC_PASS if cand.get("repair_confident") else QC_LOW
    if cand["frame_status"] == "productive":
        return QC_PASS
    return QC_LOW


# ============================================================
# Core search over one oriented sequence
# ============================================================

def search_oriented(frames, end_residue, expected_len, v_end, j_start, j_end,
                    min_len, max_len, allow_weak=False):
    """Best candidate for one orientation, or None.

    With no V alignment there is no positional anchor for Cys104, and a window
    derived from the J end alone is wide enough to catch an arbitrary upstream
    cysteine.  On the reference set that produced 34 identical 21-mers on
    germline J-only reads (v_gene = NA), which have no junction at all.  So
    when the V coordinate is missing, only a canonical Y-x-C paired with a
    canonical [WF]G-x-G is accepted.
    """
    best = None
    have_v = v_end is not None
    have_j = j_start is not None

    for pass_cfg in SEARCH_PASSES:

        requires = pass_cfg["requires"]
        if requires == "any" and not (have_v or have_j):
            continue
        if requires == "both" and not (have_v and have_j):
            continue

        # The non-productive rescue is a last resort -- but only against a
        # result whose anchors are actually canonical.  Skipping it merely
        # because *something* was found lets a junk call on degenerate anchors
        # block the correct frameshifted junction behind it.
        if (
            pass_cfg.get("relaxed", False)
            and best is not None
            and anchors_canonical(best)
        ):
            continue

        eff_v_end = v_end if pass_cfg["use_v"] else None
        eff_j_start = j_start if pass_cfg["use_j"] else None
        eff_j_end = j_end if pass_cfg["use_j"] else None

        # With no positional evidence at all, scan the whole read.
        if eff_v_end is None and eff_j_start is None:
            c_lo, c_hi = 0, frames.length
            j_lo, j_hi = 0, frames.length
        else:
            # Windows need an anchor to hang off; if the matching coordinate
            # is missing or suppressed, derive it from the other end.
            if eff_v_end is not None:
                c_lo = eff_v_end + pass_cfg["c_lo"]
                c_hi = eff_v_end + pass_cfg["c_hi"]
            else:
                c_lo = eff_j_start + pass_cfg["j_lo"] - 3 * (max_len + 1)
                c_hi = eff_j_start + pass_cfg["j_hi"] - 3 * min_len

            if eff_j_start is not None:
                j_lo = eff_j_start + pass_cfg["j_lo"]
                j_hi = eff_j_start + pass_cfg["j_hi"]
            else:
                j_lo = c_lo + 3 * (min_len + 1)
                j_hi = c_hi + 3 * (max_len + 1)

        require_canonical = eff_v_end is None
        relaxed = pass_cfg.get("relaxed", False)
        canonical_only = pass_cfg.get("canonical_only", False)
        pass_best = None

        cys = find_cys(frames, c_lo, c_hi)
        if not cys:
            continue

        motifs = find_j_motifs(frames, j_lo, j_hi, allow_weak=allow_weak)
        if not motifs:
            continue

        # If a canonical [WF]G-x-G is present in the window, it *is* the J
        # anchor -- do not let a degenerate one elsewhere stand in for it.
        # Without this, a read carrying a 1-nt indel (so the real FGQG is out
        # of frame with the Cys) gets a longer in-frame call that reads
        # straight through FR4, e.g. QQYYSTPGL -> QQYYSTPGLLARGPSWRSNI.  With
        # it, no in-frame pair survives and the junction is reported honestly
        # as out_of_frame.
        if any(tier == 2 for _, tier, _ in motifs):
            motifs = [m for m in motifs if m[1] == 2]

        for p_c in cys:
            for p_j, j_tier, end_motif in motifs:

                if p_j <= p_c:
                    continue

                nt_len = p_j - p_c - 3
                n_aa = nt_len // 3
                if n_aa < min_len or n_aa > max_len:
                    continue

                in_frame = nt_len % 3 == 0
                has_stop = in_frame and frames.has_stop_between(p_c + 3, p_j)

                if relaxed:
                    # Only a canonical Y-x-C / [WF]G-x-G pair earns a
                    # non-productive call.
                    if j_tier != 2 or frames.aa_at(p_c - 6) != "Y":
                        continue
                    if in_frame and not has_stop:
                        continue
                else:
                    # Same reading frame, and no stop inside the junction.
                    if not in_frame or has_stop:
                        continue

                if require_canonical and (
                    j_tier != 2 or frames.aa_at(p_c - 6) != "Y"
                ):
                    continue

                # Full FR3 / FR4 signature, for scans with no positional prior.
                if canonical_only:
                    t122 = frames.aa_at(p_j + 12)
                    if (
                        j_tier != 2
                        or frames.aa_at(p_c - 6) != "Y"
                        or not (t122 in ("T", "A")
                                or fr4_support(frames, p_j) >= FR4_MIN_SUPPORT)
                    ):
                        continue

                score, start_motif = score_candidate(
                    frames, p_c, p_j, j_tier, end_residue, expected_len,
                    eff_v_end, eff_j_start, eff_j_end,
                )

                nt_start = p_c + 3
                nt_end = p_j
                cdr_ntseq = frames.seq[nt_start:nt_end]
                junction_ntseq = frames.seq[p_c:p_j + 3]

                # An out-of-frame junction used to be truncated at the shift,
                # which discarded its 3' half.  Each half is now read in its
                # own frame instead; cdr_ntseq stays the exact read sequence.
                repair = "none"
                repair_confident = False

                if in_frame:
                    cdr_aa = str(Seq(cdr_ntseq).translate())
                else:
                    cdr_aa, repair, repair_confident = translate_frameshifted(
                        frames, p_c, p_j
                    )

                junction_aa = str(
                    Seq(junction_ntseq[: len(junction_ntseq) // 3 * 3]).translate()
                )

                if "X" in cdr_aa:
                    score -= 1.0

                if not in_frame:
                    frame_status = "out_of_frame"
                elif has_stop or "*" in cdr_aa:
                    frame_status = "stop_codon"
                elif "X" in cdr_aa:
                    frame_status = "ambiguous"
                else:
                    frame_status = "productive"

                cand = {
                    "p_c": p_c,
                    "p_j": p_j,
                    "frame": p_c % 3,
                    "nt_start": nt_start,
                    "nt_end": nt_end,
                    "cdr_ntseq": cdr_ntseq,
                    "cdr_aa": cdr_aa,
                    "junction_ntseq": junction_ntseq,
                    "junction_aa": junction_aa,
                    "frame_status": frame_status,
                    "productive": frame_status == "productive",
                    "score": score,
                    "start_motif": start_motif,
                    "end_motif": end_motif,
                    "search_method": pass_cfg["name"],
                    "anchor_inferred": False,
                    "repair": repair,
                    "repair_confident": repair_confident,
                }

                if pass_best is None or score > pass_best["score"] or (
                    score == pass_best["score"] and p_c < pass_best["p_c"]
                ):
                    pass_best = cand

        if pass_best is None:
            continue

        # Scores are not comparable across passes (the positional bonuses
        # differ), so rank on the signature first and only then on score.
        # Anchor quality outranks productivity: the anchors are what locate the
        # CDR3, whereas being in frame is a property of the biology.  A
        # productive-looking call on degenerate anchors is simply the wrong
        # location, and must not beat a real junction that carries an indel.
        if best is None or rank_key(pass_best) > rank_key(best):
            best = pass_best

        # A canonical call is as good as this gets; anything short of that is
        # worth letting a later pass try to beat -- which is how a read whose
        # V/J coordinates disagree with the read still gets the right answer.
        if is_canonical(best):
            return best

    return best


# ============================================================
# Confidence label
# ============================================================

def confidence_of(cand, v_end, j_start):
    """high / medium / low.

    "high" is reserved for calls carrying the complete FR3/FR4 signature.  A
    call that scores well on position alone but has lost Y-x-C or [WF]G-x-G is
    capped at "medium" however high its score -- otherwise a confident-looking
    label gets attached to exactly the calls most likely to be wrong.  Filter
    on confidence == "high" for a precision-oriented set.
    """
    if cand.get("anchor_inferred"):
        return "low"
    if not cand["productive"] or not is_canonical(cand):
        return "medium" if cand["score"] >= 6.0 else "low"

    both_coords = v_end is not None and j_start is not None

    if cand["search_method"] == "anchored" and both_coords:
        return "high"
    if cand["score"] >= 7.0:
        return "high"
    if cand["score"] >= 4.0:
        return "medium"
    return "low"


# ============================================================
# FASTA
# ============================================================

def read_fasta(path, wanted=None):
    """Load only the records named in `wanted` (all records if None).

    Returns id -> list of sequences.  Read IDs are not unique in practice
    (vdj_candidate.fa carries 57666 records under 57401 distinct IDs), and a
    plain id -> sequence dict silently keeps whichever copy came last.  When
    that copy is the shorter one, the row's alignment coordinates point past
    the end of it and the read becomes uncallable, so every copy is kept and
    the caller picks the one the coordinates actually fit.
    """
    seqs = {}
    seq_id = None
    buf = []

    def store():
        if seq_id is not None and (wanted is None or seq_id in wanted):
            seqs.setdefault(seq_id, []).append("".join(buf))

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                seq_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)

    store()
    return seqs


def select_record(records, coords, v_seq=None):
    """Pick the FASTA record a row's alignment coordinates refer to.

    Returns (index, sequence).  With a single copy there is nothing to decide;
    with several, prefer one long enough to contain the alignment, and among
    those one that actually contains the aligned V sequence.
    """
    if len(records) == 1:
        return 0, records[0]

    need = max([c for c in coords if c is not None], default=0)
    fits = [(i, r) for i, r in enumerate(records) if len(r) >= need]
    if not fits:
        fits = list(enumerate(records))

    # v_seq arrives as a float NaN when the cell is empty, so coerce rather
    # than assume a string (NaN is truthy, which a bare `if v_seq` would miss).
    probe = "" if v_seq is None else str(v_seq).strip().upper().replace("-", "")
    if probe and probe not in ("NA", "NAN", "NONE", "."):
        for i, r in fits:
            up = r.upper()
            if probe in up or probe in str(Seq(up).reverse_complement()):
                return i, r

    return max(fits, key=lambda ir: len(ir[1]))


# ============================================================
# Input adapters -- accept paths or in-memory objects
# ============================================================

REQUIRED_COLUMNS = ("Query",)

COORDINATE_COLUMNS = (
    "v_querystart", "v_queryend", "j_querystart", "j_queryend",
)


def _load_table(table, sep):
    """Accept a path, a file-like object, or a DataFrame."""
    if isinstance(table, pd.DataFrame):
        frame = table.copy()
    else:
        frame = pd.read_csv(table, sep=sep)

    frame.columns = frame.columns.str.strip()

    missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(
            "annotation table is missing required column(s): %s; got %s"
            % (", ".join(missing), ", ".join(map(str, frame.columns)))
        )

    return frame


def _load_fasta(fasta, wanted):
    """Accept a path or an already-loaded {id: seq} / {id: [seq, ...]} dict."""
    if isinstance(fasta, dict):
        loaded = {}
        for key, value in fasta.items():
            key = str(key)
            if wanted is not None and key not in wanted:
                continue
            loaded[key] = [value] if isinstance(value, str) else list(value)
        return loaded

    return read_fasta(fasta, wanted=wanted)


# ============================================================
# Main entry point
# ============================================================

def extract_cdr3(table=None, fasta=None, out_tsv=None, out_bed=None, sep="\t",
                 min_len=3, max_len=45, min_score=None, allow_weak_j=False,
                 verbose=False):
    """Extract CDR3 sequences and return them as a DataFrame.

    table   : path to the annotation TSV, or a pandas DataFrame, or None.
              Needs a Query column; v_gene/j_gene/strand and the four V/J
              coordinate columns are all optional -- each is used when present
              and worked around when absent.  None means "no annotation at
              all": every read in the FASTA is scanned on the canonical
              signature alone.
    fasta   : path to the reads, or {id: sequence} / {id: [sequence, ...]}.
    out_tsv : write the call table here (skipped when None).
    out_bed : write BED6 in query coordinates here (skipped when None).
    verbose : progress and per-category counts to stderr.  Off by default so
              importing this into a pipeline stays quiet.

    Returns the DataFrame; per-category counts are on df.attrs["stats"].
    Nothing is written unless out_tsv / out_bed are given.
    """
    if table is None:
        # No annotation at all: take every read in the FASTA and rely purely
        # on the canonical Y-x-C ... [WF]G-x-G signature.
        fasta_dict = _load_fasta(fasta, wanted=None)
        table = pd.DataFrame({"Query": sorted(fasta_dict)})
        queries = set(table["Query"])
    else:
        table = _load_table(table, sep)
        queries = set(table["Query"].astype(str))
        fasta_dict = _load_fasta(fasta, wanted=queries)

    if verbose and not any(c in table.columns for c in COORDINATE_COLUMNS):
        print(
            "[cdr3] note: no V/J coordinate columns found; falling back to a "
            "whole-read canonical scan (lower recall)",
            file=sys.stderr,
        )

    if verbose:
        print(
            f"[cdr3] {len(table)} rows, "
            f"{len(fasta_dict)}/{len(queries)} sequences found in FASTA"
            f" ({sum(len(v) for v in fasta_dict.values()) - len(fasta_dict)}"
            f" duplicate records)",
            file=sys.stderr,
        )

    cdr3_results = []
    bed_data = []
    stats = Counter()

    # Orientation is a property of (query, strand); cache the translation.
    frame_cache = {}

    def get_frames(query, seq, reverse, rec_idx=0):
        key = (query, rec_idx, reverse)
        if key not in frame_cache:
            oriented = str(Seq(seq).reverse_complement()) if reverse else seq
            frame_cache[key] = Frames(oriented)
        return frame_cache[key]

    for _, row in table.iterrows():

        query = str(row["Query"])
        v_gene = str(row.get("v_gene", "NA"))
        j_gene = str(row.get("j_gene", "NA"))
        strand = str(row.get("strand", "NA")).strip()

        if query not in fasta_dict:
            stats["no_fasta"] += 1
            if verbose:
                print(f"WARNING: {query} not found in FASTA", file=sys.stderr)
            continue

        # --- plus-strand alignment coordinates (1-based inclusive) --------
        v_qs, v_qe = to_int(row.get("v_querystart")), to_int(row.get("v_queryend"))
        j_qs, j_qe = to_int(row.get("j_querystart")), to_int(row.get("j_queryend"))

        records = fasta_dict[query]
        rec_idx, nt_seq = select_record(
            records, (v_qs, v_qe, j_qs, j_qe), row.get("v_seq"),
        )
        nt_seq = nt_seq.upper()
        seq_len = len(nt_seq)

        if len(records) > 1:
            stats["duplicate_read_id"] += 1

        # --- locus --------------------------------------------------------
        # The locus only picks which of W / F to prefer at position 118 and
        # sets the weak length prior -- neither is essential, and find_j_motifs
        # accepts both aromatics regardless.  So an unidentified V/J is no
        # reason to drop the read: search with both anchors and let the
        # canonical-signature requirement do the filtering.
        locus = locus_of(v_gene) or locus_of(j_gene)

        if locus is None:
            stats["locus_unknown"] += 1
            end_residue = None            # no basis to prefer W over F
            expected_len = DEFAULT_EXPECTED_LEN
        else:
            end_residue = LOCUS_END_RESIDUE[locus]
            expected_len = LOCUS_EXPECTED_LEN.get(locus, DEFAULT_EXPECTED_LEN)

        # --- which orientations to try -------------------------------------
        if strand == "Plus/Plus":
            orientations = [False]
        elif strand == "Plus/Minus":
            orientations = [True]
        else:
            # v4 dropped these rows; the J coordinates are still usable, so
            # score both orientations and keep the better one.
            orientations = [False, True]

        best = None
        best_meta = None

        for reverse in orientations:

            frames = get_frames(query, nt_seq, reverse, rec_idx)

            # Map plus-strand [start-1, end) onto the oriented sequence.
            if reverse:
                v_end = seq_len - v_qs + 1 if v_qs is not None else None
                j_start = seq_len - j_qe if j_qe is not None else None
                j_end = seq_len - j_qs + 1 if j_qs is not None else None
            else:
                v_end = v_qe if v_qe is not None else None
                j_start = j_qs - 1 if j_qs is not None else None
                j_end = j_qe if j_qe is not None else None

            cand = search_oriented(
                frames, end_residue, expected_len,
                v_end, j_start, j_end, min_len, max_len,
                allow_weak=allow_weak_j,
            )

            if cand is None:
                # Motif search found nothing: fall back to placing 105-117 by
                # IMGT offset from whichever alignment coordinate we have.
                cand = rescue_by_position(
                    frames, end_residue, v_end, j_start, j_end,
                    min_len, max_len,
                )

            if cand is None:
                continue

            if best is None or cand["score"] > best["score"]:
                best = cand
                best_meta = (reverse, frames, v_end, j_start)

        if best is None:
            stats["not_found"] += 1
            if verbose:
                print(f"WARNING: no CDR3 found for {query}", file=sys.stderr)
            continue

        if min_score is not None and best["score"] < min_score:
            stats["below_min_score"] += 1
            continue

        reverse, frames, v_end, j_start = best_meta
        confidence = confidence_of(best, v_end, j_start)
        qc = qc_of(best)

        stats["found"] += 1
        stats[f"qc_{qc}"] += 1
        stats[f"method_{best['search_method']}"] += 1
        stats[f"conf_{confidence}"] += 1
        stats[f"status_{best['frame_status']}"] += 1

        # --- back to query coordinates for the BED -------------------------
        if reverse:
            q_start = seq_len - best["nt_end"]
            q_end = seq_len - best["nt_start"]
            bed_strand = "-"
        else:
            q_start = best["nt_start"]
            q_end = best["nt_end"]
            bed_strand = "+"

        resolved_strand = strand
        if strand not in ("Plus/Plus", "Plus/Minus"):
            resolved_strand = "Plus/Minus" if reverse else "Plus/Plus"

        cdr3_results.append({
            "Query": query,
            "cdr_ntseq": best["cdr_ntseq"],
            "cdr_aa": best["cdr_aa"],
            "strand": resolved_strand,
            "frame": best["frame"],
            "search_method": best["search_method"],
            "locus": locus or "unknown",
            "v_gene": v_gene,
            "j_gene": j_gene,
            "cdr3_len": len(best["cdr_ntseq"]) // 3,
            "frame_status": best["frame_status"],
            "productive": best["productive"],
            "qc": qc,
            "anchor_inferred": best.get("anchor_inferred", False),
            "repair": best.get("repair", "none"),
            "junction_aa": best["junction_aa"],
            "junction_ntseq": best["junction_ntseq"],
            "score": round(best["score"], 3),
            "confidence": confidence,
            "start_motif": best["start_motif"],
            "end_motif": best["end_motif"],
            "cdr_start": q_start,
            "cdr_end": q_end,
        })

        bed_data.append({
            "chrom": query,
            "start": q_start,
            "end": q_end,
            "name": best["cdr_aa"],
            "score": int(max(0, min(1000, round(best["score"] * 50)))),
            "strand": bed_strand,
        })

    cdr3_df = pd.DataFrame(cdr3_results, columns=[
        "Query", "cdr_ntseq", "cdr_aa", "strand", "frame", "search_method",
        "locus", "v_gene", "j_gene", "cdr3_len", "frame_status", "productive",
        "qc", "anchor_inferred", "repair",
        "junction_aa", "junction_ntseq", "score", "confidence",
        "start_motif", "end_motif", "cdr_start", "cdr_end",
    ])

    bed_df = pd.DataFrame(bed_data, columns=[
        "chrom", "start", "end", "name", "score", "strand",
    ])

    cdr3_df.attrs["stats"] = dict(stats)

    if out_tsv is not None:
        cdr3_df.to_csv(out_tsv, sep="\t", index=False)

    if out_bed is not None:
        bed_df.to_csv(out_bed, sep="\t", index=False, header=False)

    if verbose:
        total = len(table)
        found = stats["found"]
        print(
            f"[cdr3] CDR3 called for {found}/{total} rows "
            f"({100.0 * found / total:.1f}%)",
            file=sys.stderr,
        )
        for key in sorted(stats):
            if key != "found":
                print(f"[cdr3]   {key}: {stats[key]}", file=sys.stderr)

    return cdr3_df


def cdr3_extract(df, fasta_file, output_bed=None, output_cdr3=None, **kwargs):
    """Backward-compatible alias for v4/v5 callers."""
    return extract_cdr3(
        table=df, fasta=fasta_file,
        out_tsv=output_cdr3, out_bed=output_bed,
        **kwargs
    )


def main():
    parser = argparse.ArgumentParser(
        description="Extract CDR3 sequences using V/J alignment anchors.",
    )
    parser.add_argument("--table", "-t", default=None,
                        help="TSV with Query and, optionally, v_gene/j_gene/"
                             "strand and v_querystart/v_queryend/j_querystart/"
                             "j_queryend.  Omit to scan every read in the FASTA "
                             "with no annotation at all.")
    parser.add_argument("--fasta", "-f", required=True)
    parser.add_argument("--out-cdr3", "-o", default=None)
    parser.add_argument("--out-bed", "-b", default=None)
    parser.add_argument("--min-len", type=int, default=3,
                        help="minimum CDR3-IMGT length in aa (default 3)")
    parser.add_argument("--max-len", type=int, default=45,
                        help="maximum CDR3-IMGT length in aa (default 45)")
    parser.add_argument("--min-score", type=float, default=None,
                        help="drop calls scoring below this (default: keep all; "
                             "filter on the score/confidence columns instead)")
    parser.add_argument("--allow-weak-j", action="store_true",
                        help="also accept J anchors with no glycine in [WF]G.G "
                             "(raises recall, lowers precision)")
    parser.add_argument("--quiet", action="store_true")

    args = parser.parse_args()

    extract_cdr3(
        table=args.table,
        fasta=args.fasta,
        out_bed=args.out_bed,
        out_tsv=args.out_cdr3,
        min_len=args.min_len,
        max_len=args.max_len,
        min_score=args.min_score,
        allow_weak_j=args.allow_weak_j,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
