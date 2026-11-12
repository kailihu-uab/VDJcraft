#!/usr/bin/env python3

import argparse
import collections
import glob
import itertools
import os
import re

from error_correct import align, signature, TRANSITIONS, is_hotspot

MIN_CLONE_READS = 3
MIN_COVERED_NT = 150          # too little V aligned to quote a rate from
BCR_V = ("IGHV", "IGKV", "IGLV")   # SHM is scored for BCR clonotypes only, not TCR

_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[i] for i, (a, b, c) in enumerate(itertools.product(_B, repeat=3))}

# IMGT unique numbering: codon ranges are fixed across all V genes.
REGIONS = [("FR1", 1, 26), ("CDR1", 27, 38), ("FR2", 39, 55),
           ("CDR2", 56, 65), ("FR3", 66, 104)]
REGION_NAMES = [r[0] for r in REGIONS]
CDR_REGIONS = ("CDR1", "CDR2")
FR_REGIONS = ("FR1", "FR2", "FR3")


def region_of(codon):
    for name, lo, hi in REGIONS:
        if lo <= codon <= hi:
            return name
    return "CDR3"


def load_germline(path, patterns=("IG*V*.fasta*", "IG*V*.fa", "IG*V*.fasta")):
 
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
    """clone_correct_v1 writes v_seq_corrected alongside v_seq_original; older files
    carried the corrected sequence in v_seq itself."""
    if rows and "v_seq_corrected" in rows[0]:
        return "v_seq_corrected"
    return "v_seq"


def analyse_clone(reads, ref, is_gapped, min_covered_nt=MIN_COVERED_NT, vcol="v_seq"):
    """Mutation profile of one clone's V consensus against its germline allele."""
    seqs = collections.Counter(r[vcol] for r in reads if r.get(vcol, ".") not in (".", ""))
    if not seqs:
        return None
    allele = collections.Counter(r["v_gene"] for r in reads).most_common(1)[0][0]
    germ_name, gapped = germline_for(ref, allele)
    if gapped is None:
        return {"skip": "no germline for this allele"}
    germ = gapped.replace(".", "").upper()
    cmap = imgt_codon_map(gapped) if is_gapped else None
    cons = seqs.most_common(1)[0][0]

    al = align(germ, cons)
    if al is None:
        return {"skip": "no alignment to germline"}
    s, e = al["span"]
    covered = e - s
    if covered < min_covered_nt:
        return {"skip": f"under {min_covered_nt}nt of V aligned"}

    def reg_at(pos):
        if cmap is None:
            return "V"
        return region_of(cmap[pos]) if pos < len(cmap) else "CDR3"

    cov_by_region = collections.Counter()
    for i in range(s, e):
        cov_by_region[reg_at(i)] += 1

    subs, ins, dels = [], 0, 0
    for v in signature(germ, al):
        if v[0] == "sub":
            subs.append(v)
        elif v[0] == "ins":
            ins += len(v[2])
        else:
            dels += v[2]

    per_region = collections.Counter()
    rs = collections.Counter()
    rs_by_region = collections.defaultdict(collections.Counter)
    ti = tv = hot = 0
    for _, pos, gb, nb in subs:
        reg = reg_at(pos)
        per_region[reg] += 1
        if (gb, nb) in TRANSITIONS:
            ti += 1
        else:
            tv += 1
        if is_hotspot(germ, pos):
            hot += 1
        k = classify_sub(germ, pos, nb)
        if k:
            rs[k] += 1
            rs_by_region[reg][k] += 1

    isotype = collections.Counter(
        r["c_gene"].split("_")[0] for r in reads if r["c_gene"] not in (".", "")
    ).most_common(1)
    return {"allele": allele, "germline": germ_name, "n_reads": len(reads), "covered": covered,
            "n_sub": len(subs), "n_ins": ins, "n_del": dels,
            "freq": len(subs) / covered,
            "per_region": per_region, "cov_by_region": cov_by_region,
            "rs": rs, "rs_by_region": rs_by_region, "ti": ti, "tv": tv, "hot": hot,
            "isotype": isotype[0][0] if isotype else ".",
            "j_gene": collections.Counter(r["j_gene"] for r in reads).most_common(1)[0][0]}


def _bar(frac, width=30):
    return "#" * int(round(frac * width)) + "." * (width - int(round(frac * width)))


def run_shm_analysis(corrected_tsv, germline, out_report=None, out_tsv=None,
                     min_clone_reads=MIN_CLONE_READS, min_covered_nt=MIN_COVERED_NT,
                     top_n=15, verbose=True):

    prefix = corrected_tsv.rsplit(".", 1)[0]
    out_report = out_report or f"{prefix}.shm_report.txt"
    out_tsv = out_tsv or f"{prefix}.shm_per_clone.tsv"
    say = print if verbose else (lambda *a, **k: None)

    ref, is_gapped = load_germline(germline)
    rows = load_rows(corrected_tsv)
    vcol = pick_v_column(rows)
    by_clone = collections.defaultdict(list)
    # BCR only (heavy and light chain): TCR reads are left out of every output
    n_other = sum(1 for r in rows if not r.get("v_gene", "").startswith(BCR_V))
    rows = [r for r in rows if r.get("v_gene", "").startswith(BCR_V)]
    for r in rows:
        # error_correct_new names the column clonotype; older tables call it clone
        cid = r.get("clonotype", r.get("clone", "."))
        if cid != ".":
            by_clone[cid].append(r)

    results, skipped = {}, collections.Counter()
    for clone, reads in by_clone.items():
        if len(reads) < min_clone_reads:
            skipped[f"clonotype under {min_clone_reads} reads"] += 1
            continue
        res = analyse_clone(reads, ref, is_gapped, min_covered_nt, vcol)
        if res is None:
            skipped["no V sequence"] += 1
        elif "skip" in res:
            skipped[res["skip"]] += 1
        else:
            results[clone] = res
    if not results:
        raise SystemExit("no BCR clonotypes analysable -- check the germline reference matches "
                         "the v_gene naming in the input")

    names = REGION_NAMES if is_gapped else ["V"]
    order = sorted(results, key=lambda c: -results[c]["n_reads"])
    tot_reads = sum(results[c]["n_reads"] for c in order)
    freqs = [results[c]["freq"] for c in order]
    wfreq = sum(results[c]["freq"] * results[c]["n_reads"] for c in order) / tot_reads
    agg = lambda key, reg: sum(results[c][key][reg] for c in order)

    with open(out_tsv, "w") as f:
        cols = (["clonotype", "n_reads", "v_allele", "j_gene", "isotype", "covered_nt",
                 "n_sub", "n_ins", "n_del", "mut_freq_pct", "ti", "tv", "R", "S",
                 "hotspot_frac"] + [f"{r}_mut" for r in names] + [f"{r}_cov" for r in names])
        f.write("\t".join(cols) + "\n")
        for c in order:
            d = results[c]
            f.write("\t".join(str(x) for x in [
                c, d["n_reads"], d["allele"], d["j_gene"], d["isotype"], d["covered"],
                d["n_sub"], d["n_ins"], d["n_del"], f"{100*d['freq']:.2f}",
                d["ti"], d["tv"], d["rs"]["R"], d["rs"]["S"],
                f"{d['hot']/d['n_sub']:.3f}" if d["n_sub"] else "."]
                + [d["per_region"][r] for r in names]
                + [d["cov_by_region"][r] for r in names]) + "\n")

    L = []
    A = L.append
    A("SOMATIC HYPERMUTATION REPORT")
    A("=" * 74)
    A(f"input            : {corrected_tsv}")
    A(f"germline         : {germline}")
    A(f"                   {len(ref)} alleles, "
      f"{'IMGT-gapped (region analysis available)' if is_gapped else 'UNGAPPED (no region analysis)'}")
    A(f"V sequence column: {vcol}")
    via = collections.Counter(f"{d['allele']} -> {d['germline']}" for d in results.values()
                              if d["germline"] != d["allele"])
    if via:
        A(f"paralog germline : {sum(via.values())} clonotypes scored against the "
          f"D-paralog of their allele (" + ", ".join(f"{k} x{v}" for k, v in via.most_common()) + ")")
    A(f"chains           : BCR only (IGH, IGK, IGL); {n_other} TCR/unassigned "
      f"reads excluded")
    A(f"unit of analysis : clonotype consensus (>= {min_clone_reads} reads, "
      f">= {min_covered_nt}nt V aligned)")
    A("")
    A(f"clonotypes used  : {len(results)}   representing {tot_reads} reads")
    if skipped:
        A("clonotypes skipped: " + ", ".join(f"{v} ({k})" for k, v in skipped.most_common()))
    A("")
    A("-" * 74)
    A("1. MUTATION LOAD")
    A("-" * 74)
    srt = sorted(freqs)
    A(f"  per clonotype: mean {100*sum(freqs)/len(freqs):.2f}%   "
      f"median {100*srt[len(srt)//2]:.2f}%   range {100*min(freqs):.2f}-{100*max(freqs):.2f}%")
    A(f"  read-weighted: {100*wfreq:.2f}%   (weights by clonal expansion)")
    A("")
    A("  distribution of clonotypes by V mutation frequency:")
    for lo, hi, lab in [(0, 1, "0-1%   (near-germline / naive)"), (1, 2, "1-2%"),
                        (2, 5, "2-5%   (moderate)"), (5, 10, "5-10%  (high)"),
                        (10, 1e9, ">10%   (very high)")]:
        n = sum(1 for x in freqs if lo <= 100 * x < hi)
        A(f"    {lab:<30} {n:>4} clonotypes  {_bar(n/len(freqs))}")
    A("")
    A("-" * 74)
    A("2. REGION TARGETING  (mutations per 100 aligned nt)")
    A("-" * 74)
    if not is_gapped:
        A("  Skipped: the germline reference is not IMGT-gapped, so standard numbering")
        A("  is unavailable. Supply a gapped IMGT V reference (sequences padded with")
        A("  '.') to place each mutation in its FR/CDR region.")
        A("")
        A(f"  whole V: {agg('per_region','V')} mutations over "
          f"{agg('cov_by_region','V')} aligned nt")
    else:
        A("  SHM under antigen selection concentrates in CDRs and spares FRs.")
        A("")
        A(f"  {'region':<8}{'mutations':>11}{'covered nt':>12}{'per 100nt':>11}")
        for reg in names:
            m, cv = agg("per_region", reg), agg("cov_by_region", reg)
            A(f"  {reg:<8}{m:>11}{cv:>12}{(100*m/cv if cv else 0):>11.2f}")
        cdr = sum(agg("per_region", r) for r in CDR_REGIONS)
        cdrc = sum(agg("cov_by_region", r) for r in CDR_REGIONS)
        fr = sum(agg("per_region", r) for r in FR_REGIONS)
        frc = sum(agg("cov_by_region", r) for r in FR_REGIONS)
        dc = 100 * cdr / cdrc if cdrc else 0
        df = 100 * fr / frc if frc else 0
        A("")
        A(f"  CDR density {dc:.2f} vs FR density {df:.2f}  ->  CDR/FR ratio {(dc/df if df else 0):.2f}")
        A(f"  {'consistent with antigen-driven selection' if dc > df else 'no CDR enrichment -- selection not evident'}")
    A("")
    A("-" * 74)
    A("3. SELECTION  (replacement / silent)")
    A("-" * 74)
    A("  R/S above the ~2.9 neutral expectation in CDRs indicates positive selection;")
    A("  R/S below it in FRs indicates purifying selection preserving the fold.")
    A("")
    A(f"  {'region':<8}{'R':>8}{'S':>8}{'R/S':>9}")
    for reg in names:
        R = sum(results[c]["rs_by_region"][reg]["R"] for c in order)
        S = sum(results[c]["rs_by_region"][reg]["S"] for c in order)
        A(f"  {reg:<8}{R:>8}{S:>8}{(R/S if S else float('inf')):>9.2f}")
    R = sum(results[c]["rs"]["R"] for c in order)
    S = sum(results[c]["rs"]["S"] for c in order)
    A(f"  {'TOTAL':<8}{R:>8}{S:>8}{(R/S if S else float('inf')):>9.2f}")
    A("")
    A("-" * 74)
    A("4. MUTATION SPECTRUM")
    A("-" * 74)
    ti = sum(results[c]["ti"] for c in order)
    tv = sum(results[c]["tv"] for c in order)
    hot = sum(results[c]["hot"] for c in order)
    nsub = sum(results[c]["n_sub"] for c in order)
    A(f"  transitions {ti}  transversions {tv}   Ti/Tv = {(ti/tv if tv else 0):.2f}"
      f"   (SHM is transition-biased, ~1.5-2.5)")
    A(f"  in AID WRC/GYW hotspots: {hot}/{nsub} = {(100*hot/nsub if nsub else 0):.0f}%")
    A(f"  indels: {sum(results[c]['n_ins'] for c in order)} inserted / "
      f"{sum(results[c]['n_del'] for c in order)} deleted nt "
      f"(SHM indels are rare; a high count suggests residual error)")
    A("")
    A("-" * 74)
    A("5. MUTATION LOAD BY ISOTYPE")
    A("-" * 74)
    A("  Class-switched cells have undergone more germinal-centre cycles than IgM.")
    A("")
    iso = collections.defaultdict(list)
    for c in order:
        iso[results[c]["isotype"]].append(results[c]["freq"])
    A(f"  {'isotype':<10}{'clonotypes':>11}{'mean mut%':>12}{'median':>10}")
    for k in sorted(iso, key=lambda x: -len(iso[x])):
        v = sorted(iso[k])
        A(f"  {k:<10}{len(v):>11}{100*sum(v)/len(v):>12.2f}{100*v[len(v)//2]:>10.2f}")
    A("")
    A("-" * 74)
    A(f"6. MOST MUTATED CLONOTYPES (top {top_n})")
    A("-" * 74)
    A(f"  {'clonotype':<11}{'reads':>7}  {'V allele':<16}{'iso':<8}{'mut%':>7}{'R/S':>7}")
    for c in sorted(order, key=lambda x: -results[x]["freq"])[:top_n]:
        d = results[c]
        rr, ss = d["rs"]["R"], d["rs"]["S"]
        A(f"  {c:<11}{d['n_reads']:>7}  {d['allele']:<16}{d['isotype']:<8}"
          f"{100*d['freq']:>7.2f}{(rr/ss if ss else float('inf')):>7.2f}")
    A("")
    A("-" * 74)
    A("CAVEATS")
    A("-" * 74)
    A("  - Mutation counts depend on the germline allele call being right. A wrong")
    A("    allele inflates mutations; check clonotypes whose load is far above the rest.")
    A("  - Only the V region is scored. CDR3 is excluded: it spans the junction, so")
    A("    differences there are rearrangement, not hypermutation.")
    A("  - Rates use aligned positions as the denominator, so partially covered V")
    A("    regions are not penalised; per-region coverage is listed above.")
    A("  - Residual sequencing error inflates apparent SHM. Clones near the read")
    A("    floor rest on the thinnest evidence.")
    A("")
    A(f"per-clonotype table written to {out_tsv}")

    open(out_report, "w").write("\n".join(L) + "\n")
    say("\n".join(L[:22]))
    say(f"\n... full report -> {out_report}")
    return {"n_clones": len(results), "n_reads": tot_reads, "is_gapped": is_gapped,
            "mean_mut_freq": sum(freqs) / len(freqs), "read_weighted_mut_freq": wfreq,
            "skipped": dict(skipped), "out_report": out_report, "out_tsv": out_tsv,
            "per_clone": results}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("corrected_tsv", help="corrected table from clone_correct_v1")
    p.add_argument("-g", "--germline", required=True,
                   help="germline V reference: a FASTA file, or a directory of them. "
                        "IMGT-gapped sequences enable the FR/CDR region analysis.")
    p.add_argument("--out-report"); p.add_argument("--out-tsv")
    p.add_argument("--min-clone-reads", type=int, default=MIN_CLONE_READS)
    p.add_argument("--min-covered-nt", type=int, default=MIN_COVERED_NT)
    p.add_argument("--top-n", type=int, default=15)
    p.add_argument("-q", "--quiet", action="store_true")
    a = p.parse_args()
    run_shm_analysis(a.corrected_tsv, a.germline, out_report=a.out_report,
                     out_tsv=a.out_tsv, min_clone_reads=a.min_clone_reads,
                     min_covered_nt=a.min_covered_nt, top_n=a.top_n, verbose=not a.quiet)


if __name__ == "__main__":
    main()
