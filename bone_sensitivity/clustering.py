"""
Supplementary Figure S5 (Leiden domains) under the second review (2026-11-02).

1. The 'unresolved' share per section and region (all bins, region of interest, outside it) and its overlap with the
   cross-captured A549 lesion: the pooled share (manuscript) mixes sections with very different shares (F20, F2).
2. One row per final Leiden domain: its manuscript name (leiden_signature_matrix: the most enriched program, which
   favours rare programs, F21) next to its section composition, depth, most frequent programs and final markers, with
   flags for domains without a significant marker or defined by one gene (F22).
3. Leiden without Harmony: Harmony corrects on the section label although each section is its own condition (F19).
   The manuscript graph and Leiden call are rerun on the Harmony and on the plain PCA embedding and the two raw
   partitions (before the unresolved and fragment merging) are compared, with their section structure.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import pipeline_config as pc
from bone_sensitivity import common

UNRESOLVED_CSV = "sensitivity_FigureS5_unresolved_by_section.csv"
ANNOTATION_CSV = "sensitivity_leiden_domain_annotation.csv"
HARMONY_CSV = "sensitivity_leiden_harmony.csv"
STABILITY_CSV = "leiden_seed_stability_adjusted_rand_index.csv"      # manuscript seed-to-seed ARI (reference)
UNRESOLVED, UNASSIGNED, AMBIGUOUS = "unresolved", "unassigned", "ambiguous"
TOP_PROGRAMS = 3                  # programs listed per domain
ONE_GENE_MAX_MARKERS = 2          # "defined by one gene": at most this many significant markers,
ONE_GENE_MIN_DOMAIN_PCT = 90.0    # the top marker detected in at least this % of the domain's bins
ONE_GENE_MAX_REST_PCT = 5.0       # and in at most this % of the other bins
SECTION_SPECIFIC = 0.9            # a cluster is section-specific when at least this share of its bins is one section
NO_MARKER, ONE_GENE = "no significant marker", "defined by one gene"


def _pct(k, n, digits=2):
    return round(100.0 * float(k) / float(n), digits) if n else np.nan


def domain_order(labels):
    """Domain labels in reading order: numbered domains ascending, then the other labels, unresolved and unassigned
    last."""
    labels = [str(x) for x in pd.unique(pd.Series(labels).astype(str))]
    num = sorted((x for x in labels if x.isdigit()), key=int)
    tail = [x for x in (UNRESOLVED, UNASSIGNED) if x in labels]
    return num + sorted(x for x in labels if not x.isdigit() and x not in tail) + tail


# ---- 1. unresolved share per section and region ---------------------------------------------------------------

def unresolved_by_section(obs, samples, lesion=None, roi_key="in_marrow_roi"):
    """Per section (and "all") and region (all bins; inside / outside the region of interest when roi_key exists):
    bins, unresolved bins and their share, unassigned bins. With a lesion flag per bin, also the lesion bins, the share
    of them that is unresolved and the share of unresolved bins that is lesion."""
    leiden = obs["leiden"].astype(str).to_numpy()
    smp = obs["sample"].astype(str).to_numpy()
    roi = common.as_bool(obs[roi_key]) if roi_key in obs else None
    lesion = None if lesion is None else np.asarray(lesion, dtype=bool)
    un_all, ua_all = leiden == UNRESOLVED, leiden == UNASSIGNED
    rows = []
    for s in list(samples) + ["all"]:
        ms = np.ones(len(leiden), dtype=bool) if s == "all" else smp == s
        regions = [("all", ms)] + ([("roi", ms & roi), ("outside_roi", ms & ~roi)] if roi is not None else [])
        for region, m in regions:
            un = m & un_all
            row = {"sample": s, "region": region, "n_bins": int(m.sum()), "n_unresolved": int(un.sum()),
                   "pct_unresolved": _pct(un.sum(), m.sum()), "n_unassigned": int((m & ua_all).sum())}
            if lesion is not None:
                ml = m & lesion
                row.update({"n_lesion": int(ml.sum()), "pct_lesion_unresolved": _pct((ml & un).sum(), ml.sum()),
                            "pct_unresolved_lesion": _pct((ml & un).sum(), un.sum())})
            rows.append(row)
    cols = ["sample", "region", "n_bins", "n_unresolved", "pct_unresolved", "n_unassigned"]
    if lesion is not None:
        cols += ["n_lesion", "pct_lesion_unresolved", "pct_unresolved_lesion"]
    return pd.DataFrame(rows, columns=cols)


def lesion_flags(a, cfg, progress=None):
    """(in-lesion flag per bin or None, where it came from): obs in_lesion when present, otherwise the shared lesion
    definition (common.lesion_signal with sensitivity_analyses.lesion), computed here and not stored."""
    if "in_lesion" in a.obs:
        return common.as_bool(a.obs["in_lesion"]), "obs in_lesion"
    les = common.settings(cfg, None)["lesion"]
    _, flag, found = common.lesion_signal(a, tuple(les["genes"]), int(les["pool_k"]), les["min_pooled_umi"])
    if not found:
        common.warn(progress, f"Leiden domains: none of the lesion genes ({', '.join(les['genes'])}) is in the data; "
                              f"no lesion overlap")
        return None, ""
    return flag, (f"{'+'.join(found)} UMIs of each bin and its {int(les['pool_k'])} nearest bins >= "
                  f"{les['min_pooled_umi']:g}")


# ---- 2. domain annotation -------------------------------------------------------------------------------------

def manuscript_names(a, cfg):
    """{domain: name} exactly as Supplementary Figure S5 names the domains (bone_analysis.leiden_signature_matrix,
    called as in _run; it only reads obs)."""
    import bone_analysis as ba
    M, _ = ba.leiden_signature_matrix(a, cfg["primary"], cfg["signatures"], lambda c: pc.signature_label(cfg, c))
    order = [str(x) for x in a.obs.groupby("leiden", observed=True).size().index]
    if len(order) != len(M):
        return {d: d for d in order}
    return dict(zip(order, M["leiden_domain"].astype(str)))


def marker_table(rgg, key="leiden"):
    """uns["rank_genes_groups"] as a long table (group, names, scores, pvals_adj), or None when it is absent or was
    computed on another grouping."""
    if not rgg or "names" not in rgg or (rgg.get("params") or {}).get("groupby") != key:
        return None
    frames = []
    for g in rgg["names"].dtype.names:
        frames.append(pd.DataFrame({"group": str(g), "names": np.asarray(rgg["names"][g]).astype(str),
                                    "scores": np.asarray(rgg["scores"][g], dtype=float),
                                    "pvals_adj": np.asarray(rgg["pvals_adj"][g], dtype=float)}))
    return pd.concat(frames, ignore_index=True) if frames else None


def summarize_markers(markers, marker_padj=0.05):
    """Per group: significant markers (BH-adjusted p < marker_padj and positive score) and the top marker (highest
    score among them; empty when there is none: the best-scoring gene of a domain without a significant marker can be
    depleted, F22)."""
    cols = ["n_significant_markers", "top_marker"]
    if markers is None or not len(markers):
        return pd.DataFrame(columns=cols)
    m = markers.sort_values(["group", "scores"], ascending=[True, False], kind="stable")
    sig = m[(m["pvals_adj"] < float(marker_padj)) & (m["scores"] > 0)]
    groups = pd.unique(m["group"])
    top_sig = sig.groupby("group", sort=False)["names"].first()
    out = pd.DataFrame(index=pd.Index(groups, name="group"))
    out["n_significant_markers"] = sig.groupby("group").size().reindex(groups).fillna(0).astype(int).to_numpy()
    out["top_marker"] = top_sig.reindex(groups).fillna("").astype(str).to_numpy()
    return out


def flag_domains(n_significant, pct_domain, pct_rest):
    """'no significant marker' (none), 'defined by one gene' (<= 2 significant markers, the top marker detected in
    >= 90% of the domain's bins and <= 5% of the others) or ''."""
    n = np.asarray(n_significant, dtype=float)
    d, r = np.asarray(pct_domain, dtype=float), np.asarray(pct_rest, dtype=float)
    one = (n <= ONE_GENE_MAX_MARKERS) & (d >= ONE_GENE_MIN_DOMAIN_PCT) & (r <= ONE_GENE_MAX_REST_PCT)
    return np.where(n == 0, NO_MARKER, np.where(one, ONE_GENE, "")).astype(object)


def top_programs(programs, n_top=TOP_PROGRAMS, ambiguous=AMBIGUOUS):
    """'osteoblast 23%; osteoclast 20%; skeletal_muscle 18%': the most frequent assigned programs with their share of
    the non-ambiguous bins."""
    p = pd.Series(programs).astype(str)
    p = p[p != ambiguous]
    if not len(p):
        return ""
    vc = p.value_counts()
    vc = vc.reset_index().set_axis(["program", "n"], axis=1).sort_values(["n", "program"], ascending=[False, True])
    return "; ".join(f"{r.program} {100 * r.n / len(p):.0f}%" for r in vc.head(n_top).itertuples())


def _detection(a, labels, genes_by_domain):
    """(% of the domain's bins, % of the other bins) detecting each domain's gene (raw counts > 0)."""
    genes = sorted({g for g in genes_by_domain.values() if g and g in a.var_names})
    out = {d: (np.nan, np.nan) for d in genes_by_domain}
    if not genes:
        return out
    X = (a.layers["counts"] if "counts" in a.layers else a.X)[:, a.var_names.get_indexer(genes)]
    det = (X.toarray() if hasattr(X, "toarray") else np.asarray(X)) > 0
    col = {g: j for j, g in enumerate(genes)}
    for d, g in genes_by_domain.items():
        if g in col:
            inside = labels == d
            v = det[:, col[g]]
            out[d] = (_pct(v[inside].sum(), inside.sum()), _pct(v[~inside].sum(), (~inside).sum()))
    return out


def domain_annotation(a, cfg, marker_padj=0.05, progress=None):
    """One row per final Leiden domain: manuscript name, size, section composition, median UMIs and genes per UMI,
    top programs and ambiguous share, significant final markers, the top marker's detection in and outside the
    domain, flags and a suggested label."""
    obs = a.obs
    ld = obs["leiden"].astype(str).to_numpy()
    domains = domain_order(ld)
    samples = common.sample_ids(a, cfg)
    names = manuscript_names(a, cfg)
    comp = pd.crosstab(pd.Series(ld, name="d"), pd.Series(obs["sample"].astype(str).to_numpy(), name="s"))
    comp = comp.reindex(index=domains, columns=samples, fill_value=0)
    umi = obs["total_counts"].to_numpy(dtype=float)
    gpu = obs["n_genes_by_counts"].to_numpy(dtype=float) / np.clip(umi, 1, None)
    med = pd.DataFrame({"d": ld, "umi": umi, "gpu": gpu}).groupby("d")[["umi", "gpu"]].median()
    prog = obs["program"].astype(str).to_numpy() if "program" in obs else None
    if prog is None:
        common.warn(progress, "Leiden domain annotation: no tissue assignment ('program'); no top programs")
    markers = marker_table(a.uns.get("rank_genes_groups"))
    if markers is None:
        common.warn(progress, "Leiden domain annotation: uns['rank_genes_groups'] does not hold the final domain "
                              "markers; no marker columns")
    summ = summarize_markers(markers, marker_padj)
    top = {d: (str(summ.at[d, "top_marker"]) if d in summ.index else "") for d in domains}
    det = _detection(a, ld, top)
    rows = []
    for d in domains:
        inside = ld == d
        n = int(inside.sum())
        row = {"domain": d, "manuscript_name": names.get(d, d), "n_bins": n}
        row.update({f"pct_{s}": _pct(comp.at[d, s], n) for s in samples})
        row.update({"median_umi": float(med.at[d, "umi"]), "median_genes_per_umi": round(float(med.at[d, "gpu"]), 3)})
        if prog is not None:
            row.update({"top_programs": top_programs(prog[inside]), "ambiguous_pct": _pct((prog[inside] == AMBIGUOUS)
                                                                                          .sum(), n)})
        else:
            row.update({"top_programs": "", "ambiguous_pct": np.nan})
        if markers is not None:
            row["n_significant_markers"] = int(summ.at[d, "n_significant_markers"]) if d in summ.index else 0
        else:
            row["n_significant_markers"] = np.nan
        row.update({"top_marker": top[d], "pct_domain_detecting_top_marker": det[d][0],
                    "pct_rest_detecting_top_marker": det[d][1]})
        rows.append(row)
    cols = (["domain", "manuscript_name", "n_bins"] + [f"pct_{s}" for s in samples]
            + ["median_umi", "median_genes_per_umi", "top_programs", "ambiguous_pct", "n_significant_markers",
               "top_marker", "pct_domain_detecting_top_marker", "pct_rest_detecting_top_marker"])
    df = pd.DataFrame(rows, columns=cols)
    df["flag"] = flag_domains(df["n_significant_markers"], df["pct_domain_detecting_top_marker"],
                              df["pct_rest_detecting_top_marker"]) if markers is not None else ""
    df["suggested_label"] = [_suggested(r) for r in df.itertuples(index=False)]
    return df


def _suggested(r):
    who = f"domain {r.domain}" if str(r.domain).isdigit() else str(r.domain)
    if r.flag == ONE_GENE:
        return f"{who}: {ONE_GENE} ({r.top_marker})"
    if r.flag:
        return f"{who}: {r.flag}"
    return f"{r.manuscript_name}: {r.top_programs}" if r.top_programs else str(r.manuscript_name)


# ---- 3. Leiden with and without Harmony -----------------------------------------------------------------------

def leiden_partition(embedding, n_neighbors=15, n_pcs=30, resolution=0.3, n_iterations=2, seed=0):
    """(labels, kNN connectivities) of the manuscript graph and Leiden call (bone_analysis.normalize_and_cluster:
    sc.pp.neighbors on the first n_pcs columns, igraph Leiden with the modularity objective) on a light AnnData object
    holding only the embedding."""
    import anndata as ad
    import scanpy as sc
    emb = np.asarray(embedding)
    light = ad.AnnData(obs=pd.DataFrame(index=pd.RangeIndex(emb.shape[0]).astype(str)), obsm={"rep": emb})
    sc.pp.neighbors(light, n_neighbors=int(n_neighbors), n_pcs=int(min(n_pcs, emb.shape[1])), use_rep="rep",
                    random_state=seed)
    sc.tl.leiden(light, resolution=resolution, key_added="leiden", flavor="igraph", objective_function="modularity",
                 n_iterations=n_iterations, random_state=seed, directed=False)
    return light.obs["leiden"].astype(str).to_numpy(), light.obsp["connectivities"]


def partition_section_stats(labels, sections, graph=None, samples=None, specific=SECTION_SPECIFIC):
    """Section structure of a partition: share of graph edges joining bins of the same section (and its value if edges
    ignored the section), clusters with >= `specific` of their bins from one section, the share of all bins in them
    and, per section, the share of its bins in clusters specific to that section."""
    labels, sections = np.asarray(labels).astype(str), np.asarray(sections).astype(str)
    samples = list(samples) if samples is not None else list(pd.unique(sections))
    ct = pd.crosstab(pd.Series(labels, name="c"), pd.Series(sections, name="s"))
    size = ct.sum(axis=1)
    spec = (ct.max(axis=1) / size) >= specific
    owner = ct.idxmax(axis=1)
    out = {"n_clusters": int(len(ct)), "same_section_edge_fraction": np.nan,
           "same_section_edge_fraction_random": round(float(((ct.sum(axis=0) / size.sum()) ** 2).sum()), 4),
           "n_section_specific_clusters": int(spec.sum()),
           "pct_bins_in_section_specific_clusters": _pct(size[spec].sum(), size.sum())}
    if graph is not None:
        g = graph.tocsr()
        codes = pd.factorize(sections)[0].astype(np.int32)       # integer codes: no large string gathers per edge
        rows = np.repeat(np.arange(g.shape[0], dtype=np.int32), np.diff(g.indptr))
        off = rows != g.indices
        if off.any():
            same = codes[rows] == codes[g.indices]
            out["same_section_edge_fraction"] = round(float(same[off].mean()), 4)
    for s in samples:
        n_s = int(ct[s].sum()) if s in ct.columns else 0
        own = ct.loc[spec & (owner == s), s].sum() if s in ct.columns else 0
        out[f"pct_{s}_in_section_specific_clusters"] = _pct(own, n_s)
    return out


def harmony_check(a, cfg, progress=None):
    """The manuscript graph and Leiden call (clustering.n_neighbors, n_pcs, resolution, leiden_n_iterations, seed) on
    X_pca_harmony and on X_pca. Raw partitions, before the unresolved and fragment merging. None without embeddings."""
    from sklearn.metrics import adjusted_rand_score
    clu, seed = cfg["clustering"], int(cfg["seed"])
    variants = [(v, k) for v, k in (("harmony", "X_pca_harmony"), ("no_harmony", "X_pca")) if k in a.obsm]
    if not variants:
        common.warn(progress, "Leiden without Harmony: no PCA embedding (Leiden step disabled); check skipped")
        return None
    if "X_pca_harmony" not in a.obsm:
        common.warn(progress, "Leiden without Harmony: no Harmony embedding (clustering.harmony off, harmonypy "
                              "missing or one section); only the PCA partition is described")
    samples = common.sample_ids(a, cfg)
    sections = a.obs["sample"].astype(str).to_numpy()
    final = a.obs["leiden"].astype(str).to_numpy() if "leiden" in a.obs else None
    used = a.uns.get("neighbors_rep")
    labels, rows = {}, []
    for variant, key in variants:
        lab, graph = leiden_partition(a.obsm[key], clu["n_neighbors"], clu["n_pcs"], clu["resolution"],
                                      clu["leiden_n_iterations"], seed)
        labels[variant] = lab
        row = {"variant": variant, "embedding": key, **partition_section_stats(lab, sections, graph, samples)}
        del graph
        if final is not None and key == used:
            kept = ~np.isin(final, [UNRESOLVED, UNASSIGNED])
            row["reproduces_manuscript_domains"] = bool((lab[kept] == final[kept]).all())
        rows.append(row)
    ref = labels.get("harmony")
    for row in rows:
        row["ari_vs_harmony"] = (round(float(adjusted_rand_score(ref, labels[row["variant"]])), 4)
                                 if ref is not None else np.nan)
    cols = (["variant", "embedding", "n_clusters", "ari_vs_harmony", "same_section_edge_fraction",
             "same_section_edge_fraction_random", "n_section_specific_clusters",
             "pct_bins_in_section_specific_clusters"] + [f"pct_{s}_in_section_specific_clusters" for s in samples]
            + ["reproduces_manuscript_domains"])
    return pd.DataFrame(rows).reindex(columns=cols)


# ---- module entry ---------------------------------------------------------------------------------------------

def run(a, cfg, progress=None):
    """The three Supplementary Figure S5 sensitivity tables (UNRESOLVED_CSV, ANNOTATION_CSV, HARMONY_CSV). Reads the
    joint object only; adds no obs column."""
    if "leiden" not in a.obs:
        common.warn(progress, "Leiden domains: no 'leiden' column (Leiden step disabled); clustering sensitivity "
                              "skipped")
        return []
    st = common.settings(cfg, "clustering")
    samples = common.sample_ids(a, cfg)
    lesion, lesion_src = lesion_flags(a, cfg, progress)
    lesion_txt = (f" Lesion bins ({lesion_src}): n_lesion, pct_lesion_unresolved (share of lesion bins that are "
                  f"unresolved) and pct_unresolved_lesion (share of unresolved bins that are lesion)."
                  if lesion is not None else "")
    outs = [common.Output(
        UNRESOLVED_CSV, "Unresolved bins per section",
        "Bins in the 'unresolved' Leiden domain (and 'unassigned' fragments) per section and pooled ('all'), over all "
        "bins, inside the region of interest (in_marrow_roi) and outside it; the manuscript quotes the pooled share "
        "only (Supplementary Figure S5)." + lesion_txt,
        unresolved_by_section(a.obs, samples, lesion))]
    padj = float(st["marker_padj"])
    outs.append(common.Output(
        ANNOTATION_CSV, "Leiden domain annotation",
        f"Per final Leiden domain: the Supplementary Figure S5 name (most enriched assigned program relative to the "
        f"whole dataset), share of its bins from each section, median UMIs and genes per UMI, the {TOP_PROGRAMS} most "
        f"frequent programs among its non-ambiguous assigned bins and the ambiguous share, final Wilcoxon markers with "
        f"BH-adjusted p < {padj:g} and positive score, the top marker (highest-scoring significant marker; empty, "
        f"with no detection, when there is none) and its detection (raw counts > 0) in and outside the domain; flag "
        f"'{NO_MARKER}' or '{ONE_GENE}' (<= {ONE_GENE_MAX_MARKERS} significant markers, top marker in >= "
        f"{ONE_GENE_MIN_DOMAIN_PCT:g}% of the domain and <= {ONE_GENE_MAX_REST_PCT:g}% elsewhere).",
        domain_annotation(a, cfg, padj, progress)))
    if common.as_bool([st["harmony_check"]])[0]:
        df = harmony_check(a, cfg, progress)
        if df is not None:
            clu = cfg["clustering"]
            outs.append(common.Output(
                HARMONY_CSV, "Leiden with and without Harmony",
                f"The manuscript neighbor graph ({clu['n_neighbors']} neighbors, {clu['n_pcs']} components) and Leiden "
                f"call (igraph, modularity, resolution {clu['resolution']}, {clu['leiden_n_iterations']} iterations, "
                f"seed {cfg['seed']}) on the Harmony-corrected (section label) and on the uncorrected PCA embedding. "
                f"Raw partitions, before the unresolved and fragment merging. ari_vs_harmony: adjusted Rand index "
                f"against the Harmony partition, to be read against the Leiden call's own variability (seed-to-seed "
                f"adjusted Rand index of the same graph and call, {STABILITY_CSV}, quoted in the figure caption); "
                f"same_section_edge_fraction: share of graph edges joining bins of the "
                f"same section (_random: if edges ignored the section); section-specific clusters: >= "
                f"{100 * SECTION_SPECIFIC:g}% of their bins from one section (per section: share of its bins in "
                f"clusters specific to it); reproduces_manuscript_domains: the partition on the manuscript embedding "
                f"equals the final domains outside 'unresolved' and 'unassigned'.", df))
    return outs
