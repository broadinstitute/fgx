# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "marimo<0.23.4",
#     "jedi<0.20.0",
#     "polars",
#     "httpx",
#     "altair",
#     "python-dotenv",
# ]
# ///

import marimo

__generated_with = "0.23.3"
app = marimo.App(width="medium")

with app.setup:
    import concurrent.futures as cf
    import os
    import sys
    import time
    from pathlib import Path

    import altair as alt
    import httpx
    import marimo as mo
    import polars as pl
    from dotenv import load_dotenv

    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
    GENEGENIE_TOKEN = os.environ.get("GENEGENIE_TOKEN")

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb01_pcsk9_walkthrough import BASE, client, fetch_tsv
    from nb02_variant_phewas import alt_alleles  # noqa: F401  (allele-fallback, see "To extend")
    from nb07_data_catalog import resource_metadata

    # Analysis constants, stated once so every count below is traceable to them.
    GW_MLOG10P = 7.30103  # p <= 5e-8, the conventional genome-wide threshold
    MAX_ABS_BETA = 1.0  # |log OR| <= 1, i.e. OR in [0.37, 2.72]
    MAX_SE = 0.5
    MIN_AAF = 0.001
    LOCUS_KB = 500

    # Statuses worth a retry; anything else is a bug in the request, not a blip.
    RETRY_STATUS = {429, 500, 502, 503, 504}


@app.function
def cancer_phenotype_codes(resource: str, prefixes=("C3_", "CD2_")) -> pl.DataFrame:
    """Every phenotype in `resource` whose FinnGen code starts with a neoplasm prefix.

    `C3_` is ICD-10 chapter C (malignant neoplasms); `CD2_` is in-situ / benign /
    uncertain-behaviour neoplasms. Returns code, label, and case counts. Note this
    is the *summary-statistics* phenotype list -- a large fraction of these have no
    fine-mapped credible sets, which `harvest_leads` discovers by probing.

    Gotcha: `/resource_metadata` is JSON, so column dtypes are inferred per resource and
    the count columns come back Int64 for some resources and String for others. Cast here
    or a `pl.concat` across resources raises SchemaError.
    """
    m = resource_metadata(resource)
    return (
        m.filter(pl.any_horizontal([pl.col("phenotype_code").str.starts_with(p) for p in prefixes]))
        .select("phenotype_code", "phenotype_string", "n_cases", "n_controls")
        .with_columns(pl.col("n_cases", "n_controls").cast(pl.Int64, strict=False))
    )


@app.function
def harvest_leads(targets: list[tuple[str, str]], max_workers: int = 12) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Fetch `credible_sets_by_phenotype_leads` for many (resource, phenotype) pairs.

    One row per credible set: the lead is the highest-PIP variant in that 95% set, so
    LD within a signal is already collapsed by fine-mapping. Returns (leads, coverage);
    coverage records the HTTP status per target so 404s are visible rather than silent.

    Gotcha: this endpoint answers JSON by default, not the API-wide TSV.
    """

    def one(t):
        res, code = t
        with client() as c:
            r = c.get(f"{BASE}/credible_sets_by_phenotype_leads/{res}/{code}", params={"format": "json"})
        payload = r.json() if r.status_code == 200 else []
        return res, code, r.status_code, payload

    rows, cov = [], []
    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        for res, code, status, payload in ex.map(one, targets):
            cov.append({"resource": res, "phenotype_code": code, "status": status, "n_leads": len(payload)})
            for d in payload:
                rows.append({**d, "resource": res, "phenotype_code": code})
    leads = pl.DataFrame(rows, infer_schema_length=None) if rows else pl.DataFrame()
    return leads, pl.DataFrame(cov)


@app.function
def cluster_loci(df: pl.DataFrame, kb: int = 500) -> pl.DataFrame:
    """Label each row with a `locus` id, merging leads within `kb` on a chromosome.

    Fine-mapping removes LD *within* a credible set, but the same signal reappears as a
    separate credible set in every phenotype definition and every resource that carries it.
    Distance clustering is the crude second pass that turns those back into one locus.
    """
    return (
        df.sort(["chr", "pos"])
        .with_columns(
            (
                (pl.col("pos") - pl.col("pos").shift(1).over("chr") > kb * 1000)
                | pl.col("pos").shift(1).over("chr").is_null()
            )
            .cum_sum()
            .over("chr")
            .alias("locus_idx")
        )
        .with_columns(pl.concat_str([pl.col("chr"), pl.col("locus_idx")], separator="_").alias("locus"))
    )


@app.function
def qc_leads(malignant: pl.DataFrame, gw: float = GW_MLOG10P, kb: int = LOCUS_KB) -> pl.DataFrame:
    """Apply the stated QC filters, dedupe, and attach `locus` labels. One code path.

    The clustering population is *all* QC-passing leads, risk and protective alike; the
    `beta < 0` filter belongs after the locus labels exist. Order matters: cluster the
    protective rows on their own and a risk lead that had been bridging two protective
    leads is gone, so the same window silently yields more loci. The headline count and
    the sensitivity grid both call this, so they cannot drift apart.
    """
    return cluster_loci(
        malignant.filter(
            (pl.col("mlog10p") >= gw)
            & (pl.col("beta").abs() <= MAX_ABS_BETA)
            & (pl.col("se") <= MAX_SE)
            & (pl.col("aaf") >= MIN_AAF)
        ).unique(subset=["chr", "pos", "ref", "alt", "phenotype_code", "beta", "se", "mlog10p"], keep="first"),
        kb,
    )


@app.function
def protective_loci_at(malignant: pl.DataFrame, gw: float, kb: int) -> int:
    """Independent protective loci at significance `gw` and clustering window `kb`."""
    return qc_leads(malignant, gw, kb).filter(pl.col("beta") < 0)["locus"].n_unique()


@app.function
def batch_credible_sets_by_variant(
    variants: list[str], chunk: int = 40, max_workers: int = 6, attempts: int = 4
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """POST /credible_sets_by_variant for many variants; the PheWAS arm of this notebook.

    **Fails closed.** A chunk that never answers 200 raises rather than returning `[]`.
    This matters more here than anywhere else in the catalog: downstream, a variant with no
    returned rows is read as "no other-disease association", so one silently dropped chunk
    would relabel up to `chunk` variants as clean. Transient statuses (429, 5xx) and
    transport errors are retried with exponential backoff first.

    Returns `(rows, coverage)`; `coverage` is one row per chunk so the caller can assert
    full coverage before computing any verdict.

    Gotcha: the `variants` body field is a single string whose separator is a NEWLINE.
    Commas, spaces and semicolons all 422 with "variant needs to contain four fields".
    """

    def one(item):
        idx, vs = item
        why = "no attempt made"
        for attempt in range(1, attempts + 1):
            try:
                with client() as c:
                    r = c.post(
                        f"{BASE}/credible_sets_by_variant",
                        json={"variants": "\n".join(vs)},
                        params={"format": "json"},
                    )
                if r.status_code == 200:
                    return {"chunk": idx, "n_variants": len(vs), "attempts": attempt}, r.json()
                why = f"HTTP {r.status_code}: {r.text[:160]}"
                if r.status_code not in RETRY_STATUS:
                    break
            except httpx.HTTPError as exc:  # timeouts, connection resets
                why = repr(exc)
            if attempt < attempts:
                time.sleep(2 ** (attempt - 1))
        raise RuntimeError(f"/credible_sets_by_variant chunk {idx} ({len(vs)} variants) failed -- {why}")

    chunks = list(enumerate([variants[i : i + chunk] for i in range(0, len(variants), chunk)]))
    rows, cov = [], []
    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        for meta, payload in ex.map(one, chunks):  # ex.map re-raises here, so a dead chunk stops the notebook
            cov.append({**meta, "n_rows": len(payload)})
            rows.extend(payload)
    return (
        pl.DataFrame(rows, infer_schema_length=None) if rows else pl.DataFrame(),
        pl.DataFrame(cov).sort("chunk"),
    )


@app.function
def cancer_site(label: str) -> str:
    """Collapse a FinnGen phenotype label to a canonical cancer site.

    FinnGen ships the same cancer under many case/control definitions (`_WIDE` adds Hilmo
    hospital records, `_EXALLC` excludes other cancers from controls, `_INCLAVO` adds
    primary-care records) plus histological subtypes. Counting those as distinct diseases
    is the single largest source of inflation in this analysis. First match wins.
    """
    rules = [
        ("basal cell", "Skin (non-melanoma)"),
        ("squamous cell neoplasms and carcinomas of the skin", "Skin (non-melanoma)"),
        ("non-melanoma", "Skin (non-melanoma)"),
        ("melanoma", "Melanoma"),
        ("neoplasm of skin", "Skin (non-melanoma)"),
        ("breast", "Breast"),
        ("prostate", "Prostate"),
        ("colorect", "Colorectal"),
        ("colon", "Colorectal"),
        ("rectum", "Colorectal"),
        ("bronchus and lung", "Lung"),
        ("lung cancer", "Lung"),
        ("ovar", "Ovary"),
        ("endometri", "Uterus/Endometrium"),
        ("corpus uteri", "Uterus/Endometrium"),
        ("uter", "Uterus/Endometrium"),
        ("cervi", "Cervix"),
        ("bladder", "Bladder"),
        ("kidney", "Kidney"),
        ("renal", "Kidney"),
        ("thyroid", "Thyroid"),
        ("stomach", "Stomach"),
        ("oesophag", "Oesophagus"),
        ("pancrea", "Pancreas"),
        ("liver", "Liver"),
        ("hepatocellular", "Liver"),
        ("testis", "Testis"),
        ("testic", "Testis"),
        ("brain", "Brain/CNS"),
        ("meningi", "Brain/CNS"),
        ("glio", "Brain/CNS"),
        ("nervous system", "Brain/CNS"),
        ("leukaemia", "Leukaemia"),
        ("leukemia", "Leukaemia"),
        ("lymphoma", "Lymphoma"),
        ("hodgkin", "Lymphoma"),
        ("myeloma", "Myeloma"),
        ("gammopathy", "Myeloma"),
        ("lymphoid", "Haematological (other)"),
        ("haematopoietic", "Haematological (other)"),
        ("myelo", "Haematological (other)"),
        ("head and neck", "Head & neck"),
        ("larynx", "Head & neck"),
        ("lip, oral", "Head & neck"),
        ("salivary", "Head & neck"),
        ("tongue", "Head & neck"),
        ("bone and articular", "Bone/soft tissue"),
        ("connective and soft tissue", "Bone/soft tissue"),
        ("small intestine", "Small intestine"),
        ("malignant neoplasm", "All cancer (umbrella)"),
    ]
    low = label.lower()
    for kw, site in rules:
        if kw in low:
            return site
    return "Other/unmapped"


@app.cell
def _():
    mo.md(r"""
    # nb11: How many protective variants are there for cancer?

    A deceptively simple two-part question, and a good stress test of the API:

    1. **How many protective variants are there for any cancer?**
    2. **Of those, which are associated with other diseases?**

    Both are answerable from the fine-mapped GWAS already in the API. The interesting part
    is not the number -- it is that the number is only meaningful once you say what a
    "variant" is, what "protective" means when allele labels are arbitrary, and which of
    FinnGen's several-hundred cancer endpoints are actually different diseases.

    This notebook answers both, states every threshold, and ends with the part that matters
    most: what these data **cannot** tell you.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## 0. What cancer data is actually here?

    Before counting anything: which resources carry fine-mapped cancer GWAS, and how many
    cancer endpoints does each define? (`nb07` is the general introspection notebook; we
    reuse its `resource_metadata` helper.)
    """)
    return


@app.cell
def _():
    RESOURCES = ["finngen", "finngen_mvp_ukbb", "finngen_ukbb"]
    catalog = pl.concat([cancer_phenotype_codes(r).with_columns(pl.lit(r).alias("resource")) for r in RESOURCES])
    catalog_summary = (
        catalog.with_columns(
            pl.when(pl.col("phenotype_code").str.starts_with("C3_"))
            .then(pl.lit("C3_ malignant"))
            .otherwise(pl.lit("CD2_ in-situ / benign"))
            .alias("class")
        )
        .group_by(["resource", "class"])
        .agg(pl.len().alias("n_endpoints"), pl.max("n_cases").alias("largest_n_cases"))
        .sort(["resource", "class"])
    )
    mo.vstack(
        [
            mo.md(
                "**Cancer endpoints defined per resource.** `finngen` is R14 (Finland only, "
                "deepest phenotyping); `finngen_ukbb` and `finngen_mvp_ukbb` are the R13 "
                "meta-analyses that add UK Biobank and the Million Veteran Program -- far fewer "
                "endpoints, far larger samples."
            ),
            catalog_summary,
        ]
    )
    return RESOURCES, catalog


@app.cell
def _(RESOURCES, catalog):
    targets = [(r["resource"], r["phenotype_code"]) for r in catalog.iter_rows(named=True)]
    raw_leads, coverage = harvest_leads(targets)
    cov_summary = (
        coverage.group_by("resource")
        .agg(
            pl.len().alias("probed"),
            (pl.col("status") == 200).sum().alias("has_credible_sets"),
            pl.sum("n_leads").alias("lead_rows"),
        )
        .sort("resource")
    )
    mo.vstack(
        [
            mo.md(
                f"### Probed all {len(targets):,} cancer endpoints across {len(RESOURCES)} resources\n\n"
                "**First real finding, and it is a data-availability one.** The phenotype list and "
                "the fine-mapping store are not the same thing: most FinnGen R14 cancer endpoints "
                "return `404` from `credible_sets_by_phenotype_leads` because they were never "
                "fine-mapped (too few cases). Any tool that answers this question off the phenotype "
                "list alone will overstate its coverage."
            ),
            cov_summary,
            mo.md(f"**Total credible-set leads retrieved: {len(raw_leads):,}**"),
        ]
    )
    return (raw_leads,)


@app.cell
def _():
    mo.md(r"""
    ## 1. The direction-of-effect problem

    Before any count, the question hides a trap.

    GWAS effect sizes are reported **per alternate allele**. "Protective" therefore means
    `beta < 0`: carrying the alt allele lowers risk. But which allele is called *alt* is a
    reference-genome convention, not biology. Every risk locus is also a protective locus
    read from the other allele. So for common variants, "how many protective variants" is
    close to a restatement of "how many cancer loci", divided by two.

    The exception is informative. Plot the fraction of significant leads with `beta < 0`
    against alt-allele frequency:
    """)
    return


@app.cell
def _(raw_leads):
    direction_probe = (
        raw_leads.filter(pl.col("mlog10p") >= GW_MLOG10P)
        .with_columns(
            pl.when(pl.col("aaf") < 0.01)
            .then(pl.lit("<1%"))
            .when(pl.col("aaf") < 0.05)
            .then(pl.lit("1-5%"))
            .when(pl.col("aaf") < 0.5)
            .then(pl.lit("5-50%"))
            .otherwise(pl.lit(">50% (alt is major)"))
            .alias("alt_allele_freq")
        )
        .group_by("alt_allele_freq")
        .agg(pl.len().alias("n_leads"), pl.col("beta").lt(0).mean().alias("frac_protective"))
        .sort("alt_allele_freq")
    )
    direction_chart = (
        alt.Chart(direction_probe)
        .mark_bar(size=42, color="#2b83ba")
        .encode(
            x=alt.X(
                "alt_allele_freq:N", sort=["<1%", "1-5%", "5-50%", ">50% (alt is major)"], title="alt allele frequency"
            ),
            y=alt.Y("frac_protective:Q", title="fraction of leads with beta < 0", scale=alt.Scale(domain=[0, 0.6])),
            tooltip=["alt_allele_freq", "n_leads", alt.Tooltip("frac_protective:Q", format=".3f")],
        )
        .properties(height=260, title="Protective fraction rises with allele frequency -- ascertainment, not biology")
    )
    rule = alt.Chart(pl.DataFrame({"y": [0.5]})).mark_rule(strokeDash=[6, 4], color="#d7191c").encode(y="y:Q")
    mo.vstack(
        [
            mo.ui.altair_chart(direction_chart + rule),
            direction_probe,
            mo.md(
                "At **alt-allele frequency above 50%** the split is ~52/48 -- a coin flip, exactly what "
                "arbitrary allele labelling predicts. It falls to **~14% at frequencies under 1%**.\n\n"
                "That gradient is a detection asymmetry, not a biological one. To reach genome-wide "
                "significance a rare allele must be *enriched* in cases; a rare allele *depleted* in "
                "cases has almost no power at the same frequency. **Rare protective variants are "
                "systematically under-ascertained in this data**, so any count below is a floor for "
                "common variants and badly incomplete for rare ones."
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(rf"""
    ## 2. Filters, stated up front

    | filter | value | why |
    |---|---|---|
    | malignant only | `C3_` prefix | `CD2_` is in-situ/benign (uterine leiomyoma, colon polyps). Not cancer. |
    | significance | `mlog10p >= {GW_MLOG10P}` (p <= 5e-8) | the leads endpoint returns credible sets regardless of significance -- the raw minimum here is p ~ 0.98 |
    | effect plausibility | `abs(beta) <= {MAX_ABS_BETA}` | log-odds beyond +/-1 at these frequencies is model separation, not biology |
    | standard error | `se <= {MAX_SE}` | same |
    | allele frequency | `aaf >= {MIN_AAF}` | below 0.1% the effect estimate is not identifiable in these sample sizes |
    | direction | `beta < 0` | alt allele lowers risk |
    | LD | credible-set leads, then {LOCUS_KB}kb clustering | fine-mapping handles LD within a signal; clustering handles the same signal recurring across endpoints |

    The plausibility filters are not cosmetic. Watch what they remove:
    """)
    return


@app.cell
def _(raw_leads):
    malignant = raw_leads.filter(pl.col("phenotype_code").str.starts_with("C3_"))
    significant = malignant.filter(pl.col("mlog10p") >= GW_MLOG10P)
    implausible = significant.filter(
        (pl.col("beta").abs() > MAX_ABS_BETA) | (pl.col("se") > MAX_SE) | (pl.col("aaf") < MIN_AAF)
    )
    mo.vstack(
        [
            mo.md(
                f"### {len(implausible):,} significant leads fail the plausibility filters\n\n"
                "Sorted by effect size, they are almost all the same thing: near-monomorphic alleles "
                "at 16q24.3 (the *MC1R* pigmentation region) in skin-cancer endpoints, where the "
                "logistic model is separating rather than estimating. A beta of -30 is an odds ratio "
                "of 1e-13. The locus is real biology; these particular estimates are not usable numbers."
            ),
            implausible.sort("beta")
            .select("phenotype_code", "resource", "chr", "pos", "beta", "se", "aaf", "mlog10p", "gene_most_severe")
            .head(10),
        ]
    )
    return (malignant,)


@app.cell
def _(malignant):
    qc = qc_leads(malignant).with_columns(
        pl.concat_str(["chr", "pos", "ref", "alt"], separator=":").alias("variant"),
        (pl.col("beta") < 0).alias("alt_protective"),
        pl.col("trait").str.replace_all("_", " ").alias("trait_label"),
    )
    qc = qc.with_columns(pl.col("trait_label").map_elements(cancer_site, return_dtype=pl.String).alias("site"))
    protective = qc.filter(pl.col("alt_protective"))
    mo.md(
        f"**After QC: {len(qc):,} cancer associations**, of which **{len(protective):,} are protective** "
        f"({len(protective) / len(qc):.1%}) and {len(qc) - len(protective):,} are risk-increasing."
    )
    return protective, qc


@app.cell
def _():
    mo.md(r"""
    ## 3. The phenotype-redundancy problem

    Now the second trap, and the one that would most distort a naive answer.

    Group the protective associations by cancer *site* rather than by phenotype code:
    """)
    return


@app.cell
def _(protective):
    by_site = (
        protective.group_by("site")
        .agg(
            pl.len().alias("associations"),
            pl.col("variant").n_unique().alias("distinct_variants"),
            pl.col("locus").n_unique().alias("independent_loci"),
            pl.col("phenotype_code").n_unique().alias("finngen_codes"),
        )
        .sort("associations", descending=True)
    )
    skin_row = by_site.filter(pl.col("site") == "Skin (non-melanoma)")
    mo.vstack(
        [
            by_site,
            mo.md(
                "**Non-melanoma skin cancer alone accounts for "
                f"{int(skin_row['associations'][0]):,} of {len(protective):,} protective associations "
                f"({int(skin_row['associations'][0]) / len(protective):.0%})** -- but only "
                f"{int(skin_row['independent_loci'][0]):,} independent loci.\n\n"
                "That is not skin cancer being unusually genetic. It is FinnGen defining it "
                f"{int(skin_row['finngen_codes'][0])} different ways (`C3_SKIN`, `C3_OTHER_SKIN_WIDE`, "
                "`C3_BASAL_CELL_CARCINOMA_WIDE`, `C3_SQUAMOUS_CELL_CARCINOMA_SKIN_WIDE`, each with "
                "`_WIDE`/`_EXALLC` control variants) and three overlapping resources each reporting "
                "it. **This is the garbage-with-the-signal problem in one number.**"
            ),
        ]
    )
    return


@app.cell
def _(protective, qc):
    counts = pl.DataFrame(
        {
            "unit": [
                "variant x cancer-endpoint associations",
                "distinct protective variants",
                "distinct (locus x cancer site) signals",
                "independent protective loci",
                "loci protective for every cancer they hit",
            ],
            "n": [
                len(protective),
                protective["variant"].n_unique(),
                protective.select(["locus", "site"]).unique().height,
                protective["locus"].n_unique(),
                qc.group_by("locus").agg(pl.col("alt_protective").mean().alias("f")).filter(pl.col("f") == 1.0).height,
            ],
        }
    )
    mixed = (
        qc.group_by("locus")
        .agg(pl.col("alt_protective").mean().alias("f"))
        .filter((pl.col("f") > 0) & (pl.col("f") < 1))
        .height
    )
    mo.vstack(
        [
            mo.md("## Answer to Q1\n\nThe count spans an order of magnitude depending on the unit:"),
            counts,
            mo.md(
                f"**The defensible headline is ~{protective['locus'].n_unique():,} independent protective loci**, "
                f"of which {counts['n'][4]:,} point protective at every cancer they reach. The other "
                f"{mixed:,} loci lower risk for one cancer while raising it for another -- so even "
                "'protective' is site-specific, not a property of the variant."
            ),
        ]
    )
    return


@app.cell
def _(malignant):
    grid = [
        {
            "p threshold": f"p<=1e-{gw:.0f}" if gw > 7.4 else "p<=5e-8",
            "locus window (kb)": kb,
            "protective loci": protective_loci_at(malignant, gw, kb),
        }
        for gw in (GW_MLOG10P, 8.0, 9.0, 12.0)
        for kb in (100, 500, 1000)
    ]
    sens = pl.DataFrame(grid).pivot(on="locus window (kb)", index="p threshold", values="protective loci")
    grid_counts = [g["protective loci"] for g in grid]
    grid_same_kb = [g["protective loci"] for g in grid if g["locus window (kb)"] == LOCUS_KB]
    mo.vstack(
        [
            mo.md(
                "### How fragile is that number?\n\n"
                "Protective loci as a function of the two choices that actually move it -- significance "
                "threshold and locus window. Columns are the clustering window in kb. The grid runs the "
                f"same `protective_loci_at` code path as the headline, so the p<=5e-8 / {LOCUS_KB}kb cell "
                "reproduces the headline exactly rather than approximating it."
            ),
            sens,
            mo.md(
                f"**A {max(grid_counts) / min(grid_counts):.1f}-fold range across the whole grid** "
                f"({max(grid_counts):,} at the loosest corner, {min(grid_counts):,} at the tightest), and "
                f"{max(grid_same_kb) / min(grid_same_kb):.1f}-fold if you hold the {LOCUS_KB}kb window and "
                "vary only the threshold. The number is threshold-dependent in the way every GWAS count "
                "is; it is not knife-edge. The plausibility filters "
                "(`beta`, `se`, `aaf`) barely move it at all -- they remove ~12% of associations but "
                "under 1% of loci, because the implausible rows pile onto loci that are already counted."
            ),
        ]
    )
    return


@app.cell
def _(protective):
    # Every distinct protective variant is queried, not one representative per locus.
    # Distance clustering says two leads are near each other; it does not say they tag the
    # same causal signal, so one lead's PheWAS is not a valid proxy for the other's. Query
    # all of them and aggregate to the locus afterwards -- that direction is safe, the
    # reverse is not.
    probes = (
        protective.sort("mlog10p", descending=True)
        .unique(subset=["variant"], keep="first")
        .select(
            "variant",
            "locus",
            "chr",
            "pos",
            "site",
            pl.col("beta").alias("cancer_beta"),
            pl.col("mlog10p").alias("cancer_mlog10p"),
            pl.col("gene_most_severe").alias("gene"),
        )
    )
    # Locus-level labels, taken from the strongest protective lead at the locus. Used for
    # display and grouping only; the PheWAS itself never collapses to these.
    loci = (
        protective.sort("mlog10p", descending=True)
        .unique(subset=["locus"], keep="first")
        .select(
            "locus",
            pl.col("gene_most_severe").alias("locus_gene"),
            pl.col("site").alias("locus_site"),
        )
    )
    n_loci = loci.height
    mo.md(
        f"### Q2 input: {len(probes):,} distinct protective variants across {n_loci:,} loci\n\n"
        f"An earlier draft sent one representative variant per locus ({n_loci:,} queries instead of "
        f"{len(probes):,}) and then labelled the whole locus from that one PheWAS. That is only valid if "
        "co-located leads are in tight LD, which 500kb distance clustering does not establish -- two "
        "independent signals 300kb apart can have entirely different pleiotropy. So every protective "
        "variant is queried and the union is taken per locus. Doing it the other way understates "
        "trade-offs, because a trade-off found at any lead in the locus is missed unless that lead "
        "happened to be the representative."
    )
    return loci, n_loci, probes


@app.cell
def _():
    mo.md(r"""
    ## 4. Q2: what else do these alleles do?

    For each protective lead we ask `credible_sets_by_variant`: which other traits have a
    credible set containing this variant? Because effect sizes are reported per alt allele
    throughout, the comparison is direct -- the same allele that gave `beta < 0` for cancer
    either lowers (`beta < 0`) or raises (`beta > 0`) the other trait.
    """)
    return


@app.cell
def _(probes):
    phewas_raw, phewas_cov = batch_credible_sets_by_variant(probes["variant"].to_list())
    # Fail closed a second time, at the analysis boundary: no verdict is computed unless the
    # chunks together account for every variant we asked about. A missing chunk would read as
    # "these variants have no other-disease association", which is the one error that silently
    # inflates the reassuring answer.
    assert phewas_cov["n_variants"].sum() == probes.height, "PheWAS coverage is incomplete"
    phewas = (
        phewas_raw.with_columns(pl.concat_str(["chr", "pos", "ref", "alt"], separator=":").alias("variant"))
        .join(probes.drop("chr", "pos"), on="variant", how="inner")
        .filter((pl.col("data_type") == "GWAS") & (pl.col("mlog10p") >= GW_MLOG10P))
        .with_columns(
            pl.any_horizontal([pl.col("trait_original").str.starts_with(p) for p in ("C3_", "CD2_")]).alias(
                "is_cancer_trait"
            )
        )
    )
    layer_mix = phewas_raw.group_by("data_type").agg(pl.len().alias("rows")).sort("rows", descending=True)
    mo.vstack(
        [
            mo.md(
                f"**{len(phewas_raw):,} credible-set rows returned** for all {probes.height:,} "
                f"protective variants, sent as {phewas_cov.height} chunks with every chunk confirmed 200 "
                "(`batch_credible_sets_by_variant` raises rather than returning an empty payload, so a "
                "throttled or expired request cannot masquerade as 'no association'). Rows span GWAS and "
                "the molecular QTL layers; we keep only genome-wide significant GWAS rows for the disease "
                "question, and leave the QTL mechanism layer to a follow-on notebook."
            ),
            layer_mix,
        ]
    )
    return (phewas,)


@app.cell
def _(phewas):
    non_cancer_all = phewas.filter(~pl.col("is_cancer_trait"))
    naive_split = (
        non_cancer_all.with_columns(
            pl.when(pl.col("beta") < 0).then(pl.lit("also lower")).otherwise(pl.lit("higher")).alias("direction")
        )
        .group_by("direction")
        .agg(pl.len().alias("n"))
    )
    mo.vstack(
        [
            mo.md(
                f"### The naive answer, and why it is wrong\n\n"
                f"{len(non_cancer_all):,} non-cancer associations. Split by direction:"
            ),
            naive_split,
            mo.md(
                "**A near-perfect 50/50 -- which should be a red flag, not a finding.** Most of these "
                "'traits' are not diseases. They are Kanta lab measurements with numeric codes "
                "(`3026361`), ATC drug-purchase endpoints (`ATC_H03AA_IRN`), and Open Targets studies "
                "(`GCST...`) that are largely quantitative. For a lab value, 'higher' has no intrinsic "
                "direction of harm, so counting sign flips over them measures nothing.\n\n"
                "To ask the question honestly we have to restrict to **binary disease endpoints**, "
                "where raising the trait unambiguously means more disease."
            ),
        ]
    )
    return (non_cancer_all,)


@app.cell
def _(RESOURCES, non_cancer_all):
    binary_codes = (
        pl.concat(
            [resource_metadata(r).select("phenotype_code", "trait_type") for r in RESOURCES],
            how="vertical_relaxed",
        )
        .filter(pl.col("trait_type") == "binary")["phenotype_code"]
        .unique()
        .to_list()
    )
    disease = non_cancer_all.filter(
        (~pl.col("trait_original").str.starts_with("GCST"))
        & (~pl.col("trait_original").str.starts_with("ATC_"))
        & (~pl.col("trait_original").str.contains("_IRN"))
        & (~pl.col("trait_original").str.contains(r"^\d+$"))
        & (pl.col("trait_original").is_in(binary_codes))
    ).with_columns(
        pl.when(pl.col("beta") < 0).then(pl.lit("also lower")).otherwise(pl.lit("higher")).alias("direction"),
        pl.col("trait").str.replace_all("_", " ").alias("other_disease"),
    )
    disease_split = disease.group_by("direction").agg(pl.len().alias("n"))
    mo.vstack(
        [
            mo.md(
                f"### Restricted to binary disease endpoints: {len(disease):,} associations "
                f"across {disease['trait_original'].n_unique():,} distinct diseases"
            ),
            disease_split,
            mo.md(
                "Now the split is real and it leans the wrong way: cancer-protective alleles raise "
                "another disease more often than they lower one."
            ),
        ]
    )
    return (disease,)


@app.cell
def _(disease, loci, n_loci, probes):
    # Counts are of distinct diseases, not association rows: the same phenotype recurs across
    # resources and across the several variants now probed at each locus, and summing rows
    # would let one disease counted three times look like three trade-offs.
    def verdicts(rows: pl.DataFrame) -> pl.DataFrame:
        return (
            rows.group_by("locus")
            .agg(
                pl.col("variant").n_unique().alias("n_variants_with_signal"),
                pl.col("trait_original").filter(pl.col("direction") == "higher").n_unique().alias("n_worse"),
                pl.col("trait_original").filter(pl.col("direction") == "also lower").n_unique().alias("n_better"),
            )
            .join(loci, on="locus", how="left")
        )

    per_locus = verdicts(disease)
    clean_loci = per_locus.filter(pl.col("n_worse") == 0)
    tradeoff_loci = per_locus.filter(pl.col("n_worse") > 0)
    silent = n_loci - per_locus.height

    # What the representative-variant shortcut would have concluded, from the same PheWAS rows
    # and no extra API calls: keep only each locus's strongest protective lead. This is the
    # cost of the shortcut measured rather than argued.
    rep_variants = probes.sort("cancer_mlog10p", descending=True).unique(subset=["locus"], keep="first")["variant"]
    rep_per_locus = verdicts(disease.filter(pl.col("variant").is_in(rep_variants)))
    rep_tradeoff = rep_per_locus.filter(pl.col("n_worse") > 0).height
    rep_silent = n_loci - rep_per_locus.height
    verdict = pl.DataFrame(
        {
            "verdict": [
                "no significant other-disease association",
                "cleanly protective (lowers cancer, raises no disease)",
                "trade-off (lowers cancer, raises another disease)",
            ],
            "loci": [silent, len(clean_loci), len(tradeoff_loci)],
        }
    )
    mo.vstack(
        [
            mo.md(f"## Answer to Q2\n\nOf the {n_loci:,} protective cancer loci:"),
            verdict,
            mo.md(
                f"**{len(tradeoff_loci) / (len(clean_loci) + len(tradeoff_loci)):.0%} of the loci that "
                "have any other-disease signal at all are trade-offs.** A locus is classified from the "
                "union over every protective variant it contains, so 'trade-off' is easy to reach and "
                "'cleanly protective' is the strictly harder claim -- which is the right way round for "
                "an error to fall. It is still the weaker of the two: absence of a significant "
                "association is mostly absence of power, not evidence of no effect."
            ),
            mo.md(
                "**What the shortcut would have cost.** Re-scoring the same PheWAS rows using only each "
                f"locus's strongest protective lead -- the representative-variant design -- gives "
                f"{rep_tradeoff:,} trade-off loci instead of {len(tradeoff_loci):,}, and {rep_silent:,} "
                f"apparently silent loci instead of {silent:,}. The single worst case is the *APOE* "
                "region: its strongest protective lead (19:44913484, liver cancer) raises three "
                "diseases, while the neighbouring protective lead 5 kb away sits on *APOE* itself and "
                "raises 46, including the entire dementia family at `mlog10p` ~1,600. Query one and "
                "call the locus described and you miss the largest trade-off in the dataset."
            ),
        ]
    )
    return (tradeoff_loci,)


@app.cell
def _(disease, tradeoff_loci):
    worse = (
        disease.filter(pl.col("direction") == "higher")
        .join(tradeoff_loci.select("locus", "locus_gene", "locus_site"), on="locus", how="inner")
        .sort("mlog10p", descending=True)
        # One row per (locus, disease), keeping the strongest. Deduping on gene instead would
        # merge two genuinely distinct loci that share a nearest-gene annotation.
        .unique(subset=["locus", "trait_original"], keep="first")
        .sort("mlog10p", descending=True)
    )
    mo.vstack(
        [
            mo.md(
                "### The trade-offs, strongest first\n\n"
                "One row per (locus, disease) pair. `gene` is the nearest/most-severe gene at the "
                "cancer locus; `other_disease` is what the same allele raises."
            ),
            worse.select(
                pl.col("locus_gene").alias("gene"),
                pl.col("locus_site").alias("site"),
                "other_disease",
                "beta",
                "mlog10p",
            ).head(30),
        ]
    )
    return


@app.cell
def _(disease, tradeoff_loci):
    # The axis says "distinct diseases", so count distinct diseases. `pl.len()` here would count
    # association rows, and the same phenotype recurs across finngen / finngen_ukbb /
    # finngen_mvp_ukbb and across the several protective variants at a locus -- enough to
    # reorder the bars, not just inflate them. `association_rows` is kept in the tooltip so the
    # gap between the two is visible rather than hidden.
    top_genes = (
        disease.join(tradeoff_loci.select("locus", "locus_gene"), on="locus", how="inner")
        .filter((pl.col("direction") == "higher") & pl.col("locus_gene").is_not_null())
        .group_by(pl.col("locus_gene").alias("gene"))
        .agg(
            pl.col("trait_original").n_unique().alias("n"),
            pl.len().alias("association_rows"),
            pl.col("mlog10p").max().alias("max_mlog10p"),
        )
        .sort("n", descending=True)
        .head(15)
    )
    tradeoff_chart = (
        alt.Chart(top_genes)
        .mark_bar(color="#d7191c", opacity=0.85)
        .encode(
            x=alt.X("n:Q", title="distinct diseases raised by the cancer-protective allele"),
            y=alt.Y("gene:N", sort="-x", title=None),
            tooltip=["gene", "n", "association_rows", alt.Tooltip("max_mlog10p:Q", format=".1f")],
        )
        .properties(height=380, title="Loci that buy cancer protection at a cost")
    )
    mo.vstack(
        [
            mo.ui.altair_chart(tradeoff_chart),
            mo.md(
                "**Read `n` as an upper bound.** Section 3 showed FinnGen defining one cancer several "
                "ways; it does the same on the other-disease side. *PHTF1*'s 64 'distinct diseases' "
                "include `Hypothyroidism, strict autoimmune`, `Hypothyroidism, drug reimbursement` and "
                "`Disorders of the thyroid gland` as three. Counting association rows instead would be "
                "worse still -- it also multiplies by resource, which is why *CDKN2B-AS1* outranks "
                "*PHTF1* on rows (99 vs 98) and loses on diseases (47 vs 64). Collapsing the "
                "non-cancer endpoints into disease families, the way `cancer_site` does for cancers, "
                "is the missing piece; it needs a mapping this notebook does not build."
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ### What the top trade-offs actually are

    These are not statistical curiosities -- they are known immunology and known pleiotropy,
    which is the best available check that the pipeline is doing something real:

    - **`PHTF1` (1p13.2)** is the *PTPN22* region, and the widest-reaching trade-off here.
      The allele that lowers skin-cancer risk raises autoimmune hypothyroidism, rheumatoid
      arthritis, type 1 diabetes, and Crohn's. Sharper immune surveillance, more
      autoimmunity -- the classic immune set-point trade.
    - **`AC011481.3` / *APOE* (19q13)** is the strongest single row in the notebook. A lead
      protective for liver cancer raises dementia and Alzheimer's at `mlog10p` ~1,600.
      The gene label is the annotation caveat in miniature: the locus is *APOE*, and the
      lead carrying the dementia signal is annotated to a neighbouring lncRNA.
    - **`CDKN2B-AS1` (9p21)** carries eight distinct protective leads spanning *CDKN2A* and
      *CDKN2B-AS1* across four cancer sites (non-melanoma skin, melanoma, brain/CNS,
      colorectal), and raises coronary atherosclerosis, ischaemic heart disease and primary
      open-angle glaucoma. The best argument in the notebook for querying every lead rather
      than one per locus.
    - **`PTCSC2` (9q22, near *FOXE1*)** lowers thyroid cancer and raises hypothyroidism and
      goitre. Same tissue, opposite ends of thyroid function.
    - **`IL2RA`, `BACH2`, `CDKAL1`, `RAB5B`, `UBE2L3`** (skin) and **`IRF5`** (kidney) --
      immune-regulation loci, all the same shape: less cancer at the site they hit, more
      autoimmune thyroid disease, type 1 diabetes, rheumatoid arthritis, Crohn's.
    - Outside the immune block: **`SPDL1`** lowers cancer overall and raises **idiopathic
      pulmonary fibrosis**, one of the cleanest known cancer/fibrosis antagonisms;
      **`VAMP8`** lowers prostate cancer and raises coronary heart disease; **`FTO`** lowers
      recorded breast-cancer risk and raises type 2 diabetes, hypertension and arthrosis,
      which is the adiposity axis rather than an immune one.

    **The recurring theme is immune tone.** A large share of what looks like "protection from
    cancer" in this data is a dial on immune activity, and turning it up costs autoimmunity.
    That is a coherent biological story, and it is the kind of answer the question was
    actually reaching for.

    **Some entries point the other way and should be treated as open questions.** *ABO* here
    is protective for pancreatic cancer while raising pulmonary embolism, and *CHEK2* is
    protective for lung cancer while raising myeloproliferative disease -- both loci are
    genuinely pleiotropic, but neither direction is the textbook one, and co-location in a
    500kb window is not colocalization. These are exactly the rows to send through
    `colocalization_by_variant` before believing them.
    """)
    return


@app.cell
def _(disease, n_loci):
    ascertainment = disease.filter(
        pl.col("gene").is_in(["KLK3", "MSMB"]) | pl.col("other_disease").str.to_lowercase().str.contains("actinic")
    )
    mo.vstack(
        [
            mo.md(
                "## 5. What these data cannot support\n\n"
                "### a. Some 'protection' is detection, not biology\n\n"
                "The strongest prostate signals sit at **`KLK3`** -- the gene encoding **PSA**. "
                "Prostate cancer in a biobank is largely ascertained *by PSA screening*. An allele "
                "that lowers PSA lowers the chance of being biopsied and therefore of being "
                "diagnosed. It may or may not lower the chance of having the disease. This data "
                "cannot distinguish those, and nothing in the API flags it."
            ),
            ascertainment.select("gene", "site", "other_disease", "beta", "mlog10p").head(8),
            mo.md(
                f"Similarly, {n_loci:,} loci were carried forward by significance alone; where a "
                "'second disease' is a precursor of the first (an allele lowering both skin cancer "
                "and **actinic keratosis**), that is one biology counted twice, not independent benefit."
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(rf"""
    ### b. The rest of the ledger

    **Ancestry.** Every resource here is European-ancestry (Finnish, UK Biobank, MVP).
    Finland is a bottlenecked population: some variants are enriched tenfold relative to
    elsewhere, which is why FinnGen finds them and also why the count does not transfer.
    This is a count of protective variants **discoverable in Europeans**, not in humans.

    **Rare protective variants are missing, and we showed why.** The frequency gradient in
    section 1 is a power artefact. A rare allele that halves cancer risk is far harder to
    detect than one that doubles it. The true rare-protective count is larger than anything
    here, and this data cannot bound it.

    **Loci are distance-defined, not LD-defined.** {LOCUS_KB}kb clustering groups leads that sit
    near each other; it does not establish that they tag one causal signal. That is why the
    PheWAS queries every protective variant and aggregates to the locus afterwards, rather than
    treating one lead as a proxy for its neighbours. The locus *count* still inherits the
    assumption, which is what the sensitivity grid is there to expose. Proper LD or
    colocalization between co-located leads would replace the assumption with a measurement.

    **Nearest gene is not causal gene.** `gene_most_severe` is a VEP annotation of the lead
    variant. *PHTF1* is a label for the *PTPN22* region, not a claim about *PHTF1*. Turning
    any of these loci into a target requires the colocalization and QTL work this notebook
    skipped (`nb01`, `nb06`).

    **Absence of a trade-off is not evidence of safety.** The "cleanly protective" loci are
    clean partly because FinnGen has fewer cases for the diseases they would affect.
    The trade-off count is a floor.

    **Effect sizes are tiny.** Most protective loci here have `abs(beta)` under 0.15, an odds
    ratio around 0.87. These are population-genetic signals, not personal risk factors.

    **What would actually answer the question better.** Multi-ancestry meta-analysis for the
    frequency floor; a case-control design not conditioned on screening for prostate and
    skin; and colocalization with pQTL/eQTL to move from locus to mechanism. In every case
    the limit is the dataset, not the query.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - Swap `C3_` for a single site (`cancer_site() == "Colorectal"`) and re-run; the
      redundancy collapse matters less and the loci become individually interpretable.
    - Feed the trade-off loci into `colocalization_by_variant` to test whether the cancer
      signal and the autoimmune signal are the *same* causal variant or two in LD. That is
      the difference between a real trade-off and a coincidence of position.
    - Keep the pQTL rows dropped in section 4 and ask which proteins the immune-tone loci
      move -- `nb06`'s `direction_consensus` is built for exactly that.
    - Use `alt_alleles` from `nb02` if you enter from an rsID rather than a phenotype; the
      credible-set index and dbSNP do not always agree on the alt allele.
    """)
    return


if __name__ == "__main__":
    app.run()
