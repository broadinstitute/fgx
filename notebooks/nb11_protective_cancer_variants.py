# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = [
#     "marimo<0.23.4",
#     "jedi<0.20.0",
#     "polars<2",
#     "httpx<0.29",
#     "altair<7",
#     "python-dotenv<2",
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
    from nb02_variant_phewas import alt_alleles  # noqa: F401
    from nb07_data_catalog import resource_metadata

    # Analysis constants used for all reported counts.
    GW_MLOG10P = 7.30103  # p <= 5e-8, the conventional genome-wide threshold
    MAX_ABS_BETA = 1.0  # |log OR| <= 1, i.e. OR in [0.37, 2.72]
    MAX_SE = 0.5
    MIN_AAF = 0.001
    LOCUS_KB = 500

    # Retry transient HTTP failures.
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

    The endpoint represents each credible set by its highest-PIP variant in the 95% set.
    This is a lead representation, not LD pruning; leads from separate sets can remain
    correlated. Returns (leads, coverage); coverage records the HTTP status per target.

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
    """Label each row with a single-linkage distance cluster on its chromosome.

    Adjacent lead positions within `kb` are grouped across phenotype definitions and
    resources. Neither the credible-set leads nor these clusters are LD-defined groups.
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
    """Apply QC filters, deduplicate rows, and attach distance-based locus labels.

    Clustering uses all QC-passing leads before filtering on effect direction. Filtering
    first could remove a lead between two protective leads and change the distance clusters.
    The primary estimate and sensitivity analysis both use this function.
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
    """Count distance-clustered protective loci at threshold `gw` and window `kb`."""
    return qc_leads(malignant, gw, kb).filter(pl.col("beta") < 0)["locus"].n_unique()


@app.function
def batch_credible_sets_by_variant(
    variants: list[str], chunk: int = 40, max_workers: int = 6, attempts: int = 4
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """POST `/credible_sets_by_variant` for multiple variants.

    A chunk that does not return HTTP 200 raises rather than returning an empty result.
    Otherwise a missing chunk could be interpreted as no other-disease association.
    Transient HTTP and transport errors are retried with exponential backoff.

    Returns `(rows, coverage)`; `coverage` is one row per chunk so the caller can assert
    full coverage before computing any verdict.

    The `variants` body field is a single string whose separator is a NEWLINE.
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

    FinnGen includes the same cancer under several case/control definitions (`_WIDE` adds Hilmo
    hospital records, `_EXALLC` excludes other cancers from controls, `_INCLAVO` adds
    primary-care records) and histological subtypes. The first matching rule is used.
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
    # nb11: Protective cancer variants

    This notebook estimates two quantities from fine-mapped cancer GWAS in GeneGenie:

    1. the number of alternate alleles associated with lower risk of at least one cancer;
    2. the number of their loci with an association to increased risk of another disease.

    Counts are reported at the association, variant, locus-site, and locus levels. Their
    interpretation depends on effect-allele coding, phenotype definitions, GWAS power, and
    the distance-based locus definition used here.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## 1. Data coverage

    The analysis queries every `C3_` or `CD2_` cancer endpoint in three resources, then
    restricts the primary estimates to malignant (`C3_`) endpoints. The endpoint metadata
    are retrieved with `resource_metadata` from `nb07`.
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
                "**Cancer endpoints defined per resource.** `finngen` is R14; `finngen_ukbb` "
                "and `finngen_mvp_ukbb` are R13 meta-analyses that include UK Biobank and the "
                "Million Veteran Program."
            ),
            catalog_summary,
        ]
    )
    return RESOURCES, catalog


@app.cell
def _(catalog):
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
                f"### Fine-mapping coverage for {len(targets):,} cancer endpoints\n\n"
                "The phenotype metadata contain more endpoints than the fine-mapping store. In "
                "particular, many FinnGen R14 endpoints return `404` from "
                "`credible_sets_by_phenotype_leads`. All estimates below use only endpoints for "
                "which the API returned credible sets."
            ),
            cov_summary,
            mo.md(f"**Total credible-set leads retrieved: {len(raw_leads):,}**"),
        ]
    )
    return (raw_leads,)


@app.cell
def _():
    mo.md(r"""
    ## 2. Effect direction and allele frequency

    GWAS effect sizes are reported per alternate allele. Here, `beta < 0` means that the
    alternate allele is associated with lower cancer risk. Reversing the effect allele
    reverses the sign, so "protective" is an allele-coded association rather than an
    intrinsic property of a locus.

    The plot shows the fraction of genome-wide significant leads with `beta < 0` by
    alternate-allele frequency.
    """)
    return


@app.cell
def _(raw_leads):
    direction_probe = (
        raw_leads.filter((pl.col("mlog10p") >= GW_MLOG10P) & pl.col("aaf").is_not_null())
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
        .agg(
            pl.len().alias("n_leads"),
            pl.col("beta").lt(0).mean().alias("frac_protective"),
        )
        .sort("alt_allele_freq")
    )
    direction_chart = (
        alt.Chart(direction_probe)
        .mark_bar(size=42, color="#2b83ba")
        .encode(
            x=alt.X(
                "alt_allele_freq:N",
                sort=["<1%", "1-5%", "5-50%", ">50% (alt is major)"],
                title="alt allele frequency",
            ),
            y=alt.Y(
                "frac_protective:Q",
                title="fraction of leads with beta < 0",
                scale=alt.Scale(domain=[0, 0.6]),
            ),
            tooltip=[
                "alt_allele_freq",
                "n_leads",
                alt.Tooltip("frac_protective:Q", format=".3f"),
            ],
        )
        .properties(
            height=260,
            title="Fraction of cancer leads with beta < 0 by alternate-allele frequency",
        )
    )
    rule = alt.Chart(pl.DataFrame({"y": [0.5]})).mark_rule(strokeDash=[6, 4], color="#d7191c").encode(y="y:Q")
    _n_significant = raw_leads.filter(pl.col("mlog10p") >= GW_MLOG10P).height
    _n_null_aaf = raw_leads.filter((pl.col("mlog10p") >= GW_MLOG10P) & pl.col("aaf").is_null()).height
    _n_below_gw = len(raw_leads) - _n_significant
    mo.vstack(
        [
            mo.ui.altair_chart(direction_chart + rule),
            direction_probe,
            mo.md(
                f"The chart includes {_n_significant - _n_null_aaf:,} genome-wide significant leads "
                f"with reported allele frequency. The {_n_below_gw:,} retrieved leads not shown are "
                f"below the p <= 5e-8 threshold; {_n_null_aaf:,} significant leads are excluded for "
                "missing allele frequency.\n\n"
                "The fraction with `beta < 0` is about 52% when the alternate allele frequency is "
                "above 50% and 14% below 1%. The discovered lead set is therefore depleted for rare "
                "alternate alleles with negative effects. These data do not estimate the number of "
                "rare protective variants that were not detected."
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(rf"""
    ## 3. Analysis definition

    | criterion | definition | purpose |
    |---|---|---|
    | malignant only | `C3_` prefix | exclude `CD2_` in-situ, benign, and uncertain-behavior endpoints |
    | significance | `mlog10p >= {GW_MLOG10P}` (p <= 5e-8) | the leads endpoint also returns sub-threshold credible sets |
    | effect estimate | `abs(beta) <= {MAX_ABS_BETA}` | exclude extreme estimates likely to be unstable |
    | standard error | `se <= {MAX_SE}` | exclude imprecise estimates |
    | allele frequency | `aaf >= {MIN_AAF}` | limit sparse estimates at very low frequency |
    | direction | `beta < 0` | alt allele lowers risk |
    | locus definition | credible-set leads, then {LOCUS_KB}kb clustering | collapse nearby leads recurring across endpoints; this does not establish LD |

    The excluded significant leads are shown below.
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
                f"### {len(implausible):,} significant leads excluded by effect-estimate filters\n\n"
                "The largest estimates occur at 16q24.3 near *MC1R* in skin-cancer endpoints. "
                "They combine very low allele frequency with extreme effect estimates or large "
                "standard errors, consistent with sparse-data separation. The associations may be "
                "real, but these estimates are not used for counting."
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
        f"After filtering, {len(qc):,} cancer associations remain: {len(protective):,} "
        f"({len(protective) / len(qc):.1%}) have `beta < 0` and "
        f"{len(qc) - len(protective):,} have `beta > 0`."
    )
    return protective, qc


@app.cell
def _():
    mo.md(r"""
    ## 4. Count by cancer site

    FinnGen includes overlapping phenotype definitions for the same cancer. Protective
    associations are therefore grouped by cancer site before loci are counted.
    """)
    return


@app.cell
def _(protective):
    _mapping_checks = {
        "Malignant neoplasm of breast": "Breast",
        "Malignant neoplasm of bronchus and lung": "Lung",
        "Malignant neoplasm": "All cancer (umbrella)",
        "Basal cell carcinoma of skin": "Skin (non-melanoma)",
    }
    assert all(cancer_site(_label) == _site for _label, _site in _mapping_checks.items())

    by_site = (
        protective.group_by("site")
        .agg(
            pl.len().alias("associations"),
            pl.col("variant").n_unique().alias("distinct_variants"),
            pl.col("locus").n_unique().alias("distance_clustered_loci"),
            pl.col("phenotype_code").n_unique().alias("finngen_codes"),
        )
        .with_columns(
            (pl.col("distinct_variants") / pl.col("distance_clustered_loci")).round(2).alias("variants_per_locus")
        )
        .sort("associations", descending=True)
    )
    skin_row = by_site.filter(pl.col("site") == "Skin (non-melanoma)")
    _prostate_row = by_site.filter(pl.col("site") == "Prostate")
    _unmapped_row = by_site.filter(pl.col("site") == "Other/unmapped")
    _all_loci = protective["locus"].n_unique()
    _skin_loci = int(skin_row["distance_clustered_loci"][0])
    _loci_without_skin = _all_loci - _skin_loci
    _multi_site_loci = (
        protective.group_by("locus")
        .agg(pl.col("site").n_unique().alias("n_sites"))
        .filter(pl.col("n_sites") > 1)
        .height
    )
    _unmapped_associations = int(_unmapped_row["associations"][0]) if _unmapped_row.height else 0
    _unmapped_loci = int(_unmapped_row["distance_clustered_loci"][0]) if _unmapped_row.height else 0
    _skin_queries_per_locus = float(skin_row["variants_per_locus"][0])
    _prostate_queries_per_locus = float(_prostate_row["variants_per_locus"][0])
    mo.vstack(
        [
            by_site,
            mo.md(
                "Non-melanoma skin cancer accounts for "
                f"{int(skin_row['associations'][0]):,} of {len(protective):,} protective associations "
                f"({int(skin_row['associations'][0]) / len(protective):.0%}) and {_skin_loci:,} of "
                f"{_all_loci:,} distance-clustered loci ({_skin_loci / _all_loci:.0%}). "
                f"The remaining {_loci_without_skin:,} loci have no non-melanoma skin association. "
                "Dermatologic ascertainment and overlapping phenotype definitions may contribute to "
                "this concentration, so the skin-excluded count is reported separately.\n\n"
                "FinnGen represents non-melanoma skin cancer with four codes: `C3_SKIN`, "
                "`C3_OTHER_SKIN_WIDE`, `C3_BASAL_CELL_CARCINOMA_WIDE`, and "
                "`C3_SQUAMOUS_CELL_CARCINOMA_SKIN_WIDE`; the three resources repeat several of these "
                "definitions.\n\n"
                f"Because {_multi_site_loci:,} loci are associated with more than one cancer site, "
                f"the site-specific counts sum to {int(by_site['distance_clustered_loci'].sum()):,} locus-site "
                f"signals rather than {_all_loci:,} loci. The mapping leaves {_unmapped_associations:,} "
                f"associations across {_unmapped_loci:,} loci as `Other/unmapped`.\n\n"
                "The number of variants queried per locus also differs by site: "
                f"{_skin_queries_per_locus:.2f} for non-melanoma skin cancer and "
                f"{_prostate_queries_per_locus:.2f} for prostate cancer. Site-specific rates of "
                "other-disease associations are therefore not compared."
            ),
        ]
    )
    return


@app.cell
def _(protective, qc):
    _all_protective_loci = protective["locus"].n_unique()
    _skin_loci = protective.filter(pl.col("site") == "Skin (non-melanoma)")["locus"].n_unique()
    _loci_without_skin = _all_protective_loci - _skin_loci
    _fully_protective_loci = (
        qc.group_by("locus").agg(pl.col("alt_protective").mean().alias("f")).filter(pl.col("f") == 1.0).height
    )
    counts = pl.DataFrame(
        {
            "unit": [
                "variant x cancer-endpoint associations",
                "distinct protective variants",
                "distinct (locus x cancer site) signals",
                f"{LOCUS_KB} kb distance-clustered loci",
                "loci with no non-melanoma skin association",
                "loci protective for every cancer they hit",
            ],
            "n": [
                len(protective),
                protective["variant"].n_unique(),
                protective.select(["locus", "site"]).unique().height,
                _all_protective_loci,
                _loci_without_skin,
                _fully_protective_loci,
            ],
        }
    )
    _mixed = (
        qc.group_by("locus")
        .agg(pl.col("alt_protective").mean().alias("f"))
        .filter((pl.col("f") > 0) & (pl.col("f") < 1))
        .height
    )
    mo.vstack(
        [
            mo.md("## 5. Number of protective cancer variants and loci\n\nCounts depend on the analysis unit:"),
            counts,
            mo.md(
                f"At the primary thresholds, the protective associations include "
                f"{_all_protective_loci:,} distance-clustered loci. {_skin_loci:,} "
                f"({_skin_loci / _all_protective_loci:.0%}) have a non-melanoma skin association; "
                f"{_loci_without_skin:,} do not. The alternate allele has `beta < 0` for every cancer "
                f"association at {_fully_protective_loci:,} loci. The remaining {_mixed:,} loci have "
                "opposite directions across cancer sites."
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
            "distance-clustered loci": protective_loci_at(malignant, gw, kb),
        }
        for gw in (GW_MLOG10P, 8.0, 9.0, 12.0)
        for kb in (100, 500, 1000)
    ]
    sens = pl.DataFrame(grid).pivot(on="locus window (kb)", index="p threshold", values="distance-clustered loci")
    grid_counts = [g["distance-clustered loci"] for g in grid]
    grid_same_kb = [g["distance-clustered loci"] for g in grid if g["locus window (kb)"] == LOCUS_KB]
    mo.vstack(
        [
            mo.md(
                "### Sensitivity to significance threshold and distance window\n\n"
                "Columns give the clustering window in kb. The p <= 5e-8, "
                f"{LOCUS_KB} kb cell is the primary estimate."
            ),
            sens,
            mo.md(
                f"Counts vary {max(grid_counts) / min(grid_counts):.1f}-fold across the full grid "
                f"({max(grid_counts):,} at the loosest corner, {min(grid_counts):,} at the tightest), and "
                f"{max(grid_same_kb) / min(grid_same_kb):.1f}-fold if you hold the {LOCUS_KB}kb window and "
                "vary only the threshold. The effect-estimate filters remove about 11% of associations "
                "but fewer than 1% of loci because most excluded rows occur at loci that remain counted."
            ),
        ]
    )
    return


@app.cell
def _(protective):
    # Query every distinct protective variant because distance clustering does not establish
    # that one lead is a proxy for another. Aggregate PheWAS results by locus afterwards.
    probes = (
        protective.sort("mlog10p", descending=True)
        .unique(subset=["variant"], keep="first", maintain_order=True)
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
        .unique(subset=["locus"], keep="first", maintain_order=True)
        .select(
            "locus",
            pl.col("gene_most_severe").alias("locus_gene"),
            pl.col("site").alias("locus_site"),
        )
    )
    n_loci = loci.height
    mo.md(
        f"### Variants tested for other-disease associations\n\n"
        f"The analysis queries all {len(probes):,} distinct protective variants across {n_loci:,} "
        "distance-clustered loci. A 500 kb window does not establish LD, so one lead cannot be used "
        "as a proxy for another. Associations are aggregated by locus only after the variant-level query."
    )
    return loci, n_loci, probes


@app.cell
def _():
    mo.md(r"""
    ## 6. Other-disease associations

    For each protective variant, `credible_sets_by_variant` returns other traits whose credible
    sets contain that variant. Effect sizes use the same alternate allele throughout. Their sign
    indicates a lower or higher coded trait value; it is interpreted as disease risk only after
    restriction to binary disease endpoints.
    """)
    return


@app.cell
def _(probes):
    phewas_raw, phewas_cov = batch_credible_sets_by_variant(probes["variant"].to_list())
    # Require complete batch coverage before classifying loci. A missing response cannot be
    # interpreted as absence of an other-disease association.
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
                f"The API returned {len(phewas_raw):,} credible-set rows for {probes.height:,} variants "
                f"in {phewas_cov.height} successful chunks. Incomplete requests raise an error rather "
                "than being treated as null results. The analysis retains genome-wide significant GWAS "
                "rows; molecular QTL rows are not used here."
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
        .sort("direction")
    )
    mo.vstack(
        [
            mo.md(
                "### Restriction to disease endpoints\n\n"
                f"All {len(non_cancer_all):,} non-cancer associations, grouped by effect direction:"
            ),
            naive_split,
            mo.md(
                "These rows include laboratory measurements, drug-purchase endpoints, and quantitative "
                "Open Targets studies. For those traits, `beta > 0` does not consistently mean harm. "
                "The locus classification therefore uses only binary disease endpoints."
            ),
        ]
    )
    return (non_cancer_all,)


@app.cell
def _(RESOURCES, non_cancer_all):
    _binary_endpoint_catalog = pl.concat(
        [
            resource_metadata(r)
            .filter(pl.col("trait_type") == "binary")
            .select(pl.lit(r).alias("resource"), "phenotype_code")
            for r in RESOURCES
        ],
        how="vertical_relaxed",
    )
    binary_codes = _binary_endpoint_catalog["phenotype_code"].unique().to_list()
    n_binary_endpoint_instances = _binary_endpoint_catalog.select("resource", "phenotype_code").unique().height
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
    disease_split = disease.group_by("direction").agg(pl.len().alias("n")).sort("direction")
    _n_higher = int(disease_split.filter(pl.col("direction") == "higher")["n"][0])
    _n_lower = int(disease_split.filter(pl.col("direction") == "also lower")["n"][0])
    mo.vstack(
        [
            mo.md(
                f"### Binary disease endpoints: {len(disease):,} associations "
                f"across {disease['trait_original'].n_unique():,} distinct disease codes"
            ),
            disease_split,
            mo.md(
                f"At the association-row level, {_n_higher:,} rows have `beta > 0` and {_n_lower:,} "
                "have `beta < 0`. Rows recur across overlapping resources, correlated disease "
                "definitions, and variants at the same locus, so they are not independent observations. "
                "The analysis unit below is the locus."
            ),
        ]
    )
    return disease, n_binary_endpoint_instances


@app.cell
def _(disease, loci, n_binary_endpoint_instances, n_loci, probes):
    # Count distinct disease codes because rows recur across resources and variants.
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

    # Recompute classifications using only the strongest protective lead at each locus.
    # `.implode()` because `is_in` against a bare Series of the same dtype is deprecated in
    # polars (pola-rs/polars#22149) -- element-wise vs collection membership is ambiguous
    # there. Imploding to a single list cell says "membership in this collection" explicitly.
    rep_variants = probes.sort("cancer_mlog10p", descending=True).unique(
        subset=["locus"], keep="first", maintain_order=True
    )["variant"]
    rep_per_locus = verdicts(disease.filter(pl.col("variant").is_in(rep_variants.implode())))
    rep_tradeoff = rep_per_locus.filter(pl.col("n_worse") > 0).height
    rep_silent = n_loci - rep_per_locus.height
    verdict = pl.DataFrame(
        {
            "verdict": [
                "no significant other-disease association",
                "disease associations detected, none with beta > 0",
                "at least one disease association with beta > 0",
            ],
            "loci": [silent, len(clean_loci), len(tradeoff_loci)],
        }
    )
    _expected_null_hits = probes.height * n_binary_endpoint_instances * 5e-8
    _tradeoff_denominator = len(clean_loci) + len(tradeoff_loci)
    mo.vstack(
        [
            mo.md(f"## 7. Locus-level other-disease results\n\nOf the {n_loci:,} protective cancer loci:"),
            verdict,
            mo.md(
                f"Among loci with a significant other-disease association, {len(tradeoff_loci):,} "
                f"have at least one association with `beta > 0` and {len(clean_loci):,} have only "
                f"associations with `beta < 0` ({len(tradeoff_loci) / _tradeoff_denominator:.0%} with "
                "`beta > 0`). Each locus is classified from the union of results for all protective "
                "variants in that locus.\n\n"
                "This percentage is conditional on the variants queried, available disease GWAS, and "
                "endpoint power. Loci with more queried variants have more opportunities to acquire a "
                "`beta > 0` association, and a null result is weaker evidence than a detected association.\n\n"
                f"A conservative scale calculation treats all {probes.height:,} variants x "
                f"{n_binary_endpoint_instances:,} binary resource-endpoints as tests. At 5e-8 this gives "
                f"{_expected_null_hits:.2f} expected threshold crossings. This calculation does not "
                "address phenotype dependence, selection, winner's curse, or unequal query density."
            ),
            mo.md(
                f"Using only the strongest cancer lead at each locus gives {rep_tradeoff:,} loci with "
                f"a `beta > 0` disease association, compared with {len(tradeoff_loci):,} when all "
                f"{probes.height:,} protective variants are queried. The corresponding number with no "
                f"other-disease association is {rep_silent:,} rather than {silent:,}. A representative "
                "variant therefore misses associations carried by other leads in the same distance cluster."
            ),
        ]
    )
    return (tradeoff_loci,)


@app.cell
def _(disease, tradeoff_loci):
    worse = (
        disease.filter(pl.col("direction") == "higher")
        .join(
            tradeoff_loci.select("locus", "locus_gene", "locus_site"),
            on="locus",
            how="inner",
        )
        .sort("mlog10p", descending=True)
        # One row per (locus, disease), keeping the strongest. Deduping on gene instead would
        # merge two genuinely distinct loci that share a nearest-gene annotation.
        .unique(subset=["locus", "trait_original"], keep="first", maintain_order=True)
        .sort("mlog10p", descending=True)
    )
    mo.vstack(
        [
            mo.md(
                "### Disease-increasing associations by locus\n\n"
                "The table contains one row per locus and disease code, retaining the strongest "
                "association. `gene` is the VEP nearest/most-severe annotation, not an assigned causal gene."
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
        .sort(["n", "gene"], descending=[True, False])
        .head(15)
    )
    tradeoff_chart = (
        alt.Chart(top_genes)
        .mark_bar(color="#d7191c", opacity=0.85)
        .encode(
            x=alt.X("n:Q", title="disease codes with beta > 0"),
            y=alt.Y("gene:N", sort="-x", title=None),
            tooltip=["gene", "n", "association_rows", alt.Tooltip("max_mlog10p:Q", format=".1f")],
        )
        .properties(height=380, title="Cancer-protective loci with disease associations at beta > 0")
    )
    mo.vstack(
        [
            mo.ui.altair_chart(tradeoff_chart),
            mo.md(
                "`n` is an upper bound because related disease codes are counted separately. For "
                "example, the *PHTF1*-labelled locus includes strict autoimmune hypothyroidism, "
                "hypothyroidism with drug reimbursement, and disorders of the thyroid gland as three "
                "codes. No disease-family mapping is applied to the non-cancer endpoints."
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ### Selected loci

    Several high-count loci have associations consistent with established biology. Positional
    overlap within a 500 kb cluster does not show that the cancer and other-disease associations
    share a causal variant.

    - The *PTPN22* region, labelled `PHTF1` by VEP, includes lower skin-cancer risk and higher
      risk of autoimmune hypothyroidism, rheumatoid arthritis, type 1 diabetes, and Crohn's disease.
    - The 9p21 `CDKN2B-AS1` cluster includes protective associations across several cancer sites
      and associations with coronary atherosclerosis, ischaemic heart disease, and glaucoma.
    - The `PTCSC2` region near *FOXE1* includes lower thyroid-cancer risk and higher risk of
      hypothyroidism and goitre.
    - `IL2RA`, `BACH2`, `CDKAL1`, `RAB5B`, `UBE2L3`, and `IRF5` also pair lower cancer risk with
      higher risk of one or more autoimmune diseases.

    Several of the highest-count loci are annotated near immune-regulation genes. The correlated
    endpoints and distance-defined loci do not establish a shared mechanism.

    A Polish case-control study of four *CHEK2* founder alleles reported lower lung-cancer odds
    among carriers (895 cases, 6,391 controls; OR 0.3, 95% CI 0.2-0.5, p = 3e-8;
    [PMID 18281249](https://pubmed.ncbi.nlm.nih.gov/18281249/)). This supports the cancer-effect
    direction but does not show that the cancer and myeloproliferative associations share a
    causal signal.
    """)
    return


@app.cell
def _(disease, n_loci, protective):
    psa_signals = (
        protective.filter(pl.col("gene_most_severe").is_in(["KLK3", "MSMB"]))
        .sort("mlog10p", descending=True)
        .unique(
            subset=["gene_most_severe", "variant"],
            keep="first",
            maintain_order=True,
        )
    )
    actinic_signals = (
        disease.filter(pl.col("other_disease").str.to_lowercase().str.contains("actinic"))
        .sort("mlog10p", descending=True)
        .unique(subset=["variant", "trait_original"], keep="first", maintain_order=True)
    )
    apoe_cancer_signals = (
        protective.filter(pl.col("gene_most_severe").is_in(["APOE", "AC011481.3"]))
        .sort("mlog10p", descending=True)
        .unique(
            subset=["gene_most_severe", "variant"],
            keep="first",
            maintain_order=True,
        )
    )
    mo.vstack(
        [
            mo.md(
                "## 8. Interpretation limits\n\n"
                "### Ascertainment and selection\n\n"
                "The strongest prostate signals include `KLK3`, which encodes PSA, and `MSMB`. "
                "An allele that lowers PSA can reduce biopsy and diagnosis without reducing disease "
                "incidence. These data do not distinguish those pathways."
            ),
            mo.md("Cancer associations at the PSA-related loci:"),
            psa_signals.select(
                pl.col("gene_most_severe").alias("gene"),
                "site",
                pl.col("trait_label").alias("cancer_trait"),
                "variant",
                pl.col("beta").alias("cancer_beta"),
                pl.col("mlog10p").alias("cancer_mlog10p"),
            ),
            mo.md(
                f"Among the {n_loci:,} loci, several alleles are associated with lower risk of both "
                "skin cancer and actinic keratosis. These correlated endpoints do not represent "
                "independent protective effects."
            ),
            actinic_signals.select("gene", "site", "other_disease", "variant", "beta", "mlog10p").head(8),
            mo.md(
                "The *APOE* liver-cancer association may also reflect selection. Simulations in a "
                "case-only stroke GWAS showed that an *APOE* age-at-onset signal could arise from survival "
                "bias ([PMID 39286457](https://pubmed.ncbi.nlm.nih.gov/39286457/)). A 38-study meta-analysis "
                "found no major overall association between common *APOE* genotype and cancer, while "
                "allowing small site-specific effects "
                "([PMID 41197273](https://pubmed.ncbi.nlm.nih.gov/41197273/)). This notebook has no "
                "survival model and cannot distinguish competing mortality from a hepatic effect."
            ),
            apoe_cancer_signals.select(
                pl.col("gene_most_severe").alias("gene"),
                "site",
                pl.col("trait_label").alias("cancer_trait"),
                "variant",
                pl.col("beta").alias("cancer_beta"),
                pl.col("mlog10p").alias("cancer_mlog10p"),
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(rf"""
    ### Generalizability and interpretation

    **Population.** FinnGen is a Finnish founder-population resource, and the combined resources
    do not provide a balanced multi-ancestry estimate. The result is not a human-wide count.

    **Rare variants.** Only 14% of significant leads below 1% alternate-allele frequency have
    `beta < 0`. The analysis does not bound the number of rare protective variants that were
    not detected.

    **Locus definition.** The {LOCUS_KB} kb clusters are based on physical distance, not LD.
    The sensitivity analysis varies this window, but LD-aware grouping or colocalization is
    required to establish whether nearby leads represent the same signal.

    **Gene labels.** `gene_most_severe` is a VEP annotation of the lead variant, not a causal-gene
    assignment. For example, `PHTF1` labels the *PTPN22* region.

    **Null results.** No detected disease-increasing association does not imply safety. Querying
    one representative variant per locus identifies 76 loci with such an association; querying
    all protective variants identifies 111. Detection depends on query density, endpoint coverage,
    and power.

    **Chromosome labels.** The API uses `chr == 23` for chromosome X. Any extension that joins a
    source labelled `X` must normalize this field.

    **Effect size.** About 80% of protective association rows have `abs(beta) < 0.15`, corresponding
    to an odds ratio near 0.87. These are association estimates, not individual risk predictions.

    Multi-ancestry analysis, designs less sensitive to screening-related ascertainment, and
    colocalization with eQTL or pQTL data would address these limitations.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - Repeat the analysis within one cancer site to reduce phenotype overlap.
    - Test whether cancer and other-disease associations colocalize rather than relying on
      physical proximity.
    - Use the excluded pQTL rows with `nb06` to examine protein-level direction of effect.
    """)
    return


if __name__ == "__main__":
    app.run()
