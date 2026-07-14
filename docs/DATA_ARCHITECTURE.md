# Data architecture — gnomAD · AlphaMissense · ClinVar

How the three frequency/pathogenicity/clinical databases are stored and joined, and
why they are kept as **separate per-source Parquet stores** rather than one unified table.

## Q1 — Is the data accessed through MCP? No.

There are two independent layers:

```
Claude (Desktop / Code)  ──MCP or Bash──▶  vcf2report engine  ──direct──▶  local data
                                            (Python package)               gnomAD Parquet (DuckDB)
                                                                           AlphaMissense (tabix / Parquet)
                                                                           ClinVar (tabix / Parquet)
```

`mcp_server.py` is a **thin FastMCP adapter** that wraps the pipeline/lookup functions as
tools; `scripts/*.py` are the equivalent Bash adapters. **No variant data ever flows through
MCP** — the engine opens DuckDB / pysam directly on the machine. MCP is only the *Claude ↔ engine*
interface, never *engine ↔ data*.

## Q2 — One unified Parquet, or separate? Separate, joined by DuckDB.

Keep **three per-source Parquet stores**, each Hive-partitioned by `chrom` and locus-sorted,
and join them at query time with **one DuckDB `LEFT JOIN`** over a temp table of candidate
variants. Do **not** materialise a unified wide table.

| factor | unified Parquet | **separate + DuckDB join (chosen)** |
|---|---|---|
| Update cadence | ClinVar changes **weekly** → rebuild the whole 1.3 GB store every week | each source rebuilt independently; ClinVar weekly touches only its 69 MB store |
| Key-space | gnomAD ~70M (all classes), AM ~71M (missense-SNV only), ClinVar ~4M → a wide outer join is **mostly NULL** | each `LEFT JOIN` pulls only the columns a source has for the ~10²–10⁴ candidate loci |
| Coverage semantics | blurs three different "absence" rules into one row | each source keeps its own presence flag + coverage rule (invariant stays per-source) |
| Query speed | one scan | predicate pushdown + row-group min/max + locus-sort prune to a few row groups per partition — **no measurable penalty** |
| Blast radius | a failed rebuild can corrupt gnomAD frequencies | a ClinVar rebuild failure can't touch gnomAD/AM |

### Per-store schema (all sorted by `pos, ref, alt`, partitioned by `chrom`, chr-prefixed)

- **gnomAD** (unchanged, already Parquet): `chrom, pos, ref, alt, filter, af, af_grpmax, ac, an,
  nhomalt, faf95, grpmax_pop, af_<pop>…` + sidecar `_meta.json` (coverage mode: full/partial/bed).
- **AlphaMissense** (`scripts/build_alphamissense_parquet.py`, MAX-per-locus so the join is 1:1):
  `chrom, pos, ref, alt, am_pathogenicity, am_class`. Missense-SNV only. No coverage sidecar —
  AM absence is *never* asserted (no prediction → None, fires nothing).
- **ClinVar** (`scripts/build_clinvar_parquet.py`, weekly): `chrom, pos, ref, alt, significance,
  review_status, review_stars (0–4, precomputed for the ≥2★ safety flag), accession, condition`.
  No coverage sidecar — a missing row is simply "not in ClinVar".

### The batch-annotate query (one pass, invariant-preserving)

```sql
SELECT q.chrom, q.pos, q.key,
       g.pos AS g_pos, g.filter AS g_filter, g.af, g.af_grpmax, g.faf95, g.grpmax_pop, g.nhomalt,
       a.am_pathogenicity, a.am_class,
       c.significance AS clinvar_sig, c.review_status, c.review_stars, c.accession, c.condition
FROM q
LEFT JOIN read_parquet('data/gnomad/gnomad_parquet/**/*.parquet')  g
       ON g.chrom=q.chrom AND g.pos=q.pos AND upper(g.ref)=q.ref AND upper(g.alt)=q.alt
LEFT JOIN read_parquet('data/alphamissense/am_parquet/**/*.parquet') a
       ON a.chrom=q.chrom AND a.pos=q.pos AND upper(a.ref)=q.ref AND upper(a.alt)=q.alt
LEFT JOIN read_parquet('data/clinvar/clinvar_parquet/**/*.parquet') c
       ON c.chrom=q.chrom AND c.pos=q.pos AND upper(c.ref)=q.ref AND upper(c.alt)=q.alt;
```

**Never fabricate a false absence** (the core invariant), enforced *per source* after the join —
`filter` is `SELECT`ed, never a join predicate:
- gnomAD `g_pos` NOT NULL + `g_filter` ≠ PASS → present-but-filtered (AF unavailable, **no PM2**);
  `g_pos` NULL → candidate absence, asserted `af=0.0` **only** where `_meta` mode = full/bed vouches
  for the locus, else left None.
- AM `am_pathogenicity` NULL → no missense prediction → fires nothing (never an assertion).
- ClinVar `clinvar_sig` NULL → "no record" (never a benign call).

## Build / refresh

| store | script | cadence | size |
|---|---|---|---|
| gnomAD | `scripts/build_gnomad_parquet.py` | frozen (v4.1) | ~1.3 GB |
| AlphaMissense | `scripts/build_alphamissense_parquet.py` | frozen | ~0.4–0.6 GB |
| ClinVar | `scripts/build_clinvar_parquet.py` | **weekly** | ~69 MB |

Each builder writes to a `.building` temp dir and swaps atomically, so a mid-rebuild reader
never sees a partial store. All three are git-ignored (data, not code).

## Performance (measured on this machine)

The 3-way join was validated on real NA12878 candidate loci (POGZ spike + an AlphaMissense-only
missense + a ≥2★ ClinVar variant), and the invariant held (AM `NULL` for the stop-gain — absence
never asserted; gnomAD `NULL` for the AM-only variant — genuine absence; ClinVar `NULL` = no record):

| query | latency |
|---|---|
| naive `read_parquet('**/*.parquet')` LEFT JOIN, no pruning | **~29 s** (scans every partition) |
| **chr-pruned** (`WHERE chrom IN (candidate chroms)`) | **321 ms cold / 162 ms warm** |

**The production wiring must prune each store to the candidate chromosomes** (as gnomAD's `prime()`
effectively already does) — then the whole 3-way annotate is sub-second. Store sizes as built:
ClinVar 69 MB (4,195,020 rows), AlphaMissense 509 MB (71,034,269 loci, 0 duplicates after MAX-agg).

## Status

- ✅ Build scripts landed and validated (ClinVar 4.2M variants in ~5 s; AlphaMissense 71M loci,
  MAX-per-locus, 0 duplicates). 3-way join proven correct + fast (chr-pruned) on real loci.
- ✅ **Wired into the annotate stage.** `annotate/alphamissense_parquet.py` and
  `annotate/clinvar_parquet.py` each run a chr-pruned `prime()` mirroring `gnomad_parquet`;
  `alphamissense.prime()` prefers the Parquet store (tabix fallback preserved), `clinvar.lookup()`
  reads the Parquet cache before the per-variant tabix (a miss falls through unchanged), and the
  pipeline primes ClinVar over the whole post-QC set — its **first-ever batch path**. AlphaMissense
  is stored as **DOUBLE** so scores are bit-for-bit equal to the tabix `float()` parse.
- ✅ **Validated byte-identical** on real NA12878 (2,395 candidates): parquet path == tabix path on
  every gene / tier / ACMG criterion / AlphaMissense score / ClinVar field (0 differences). ~12%
  faster offline (7.98 s → 7.05 s); ClinVar resolves in one 0.24 s join instead of ~24k per-variant
  lookups (a larger win when network ClinVar would otherwise be hit). The tabix clients remain the
  offline / no-duckdb fallback and are not retired.

> ⚠️ **AlphaMissense license:** the source file header declares **CC BY-NC-SA 4.0** (non-commercial,
> share-alike), while `scripts/fetch_alphamissense.sh` states CC BY 4.0. Resolve which is authoritative
> before redistributing any derived store. The `am_parquet` store is git-ignored regardless.
