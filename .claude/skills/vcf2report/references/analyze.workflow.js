export const meta = {
  name: 'vcf2report-analyze',
  description: 'Analyze one exome as 8 visible stages — each labeled LOCAL (your machine) or CLAUDE',
  // Phase titles carry who does the work: 🖥️ = local (vcf2report on the machine, deterministic,
  // offline) · 🤖 = Claude (reasoning). They render as the step boxes in Background Tasks.
  phases: [
    { title: '🖥️ Local · Dependency check', detail: 'probe python + tools + local stores, explain what each enables' },
    { title: '🖥️ Local · Inspect VCF', detail: 'build, sample, counts, is it annotated (VEP/SnpEff/consequence)' },
    { title: '🖥️ Local · Annotate', detail: 'annotate a raw VCF (bcftools+SnpEff+vcfanno) or explain coordinate-only limits' },
    { title: '🤖 Claude · Analysis triage', detail: 'which ACMG criteria are computable from this VCF + stores' },
    { title: '🤖 Claude · Phenotype → HPO', detail: 'map free-text phenotype to HPO terms' },
    { title: '🖥️ Local · Prioritize (gnomAD+AM+ClinVar+HPO)', detail: 'run the engine, rank candidates against the population background' },
    { title: '🖥️ Local · QC', detail: 'funnel, sequencing quality, gnomAD-store safety net' },
    { title: '🤖 Claude · Laudo', detail: 'synthesize the auditable report from the facts' },
  ],
}

// args (from the vcf2report skill): { repo, vcf, sample, hpo?, phenotypeText?, out?, lift?,
// chain?, reference? }. The skill fills these after locating the repo and the VCF.
const a = (typeof args === 'string' ? JSON.parse(args) : args) || {}
const REPO = a.repo || '.'
const SAMPLE = a.sample || 'sample'
const OUT = a.out || `out/${SAMPLE}`

if (!a.vcf) {
  log('No VCF path in args — nothing to analyze.')
  return { error: 'no vcf' }
}
let vcf = a.vcf
let hpoFile = a.hpo || ''
const REF = a.reference || ''

// 1 — 🖥️ DEPENDENCY CHECK — opens Background Tasks; the Parquet-store GATE runs here, before
// anything patient-specific. The analysis proceeds ONLY if the stores are available + intact.
phase('🖥️ Local · Dependency check')
const deps = await agent(
  `cd ${REPO} && python3 scripts/preflight.py. Return, in plain language: the python version and which of ` +
  `bcftools/snpEff/vcfanno are on PATH. (Parquet-store details come from the gate step next.)`,
  { label: '🖥️ preflight', phase: '🖥️ Local · Dependency check' })
// HARD GATE — availability + version + build date + integrity of the 3 Parquet stores.
const gate = await agent(
  `cd ${REPO} && python3 scripts/check_stores.py --gate; echo "STORES_EXIT=$?". Show the user the store ` +
  `table VERBATIM — for EACH parquet (gnomAD · AlphaMissense · ClinVar): available? · version (source) · ` +
  `BUILD DATE · integrity · complete — then the READY/BLOCKED banner and the STORES_EXIT line (do not omit it).`,
  { label: '🖥️ parquet stores — availability + integrity GATE', phase: '🖥️ Local · Dependency check' })
if (!/STORES_EXIT=0/.test(gate || '')) {
  log('⛔ ANALYSIS BLOCKED — the gnomAD / AlphaMissense / ClinVar Parquet stores are not all available and ' +
      'intact. Fix the flagged store(s) — build_gnomad_parquet.py / build_alphamissense_parquet.py / ' +
      'build_clinvar_parquet.py (or stamp_store_manifest.py) — then re-run. No analysis was performed.')
  return { error: 'stores_unavailable', deps, gate }
}
log('✅ Parquet stores available + intact — proceeding with the analysis.')

// 2 — 🖥️ INSPECT VCF — liftover first if needed, then detect build + whether it is annotated
phase('🖥️ Local · Inspect VCF')
if (a.lift) {
  const chain = a.chain ? `--chain ${a.chain}` : ''
  const lifted = `${OUT}/${SAMPLE}.hg38.vcf`
  await agent(
    `cd ${REPO} && mkdir -p ${OUT} && VCF2REPORT_ALLOW_NETWORK=1 python3 ` +
    `scripts/liftover_to_grch38.py ${vcf} ${lifted} ${chain} 2>&1 | tail -3. ` +
    `Return the [liftover] summary line (records/lifted/skipped).`,
    { label: '🖥️ liftover 37→38', phase: '🖥️ Local · Inspect VCF' })
  vcf = lifted
}
const hpoFlag = (hpoFile || a.phenotypeText) ? '--hpo' : ''
const inspect = await agent(
  `cd ${REPO} && python3 scripts/inspect_vcf.py ${vcf} ${hpoFlag}. From the JSON, output the VERY FIRST line as ` +
  `exactly \`ANNOTATED=yes\` or \`ANNOTATED=no\` (from inspect.annotated), then summarize: build, sample id, total + ` +
  `PASS counts, annotation source (VEP CSQ / SnpEff ANN / consequence / none), and the capabilities.criteria map ` +
  `(each ACMG criterion → available|limited|na + reason). Do NOT classify yet.`,
  { label: '🖥️ inspect + capabilities', phase: '🖥️ Local · Inspect VCF' })
const annotated = /ANNOTATED=yes/i.test(inspect || '')

// 3 — 🖥️ ANNOTATE — only if the VCF isn't annotated; otherwise skip (visibly). No reference → explain the limits.
phase('🖥️ Local · Annotate')
if (annotated) {
  log('🖥️ VCF already annotated — skipping annotation (consequence terms present).')
} else if (REF) {
  await agent(
    `cd ${REPO} && bash scripts/annotate_vcf.sh ${vcf} ${REF} ${OUT}/${SAMPLE}.annotated.vcf.gz 2>&1 | tail -8. ` +
    `Return whether annotation succeeded and the steps run (bcftools norm + SnpEff + vcfanno).`,
    { label: '🖥️ annotate (SnpEff+vcfanno)', phase: '🖥️ Local · Annotate' })
  vcf = `${OUT}/${SAMPLE}.annotated.vcf.gz`
} else {
  await agent(
    `The VCF is NOT annotated and no reference FASTA was provided, so local annotation (bcftools+SnpEff+vcfanno) ` +
    `can't run. In 2-3 lines tell the user: classification will be COORDINATE-ONLY — PVS1/PM4/PP3/BP4 and HGVS ` +
    `c./p. are unavailable; gnomAD/ClinVar coordinate lookups + the ≥2★ ClinVar safety flag still work. ` +
    `Recommend annotating first (docs/ANNOTATION.md) for full ACMG.`,
    { label: '🤖 annotation options', phase: '🖥️ Local · Annotate' })
}

// 4 — 🤖 ANALYSIS TRIAGE — the honesty gate: say what this run can and cannot conclude
phase('🤖 Claude · Analysis triage')
const triage = await agent(
  `From the inspection + capabilities already gathered for ${SAMPLE} (annotated=${annotated}, ` +
  `reference=${REF ? 'provided' : 'none'}), give a short capability list: which ACMG criteria this run CAN ` +
  `evaluate vs limited/NA, each with a one-line consequence for the laudo (e.g. 'no gnomAD store → PM2 disabled ` +
  `→ over-call risk', 'single-proband → segregation N/A'). This gate is stated before running, not after.`,
  { label: '🤖 what can we conclude', phase: '🤖 Claude · Analysis triage' })

// 5 — 🤖 PHENOTYPE → HPO — turn the clinician's free-text phenotype into HPO ids (only if no file given)
phase('🤖 Claude · Phenotype → HPO')
if (!hpoFile && a.phenotypeText) {
  const mapped = await agent(
    `The patient's phenotype, free text: "${a.phenotypeText}". Map it to Human Phenotype ` +
    `Ontology (HPO) ids. Return ONLY the HP: ids, comma-separated (e.g. HP:0001250,HP:0001263). ` +
    `Be precise and conservative — only terms clearly implied by the text.`,
    { label: '🤖 map phenotype→HPO', phase: '🤖 Claude · Phenotype → HPO' })
  const ids = (mapped || '').match(/HP:\d+/gi) || []
  if (ids.length) {
    hpoFile = `${OUT}/hpo.txt`
    await agent(
      `cd ${REPO} && mkdir -p ${OUT} && printf '${ids.join('\\n')}\\n' > ${hpoFile} && ` +
      `echo wrote ${ids.length} HPO terms. Return the count.`,
      { label: '🖥️ write hpo.txt', phase: '🤖 Claude · Phenotype → HPO' })
  }
} else {
  log(hpoFile ? '🖥️ HPO terms supplied as a file.' : '🤖 No phenotype given — running genotype-only.')
}
const HPO = hpoFile ? `--hpo ${hpoFile}` : ''

// 6 — 🖥️ PRIORITIZE — the deterministic engine over gnomAD + AlphaMissense + ClinVar + HPO
phase('🖥️ Local · Prioritize (gnomAD+AM+ClinVar+HPO)')
const classify = await agent(
  `cd ${REPO} && python3 scripts/run_headless.py ${vcf} ${HPO} --sample-id ${SAMPLE} ` +
  `--out ${OUT} --timing 2>&1 | tail -40. Return: the RANKED candidate list (gene → tier, the phenotype- and ` +
  `tier-topped rows first), the funnel counts, and the per-stage timings. This runs the engine over gnomAD + ` +
  `AlphaMissense + ClinVar + HPO — all local + deterministic, no network, no LLM.`,
  { label: '🖥️ run engine (gnomAD+AM+ClinVar+HPO)', phase: '🖥️ Local · Prioritize (gnomAD+AM+ClinVar+HPO)' })

// 7 — 🖥️ QC — surface the funnel + sequencing quality + the gnomAD-store safety net as its own gate
phase('🖥️ Local · QC')
const qc = await agent(
  `Read ${REPO}/${OUT}/${SAMPLE}_report.md and return the QC section only: the funnel (total → PASS → ` +
  `QC-passing → candidates), the Sequencing-quality panel numbers (depth/GQ, Ti/Tv, het:hom, indel:SNV, ` +
  `multiallelic, novelty), and ANY gnomAD-store / coverage warning. This is the QC gate before the laudo.`,
  { label: '🖥️ QC funnel + safety net', phase: '🖥️ Local · QC' })

// 8 — 🤖 LAUDO — read the machine's facts and synthesize the auditable narrative
phase('🤖 Claude · Laudo')
const report = await agent(
  `Read ${REPO}/${OUT}/${SAMPLE}_report.md and return, as compact structured text for ` +
  `rendering a laudo: the Conclusion bullets; the Sequencing-quality numbers; the primary ` +
  `and secondary findings rows (gene, transcript, HGVS, consequence, tier, ACMG criteria); ` +
  `and the Performance/timing lines. Do NOT re-classify — report the engine's calls faithfully.`,
  { label: '🤖 synthesize laudo', phase: '🤖 Claude · Laudo' })

return { sample: SAMPLE, out: OUT, hpo: hpoFile, annotated, deps, inspect, triage, classify, qc, report }
