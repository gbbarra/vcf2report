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
  `exactly \`ANNOTATED=yes\` or \`ANNOTATED=no\` (from inspect.annotated), then show the user, clearly:\n` +
  `• build, sample id, total + PASS counts;\n` +
  `• **is it annotated?** — annotation_source (VEP CSQ / SnpEff ANN / consequence / none) AND the coverage: ` +
  `\`variants_with_consequence\`/\`total_variants\` (consequence_coverage). A raw caller's VCF is ~0% — say so plainly;\n` +
  `• **WHICH FIELDS the VCF actually carries** — list \`info_fields\` (e.g. AC/AF/DP/MQ… from the caller) and ` +
  `highlight \`info_fields_annotation\` (ANN/CSQ/CLNSIG/gnomad_*…) if any;\n` +
  `• the capabilities.criteria map (each ACMG criterion → available|limited|na + reason).\n` +
  `If coverage is ~0 the file is NOT annotated: PVS1/PM4/PP3 and HGVS are unavailable for the whole callset — ` +
  `flag that Stage 4 must annotate it, because the laudo depends on it. Do NOT classify yet.`,
  { label: '🖥️ inspect: anotado? quais campos?', phase: '🖥️ Local · Inspect VCF' })
const annotated = /ANNOTATED=yes/i.test(inspect || '')

// 3 — 🖥️ ANNOTATE — only if the VCF isn't annotated; otherwise skip (visibly). A reference FASTA is
// OPTIONAL (it only adds indel left-alignment), so a raw VCF gets annotated whenever SnpEff is
// installed. Annotation is what gives the laudo its gene/consequence/HGVS — only if it genuinely
// cannot run do we fall back to explaining the coordinate-only limits.
phase('🖥️ Local · Annotate')
let isAnnotated = annotated   // updated below when Stage 3 annotates, so the triage reflects reality
if (annotated) {
  log('🖥️ VCF already annotated — skipping annotation (consequence terms present).')
} else {
  const out = `${OUT}/${SAMPLE}.annotated.vcf.gz`
  const ann = await agent(
    `cd ${REPO} && mkdir -p ${OUT} && bash scripts/annotate_vcf.sh ${vcf} ${out} ${REF} > ${OUT}/annotate.log 2>&1; ` +
    `echo "ANNOTATE_EXIT=$?"; tail -12 ${OUT}/annotate.log. ` +
    `Output the ANNOTATE_EXIT line VERBATIM as your FIRST line, then report to the user: whether annotation ` +
    `succeeded, and the "annotated N/M records" count (that is the share of the callset that now has gene + ` +
    `consequence + HGVS). If ANNOTATE_EXIT is non-zero, say which step failed and quote the error.`,
    { label: '🖥️ annotate (bcftools norm + SnpEff MANE)', phase: '🖥️ Local · Annotate' })
  if (/ANNOTATE_EXIT=0/.test(ann || '')) {
    vcf = out
    isAnnotated = true
    log('🖥️ Annotated — the laudo will carry gene, consequence and HGVS c./p.')
  } else {
    await agent(
      `Local annotation did NOT run (SnpEff missing, or the script failed). In 2-3 lines tell the user: ` +
      `classification will be COORDINATE-ONLY — PVS1/PM4/PP3/BP4 and HGVS c./p. are unavailable for the whole ` +
      `callset; gnomAD/ClinVar coordinate lookups + the ≥2★ ClinVar safety flag still work. To fix it, install ` +
      `the annotator: \`bash scripts/setup_snpeff.sh\` (see docs/ANNOTATION.md).`,
      { label: '🤖 annotation unavailable — limits', phase: '🖥️ Local · Annotate' })
  }
}

// 4 — 🤖 ANALYSIS TRIAGE — the honesty gate: say what this run can and cannot conclude
phase('🤖 Claude · Analysis triage')
const triage = await agent(
  `State the capability gate for ${SAMPLE} as it stands NOW, after annotation: the VCF ` +
  `${isAnnotated ? 'IS functionally annotated (SnpEff MANE ran in Stage 4) — so PVS1/PM4/PP3/BP4 and HGVS ARE ' +
    'evaluable' : 'is NOT annotated — PVS1/PM4/PP3/BP4 and HGVS are unavailable (coordinate-only)'}. ` +
  `The three required Parquet stores (gnomAD · AlphaMissense · ClinVar) PASSED the Stage-1 gate, so PM2/BA1/BS1, ` +
  `PP3/BP4 and PS1/PM5/PP5 ARE available — do NOT say 'no stores'. (A reference FASTA is optional and only adds ` +
  `indel left-alignment; its absence does not disable any criterion.) Give a short list: which ACMG criteria ` +
  `this run CAN evaluate vs limited/NA, each with a one-line consequence (e.g. 'single-proband → PS2/PM3 N/A'). ` +
  `Be consistent with what the engine will actually do next — do not claim '0 criteria evaluable' when annotation ` +
  `and the stores are both present.`,
  { label: '🤖 what can we conclude', phase: '🤖 Claude · Analysis triage' })

// 5 — 🤖 PHENOTYPE → HPO — turn the clinician's free-text phenotype into HPO ids (only if no file given).
// Claude MAPS the text to ids here; the FILE is written by the deterministic prioritize command below
// (Stage 6), never by a separate agent — delegating a file write to an LLM lost the file silently before.
phase('🤖 Claude · Phenotype → HPO')
let hpoIds = []
if (!hpoFile && a.phenotypeText) {
  const mapped = await agent(
    `The patient's phenotype, free text: "${a.phenotypeText}". Map it to Human Phenotype ` +
    `Ontology (HPO) ids. Return ONLY the HP: ids, comma-separated (e.g. HP:0001250,HP:0001263). ` +
    `Be precise and conservative — only terms clearly implied by the text.`,
    { label: '🤖 map phenotype→HPO', phase: '🤖 Claude · Phenotype → HPO' })
  hpoIds = (mapped || '').match(/HP:\d+/gi) || []
  if (hpoIds.length) { hpoFile = `${OUT}/hpo.txt`; log(`🤖 Mapped the phenotype to ${hpoIds.length} HPO term(s).`) }
  else { log('⚠️ The phenotype text did not map to any HPO id — the run will be genotype-only (no PP4).') }
} else {
  log(hpoFile ? '🖥️ HPO terms supplied as a file.' : '🤖 No phenotype given — running genotype-only.')
}
const HPO = hpoFile ? `--hpo ${hpoFile}` : ''
// Write the HPO file as part of the SAME command that runs the engine, and `test -s` it, so the
// phenotype provably reaches the classifier (or the command fails loudly) rather than being dropped.
const writeHpo = hpoIds.length ? `mkdir -p ${OUT} && printf '${hpoIds.join('\\n')}\\n' > ${hpoFile} && test -s ${hpoFile} && ` : ''

// 6 — 🖥️ PRIORITIZE — the deterministic engine over gnomAD + AlphaMissense + ClinVar + HPO
phase('🖥️ Local · Prioritize (gnomAD+AM+ClinVar+HPO)')
const classify = await agent(
  `cd ${REPO} && ${writeHpo}python3 scripts/run_headless.py ${vcf} ${HPO} --sample-id ${SAMPLE} ` +
  `--out ${OUT} --timing 2>&1 | tail -40. Return: the RANKED candidate list (gene → tier, the phenotype- and ` +
  `tier-topped rows first), the funnel counts, and the per-stage timings. This runs the engine over gnomAD + ` +
  `AlphaMissense + ClinVar + HPO — all local + deterministic, no network, no LLM. ` +
  (a.phenotypeText
    ? `A phenotype WAS provided, so the run must show phenotype terms loaded and phenotype-matched ranking — if ` +
      `it instead reports "no phenotype terms / none provided", FLAG that loudly as a failure: the phenotype did ` +
      `not reach the engine.`
    : `No phenotype was provided — this is a genotype-only run.`),
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
  `rendering a laudo: the Conclusion bullets; the Sequencing-quality numbers; and the findings rows ` +
  `(gene, transcript, HGVS, consequence, tier, ACMG criteria) grouped by BUCKET in this order — ` +
  `PRIMARY (phenotype-matched), SECONDARY (ACMG SF), **CARRIER** (recessive het — the tier is real but ` +
  `it is NOT the diagnosis, reproductive relevance only), **PROBABLE-PATHOGENIC VUS** (a VUS the engine ` +
  `held but that is phenotype-relevant + molecularly suggestive — list its signals, note it is ` +
  `prioritised for expert+Claude exploration, tier UNCHANGED), then OTHER; plus the Performance/timing ` +
  `lines. For each met PVS1, carry its mechanism basis (constraint / ClinGen HI=3 / AR phenotype). ` +
  `Do NOT re-classify — report the engine's calls and routing faithfully.`,
  { label: '🤖 synthesize laudo', phase: '🤖 Claude · Laudo' })

return { sample: SAMPLE, out: OUT, hpo: hpoFile, annotated, deps, inspect, triage, classify, qc, report }
