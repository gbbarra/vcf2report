export const meta = {
  name: 'vcf2report-analyze',
  description: 'Analyze one exome as visible phases — each labeled LOCAL (your machine) or CLAUDE',
  // Phase titles carry who does the work: 🖥️ = local (vcf2report on the machine, deterministic,
  // offline) · 🤖 = Claude (reasoning). They render as the step boxes in Background Tasks.
  phases: [
    { title: '🖥️ Local · Setup', detail: 'liftover GRCh37→38 if needed (bcftools/pyliftover)' },
    { title: '🤖 Claude · Phenotype → HPO', detail: 'map free-text phenotype to HPO terms' },
    { title: '🖥️ Local · Frequencies', detail: 'local gnomAD DuckDB/Parquet store' },
    { title: '🖥️ Local · Classify (ACMG)', detail: 'parse → QC → annotate → ACMG engine' },
    { title: '🤖 Claude · Laudo', detail: 'synthesize the auditable report from the facts' },
  ],
}

// args (from the analyze-vcf skill): { repo, vcf, sample, hpo?, phenotypeText?, out?, lift?,
// chain? }. The skill fills these after locating the repo and the VCF.
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

// 🖥️ LOCAL — build conversion on the machine
phase('🖥️ Local · Setup')
if (a.lift) {
  const chain = a.chain ? `--chain ${a.chain}` : ''
  const lifted = `${OUT}/${SAMPLE}.hg38.vcf`
  await agent(
    `cd ${REPO} && mkdir -p ${OUT} && VCF2REPORT_ALLOW_NETWORK=1 python3 ` +
    `scripts/liftover_to_grch38.py ${vcf} ${lifted} ${chain} 2>&1 | tail -3. ` +
    `Return the [liftover] summary line (records/lifted/skipped).`,
    { label: '🖥️ liftover 37→38', phase: '🖥️ Local · Setup' })
  vcf = lifted
} else {
  log('🖥️ GRCh38 input — no liftover needed.')
}

// 🤖 CLAUDE — turn the clinician's free-text phenotype into HPO ids (only if not given a file)
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

// 🖥️ LOCAL — the population-frequency store lives on the machine, queried offline
phase('🖥️ Local · Frequencies')
await agent(
  `cd ${REPO} && python3 -c "from vcf2report import config; import os; ` +
  `p=config._resolve_gnomad_parquet(); print('gnomAD parquet:', p or 'bundled slice only'); ` +
  `print('exists:', bool(p and os.path.exists(p)))". ` +
  `Return whether a local gnomAD Parquet store is auto-detected (offline frequencies).`,
  { label: '🖥️ gnomAD store', phase: '🖥️ Local · Frequencies' })

// 🖥️ LOCAL — the whole deterministic pipeline runs on the machine, no network, no LLM
phase('🖥️ Local · Classify (ACMG)')
const classify = await agent(
  `cd ${REPO} && python3 scripts/run_headless.py ${vcf} ${HPO} --sample-id ${SAMPLE} ` +
  `--out ${OUT} --timing 2>&1 | tail -40. Return: the candidate list (gene → tier), ` +
  `the QC funnel counts, and the per-stage timings. This is all local + deterministic.`,
  { label: '🖥️ parse→QC→annotate→ACMG', phase: '🖥️ Local · Classify (ACMG)' })

// 🤖 CLAUDE — read the machine's facts and synthesize the auditable narrative
phase('🤖 Claude · Laudo')
const report = await agent(
  `Read ${REPO}/${OUT}/${SAMPLE}_report.md and return, as compact structured text for ` +
  `rendering a laudo: the Conclusion bullets; the Sequencing-quality numbers; the primary ` +
  `and secondary findings rows (gene, transcript, HGVS, consequence, tier, ACMG criteria); ` +
  `and the Performance/timing lines. Do NOT re-classify — report the engine's calls faithfully.`,
  { label: '🤖 synthesize laudo', phase: '🤖 Claude · Laudo' })

return { sample: SAMPLE, out: OUT, hpo: hpoFile, classify, report }
