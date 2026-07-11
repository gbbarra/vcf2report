export const meta = {
  name: 'vcf2report-analyze',
  description: 'Run the vcf2report pipeline for one sample as visible, named phases',
  phases: [
    { title: 'Setup', detail: 'liftover GRCh37->38 if needed' },
    { title: 'Frequencies', detail: 'build the local gnomAD table' },
    { title: 'Classify', detail: 'parse -> QC -> annotate -> ACMG -> report' },
    { title: 'Report', detail: 'collect the laudo facts for rendering' },
  ],
}

// args (from the analyze-vcf skill): { repo, vcf, sample, hpo?, out?, lift?, chain?,
// buildGnomad?, jobs? }. The skill fills these after locating the repo and the VCF.
const a = args || {}
const REPO = a.repo || '.'
const SAMPLE = a.sample || 'sample'
const OUT = a.out || `out/${SAMPLE}`
const HPO = a.hpo ? `--hpo ${a.hpo}` : ''
const JOBS = a.jobs || 24

if (!a.vcf) {
  log('No VCF path in args — nothing to analyze.')
  return { error: 'no vcf' }
}
let vcf = a.vcf

phase('Setup')
if (a.lift) {
  const chain = a.chain ? `--chain ${a.chain}` : ''
  const lifted = `${OUT}/${SAMPLE}.hg38.vcf`
  await agent(
    `cd ${REPO} && mkdir -p ${OUT} && VCF2REPORT_ALLOW_NETWORK=1 python3 ` +
    `scripts/liftover_to_grch38.py ${vcf} ${lifted} ${chain} 2>&1 | tail -3. ` +
    `Return the [liftover] summary line (records/lifted/skipped).`,
    { label: 'liftover', phase: 'Setup' })
  vcf = lifted
} else {
  log('GRCh38 input — no liftover needed.')
}

phase('Frequencies')
if (a.buildGnomad) {
  await agent(
    `cd ${REPO} && VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_gnomad_local.py ` +
    `--from-vcf ${vcf} --jobs ${JOBS} 2>&1 | tail -6. ` +
    `Return rows written, skipped, and the wall time if present.`,
    { label: 'gnomad-local', phase: 'Frequencies' })
} else {
  log('Using the existing local gnomAD table / bundled data.')
}

phase('Classify')
const classify = await agent(
  `cd ${REPO} && python3 scripts/run_headless.py ${vcf} ${HPO} --sample-id ${SAMPLE} ` +
  `--out ${OUT} --timing 2>&1 | tail -40. Return: the candidate list (gene -> tier), ` +
  `the QC funnel counts, and the per-stage timings.`,
  { label: 'run_headless', phase: 'Classify' })

phase('Report')
const report = await agent(
  `Read ${REPO}/${OUT}/${SAMPLE}_report.md and return, as compact structured text for ` +
  `rendering a laudo: the Conclusion bullets; the Sequencing-quality numbers; the primary ` +
  `and secondary findings rows (gene, transcript, HGVS, consequence, tier); and the ` +
  `Performance/timing lines.`,
  { label: 'report', phase: 'Report' })

return { sample: SAMPLE, out: OUT, classify, report }
