export const meta = {
  name: 'vcf2report-preflight',
  description: 'vcf2report Stage 1 — dependency check + Parquet store GATE, rendered in Background Tasks BEFORE any analysis (no VCF needed)',
  phases: [
    { title: '🖥️ Local · Dependency check', detail: 'python + annotation tools on PATH, network flag' },
    { title: '🖥️ Local · Parquet store GATE', detail: 'gnomAD · AlphaMissense · ClinVar — availability, version, build date, integrity' },
  ],
}

// args: { repo }. Launch this the moment /vcf2report is invoked — it needs NO VCF, so Stage 1
// (deps + the store gate) shows in the right-side Background Tasks pane before anything is asked.
const a = (typeof args === 'string' ? JSON.parse(args) : args) || {}
const REPO = a.repo || '/Users/gbbarra/vcf2report'

phase('🖥️ Local · Dependency check')
const deps = await agent(
  `cd ${REPO} && python3 scripts/preflight.py. Report in ONE line: python version; which of ` +
  `bcftools/snpEff/vcfanno are on PATH; network on/off. Nothing patient-specific.`,
  { label: '🖥️ preflight (python + tools)', phase: '🖥️ Local · Dependency check' })

phase('🖥️ Local · Parquet store GATE')
const gate = await agent(
  `cd ${REPO} && python3 scripts/check_stores.py --gate; echo "STORES_EXIT=$?". Show the store table ` +
  `VERBATIM — for EACH parquet (gnomAD · AlphaMissense · ClinVar): available? · version · BUILD DATE · ` +
  `integrity · complete — then the READY/BLOCKED banner and the STORES_EXIT line (do not omit it).`,
  { label: '🖥️ gnomAD · AlphaMissense · ClinVar (gate)', phase: '🖥️ Local · Parquet store GATE' })

const ready = /STORES_EXIT=0/.test(gate || '')
log(ready ? '✅ Stores available + intact — ready to analyze (waiting for VCF).'
          : '⛔ Stores not ready — analysis blocked; fix the flagged store(s) and re-run.')
return { ready, deps, gate }
