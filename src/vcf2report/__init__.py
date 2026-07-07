"""vcf2report — raw exome VCF to auditable ACMG variant report.

The package is designed as the *engine* behind a Claude Desktop workflow:
Agent Skills encode the clinical SOP, an MCP server exposes thin tool wrappers,
and everything here does the actual work and is importable/testable headless.
"""

__version__ = "0.1.0"

GENOME_BUILD = "GRCh38"
