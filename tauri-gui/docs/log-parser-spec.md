# Log Parser Spec (v1)

## Target lines
The GUI parser reads comment lines emitted by TaxonDBBuilder:

- `# query count taxid=...: N`
- `# total records: N`
- `# matched records: N`
- `# kept records before post_prep: N`
- `# post_prep primer trim: ... removed=N ...`
- `# post_prep length filter: ... removed=N ...`
- `# post_prep duplicate_acc_report: ... groups=N ... cross_organism_groups=N ...`
- `# output: ...`
- `# finished: ...`

## Progress mapping
- Phase 1 Query count: `0-10%`
- Phase 2 Fetch/Parse: `10-80%`
- Phase 3 Post-pPrep: `80-95%`
- Phase 4 Finalize: `95-100%`

## Metrics
- query count per taxid
- matched records
- kept records before post_prep
- primer_trim removed
- length_filter removed
- duplicate groups
- cross_organism_groups
