# Input provenance and redistribution rules

The original data-source artifacts belong to their respective publishers. This public code-only repository deliberately excludes raw source bytes and certain derived records until a file-by-file rights decision is documented.

| Source | Use in research | Public release treatment |
|---|---|---|
| MovingAI street maps | closed-cell topology and engineering clearances | link to official source; preserve exact hash and attribution; do not mirror without permissions confirmation |
| Copernicus EMSR659 | source-derived flood polygon certificates | use official product reference; honour product-specific conditions; do not include original archive in this GitHub code-only stage |
| USGS PAD-US | protected-area transfer tests | attribute official publication and release version; raw GeoJSON omitted at this stage |
| sakayph Metro Manila GTFS | threshold/window proof case | custom DOTC developer agreement: no source rows or feed in this stage |
| TLLVerifyBench | neural benchmark slice | source Git blobs/hashes only; no model coefficients redistributed pending licence determination |
| factory scheduling workbook | real-data scheduling core | no original workbooks here; reacquire where permitted |
| UC Davis battery power data | cleaned measured-power experiment | source attribution and dataset DOI required; workbook bytes omitted |

The R8 `data/map_source_lock.json` contains source identities and hashes, not the benchmark map bytes.

A reviewer granted a confidential Zenodo draft link may inspect deposited code-only materials. **Do not place that secret reviewer URL in a public GitHub README.**
