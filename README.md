# PSS — Source-Aware Boundary-Span Certificates

**Research title:** *Source-Aware Boundary-Span Certificates for Polyhedral Unions under Partial Intersection Information*  
**Corresponding author:** Ghassan Malkawi, Higher Colleges of Technology, Al Ain, UAE  
**Correspondence:** gmalkawi@hct.ac.ae · [ORCID 0000-0002-5320-7561](https://orcid.org/0000-0002-5320-7561)

**Status:** Public, code-only research companion for an unpublished author-review manuscript. This repository is NOT the complete journal reproducibility deposit, and it excludes original third-party datasets and private research archives.

## Scope

This code implements a solver-free rational certificate checker for the declared polyhedral-branch model and an engineering source adapter that converts a binary occupancy grid into closed-square topology and a fixed triangular polyhedral cover. It supports complete and boundary-span certificates whose acceptance is checked independently from the generation routine.

The mathematical contribution is a **conditional proof-information sufficiency rule**. The software does not prove universal runtime gains, safe robot navigation, physical obstacle clearance, or topological properties of raw unsimplified source GIS geometry.

## Built-in source-free smoke test

Requirements: Python 3.11+; the exact certificate engine itself is standard-library only. For the optional generated plot, the packaged R8 engine uses NumPy and Matplotlib (install `pip install -r requirements-optional.txt`).

```bash
python tests/run_synthetic_smoke.py
```

This converts an openly generated **synthetic** 8x8 CSV ring to a strict `.map`, calculates modeled `(beta0,beta1)=(1,1)`, builds both certificates, and requires both to be accepted. This smoke test is NOT an independent field-data validation or a substitute for reproducing the manuscript.

## Engineer-provided data

```bash
python code/csv_to_map.py --input YOUR_STRICT_BINARY_0_1.csv --output output/site.map
python code/pss_engineer.py window --map output/site.map --allow-other --x 0 --y 0 --factor 1 --window 8 --clearance-cells 0 --output output/site_case
```

The site CSV must have uniform row length and only exact values `0` or `1`, with `1` free and `0` blocked. Adjust the region of interest and document the source. `--allow-other` allows an external file with a reported hash; it **does not certify that the source faithfully represents physical space**.

## What this staging package does not contain

It does not redistribute the MovingAI `.map` files, the Metro Manila GTFS rows, third-party TLLVerifyBench coefficients, original factory workbooks, Copernicus/USGS raw datasets, the private full research archive, unpublished reviewer links, passwords, or access tokens.

Full reproduction of the published empirical claims requires the separately frozen, source-hashed evidence archives and rights-permitted re-acquisition. See `SOURCES_AND_RIGHTS.md` and `RELEASE_GATE.md`. The inherited `code/test_engineering.py` requires the original benchmark map inputs and therefore is **not** run as part of the synthetic smoke test.

## Version and provenance

This code snapshot is derived, without modifying the checker algorithms, from the frozen R8 engineering kit inside the PSS V25 reproducibility archive. The latest manuscript author-review successor is **PSS V26R1**, a scientific-result-preserving editorial derivative of PSS V26. R2/R3 adaptive-query investigations are **separate, non-integrated research extensions**, not a claimed part of the accepted manuscript.

No open-source licence or publication DOI has been assigned as of this upload. Public repository visibility does not grant redistribution rights in source data or permission to reuse code under an unstated licence.
