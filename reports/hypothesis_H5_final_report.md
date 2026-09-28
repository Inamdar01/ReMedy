# H5 — ReMedy Mechanism-Path Overlap Validation

## Hypothesis

**H5:** ReMedy treatment-use drug–disease pairs show greater overlap with DrugMechDB mechanism paths than expected under a random-pair null model.

### Null hypothesis (H0)

The observed ReMedy treatment-use pairs do not have greater mechanism-path overlap than would be expected from random selection of the same number of drug–disease pairs from the eligible DrugBank × target-disease universe.

### Alternative hypothesis (H1)

The observed ReMedy treatment-use pairs have greater mechanism-path overlap than expected under the random-pair null model.

## Operational definition

A pair is classified as having **mechanism overlap** when:

1. the DrugBank drug has a validated DrugBank→ChEMBL mapping;
2. the disease is one of the 10 ReMedy diseases with a direct MONDO→MeSH cross-reference;
3. the DrugMechDB mechanism path contains at least one UniProt protein; and
4. at least one of those proteins is also present in the ReMedy drug→disease biological connectivity.

ReMedy drug→protein evidence uses the union of the ChEMBL and PrimeKG Drug→Protein canonical relations.

Only ReMedy `indication` and `off-label use` relationships are treated as observed treatment-use pairs. `contraindication` relationships are excluded.

The two ReMedy diseases without a direct MeSH cross-reference in the downloaded MONDO ontology (`MONDO:0005090` schizophrenia and `MONDO:0007743` attention deficit-hyperactivity disorder) are excluded from this direct-ID validation rather than matched by name.

## Statistical test

A one-sided **hypergeometric enrichment test** is used.

- Universe size, **N** = 1180
- Mechanism-overlap pairs in the universe, **K** = 61
- Observed ReMedy treatment-use pairs, **n** = 169
- Observed mechanism-overlap pairs, **x** = 60
- Expected overlap under the null = **8.7364**
- Observed overlap rate = **35.50%**
- Universe overlap rate = **5.17%**
- Enrichment ratio = **6.8678×**
- Hypergeometric p-value = **3.64916e-54**

### 95% confidence interval

A 95% Wilson interval for the observed overlap proportion is:

**28.68% to 42.97%**

This interval describes uncertainty around the observed pair-level overlap proportion; it is not a confidence interval for clinical efficacy.

## Result

The observed overlap was **60 of 169 pairs (35.50%)**, compared with **8.74 pairs expected under the random-pair null**. The observed overlap was **6.87×** the random expectation. The one-sided hypergeometric p-value was **3.649e-54**.

Under this operational definition, the observed ReMedy treatment-use pairs show substantially greater exact UniProt-level mechanism overlap with DrugMechDB than expected under the specified random-pair null.

## Interpretation

This result supports **H5 as a dataset-level consistency/validation finding**: the ReMedy treatment-use pairs overlap curated mechanism paths more often than random pairs in the defined universe.

It does **not** establish clinical efficacy, causal mechanism, or therapeutic suitability. DrugMechDB overlap is evidence of mechanistic agreement with a curated knowledge source, not proof that a drug will work for a patient.

## Limitations

1. H5 validates **exact UniProt overlap**, so mechanisms expressed through other node types are not counted unless a shared UniProt protein is present.
2. Only 10 of the 12 ReMedy diseases have direct MONDO→MeSH mappings and are therefore included.
3. The mechanism-overlap universe is constructed from the ReMedy biological graph and DrugMechDB; the test evaluates enrichment of observed treatment labels within that predefined universe.
4. A hypergeometric test treats candidate pairs as the sampling units; it does not model dependence among pairs sharing drugs, diseases, or proteins.
5. The mechanism-overlap definition is an operational proxy for validation and should not be described as clinical validation.

## Reproducibility

DrugMechDB input:

`data/raw/drugmechdb/indication_paths.yaml`

SHA-256:

`b973118ff3a43450aac9facc9957ec5ce1310efa049ed11b6e94b42807050ee7`

The executable analysis is:

`src/validation/hypothesis_H5_final.py`

The pair-level output is:

`reports/hypothesis_H5_pair_results.csv`

The report is:

`reports/hypothesis_H5_final_report.md`
