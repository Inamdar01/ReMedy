# ReMedy — Explainable AI Candidate Report

## Purpose

This report summarizes the top 20 candidates from the existing R-GCN ranking and connects each candidate to available knowledge-graph and ChEMBL drug-target evidence.

## Important interpretation limits

- R-GCN scores are ranking scores, not calibrated probabilities.
- Shared proteins and pathways indicate graph connections, not proven causal mechanisms.
- Drug-target assay evidence does not establish efficacy for the candidate disease.
- Missing evidence in this export is not proof that no evidence exists elsewhere.
- All candidates are research hypotheses, not clinical advice.

## Candidate overview

| Rank | Drug | Disease | R-GCN score | Shared proteins | Shared pathways | ChEMBL assay records |
|---:|---|---|---:|---:|---:|---:|
| 1 | CAPTODIAME | schizophrenia | 0.487865 | 3 | 15 | 0 |
| 2 | GEPIRONE | schizophrenia | 0.487848 | 3 | 9 | 8 |
| 3 | CAPTODIAME | Parkinson disease | 0.487440 | 0 | 11 | 0 |
| 4 | GEPIRONE | Parkinson disease | 0.487399 | 1 | 7 | 8 |
| 5 | CAPTODIAME | major depressive disorder | 0.487360 | 2 | 15 | 0 |
| 6 | FROVATRIPTAN | schizophrenia | 0.487194 | 3 | 56 | 7 |
| 7 | GEPIRONE | migraine disorder | 0.487164 | 1 | 8 | 8 |
| 8 | CAPTODIAME | migraine disorder | 0.487142 | 0 | 9 | 0 |
| 9 | CITALOPRAM | schizophrenia | 0.486904 | 6 | 44 | 0 |
| 10 | FROVATRIPTAN | Parkinson disease | 0.486795 | 0 | 29 | 7 |
| 11 | CAPTODIAME | epilepsy | 0.486791 | 2 | 15 | 0 |
| 12 | GEPIRONE | epilepsy | 0.486782 | 2 | 9 | 8 |
| 13 | METHYLPHENIDATE | schizophrenia | 0.486704 | 4 | 20 | 10 |
| 14 | FROVATRIPTAN | major depressive disorder | 0.486645 | 2 | 54 | 7 |
| 15 | VILAZODONE | schizophrenia | 0.486562 | 12 | 109 | 24 |
| 16 | CITALOPRAM | Parkinson disease | 0.486516 | 2 | 35 | 0 |
| 17 | PHENYLTOLOXAMINE | schizophrenia | 0.486469 | 5 | 9 | 5 |
| 18 | EPICRIPTINE | schizophrenia | 0.486335 | 5 | 9 | 0 |
| 19 | DEXMETHYLPHENIDATE | schizophrenia | 0.486315 | 3 | 14 | 4 |
| 20 | FROVATRIPTAN | anxiety disorder | 0.486281 | 4 | 59 | 7 |

## Candidate explanations

### 1. CAPTODIAME — schizophrenia

The R-GCN ranked this drug-disease pair at position 1 with a score of 0.487865. This score is an uncalibrated ranking score, not a probability. The graph contains 3 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 15 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR2C;DRD3;SIGMAR1

**Shared pathways:** Signal Transduction;Disease;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;Serotonin receptors;G alpha (q) signalling events;G alpha (i) signalling events;GPCR ligand binding;Infectious disease;Potential therapeutics for SARS;SARS-CoV Infections;Viral Infection Pathways

---

### 2. GEPIRONE — schizophrenia

The R-GCN ranked this drug-disease pair at position 2 with a score of 0.487848. This score is an uncalibrated ranking score, not a probability. The graph contains 3 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 9 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 8 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR1A;DRD2;HTR2A

**Shared pathways:** Signal Transduction;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;Serotonin receptors;G alpha (q) signalling events;GPCR ligand binding

**Targets with ChEMBL assay records:** DRD2; HTR1A; HTR2A

**Reported activity types:** Activity; EC50; IC50; Ki; pKd

**Available assay details:**

- DRD2: Ki = 4320.0 nM | assay=Assay:CHEMBL4037101 | document=Document:CHEMBL4033775
- DRD2: Ki = 58.0 nM | assay=Assay:CHEMBL870200 | document=Document:CHEMBL1147799
- HTR1A: Activity = 80.0 % | assay=Assay:CHEMBL870203 | document=Document:CHEMBL1147799
- HTR1A: EC50 = 288.0 nM | assay=Assay:CHEMBL870204 | document=Document:CHEMBL1147799
- HTR1A: IC50 = 114.0 nM | assay=Assay:CHEMBL616273 | document=Document:CHEMBL1148245
- Plus 3 additional assay record(s); see the CSV report.

---

### 3. CAPTODIAME — Parkinson disease

The R-GCN ranked this drug-disease pair at position 3 with a score of 0.487440. This score is an uncalibrated ranking score, not a probability. No shared protein connections were found by the current neighborhood-tracing method. The two graph neighborhoods share 11 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** nan

**Shared pathways:** Signal Transduction;Disease;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;G alpha (i) signalling events;GPCR ligand binding;Infectious disease;Viral Infection Pathways

---

### 4. GEPIRONE — Parkinson disease

The R-GCN ranked this drug-disease pair at position 4 with a score of 0.487399. This score is an uncalibrated ranking score, not a probability. The graph contains 1 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 7 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 8 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** DRD2

**Shared pathways:** Signal Transduction;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;GPCR ligand binding

**Targets with ChEMBL assay records:** DRD2; HTR1A; HTR2A

**Reported activity types:** Activity; EC50; IC50; Ki; pKd

**Available assay details:**

- DRD2: Ki = 4320.0 nM | assay=Assay:CHEMBL4037101 | document=Document:CHEMBL4033775
- DRD2: Ki = 58.0 nM | assay=Assay:CHEMBL870200 | document=Document:CHEMBL1147799
- HTR1A: Activity = 80.0 % | assay=Assay:CHEMBL870203 | document=Document:CHEMBL1147799
- HTR1A: EC50 = 288.0 nM | assay=Assay:CHEMBL870204 | document=Document:CHEMBL1147799
- HTR1A: IC50 = 114.0 nM | assay=Assay:CHEMBL616273 | document=Document:CHEMBL1148245
- Plus 3 additional assay record(s); see the CSV report.

---

### 5. CAPTODIAME — major depressive disorder

The R-GCN ranked this drug-disease pair at position 5 with a score of 0.487360. This score is an uncalibrated ranking score, not a probability. The graph contains 2 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 15 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR2C;DRD3

**Shared pathways:** Signal Transduction;Disease;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;Serotonin receptors;G alpha (q) signalling events;G alpha (i) signalling events;GPCR ligand binding;Infectious disease;Potential therapeutics for SARS;SARS-CoV Infections;Viral Infection Pathways

---

### 6. FROVATRIPTAN — schizophrenia

The R-GCN ranked this drug-disease pair at position 6 with a score of 0.487194. This score is an uncalibrated ranking score, not a probability. The graph contains 3 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 56 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 7 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** CYP1A2;HTR1A;HTR1B

**Shared pathways:** Developmental Biology;Metabolism;Phase II - Conjugation of compounds;Methylation;Signaling by NOTCH;Signal Transduction;Macroautophagy;Disease;Organelle biogenesis and maintenance;Signaling by NOTCH1;Biological oxidations;Cytochrome P450 - arranged by substrate type;Phase I - Functionalization of compounds;Aromatic amines can be N-hydroxylated or N-dealkylated by CYP1A2;Xenobiotics;NOTCH1 Intracellular Domain Regulates Transcription;Generic Transcription Pathway;Synthesis of epoxy (EET) and dihydroxyeicosatrienoic acids (DHET);Arachidonate metabolism;Synthesis of (16-20)-hydroxyeicosatetraenoic acids (HETE);Cellular responses to stress;Signaling by NOTCH1 PEST Domain Mutants in Cancer;Signaling by NOTCH1 in Cancer;Constitutive Signaling by NOTCH1 PEST Domain Mutants;Signaling by NOTCH1 HD+PEST Domain Mutants in Cancer;Constitutive Signaling by NOTCH1 HD+PEST Domain Mutants;Cellular response to heat stress;Notch-HLH transcription pathway;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Serotonin receptors;G alpha (i) signalling events;GPCR ligand binding;Aflatoxin activation and detoxification;Metabolism of lipids;Cilium Assembly;Diseases of signal transduction by growth factor receptors and second messengers;RNA Polymerase II Transcription;Gene expression (Transcription);Transcriptional regulation by RUNX2;RUNX2 regulates osteoblast differentiation;RUNX2 regulates bone development;Cellular responses to stimuli;Fatty acid metabolism;Biosynthesis of DHA-derived SPMs;Biosynthesis of specialized proresolving mediators (SPMs);Biosynthesis of protectins;Biosynthesis of maresins;Biosynthesis of maresin-like SPMs;Autophagy;Selective autophagy;Differentiation of T cells;Assembly of the 9+0 primary cilium;Differentiation of naive CD4+ T cells to T helper 2 cells (Th2 cells)

**Targets with ChEMBL assay records:** HDAC6; HTR1A; HTR1B; HTR1D

**Reported activity types:** Ki; inhibition

**Available assay details:**

- HDAC6: Inhibition = -10.95 % | assay=Assay:CHEMBL4808150 | document=Document:CHEMBL4808148
- HDAC6: Inhibition = 9.03 % | assay=Assay:CHEMBL4808149 | document=Document:CHEMBL4808148
- HTR1A: Ki = 62.0 nM | assay=Assay:CHEMBL616130 | document=Document:CHEMBL1132575
- HTR1B: Ki = 10.3 nM | assay=Assay:CHEMBL616430 | document=Document:CHEMBL1132575
- HTR1B: Ki = 4.4 nM | assay=Assay:CHEMBL5346609 | document=Document:CHEMBL5344475
- Plus 2 additional assay record(s); see the CSV report.

---

### 7. GEPIRONE — migraine disorder

The R-GCN ranked this drug-disease pair at position 7 with a score of 0.487164. This score is an uncalibrated ranking score, not a probability. The graph contains 1 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 8 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 8 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR2A

**Shared pathways:** Signal Transduction;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Serotonin receptors;G alpha (q) signalling events;GPCR ligand binding

**Targets with ChEMBL assay records:** DRD2; HTR1A; HTR2A

**Reported activity types:** Activity; EC50; IC50; Ki; pKd

**Available assay details:**

- DRD2: Ki = 4320.0 nM | assay=Assay:CHEMBL4037101 | document=Document:CHEMBL4033775
- DRD2: Ki = 58.0 nM | assay=Assay:CHEMBL870200 | document=Document:CHEMBL1147799
- HTR1A: Activity = 80.0 % | assay=Assay:CHEMBL870203 | document=Document:CHEMBL1147799
- HTR1A: EC50 = 288.0 nM | assay=Assay:CHEMBL870204 | document=Document:CHEMBL1147799
- HTR1A: IC50 = 114.0 nM | assay=Assay:CHEMBL616273 | document=Document:CHEMBL1148245
- Plus 3 additional assay record(s); see the CSV report.

---

### 8. CAPTODIAME — migraine disorder

The R-GCN ranked this drug-disease pair at position 8 with a score of 0.487142. This score is an uncalibrated ranking score, not a probability. No shared protein connections were found by the current neighborhood-tracing method. The two graph neighborhoods share 9 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** nan

**Shared pathways:** Signal Transduction;Disease;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Serotonin receptors;G alpha (q) signalling events;GPCR ligand binding

---

### 9. CITALOPRAM — schizophrenia

The R-GCN ranked this drug-disease pair at position 9 with a score of 0.486904. This score is an uncalibrated ranking score, not a probability. The graph contains 6 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 44 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** CYP1A2;ABCB1;CYP3A4;CYP2D6;SLC6A4;HRH1

**Shared pathways:** Neurotransmitter clearance;Transmission across Chemical Synapses;Neuronal System;Metabolism;Phase II - Conjugation of compounds;Methylation;Signal Transduction;Biological oxidations;Cytochrome P450 - arranged by substrate type;Fatty acids;Phase I - Functionalization of compounds;Aromatic amines can be N-hydroxylated or N-dealkylated by CYP1A2;Miscellaneous substrates;Xenobiotics;CYP2E1 reactions;Synthesis of epoxy (EET) and dihydroxyeicosatrienoic acids (DHET);Arachidonate metabolism;Synthesis of (16-20)-hydroxyeicosatetraenoic acids (HETE);Abacavir transmembrane transport;Abacavir ADME;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;Serotonin clearance from the synaptic cleft;Transport of small molecules;ABC-family protein mediated transport;GPCR downstream signalling;Histamine receptors;G alpha (q) signalling events;SLC-mediated transmembrane transport;SLC-mediated transport of neurotransmitters;GPCR ligand binding;Aflatoxin activation and detoxification;Metabolism of lipids;Fatty acid metabolism;Biosynthesis of DHA-derived SPMs;Biosynthesis of specialized proresolving mediators (SPMs);Biosynthesis of protectins;Biosynthesis of maresins;Biosynthesis of maresin-like SPMs;Drug ADME;Aspirin ADME;Atorvastatin ADME;Prednisone ADME

---

### 10. FROVATRIPTAN — Parkinson disease

The R-GCN ranked this drug-disease pair at position 10 with a score of 0.486795. This score is an uncalibrated ranking score, not a probability. No shared protein connections were found by the current neighborhood-tracing method. The two graph neighborhoods share 29 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 7 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** nan

**Shared pathways:** Developmental Biology;Metabolism;Signal Transduction;Disease;Organelle biogenesis and maintenance;Biological oxidations;Cytochrome P450 - arranged by substrate type;Phase I - Functionalization of compounds;Xenobiotics;Generic Transcription Pathway;Arachidonate metabolism;Cellular responses to stress;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;G alpha (i) signalling events;GPCR ligand binding;Metabolism of lipids;RNA Polymerase II Transcription;Gene expression (Transcription);Transcriptional regulation by RUNX2;Cellular responses to stimuli;Fatty acid metabolism;Biosynthesis of DHA-derived SPMs;Biosynthesis of specialized proresolving mediators (SPMs);Biosynthesis of maresins;Biosynthesis of maresin-like SPMs;Differentiation of T cells

**Targets with ChEMBL assay records:** HDAC6; HTR1A; HTR1B; HTR1D

**Reported activity types:** Ki; inhibition

**Available assay details:**

- HDAC6: Inhibition = -10.95 % | assay=Assay:CHEMBL4808150 | document=Document:CHEMBL4808148
- HDAC6: Inhibition = 9.03 % | assay=Assay:CHEMBL4808149 | document=Document:CHEMBL4808148
- HTR1A: Ki = 62.0 nM | assay=Assay:CHEMBL616130 | document=Document:CHEMBL1132575
- HTR1B: Ki = 10.3 nM | assay=Assay:CHEMBL616430 | document=Document:CHEMBL1132575
- HTR1B: Ki = 4.4 nM | assay=Assay:CHEMBL5346609 | document=Document:CHEMBL5344475
- Plus 2 additional assay record(s); see the CSV report.

---

### 11. CAPTODIAME — epilepsy

The R-GCN ranked this drug-disease pair at position 11 with a score of 0.486791. This score is an uncalibrated ranking score, not a probability. The graph contains 2 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 15 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** DRD3;SIGMAR1

**Shared pathways:** Signal Transduction;Disease;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;Serotonin receptors;G alpha (q) signalling events;G alpha (i) signalling events;GPCR ligand binding;Infectious disease;Potential therapeutics for SARS;SARS-CoV Infections;Viral Infection Pathways

---

### 12. GEPIRONE — epilepsy

The R-GCN ranked this drug-disease pair at position 12 with a score of 0.486782. This score is an uncalibrated ranking score, not a probability. The graph contains 2 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 9 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 8 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR1A;DRD2

**Shared pathways:** Signal Transduction;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;Serotonin receptors;G alpha (q) signalling events;GPCR ligand binding

**Targets with ChEMBL assay records:** DRD2; HTR1A; HTR2A

**Reported activity types:** Activity; EC50; IC50; Ki; pKd

**Available assay details:**

- DRD2: Ki = 4320.0 nM | assay=Assay:CHEMBL4037101 | document=Document:CHEMBL4033775
- DRD2: Ki = 58.0 nM | assay=Assay:CHEMBL870200 | document=Document:CHEMBL1147799
- HTR1A: Activity = 80.0 % | assay=Assay:CHEMBL870203 | document=Document:CHEMBL1147799
- HTR1A: EC50 = 288.0 nM | assay=Assay:CHEMBL870204 | document=Document:CHEMBL1147799
- HTR1A: IC50 = 114.0 nM | assay=Assay:CHEMBL616273 | document=Document:CHEMBL1148245
- Plus 3 additional assay record(s); see the CSV report.

---

### 13. METHYLPHENIDATE — schizophrenia

The R-GCN ranked this drug-disease pair at position 13 with a score of 0.486704. This score is an uncalibrated ranking score, not a probability. The graph contains 4 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 20 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 10 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR1A;SLC6A2;SLC6A4;SLC6A3

**Shared pathways:** Neurotransmitter clearance;Transmission across Chemical Synapses;Neuronal System;Signal Transduction;Disease;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;Dopamine clearance from the synaptic cleft;Serotonin clearance from the synaptic cleft;Transport of small molecules;Serotonin receptors;SLC-mediated transmembrane transport;SLC-mediated transport of neurotransmitters;GPCR ligand binding;Defective neurotransmitter clearance by SLC6A3 causes Parkinsonism-dystonia infantile (PKDYS);SLC transporter disorders;Defective SLC6A2 causes orthostatic intolerance (OI);Disorders of transmembrane transporters;Defective transport of neurotransmitters by SLC6A3 causes Parkinsonism-dystonia infantile (PKDYS)

**Targets with ChEMBL assay records:** SLC6A2; SLC6A3; SLC6A4

**Reported activity types:** IC50; Ki

**Available assay details:**

- SLC6A2: IC50 = 61.0 nM | assay=Assay:CHEMBL911129 | document=Document:CHEMBL1139231
- SLC6A2: Ki = 340.0 nM | assay=Assay:CHEMBL3379356 | document=Document:CHEMBL3351632
- SLC6A2: Ki = 660.0 nM | assay=Assay:CHEMBL911128 | document=Document:CHEMBL1139231
- SLC6A3: IC50 = 17.0 nM | assay=Assay:CHEMBL673257 | document=Document:CHEMBL1136769
- SLC6A3: IC50 = 79.0 nM | assay=Assay:CHEMBL911125 | document=Document:CHEMBL1139231
- Plus 5 additional assay record(s); see the CSV report.

---

### 14. FROVATRIPTAN — major depressive disorder

The R-GCN ranked this drug-disease pair at position 14 with a score of 0.486645. This score is an uncalibrated ranking score, not a probability. The graph contains 2 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 54 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 7 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR1A;HTR1B

**Shared pathways:** Developmental Biology;Metabolism;Phase II - Conjugation of compounds;Methylation;Signaling by NOTCH;Signal Transduction;Macroautophagy;Disease;Organelle biogenesis and maintenance;Signaling by NOTCH1;Biological oxidations;Cytochrome P450 - arranged by substrate type;Phase I - Functionalization of compounds;Xenobiotics;NOTCH1 Intracellular Domain Regulates Transcription;Generic Transcription Pathway;Synthesis of epoxy (EET) and dihydroxyeicosatrienoic acids (DHET);Arachidonate metabolism;Synthesis of (16-20)-hydroxyeicosatetraenoic acids (HETE);Cellular responses to stress;Signaling by NOTCH1 PEST Domain Mutants in Cancer;Signaling by NOTCH1 in Cancer;Constitutive Signaling by NOTCH1 PEST Domain Mutants;Signaling by NOTCH1 HD+PEST Domain Mutants in Cancer;Constitutive Signaling by NOTCH1 HD+PEST Domain Mutants;HSF1 activation;Cellular response to heat stress;Notch-HLH transcription pathway;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Serotonin receptors;G alpha (i) signalling events;GPCR ligand binding;Metabolism of lipids;Cilium Assembly;Diseases of signal transduction by growth factor receptors and second messengers;RNA Polymerase II Transcription;Gene expression (Transcription);Transcriptional regulation by RUNX2;Cellular responses to stimuli;Fatty acid metabolism;Biosynthesis of DHA-derived SPMs;Biosynthesis of specialized proresolving mediators (SPMs);Biosynthesis of maresins;Biosynthesis of maresin-like SPMs;Autophagy;Chaperone Mediated Autophagy;Aggrephagy;Selective autophagy;Differentiation of T cells;Assembly of the 9+0 primary cilium;Differentiation of naive CD4+ T cells to T helper 2 cells (Th2 cells)

**Targets with ChEMBL assay records:** HDAC6; HTR1A; HTR1B; HTR1D

**Reported activity types:** Ki; inhibition

**Available assay details:**

- HDAC6: Inhibition = -10.95 % | assay=Assay:CHEMBL4808150 | document=Document:CHEMBL4808148
- HDAC6: Inhibition = 9.03 % | assay=Assay:CHEMBL4808149 | document=Document:CHEMBL4808148
- HTR1A: Ki = 62.0 nM | assay=Assay:CHEMBL616130 | document=Document:CHEMBL1132575
- HTR1B: Ki = 10.3 nM | assay=Assay:CHEMBL616430 | document=Document:CHEMBL1132575
- HTR1B: Ki = 4.4 nM | assay=Assay:CHEMBL5346609 | document=Document:CHEMBL5344475
- Plus 2 additional assay record(s); see the CSV report.

---

### 15. VILAZODONE — schizophrenia

The R-GCN ranked this drug-disease pair at position 15 with a score of 0.486562. This score is an uncalibrated ranking score, not a probability. The graph contains 12 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 109 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 24 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** CYP1A2;CYP3A4;HTR1A;CYP2D6;DRD2;SLC6A2;HTR1B;SLC6A4;HTR7;SLC6A3;KCNH2;AGER

**Shared pathways:** Neurotransmitter clearance;Transmission across Chemical Synapses;Neuronal System;Developmental Biology;Cytokine Signaling in Immune system;Potassium Channels;Voltage gated Potassium channels;Metabolism;Phase II - Conjugation of compounds;Methylation;Signaling by NOTCH;Signal Transduction;Macroautophagy;Disease;Toll Like Receptor 4 (TLR4) Cascade;MyD88:MAL(TIRAP) cascade initiated on plasma membrane;MyD88-independent TLR4 cascade;Toll Like Receptor 9 (TLR9) Cascade;Toll Like Receptor 10 (TLR10) Cascade;Toll Like Receptor 3 (TLR3) Cascade;Toll Like Receptor 5 (TLR5) Cascade;Toll Like Receptor TLR1:TLR2 Cascade;Toll Like Receptor 7/8 (TLR7/8) Cascade;Toll Like Receptor TLR6:TLR2 Cascade;Innate Immune System;Immune System;Toll-like Receptor Cascades;DDX58/IFIH1-mediated induction of interferon-alpha/beta;Toll Like Receptor 2 (TLR2) Cascade;Organelle biogenesis and maintenance;Signaling by NOTCH1;Biological oxidations;Cytochrome P450 - arranged by substrate type;Fatty acids;Phase I - Functionalization of compounds;Aromatic amines can be N-hydroxylated or N-dealkylated by CYP1A2;Miscellaneous substrates;Xenobiotics;CYP2E1 reactions;NOTCH1 Intracellular Domain Regulates Transcription;Generic Transcription Pathway;Synthesis of epoxy (EET) and dihydroxyeicosatrienoic acids (DHET);Arachidonate metabolism;Synthesis of (16-20)-hydroxyeicosatetraenoic acids (HETE);Cellular responses to stress;Signaling by NOTCH1 PEST Domain Mutants in Cancer;Signaling by NOTCH1 in Cancer;Constitutive Signaling by NOTCH1 PEST Domain Mutants;Signaling by NOTCH1 HD+PEST Domain Mutants in Cancer;Constitutive Signaling by NOTCH1 HD+PEST Domain Mutants;Cellular response to heat stress;Notch-HLH transcription pathway;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;Dopamine clearance from the synaptic cleft;Serotonin clearance from the synaptic cleft;Transport of small molecules;GPCR downstream signalling;Dopamine receptors;Serotonin receptors;Muscle contraction;G alpha (s) signalling events;G alpha (i) signalling events;SLC-mediated transmembrane transport;SLC-mediated transport of neurotransmitters;TAK1-dependent IKK and NF-kappa-B activation;Interleukin-1 family signaling;Signaling by Interleukins;GPCR ligand binding;Aflatoxin activation and detoxification;Metabolism of lipids;Phase 3 - rapid repolarisation;Cardiac conduction;Cilium Assembly;Defective neurotransmitter clearance by SLC6A3 causes Parkinsonism-dystonia infantile (PKDYS);SLC transporter disorders;Defective SLC6A2 causes orthostatic intolerance (OI);Disorders of transmembrane transporters;Defective transport of neurotransmitters by SLC6A3 causes Parkinsonism-dystonia infantile (PKDYS);Diseases of signal transduction by growth factor receptors and second messengers;RNA Polymerase II Transcription;Gene expression (Transcription);Advanced glycosylation endproduct receptor signaling;Transcriptional regulation by RUNX2;RUNX2 regulates osteoblast differentiation;RUNX2 regulates bone development;Cellular responses to stimuli;Fatty acid metabolism;Biosynthesis of DHA-derived SPMs;Biosynthesis of specialized proresolving mediators (SPMs);Biosynthesis of protectins;Biosynthesis of maresins;Interleukin-1 signaling;Biosynthesis of maresin-like SPMs;TRAF6 mediated NF-kB activation;TRIF (TICAM1)-mediated TLR4 signaling;Autophagy;Selective autophagy;Drug ADME;Aspirin ADME;TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation;MyD88 dependent cascade initiated on endosome;Atorvastatin ADME;Prednisone ADME;MyD88 cascade initiated on plasma membrane;Differentiation of T cells;Assembly of the 9+0 primary cilium;Differentiation of naive CD4+ T cells to T helper 2 cells (Th2 cells)

**Targets with ChEMBL assay records:** AGER; CYP1A2; CYP2C19; CYP2C9; CYP3A4; DRD2; HDAC6; HTR1A; HTR1B; HTR7; KCNH2; SLC6A2; SLC6A3; SLC6A4

**Reported activity types:** EC50; IC50; INH; Ki; inhibition

**Available assay details:**

- AGER: IC50 > 200000.0 nM | assay=Assay:CHEMBL5165149 | document=Document:CHEMBL5154874
- CYP1A2: Inhibition = 14.4 % | assay=Assay:CHEMBL5165165 | document=Document:CHEMBL5154874
- CYP2C19: Inhibition = 51.3 % | assay=Assay:CHEMBL5165167 | document=Document:CHEMBL5154874
- CYP2C9: Inhibition = 49.2 % | assay=Assay:CHEMBL5165166 | document=Document:CHEMBL5154874
- CYP3A4: Inhibition = 29.0 % | assay=Assay:CHEMBL5165168 | document=Document:CHEMBL5154874
- Plus 19 additional assay record(s); see the CSV report.

---

### 16. CITALOPRAM — Parkinson disease

The R-GCN ranked this drug-disease pair at position 16 with a score of 0.486516. This score is an uncalibrated ranking score, not a probability. The graph contains 2 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 35 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** ABCB1;CYP2D6

**Shared pathways:** Neurotransmitter clearance;Transmission across Chemical Synapses;Neuronal System;Metabolism;Signal Transduction;Biological oxidations;Cytochrome P450 - arranged by substrate type;Fatty acids;Phase I - Functionalization of compounds;Miscellaneous substrates;Xenobiotics;CYP2E1 reactions;Arachidonate metabolism;Abacavir transmembrane transport;Abacavir ADME;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;Serotonin clearance from the synaptic cleft;Transport of small molecules;ABC-family protein mediated transport;GPCR downstream signalling;SLC-mediated transmembrane transport;SLC-mediated transport of neurotransmitters;GPCR ligand binding;Metabolism of lipids;Fatty acid metabolism;Biosynthesis of DHA-derived SPMs;Biosynthesis of specialized proresolving mediators (SPMs);Biosynthesis of maresins;Biosynthesis of maresin-like SPMs;Drug ADME;Aspirin ADME;Atorvastatin ADME;Prednisone ADME

---

### 17. PHENYLTOLOXAMINE — schizophrenia

The R-GCN ranked this drug-disease pair at position 17 with a score of 0.486469. This score is an uncalibrated ranking score, not a probability. The graph contains 5 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 9 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 5 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** DRD2;DRD1;DRD4;DRD5;DRD3

**Shared pathways:** Signal Transduction;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;G alpha (s) signalling events;G alpha (i) signalling events;GPCR ligand binding

**Targets with ChEMBL assay records:** DRD1; DRD2; DRD3; DRD4; DRD5

**Reported activity types:** Ki

**Available assay details:**

- DRD1: Ki = 343.0 nM | assay=Assay:CHEMBL839544 | document=Document:CHEMBL1138672
- DRD2: Ki = 2740.0 nM | assay=Assay:CHEMBL838569 | document=Document:CHEMBL1138672
- DRD3: Ki > 10000.0 nM | assay=Assay:CHEMBL839545 | document=Document:CHEMBL1138672
- DRD4: Ki = 2765.0 nM | assay=Assay:CHEMBL876688 | document=Document:CHEMBL1138672
- DRD5: Ki = 414.0 nM | assay=Assay:CHEMBL839546 | document=Document:CHEMBL1138672

---

### 18. EPICRIPTINE — schizophrenia

The R-GCN ranked this drug-disease pair at position 18 with a score of 0.486335. This score is an uncalibrated ranking score, not a probability. The graph contains 5 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 9 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. No ChEMBL drug-target assay records were found for this pair in the current top-20 evidence export. This does not mean that no such evidence exists elsewhere. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** DRD2;DRD1;DRD4;DRD5;DRD3

**Shared pathways:** Signal Transduction;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Dopamine receptors;G alpha (s) signalling events;G alpha (i) signalling events;GPCR ligand binding

---

### 19. DEXMETHYLPHENIDATE — schizophrenia

The R-GCN ranked this drug-disease pair at position 19 with a score of 0.486315. This score is an uncalibrated ranking score, not a probability. The graph contains 3 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 14 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 4 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** SLC6A2;SLC6A4;SLC6A3

**Shared pathways:** Neurotransmitter clearance;Transmission across Chemical Synapses;Neuronal System;Disease;Dopamine clearance from the synaptic cleft;Serotonin clearance from the synaptic cleft;Transport of small molecules;SLC-mediated transmembrane transport;SLC-mediated transport of neurotransmitters;Defective neurotransmitter clearance by SLC6A3 causes Parkinsonism-dystonia infantile (PKDYS);SLC transporter disorders;Defective SLC6A2 causes orthostatic intolerance (OI);Disorders of transmembrane transporters;Defective transport of neurotransmitters by SLC6A3 causes Parkinsonism-dystonia infantile (PKDYS)

**Targets with ChEMBL assay records:** SLC6A3

**Reported activity types:** IC50; Ki

**Available assay details:**

- SLC6A3: IC50 = 156.0 nM | assay=Assay:CHEMBL1647747 | document=Document:CHEMBL1641574
- SLC6A3: IC50 = 88.0 nM | assay=Assay:CHEMBL673257 | document=Document:CHEMBL1136769
- SLC6A3: Ki = 16.0 nM | assay=Assay:CHEMBL2065334 | document=Document:CHEMBL2062423
- SLC6A3: Ki = 25.0 nM | assay=Assay:CHEMBL1647746 | document=Document:CHEMBL1641574

---

### 20. FROVATRIPTAN — anxiety disorder

The R-GCN ranked this drug-disease pair at position 20 with a score of 0.486281. This score is an uncalibrated ranking score, not a probability. The graph contains 4 shared protein connection(s) in the drug and disease neighborhoods. The two graph neighborhoods share 59 pathway connection(s). These are indirect graph connections, not proof of a causal mechanism. The current evidence export contains 7 ChEMBL drug-target assay record(s) associated with this drug-disease candidate. These support drug-target relationships, not disease-treatment efficacy. This candidate is a research hypothesis requiring independent validation, not a clinical recommendation.

**Shared proteins:** HTR1A;HTR1D;HTR1B;HDAC6

**Shared pathways:** Developmental Biology;Metabolism;Phase II - Conjugation of compounds;Methylation;Signaling by NOTCH;Signal Transduction;Macroautophagy;Disease;Organelle biogenesis and maintenance;Signaling by NOTCH1;Biological oxidations;Cytochrome P450 - arranged by substrate type;Phase I - Functionalization of compounds;Xenobiotics;NOTCH1 Intracellular Domain Regulates Transcription;Generic Transcription Pathway;Synthesis of epoxy (EET) and dihydroxyeicosatrienoic acids (DHET);Arachidonate metabolism;Synthesis of (16-20)-hydroxyeicosatetraenoic acids (HETE);Cellular responses to stress;Signaling by NOTCH1 PEST Domain Mutants in Cancer;Signaling by NOTCH1 in Cancer;Constitutive Signaling by NOTCH1 PEST Domain Mutants;Signaling by NOTCH1 HD+PEST Domain Mutants in Cancer;Constitutive Signaling by NOTCH1 HD+PEST Domain Mutants;HSF1 activation;Cellular response to heat stress;Notch-HLH transcription pathway;Signaling by GPCR;Class A/1 (Rhodopsin-like receptors);Amine ligand-binding receptors;GPCR downstream signalling;Serotonin receptors;G alpha (i) signalling events;GPCR ligand binding;Metabolism of lipids;Cilium Assembly;Cargo trafficking to the periciliary membrane;Diseases of signal transduction by growth factor receptors and second messengers;RNA Polymerase II Transcription;Gene expression (Transcription);Transcriptional regulation by RUNX2;RUNX2 regulates osteoblast differentiation;RUNX2 regulates bone development;Cellular responses to stimuli;Fatty acid metabolism;Biosynthesis of DHA-derived SPMs;Biosynthesis of specialized proresolving mediators (SPMs);Biosynthesis of protectins;Biosynthesis of maresins;Biosynthesis of maresin-like SPMs;Autophagy;Chaperone Mediated Autophagy;Late endosomal microautophagy;Aggrephagy;Selective autophagy;Differentiation of T cells;Assembly of the 9+0 primary cilium;Differentiation of naive CD4+ T cells to T helper 2 cells (Th2 cells)

**Targets with ChEMBL assay records:** HDAC6; HTR1A; HTR1B; HTR1D

**Reported activity types:** Ki; inhibition

**Available assay details:**

- HDAC6: Inhibition = -10.95 % | assay=Assay:CHEMBL4808150 | document=Document:CHEMBL4808148
- HDAC6: Inhibition = 9.03 % | assay=Assay:CHEMBL4808149 | document=Document:CHEMBL4808148
- HTR1A: Ki = 62.0 nM | assay=Assay:CHEMBL616130 | document=Document:CHEMBL1132575
- HTR1B: Ki = 10.3 nM | assay=Assay:CHEMBL616430 | document=Document:CHEMBL1132575
- HTR1B: Ki = 4.4 nM | assay=Assay:CHEMBL5346609 | document=Document:CHEMBL5344475
- Plus 2 additional assay record(s); see the CSV report.

---
