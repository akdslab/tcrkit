# Data

| file | what | source |
|---|---|---|
| `mhc_uniprot.fasta` | MHC protein sequences, exactly as UniProt serves them | UniProt (CC BY 4.0) |
| `mhc_annotations.csv` | one row per MHC allele: its UniProt entry, membrane topology, and what `annotate_mhc` needs | DeepTMHMM, MHC Motif Atlas, tcrkit |
| `imgt_gdomains.csv` | IMGT-numbered G-domain templates, for `annotate_mhc`'s IMGT numbering | IMGT (CC BY 4.0) |
| `cdr3_junction_lookup.csv` | per-gene CDR3 junction residues, for `add_cdr3_junctions` | OTS (CC BY 4.0) |

`tcrkit.load_mhc_table()` joins the first two.

## Rebuilding

```bash
python -m tcrkit._build.fetch_mhc_sequences        # the FASTA, from the accessions in the CSV
python -m tcrkit._build.annotate_mhc_deeptmhmm     # DeepTMHMM topology into the CSV (needs pybiolib)
python -m tcrkit._build.fetch_imgt_gdomains        # imgt_gdomains.csv
```

## mhc_uniprot.fasta, mhc_annotations.csv

186 alleles (160 human, 24 mouse, 2 B2M), each with its own UniProt entry identical to
the allele (B2M_mouse differs at one residue, the allotype site). H2-AAs, H2-AAu and
H2-Dq are mature chains without a signal peptide; H2-Dq comes from UniParc
(`UPI00005FF33A`), as UniProtKB has no H-2D^q.

| column | meaning |
|---|---|
| `allele`, `species`, `gene`, `mhc_class`, `mhc_prefix` | the allele as `normalize_mhc` writes it, and descriptors |
| `source`, `accession`, `source_release`, `sequence_version` | `uniprot` or `uniparc`, the entry, the UniProt release, its `SV=` |
| `identity_to_previous`, `coverage_of_previous`, `note` | how the entry was chosen |
| `deeptmhmm_prediction`, `deeptmhmm_version`, `topology_length` | DeepTMHMM 1.0.24: `SP+TM`, `TM`, `SP` or `GLOB` |
| `signal_spans`, `outside_spans`, `tmhelix_spans`, `inside_spans` | 1-based residue ranges per region, into the UniProt sequence |
| `mature_start` | signal-peptide length; IMGT position 1 is residue `mature_start + 1` |
| `in_detect`, `detect_rank`, `detect_class` | the alleles `annotate_mhc` matches genes on, and in which order |
| `helix_ref_span`, `helix_ind` | MHC Motif Atlas helix ranges on 18 gene representatives, for `helix='tpm-bench'` |

Credit: The UniProt Consortium, *Nucleic Acids Res.* 53, D609-D617 (2025). Hallgren J.
et al., DeepTMHMM, *bioRxiv* (2022), https://doi.org/10.1101/2022.04.08.487609. Helix
ranges: [MHC Motif Atlas](http://mhcmotifatlas.org).

## DecoderTCR (optional, not shipped)

`source='decodertcr'` uses [DecoderTCR](https://github.com/Biohub/DecoderTCR)'s HLA
reference (MIT): `hla_reference.json` at commit `3e3d9889`, downloaded on first use,
checked against its SHA-256 and cached in `~/.cache/tcrkit/decodertcr/`. Its 88 pairs
become 87 chains; each gets the topology of the matching bundled allele, so regions
mean the same stretch of the chain from either source.

Cite: Lai B, Englund M, Bharanikumar R, Nocedal I, Davariashtiyani A, Perera J, Khan AA.
DecoderTCR: Compositional Pretraining and Entropy-Guided Decoding for TCR-pMHC
Interactions. *ICML* (2026). https://openreview.net/pdf?id=yzes8qBM70

## imgt_gdomains.csv

Parsed from IMGT's G-DOMAIN page
(https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/G-DOMAIN/Gdomains.html): the MH1
G-ALPHA1/G-ALPHA2 and MH2 G-ALPHA/G-BETA sections, human, mouse and rat - 42 entries.
Residues upper-cased, gaps removed; positions before 1 are `'-10'`..`'-1'`, those
between strand D and the helix `'49.1'`..`'49.7'`.

Cite: Lefranc M.-P. et al., IMGT unique numbering for the G-DOMAIN and G-LIKE-DOMAIN,
*Dev. Comp. Immunol.* 29, 917-938 (2005); IMGT®, https://www.imgt.org.
