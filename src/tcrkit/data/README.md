# Data

| file | what | source |
|---|---|---|
| `mhc_uniprot.fasta` | MHC protein sequences, exactly as UniProt serves them | UniProt (CC BY 4.0) |
| `mhc_annotations.csv` | one row per MHC allele: its UniProt entry, membrane topology, and what `annotate_mhc` needs | DeepTMHMM, MHC Motif Atlas, tcrkit |
| `imgt_gdomains.csv` | IMGT-numbered G-domain templates, for `annotate_mhc`'s IMGT numbering | IMGT (CC BY 4.0) |
| `imgt_cdomains.csv` | IMGT-numbered C / C-like domain templates (B2M, TR/IG constant, MHC α3/β2), for `domains=True` on `anarci` / `annotate_mhc` | IMGT (CC BY 4.0) |
| `cdr3_junction_lookup.csv` | per-gene CDR3 junction residues, for `add_cdr3_junctions` | OTS (CC BY 4.0) |

`tcrkit.load_mhc_table()` joins the first two.

## Rebuilding

```bash
python -m tcrkit._build.fetch_mhc_sequences        # the FASTA, from the accessions in the CSV
python -m tcrkit._build.annotate_mhc_deeptmhmm     # DeepTMHMM topology into the CSV (needs pybiolib)
python -m tcrkit._build.fetch_imgt_gdomains        # imgt_gdomains.csv
python -m tcrkit._build.fetch_imgt_cdomains        # imgt_cdomains.csv
```

## mhc_uniprot.fasta, mhc_annotations.csv

186 alleles (160 human, 24 mouse, 2 B2M), each with its own UniProt entry identical to
the allele (B2M_mouse differs at one residue, the allotype site). H2-AAs, H2-AAu and
H2-Dq are mature chains without a signal peptide; H2-Dq comes from UniParc
(`UPI00005FF33A`), as UniProtKB has no H-2D^q.

Then 76 MHC and B2M chains of 20 other species, for gene- and species-level matching:
the UniProt entries RCSB (SIFTS) maps the MHC and B2M chains of TCR3d's other-species
pMHC structures to, one per entry (near-identical ones left out). `note` names the
structures and the identity over their chain. Two have no UniProtKB entry and come from
UniParc, as deposited: the black flying fox B2M (`UPI00116155BB`; its UniProtKB entry is
a 331-residue gene model) and the big brown bat DR beta (`UPI0032B08EC0`). Per species:
rhesus macaque (*Macaca mulatta*) 6; pig (*Sus scrofa*) 13; cattle (*Bos taurus*) 3; horse (*Equus caballus*) 1; dog (*Canis lupus familiaris*) 3; cat (*Felis catus*) 4; giant panda (*Ailuropoda melanoleuca*) 2; rabbit (*Oryctolagus cuniculus*) 1; rat (*Rattus norvegicus*) 3; black flying fox (*Pteropus alecto*) 2; little brown bat (*Myotis lucifugus*) 1; big brown bat (*Eptesicus fuscus*) 2; brushtail possum (*Trichosurus vulpecula*) 1; chicken (*Gallus gallus*) 13; duck (*Anas platyrhynchos*) 9; great reed warbler (*Acrocephalus arundinaceus*) 2; grass carp (*Ctenopharyngodon idella*) 4; clawed frog (*Xenopus laevis*) 2; green anole (*Anolis carolinensis*) 2; nurse shark (*Ginglymostoma cirratum*) 2.

| column | meaning |
|---|---|
| `allele`, `species`, `organism`, `gene`, `mhc_class`, `mhc_prefix` | the allele as `normalize_mhc` writes it (other species: TCR3d's name), and descriptors |
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

## MHC Motif Atlas (not shipped)

`annotate_mhc`'s alleles (by default), `source='motifatlas'`, `peptide_anchors` and
`peptide_motif` use the [MHC Motif Atlas](http://mhcmotifatlas.org)'s four data files,
downloaded on first use into `~/.cache/tcrkit/motifatlas/` (`manifest.json` there records
the URL, date and SHA-256 of each): the class I and class II sequences (13,966 class I
grooves, 865 class II chains) and their mass-spectrometry ligands (with the class II
binding cores). Names are normalized once into `table.csv`; class II chains get the
topology of the matching bundled allele.

Cite: Tadros DM, Eggenschwiler S, Racle J, Gfeller D. The MHC Motif Atlas: a database of
MHC binding specificities and ligands. *Nucleic Acids Research* 51(D1), D428-D437 (2023).
https://doi.org/10.1093/nar/gkac965

## imgt_gdomains.csv

Parsed from IMGT's G-DOMAIN page
(https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/G-DOMAIN/Gdomains.html): the MH1
G-ALPHA1/G-ALPHA2 and MH2 G-ALPHA/G-BETA sections, human, mouse and rat - 42 entries.
Residues upper-cased, gaps removed; positions before 1 are `'-10'`..`'-1'`, those
between strand D and the helix `'49.1'`..`'49.7'`.

Cite: Lefranc M.-P. et al., IMGT unique numbering for the G-DOMAIN and G-LIKE-DOMAIN,
*Dev. Comp. Immunol.* 29, 917-938 (2005); IMGT®, https://www.imgt.org.

## imgt_cdomains.csv

Parsed from IMGT's C-DOMAIN / C-LIKE-DOMAIN protein display
(https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/C-DOMAIN/CLdomains.html), human
and mouse rows. Positions 1-128 with 46-76 absent; tcrkit adds 1000 when a chain also
has a V or G domain.

Cite: Lefranc M.-P. et al., IMGT unique numbering for immunoglobulin and T cell receptor
constant domains and Ig superfamily C-like domains, *Dev. Comp. Immunol.* 29, 185-203
(2005); IMGT®, https://www.imgt.org.

## cdr3_junction_lookup.csv

Per-V/J-gene CDR3 junction residues aggregated from OTS (CC BY 4.0; aggregates only, no
individual sequences).

Cite: Raybould M.I.J. et al., The Observed T Cell Receptor Space database, *Cell Reports*
43(9), 114704 (2024), https://doi.org/10.1016/j.celrep.2024.114704.
