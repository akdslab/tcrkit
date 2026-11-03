"""Run-once scripts that rebuild tcrkit's bundled data. Not part of the API.

Each is run as a module and writes into ``tcrkit/data/``::

    python -m tcrkit._build.fetch_mhc_sequences       # MHC sequences from UniProt
    python -m tcrkit._build.annotate_mhc_deeptmhmm    # their membrane topology (DeepTMHMM)
    python -m tcrkit._build.fetch_imgt_gdomains       # IMGT G-DOMAIN numbering templates

Run ``fetch_mhc_sequences`` before ``annotate_mhc_deeptmhmm``: the topology belongs to
the sequences. See ``tcrkit/data/README.md`` for what each file holds.
"""
