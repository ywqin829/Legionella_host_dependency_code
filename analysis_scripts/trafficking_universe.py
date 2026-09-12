"""Curated human membrane-trafficking gene universe (A3/A4 shared).

Families are assigned manually from HGNC nomenclature. Genes are intersected
with the atlas gene symbols downstream, so imperfect exhaustiveness is safe.
"""

TRAFFICKING = {
    "RAB": [
        "RAB1A", "RAB1B", "RAB2A", "RAB2B", "RAB3A", "RAB3B", "RAB3C", "RAB3D",
        "RAB4A", "RAB4B", "RAB5A", "RAB5B", "RAB5C", "RAB6A", "RAB6B", "RAB6C",
        "RAB7A", "RAB7B", "RAB8A", "RAB8B", "RAB9A", "RAB9B", "RAB10", "RAB11A",
        "RAB11B", "RAB12", "RAB13", "RAB14", "RAB15", "RAB17", "RAB18", "RAB19",
        "RAB20", "RAB21", "RAB22A", "RAB23", "RAB24", "RAB25", "RAB26", "RAB27A",
        "RAB27B", "RAB28", "RAB29", "RAB30", "RAB31", "RAB32", "RAB33A", "RAB33B",
        "RAB34", "RAB35", "RAB36", "RAB37", "RAB38", "RAB39A", "RAB39B", "RAB40A",
        "RAB40B", "RAB40C", "RAB41", "RAB43", "RAB44",
    ],
    "RAB_effector_GEF": ["RABIF", "CHM", "CHML", "GDI1", "GDI2", "DENND1A", "DENND1B", "DENND1C"],
    "RAB_GAP": ["TBC1D20", "TBC1D23", "TBC1D2", "TBC1D4", "TBC1D5", "TBC1D9B", "TBC1D15", "TBC1D17"],
    "ARF": ["ARF1", "ARF3", "ARF4", "ARF5", "ARF6"],
    "ARL": [
        "ARL1", "ARL2", "ARL3", "ARL4A", "ARL4C", "ARL4D", "ARL5A", "ARL5B", "ARL5C",
        "ARL6", "ARL8A", "ARL8B", "ARL9", "ARL10", "ARL11", "ARL13A", "ARL13B", "ARL14",
        "ARL15", "ARL16", "ARL17A", "ARL17B", "ARFRP1",
    ],
    "COPII": ["SAR1A", "SAR1B", "SEC13", "SEC23A", "SEC23B", "SEC23IP", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC31A", "SEC31B"],
    "COPI": ["COPA", "COPB1", "COPB2", "COPE", "COPG1", "COPG2", "COPZ1", "COPZ2", "ARCN1"],
    "SNARE": [
        "VAMP1", "VAMP2", "VAMP3", "VAMP4", "VAMP5", "VAMP7", "VAMP8",
        "SEC22A", "SEC22B", "SEC22C", "STX1A", "STX1B", "STX2", "STX3", "STX4",
        "STX5", "STX6", "STX7", "STX8", "STX10", "STX11", "STX12", "STX16", "STX17",
        "STX18", "STX19", "SNAP23", "SNAP25", "SNAP29", "SNAP47", "YKT6",
        "BET1", "BET1L", "GOSR1", "GOSR2", "USE1", "VTI1A", "VTI1B", "SCFD1",
        "STXBP1", "STXBP2", "STXBP3", "STXBP4", "STXBP5", "STXBP6", "NAPA", "NAPB",
    ],
    "TRAPP": ["TRAPPC1", "TRAPPC2", "TRAPPC2L", "TRAPPC3", "TRAPPC4", "TRAPPC5", "TRAPPC6A", "TRAPPC6B", "TRAPPC8", "TRAPPC9", "TRAPPC10", "TRAPPC11", "TRAPPC12", "TRAPPC13", "TRAPPC14"],
    "COG_tether": ["COG1", "COG2", "COG3", "COG4", "COG5", "COG6", "COG7", "COG8"],
    "GARP_EARP": ["VPS51", "VPS52", "VPS53", "VPS54", "VPS50"],
    "HOPS": ["VPS11", "VPS16", "VPS18", "VPS33A", "VPS33B", "VPS39", "VPS41", "VPS8"],
    "Exocyst": ["EXOC1", "EXOC2", "EXOC3", "EXOC4", "EXOC5", "EXOC6", "EXOC6B", "EXOC7", "EXOC8"],
    "Golgin": ["GOLGA1", "GOLGA2", "GOLGA3", "GOLGA4", "GOLGA5", "GOLGA7", "GOLGA8A", "GMAP", "GOLPH3", "GOLPH3L", "GOLM1", "GORASP1", "GORASP2"],
    "ER_contact": ["VAPA", "VAPB", "OSBPL1A", "OSBPL2", "OSBPL5", "OSBPL6", "OSBPL8", "OSBPL9", "OSBPL10", "OSBPL11", "ESYT1", "ESYT2", "ESYT3", "RTN1", "RTN3", "RTN4", "ATL1", "ATL2", "ATL3"],
    "ER_retention": ["KDELR1", "KDELR2", "KDELR3", "ERP44", "RER1", "SURF4"],
    "Arf_effectors_GEF": ["GBF1", "ARFGEF1", "ARFGEF2", "ARFGAP1", "ARFGAP2", "ARFGAP3", "ARFIP1", "ARFIP2", "PICALM", "AP1G1", "AP1B1", "AP1M1", "AP1S1", "AP1S2", "AP3B1", "AP3D1", "AP3M1", "AP3M2", "AP3S1", "AP3S2"],
}

# flatten with family
GENE_TO_FAMILY = {}
for fam, genes in TRAFFICKING.items():
    for g in genes:
        GENE_TO_FAMILY.setdefault(g, fam)

ALL_GENES = sorted(GENE_TO_FAMILY)
RAB_GENES = sorted(TRAFFICKING["RAB"] + TRAFFICKING["RAB_effector_GEF"] + TRAFFICKING["RAB_GAP"])
