// Canonical reference variants shown in the "Known variants" tab. Each is a
// well-established substitution on a bundled transcript, so it resolves offline.
// Every entry is checked by src/tests/knownVariants.test.mjs: the resolver must
// reproduce `expected` with the reference base matching. If you add one and the
// test fails, the c. position or transcript version is wrong — fix it, don't ship it.

const RAW = [
  // --- Canonical haem-onc hotspots (VV-verified; genomic coords, exon and dbSNP from Chris's reference table) ---
  { gene: 'JAK2', hgvsC: 'NM_004972.4:c.1849G>T', expected: 'V617F', hgvsP: 'p.Val617Phe', consequence: 'missense', exon: 14, dbSNP: 'rs77375493', chrPos: '9:5,073,770', genomicRefAlt: 'G>T', category: 'Haem-onc hotspots', note: 'MPN driver (PV, ET, PMF).' },
  { gene: 'MYD88', hgvsC: 'NM_002468.5:c.755T>C', expected: 'L252P', hgvsP: 'p.Leu252Pro', consequence: 'missense', exon: 5, dbSNP: 'rs387907272', chrPos: '3:38,141,150', genomicRefAlt: 'T>C', category: 'Haem-onc hotspots', note: 'WM, ABC-DLBCL. Equals legacy L265P on NM_002468.4.' },
  { gene: 'BRAF', hgvsC: 'NM_004333.6:c.1799T>A', expected: 'V600E', hgvsP: 'p.Val600Glu', consequence: 'missense', exon: 15, dbSNP: 'rs113488022', chrPos: '7:140,753,336', genomicRefAlt: 'A>T', category: 'Haem-onc hotspots', note: 'HCL; melanoma (minus strand).' },
  { gene: 'IDH1', hgvsC: 'NM_005896.4:c.395G>A', expected: 'R132H', hgvsP: 'p.Arg132His', consequence: 'missense', exon: 4, dbSNP: 'rs121913500', chrPos: '2:208,248,388', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'AML (CIMP+); gliomas (minus strand).' },
  { gene: 'IDH2', hgvsC: 'NM_002168.4:c.419G>A', expected: 'R140Q', hgvsP: 'p.Arg140Gln', consequence: 'missense', exon: 4, dbSNP: 'rs121913502', chrPos: '15:90,088,702', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'AML (minus strand).' },
  { gene: 'IDH2', hgvsC: 'NM_002168.4:c.515G>A', expected: 'R172K', hgvsP: 'p.Arg172Lys', consequence: 'missense', exon: 4, dbSNP: 'rs121913503', chrPos: '15:90,088,606', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'AML / AITL (worse prognosis); T-cell BAM check (minus strand).' },
  { gene: 'KIT', hgvsC: 'NM_000222.3:c.2447A>T', expected: 'D816V', hgvsP: 'p.Asp816Val', consequence: 'missense', exon: 17, dbSNP: 'rs121913507', chrPos: '4:54,733,155', genomicRefAlt: 'A>T', category: 'Haem-onc hotspots', note: 'Systemic mastocytosis; AML.' },
  { gene: 'NPM1', hgvsC: 'NM_002520.7:c.860_863dup', expected: 'W288Cfs*12', hgvsP: 'p.Trp288CysfsTer12', consequence: 'frameshift', frameshift: true, exon: 11, dbSNP: 'rs587776806', chrPos: '5:171,410,539', genomicRefAlt: 'C>CTCTG', category: 'Haem-onc hotspots', note: 'AML — type A 4-bp TCTG duplication.' },
  { gene: 'FLT3', hgvsC: 'NM_004119.3:c.2503G>T', expected: 'D835Y', hgvsP: 'p.Asp835Tyr', consequence: 'missense', exon: 20, dbSNP: null, chrPos: '13:28,018,505', genomicRefAlt: 'C>A', category: 'Haem-onc hotspots', note: 'AML (TKD) (minus strand).' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.524G>A', expected: 'R175H', hgvsP: 'p.Arg175His', consequence: 'missense', exon: 5, dbSNP: 'rs28934578', chrPos: '17:7,675,088', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'Pan-cancer (minus strand).' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.743G>A', expected: 'R248Q', hgvsP: 'p.Arg248Gln', consequence: 'missense', exon: 7, dbSNP: 'rs11540652', chrPos: '17:7,674,220', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'Pan-cancer (minus strand).' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.818G>A', expected: 'R273H', hgvsP: 'p.Arg273His', consequence: 'missense', exon: 8, dbSNP: 'rs28934576', chrPos: '17:7,673,802', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'Pan-cancer (minus strand).' },
  { gene: 'DNMT3A', hgvsC: 'NM_022552.5:c.2645G>A', expected: 'R882H', hgvsP: 'p.Arg882His', consequence: 'missense', exon: 23, dbSNP: 'rs147001633', chrPos: '2:25,234,373', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'AML (minus strand).' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.35G>A', expected: 'G12D', hgvsP: 'p.Gly12Asp', consequence: 'missense', exon: 2, dbSNP: 'rs121913529', chrPos: '12:25,245,350', genomicRefAlt: 'C>T', category: 'Haem-onc hotspots', note: 'Pan-cancer (minus strand).' },
  { gene: 'NRAS', hgvsC: 'NM_002524.5:c.182A>G', expected: 'Q61R', hgvsP: 'p.Gln61Arg', consequence: 'missense', exon: 3, dbSNP: null, chrPos: '1:114,713,908', genomicRefAlt: 'T>C', category: 'Haem-onc hotspots', note: 'AML, T-ALL (minus strand).' },
  { gene: 'SF3B1', hgvsC: 'NM_012433.4:c.2098A>G', expected: 'K700E', hgvsP: 'p.Lys700Glu', consequence: 'missense', exon: 15, dbSNP: 'rs559063155', chrPos: '2:197,402,110', genomicRefAlt: 'T>C', category: 'Haem-onc hotspots', note: 'MDS-RS, CLL (minus strand).' },
  { gene: 'SRSF2', hgvsC: 'NM_003016.5:c.284C>A', expected: 'P95H', hgvsP: 'p.Pro95His', consequence: 'missense', exon: 1, dbSNP: 'rs751713049', chrPos: '17:76,736,877', genomicRefAlt: 'G>T', category: 'Haem-onc hotspots', note: 'MDS (minus strand).' },
  { gene: 'U2AF1', hgvsC: 'NM_006758.3:c.101C>T', expected: 'S34F', hgvsP: 'p.Ser34Phe', consequence: 'missense', exon: 2, dbSNP: 'rs371769427', chrPos: '21:43,104,346', genomicRefAlt: 'G>A', category: 'Haem-onc hotspots', note: 'MDS (minus strand).' },
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.1936T>A', expected: 'Y646N', hgvsP: 'p.Tyr646Asn', consequence: 'missense', exon: 16, dbSNP: null, chrPos: '7:148,811,636', genomicRefAlt: 'A>T', category: 'Haem-onc hotspots', note: 'FL, DLBCL (minus strand).' },
  { gene: 'ASXL1', hgvsC: 'NM_015338.5:c.1934dup', expected: 'G646Wfs*12', hgvsP: 'p.Gly646TrpfsTer12', consequence: 'frameshift', frameshift: true, exon: 13, dbSNP: 'rs750318549', chrPos: '20:32,434,638', genomicRefAlt: 'A>AG', category: 'Haem-onc hotspots', note: 'AML/MDS frameshift hotspot (artefact-prone <10% VAF).' },
  { gene: 'CSF3R', hgvsC: 'NM_000760.4:c.1853C>T', expected: 'T618I', hgvsP: 'p.Thr618Ile', consequence: 'missense', exon: 14, dbSNP: null, chrPos: '1:36,467,833', genomicRefAlt: 'G>A', category: 'Haem-onc hotspots', note: 'CNL (minus strand).' },
  { gene: 'CBL', hgvsC: 'NM_005188.4:c.1111T>C', expected: 'Y371H', hgvsP: 'p.Tyr371His', consequence: 'missense', exon: 8, dbSNP: 'rs267606706', chrPos: '11:119,278,181', genomicRefAlt: 'T>C', category: 'Haem-onc hotspots', note: 'JMML, CMML.' },
  { gene: 'RHOA', hgvsC: 'NM_001664.4:c.50G>T', expected: 'G17V', hgvsP: 'p.Gly17Val', consequence: 'missense', exon: 2, dbSNP: null, chrPos: '3:49,375,540', genomicRefAlt: 'C>A', category: 'Haem-onc hotspots', note: 'AITL / T-NHL driver (minus strand); mandatory T-cell BAM check, pairs with IDH2 R172.' },

  // --- Lymphoma / lymphoid ---
  { gene: 'MYD88', hgvsC: 'NM_002468.4:c.794T>C', expected: 'L265P', category: 'Lymphoma', note: 'ABC-DLBCL and Waldenström macroglobulinaemia; drives NF-κB signalling.' },
  { gene: 'CD79B', hgvsC: 'NM_000626.4:c.587A>G', expected: 'Y196C', category: 'Lymphoma', note: 'ABC-DLBCL; disrupts the B-cell receptor ITAM.' },
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.1936T>A', expected: 'Y646N', category: 'Lymphoma', note: 'GCB-DLBCL and follicular lymphoma; SET-domain gain of function.' },
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.1937A>T', expected: 'Y646F', category: 'Lymphoma', note: 'GCB-DLBCL and follicular lymphoma; same hotspot codon as Y646N.' },
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.1937A>C', expected: 'Y646S', category: 'Lymphoma', note: 'GCB-DLBCL and follicular lymphoma; same hotspot codon.' },
  { gene: 'BRAF', hgvsC: 'NM_004333.6:c.1799T>A', expected: 'V600E', category: 'Lymphoma', note: 'Hairy cell leukaemia (near-universal); also many solid tumours.' },
  { gene: 'XPO1', hgvsC: 'NM_003400.4:c.1711G>A', expected: 'E571K', category: 'Lymphoma', note: 'Primary mediastinal B-cell and Hodgkin lymphoma; nuclear export.' },
  { gene: 'SF3B1', hgvsC: 'NM_012433.4:c.2098A>G', expected: 'K700E', category: 'Lymphoma', note: 'Chronic lymphocytic leukaemia and MDS; splicing-factor hotspot.' },

  // --- TP53 hotspots ---
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.818G>A', expected: 'R273H', category: 'TP53 hotspots', note: 'DNA-contact hotspot.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.817C>T', expected: 'R273C', category: 'TP53 hotspots', note: 'DNA-contact hotspot (same codon as R273H).' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.524G>A', expected: 'R175H', category: 'TP53 hotspots', note: 'Structural hotspot.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.742C>T', expected: 'R248W', category: 'TP53 hotspots', note: 'DNA-contact hotspot.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.743G>A', expected: 'R248Q', category: 'TP53 hotspots', note: 'DNA-contact hotspot (same codon as R248W).' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.844C>T', expected: 'R282W', category: 'TP53 hotspots', note: 'Structural/contact hotspot.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.659A>G', expected: 'Y220C', category: 'TP53 hotspots', note: 'Structural hotspot; creates a druggable pocket.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.733G>A', expected: 'G245S', category: 'TP53 hotspots', note: 'Structural hotspot.' },

  // --- General oncogenes ---
  // KRAS is on NM_004985.5 (MANE Select) throughout, matching the haem-onc block above and the
  // imported whitelist. These sat on NM_033360.4 (MANE Plus Clinical) with no coordinate, so the
  // merged list could not tell they were the same molecules as the whitelist rows and showed each
  // hotspot twice. Coordinates are computed from the bundled coding-exon table and checked by
  // src/tests/knownVariants.test.mjs.
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.35G>A', expected: 'G12D', chrPos: '12:25,245,350', genomicRefAlt: 'C>T', category: 'General oncogene', note: 'Codon-12 RAS driver.' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.35G>T', expected: 'G12V', chrPos: '12:25,245,350', genomicRefAlt: 'C>A', category: 'General oncogene', note: 'Codon-12 RAS driver.' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.34G>T', expected: 'G12C', chrPos: '12:25,245,351', genomicRefAlt: 'C>A', category: 'General oncogene', note: 'Codon-12 driver; targetable in lung cancer.' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.38G>A', expected: 'G13D', chrPos: '12:25,245,347', genomicRefAlt: 'C>T', category: 'General oncogene', note: 'Codon-13 RAS driver.' },
  { gene: 'NRAS', hgvsC: 'NM_002524.5:c.182A>G', expected: 'Q61R', category: 'General oncogene', note: 'Codon-61 RAS driver; melanoma and AML.' },
  { gene: 'NRAS', hgvsC: 'NM_002524.5:c.181C>A', expected: 'Q61K', category: 'General oncogene', note: 'Codon-61 RAS driver.' },
  { gene: 'NRAS', hgvsC: 'NM_002524.5:c.35G>A', expected: 'G12D', category: 'General oncogene', note: 'Codon-12 RAS driver.' },
  { gene: 'IDH1', hgvsC: 'NM_005896.4:c.395G>A', expected: 'R132H', category: 'General oncogene', note: 'Neomorphic (2-HG) hotspot; glioma and AML.' },
  { gene: 'IDH1', hgvsC: 'NM_005896.4:c.394C>T', expected: 'R132C', category: 'General oncogene', note: 'Neomorphic hotspot (same codon).' },
  { gene: 'JAK2', hgvsC: 'NM_004972.4:c.1849G>T', expected: 'V617F', category: 'General oncogene', note: 'Myeloproliferative-neoplasm driver.' },
  { gene: 'PTEN', hgvsC: 'NM_000314.8:c.389G>A', expected: 'R130Q', category: 'General oncogene', note: 'Phosphatase-domain hotspot (tumour suppressor).' },
  { gene: 'PTEN', hgvsC: 'NM_000314.8:c.388C>T', expected: 'R130*', category: 'General oncogene', note: 'Nonsense at the same codon; loss of function.' },
  { gene: 'ARID1A', hgvsC: 'NM_006015.6:c.461A>C', expected: 'Y154S', category: 'General oncogene', note: 'Chromatin remodeller (tumour suppressor); the app’s worked example.' },

  // --- expansion (each verified by the test) ---
  // Lymphoma
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.1936T>C', expected: 'Y646H', category: 'Lymphoma', note: 'GCB-DLBCL / FL; Y646 hotspot allele.' },
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.1937A>G', expected: 'Y646C', category: 'Lymphoma', note: 'GCB-DLBCL / FL; Y646 hotspot allele.' },
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.2045C>G', expected: 'A682G', category: 'Lymphoma', note: 'SET-domain gain-of-function; second EZH2 hotspot.' },
  { gene: 'EZH2', hgvsC: 'NM_004456.5:c.2075C>T', expected: 'A692V', category: 'Lymphoma', note: 'SET-domain gain-of-function.' },
  { gene: 'CD79B', hgvsC: 'NM_000626.4:c.586T>C', expected: 'Y196H', category: 'Lymphoma', note: 'ABC-DLBCL ITAM hotspot allele.' },
  { gene: 'CD79B', hgvsC: 'NM_000626.4:c.586T>A', expected: 'Y196N', category: 'Lymphoma', note: 'ABC-DLBCL ITAM hotspot allele.' },
  { gene: 'CD79B', hgvsC: 'NM_000626.4:c.587A>C', expected: 'Y196S', category: 'Lymphoma', note: 'ABC-DLBCL ITAM hotspot allele.' },
  { gene: 'XPO1', hgvsC: 'NM_003400.4:c.1712A>G', expected: 'E571G', category: 'Lymphoma', note: 'Alternative XPO1 E571 allele.' },
  { gene: 'SF3B1', hgvsC: 'NM_012433.4:c.1874G>A', expected: 'R625H', category: 'Lymphoma', note: 'Splicing-factor hotspot; CLL / MDS.' },
  { gene: 'SF3B1', hgvsC: 'NM_012433.4:c.1984C>T', expected: 'H662Y', category: 'Lymphoma', note: 'Splicing-factor hotspot.' },

  // TP53 hotspots
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.747G>T', expected: 'R249S', category: 'TP53 hotspots', note: 'DNA-contact hotspot.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.734G>A', expected: 'G245D', category: 'TP53 hotspots', note: 'Structural hotspot.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.733G>T', expected: 'G245C', category: 'TP53 hotspots', note: 'Structural hotspot (same codon as G245S/D).' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.637C>T', expected: 'R213*', category: 'TP53 hotspots', note: 'Recurrent nonsense; loss of function.' },
  { gene: 'TP53', hgvsC: 'NM_000546.6:c.586C>T', expected: 'R196*', category: 'TP53 hotspots', note: 'Recurrent nonsense; loss of function.' },

  // General oncogenes
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.35G>C', expected: 'G12A', chrPos: '12:25,245,350', genomicRefAlt: 'C>G', category: 'General oncogene', note: 'Codon-12 driver.' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.34G>A', expected: 'G12S', chrPos: '12:25,245,351', genomicRefAlt: 'C>T', category: 'General oncogene', note: 'Codon-12 driver.' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.34G>C', expected: 'G12R', chrPos: '12:25,245,351', genomicRefAlt: 'C>G', category: 'General oncogene', note: 'Codon-12 driver.' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.183A>T', expected: 'Q61H', chrPos: '12:25,227,341', genomicRefAlt: 'T>A', category: 'General oncogene', note: 'Codon-61 driver.' },
  { gene: 'KRAS', hgvsC: 'NM_004985.5:c.436G>A', expected: 'A146T', chrPos: '12:25,225,628', genomicRefAlt: 'C>T', category: 'General oncogene', note: 'Codon-146 driver.' },
  { gene: 'NRAS', hgvsC: 'NM_002524.5:c.34G>A', expected: 'G12S', category: 'General oncogene', note: 'Codon-12 driver.' },
  { gene: 'NRAS', hgvsC: 'NM_002524.5:c.182A>T', expected: 'Q61L', category: 'General oncogene', note: 'Codon-61 driver.' },
  { gene: 'IDH1', hgvsC: 'NM_005896.4:c.394C>G', expected: 'R132G', category: 'General oncogene', note: 'Neomorphic hotspot allele.' },
  { gene: 'IDH1', hgvsC: 'NM_005896.4:c.395G>T', expected: 'R132L', category: 'General oncogene', note: 'Neomorphic hotspot allele.' },
  { gene: 'IDH1', hgvsC: 'NM_005896.4:c.394C>A', expected: 'R132S', category: 'General oncogene', note: 'Neomorphic hotspot allele.' },
  { gene: 'BRAF', hgvsC: 'NM_004333.6:c.1801A>G', expected: 'K601E', category: 'General oncogene', note: 'Activating BRAF variant outside V600.' },
  { gene: 'BRAF', hgvsC: 'NM_004333.6:c.1406G>C', expected: 'G469A', category: 'General oncogene', note: 'Activating BRAF kinase-domain variant.' },
  { gene: 'PTEN', hgvsC: 'NM_000314.8:c.697C>T', expected: 'R233*', category: 'General oncogene', note: 'Recurrent nonsense; loss of function.' },
  { gene: 'PTEN', hgvsC: 'NM_000314.8:c.517C>T', expected: 'R173C', category: 'General oncogene', note: 'Phosphatase-domain missense.' },
];

// Deduplicate: the haem-onc set is listed first, so it wins. A later entry for the
// same gene + protein change (e.g. BRAF V600E already in the haem-onc set) is dropped
// so each variant appears once. Same-gene entries with a different protein change,
// or the same change on a different transcript version (MYD88 L252P vs L265P), are kept.
const seen = new Set();
export const KNOWN_VARIANTS = RAW.filter((v) => {
  const key = `${v.gene}|${v.expected}`;
  if (seen.has(key)) return false;
  seen.add(key);
  return true;
});

/** Distinct categories in list order. */
export function knownCategories() {
  const seen = [];
  for (const v of KNOWN_VARIANTS) if (!seen.includes(v.category)) seen.push(v.category);
  return seen;
}
