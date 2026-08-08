// Amino-acid reference data for Codon Compass.
//
// Polarity uses the common textbook four-way grouping:
//   nonpolar | polar uncharged | acidic (negatively charged) | basic (positively charged)
// Aromatic means the side chain contains an aromatic ring — a flat, fully conjugated
// ring. That is histidine, phenylalanine, tryptophan and tyrosine. Many teaching tables
// name only phenylalanine, tryptophan and tyrosine as "the aromatic amino acids" because
// histidine's imidazole ring is also basic and far less hydrophobic, but the ring is
// aromatic all the same, and the shipped structures say so: the generated structure file
// writes histidine's ring atoms as aromatic. Proline's ring is the only other side-chain
// ring, and it is saturated, so it is not aromatic. AROMATICITY_NOTE below states this
// for the reader; the control test checks the flags against the shipped structures.
// Charge is the net side-chain charge at physiological pH (~7.4).
// Hydropathy is the Kyte-Doolittle index. avgResidueMass is the average residue
// mass in Daltons (a residue in a chain, i.e. free amino acid minus one water).
// sideChainPka is given only for ionizable side chains.
//
// These property values are cross-checked by an independent verification workflow;
// the classification scheme is named in the UI so the reader knows the convention.

import { STANDARD_CODE } from './geneticCode.mjs';

const POLARITY_COLOUR = Object.freeze({
  nonpolar: '#d99a4e',        // sand / amber
  'polar uncharged': '#37a89a', // teal
  acidic: '#e05561',          // red
  basic: '#4f83e2',           // blue
});

// name/property table. Codons are derived from the genetic code below so they
// can never drift out of sync with the translation table.
const AA = [
  { one: 'A', three: 'Ala', name: 'Alanine',       polarity: 'nonpolar',        charge: 0,   aromatic: false, hydropathy: 1.8,  mass: 71.0788,  pka: null,  essential: false, sideChain: 'Methyl group; small and hydrophobic.' },
  { one: 'R', three: 'Arg', name: 'Arginine',      polarity: 'basic',           charge: 1,   aromatic: false, hydropathy: -4.5, mass: 156.1875, pka: 12.48, essential: false, sideChain: 'Guanidinium group; positively charged, strong base.' },
  { one: 'N', three: 'Asn', name: 'Asparagine',    polarity: 'polar uncharged', charge: 0,   aromatic: false, hydropathy: -3.5, mass: 114.1038, pka: null,  essential: false, sideChain: 'Amide group; polar, forms hydrogen bonds.' },
  { one: 'D', three: 'Asp', name: 'Aspartate',     polarity: 'acidic',          charge: -1,  aromatic: false, hydropathy: -3.5, mass: 115.0886, pka: 3.65,  essential: false, sideChain: 'Carboxylate group; negatively charged.' },
  { one: 'C', three: 'Cys', name: 'Cysteine',      polarity: 'polar uncharged', charge: 0,   aromatic: false, hydropathy: 2.5,  mass: 103.1388, pka: 8.3,   essential: false, sideChain: 'Thiol group; forms disulfide bonds.' },
  { one: 'Q', three: 'Gln', name: 'Glutamine',     polarity: 'polar uncharged', charge: 0,   aromatic: false, hydropathy: -3.5, mass: 128.1307, pka: null,  essential: false, sideChain: 'Amide group; polar, forms hydrogen bonds.' },
  { one: 'E', three: 'Glu', name: 'Glutamate',     polarity: 'acidic',          charge: -1,  aromatic: false, hydropathy: -3.5, mass: 129.1155, pka: 4.25,  essential: false, sideChain: 'Carboxylate group; negatively charged.' },
  { one: 'G', three: 'Gly', name: 'Glycine',       polarity: 'nonpolar',        charge: 0,   aromatic: false, hydropathy: -0.4, mass: 57.0519,  pka: null,  essential: false, sideChain: 'Single hydrogen; smallest, very flexible.' },
  { one: 'H', three: 'His', name: 'Histidine',     polarity: 'basic',           charge: 0.1, aromatic: true,  hydropathy: -3.2, mass: 137.1411, pka: 6.0,   essential: true,  sideChain: 'Imidazole ring; aromatic; weakly basic, near neutral at pH 7.4.' },
  { one: 'I', three: 'Ile', name: 'Isoleucine',    polarity: 'nonpolar',        charge: 0,   aromatic: false, hydropathy: 4.5,  mass: 113.1594, pka: null,  essential: true,  sideChain: 'Branched hydrocarbon; strongly hydrophobic.' },
  { one: 'L', three: 'Leu', name: 'Leucine',       polarity: 'nonpolar',        charge: 0,   aromatic: false, hydropathy: 3.8,  mass: 113.1594, pka: null,  essential: true,  sideChain: 'Branched hydrocarbon; strongly hydrophobic.' },
  { one: 'K', three: 'Lys', name: 'Lysine',        polarity: 'basic',           charge: 1,   aromatic: false, hydropathy: -3.9, mass: 128.1741, pka: 10.53, essential: true,  sideChain: 'Amino group; positively charged.' },
  { one: 'M', three: 'Met', name: 'Methionine',    polarity: 'nonpolar',        charge: 0,   aromatic: false, hydropathy: 1.9,  mass: 131.1926, pka: null,  essential: true,  sideChain: 'Thioether; hydrophobic; the start codon residue.' },
  { one: 'F', three: 'Phe', name: 'Phenylalanine', polarity: 'nonpolar',        charge: 0,   aromatic: true,  hydropathy: 2.8,  mass: 147.1766, pka: null,  essential: true,  sideChain: 'Benzyl ring; aromatic and hydrophobic.' },
  { one: 'P', three: 'Pro', name: 'Proline',       polarity: 'nonpolar',        charge: 0,   aromatic: false, hydropathy: -1.6, mass: 97.1167,  pka: null,  essential: false, sideChain: 'Ring to the backbone; rigid, breaks helices.' },
  { one: 'S', three: 'Ser', name: 'Serine',        polarity: 'polar uncharged', charge: 0,   aromatic: false, hydropathy: -0.8, mass: 87.0782,  pka: null,  essential: false, sideChain: 'Hydroxyl group; polar; common phosphorylation site.' },
  { one: 'T', three: 'Thr', name: 'Threonine',     polarity: 'polar uncharged', charge: 0,   aromatic: false, hydropathy: -0.7, mass: 101.1051, pka: null,  essential: true,  sideChain: 'Hydroxyl group; polar; common phosphorylation site.' },
  { one: 'W', three: 'Trp', name: 'Tryptophan',    polarity: 'nonpolar',        charge: 0,   aromatic: true,  hydropathy: -0.9, mass: 186.2132, pka: null,  essential: true,  sideChain: 'Indole ring; largest; aromatic.' },
  { one: 'Y', three: 'Tyr', name: 'Tyrosine',      polarity: 'polar uncharged', charge: 0,   aromatic: true,  hydropathy: -1.3, mass: 163.1760, pka: 10.07, essential: false, sideChain: 'Phenol ring; aromatic; common phosphorylation site.' },
  { one: 'V', three: 'Val', name: 'Valine',        polarity: 'nonpolar',        charge: 0,   aromatic: false, hydropathy: 4.2,  mass: 99.1326,  pka: null,  essential: true,  sideChain: 'Branched hydrocarbon; strongly hydrophobic.' },
];

// Shown next to the Aromatic column so the reader knows which grouping is in use.
export const AROMATICITY_NOTE = 'Aromatic means the side chain carries a flat ring whose electrons are shared right around it: '
  + 'histidine, phenylalanine, tryptophan and tyrosine. Many tables list only phenylalanine, tryptophan '
  + 'and tyrosine, because histidine’s imidazole ring is also basic and much less hydrophobic — but it is '
  + 'an aromatic ring.';

// Codons per amino acid, derived from the genetic code (single source of truth).
const codonsByAa = {};
for (const [codon, one] of Object.entries(STANDARD_CODE)) {
  (codonsByAa[one] ||= []).push(codon);
}
for (const list of Object.values(codonsByAa)) list.sort();

export const AMINO_ACIDS = Object.freeze(
  AA.map((a) => Object.freeze({
    ...a,
    colour: POLARITY_COLOUR[a.polarity],
    codons: Object.freeze(codonsByAa[a.one] || []),
  }))
);

export const STOP = Object.freeze({
  one: '*', three: 'Ter', name: 'Stop (termination)',
  polarity: null, colour: '#8a8f99', codons: Object.freeze(codonsByAa['*'] || []),
});

const byOne = new Map(AMINO_ACIDS.map((a) => [a.one, a]));
byOne.set('*', STOP);

/** Look up an amino acid record by its one-letter code ('*' returns the stop record). */
export function aaByOne(one) {
  return byOne.get(String(one || '').toUpperCase()) || null;
}

/** "Tyr", "Methionine", etc. Helpers for display. */
export function three(one) { const a = aaByOne(one); return a ? a.three : '?'; }
export function fullName(one) { const a = aaByOne(one); return a ? a.name : 'Unknown'; }
export function colourOf(one) { const a = aaByOne(one); return a ? a.colour : '#8a8f99'; }

export { POLARITY_COLOUR };
