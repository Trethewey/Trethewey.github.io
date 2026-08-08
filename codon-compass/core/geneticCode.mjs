// Standard genetic code (NCBI translation table 1), DNA sense strand (T not U).
// This table is the single source of truth for translation. The control test in
// src/tests translates a whole RefSeq coding sequence with it and requires the
// result to equal the published protein exactly, so any error here fails a test.

export const STANDARD_CODE = Object.freeze({
  TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L',
  CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
  ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M',
  GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
  TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S',
  CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*',
  CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
  AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K',
  GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
  TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W',
  CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
  AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
  GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G',
});

export const STOP_CODONS = Object.freeze(['TAA', 'TAG', 'TGA']);
export const START_CODON = 'ATG';
export const BASES = Object.freeze(['T', 'C', 'A', 'G']);

/** Return the one-letter amino acid for a codon, or null if the codon is not three valid bases. */
export function translateCodon(codon) {
  if (typeof codon !== 'string') return null;
  const c = codon.toUpperCase();
  return Object.prototype.hasOwnProperty.call(STANDARD_CODE, c) ? STANDARD_CODE[c] : null;
}

/**
 * Translate a coding sequence into a one-letter protein string.
 * @param {string} cds  coding sequence (should start at the A of ATG)
 * @param {object} [opts]
 * @param {boolean} [opts.stopAtStop=true]  stop at the first stop codon (do not include '*')
 * @returns {{protein:string, stoppedAtStop:boolean, stopIndex:number, ambiguousCodons:number}}
 */
export function translate(cds, opts = {}) {
  const stopAtStop = opts.stopAtStop !== false;
  const seq = String(cds || '').toUpperCase().replace(/\s/g, '');
  let protein = '';
  let stoppedAtStop = false;
  let stopIndex = -1;
  let ambiguousCodons = 0;
  for (let i = 0; i + 3 <= seq.length; i += 3) {
    const codon = seq.slice(i, i + 3);
    const aa = translateCodon(codon);
    if (aa === null) { protein += 'X'; ambiguousCodons += 1; continue; }
    if (aa === '*') {
      stoppedAtStop = true;
      stopIndex = i / 3;
      if (stopAtStop) break;
      protein += '*';
      continue;
    }
    protein += aa;
  }
  return { protein, stoppedAtStop, stopIndex, ambiguousCodons };
}

/** Reverse complement of a DNA string (non-ACGT characters are passed through unchanged, uppercased). */
export function reverseComplement(dna) {
  const comp = { A: 'T', T: 'A', C: 'G', G: 'C', N: 'N' };
  return String(dna || '')
    .toUpperCase()
    .split('')
    .reverse()
    .map((b) => comp[b] || b)
    .join('');
}

/** DNA sense strand -> mRNA (T -> U). */
export function toMrna(dna) {
  return String(dna || '').toUpperCase().replace(/T/g, 'U');
}
