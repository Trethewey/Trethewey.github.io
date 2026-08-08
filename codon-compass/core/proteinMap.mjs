// Map a residue number between the MANE/RefSeq protein and the UniProt canonical
// protein. Domains are stored in MANE coordinates, so the lollipop needs no
// mapping; but the AlphaFold structure is numbered by UniProt, so the 3D highlight
// maps the MANE residue into UniProt numbering.
//
// What the map holds. `uniprotMap` on each transcript records only four numbers,
// worked out at build time by comparing the two sequences from each end:
//   uLen, mLen  the length of the UniProt protein and of the MANE protein
//   prefix      residues 1..prefix are the same residue in both proteins
//   suffix      the last `suffix` residues are the same residue in both proteins
// Anything between those two matching ends is a block where the sequences differ.
// The map does not record the alignment inside that block.
//
// So this module maps a residue exactly in the two matching ends, and inside the
// differing block only when the block is the same length in both proteins (the
// numbering then runs straight through it). In every other case it returns
// `pos: null` with a reason. That is deliberate: the only caller is the 3D
// structure highlight, and a guessed number does not mean "about right" there — it
// lights up a specific, different residue on the model.
//
// Result shape: { pos, exact, offset, reason }
//   pos     UniProt residue number, or null when it cannot be mapped
//   exact   true only where the two sequences are known to match residue for residue
//   offset  UniProt minus MANE numbering shift where one applies, else null
//   reason  plain-English reason when pos is null (or when pos rests on an
//           assumption), else null

const cannotMap = (reason) => ({ pos: null, exact: false, offset: null, reason });

/**
 * Map a residue number from MANE/RefSeq numbering into UniProt numbering.
 *
 * THIS OFTEN RETURNS NO NUMBER. `pos` is null whenever the alignment is not
 * recorded: inside a stretch where the two sequences differ and are of different
 * lengths (KRAS 151-186, EZH2 298-303), for a residue that exists only in the MANE
 * protein (MYD88 2-14), for a number that is missing, not whole, or past the end of
 * the protein, and for a broken map. There is no fallback number by design — the
 * caller must test `pos === null` and show `reason` instead, because a guessed
 * number is not "about right" on a 3D model: it lights up a different residue.
 *
 * When `pos` is a number: `exact: true` means the two sequences are known to match
 * residue for residue there, and `reason` is null. `exact: false` means the number
 * rests on a stated assumption (no map recorded for this transcript, or a differing
 * stretch that is the same length in both) — `reason` says which, and the caller
 * must pass that caveat on rather than present the number as checked.
 *
 * @param {{uLen:number,mLen:number,prefix:number,suffix:number,identical:boolean}|null|undefined} map
 * @param {number} manePos residue number in MANE/RefSeq numbering, 1-based
 * @returns {{pos:number|null, exact:boolean, offset:number|null, reason:string|null}}
 */
export function maneToUniprot(map, manePos) {
  const pos = typeof manePos === 'number' ? manePos : Number(manePos);
  if (!Number.isInteger(pos) || pos < 1) {
    return cannotMap('no residue number to map (it must be a whole number of 1 or more)');
  }
  if (!map) {
    // Nothing was recorded for this transcript. Fall back to the same numbering,
    // but say that it is an assumption rather than a checked fact.
    return { pos, exact: false, offset: 0, reason: 'this transcript has no numbering map, so the two proteins are assumed to be numbered the same' };
  }
  if (map.identical) return { pos, exact: true, offset: 0, reason: null };

  const { uLen, mLen, prefix, suffix } = map;
  const sound = [uLen, mLen, prefix, suffix].every((n) => Number.isInteger(n) && n >= 0) && uLen >= 1 && mLen >= 1;
  if (!sound) return cannotMap('the numbering map for this transcript is incomplete');
  if (pos > mLen) return cannotMap(`the MANE protein is only ${mLen} residues long`);

  const sameLength = uLen === mLen;
  // Matching ends that overlap can only be read one way when the proteins are the
  // same length; otherwise the two ends disagree about the shift and neither wins.
  if (!sameLength && (prefix + suffix > mLen || prefix + suffix > uLen)) {
    return cannotMap('the numbering map is inconsistent (its matching ends overlap)');
  }

  // Matching start: same residue, same number.
  if (pos <= prefix) return { pos, exact: true, offset: 0, reason: null };

  // Matching end: same residue, shifted by the difference in length.
  if (pos > mLen - suffix) {
    const uniprotPos = pos - (mLen - uLen);
    if (uniprotPos < 1 || uniprotPos > uLen) return cannotMap('the numbering map is inconsistent (the matching end falls outside the UniProt protein)');
    return { pos: uniprotPos, exact: true, offset: uLen - mLen, reason: null };
  }

  // The differing block: MANE prefix+1 .. mLen-suffix against UniProt prefix+1 .. uLen-suffix.
  const uniprotBlock = uLen - suffix - prefix;
  const maneBlock = mLen - suffix - prefix;
  if (uniprotBlock <= 0) {
    return cannotMap('these residues are extra in the MANE protein — the UniProt protein has no residue here');
  }
  if (uniprotBlock === maneBlock) {
    // Same number of residues in both, so the numbering runs straight through.
    // Still not called exact: the residues themselves differ, and an insertion and
    // a deletion of the same size inside the block would cancel out unseen.
    return { pos, exact: false, offset: 0, reason: 'inside a stretch where the two sequences differ; it is the same length in both, so the numbering is taken to run straight through' };
  }
  return cannotMap(`inside a stretch where the two sequences differ and are of different lengths (${maneBlock} residues in MANE, ${uniprotBlock} in UniProt) — which residue matches which is not recorded`);
}

/**
 * True when MANE and UniProt numbering coincide for this transcript. Same-length
 * proteins share their numbering unless an insertion and a deletion of the same
 * size cancel out inside the differing block, which this map cannot show. A
 * missing map is unknown, not a match.
 * (No caller today; kept as the counterpart to maneToUniprot.)
 */
export function numberingMatches(map) {
  if (!map) return false;
  if (map.identical) return true;
  return map.uLen === map.mLen;
}
