// Resolve a coding single-nucleotide substitution to its amino-acid consequence,
// given the transcript coding sequence (CDS). CDS-relative position N equals HGVS
// c.N for coding substitutions, so no cross-UTR coordinate maths is needed.

import { translateCodon, translate } from './geneticCode.mjs';
import { three } from './aminoAcids.mjs';

/**
 * @param {string} cds  the coding sequence, starting at the A of the start codon
 * @param {{position:number, ref:string, alt:string}} sub  a parsed CDS substitution
 * @param {object} [opts] { flank:number } bases of context each side for the base track
 */
export function resolveSubstitution(cds, sub, opts = {}) {
  const flank = Number.isInteger(opts.flank) ? opts.flank : 12;
  const seq = String(cds || '').toUpperCase().replace(/\s/g, '');
  const pos = sub.position; // 1-based
  const ref = String(sub.ref || '').toUpperCase();
  const alt = String(sub.alt || '').toUpperCase();

  const result = {
    ok: false, message: '',
    position: pos, ref, alt,
    cdsLength: seq.length,
    refBaseInTranscript: null, refMatches: null,
    codonNumber: null, posInCodon: null,
    refCodon: null, altCodon: null,
    refAa: null, altAa: null,
    refAaThree: null, altAaThree: null,
    consequence: null, proteinHgvs: null, proteinShort: null,
    window: null,
  };

  if (!seq) { result.message = 'No coding sequence available for this transcript.'; return result; }
  if (!Number.isInteger(pos) || pos < 1 || pos > seq.length) {
    result.message = `Position c.${pos} is outside the coding sequence (length ${seq.length}).`;
    return result;
  }

  const refBaseInTranscript = seq[pos - 1];
  const refMatches = refBaseInTranscript === ref;

  const codonNumber = Math.floor((pos - 1) / 3) + 1;
  const codonStart0 = (codonNumber - 1) * 3;
  const posInCodon = (pos - 1) % 3;
  const refCodon = seq.slice(codonStart0, codonStart0 + 3);
  if (refCodon.length < 3) {
    result.message = 'The affected codon runs past the end of the coding sequence.';
    return result;
  }
  const altArr = refCodon.split('');
  altArr[posInCodon] = alt;
  const altCodon = altArr.join('');

  const refAa = translateCodon(refCodon);
  const altAa = translateCodon(altCodon);

  // Consequence classification.
  let consequence;
  if (refAa === null || altAa === null) {
    consequence = 'unknown';
  } else if (codonNumber === 1 && refAa === 'M' && altAa !== 'M') {
    consequence = 'start-loss';
  } else if (refAa === '*') {
    consequence = altAa === '*' ? 'stop-retained' : 'stop-loss';
  } else if (altAa === '*') {
    consequence = 'nonsense';
  } else if (refAa === altAa) {
    consequence = 'synonymous';
  } else {
    consequence = 'missense';
  }

  const n = codonNumber;
  let proteinHgvs, proteinShort;
  switch (consequence) {
    case 'synonymous':
      proteinHgvs = `p.(${three(refAa)}${n}=)`; proteinShort = `${refAa}${n}=`; break;
    case 'missense':
      proteinHgvs = `p.(${three(refAa)}${n}${three(altAa)})`; proteinShort = `${refAa}${n}${altAa}`; break;
    case 'nonsense':
      proteinHgvs = `p.(${three(refAa)}${n}Ter)`; proteinShort = `${refAa}${n}*`; break;
    case 'start-loss':
      proteinHgvs = `p.(Met1?)`; proteinShort = `M1?`; break;
    case 'stop-loss':
      proteinHgvs = `p.(Ter${n}${three(altAa)}ext*?)`; proteinShort = `*${n}${altAa}ext*?`; break;
    case 'stop-retained':
      proteinHgvs = `p.(Ter${n}=)`; proteinShort = `*${n}=`; break;
    default:
      proteinHgvs = 'p.(?)'; proteinShort = '?';
  }

  // Context window for the base track.
  const from = Math.max(0, pos - 1 - flank);
  const to = Math.min(seq.length, pos - 1 + flank + 1);
  const refWindow = seq.slice(from, to);
  const changeIndex = (pos - 1) - from; // index of the changed base within the window
  const altWindow = refWindow.slice(0, changeIndex) + alt + refWindow.slice(changeIndex + 1);

  Object.assign(result, {
    ok: true,
    message: refMatches ? 'Resolved.' : 'Resolved, but the stated reference base does not match the transcript (check transcript version and strand).',
    refBaseInTranscript, refMatches,
    codonNumber, posInCodon,
    refCodon, altCodon,
    refAa, altAa,
    refAaThree: three(refAa), altAaThree: three(altAa),
    consequence, proteinHgvs, proteinShort,
    window: {
      ref: refWindow, alt: altWindow, changeIndex,
      startCds: from + 1, endCds: to, // 1-based inclusive-ish for display
    },
  });
  return result;
}

/**
 * Resolve a coding indel (deletion, duplication, insertion, delins) to its protein
 * consequence, by applying the edit to the CDS and re-translating.
 *
 * Frameshifts are reported in full HGVS form, e.g. p.(Gly646TrpfsTer12): the first
 * amino acid that differs, the new amino acid there, and the position of the new stop
 * counted from that residue. In-frame events use del / dup / ins / delins with the
 * 3'-most alignment (maximising the common prefix, as HGVS requires).
 *
 * @param {string} cds coding sequence
 * @param {{kind:string,start:number,end:number,insSeq?:string,delSeq?:string,dupSeq?:string}} parsed
 */
export function resolveIndel(cds, parsed, opts = {}) {
  const flank = Number.isInteger(opts.flank) ? opts.flank : 12;
  const seq = String(cds || '').toUpperCase().replace(/\s/g, '');
  const { kind, start, end } = parsed;

  const result = {
    ok: false, message: '', kind, start, end, cdsLength: seq.length,
    refMatches: null, netNt: 0, frameshift: false,
    firstAffected: null, refAa: null, altAa: null, refAaThree: null, altAaThree: null,
    consequence: null, proteinHgvs: null, proteinShort: null,
    insertedNt: '', deletedNt: '', window: null,
  };

  if (!seq) { result.message = 'No coding sequence available for this transcript.'; return result; }
  if (!Number.isInteger(start) || start < 1 || start > seq.length) {
    result.message = `Position c.${start} is outside the coding sequence (length ${seq.length}).`; return result;
  }
  if (!Number.isInteger(end) || end < start || end > seq.length) {
    result.message = `Range c.${start}_${end} is outside the coding sequence (length ${seq.length}).`; return result;
  }

  let mut; let insertedNt = ''; let deletedNt = '';
  if (kind === 'deletion') {
    deletedNt = seq.slice(start - 1, end);
    mut = seq.slice(0, start - 1) + seq.slice(end);
    if (parsed.delSeq) result.refMatches = deletedNt === parsed.delSeq;
  } else if (kind === 'duplication') {
    insertedNt = seq.slice(start - 1, end);
    mut = seq.slice(0, end) + insertedNt + seq.slice(end);
    if (parsed.dupSeq) result.refMatches = insertedNt === parsed.dupSeq;
  } else if (kind === 'insertion') {
    insertedNt = String(parsed.insSeq || '').toUpperCase();
    if (!insertedNt) { result.message = 'Insertion has no inserted sequence.'; return result; }
    mut = seq.slice(0, start) + insertedNt + seq.slice(start);
  } else if (kind === 'delins') {
    deletedNt = seq.slice(start - 1, end);
    insertedNt = String(parsed.insSeq || '').toUpperCase();
    mut = seq.slice(0, start - 1) + insertedNt + seq.slice(end);
  } else {
    result.message = `Event type "${kind}" is not resolved here.`; return result;
  }

  const refProt = translate(seq).protein;
  // A frameshift can run past the original stop codon, so the new stop may lie in the
  // 3'UTR (NPM1 type A is the classic case). Translate the mutant with the UTR appended,
  // otherwise the fsTer count comes out short.
  const utr3 = String(opts.utr3 || '').toUpperCase();
  const mutT = translate(mut + utr3);
  const mutProt = mutT.protein;
  const netNt = mut.length - seq.length;
  const frameshift = (netNt % 3) !== 0;

  // First amino acid that differs (0-based). This also normalises the event to the
  // 3'-most position, matching HGVS.
  let f = 0;
  const minLen = Math.min(refProt.length, mutProt.length);
  while (f < minLen && refProt[f] === mutProt[f]) f += 1;
  const firstAffected = f + 1;

  const three3 = (s) => s.split('').map((c) => three(c)).join('');
  let consequence; let proteinHgvs; let proteinShort;

  if (frameshift) {
    const noStopFound = !mutT.stoppedAtStop; // ran out of sequence before a new stop
    consequence = 'frameshift';
    if (mutProt.length === f && !noStopFound) {
      // the new reading frame hits a stop immediately at this residue
      proteinHgvs = `p.(${three(refProt[f])}${firstAffected}Ter)`;
      proteinShort = `${refProt[f]}${firstAffected}*`;
    } else {
      const newAa = mutProt[f] || '?';
      const terN = mutProt.length - f + 1; // residues in the new frame up to and including the stop
      proteinHgvs = `p.(${three(refProt[f])}${firstAffected}${three(newAa)}fs${noStopFound ? 'Ter?' : `Ter${terN}`})`;
      proteinShort = `${refProt[f]}${firstAffected}${newAa}fs${noStopFound ? '*?' : `*${terN}`}`;
      result.altAa = newAa; result.altAaThree = three(newAa);
      if (noStopFound) result.stopNotFound = true;
    }
    result.refAa = refProt[f] || null;
    result.refAaThree = result.refAa ? three(result.refAa) : null;
  } else {
    // in-frame: longest common prefix (already f), then longest common suffix
    let s = 0;
    const maxS = Math.min(refProt.length - f, mutProt.length - f);
    while (s < maxS && refProt[refProt.length - 1 - s] === mutProt[mutProt.length - 1 - s]) s += 1;
    const refMid = refProt.slice(f, refProt.length - s);
    const mutMid = mutProt.slice(f, mutProt.length - s);
    const startRes = f + 1;
    const endRes = f + refMid.length;

    if (!refMid.length && mutMid.length) {
      const n = mutMid.length;
      const preceding = refProt.slice(f - n, f);
      if (preceding === mutMid) {
        consequence = 'inframe-duplication';
        proteinHgvs = n === 1
          ? `p.(${three(preceding[0])}${f}dup)`
          : `p.(${three(preceding[0])}${f - n + 1}_${three(preceding[n - 1])}${f}dup)`;
        proteinShort = n === 1 ? `${preceding[0]}${f}dup` : `${preceding[0]}${f - n + 1}_${preceding[n - 1]}${f}dup`;
      } else {
        consequence = 'inframe-insertion';
        proteinHgvs = `p.(${three(refProt[f - 1])}${f}_${three(refProt[f])}${f + 1}ins${three3(mutMid)})`;
        proteinShort = `${refProt[f - 1]}${f}_${refProt[f]}${f + 1}ins${mutMid}`;
      }
    } else if (refMid.length && !mutMid.length) {
      consequence = 'inframe-deletion';
      proteinHgvs = refMid.length === 1
        ? `p.(${three(refMid[0])}${startRes}del)`
        : `p.(${three(refMid[0])}${startRes}_${three(refMid[refMid.length - 1])}${endRes}del)`;
      proteinShort = refMid.length === 1 ? `${refMid[0]}${startRes}del` : `${refMid[0]}${startRes}_${refMid[refMid.length - 1]}${endRes}del`;
    } else if (refMid.length === 1 && mutMid.length === 1) {
      // One residue replaced by one residue. At the coding level this is a deletion-insertion — two
      // adjacent bases changed at once, which is how a caller reports it — but at the protein level it
      // is an ordinary substitution, and HGVS says to describe it as one: p.Gly13Val, never
      // p.Gly13delinsVal. Reporting it as a deletion-insertion gave six whitelisted missense hotspots
      // (KRAS G13V and Q61K, TP53 G105V, NOTCH1 A6V, CHEK2 K373E, KMT2D P2717S) the wrong consequence
      // badge and the wrong marker colour on the lollipop.
      consequence = 'missense';
      proteinHgvs = `p.(${three(refMid[0])}${startRes}${three(mutMid[0])})`;
      proteinShort = `${refMid[0]}${startRes}${mutMid[0]}`;
      result.altAa = mutMid[0];
      result.altAaThree = three(mutMid[0]);
    } else if (refMid.length && mutMid.length) {
      consequence = 'inframe-delins';
      proteinHgvs = refMid.length === 1
        ? `p.(${three(refMid[0])}${startRes}delins${three3(mutMid)})`
        : `p.(${three(refMid[0])}${startRes}_${three(refMid[refMid.length - 1])}${endRes}delins${three3(mutMid)})`;
      proteinShort = refMid.length === 1
        ? `${refMid[0]}${startRes}delins${mutMid}`
        : `${refMid[0]}${startRes}_${refMid[refMid.length - 1]}${endRes}delins${mutMid}`;
    } else {
      consequence = 'synonymous';
      proteinHgvs = 'p.(=)';
      proteinShort = '=';
    }
    result.refAa = refProt[f] || null;
    result.refAaThree = result.refAa ? three(result.refAa) : null;
  }

  // context window for the base track
  const from = Math.max(0, start - 1 - flank);
  const refTo = Math.min(seq.length, end + flank);
  const refWin = seq.slice(from, refTo);
  const mutWin = mut.slice(from, Math.min(mut.length, refTo + netNt));

  Object.assign(result, {
    ok: true,
    message: result.refMatches === false
      ? 'Resolved, but the stated reference bases do not match the transcript (check the transcript version).'
      : 'Resolved.',
    netNt, frameshift, firstAffected, consequence, proteinHgvs, proteinShort,
    insertedNt, deletedNt,
    proteinLength: refProt.length, mutantProteinLength: mutProt.length,
    window: { ref: refWin, alt: mutWin, startCds: from + 1, endCds: refTo, eventIndex: start - 1 - from, eventLength: end - start + 1 },
  });
  return result;
}

/** Human-readable label for a consequence class. */
export function consequenceLabel(consequence) {
  return ({
    synonymous: 'Synonymous (silent) — same amino acid',
    missense: 'Missense — amino acid changed',
    nonsense: 'Nonsense — new premature stop',
    'stop-loss': 'Stop-loss — stop codon lost, protein extended',
    'stop-retained': 'Stop-retained — still a stop',
    'start-loss': 'Start-loss — start codon disrupted',
    frameshift: 'Frameshift — reading frame shifted, new stop downstream',
    'inframe-deletion': 'In-frame deletion — residue(s) removed, reading frame kept',
    'inframe-insertion': 'In-frame insertion — residue(s) added, reading frame kept',
    'inframe-duplication': 'In-frame duplication — residue(s) duplicated, reading frame kept',
    'inframe-delins': 'In-frame deletion-insertion — residues replaced, reading frame kept',
    unknown: 'Unknown — could not translate a codon',
  })[consequence] || consequence;
}
