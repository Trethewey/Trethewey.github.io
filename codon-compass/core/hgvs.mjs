// Parser for HGVS coding (c.) variant strings.
//
// Fully supported: single-nucleotide substitutions inside the coding region,
//   e.g. NM_006015.6:c.461A>C
// Recognised but not translated (clearly flagged): UTR positions (c.-12, c.*34),
//   intronic/splice offsets (c.461+3), and del/dup/ins/delins events. The parser
//   still returns what it understood so the UI can explain why it stopped.

const COORD_RE = /^([cngmr])\.(.+)$/i;

/** Classify a position token into a transcript region. */
function classifyRegion(token) {
  if (/^-\d+$/.test(token)) return 'utr5';            // 5' UTR, e.g. -12
  if (/^\*\d+$/.test(token)) return 'utr3';           // 3' UTR, e.g. *34
  if (/[+-]\d+$/.test(token) && /^\d/.test(token)) return 'intronic'; // e.g. 461+3
  if (/^\d+$/.test(token)) return 'cds';
  return 'unknown';
}

/**
 * Split a coordinate token into the coding position it is anchored to and its signed offset, so a
 * non-coding coordinate can be placed on the chromosome rather than only recognised.
 *
 *   '803'    -> { position: 803,  offset: 0 }
 *   '803+6'  -> { position: 803,  offset: 6 }    6 bases into the intron after 803 (donor side)
 *   '803-6'  -> { position: 803,  offset: -6 }   6 bases into the intron before 803 (acceptor side)
 *   '-28'    -> { position: 28,   offset: 0 }    28 bases before the start codon
 *   '*15'    -> { position: 15,   offset: 0 }    15 bases after the stop codon
 *   '*15+3'  -> { position: 15,   offset: 3 }
 *
 * `position` is always the magnitude; which end of the transcript it counts from is carried by the
 * region, not by a sign, so a caller cannot silently read a 5′ offset as a coding position.
 */
export function splitCoordinateToken(token) {
  const m = /^(\*?)(-?\d+)([+-]\d+)?$/.exec(String(token || '').trim());
  if (!m) return { position: null, offset: 0 };
  return { position: Math.abs(parseInt(m[2], 10)), offset: m[3] ? parseInt(m[3], 10) : 0 };
}

/**
 * Parse an HGVS c. string.
 * @returns {{
 *   raw:string, ok:boolean, transcript:(string|null), coordinate:(string|null),
 *   kind:('substitution'|'deletion'|'duplication'|'insertion'|'delins'|'unsupported'),
 *   region:(string|null), position:(number|null), ref:(string|null), alt:(string|null),
 *   supported:boolean, message:string
 * }}
 */
export function parseHgvs(input) {
  const raw = String(input || '').trim();
  const base = {
    raw, ok: false, transcript: null, coordinate: null, kind: 'unsupported',
    region: null, position: null, ref: null, alt: null, supported: false, message: '',
  };
  if (!raw) return { ...base, message: 'Enter a variant, e.g. NM_006015.6:c.461A>C' };

  let transcript = null;
  let rest = raw;
  if (raw.includes(':')) {
    const idx = raw.indexOf(':');
    transcript = raw.slice(0, idx).trim();
    rest = raw.slice(idx + 1).trim();
  }

  const cm = rest.match(COORD_RE);
  if (!cm) {
    return { ...base, transcript, message: 'Could not find a coordinate prefix such as "c.". Expected e.g. c.461A>C' };
  }
  const coordinate = cm[1].toLowerCase();
  const body = cm[2].trim();
  const out = { ...base, transcript, coordinate };

  if (coordinate !== 'c') {
    return { ...out, message: `Only coding "c." variants are assessed here; got "${coordinate}."` };
  }

  // Single-nucleotide substitution inside the CDS: digits, one base, ">", one base.
  const sub = body.match(/^(\d+)([ACGT])>([ACGT])$/i);
  if (sub) {
    return {
      ...out, ok: true, kind: 'substitution', region: 'cds',
      position: parseInt(sub[1], 10), ref: sub[2].toUpperCase(), alt: sub[3].toUpperCase(),
      supported: true, message: 'Coding single-nucleotide substitution.',
    };
  }

  // Substitution outside the plain CDS (UTR or intronic) — recognised, not translated.
  const subAny = body.match(/^(\*?-?\d+(?:[+-]\d+)?)([ACGT]+)>([ACGT]+)$/i);
  if (subAny) {
    const region = classifyRegion(subAny[1]);
    const where = region === 'utr5' ? "the 5' untranslated region"
      : region === 'utr3' ? "the 3' untranslated region"
      : region === 'intronic' ? 'an intron / splice region'
      : 'a non-coding position';
    // `supported` stays false: the protein effect genuinely cannot be computed from the coding
    // sequence. But the coordinate itself IS placeable, and `position`/`offset` are what
    // src/core/codingCoordinate.mjs needs to put it on the chromosome — which is what then lets
    // gnomAD, ClinVar and SpliceAI answer for it. Without these fields the variant parsed and
    // then stopped, which is the whole reason intronic variants dead-ended.
    const { position, offset } = splitCoordinateToken(subAny[1]);
    return {
      ...out, ok: true, kind: 'substitution', region,
      position, offset, placeable: position != null && region !== 'unknown',
      ref: subAny[2].toUpperCase(), alt: subAny[3].toUpperCase(),
      supported: false,
      message: `This substitution is in ${where}; its protein effect is not computed from the coding sequence alone.`,
    };
  }

  // Coding-region indels: delins must be tested before del. Plain digit positions
  // only — anything with +/-/* offsets is intronic/UTR and cannot be resolved here.
  let m;
  if ((m = body.match(/^(\d+)(?:_(\d+))?delins([ACGT]+)$/i))) {
    return { ...out, ok: true, kind: 'delins', region: 'cds', start: parseInt(m[1], 10), end: m[2] ? parseInt(m[2], 10) : parseInt(m[1], 10), insSeq: m[3].toUpperCase(), supported: true, message: 'Coding deletion-insertion.' };
  }
  if ((m = body.match(/^(\d+)(?:_(\d+))?del([ACGT]*)$/i))) {
    return { ...out, ok: true, kind: 'deletion', region: 'cds', start: parseInt(m[1], 10), end: m[2] ? parseInt(m[2], 10) : parseInt(m[1], 10), delSeq: (m[3] || '').toUpperCase(), supported: true, message: 'Coding deletion.' };
  }
  if ((m = body.match(/^(\d+)(?:_(\d+))?dup([ACGT]*)$/i))) {
    return { ...out, ok: true, kind: 'duplication', region: 'cds', start: parseInt(m[1], 10), end: m[2] ? parseInt(m[2], 10) : parseInt(m[1], 10), dupSeq: (m[3] || '').toUpperCase(), supported: true, message: 'Coding duplication.' };
  }
  if ((m = body.match(/^(\d+)_(\d+)ins([ACGT]+)$/i))) {
    return { ...out, ok: true, kind: 'insertion', region: 'cds', start: parseInt(m[1], 10), end: parseInt(m[2], 10), insSeq: m[3].toUpperCase(), supported: true, message: 'Coding insertion.' };
  }

  // Same events outside the plain coding numbering (intronic / UTR) — recognised, not resolvable.
  let kind = 'unsupported';
  if (/delins/i.test(body)) kind = 'delins';
  else if (/del/i.test(body)) kind = 'deletion';
  else if (/dup/i.test(body)) kind = 'duplication';
  else if (/ins/i.test(body)) kind = 'insertion';
  if (kind !== 'unsupported') {
    return {
      ...out, ok: true, kind, region: classifyRegion(body.match(/^[*\-\d+]+/)?.[0] || ''),
      supported: false,
      message: `Recognised a ${kind} event outside the plain coding numbering (intronic or untranslated). Its protein effect cannot be computed from the coding sequence alone.`,
    };
  }

  return { ...out, message: `Could not parse the change "${body}". Supported example: c.461A>C` };
}
