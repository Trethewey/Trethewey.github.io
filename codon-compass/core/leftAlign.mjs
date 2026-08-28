// Left-aligning an insertion or deletion to the position gnomAD, ClinVar's search index and variant
// callers file it under.
//
// One extra CAG in a run of thirteen can be written at any of thirteen places. HGVS writes a
// duplication at the 3′-most one (SMARCA2 c.705_707dup, which Ensembl renders as 9-2039814-A-AGCA);
// the VCF convention puts it at the left-most (9-2039776-A-ACAG). Same allele, 38 bases apart, and a
// lookup keyed on the first came back "not found" — which the app then showed as "not observed in
// gnomAD, consistent with a somatic or very rare variant" for an allele carried by about 1 in 30
// people (rs113070757). CFTR p.Phe508del (7-117559591-TCTT-T against gnomAD's 7-117559590-ATCT-A)
// and ASXL1 c.1934dup (20-32434645-G-GG against 20-32434638-A-AG) fail the same way.
//
// The rule is the one bcftools norm and vt normalize apply (Tan, Abecasis and Kang, 2015):
//   1. while both alleles end in the same base, drop it from both; whenever either allele becomes
//      empty, put the reference base to the left on the front of both and step the position back one;
//   2. then drop shared leading bases while both alleles keep at least one.
// A single-base substitution has one place to sit and is returned untouched, without a window.
//
// Pure: the reference bases arrive as a window of plus-strand sequence, so the controls run offline
// and nothing here reaches the network. The window is read only over the reference allele itself —
// which is checked against the genome before anything moves — and to the left of it.

/** True for any change other than a one-base substitution: only those can sit in more than one place. */
export function needsLeftAlignment(ref, alt) {
  const r = String(ref || '').toUpperCase();
  const a = String(alt || '').toUpperCase();
  return /^[ACGT]+$/.test(r) && /^[ACGT]+$/.test(a) && (r.length > 1 || a.length > 1);
}

/** "9-2039776-A-ACAG" — the form gnomAD keys on. */
export function variantId(c) {
  return `${String(c.chrom ?? c.chr).replace(/^chr/i, '')}-${c.pos}-${c.ref}-${c.alt}`;
}

/**
 * @param {{chrom?:string,chr?:string,pos:number,ref:string,alt:string}} coord  VCF-style: 1-based, plus strand, anchor base included
 * @param {{start:number,seq:string}|null} window  plus-strand reference; `start` is the 1-based position of seq[0]
 * @returns {{ok:true,chrom:string,pos:number,ref:string,alt:string,moved:boolean,from:{pos:number,ref:string,alt:string}}
 *         |{ok:false,code:'window'|'uncovered'|'mismatch'|'unreadable',reason:string,neededBefore?:number}}
 */
export function leftAlign(coord, window) {
  const chrom = String(coord.chrom ?? coord.chr).replace(/^chr/i, '');
  const pos = Number(coord.pos);
  const ref = String(coord.ref || '').toUpperCase();
  const alt = String(coord.alt || '').toUpperCase();
  const from = { pos, ref, alt };
  if (!needsLeftAlignment(ref, alt)) return { ok: true, chrom, pos, ref, alt, moved: false, from };
  if (!window || !Number.isInteger(window.start) || typeof window.seq !== 'string') {
    return { ok: false, code: 'window', reason: 'no reference window', neededBefore: pos };
  }
  const seq = window.seq.toUpperCase();
  const end = window.start + seq.length - 1;
  const baseAt = (p) => (p >= window.start && p <= end ? seq[p - window.start] : null);

  // The reference allele must be what the genome carries here, or the walk below would be shifting
  // a description of the wrong sequence.
  if (pos < window.start || pos + ref.length - 1 > end) {
    return { ok: false, code: 'uncovered', reason: `window ${chrom}:${window.start}-${end} does not cover ${chrom}:${pos}-${pos + ref.length - 1}` };
  }
  const onGenome = seq.slice(pos - window.start, pos - window.start + ref.length);
  if (onGenome !== ref) {
    return { ok: false, code: 'mismatch', reason: `reference mismatch at ${chrom}:${pos}: genome ${onGenome}, allele ${ref}` };
  }

  let p = pos; let r = ref; let a = alt;
  for (;;) {
    if (r.length && a.length && r[r.length - 1] === a[a.length - 1]) { r = r.slice(0, -1); a = a.slice(0, -1); continue; }
    if (!r.length || !a.length) {
      const b = baseAt(p - 1);
      if (b === null) return { ok: false, code: 'window', reason: `window starts at ${chrom}:${window.start}; the walk needs the base at ${chrom}:${p - 1}`, neededBefore: p };
      if (!/^[ACGT]$/.test(b)) return { ok: false, code: 'unreadable', reason: `reference base ${b} at ${chrom}:${p - 1}` };
      r = b + r; a = b + a; p -= 1; continue;
    }
    break;
  }
  while (r.length > 1 && a.length > 1 && r[0] === a[0]) { r = r.slice(1); a = a.slice(1); p += 1; }

  const moved = p !== pos || r !== ref || a !== alt;
  return { ok: true, chrom, pos: p, ref: r, alt: a, moved, from };
}
