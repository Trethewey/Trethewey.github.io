// Turning an HGVS coding coordinate into a genomic one — including the ones that are not in the
// coding sequence at all.
//
// This is what let intronic and untranslated-region variants dead-end: `NM_001614.5:c.803-6C>T`
// parsed correctly and then stopped, because nothing converted `c.803-6` into a position on the
// chromosome. Everything downstream — gnomAD, ClinVar, SpliceAI — is keyed on that position, and
// SpliceAI in particular exists FOR these variants, so the one missing step blocked the tool's whole
// answer for the class of variant it matters most for.
//
// The four forms, all relative to the coding sequence of one transcript:
//   c.803       inside the coding sequence
//   c.803+6     6 bases into the intron AFTER coding position 803   (donor side)
//   c.803-6     6 bases into the intron BEFORE coding position 803  (acceptor side)
//   c.-28       28 bases before the A of the start codon (5' untranslated region)
//   c.*15       15 bases after the last base of the stop codon (3' untranslated region)
//
// Everything here works in TRANSCRIPT order, then converts to genomic at the last step. For a
// minus-strand gene transcript order runs downwards through the chromosome, so "the next base along
// the transcript" is genomic position MINUS one. Getting that backwards puts the variant in the wrong
// intron — on the other side of the exon — and every downstream answer inherits the error silently,
// which is why the strand conversion is done once, here, and controlled both ways.
//
// Coordinates in data/gene-exons.json are genomic, 1-based, inclusive.

/**
 * The alleles as the CHROMOSOME carries them.
 *
 * An HGVS `c.` change is written in transcript orientation. For a minus-strand gene that is the
 * reverse complement of the genome, so KRAS `c.35G>A` is `C>T` at chr12:25,245,350 — which is what
 * every coordinate-keyed service expects. Sending the transcript-orientation alleles instead asks
 * about an allele that does not exist at that position, and the answer comes back empty rather than
 * wrong, which is why it went unnoticed the last time it happened here: `hgvsToGenomic` silently
 * broke gnomAD, SpliceAI and ClinVar together.
 */
const COMPLEMENT = { A: 'T', T: 'A', C: 'G', G: 'C', N: 'N' };
export function toGenomicAlleles(strand, ref, alt) {
  const flip = (s) => String(s || '').toUpperCase().split('').reverse().map((c) => COMPLEMENT[c] || c).join('');
  return strand === '-'
    ? { ref: flip(ref), alt: flip(alt), complemented: true }
    : { ref: String(ref || '').toUpperCase(), alt: String(alt || '').toUpperCase(), complemented: false };
}

/** Pairs of a flat [start, end, start, end, …] array, ascending genomically. */
function pairs(flat) {
  const out = [];
  for (let i = 0; i + 1 < (flat || []).length; i += 2) out.push([flat[i], flat[i + 1]]);
  return out;
}

/** Exon blocks in transcript 5′→3′ order. For a minus-strand gene that is genomic order reversed. */
function blocksInTranscriptOrder(flat, strand) {
  const b = pairs(flat);
  return strand === '-' ? b.reverse() : b;
}

/**
 * The genomic position of coding base `n` (c.n, 1-based), or null when the transcript has no coding
 * sequence or n lies outside it.
 */
export function codingBaseToGenomic(gene, n) {
  const blocks = blocksInTranscriptOrder(gene.cdsExons, gene.strand);
  let seen = 0;
  for (const [lo, hi] of blocks) {
    const len = hi - lo + 1;
    if (seen + len >= n) {
      const offset = n - seen - 1;                 // 0-based within this block, in transcript order
      return gene.strand === '-' ? hi - offset : lo + offset;
    }
    seen += len;
  }
  return null;
}

/** Total coding bases of the transcript. */
export function codingLength(gene) {
  return pairs(gene.cdsExons).reduce((a, [lo, hi]) => a + hi - lo + 1, 0);
}

/** One step along the transcript, in genomic terms: +1 on the plus strand, −1 on the minus. */
const step = (gene) => (gene.strand === '-' ? -1 : 1);

/**
 * Which exon (1-based, transcript order) contains a genomic position, and where it sits in it.
 * Returns null when the position is not inside an exon.
 */
export function exonAt(gene, pos) {
  const ex = blocksInTranscriptOrder(gene.exons, gene.strand);
  for (let i = 0; i < ex.length; i += 1) {
    const [lo, hi] = ex[i];
    if (pos >= lo && pos <= hi) {
      // distance to each end of the exon, in transcript order
      const fromStart = gene.strand === '-' ? hi - pos : pos - lo;
      const toEnd = gene.strand === '-' ? pos - lo : hi - pos;
      return { number: i + 1, of: ex.length, fromStart, toEnd };
    }
  }
  return null;
}

/**
 * Resolve a parsed HGVS coding coordinate to a genomic position.
 *
 * `parsed` is what src/core/hgvs.mjs produces: { region, position, offset } where `region` is
 * 'cds' | 'intronic' | 'utr5' | 'utr3'. An intronic coordinate carries the anchor coding position and
 * a signed offset; c.-N and c.*N carry the distance beyond each end of the coding sequence.
 *
 * Returns { ok, pos, chr, strand, region, anchorPos, offset, exon, intron, note } or
 * { ok: false, reason } — never a guessed position.
 */
export function resolveCodingCoordinate(gene, parsed) {
  const out = resolvePosition(gene, parsed);
  if (!out.ok) return out;
  // Carry the alleles in CHROMOSOME orientation alongside the position, so a caller cannot pair a
  // genomic coordinate with transcript-orientation alleles — which is the shape of the bug that broke
  // gnomAD, SpliceAI and ClinVar here once before.
  if (parsed && parsed.ref && parsed.alt) {
    const g = toGenomicAlleles(gene.strand, parsed.ref, parsed.alt);
    out.ref = g.ref; out.alt = g.alt; out.complemented = g.complemented;
    out.id = `chr${out.chr}-${out.pos}-${g.ref}-${g.alt}`;
  }
  return out;
}

function resolvePosition(gene, parsed) {
  if (!gene || !gene.chr) return { ok: false, reason: 'no exon record for this gene' };
  if (!(gene.cdsExons || []).length) return { ok: false, reason: 'the transcript has no coding sequence in the bundled table' };
  const chr = gene.chr; const strand = gene.strand;
  const cdsLen = codingLength(gene);
  const d = step(gene);

  const region = parsed.region || 'cds';
  const offset = Number(parsed.offset || 0);
  const at = Number(parsed.position);

  // ---- inside the coding sequence -------------------------------------------------------------
  if (region === 'cds' && !offset) {
    if (!(at >= 1 && at <= cdsLen)) return { ok: false, reason: `c.${at} is outside the coding sequence (${cdsLen} bases)` };
    const pos = codingBaseToGenomic(gene, at);
    return { ok: true, pos, chr, strand, region: 'cds', anchorPos: pos, offset: 0, exon: exonAt(gene, pos) };
  }

  // ---- 5′ untranslated region: c.-N, N bases before the first coding base ----------------------
  if (region === 'utr5') {
    const first = codingBaseToGenomic(gene, 1);
    const n = Math.abs(at);
    const pos = first - d * n;
    return withinTranscript(gene, { ok: true, pos, chr, strand, region: 'utr5', anchorPos: first, offset: -n,
      note: `${n} base${n === 1 ? '' : 's'} before the start codon` });
  }

  // ---- 3′ untranslated region: c.*N, N bases after the last coding base ------------------------
  if (region === 'utr3') {
    const last = codingBaseToGenomic(gene, cdsLen);
    const n = Math.abs(at);
    const pos = last + d * n;
    return withinTranscript(gene, { ok: true, pos, chr, strand, region: 'utr3', anchorPos: last, offset: n,
      note: `${n} base${n === 1 ? '' : 's'} after the stop codon` });
  }

  // ---- intronic: c.N+M (donor side) or c.N-M (acceptor side) -----------------------------------
  if (region === 'intronic' || offset) {
    if (!(at >= 1 && at <= cdsLen)) return { ok: false, reason: `c.${at} is outside the coding sequence (${cdsLen} bases)` };
    if (!offset) return { ok: false, reason: 'an intronic coordinate needs an offset' };
    const anchor = codingBaseToGenomic(gene, at);
    const anchorExon = exonAt(gene, anchor);
    // An intronic coordinate is only meaningful when its anchor is ON an exon boundary: c.N+M counts
    // out of the exon that ENDS at c.N, and c.N−M counts back out of the exon that BEGINS at c.N.
    // Without this check the arithmetic still produces a position — it just walks M bases along the
    // exon and reports a coding base as though it were intronic. BRAF c.1799 is the V600E base in the
    // middle of exon 15, so "c.1799+1" is not a coordinate on this transcript at all, and saying
    // chr7:140,753,335 for it would be a fabricated answer that looks entirely plausible.
    if (!anchorExon) return { ok: false, reason: `c.${at} does not fall in an exon of this transcript` };
    const atBoundary = offset > 0 ? anchorExon.toEnd === 0 : anchorExon.fromStart === 0;
    if (!atBoundary) {
      const side = offset > 0 ? 'end' : 'start';
      return { ok: false,
        reason: `c.${at} is not at the ${side} of exon ${anchorExon.number}, so c.${at}${offset > 0 ? '+' : ''}${offset} is not a position on this transcript. `
          + `An intronic coordinate counts outwards from an exon boundary.` };
    }
    // The offset counts along the TRANSCRIPT from the anchor base, so it steps by the strand.
    const pos = anchor + d * offset;
    // Which intron: the one after the anchor's exon for +N, before it for −N.
    const intron = anchorExon
      ? (offset > 0
        ? { after: anchorExon.number, before: anchorExon.number + 1 }
        : { after: anchorExon.number - 1, before: anchorExon.number })
      : null;
    return withinTranscript(gene, {
      ok: true, pos, chr, strand, region: 'intronic', anchorPos: anchor, offset, exon: anchorExon, intron,
    });
  }

  return { ok: false, reason: `cannot place a ${region} coordinate` };
}

/** Guard: a resolved position that falls outside the transcript's own span is not reported as fact. */
function withinTranscript(gene, res) {
  const all = pairs(gene.exons);
  const lo = Math.min(...all.map((p) => p[0]));
  const hi = Math.max(...all.map((p) => p[1]));
  if (res.pos < lo || res.pos > hi) {
    return { ...res, ok: true, beyondTranscript: true,
      note: `${res.note ? `${res.note}. ` : ''}This lands outside the transcript's annotated span (${lo.toLocaleString()}–${hi.toLocaleString()}), so it is past the end of the record rather than inside it.` };
  }
  return res;
}

// ---------------------------------------------------------------------------------------------
//  Splice-region classification
// ---------------------------------------------------------------------------------------------
// The convention is stated rather than assumed, because published definitions differ at the edges:
//   * the ESSENTIAL splice site is the first two intronic bases at each end of an intron — the GT of
//     the donor (+1, +2) and the AG of the acceptor (−1, −2). Almost all variants here abolish the site.
//   * the SPLICE REGION, as the Sequence Ontology defines it, is 1–3 bases INTO the exon and 3–8 bases
//     into the intron at either end.
// Anything further into an intron than that is reported as deep intronic, with its distance, and no
// consequence is claimed for it — that is exactly what SpliceAI is for.
export const SPLICE_CONVENTION =
  'Essential site: the first two intronic bases at either end of an intron (the donor GT and acceptor AG). '
  + 'Splice region: 1–3 bases into the exon and 3–8 bases into the intron, as the Sequence Ontology defines it. '
  + 'Further in than that is reported as deep intronic, with the distance, and no effect is claimed.';

/**
 * What a resolved coordinate means for splicing. Returns { kind, distance, side, exon, label } where
 * kind is 'essential-site' | 'splice-region' | 'deep-intronic' | 'exonic-splice-region' | 'none'.
 */
export function spliceRegion(gene, res) {
  if (!res || !res.ok) return { kind: 'none', label: 'not placed' };

  if (res.region === 'intronic') {
    const n = Math.abs(res.offset);
    const side = res.offset > 0 ? 'donor' : 'acceptor';
    const where = res.intron ? ` of the intron between exons ${res.intron.after} and ${res.intron.before}` : '';
    if (n <= 2) {
      return { kind: 'essential-site', distance: n, side,
        label: `Essential ${side} splice site — ${side === 'donor' ? '+' : '−'}${n}${where}. These two bases are the ${side === 'donor' ? 'GT' : 'AG'} the spliceosome binds.` };
    }
    if (n <= 8) {
      return { kind: 'splice-region', distance: n, side,
        label: `Splice region — ${side === 'donor' ? '+' : '−'}${n}${where}, within the 3–8 base window either side of the intron.` };
    }
    return { kind: 'deep-intronic', distance: n, side,
      label: `Deep intronic — ${n} bases into the intron on the ${side} side${where}. Too far in for a rule to say anything; a splice predictor is the tool for this.` };
  }

  if (res.region === 'cds' && res.exon) {
    // 1–3 bases from either end of an exon is inside the splice region too
    const toDonor = res.exon.toEnd;      // bases to the exon's 3′ end, transcript order
    const toAcceptor = res.exon.fromStart;
    if (toDonor <= 2) {
      return { kind: 'exonic-splice-region', distance: toDonor + 1, side: 'donor',
        label: `Last ${toDonor + 1} base${toDonor ? 's' : ''} of exon ${res.exon.number} — inside the splice region, so this may affect splicing as well as the protein.` };
    }
    if (toAcceptor <= 2) {
      return { kind: 'exonic-splice-region', distance: toAcceptor + 1, side: 'acceptor',
        label: `First ${toAcceptor + 1} base${toAcceptor ? 's' : ''} of exon ${res.exon.number} — inside the splice region, so this may affect splicing as well as the protein.` };
    }
  }

  if (res.region === 'utr5') return { kind: 'none', label: "5′ untranslated region — before the start codon, so no amino acid changes." };
  if (res.region === 'utr3') return { kind: 'none', label: "3′ untranslated region — after the stop codon, so no amino acid changes." };
  return { kind: 'none', label: 'Within an exon, away from either splice site.' };
}
