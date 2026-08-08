// Exon-level reasoning about a breakpoint: which exon or intron it falls in, how much coding sequence
// survives on the kept side, and whether two pieces join in frame.
//
// Everything here works in TRANSCRIPT order. For a minus-strand gene that is the reverse of genomic
// order, and getting that backwards is the single easiest way to produce a confident wrong answer, so
// every function converts to transcript order first and never mixes the two.
//
// Coordinates come from data/gene-exons.json (MANE Select, genomic, 1-based inclusive):
//   { tx, chr, strand, exons: [s,e,s,e,…], cds: [start,end]|null, cdsExons: [s,e,s,e,…] }

/** Pairs of a flat [s,e,s,e,…] array, in ascending genomic order. */
function pairs(flat) {
  const out = [];
  for (let i = 0; i + 1 < flat.length; i += 2) out.push([flat[i], flat[i + 1]]);
  return out;
}

/** The gene's exons in transcript 5′→3′ order, numbered from 1. */
export function exonsInTranscriptOrder(g) {
  const list = pairs(g.exons || []).map(([start, end], i) => ({ start, end, genomicIndex: i }));
  if (g.strand === '-') list.reverse();
  return list.map((e, i) => ({ ...e, number: i + 1 }));
}

/** Total exon count. */
export function exonCount(g) { return Math.floor((g.exons || []).length / 2); }

/**
 * Where a genomic position falls in the gene's structure.
 * Returns { kind: 'exon'|'intron'|'outside', exon, afterExon, beforeExon } with exon numbers in
 * transcript order. An intronic position reports the exons on either side of it.
 */
export function featureAt(g, pos) {
  const ex = exonsInTranscriptOrder(g);
  if (!ex.length) return { kind: 'outside' };
  for (const e of ex) if (pos >= e.start && pos <= e.end) return { kind: 'exon', exon: e.number, of: ex.length };
  // between two exons? compare in genomic space, then report in transcript numbering
  const genomic = pairs(g.exons);
  for (let i = 0; i + 1 < genomic.length; i += 1) {
    if (pos > genomic[i][1] && pos < genomic[i + 1][0]) {
      const a = ex.find((e) => e.genomicIndex === i).number;
      const b = ex.find((e) => e.genomicIndex === i + 1).number;
      return { kind: 'intron', afterExon: Math.min(a, b), beforeExon: Math.max(a, b), of: ex.length };
    }
  }
  return { kind: 'outside', of: ex.length };
}

/**
 * Coding bases of this gene that survive on the kept side of a breakpoint.
 *
 * `keepLeft` is genomic: true keeps everything at or below `pos`. Whether that is the gene's 5′ or 3′
 * portion depends on its strand, which is exactly what this resolves.
 *
 * A break inside an intron is resolved by splicing: the kept side keeps whole exons only, because a
 * transcript spliced across the junction cannot include half an intron.
 *
 * Returns { bases, side: '5prime'|'3prime', lastWholeExon, splitExon } where `bases` counts coding
 * bases in transcript order from the CDS start (for a 5′ piece) or to the CDS end (for a 3′ piece).
 */
export function codingRetained(g, pos, keepLeft) {
  const cds = pairs(g.cdsExons || []);
  if (!cds.length) return { bases: 0, side: null, coding: false };
  const forward = g.strand !== '-';
  // the kept side in transcript terms: a plus-strand gene keeping the low side keeps its 5′ portion
  const side = (forward === keepLeft) ? '5prime' : '3prime';

  let bases = 0;
  for (const [s, e] of cds) {
    if (keepLeft) {
      if (e <= pos) bases += e - s + 1;
      else if (s <= pos) bases += pos - s + 1;      // the breakpoint splits this coding block
    } else {
      if (s >= pos) bases += e - s + 1;
      else if (e >= pos) bases += e - pos + 1;
    }
  }
  return { bases, side, coding: bases > 0 };
}

/**
 * Coding offset, in transcript order, of the first coding base retained on a 3′ piece — i.e. how many
 * coding bases of this gene were LOST at its 5′ end. This is what decides the acceptor's own frame.
 */
export function codingOffsetOfAcceptor(g, pos, keepLeft) {
  const total = pairs(g.cdsExons || []).reduce((n, [s, e]) => n + e - s + 1, 0);
  const kept = codingRetained(g, pos, keepLeft).bases;
  return total - kept;
}

/** Total coding bases of a gene, from its coding exons. */
export function codingTotal(g) {
  return pairs(g.cdsExons || []).reduce((n, [s, e]) => n + e - s + 1, 0);
}

/**
 * Whether a chimeric junction joins the two coding sequences in frame.
 *
 * The fused coding sequence is the donor's retained coding bases followed by the acceptor's. The
 * acceptor's surviving sequence only reads in its own native frame if the number of donor bases in
 * front of it is congruent, modulo three, to the number of acceptor bases that were removed.
 *
 * `shift` and `inFrame` are that arithmetic and nothing more. They say the two coding sequences line
 * up; they do not say a ribosome ever gets from one to the other.
 *
 * ONE CASE HAS NO FUSED READING FRAME AT ALL. When the donor keeps EVERY coding base on its 5′ side
 * the breakpoint lies 3′ of its stop codon, so that stop travels into the fused message ahead of the
 * junction and translation ends there — no fusion protein forms however neatly the two frames line
 * up. `donorStopRetained` says so, and `inFrame`, `shift` and `fusedCodingBases` are then all null:
 * there is nothing to be in or out of, and nothing fused to measure. They are null rather than
 * absent so that a caller reading them gets no answer instead of a wrong one — the mod-three sum is
 * zero often enough that reporting it is how this used to promise a protein that cannot be made.
 *
 * Returns null when either side has no coding sequence retained — there is nothing to be in frame.
 * That null is an ANSWER and never means "no data": this function is only ever reached with both
 * exon tables in hand. A caller that reports it as missing exon coordinates states the opposite of
 * what was computed, which is what the fusion classifier used to do.
 */
export function junctionFrame(donor, donorPos, donorKeepLeft, acceptor, accPos, accKeepLeft) {
  const d = codingRetained(donor, donorPos, donorKeepLeft);
  const a = codingRetained(acceptor, accPos, accKeepLeft);
  if (!d.coding || !a.coding) return null;
  const lost = codingOffsetOfAcceptor(acceptor, accPos, accKeepLeft);
  const shift = ((d.bases - lost) % 3 + 3) % 3;
  const donorTotal = codingTotal(donor);
  const donorStopRetained = donorTotal > 0 && d.side === '5prime' && d.bases === donorTotal;
  return {
    inFrame: donorStopRetained ? null : shift === 0,
    shift: donorStopRetained ? null : shift,          // 0, 1 or 2 bases out; null = no fused frame
    donorCodingBases: d.bases,
    donorCodingTotal: donorTotal,
    donorCodons: d.bases / 3,
    donorStopRetained,                  // true = translation ends before the junction
    acceptorCodingBases: a.bases,
    acceptorCodingLost: lost,
    // The donor's own coding sequence, stop codon and all, is not part of a fusion protein, so it is
    // not summed with the acceptor's as if it were.
    fusedCodingBases: donorStopRetained ? null : d.bases + a.bases,
  };
}

/**
 * The exons this piece contributes to the fusion transcript, in transcript 5′→3′ order.
 *
 * This is the fusion PRODUCT, not the gene: exons the breakpoint removes are absent, and an exon the
 * breakpoint cuts through is truncated to the surviving part. Introns are gone, because they are
 * spliced out — which is why a fusion diagram is drawn in transcript space, not genomic space.
 *
 * Each entry: { number, start, end, length, coding, partial }
 *   number  — the exon's number in its own gene, in transcript order
 *   coding  — 'full' | 'none' | 'partial' relative to the coding sequence
 *   partial — true when the breakpoint cut this exon short
 */
export function retainedExons(g, pos, keepLeft) {
  const all = exonsInTranscriptOrder(g);
  const cds = g.cds || null;
  const out = [];
  for (const e of all) {
    // does any of this exon survive?
    const survives = keepLeft ? e.start <= pos : e.end >= pos;
    if (!survives) continue;
    // clip it to the surviving side
    const start = keepLeft ? e.start : Math.max(e.start, pos);
    const end = keepLeft ? Math.min(e.end, pos) : e.end;
    if (end < start) continue;
    const codingBases = cds ? Math.max(0, Math.min(end, cds[1]) - Math.max(start, cds[0]) + 1) : 0;
    const length = end - start + 1;
    out.push({
      number: e.number,
      start,
      end,
      length,
      coding: !cds || codingBases === 0 ? 'none' : (codingBases === length ? 'full' : 'partial'),
      codingBases,
      partial: start !== e.start || end !== e.end,
    });
  }
  return out;
}

/** A plain-language description of where a breakpoint sits, for display. */
export function featureLabel(g, pos) {
  const f = featureAt(g, pos);
  if (f.kind === 'exon') return `in exon ${f.exon} of ${f.of}`;
  if (f.kind === 'intron') return `in the intron between exons ${f.afterExon} and ${f.beforeExon} (of ${f.of})`;
  return 'outside the transcribed region';
}
