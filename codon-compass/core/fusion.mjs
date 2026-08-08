// Breakend (BND) interpretation for translocation / gene-fusion assessment.
//
// Pure logic — no DOM, no network — so it is imported by both the renderer and `node --test`.
// Scope: BND records only (translocations/fusions). CNV and intra-chromosomal DEL/DUP/INV rows from the
// SV report are deliberately out of scope for now.
//
// The SV report gives each breakend as a local locus glued to a VCF-style ALT string, e.g.
//   "3:187746165T]12:25056225]"            local 3:187746165, ALT "T]12:25056225]"
//   "14:105864633ACCC[18:63126216["         inserted bases ACCC before the bracket (t-first)
//   "14:105891702]18:63126209]GGG(+4)CTC"   inserted sequence after the brackets (t-last)
//
// VCF 4.x breakend ALT grammar (t = local ref/inserted bases, p = mate locus). The spec words the
// four forms as "which mate piece is joined, and on which side of t":
//   t[p[   local LEFT (<=pos) first,  then mate RIGHT (>=p)          -- mate as-is
//   t]p]   local LEFT (<=pos) first,  then mate LEFT  (<=p) reversed
//   ]p]t   mate LEFT  (<=p)   first,  then local RIGHT (>=pos)       -- mate as-is
//   [p[t   mate RIGHT (>=p)   first reversed, then local RIGHT (>=pos)
// So: t before the brackets => local LEFT kept and written first;  t after => local RIGHT kept and
//     written second.  Bracket "[" => mate RIGHT kept;  "]" => mate LEFT kept.
//
// The RELATIVE ORIENTATION of the two pieces is NOT the bracket on its own. A junction joins the
// end of one piece to the start of the next, so two pieces that are both "LEFT of their breakpoint"
// (or both "RIGHT of it") can only meet if one of them is read backwards:
//     local LEFT + mate RIGHT  -> same orientation      local LEFT + mate LEFT  -> reverse
//     local RIGHT + mate LEFT  -> same orientation      local RIGHT + mate RIGHT -> reverse
// i.e. reverse-complement exactly when the two kept sides are the same word. Reading it off the
// bracket alone is right for t[p[ / t]p] and inverted for ]p]t / [p[t.

import {
  featureAt, featureLabel, codingRetained, codingOffsetOfAcceptor, exonCount, junctionFrame,
} from './exons.mjs';

/** Normalise a chromosome token: drop a leading "chr", upper-case letters. */
export function normChr(c) {
  return String(c || '').replace(/^chr/i, '').toUpperCase();
}

const LEAD_LOCUS = /^([\w.]+):(\d+)/; // leading local chr:pos
// t = the reference base(s) at the record's own position plus any inserted sequence. Some callers
// elide a run of bases as "(+N)" in the middle of it.
const T_FIELD = /^[ACGTN]+(?:\(\+\d+\)[ACGTN]+)*$/i;

/** Parse one BND cell into a structured breakend, or null if it is not a parseable BND.
 *  Every part of the cell must validate: a cell that only half-parses is rejected, never
 *  silently repaired, because a half-parsed breakpoint is a wrong breakpoint. */
export function parseBnd(cell) {
  const raw = String(cell || '').trim();
  // Coordinates are sometimes written with thousands separators. Strip them deliberately, or the
  // position regex stops at the first comma and 187,746,165 is read as position 187.
  const s = raw.replace(/,(?=\d{3}(?!\d))/g, '');
  const lead = LEAD_LOCUS.exec(s);
  if (!lead) return null;
  const alt = s.slice(lead[0].length);
  const a = parseAlt(alt);
  if (!a) return null;
  return {
    raw,
    localChr: normChr(lead[1]),
    localPos: Number(lead[2]),
    mateChr: a.mateChr,
    matePos: a.matePos,
    bases: a.bases,             // local ref + any inserted sequence, verbatim (e.g. "T", "ACCC", "GGG(+4)CTC")
    bracket: a.bracket,         // "[" or "]"
    tPosition: a.tPosition,     // "first" (t before brackets) or "last"
    localSide: a.localSide,     // "left" = pter..breakpoint (<=pos), "right" = breakpoint..qter (>=pos)
    mateSide: a.mateSide,       // "left" (<=matePos) or "right" (>=matePos)
    strand: a.strand,           // "same" or "reverse" (reverse-complement join)
    order: a.order,             // "local-first" or "mate-first" — the order the ALT form writes the pieces
  };
}

/** Parse just the ALT portion (after the local chr:pos). */
function parseAlt(altRaw) {
  const s = String(altRaw || '').trim();
  if (!s) return null;
  // Both brackets must match, and t must be sequence. Accepting an empty t would silently reinterpret
  // "3:100[5:200[" as the t-last form, flipping which side is kept AND the join orientation; letting
  // t be ".*" would swallow the rest of a space-separated report row and print it as inserted bases.
  const ok = (b1, b2, t) => b1 === b2 && (b1 === '[' || b1 === ']') && T_FIELD.test(t);
  let m;
  if (s[0] === '[' || s[0] === ']') {
    // t-last:  [p[t  or  ]p]t
    m = /^([[\]])([\w.]+):(\d+)([[\]])(.*)$/.exec(s);
    if (!m || !ok(m[1], m[4], m[5])) return null;
    return sides(m[1], normChr(m[2]), Number(m[3]), m[5], 'last');
  }
  // t-first:  t[p[  or  t]p]
  m = /^(.*?)([[\]])([\w.]+):(\d+)([[\]])$/.exec(s);
  if (!m || !ok(m[2], m[5], m[1])) return null;
  return sides(m[2], normChr(m[3]), Number(m[4]), m[1], 'first');
}

function sides(bracket, mateChr, matePos, bases, tPosition) {
  const localSide = tPosition === 'first' ? 'left' : 'right';
  const mateSide = bracket === '[' ? 'right' : 'left';
  return {
    mateChr, matePos, bases, bracket, tPosition, localSide, mateSide,
    // Same kept side on both pieces => they can only be joined by reading one backwards.
    strand: localSide === mateSide ? 'reverse' : 'same',
    order: tPosition === 'first' ? 'local-first' : 'mate-first',
  };
}

/** Parse a support token like "PR-40/194;SR-23/159" (either part optional). The tag must stand on
 *  its own — without a boundary, "SPR-1/2" reads as a PR count. When a row carries more than one
 *  PR or SR figure (a normal column beside a tumour one) the first is taken and `ambiguous` is set,
 *  because picking silently would let a normal sample's zero counts stand for the tumour's. */
export function parseSupport(text) {
  const s = String(text || '');
  const grab = (tag) => {
    const all = [...s.matchAll(new RegExp(`(?:^|[^A-Za-z])${tag}-(\\d+)/(\\d+)`, 'g'))];
    return { value: all.length ? { alt: Number(all[0][1]), total: Number(all[0][2]) } : null, count: all.length };
  };
  const pr = grab('PR'); const sr = grab('SR');
  if (!pr.value && !sr.value) return null;
  return { pr: pr.value, sr: sr.value, ambiguous: pr.count > 1 || sr.count > 1 };
}

/** Parse the germline population allele-frequency field (e.g. "0.7116 (GE)", "0.0001"); null for
 *  "N/A"/"-"/absent. A high value flags a common germline SV rather than a somatic driver, so a
 *  wrong read here dismisses a real driver.
 *
 *  Two guards: a frequency cannot exceed 1, and a field tagged "(GE)" is preferred over a bare
 *  decimal. A bare decimal that is the only candidate is still used — that is how the example
 *  reports write it — but is marked `assumed`, because a somatic variant-allele-fraction column
 *  looks identical and the caller must say so rather than present it as fact. */
export function parseGermlineAf(text) {
  const found = [];
  for (const f of String(text || '').split(/\t|\s{2,}/)) {
    // Reports write small frequencies either way: "0.0001" and "1E-4" both appear. A bare integer is
    // deliberately not accepted — that shape is a count, not a frequency.
    const m = /^(\d*\.\d+(?:[eE][-+]?\d+)?|\d+[eE][-+]?\d+)\s*(\(GE\))?$/.exec(f.trim());
    if (!m) continue;
    const af = Number(m[1]);
    if (!(af >= 0 && af <= 1)) continue;
    found.push({ af, populationDb: Boolean(m[2]) });
  }
  if (!found.length) return null;
  const tagged = found.find((x) => x.populationDb);
  if (tagged) return { ...tagged, assumed: false };
  return { ...found[0], assumed: true, otherCandidates: found.length - 1 };
}

const CHR_NAMES = new Set([...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y']);

/** Every reading of a "<chromosome><p|q><band>" token in `s` for the chromosomes in `want`.
 *  Reports are inconsistent about separating the two bands ("3q27.3;12p12.1" but also
 *  "3q27.312p12.1"), so each possible sub-band length is kept as a separate candidate and the
 *  ambiguity is resolved later by requiring the two chosen tokens not to overlap. */
function bandCandidates(s, want) {
  const out = [];
  for (let i = 0; i < s.length; i += 1) {
    if (i > 0 && /[A-Za-z]/.test(s[i - 1])) continue;         // mid-word: not the start of a band
    const j = /^chr/i.test(s.slice(i, i + 3)) ? i + 3 : i;
    for (const len of [2, 1]) {
      const chr = normChr(s.slice(j, j + len));
      if (!CHR_NAMES.has(chr) || !want.has(chr)) continue;
      const arm = (s[j + len] || '').toLowerCase();
      if (arm !== 'p' && arm !== 'q') continue;
      const first = j + len + 1;
      let k = first;
      while (k < s.length && s[k] >= '0' && s[k] <= '9') k += 1;
      if (k === first) continue;                              // a band needs at least one digit
      // Offer every prefix of the digit run, not just the whole of it: with no separator between
      // the two band fields ("1q4119p13.2") the run swallows the next chromosome's number, and only
      // the shorter reading (1q41 + 19p13.2) leaves room for both.
      for (let e = first + 1; e <= k; e += 1) {
        out.push({ chr, arm, start: i, end: e });             // band with no sub-band
        if (s[e] === '.') for (let d = 1; d <= 3; d += 1) {   // ".3", ".31", ".312" all possible
          const c = s[e + d];
          if (!(c >= '0' && c <= '9')) break;
          out.push({ chr, arm, start: i, end: e + d + 1 });
        }
      }
    }
  }
  return out;
}

/** The cytoband of each of a breakend's two chromosomes, read off the report's cytoband column
 *  (e.g. "3q27.3;12p12.1"), as {chr: {arm, band}}. Returns {} when the column is absent — the caller
 *  must then simply not name the derivative rather than guess. */
export function parseCytobands(text, chrA, chrB) {
  const s = String(text || '');
  const a = normChr(chrA); const b = normChr(chrB);
  if (a === b) return {};                                     // one chromosome: bands are not separable
  const cands = bandCandidates(s, new Set([a, b]));
  const forA = cands.filter((c) => c.chr === a);
  const forB = cands.filter((c) => c.chr === b);
  // `label` is the band as written, chromosome included (e.g. "14q32.33"); `arm` is just p or q.
  const entry = (c) => ({ arm: c.arm, label: s.slice(c.start, c.end).replace(/^chr/i, '') });
  let best = null;
  for (const x of forA) for (const y of forB) {
    if (x.end > y.start && y.end > x.start) continue;         // overlapping readings cannot both be right
    // prefer the leftmost pair, and within that the fullest sub-band reading
    const score = Math.min(x.start, y.start) * 1e6 + Math.abs(x.start - y.start) * 1000
      - (x.end - x.start) - (y.end - y.start);
    if (!best || score < best.score) best = { score, [a]: entry(x), [b]: entry(y) };
  }
  if (best) { const { score, ...bands } = best; return bands; }
  if (forA.length && !forB.length) return { [a]: entry(forA[forA.length - 1]) };
  if (forB.length && !forA.length) return { [b]: entry(forB[forB.length - 1]) };
  return {};
}

/** Just the chromosome arm (p or q) per chromosome. The arm is what tells us which side of a
 *  breakpoint carries the centromere, and hence which derivative chromosome a junction builds. */
export function parseCytobandArms(text, chrA, chrB) {
  const bands = parseCytobands(text, chrA, chrB);
  return Object.fromEntries(Object.entries(bands).map(([c, v]) => [c, v.arm]));
}

/** Split a pasted report into breakends. Each non-empty line may be a bare BND cell or a full
 *  tab/multi-space-separated row; we pull the first BND-looking token and any PR/SR support.
 *
 *  Lines that do not yield a breakend are returned in `unparsed` rather than dropped: silence would
 *  otherwise stand for "nothing there" when it actually means "this row was not understood". */
export function parseReport(text) {
  const out = [];
  out.unparsed = [];
  for (const line of String(text || '').split(/\r?\n/)) {
    if (!line.trim()) continue;
    const fields = line.split(/\t|\s{2,}/).map((f) => f.trim()).filter(Boolean);
    // Also try a plain-whitespace split, so a single-space-separated paste still finds its cell.
    const tokens = [...new Set([...fields, ...line.trim().split(/\s+/), line.trim()])];
    let bnd = null; let cellRaw = null;
    for (const f of tokens) {
      const p = parseBnd(f);
      if (p) { bnd = p; cellRaw = f; break; }
    }
    if (!bnd) { out.unparsed.push(line.trim()); continue; }
    // Look for the cytoband only in the fields that are not the BND cell itself, so a coordinate
    // such as "12:25056225" inside the notation can never be read as a band.
    const rest = fields.filter((f) => f !== cellRaw).join(' ');
    const bands = parseCytobands(rest, bnd.localChr, bnd.mateChr);
    out.push({
      ...bnd,
      cell: cellRaw,
      support: parseSupport(line),
      germline: parseGermlineAf(line),
      bands,
      arms: Object.fromEntries(Object.entries(bands).map(([c, v]) => [c, v.arm])),
    });
  }
  return out;
}

/** The unordered chromosome pair of a breakend, as a stable key "A|B" (sorted). */
export function chrPairKey(b) {
  return [b.localChr, b.mateChr].sort().join('|');
}

/** The two ends of a breakend, each as {chr, pos, side}. */
function endsOf(b) {
  return [
    { chr: b.localChr, pos: b.localPos, side: b.localSide },
    { chr: b.mateChr, pos: b.matePos, side: b.mateSide },
  ];
}

/** The end of a breakend that sits on `chr`; null if the breakend is intra-chromosomal (both ends on
 *  the same chromosome), where "the end on chr X" is not a single thing. */
function endOnChr(b, chr) {
  const [l, m] = endsOf(b);
  if (l.chr === m.chr) return null;
  return l.chr === chr ? l : (m.chr === chr ? m : null);
}

/**
 * A canonical, order-independent key for the *adjacency* a breakend describes. Callers normally emit
 * both mate records of one junction ("3:100T]12:200]" and "12:200T]3:100]"); those are one junction
 * seen from its two ends, not two junctions, and they share this key.
 */
export function adjacencyKey(b) {
  return endsOf(b).map((e) => `${e.chr}:${e.pos}:${e.side}`).sort().join('~');
}

const OPPOSITE_SIDE = { left: 'right', right: 'left' };

/** Fill in the chromosomes a per-chromosome map is missing from another map, without overwriting
 *  anything the first record already gave. A new object is returned every time: the merged record
 *  shares its maps with the record it was copied from, so filling them in place would edit the
 *  caller's parsed rows. */
function fillByChr(seen, extra) {
  const out = { ...(seen || {}) };
  for (const [chr, v] of Object.entries(extra || {})) if (out[chr] === undefined) out[chr] = v;
  return out;
}

/** Collapse mate records of the same adjacency into one breakend, keeping whichever fields are
 *  populated (a caller sometimes puts the support figures on only one of the two rows). */
export function dedupeMateRecords(breakends) {
  const byKey = new Map();
  const out = [];
  for (const b of breakends) {
    const k = adjacencyKey(b);
    const seen = byKey.get(k);
    if (!seen) { const copy = { ...b, mateRecords: [b.raw] }; byKey.set(k, copy); out.push(copy); continue; }
    seen.mateRecords.push(b.raw);
    seen.support = seen.support || b.support;
    seen.germline = seen.germline || b.germline;
    // Bands and arms are merged per chromosome, not all-or-nothing. Bands used to be dropped
    // altogether: the arms came across from the mate row and named the derivative while every band
    // label — headline, evidence row, the "this partner is at 14q32" prompt — silently vanished.
    // Per chromosome because a row can resolve one side's band and not the other's.
    seen.bands = fillByChr(seen.bands, b.bands);
    seen.arms = fillByChr(seen.arms, b.arms);
  }
  return out;
}

/**
 * Detect reciprocal junction pairs among the breakends.
 *
 * A reciprocal translocation cuts each chromosome once and swaps the pieces, so its two junctions
 * use *opposite* sides of the (near-identical) cut on BOTH chromosomes: if junction 1 keeps
 * chr3-left + chr12-left, junction 2 must keep chr3-right + chr12-right. That side-opposition test
 * is the definition; testing the bracket characters instead (as this used to) both misses cases and
 * happily pairs a record with its own mate record, which is the same junction twice.
 *
 * Returns { events } where each event is either a reciprocal pair or a lone junction. Offsets
 * between the two junctions are reported so the tightness of the call is visible (junction-associated
 * deletions can offset one side by kilobases).
 */
export function detectEvents(breakends, opts = {}) {
  const window = opts.window ?? 50000;
  const tight = opts.tight ?? 1000;
  const items = dedupeMateRecords(breakends).map((b, i) => ({ b, i, used: false }));
  const events = [];
  const candidates = [];
  const pairs = new Map();      // index of the first breakend of a pair -> the pairing
  const claimed = new Set();    // indexes that are the SECOND breakend of a pair

  for (let a = 0; a < items.length; a += 1) {
    const A = items[a].b;
    for (let c = a + 1; c < items.length; c += 1) {
      const B = items[c].b;
      if (chrPairKey(A) !== chrPairKey(B)) continue;
      const [x, y] = [A.localChr, A.mateChr];
      const ax = endOnChr(A, x); const ay = endOnChr(A, y);
      const bx = endOnChr(B, x); const by = endOnChr(B, y);
      if (!ax || !ay || !bx || !by) continue;                              // intra-chromosomal: skip
      if (bx.side !== OPPOSITE_SIDE[ax.side]) continue;                    // not the reciprocal cut
      if (by.side !== OPPOSITE_SIDE[ay.side]) continue;
      const dx = Math.abs(ax.pos - bx.pos);
      const dy = Math.abs(ay.pos - by.pos);
      // BOTH offsets must be inside the window. Admitting a pair because one side happens to be
      // tight lets a junction 127 Mb away on the same chromosome pair pass as "reciprocal".
      if (dx > window || dy > window) continue;
      candidates.push({ a, c, A, B, dx, dy, x, y, score: dx + dy });
    }
  }
  // Best-first, not first-row-wins: a loose but earlier candidate must not consume a breakend that
  // a tighter pair needs.
  candidates.sort((p, q) => p.score - q.score || p.a - q.a || p.c - q.c);
  for (const k of candidates) {
    if (items[k.a].used || items[k.c].used) continue;
    items[k.a].used = true; items[k.c].used = true;
    pairs.set(k.a, k); claimed.add(k.c);
  }
  for (let i = 0; i < items.length; i += 1) {
    const k = pairs.get(i);
    if (k) {
      events.push({
        reciprocal: true,
        chrs: chrPairKey(k.A).split('|'),
        breakends: [k.A, k.B],
        offsets: { [k.x]: k.dx, [k.y]: k.dy },
        confidence: Math.max(k.dx, k.dy) <= tight ? 'tight' : 'loose',
      });
    } else if (!claimed.has(i)) {
      const A = items[i].b;
      events.push({ reciprocal: false, chrs: chrPairKey(A).split('|'), breakends: [A], offsets: null, confidence: null });
    }
  }
  return { events };
}

/**
 * Name the gene(s) at a breakpoint. Returns any gene whose body spans the position, and — when the
 * breakpoint falls just outside a gene — the nearest gene with its distance and side (5' upstream /
 * 3' downstream, strand-aware). A breakpoint 5' of a gene is how promoter-capture rearrangements
 * present (e.g. a breakpoint ~0.7 kb 5' of BCL6). `genes` is the map from gene-loci.json.
 */
export function annotateLocus(genes, chr, pos, opts = {}) {
  const window = opts.window ?? 1000000;
  const c = normChr(chr);
  const inGene = [];
  let nearest = null;
  for (const sym of Object.keys(genes || {})) {
    const e = genes[sym];
    if (!e || e.chr !== c) continue;
    if (pos >= e.start && pos <= e.end) {
      inGene.push({ gene: sym, strand: e.strand, tx: e.tx, span: e.end - e.start });
      continue;
    }
    const gap = pos < e.start ? e.start - pos : pos - e.end;
    if (gap > window) continue;
    // 5' of a + gene = below its start; 5' of a - gene = above its end.
    const side = pos > e.end ? (e.strand === '-' ? '5prime' : '3prime') : (e.strand === '-' ? '3prime' : '5prime');
    if (!nearest || gap < nearest.distance || (gap === nearest.distance && sym < nearest.gene)) {
      nearest = { gene: sym, strand: e.strand, tx: e.tx, distance: gap, side };
    }
  }
  // Overlapping gene bodies are common. Order them deterministically — largest span first, then
  // alphabetically — rather than by whatever order the loci file happened to be written in. Largest
  // first because the small genes tucked inside a driver (DEXI inside CIITA, for one) would
  // otherwise take the name and invert the donor/acceptor call.
  inGene.sort((p, q) => q.span - p.span || (p.gene < q.gene ? -1 : 1));
  return { inGene, nearest };
}

/**
 * Every gene near a breakpoint, with the part of it that survives on the kept side of the junction,
 * best candidate first. `keepLeft` says which side survives.
 *
 * Ordering is the point of this function. A gene whose body is carried across on the kept side is
 * far more relevant than a nearer gene that is left behind entirely — a breakpoint 5 kb upstream of
 * BCL2 that carries the whole of BCL2 to an immunoglobulin locus must name BCL2, not whichever gene
 * happens to sit 3 kb away on the discarded side.
 */
export function geneCandidates(genes, chr, pos, keepLeft, opts = {}) {
  const window = opts.window ?? 50000;
  const c = normChr(chr);
  const out = [];
  for (const sym of Object.keys(genes || {})) {
    const g = genes[sym];
    if (!g || g.chr !== c) continue;
    const inside = pos >= g.start && pos <= g.end;
    const gap = inside ? 0 : (pos < g.start ? g.start - pos : pos - g.end);
    if (gap > window) continue;
    // Gene spans are 1-based inclusive, so the body is end - start + 1 bases and the kept side has
    // to be counted inclusively too. Counting interval lengths instead left the two ends of a gene
    // inconsistent: a breakpoint exactly on the last base kept the whole gene ("intact"), while one
    // exactly on the first base kept none of it and was labelled 5′ flank at a distance of 0 —
    // a flank the breakpoint is sitting inside.
    const span = g.end - g.start + 1;
    // bases of the gene body that survive on the kept side
    const kept = keepLeft
      ? Math.max(0, Math.min(pos, g.end) - g.start + 1)
      : Math.max(0, g.end - Math.max(pos, g.start) + 1);

    // A breakpoint sitting exactly ON the gene's first or last base, with the rest of the gene on the
    // discarded side, keeps one base of it. That is the gene's boundary, not a cut through it: the
    // kept side is flanking sequence and the honest thing to say is where the breakpoint sits, not
    // that the gene contributes a 1-base "5′ part". Calling it a cut also let such a gene outrank a
    // gene carried across whole — the inverse of the rule this ranking exists to enforce.
    const atBoundaryBase = span > 1 && kept === 1 && (pos === g.start || pos === g.end);

    let part; let extra = {};
    if (kept <= 0 || atBoundaryBase) {
      // Nothing of the gene survives; the kept side is flanking sequence. Which flank it is follows
      // from which side is kept, not from comparing the breakpoint to the gene's end — at a position
      // exactly on the boundary those two disagree.
      const keptSideIsAbove = !keepLeft;
      const fivePrimeFlank = keptSideIsAbove ? g.strand === '-' : g.strand === '+';
      part = fivePrimeFlank ? 'upstream' : 'downstream';
      // `onGeneBoundary` says the breakpoint is on the gene's own first or last base rather than some
      // distance outside it — which is what "distance 0" used to mean and could not express.
      extra = { distance: gap, onGeneBoundary: atBoundaryBase || gap === 0 };
    } else if (kept >= span) {
      // the whole gene body survives; which end the junction sits beyond is pure geometry
      part = 'intact';
      extra = {
        junctionSide: keepLeft ? (g.strand === '-' ? '5prime' : '3prime') : (g.strand === '-' ? '3prime' : '5prime'),
        distance: gap,
      };
    } else {
      const fivePrime = g.strand === '+' ? keepLeft : !keepLeft;
      part = fivePrime ? '5prime' : '3prime';
      extra = { span, kept, keptFraction: span > 0 ? kept / span : null };
    }
    // 0 = the breakpoint cuts this gene, 1 = the gene is carried across whole, 2 = only its flank.
    // A gene the breakpoint merely borders ranks with the flanks, so it can never outrank a gene the
    // derivative actually carries: rank by what survives on the kept side, not by distance.
    const rank = (kept <= 0 || atBoundaryBase) ? 2 : (kept >= span ? 1 : 0);
    out.push({ gene: sym, strand: g.strand, tx: g.tx, span, kept, gap, rank, part, ...extra });
  }
  out.sort((p, q) => p.rank - q.rank
    || (p.rank === 0 ? q.span - p.span : p.gap - q.gap)
    || (p.gene < q.gene ? -1 : 1));
  return out;
}

/** A short human label for a breakpoint's gene context. Only names a nearby gene within 50 kb, so an
 *  IG/TR-locus or intergenic breakpoint is not mislabelled after a distant MANE gene. */
export function locusLabel(ann) {
  if (ann.inGene.length) return { gene: ann.inGene[0].gene, note: 'within gene', promoter: false };
  if (ann.nearest && ann.nearest.distance <= 50000) {
    const n = ann.nearest;
    const promoter = n.side === '5prime' && n.distance <= 20000;
    const where = `${n.distance.toLocaleString()} bp ${n.side === '5prime' ? '5′' : '3′'}`;
    return { gene: n.gene, note: `${where}${promoter ? ' (promoter/regulatory)' : ''}`, promoter };
  }
  return { gene: null, note: 'no MANE gene', promoter: false };
}

/**
 * What a kept derivative segment actually carries, relative to the gene at its breakpoint. This is the
 * clinically meaningful label: a breakpoint 697 bp 5' of BCL6 whose gene-body side is kept carries the
 * INTACT gene, not a "promoter fragment"; a breakpoint inside IRAG2 whose 5' side is kept carries that
 * gene's promoter + early exons. seg = {chr, pos, side:'left'|'right'}.
 *   part: 'intact' | '5prime' | '3prime' | 'upstream' (5'-flank/promoter) | 'downstream' | 'none'
 */
export function segmentContent(genes, seg, opts = {}) {
  const c = normChr(seg.chr); const keepLeft = seg.side === 'left';
  const cands = geneCandidates(genes, c, seg.pos, keepLeft, { window: opts.nameWindow ?? 50000 });
  if (!cands.length) {
    // Say how far the nearest gene is, so a true gene desert is distinguishable from a breakpoint
    // that fell just outside the naming window. Never fall back to that gene's name.
    const wider = geneCandidates(genes, c, seg.pos, keepLeft, { window: 2000000 });
    const nearest = wider.slice().sort((p, q) => p.gap - q.gap)[0] || null;
    // An immunoglobulin or T-cell-receptor locus is not in MANE Select yet is the clinically
    // important partner in most lymphoma translocations; when the caller supplies the locus table
    // (data/ig-loci.json), a breakpoint inside one is named instead of shrugged at. TRD nests inside
    // TRA, so the SMALLEST containing locus wins.
    const igLocus = (opts.igLoci || [])
      .filter((l) => normChr(l.chr) === c && seg.pos >= l.start && seg.pos <= l.end)
      .sort((p1, q1) => (p1.end - p1.start) - (q1.end - q1.start))[0] || null;
    return {
      gene: null,
      part: 'none',
      igLocus,
      alternatives: [],
      nearestGene: nearest ? { gene: nearest.gene, distance: nearest.gap } : null,
      label: igLocus ? `${igLocus.symbol} locus` : `chr${c}`,
    };
  }
  const [best, ...rest] = cands;
  // A real gene embedded wholly inside an IG/TCR locus (IGLL5 inside IGL) must not out-name the
  // locus: V(D)J-mediated breaks belong to the locus, and naming the bystander gene misleads. The
  // embedded gene is surfaced as context instead.
  const igAll = (opts.igLoci || [])
    .filter((l) => normChr(l.chr) === c && seg.pos >= l.start && seg.pos <= l.end)
    .sort((p1, q1) => (p1.end - p1.start) - (q1.end - q1.start));
  if (igAll.length) {
    const locus = igAll[0];
    const g = genes[best.gene];
    if (g && g.start >= locus.start && g.end <= locus.end) {
      return {
        gene: null, part: 'none', igLocus: locus, embeddedGene: best.gene,
        alternatives: [], nearestGene: null, label: `${locus.symbol} locus`,
      };
    }
  }
  const label = {
    intact: `intact ${best.gene}`,
    '5prime': `${best.gene} 5′ part (promoter + 5′ exons)`,
    '3prime': `${best.gene} 3′ part (3′ exons)`,
    upstream: `5′-flank/promoter of ${best.gene}`,
    downstream: `3′ of ${best.gene}`,
  }[best.part];
  // An alternative is a gene that genuinely competes to be THE gene at this junction: one the
  // breakpoint falls inside, or one carried across whole from within 20 kb. A gene 40 kb away that
  // happens to ride along on the same fragment is not an alternative reading, it is noise.
  const alternatives = rest.filter((r) => r.rank === 0 || (r.rank === 1 && r.gap <= 20000)).slice(0, 3)
    .map((r) => ({ gene: r.gene, part: r.part, strand: r.strand, gap: r.gap, inside: r.rank === 0 }));
  return { ...best, gene: best.gene, alternatives, alsoIn: alternatives.map((a) => a.gene), label };
}

/** Short part label for a segment content, for compact display. */
export function partShort(part) {
  return ({ intact: 'intact', '5prime': '5′ part', '3prime': '3′ part', upstream: '5′ promoter', downstream: '3′ flank', none: 'no MANE gene' })[part] || part;
}

/** The novel sequence at the junction. The t field is the reference base at the record's own
 *  position PLUS any insertion, so one base has to come off: the reference base leads in the t-first
 *  forms and trails in the t-last ones. Reporting t verbatim overstates a 3 bp insert as 4 bp.
 *  "(+N)" is the caller's shorthand for N elided bases and is passed through as written. */
function cleanInsertion(bases, tPosition) {
  const s = String(bases || '').trim();
  if (s.length <= 1) return '';
  return tPosition === 'last' ? s.slice(0, -1) : s.slice(1);
}

// =====================================================================
//  READING A DERIVATIVE IN TRANSCRIPT ORDER
// =====================================================================
// A derivative chromosome is double-stranded: it has no inherent left or right, and the order the
// VCF happens to write the two pieces in is an artefact of which chromosome the record was filed
// under. What a clinician needs is the order a *transcript* would read them in, 5′ to 3′. That is
// what the functions below compute: lay the two pieces out, work out which way each gene points in
// the rearranged molecule, and — if they both point the same way — read the derivative in that
// direction, promoter donor first.

const SIGN = { '+': 1, '-': -1 };

/**
 * The two pieces this breakend joins, in the order the ALT form writes them, each with the direction
 * its reference sequence is traversed ('+' = ascending coordinate, '-' = descending).
 *
 * A junction glues the END of the first piece to the START of the second, and that alone fixes the
 * directions: a "left" piece (pter..breakpoint) ends at its breakpoint when read forwards, so as the
 * first piece it is traversed '+' and as the second piece '-'; a "right" piece is the mirror image.
 */
export function derivativeSegments(b) {
  const local = { chr: b.localChr, pos: b.localPos, side: b.localSide, from: 'local' };
  const mate = { chr: b.mateChr, pos: b.matePos, side: b.mateSide, from: 'mate' };
  const [first, second] = b.order === 'mate-first' ? [mate, local] : [local, mate];
  return [
    { ...first, orient: first.side === 'left' ? '+' : '-' },
    { ...second, orient: second.side === 'right' ? '+' : '-' },
  ];
}

/**
 * Which way the gene on a piece points once the piece sits in the derivative: +1 = left to right in
 * the reading, -1 = right to left, null = no gene. A gene on the reference minus strand carried on a
 * piece that is itself traversed backwards ends up pointing forwards again, which is exactly the
 * case that made the old tool read BCL6 the wrong way round.
 */
export function segmentTxDirection(seg) {
  if (!seg || !seg.content || !seg.content.gene) return null;
  return SIGN[seg.orient] * SIGN[seg.content.strand];
}

/** Read the same molecule from its other end: reverse the pieces and flip every traversal. */
function flipReading(segs) {
  return segs.slice().reverse().map((s) => ({ ...s, orient: s.orient === '+' ? '-' : '+' }));
}

/** Does this piece carry its chromosome's centromere? The centromere lies between the arms, so it is
 *  below a q-arm breakpoint and above a p-arm one. `arm` comes from the report's cytoband column;
 *  without it the answer is unknown and the derivative simply goes unnamed. */
function carriesCentromere(seg, arm) {
  if (arm !== 'p' && arm !== 'q') return null;
  return seg.side === 'left' ? arm === 'q' : arm === 'p';
}

/**
 * Name the derivative chromosome a junction builds — der(3), der(12) — from which piece carries
 * which centromere. Every normal derivative has exactly one; none means an acentric fragment (which
 * would be lost at mitosis) and two means a dicentric, both of which are worth saying out loud.
 */
export function nameDerivative(segs, arms = {}) {
  const carried = segs.map((s) => ({ chr: s.chr, has: carriesCentromere(s, (arms || {})[s.chr]) }));
  if (carried.some((c) => c.has === null)) {
    return { name: null, basis: 'the pasted row carries no cytoband, so the centromere-bearing piece cannot be identified' };
  }
  const withCentromere = carried.filter((c) => c.has).map((c) => c.chr);
  if (withCentromere.length === 1) {
    return { name: `der(${withCentromere[0]})`, basis: `this piece carries the chr${withCentromere[0]} centromere` };
  }
  if (withCentromere.length === 0) {
    return { name: null, warning: 'acentric',
      basis: 'neither piece carries a centromere, so this junction on its own gives an acentric fragment. In a balanced exchange that is normal: the centromere-bearing product is the other, reciprocal junction' };
  }
  return { name: null, warning: 'dicentric', basis: `both pieces carry a centromere (chr${withCentromere.join(' and chr')}) — a dicentric product` };
}

const READ_THROUGH_DONOR = new Set(['5prime', 'upstream']);

/**
 * What a donor's retained stop codon means for the gene behind it, with the assumption stated.
 *
 * The computed part is the geometry: the donor keeps every coding base, so its stop codon sits in
 * the fused message ahead of the junction. What follows from that is an assumption about the
 * ribosome — that it stops there and does not start again further down the same message — so the
 * sentence says "unless a ribosome restarts" rather than presenting it as a settled fact.
 */
const stopEndsTranslation = (donor, acceptor) => `Translation from ${donor}'s promoter ends at that `
  + `stop, so ${acceptor} is not translated unless a ribosome restarts.`;

/**
 * Classify what the read-through, if any, produces. `five` and `three` are the 5′ and 3′ pieces of
 * the derivative *as read*, with their `content` from segmentContent().
 *
 * The tool holds gene spans only — no exon or coding-sequence coordinates — so it never claims a
 * reading frame. Anything frame-dependent is reported as a candidate with that limit stated.
 */
/**
 * Exon-level facts about one segment, when the exon table is supplied. Absent it, every consumer
 * falls back to the span-only reasoning — the tool degrades to what it used to say rather than
 * guessing. `keepLeft` is genomic; codingRetained resolves it against the gene's strand.
 */
export function exonFactsFor(seg, exonTable) {
  const sym = seg.content && seg.content.gene;
  const g = sym && exonTable ? exonTable[sym] : null;
  if (!g) return null;
  const keepLeft = seg.side === 'left';
  const where = featureAt(g, seg.pos);
  const kept = codingRetained(g, seg.pos, keepLeft);
  const lost = codingOffsetOfAcceptor(g, seg.pos, keepLeft);
  let cdsTotal = 0;
  for (let i = 0; i + 1 < (g.cdsExons || []).length; i += 2) cdsTotal += g.cdsExons[i + 1] - g.cdsExons[i] + 1;
  return {
    tx: g.tx,
    exons: exonCount(g),
    where,                                   // {kind:'exon'|'intron'|'outside', …}
    label: featureLabel(g, seg.pos),
    codingKept: kept.bases,
    codingLost: lost,
    codingTotal: cdsTotal,
    codingSide: kept.side,
    // The break spares the whole coding sequence, so nothing of the protein is lost on this side.
    // Which side of the coding sequence the break lies on is `codingSide`, and `where` says whether
    // it is exonic or intronic — so there is no separate "untranslated region" flag. There used to
    // be one, and it read true for intronic breaks as well, which is not what its name said. It was
    // read nowhere, and a wrong flag nobody reads is still a trap for whoever reads it next.
    codingIntact: cdsTotal > 0 && kept.bases === cdsTotal,
  };
}

/**
 * Does the acceptor piece keep its entire protein-coding sequence?
 *
 * This is the question the span-based part label ('3prime' = the breakpoint fell somewhere inside the
 * gene) cannot answer. It is true only when the exon table is loaded AND every coding base survives on
 * the kept side AND that kept side is the gene's 3′ portion in transcript terms — which places the
 * breakpoint 5′ of the start codon, in the 5′ untranslated region or an intron ahead of it. The gene
 * then loses its promoter and gains the donor's, with the protein unchanged.
 *
 * Returns false whenever exon facts are absent, so the classifier degrades to its span-only wording
 * rather than guessing.
 */
export function acceptorCodingIntact(bx) {
  return Boolean(bx && bx.codingIntact && bx.codingSide === '3prime' && bx.codingTotal > 0);
}

/**
 * The donor's mirror image: does the 5′ piece keep the gene's entire protein-coding sequence?
 *
 * True only with the exon table loaded, every coding base on the kept side, and that kept side the
 * gene's 5′ portion — which puts the breakpoint 3′ of the stop codon, in the 3′ untranslated region
 * or beyond it. Two things follow, and both are facts rather than proxies: the donor's protein is
 * unchanged, and its stop codon sits in front of the junction, so nothing downstream is translated.
 *
 * False whenever exon facts are absent, so the classifier degrades to its span-only wording.
 */
export function donorCodingIntact(ax) {
  return Boolean(ax && ax.codingIntact && ax.codingSide === '5prime' && ax.codingTotal > 0);
}

/**
 * Can the donor's protein fate be settled from exon facts at all? It needs the exon table, a coding
 * sequence in it, and the kept side to be the gene's 5′ portion. When this is false the classifier
 * falls back to the genomic-span proxy and says so.
 */
function donorCodingDecidable(ax) {
  return Boolean(ax && ax.codingTotal > 0 && ax.codingSide === '5prime');
}

export function classifyProduct(five, three, ctx = {}) {
  const bands = ctx.bands || {};
  const a = five.content; const b = three.content;
  // Exon-level facts, when the exon table was supplied. Everything below tests for these before
  // using them, so the classifier still works span-only when they are absent.
  const ax = five.exon || null; const bx = three.exon || null;
  // junctionFrame returns null for TWO different reasons, and they must never be reported as one:
  // the exon table is missing for a gene, or it is present and shows that one side retains no coding
  // bases at all. The second is an answer, not a gap — saying "no exon coordinates are loaded" there
  // is a false statement about the tool's own data. This flag separates them.
  const exonTablesLoaded = Boolean(ctx.exons && a.gene && b.gene && ctx.exons[a.gene] && ctx.exons[b.gene]);
  const frame = exonTablesLoaded
    ? junctionFrame(ctx.exons[a.gene], five.pos, five.side === 'left', ctx.exons[b.gene], three.pos, three.side === 'left')
    : null;
  // "This gene has a coding sequence in the table and none of it survives on the kept side."
  const contributesNoCoding = (x) => Boolean(x && x.codingTotal > 0 && x.codingKept === 0);
  const nonCodingInTable = (x) => Boolean(x && x.codingTotal === 0);
  const ap = a.gene ? a.part : 'none';
  const bp = b.gene ? b.part : 'none';
  const caveats = [];
  const name5 = a.gene || `chr${five.chr}`;
  const name3 = b.gene || `chr${three.chr}`;
  // How much of a split gene's genomic span survives on this side, in plain words.
  const keptPhrase = (c) => (c.keptFraction == null ? ''
    : `${Math.round(c.keptFraction * 100)}% of its ${Math.round(c.span / 1000)} kb genomic span`);
  const keptText = (c) => (c.keptFraction == null ? '' : ` (${keptPhrase(c)})`);

  // Overlapping genes are common and the 5′/3′ call is strand-dependent, so a breakpoint that sits
  // in more than one gene has more than one reading. Say so every time rather than picking silently.
  for (const [seg, content] of [[five, a], [three, b]]) {
    if (!content.alternatives || !content.alternatives.length) continue;
    const where = (x) => (x.inside
      ? `${x.gene} (${x.strand} strand), which the breakpoint also falls inside — ${partShort(x.part)} here`
      : `${x.gene} (${x.strand} strand) ${(x.gap / 1000).toFixed(1)} kb away, also carried across ${partShort(x.part)}`);
    caveats.push(`At chr${seg.chr}:${seg.pos.toLocaleString()} this piece could equally be read as ${content.alternatives.map(where).join(', or ')}. ${content.gene} was taken as the primary gene because it has the widest span at this position; on the other reading the donor and acceptor roles can swap.`);
  }

  if (a.gene && a.gene === b.gene) {
    return { kind: 'same-gene', productive: false,
      headline: `Both breakpoints fall in ${a.gene}`,
      summary: `Internal rearrangement of ${a.gene}, not a fusion with a partner.`,
      caveats };
  }

  if (READ_THROUGH_DONOR.has(ap) && bp === 'intact') {
    // Same geometry as the donor-stop-codon-retained branch below, with an intact acceptor instead
    // of a cut one: the donor keeps every coding base, so the breakpoint is 3′ of its stop codon and
    // that stop sits in the fused message in front of the acceptor. Two things the span-based
    // wording used to say here are then both false — the donor's allele is NOT disrupted (its whole
    // reading frame and stop codon survive) and the acceptor is NOT translated from this message.
    if (donorCodingIntact(ax)) {
      return { kind: 'donor-stop-codon-retained', productive: false,
        headline: `${name5} → intact ${name3} — translation ends at ${name5}'s own stop codon`,
        summary: [
          `No fusion protein: ${name5} keeps its stop codon, cut ${ax.label}, past it.`,
          `${name5} protein unchanged — all ${ax.codingTotal.toLocaleString()} coding bases kept.`,
          `${name3} is intact${bx && bx.codingIntact ? ` — all ${bx.codingTotal.toLocaleString()} coding bases across ${bx.exons} exons` : ''}.`,
          stopEndsTranslation(name5, name3),
        ].join(' '),
        caveats };
    }
    const donorText = ap === '5prime'
      ? `${name5}'s promoter and 5′ region${keptText(a)}`
      : `${name5}'s 5′ regulatory region (promoter side, no exons of ${name5})`;
    // Only genuine limitations belong in caveats. That this allele of the donor is disrupted is a
    // FINDING, so it goes in the summary; and with exon coordinates loaded the donor-exon question is
    // answered rather than deferred, so the caveat is raised only when they are absent.
    if (!ax) caveats.push(`Whether any ${name5} exons are transcribed into the ${name3} message cannot be settled without exon coordinates.`);
    // "Promoter replaced" was printed unconditionally, for an acceptor up to 50 kb from the junction.
    // The acceptor is intact and the junction lies 5′ of it, so everything between the two — its own
    // upstream region, proximal promoter and all — rides across on the same derivative. What is
    // computed is that distance, so that is what gets said; whether the acceptor's own promoter is
    // inside it cannot be known here, because the tool holds no promoter annotation.
    //
    // What the number measures matters: it is the gap from the junction to the acceptor's gene body,
    // i.e. how much of the ACCEPTOR's own upstream sequence came across. It is not the distance to the
    // donor's promoter, which lies somewhere beyond the junction and is not computed here at all. The
    // sentence used to attach the first number to the second thing, which reads as a measurement the
    // tool never made.
    const upstreamKept = Number.isFinite(b.distance) ? b.distance : null;
    const promoterSentence = upstreamKept
      ? `The junction lies ${upstreamKept.toLocaleString()} bp 5′ of ${name3}; ${name5}'s promoter is further 5′ again, at a distance not computed here.`
      : `${name3} keeps none of its own upstream sequence; ${name5}'s promoter now leads it.`;
    if (upstreamKept) {
      caveats.push(`${name3} carries ${upstreamKept.toLocaleString()} bp of its own upstream sequence across too, so whether ${name5}'s promoter replaces ${name3}'s or is merely added in front of it is not settled here — this tool holds no promoter annotation.`);
    }
    // The donor's fate follows its own numbers. With the breakpoint in the donor's 5′ flank the gene
    // body is not on this piece at all, so it is not "cut" anywhere and saying so would be false.
    const donorFate = ap === 'upstream'
      ? `${name5}'s own gene body stays behind — this piece carries only its 5′ regulatory sequence.`
      : ax
        ? (ax.codingTotal > 0
          // The denominator matters: "501 coding bases carried across" with a total of 501 reads as
          // a loss when nothing was lost. That case now returns above, and this one prints both.
          ? `${name5} allele disrupted: cut ${ax.label}, ${ax.codingKept.toLocaleString()} of ${ax.codingTotal.toLocaleString()} coding bases carried across.`
          : `${name5} allele disrupted: cut ${ax.label}; the exon table lists no coding sequence for it.`)
        : `${name5} allele disrupted at the breakpoint.`;
    return { kind: 'promoter-substitution', productive: true,
      headline: `${donorText} → intact ${name3}`,
      summary: [
        `${name3} protein unchanged${bx && bx.codingIntact ? ` — ${bx.codingTotal.toLocaleString()} coding bases across ${bx.exons} exons intact` : ''}.`,
        promoterSentence,
        'No chimeric protein.',
        donorFate,
      ].join(' '),
      caveats };
  }

  if (ap === '5prime' && bp === '3prime') {
    // A donor cut 3′ of its own stop codon carries that stop into the fused message, in front of the
    // junction. Translation starts at the donor's start codon and ends there, so no fusion protein
    // can form however neatly the two coding sequences line up modulo three. The frame arithmetic is
    // deliberately not attached to this product: there is nothing fused for it to describe.
    //
    // junctionFrame reaches the same conclusion from the same tables, so either flag is enough. The
    // exon-fact one is preferred because it carries the wording (which exon the cut falls in); the
    // frame flag is the backstop for a caller that passes ctx.exons without pre-computed exon facts,
    // and it keeps a frame object with null fields from ever reaching the chimeric wording below.
    const donorStop = donorCodingIntact(ax) || Boolean(frame && frame.donorStopRetained);
    if (donorStop) {
      const donorFacts = ax
        ? [`No fusion protein: ${name5} keeps its stop codon, cut ${ax.label}, past it.`,
          `${name5} protein unchanged — all ${ax.codingTotal.toLocaleString()} coding bases kept.`]
        : [`No fusion protein: ${name5} keeps its stop codon, and the breakpoint lies past it.`,
          `${name5} protein unchanged — all ${frame.donorCodingTotal.toLocaleString()} coding bases kept.`];
      return { kind: 'donor-stop-codon-retained', productive: false,
        headline: `${name5} 5′ part → ${name3} 3′ part — translation ends at ${name5}'s own stop codon`,
        summary: [
          ...donorFacts,
          acceptorCodingIntact(bx)
            ? `${name3} keeps its whole coding sequence but loses its own promoter.`
            : `${name3} has lost its 5′ end and promoter${bx ? `, cut ${bx.label}` : ''}.`,
          stopEndsTranslation(name5, name3),
        ].join(' '),
        caveats };
    }
    // With both exon tables loaded, a null frame is an ANSWER, not a gap: one side retains no coding
    // bases, so there is no fused reading frame to be in or out of. The branch below used to blame
    // missing exon coordinates for it — the opposite of what the tool knows — and then offered a
    // "chimeric transcript expected" product on numbers that rule the chimeric protein out.
    if (!frame && exonTablesLoaded && ax && bx) {
      // The donor's start codon is on the discarded side: its promoter and 5′ untranslated exons
      // ride across, none of its protein can.
      if (contributesNoCoding(ax) || nonCodingInTable(ax)) {
        const donorNote = nonCodingInTable(ax)
          ? `${name5} has no coding sequence in the exon table, so it contributes none.`
          : `${name5} contributes no coding bases — cut ${ax.label}, 5′ of its own start codon.`;
        // Every acceptor coding base surviving makes this the promoter-capture event, the same one
        // the span label misses: an unchanged protein read from a new promoter.
        if (acceptorCodingIntact(bx)) {
          return { kind: 'promoter-substitution', productive: true,
            headline: `${name5}'s promoter and 5′ exons → ${name3}, coding sequence complete`,
            summary: [
              `${name3} protein unchanged — all ${bx.codingTotal.toLocaleString()} coding bases across ${bx.exons} exons kept.`,
              `Cut ${bx.label}, 5′ of the start codon: only ${name3}'s own promoter and 5′ untranslated region are lost.`,
              donorNote,
              'No chimeric protein.',
            ].join(' '),
            caveats };
        }
        return { kind: 'no-fusion-protein', productive: false,
          headline: `${name5} 5′ part → ${name3} 3′ part — neither side brings a start codon`,
          summary: [
            donorNote,
            nonCodingInTable(bx)
              ? `${name3} has no coding sequence in the exon table either.`
              : `${name3} has lost its own start codon — cut ${bx.label}, ${bx.codingLost.toLocaleString()} of ${bx.codingTotal.toLocaleString()} coding bases lost.`,
            'A chimeric transcript still forms; no start codon this tool can see lies in it.',
          ].join(' '),
          caveats };
      }
      // The mirror: the acceptor keeps no coding bases, so the break lies 3′ of its stop codon and it
      // contributes untranslated sequence only. The donor is cut inside its own coding sequence — the
      // intact-donor case was settled above — so a truncated donor protein is what this can make.
      if (contributesNoCoding(bx) || nonCodingInTable(bx)) {
        return { kind: 'truncated-donor', productive: 'possible',
          headline: `${name5} 5′ part → ${name3}, which contributes no coding sequence`,
          summary: [
            `${name5} truncated: cut ${ax.label}, ${ax.codingKept.toLocaleString()} of ${ax.codingTotal.toLocaleString()} coding bases kept.`,
            nonCodingInTable(bx)
              ? `${name3} has no coding sequence in the exon table.`
              : `${name3} contributes no coding bases — the break lies 3′ of its stop codon.`,
            // "No fusion protein" is this file's flat verdict, and it belongs to the branches that
            // return productive:false. Here the donor's start codon and part of its coding sequence
            // are on the message, so a truncated donor protein can still be made — which is what
            // productive:'possible' says. Deny the chimera only, and say what can form.
            `No chimeric protein, but a truncated ${name5} protein is possible: translation runs off its kept coding bases into ${name3} sequence. Where it stops is not computed here.`,
          ].join(' '),
          caveats };
      }
    }
    if (!frame) {
      caveats.push(exonTablesLoaded
        ? 'Whether the join is in frame cannot be determined here: one side retains no coding sequence, so there is no fused frame to measure.'
        : 'Whether the join is in frame cannot be determined here: no exon coordinates are loaded for one or both genes.');
    } else {
      // Splicing is assumed to be canonical: an intronic break joins to the next exon boundary. That
      // assumption is stated rather than buried, because a cryptic splice site would break it.
      caveats.push(`The frame call assumes canonical splicing across the junction. ${ax && ax.where.kind === 'intron' ? 'This breakpoint is intronic, so the transcript is taken to splice from the last complete donor exon; ' : ''}a cryptic splice site would change it.`);
    }
    const where = [ax ? `${name5} is cut ${ax.label}` : null, bx ? `${name3} is cut ${bx.label}` : null].filter(Boolean).join(', and ');
    return { kind: 'chimeric-candidate', productive: frame ? (frame.inFrame ? true : 'possible') : 'possible',
      headline: frame
        ? `${name5} 5′ part → ${name3} 3′ part — ${frame.inFrame ? 'in frame' : `out of frame by ${frame.shift}`}`
        : `${name5} 5′ part → ${name3} 3′ part`,
      summary: frame
        ? [
          frame.inFrame ? 'In frame: fusion protein expected.' : `Out of frame by ${frame.shift}: premature stop, truncated product.`,
          `Fused coding sequence ${frame.donorCodingBases.toLocaleString()} + ${frame.acceptorCodingBases.toLocaleString()} = ${frame.fusedCodingBases.toLocaleString()} bases.`,
          where ? `${where.charAt(0).toUpperCase()}${where.slice(1)}.` : '',
        ].filter(Boolean).join(' ')
        : [
          exonTablesLoaded
            ? 'Chimeric transcript expected; one side retains no coding sequence, so there is no fused frame to measure.'
            : 'Chimeric transcript expected; frame not determinable without exon coordinates.',
          `${name5} keeps ${keptPhrase(a)}; ${name3} keeps ${keptPhrase(b)}.`,
        ].join(' '),
      frame,
      caveats };
  }

  if (ap === 'upstream' && bp === '3prime') {
    // "3′ part" is a SPAN label: it only says the breakpoint fell inside the gene's genomic span. It
    // does not say the start codon went with it. When the break lands 5′ of the coding sequence — in
    // the 5′ untranslated region or an intron ahead of the start codon — the whole open reading frame
    // survives and only the gene's own promoter is lost. That is promoter substitution, and it is the
    // t(3;14) IGH-BCL6 geometry. Decide it from the exon facts whenever they are loaded.
    if (acceptorCodingIntact(bx)) {
      return { kind: 'promoter-substitution', productive: true,
        headline: `${name5}'s 5′ regulatory region (promoter side, no exons of ${name5}) → ${name3}, coding sequence complete`,
        summary: [
          `${name3} protein unchanged — all ${bx.codingTotal.toLocaleString()} coding bases across ${bx.exons} exons kept.`,
          `Cut ${bx.label}, 5′ of the start codon: only ${name3}'s own promoter and 5′ untranslated region are lost.`,
          `Promoter replaced by ${name5}'s. No chimeric protein.`,
        ].join(' '),
        caveats };
    }
    return { kind: 'headless-acceptor', productive: false,
      headline: `5′ regulatory region of ${name5} → ${name3} 3′ part`,
      // Without exon coordinates the start codon's fate is unknown, so it is not asserted: the loss
      // of the 5′ end and promoter is what the spans actually show.
      summary: `No product. ${name5} promoter carries none of its exons; ${name3} has lost its 5′ end and promoter${bx ? `, cut ${bx.label} — ${bx.codingLost.toLocaleString()} of ${bx.codingTotal.toLocaleString()} coding bases lost, including the start codon` : ''}.`,
      caveats };
  }

  if (READ_THROUGH_DONOR.has(ap)) {                       // reads into no gene, or into a 3′ flank
    const partnerBand = bands[three.chr] ? bands[three.chr].label : null;
    const igl3 = b.igLocus || null;
    const into = bp === 'downstream' ? `sequence 3′ of ${name3}`
      : igl3 ? `the ${igl3.symbol} locus (${igl3.name}${partnerBand ? `, ${partnerBand}` : ''})`
        : `chr${three.chr}${partnerBand ? ` (${partnerBand})` : ''}, where there is no MANE Select gene`;
    // With the locus table loaded the partner IS identified, so the go-and-check caveat would be
    // noise; it survives only when the tool genuinely cannot say.
    if (bp === 'none' && !igl3) caveats.push(`Immunoglobulin and T-cell-receptor loci are absent from MANE Select, so they always show as "no MANE gene" here even when they are the clinically important side.${partnerBand ? ` This partner is at ${partnerBand} — worth checking against the immunoglobulin heavy chain (14q32), kappa (2p11), lambda (22q11) and the receptor loci (7q34, 14q11) before dismissing it.` : ''}`);
    if (ap !== '5prime') {
      return { kind: 'regulatory-into-nothing', productive: false,
        headline: `5′ regulatory region of ${name5} → ${into}`,
        summary: 'No product. Regulatory sequence only; no gene body follows it.',
        caveats };
    }
    // The gene keeps its own promoter and reads across the junction. Whether its PROTEIN survives is
    // an exon question, and with the exon table loaded it is answered rather than guessed: every
    // coding base on the kept side means the break is 3′ of the stop codon, so the protein is whole.
    //
    // The ≥90%-of-the-genomic-span rule below is what the tool used before it had exon coordinates.
    // It is a proxy for the same question and it disagrees with the answer often — a gene whose
    // coding sequence ends early in a long span reads as "truncated" on span while every coding base
    // survives, and a gene with a big last coding exon reads as "whole" on span while a third of its
    // protein is gone. So the proxy is now the FALLBACK only, used when the exon table has nothing to
    // say about this gene, and it still declares itself a proxy when it is the thing deciding.
    const spanNearlyWhole = a.keptFraction != null && a.keptFraction >= 0.9;
    const codingDecidable = donorCodingDecidable(ax);
    const donorWhole = codingDecidable ? donorCodingIntact(ax) : spanNearlyWhole;
    if (!ax) caveats.push(`This tool has no exon coordinates loaded for ${name5}, so it reads "near the 3′ end" from the genomic span: it says likely-3′-end at 90% or more retained. That is a proxy, not an exon-level call.`);
    else if (!codingDecidable) caveats.push(`The exon table holds no coding sequence for ${name5}, so "near the 3′ end" is read from the genomic span at 90% or more retained. That is a proxy, not an exon-level call.`);
    // Name the partner the same way the headline does. Saying "no annotated gene" while the headline
    // names an immunoglobulin locus contradicts the tool's own sentence above.
    const intoShort = bp === 'downstream' ? `sequence 3′ of ${name3}`
      : igl3 ? `the ${igl3.symbol} locus`
        : 'sequence with no annotated gene';
    // What is computed here is that no gene body survives on the partner piece, so there is nothing
    // for the donor to fuse with. Enhancer activity is NOT computed — this tool holds no enhancer
    // data at all — so it is named as the recognised mechanism for an immunoglobulin or receptor
    // locus and left open for a stretch of chromosome that is merely unannotated.
    const partnerNote = igl3
      ? `Partner: the ${igl3.symbol} locus — no gene there to fuse with; deregulation by its enhancers is the recognised mechanism.`
      : `No gene on the partner side to fuse with; what that sequence does to ${name5} is not computed here.`;
    return {
      kind: donorWhole ? 'deregulation-candidate' : 'truncated-donor',
      productive: 'possible',
      headline: donorWhole
        ? `${name5} keeps its own promoter and ${codingDecidable ? 'its whole coding sequence' : keptPhrase(a)} → ${into}`
        : `${name5} 5′ part (${codingDecidable ? `${ax.codingKept.toLocaleString()} of ${ax.codingTotal.toLocaleString()} coding bases` : keptPhrase(a)}) → ${into}`,
      summary: donorWhole
        ? [
          codingDecidable
            ? `${name5} protein unchanged — all ${ax.codingTotal.toLocaleString()} coding bases intact; cut ${ax.label}, past the stop codon.`
            : `${name5} cut near its 3′ end, keeping ${keptPhrase(a)}.`,
          'Own promoter retained; 3′ end replaced.',
          partnerNote,
          `Transcript order: ${name5} first, partner second.`,
        ].join(' ')
        : [
          `${name5} truncated: ${codingDecidable ? `cut ${ax.label}, ${ax.codingKept.toLocaleString()} of ${ax.codingTotal.toLocaleString()} coding bases kept` : `${keptPhrase(a)} kept`}.`,
          `Own promoter retained; reads into ${intoShort}.`,
          `Likely consequence: disruption of ${name5}.`,
        ].join(' '),
      caveats,
    };
  }

  if (ap === 'intact') {
    // An intact gene with an immunoglobulin or T-cell-receptor locus brought in 3′ of it is the
    // VARIANT translocation pattern — t(2;8) and t(8;22) in Burkitt lymphoma. The enhancers act at a
    // distance, so no read-through is needed for the junction to matter.
    const igl3i = b.igLocus || null;
    if (bp === 'none' && igl3i) {
      return { kind: 'enhancer-adoption', productive: 'possible',
        headline: `intact ${name5} · ${igl3i.symbol} locus enhancers brought 3′`,
        summary: `${name5} complete${ax && ax.codingIntact ? ` — all ${ax.codingTotal.toLocaleString()} coding bases` : ''}; the junction lies past its 3′ end. The ${igl3i.symbol} locus (${igl3i.name}) now sits downstream. Enhancer adoption without read-through — the variant immunoglobulin-translocation pattern.`,
        caveats };
    }
    // Mirror of the unnamed-donor branches further down. With the locus table missing, an unnamed
    // partner 3′ of an intact gene is exactly the geometry the table would have called enhancer
    // adoption, so the lead has to be offered rather than the junction dismissed in silence. This
    // side used to be the only immunoglobulin-shaped branch that said nothing at all.
    const partnerBand3 = bands[three.chr] ? bands[three.chr].label : null;
    if (bp === 'none' && !igl3i) {
      caveats.push(`Immunoglobulin and T-cell-receptor loci are absent from MANE Select, so this unnamed partner 3′ of ${name5} may well be one.${partnerBand3 ? ` It is at ${partnerBand3} — check that against 14q32 (heavy chain), 2p11 (kappa), 22q11 (lambda) and the receptor loci at 7q34 and 14q11.` : ' Check its cytoband against 14q32, 2p11, 22q11, 7q34 and 14q11.'}`);
      caveats.push(`If it is one, an intact ${name5} with those enhancers 3′ of it is the variant immunoglobulin-translocation pattern — the t(2;8) and t(8;22) geometry. Load the locus table to settle it; treat this as a lead, not a conclusion.`);
    }
    const shared = `The junction lies 3′ of ${name5}, past its end, so ${name5} itself is complete and reads out normally into its own terminator — nothing transcribes through the junction from this side.`;
    if (bp === 'intact') {
      return { kind: 'no-disruption', productive: false,
        headline: `intact ${name5} · intact ${name3}`,
        summary: `${shared} ${name3} is also complete. Neither gene is broken; the junction rearranges the sequence between them.`,
        caveats };
    }
    if (bp === '3prime') {
      return { kind: 'disrupted-acceptor', productive: false,
        headline: `intact ${name5} → ${name3} 3′ part (${name3} disrupted)`,
        summary: `${shared} ${name3} has lost its 5′ end and promoter at the breakpoint, so this derivative disrupts ${name3} without making a fusion transcript.`,
        caveats };
    }
    return { kind: 'no-fusion-transcript', productive: false,
      headline: `intact ${name5} → ${bp === 'none' ? `chr${three.chr} (no MANE gene)` : b.label}`,
      summary: shared, caveats };
  }

  if (ap === 'downstream') {
    return { kind: 'no-promoter', productive: false,
      headline: `sequence 3′ of ${name5} → ${bp === 'none' ? `chr${three.chr} (no MANE gene)` : b.label}`,
      summary: `The 5′ piece is flanking sequence downstream of ${name5} with no promoter, so nothing transcribes into the junction.`,
      caveats };
  }

  // 5′ side has no named gene
  const donorBand = bands[five.chr] ? bands[five.chr].label : null;
  const igl5 = a.igLocus || null;
  if (bp === 'intact') {
    if (!igl5) {
      caveats.push(`Immunoglobulin and T-cell-receptor loci are absent from MANE Select, so this unnamed partner may well be one.${donorBand ? ` It is at ${donorBand} — check that against 14q32 (heavy chain), 2p11 (kappa), 22q11 (lambda) and the receptor loci at 7q34 and 14q11.` : ' Check its cytoband against 14q32, 2p11, 22q11, 7q34 and 14q11.'}`);
      caveats.push('Without a named donor the tool cannot say what drives the gene; treat this as a lead, not a conclusion.');
    }
    return { kind: 'unnamed-donor-intact', productive: 'possible',
      headline: igl5
        ? `${igl5.symbol} locus → intact ${name3}`
        : `chr${five.chr}${donorBand ? ` (${donorBand})` : ''} — no MANE gene → intact ${name3}`,
      summary: igl5
        ? `${name3} intact${bx && bx.codingIntact ? ` — ${bx.codingTotal.toLocaleString()} coding bases across ${bx.exons} exons` : ''}; upstream is now the ${igl5.symbol} locus (${igl5.name}). Deregulation by its enhancers is the recognised mechanism.`
        : `${name3} intact${bx && bx.codingIntact ? ` — ${bx.codingTotal.toLocaleString()} coding bases across ${bx.exons} exons` : ''}; new upstream region acquired from chr${five.chr}. Deregulation if that region carries a promoter or enhancer, as an immunoglobulin locus does.`,
      caveats };
  }
  if (bp === '3prime') {
    // Same rule as headless-acceptor above: the span label says the breakpoint is inside the gene,
    // not that the reading frame is gone. An acceptor cut 5′ of its start codon keeps its whole
    // protein, and an immunoglobulin locus upstream then supplies the promoter and enhancers instead
    // of "cannot restore them". This is t(3;14) IGH-BCL6 with the usual intron-1 break.
    if (acceptorCodingIntact(bx)) {
      if (!igl5) {
        caveats.push(`Immunoglobulin and T-cell-receptor loci are absent from MANE Select, so this unnamed partner may well be one.${donorBand ? ` It is at ${donorBand} — check that against 14q32 (heavy chain), 2p11 (kappa), 22q11 (lambda) and the receptor loci at 7q34 and 14q11.` : ' Check its cytoband against 14q32, 2p11, 22q11, 7q34 and 14q11.'}`);
        caveats.push('Without a named donor the tool cannot say what drives the gene; treat this as a lead, not a conclusion.');
      }
      return { kind: 'promoter-substitution-unnamed-donor', productive: 'possible',
        headline: `${igl5 ? `${igl5.symbol} locus` : `chr${five.chr}${donorBand ? ` (${donorBand})` : ''} — no MANE gene`} → ${name3}, coding sequence complete`,
        summary: [
          `${name3} protein unchanged — all ${bx.codingTotal.toLocaleString()} coding bases across ${bx.exons} exons kept.`,
          `Cut ${bx.label}, 5′ of the start codon: only its own promoter and 5′ untranslated region are lost.`,
          igl5
            ? `The ${igl5.symbol} locus now supplies both — deregulation of an unchanged protein, no fusion protein.`
            : `Deregulation if the new upstream region on chr${five.chr} carries a promoter or enhancer, as an immunoglobulin locus does.`,
        ].join(' '),
        caveats };
    }
    return { kind: 'unnamed-donor-body', productive: false,
      headline: `${igl5 ? `${igl5.symbol} locus` : `chr${five.chr} (no MANE gene)`} → ${name3} 3′ part`,
      summary: `No product. ${name3} has lost its 5′ end and promoter${bx ? `, cut ${bx.label} — ${bx.codingLost.toLocaleString()} of ${bx.codingTotal.toLocaleString()} coding bases lost, including the start codon` : ''}; ${igl5 ? `the ${igl5.symbol} locus upstream cannot restore them` : 'no annotated gene upstream to drive it'}.`,
      caveats };
  }
  // No gene BODY survives on either piece — that is what this branch means, and it is what makes the
  // junction non-productive. A piece can still NAME a gene: its flank rode across while the gene
  // itself stayed on the discarded side. Saying "no MANE gene within 2 Mb" for such a piece
  // contradicts the tool's own segment content, which named the gene 3 kb away.
  const flankOf = (seg) => {
    const c = seg.content;
    if (!c.gene || (c.part !== 'upstream' && c.part !== 'downstream')) return null;
    return { gene: c.gene, side: c.part === 'upstream' ? '5′' : '3′', bp: c.distance || 0 };
  };
  const howFar = (seg) => {
    const f = flankOf(seg);
    if (f) return `chr${seg.chr} — ${f.bp.toLocaleString()} bp ${f.side} of ${f.gene}, which is not on this piece`;
    const igp = seg.content.igLocus;
    if (igp) return `chr${seg.chr} — inside the ${igp.symbol} locus`;
    const n = seg.content.nearestGene;
    return n ? `chr${seg.chr} — nearest MANE gene is ${n.gene}, ${Math.round(n.distance / 1000).toLocaleString()} kb away`
      : `chr${seg.chr} — no MANE gene within 2 Mb`;
  };
  const sentences = ['No gene body on either piece — nothing to call a fusion.'];
  for (const seg of [five, three]) {
    const f = flankOf(seg);
    if (f) sentences.push(`The chr${seg.chr} side carries only sequence ${f.bp.toLocaleString()} bp ${f.side} of ${f.gene}; the gene itself went the other way.`);
    else if (seg.content.igLocus) sentences.push(`The chr${seg.chr} side lies inside the ${seg.content.igLocus.symbol} locus (${seg.content.igLocus.name}).`);
  }
  if (sentences.length === 1) {
    sentences.push('No MANE gene within 50 kb of either breakpoint.');
    sentences.push('An immunoglobulin or receptor locus looks identical here — use the distances above to tell them apart.');
  }
  return { kind: 'no-gene', productive: false,
    headline: `${howFar(five)}; ${howFar(three)}`,
    summary: sentences.join(' '),
    caveats };
}

/**
 * The full interpretation of one junction: its two pieces read in transcript 5′→3′ order, which
 * derivative chromosome it builds, and what — if anything — the read-through produces.
 */
export function orientDerivative(genes, b, opts = {}) {
  const grammarOrder = derivativeSegments(b).map((s) => ({ ...s, content: segmentContent(genes, s, opts) }))
    .map((s) => ({ ...s, exon: exonFactsFor(s, opts.exons) }));
  const dirOf = (s) => (s.content.gene ? SIGN[s.orient] * SIGN[s.content.strand] : null);
  const d0 = dirOf(grammarOrder[0]); const d1 = dirOf(grammarOrder[1]);

  let coherence; let flipped = false;
  if (d0 !== null && d1 !== null) {
    if (d0 === d1) { coherence = 'coherent'; flipped = d0 === -1; }
    // The genes disagree: pointing towards each other across the junction (tail to tail) or away
    // from it (head to head). Either way nothing can transcribe through. Both labels survive a flip
    // of the whole molecule, so the grammar order is enough to tell them apart.
    else coherence = d0 === 1 ? 'convergent' : 'divergent';
  } else if (d0 !== null || d1 !== null) {
    coherence = 'coherent';
    flipped = (d0 !== null ? d0 : d1) === -1;
  } else {
    coherence = 'no-gene';
  }

  const reading = flipped ? flipReading(grammarOrder) : grammarOrder;
  // "Promoters now back to back" used to be printed for every divergent junction. It describes
  // sequence that is often not there: with both genes read away from the junction, the piece keeps
  // the 3′ portion of each and both promoters were discarded. Which 5′ ends are actually on this
  // derivative is computable from the part label, so say that instead of assuming.
  const carriesFivePrimeEnd = (s) => ['intact', '5prime', 'upstream'].includes(s.content.part);
  const withPromoter = reading.filter(carriesFivePrimeEnd).map((s) => s.content.gene);
  const divergentTail = withPromoter.length === 2
    ? 'Both promoters are on this derivative, now back to back — expression may change, but no fusion.'
    : withPromoter.length === 1
      ? `Only ${withPromoter[0]}'s promoter is on this derivative; the other gene's 5′ end went the other way.`
      : 'Neither promoter is on this derivative — both 5′ ends went the other way.';
  // With no gene on either side there is no reading direction to get wrong, so the classifier still
  // runs — it is what reports how far away the nearest gene is. Only a genuine direction clash
  // (the genes point at or away from each other) short-circuits it.
  const product = coherence === 'convergent' || coherence === 'divergent'
    ? {
      kind: 'no-read-through',
      productive: false,
      headline: `${reading[0].content.gene} and ${reading[1].content.gene} point ${coherence === 'convergent' ? 'towards each other' : 'away from each other'} across the junction`,
      summary: coherence === 'convergent'
        ? 'No transcript. Genes transcribed towards the junction, tail to tail. Each affected only by what it lost at its own breakpoint.'
        : `No transcript. Genes transcribed away from the junction, head to head. ${divergentTail}`,
      caveats: [],
    }
    : classifyProduct(reading[0], reading[1], { bands: b.bands, exons: opts.exons });

  return {
    notation: b.raw,
    mateRecords: b.mateRecords || [b.raw],
    support: b.support || null,
    germline: b.germline || null,
    insertion: cleanInsertion(b.bases, b.tPosition),
    arms: b.arms || {},
    bands: b.bands || {},
    derivative: nameDerivative(reading, b.arms),
    reading,                 // the two pieces in transcript 5′→3′ order
    grammarOrder,            // the same pieces as the VCF record writes them
    flipped,                 // true when the record's own order is the reverse of the reading order
    reverseComplement: b.strand === 'reverse',
    coherence,               // coherent | convergent | divergent | no-gene
    product,
  };
}

/** Interpret every junction of an event, productive derivative first. */
export function describeDerivatives(event, genes = {}, opts = {}) {
  const out = event.breakends.map((b, i) => {
    const d = orientDerivative(genes, b, opts);
    return { ...d, index: i, label: d.derivative.name || `Junction ${i + 1}` };
  });
  const rank = (d) => (d.product.productive === true ? 0 : (d.product.productive === 'possible' ? 1 : 2));
  return out.slice().sort((p, q) => rank(p) - rank(q) || p.index - q.index);
}
