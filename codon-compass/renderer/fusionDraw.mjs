// The Fusions tab's figures, as SVG strings: the chromosome ideograms, the derivative molecule and
// the reference loci. Ported from the reviewed design (scripts/build_fusion_versions.mjs) so the app
// shows exactly what was signed off.
//
// Pure string builders — no DOM access — so the same code serves the desktop renderer, the browser
// builds and any generator script. Every function takes a `ctx` carrying the reference tables:
//   { genes, exons, cytobands, igLoci, igGenes }
// and every drawn element derives from those tables or the parsed row; nothing is invented here.

import { partShort } from '../core/fusion.mjs';

const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const num = (x) => Number(x).toLocaleString('en-GB');

// Giemsa stain shades: light body, dark positive bands — how a karyotype is read.
export const STAIN = {
  gneg: '#eef2f9', gpos25: '#c2cadb', gpos33: '#b1bacd', gpos50: '#939db6',
  gpos66: '#7e879f', gpos75: '#6d7791', gpos100: '#4a536b', gpos: '#6d7791',
  // The centromere is read against a pale chromosome body, so it needs to be light enough to show as a
  // constriction rather than a dark smudge.
  acen: '#e08a8f', gvar: '#a8b0c6', stalk: '#8f97ad',
};

/**
 * The immunoglobulin or T-cell-receptor locus a breakpoint falls in.
 *
 * Loci nest: the T-cell receptor delta locus (TRD, 44 kb) sits wholly inside the alpha locus (TRA,
 * 930 kb). The SMALLEST containing locus wins, which is the rule `segmentContent` in core/fusion.mjs
 * already uses for the wording. Taking whichever locus the table happened to list first named a TRD
 * breakpoint "TRA locus" in the figures while the text beside them said "TRD locus", and drew the
 * fragment box round the whole 930 kb instead of the 44 kb the event is in.
 */
export function igLocusAt(ctx, chr, pos) {
  const loci = (ctx.igLoci && ctx.igLoci.loci) || [];
  return loci
    .filter((l) => l.chr === String(chr) && pos >= l.start && pos <= l.end)
    .sort((p, q) => (p.end - p.start) - (q.end - q.start))[0] || null;
}

// The two panel colours. A chromosome keeps its colour across every figure of one event card.
export const CHROMOSOME_COLOURS = ['#c8823a', '#3f7fbf'];

/**
 * Display order for chromosome names: numeric, then X, then Y, with any other name last and
 * alphabetical among its equals. This is the order the Fusions tab uses everywhere, so t(3;12) is
 * never shown as 12/3. Kept identical to `fusionChrOrder` in renderer/app.mjs.
 */
export function chromosomeOrder(a, b) {
  if (a === b) return 0;
  const rank = (c) => (c === 'X' ? 23 : (c === 'Y' ? 24 : (Number(c) || 99)));
  return rank(a) - rank(b) || (a < b ? -1 : 1);
}

/**
 * One colour per chromosome: the first in display order takes the first palette colour, everything
 * else the second.
 *
 * Colour has to mean the same thing in every figure of one event card, because that is what lets a
 * reader carry a chromosome from the miniature circos down into the panels. Assigning by READING
 * order instead made a reciprocal pair contradict itself: the two derivatives read their pieces in
 * opposite orders, so chr3 was orange in the circos and in the der(12) panel but blue in the der(3)
 * panel directly beneath it.
 *
 * On an intra-chromosomal event — an inversion, say — both pieces come from the same chromosome and
 * so share its colour. The junction mark and the fragment labels separate them; colour cannot, since
 * there is only one chromosome to name.
 */
export function chromosomeColours(chrs) {
  const out = {};
  [...new Set([...chrs].map(String))].sort(chromosomeOrder)
    .forEach((c, i) => { out[c] = CHROMOSOME_COLOURS[i === 0 ? 0 : 1]; });
  return out;
}

// Every drawn ideogram needs its own clip identifier. Identifiers are resolved across the whole
// page, not within one figure, and two ideograms on a page can easily share a chromosome and a first
// breakpoint — several pasted rows all breaking at 14:105,864,327, say. An identifier built from
// those alone repeats, and the second figure then clips to the first figure's rectangle, which is
// the wrong width whenever the two were drawn at different sizes. A counter cannot repeat.
let clipSeq = 0;

/** One chromosome to scale: real bands, real centromere, breakpoints ticked, band label above them. */
export function ideogram(ctx, chr, positions, { width = 1080, height = 72, label = true } = {}) {
  const c = ctx.cytobands && ctx.cytobands.chromosomes && ctx.cytobands.chromosomes[chr];
  if (!c) return `<svg class="ideo" viewBox="0 0 ${width} ${height}"><text x="8" y="30" class="ideo-sub">chr${esc(chr)} — no banding table</text></svg>`;
  const X0 = label ? 108 : 8;
  const W = width - X0 - 12; const H = 17; const Y = 26; const rr = H / 2;
  const at = (bp) => X0 + (Math.max(0, Math.min(bp, c.length)) / c.length) * W;

  const parts = [];
  for (const b of c.bands) {
    if (b.stain === 'acen') continue;
    const x = at(b.start); const w = Math.max(0.6, at(b.end) - x);
    parts.push(`<rect x="${x.toFixed(2)}" y="${Y}" width="${w.toFixed(2)}" height="${H}" fill="${STAIN[b.stain] || STAIN.gneg}"/>`);
  }
  // The centromere is a CONSTRICTION along the chromosome, so the two triangles meet point to point
  // across the bar: each has its flat base against the arm it belongs to and narrows to the waist in
  // the middle. Drawn the other way up — bases on the top and bottom edges — it reads as a cross
  // sitting on the chromosome rather than a pinch in it.
  const cenX = at(c.centromere); const cw = 7;
  parts.push(`<path d="M${cenX - cw} ${Y} L${cenX} ${Y + H / 2} L${cenX - cw} ${Y + H} Z`
    + ` M${cenX + cw} ${Y} L${cenX} ${Y + H / 2} L${cenX + cw} ${Y + H} Z" fill="${STAIN.acen}"/>`);

  const uniq = [...new Set(positions)].sort((a, b) => a - b);
  const ticks = uniq.map((p) => {
    const x = at(p);
    return `<line x1="${x.toFixed(2)}" y1="${Y - 7}" x2="${x.toFixed(2)}" y2="${Y + H + 7}" stroke="#ff5f57" stroke-width="2"/>`
      + `<polygon points="${(x - 4).toFixed(2)},${Y - 12} ${(x + 4).toFixed(2)},${Y - 12} ${x.toFixed(2)},${Y - 5}" fill="#ff5f57"/>`;
  }).join('');

  const bandAt = (p) => { const b = c.bands.find((x) => p >= x.start && p < x.end); return b ? b.band : null; };
  const called = [...new Set(uniq.map((p) => `${chr}${bandAt(p) || ''}`))].join(', ');
  const xs = uniq.map((p) => at(p));
  const cx = xs.reduce((q, r) => q + r, 0) / (xs.length || 1);
  const anchor = cx < X0 + 50 ? 'start' : (cx > X0 + W - 50 ? 'end' : 'middle');

  clipSeq += 1;
  const clipId = `clip-${esc(chr)}-${uniq[0]}-${clipSeq}`;
  return `<svg class="ideo" viewBox="0 0 ${width} ${height}" role="img" aria-label="Chromosome ${esc(chr)}, breakpoint in ${esc(called)}">
  ${label ? `<text x="8" y="${Y + 13}" class="ideo-name">chr${esc(chr)}</text>
  <text x="8" y="${Y + 27}" class="ideo-sub">${(c.length / 1e6).toFixed(0)} Mb</text>` : ''}
  <clipPath id="${clipId}"><rect x="${X0}" y="${Y}" width="${W}" height="${H}" rx="${rr}"/></clipPath>
  <g clip-path="url(#${clipId})">${parts.join('')}</g>
  <rect x="${X0}" y="${Y}" width="${W}" height="${H}" rx="${rr}" fill="none" stroke="#0d1220" stroke-width="1.2"/>
  ${ticks}
  <text x="${X0 - 4}" y="${Y + H + 16}" text-anchor="end" class="ideo-arm">p</text>
  <text x="${X0 + W + 4}" y="${Y + H + 16}" class="ideo-arm">q</text>
  <text x="${cx.toFixed(2)}" y="${Y - 15}" text-anchor="${anchor}" class="ideo-band">${esc(called)}</text>
</svg>`;
}

/** A pentagon pointing the way the gene is transcribed. */
function geneArrow(x, w, y, h, colour, pointsRight, opacity = 1) {
  const tip = Math.min(16, Math.max(6, w * 0.18));
  const pts = pointsRight
    ? `${x},${y} ${x + w - tip},${y} ${x + w},${y + h / 2} ${x + w - tip},${y + h} ${x},${y + h}`
    : `${x + w},${y} ${x + tip},${y} ${x},${y + h / 2} ${x + tip},${y + h} ${x + w},${y + h}`;
  return `<polygon points="${pts}" fill="${colour}" fill-opacity="${opacity}" stroke="${colour}" stroke-width="1"/>`;
}

function breakMark(x, yTop, yBottom, label) {
  return `<line x1="${x.toFixed(2)}" y1="${yTop}" x2="${x.toFixed(2)}" y2="${yBottom}" stroke="#d9534f" stroke-width="1.6" stroke-dasharray="4 3"/>`
    + `<circle cx="${x.toFixed(2)}" cy="${(yTop + yBottom) / 2}" r="3.5" fill="#d9534f"/>`
    + (label ? `<text x="${x.toFixed(2)}" y="${yBottom + 14}" text-anchor="middle" class="rx-bp">${label}</text>` : '');
}

/** The V/D/J/constant segment genes inside an IG or TCR locus, as ticks in the band. */
function igSegments(ctx, out, symbol, at, yMid, lo, hi) {
  const seg = (ctx.igGenes && ctx.igGenes.loci && ctx.igGenes.loci[symbol]) || [];
  for (const sg of seg) {
    if (sg.end < lo || sg.start > hi) continue;
    const sa = at(Math.max(sg.start, lo)); const sb = at(Math.min(sg.end, hi));
    const x = Math.min(sa, sb); const w = Math.max(1, Math.abs(sb - sa));
    const h = sg.cls === 'C' ? 18 : (sg.cls === 'V' ? 9 : 12);
    out.push(`<rect x="${x.toFixed(2)}" y="${(yMid - h / 2).toFixed(1)}" width="${w.toFixed(2)}" height="${h}" fill="#9aa6bc" fill-opacity="${sg.cls === 'C' ? 0.95 : sg.cls === 'V' ? 0.4 : 0.65}"><title>${esc(sg.symbol)}</title></rect>`);
  }
}

/**
 * A gene with its structure on top: body, dark intron line, dark exon blocks, and — where it is true —
 * a 5′ mark at the transcription start.
 *
 * `shape` is 'arrow' or 'block'. The reference figure draws each gene as an arrow, because there it is
 * showing the genome as it is and each gene's own direction is the point. The derivative draws plain
 * blocks that butt flush at the junction: the whole molecule is read left to right and says so once,
 * at the top, so a per-gene arrowhead in the middle of the fused product is noise.
 *
 * `fiveMark` is false on the derivative. A gene cut at its 5′ end has LOST its transcription start, so
 * drawing a 5′ mark at the breakpoint would put a label on sequence the derivative does not carry.
 */
function geneStructure(out, gEx, at, yMid, colour, lo, hi, { shape = 'arrow', fiveMark = true } = {}) {
  const flat = gEx.exons;
  const gs = flat[0]; const ge = flat[flat.length - 1];
  const vs = Math.max(gs, lo); const ve = Math.min(ge, hi);
  if (ve <= vs) return;
  const xa = at(vs); const xb = at(ve);
  const xl = Math.min(xa, xb); const xr = Math.max(xa, xb);
  const p5 = gEx.strand === '+' ? gs : ge;
  const p3 = gEx.strand === '+' ? ge : gs;
  const sd = Math.sign(at(p3) - at(p5)) || 1;
  const bodyW = Math.max(10, xr - xl);
  // A gene arrow's TIP is its 3′ end. On the derivative that end is either the breakpoint — where the
  // gene was cut, so the piece must finish flush against the junction — or the gene's own end, which
  // is still there and still points the way it is transcribed. So the shape follows the biology: keep
  // the point where the 3′ end survives, square it off where the breakpoint took it.
  const threePrimeEndKept = p3 >= lo && p3 <= hi;
  out.push(shape === 'arrow' || threePrimeEndKept
    ? geneArrow(xl, bodyW, yMid - 13, 26, colour, sd > 0, 0.92)
    : `<rect x="${xl.toFixed(2)}" y="${yMid - 13}" width="${bodyW.toFixed(2)}" height="26" fill="${colour}" fill-opacity="0.92" stroke="${colour}" stroke-width="1"/>`);
  out.push(`<line x1="${xl.toFixed(2)}" y1="${yMid}" x2="${xr.toFixed(2)}" y2="${yMid}" stroke="#0f1420" stroke-width="1.4" stroke-opacity="0.5"/>`);
  const cds = gEx.cds;
  for (let k = 0; k + 1 < flat.length; k += 2) {
    const es = Math.max(flat[k], lo); const ee = Math.min(flat[k + 1], hi);
    if (ee < es) continue;
    const segs = [];
    if (!cds) segs.push([es, ee, false]);
    else {
      if (es < cds[0]) segs.push([es, Math.min(ee, cds[0] - 1), false]);
      const cs = Math.max(es, cds[0]); const ce = Math.min(ee, cds[1]);
      if (cs <= ce) segs.push([cs, ce, true]);
      if (ee > cds[1]) segs.push([Math.max(es, cds[1] + 1), ee, false]);
    }
    for (const [ps, pe, coding] of segs) {
      const pa = at(ps); const pb = at(pe);
      const ex = Math.min(pa, pb); const ew = Math.max(2.2, Math.abs(pb - pa));
      const h = coding ? 22 : 11;
      out.push(`<rect x="${ex.toFixed(2)}" y="${(yMid - h / 2).toFixed(1)}" width="${ew.toFixed(2)}" height="${h}" fill="#0f1420" fill-opacity="${coding ? 0.7 : 0.4}"/>`);
    }
  }
  // Only mark the transcription start when it is actually on the drawn piece. Clamping it into the
  // window puts the mark on a breakpoint, which says the gene starts where it was in fact cut.
  if (fiveMark && p5 >= lo && p5 <= hi) {
    const cx5 = at(p5);
    out.push(`<text x="${(cx5 - sd * 8).toFixed(2)}" y="${yMid + 4}" text-anchor="${sd > 0 ? 'end' : 'start'}" class="rx-five">5′</text>`);
  }
}

// How much flanking sequence to draw for a fragment with no gene and no immunoglobulin locus at its
// breakpoint. There is nothing on the fragment to set a natural window, so a fixed one is used. It is
// a viewing window only — the fragment itself runs on to the telomere — and it is laid entirely on
// the RETAINED side of the breakpoint, never across it.
const NO_GENE_WINDOW = 30000;

/**
 * The drawn window for each fragment of a derivative, and where a genomic position lands inside it.
 *
 * Split out from the drawing so the controls can check the geometry that carries the meaning: which
 * side of the junction a piece sits on, whether it is drawn mirrored, and — the invariant that makes
 * the figure honest — that every window ends ON the breakpoint, so only retained sequence is drawn
 * and the breakpoint always lands on the junction edge.
 *
 * Returns { width, pad, junctionX, toScale, sides: [{ seg, g, ig, keepLeft, lo, hi, bp, mirrored,
 * x0, w, at, junctionEdge }] }, where `at(bp)` is the genomic-to-screen mapping for that fragment.
 */
export function derivativeLayout(ctx, d, width = 1080) {
  // The margin carries the 5′ and 3′ marks for the whole molecule, so it is wider than a plain inset.
  const pad = 30;
  const EX = ctx.exons && ctx.exons.genes ? ctx.exons.genes : {};
  // Keyed by chromosome, not by position in the reading, so the same chromosome keeps its colour in
  // the circos above and in the other derivative of a reciprocal pair.
  const COL = chromosomeColours(d.reading.map((s) => s.chr));

  const sides = d.reading.map((seg) => {
    const g = seg.content.gene ? EX[seg.content.gene] : null;
    const ig = !g ? igLocusAt(ctx, seg.chr, seg.pos) : null;
    // `left` means the sequence below the breakpoint is the piece kept on this derivative, so the
    // window must run up to the breakpoint and stop; `right` means it must start there.
    const keepLeft = seg.side === 'left';
    let lo; let hi;
    if (g) {
      const gs = g.exons[0]; const ge = g.exons[g.exons.length - 1];
      lo = keepLeft ? Math.min(gs, seg.pos) : seg.pos;
      hi = keepLeft ? seg.pos : Math.max(ge, seg.pos);
      // a promoter-side fragment holds none of the gene body; show its flanking DNA, not a sliver
      if (hi - lo < 2000) { if (keepLeft) lo = seg.pos - 12000; else hi = seg.pos + 12000; }
    } else if (ig) {
      lo = keepLeft ? ig.start : seg.pos;
      hi = keepLeft ? seg.pos : ig.end;
    } else {
      lo = keepLeft ? seg.pos - NO_GENE_WINDOW : seg.pos;
      hi = keepLeft ? seg.pos : seg.pos + NO_GENE_WINDOW;
    }
    return { seg, g, ig, keepLeft, lo, hi, bp: Math.max(1, hi - lo), mirrored: seg.orient === '-', col: COL[String(seg.chr)] };
  });

  // No gap between the fragments: a junction joins the end of one piece to the start of the next, so
  // the two are drawn flush and the junction mark sits on the seam they share.
  const avail = width - pad * 2;
  const total = sides.reduce((a, s) => a + s.bp, 0);
  let w = sides.map((s) => (s.bp / total) * avail);
  const MINW = 170; let toScale = true;
  if (w.some((v) => v < MINW)) {
    toScale = false;
    const shortIdx = w.map((v, i) => (v < MINW ? i : -1)).filter((i) => i >= 0);
    const rest = avail - shortIdx.length * MINW;
    const restBp = sides.filter((_, i) => !shortIdx.includes(i)).reduce((a, s) => a + s.bp, 0) || 1;
    w = w.map((v, i) => (shortIdx.includes(i) ? MINW : (sides[i].bp / restBp) * rest));
  }

  let x = pad; let junctionX = null;
  const laid = sides.map((s, i) => {
    const x0 = x; const wi = w[i];
    const at = (bp) => {
      const f = (Math.max(s.lo, Math.min(bp, s.hi)) - s.lo) / s.bp;
      return x0 + (s.mirrored ? 1 - f : f) * wi;
    };
    // The first piece meets the junction at its right-hand edge, the second at its left-hand edge.
    const out = { ...s, x0, w: wi, at, junctionEdge: i === 0 ? x0 + wi : x0 };
    x += wi;
    if (i === 0) junctionX = x;
    return out;
  });

  return { width, pad, junctionX, toScale, sides: laid };
}

/** The derivative: both retained fragments joined at the junction, reading 5′ → 3′ left to right. */
export function derivativeMolecule(ctx, d, width = 1080) {
  const H = 142; const y = 68;
  const L = derivativeLayout(ctx, d, width);
  const pad = L.pad; const jx = L.junctionX; const toScale = L.toScale;

  let out = '';
  L.sides.forEach((s) => {
    const mirrored = s.mirrored;
    const x0 = s.x0; const wi = s.w; const at = s.at;
    out += `<line x1="${x0.toFixed(2)}" y1="${y}" x2="${(x0 + wi).toFixed(2)}" y2="${y}" stroke="${s.col}" stroke-width="2" stroke-opacity="0.45"/>`;
    const parts = [];
    let name; let sub;
    if (s.g) {
      geneStructure(parts, s.g, at, y, s.col, s.lo, s.hi, { shape: 'block', fiveMark: false });
      name = s.seg.content.gene;
      const bodyIn = !(s.g.exons[s.g.exons.length - 1] < s.lo || s.g.exons[0] > s.hi);
      sub = bodyIn
        ? `${partShort(s.seg.content.part)} · ${(s.bp / 1000).toFixed(1)} kb${mirrored ? ' · reverse-complemented' : ''}`
        : `flanking DNA only${mirrored ? ' · reverse-complemented' : ''}`;
    } else if (s.ig) {
      const ba = at(s.lo); const bb = at(s.hi);
      const xl = Math.min(ba, bb); const xr = Math.max(ba, bb);
      parts.push(`<rect x="${xl.toFixed(2)}" y="${y - 11}" width="${(xr - xl).toFixed(2)}" height="22" fill="${s.col}" fill-opacity="0.2" stroke="${s.col}" stroke-dasharray="4 3"/>`);
      igSegments(ctx, parts, s.ig.symbol, at, y, s.lo, s.hi);
      name = `${s.ig.symbol} locus`;
      sub = `${(s.bp / 1000).toFixed(0)} kb of the ${s.ig.name}${mirrored ? ' · reverse-complemented' : ''}`;
    } else {
      name = `chr${s.seg.chr}`;
      const ng = s.seg.content.nearestGene;
      // Nothing on this fragment fixes its drawn extent, so say how much of it is shown: the piece
      // itself carries on past the window, away from the junction.
      const shown = `${(s.bp / 1000).toFixed(0)} kb shown, retained side only`;
      sub = `${ng ? `nearest gene ${ng.gene}, ${(ng.distance / 1000).toFixed(0)} kb` : 'no annotated gene'} · ${shown}${mirrored ? ' · reverse-complemented' : ''}`;
    }
    out += parts.join('');
    out += `<text x="${(x0 + wi / 2).toFixed(2)}" y="${y - 26}" text-anchor="middle" class="rx-gfocus" fill="${s.col}">${esc(name)}</text>`;
    out += `<text x="${(x0 + wi / 2).toFixed(2)}" y="${y + 32}" text-anchor="middle" class="rx-dsub">${esc(sub)}</text>`;
    out += `<text x="${(x0 + wi / 2).toFixed(2)}" y="${y + 45}" text-anchor="middle" class="rx-dpos">chr${esc(s.seg.chr)}:${num(s.seg.pos)}</text>`;
  });
  let out2 = out + breakMark(jx, y - 30, y + 18, null);
  // The molecule's own ends. Only two marks, and they are the only ones that are true: the fused
  // transcript starts at the left end of the drawing and finishes at the right.
  out2 += `<text x="${(pad - 6).toFixed(2)}" y="${y + 4}" text-anchor="end" class="rx-five">5′</text>`;
  out2 += `<text x="${(width - pad + 6).toFixed(2)}" y="${y + 4}" text-anchor="start" class="rx-five">3′</text>`;
  // The direction rail sits above everything, with its caption at the opposite end from its arrowhead
  // so neither can collide with a fragment name.
  out2 += `<text x="${pad}" y="20" class="rx-read">read 5′ → 3′</text>`;
  out2 += `<line x1="${pad + 84}" y1="16" x2="${width - pad}" y2="16" stroke="#5a6782" stroke-width="1"/>`;
  out2 += `<polygon points="${width - pad},16 ${width - pad - 8},13 ${width - pad - 8},19" fill="#5a6782"/>`;
  if (!toScale) out2 += `<text x="${width - pad}" y="${H - 5}" text-anchor="end" class="rx-note">fragments not to a common scale</text>`;
  return `<svg class="rx" viewBox="0 0 ${width} ${H}" role="img" aria-label="Derivative ${esc(d.derivative.name || '')}">${out2}</svg>`;
}

/** Reference loci: one row per chromosome — the focus gene as structure, neighbours grey, an IG
 *  locus as a populated band, breakpoints as dashed rules. */
export function referenceLoci(ctx, d, width = 1080) {
  const rowH = 92; const pad = 14; const labelW = 96;
  const chrs = [...new Set(d.reading.map((s) => s.chr))];
  const H = 30 + chrs.length * rowH;
  const GEN = ctx.genes || {};
  const EX = ctx.exons && ctx.exons.genes ? ctx.exons.genes : {};
  // Same colour rule as the derivative and the circos: by chromosome, in display order.
  const COL = chromosomeColours(chrs);

  let outAll = `<text x="${pad}" y="16" class="rx-head">Reference loci</text>`;
  chrs.forEach((chr, ri) => {
    const y = 48 + ri * rowH;
    const out = [];
    const segs = d.reading.filter((s) => s.chr === chr);
    const focus = segs.find((s) => s.content.gene) || segs[0];
    const fg = focus.content.gene ? EX[focus.content.gene] : null;
    const ig = !fg ? igLocusAt(ctx, chr, focus.pos) : null;
    const positions = segs.map((s) => s.pos);
    let lo0 = Math.min(...positions, fg ? fg.exons[0] : Infinity);
    let hi0 = Math.max(...positions, fg ? fg.exons[fg.exons.length - 1] : -Infinity);
    if (ig) { lo0 = Math.min(lo0, Math.max(ig.start, focus.pos - 500000)); hi0 = Math.max(hi0, Math.min(ig.end, focus.pos + 500000)); }
    const padBp = Math.max(4000, (hi0 - lo0) * 0.4);
    const lo = lo0 - padBp; const hi = hi0 + padBp; const span = Math.max(1, hi - lo);
    const x0 = pad + labelW; const inner = width - x0 - pad;
    const at = (bp) => x0 + ((bp - lo) / span) * inner;

    out.push(`<text x="${pad}" y="${y + 4}" class="rx-chr">chr${esc(chr)}</text>`);
    out.push(`<text x="${pad}" y="${y + 18}" class="rx-chrsub">${d.bands && d.bands[chr] ? esc(d.bands[chr].label) : ''}</text>`);
    out.push(`<line x1="${x0}" y1="${y}" x2="${x0 + inner}" y2="${y}" stroke="#4a5570" stroke-width="1"/>`);

    if (ig) {
      const ba = at(Math.max(ig.start, lo)); const bb = at(Math.min(ig.end, hi));
      out.push(`<rect x="${ba.toFixed(2)}" y="${y - 11}" width="${(bb - ba).toFixed(2)}" height="22" fill="#8a94a8" fill-opacity="0.22" stroke="#8a94a8" stroke-dasharray="4 3"/>`);
      out.push(`<text x="${((ba + bb) / 2).toFixed(2)}" y="${y - 18}" text-anchor="middle" class="rx-igname">${esc(ig.symbol)} — ${esc(ig.name)}</text>`);
      igSegments(ctx, out, ig.symbol, at, y, lo, hi);
      out.push(`<text x="${((ba + bb) / 2).toFixed(2)}" y="${y + 26}" text-anchor="middle" class="rx-bound">V segments thin · D/J mid · constant genes tall</text>`);
      if (ig.start >= lo) out.push(`<text x="${ba.toFixed(2)}" y="${y + 26}" text-anchor="start" class="rx-bound">locus start</text>`);
      if (ig.end <= hi) out.push(`<text x="${bb.toFixed(2)}" y="${y + 26}" text-anchor="middle" class="rx-bound">locus end</text>`);
    }

    const inWin = Object.entries(GEN).filter(([, g]) => g.chr === chr && g.end > lo && g.start < hi);
    for (const [sym, g] of inWin) {
      if (sym === (focus.content.gene || '')) continue;
      const gx = at(Math.max(g.start, lo)); const gw = Math.max(14, at(Math.min(g.end, hi)) - gx);
      out.push(geneArrow(gx, gw, y - 9, 18, '#6c7688', g.strand === '+', 0.5));
      out.push(`<text x="${(gx + gw / 2).toFixed(2)}" y="${y + 4}" text-anchor="middle" class="rx-gname-sm">${esc(sym)}</text>`);
    }
    if (fg) {
      geneStructure(out, fg, at, y, COL[chr], lo, hi);
      out.push(`<text x="${at((Math.max(fg.exons[0], lo) + Math.min(fg.exons[fg.exons.length - 1], hi)) / 2).toFixed(2)}" y="${y - 20}" text-anchor="middle" class="rx-gfocus" fill="${COL[chr]}">${esc(focus.content.gene)}</text>`);
    }

    // Dim everything this derivative does NOT carry, so the contributing stretch reads at a glance.
    // A segment kept on the "left" keeps the sequence at or below its breakpoint, so the region above
    // it is discarded, and the mirror for a "right" segment. With two breakpoints on one chromosome
    // only the sequence kept by BOTH survives, so the discarded regions are drawn one per segment.
    // This goes on after the genes and before the breakpoint rules, so the rules and their labels stay
    // legible on the boundary they mark.
    for (const seg of segs) {
      const bx = at(seg.pos);
      const from = seg.side === 'left' ? bx : x0;
      const to = seg.side === 'left' ? x0 + inner : bx;
      if (to - from <= 0.5) continue;
      out.push(`<rect x="${from.toFixed(2)}" y="${y - 34}" width="${(to - from).toFixed(2)}" height="68" fill="#0f1420" fill-opacity="0.68"/>`);
    }
    // Said once for the whole figure, not once per row.
    if (ri === 0) out.push(`<text x="${(x0 + inner).toFixed(2)}" y="${y - 34}" text-anchor="end" class="rx-bound">dimmed: not carried on this derivative</text>`);

    for (const seg of segs) out.push(breakMark(at(seg.pos), y - 30, y + 30, `chr${esc(chr)}:${num(seg.pos)}`));
    outAll += out.join('');
  });
  return `<svg class="rx" viewBox="0 0 ${width} ${H}" role="img" aria-label="Reference loci">${outAll}</svg>`;
}

/**
 * A miniature circos: every chromosome as an arc to genomic scale, the involved pair highlighted in
 * the panel colours, and a chord joining the two breakpoints. Pure geometry from the banding table.
 */
export function miniCircos(ctx, junctions, colours, size = 200) {
  const chroms = ctx.cytobands && ctx.cytobands.chromosomes;
  if (!chroms || !junctions.length) return '';
  const order = [...Array(22)].map((_, i) => String(i + 1)).concat(['X', 'Y']).filter((c) => chroms[c]);
  const total = order.reduce((a, c) => a + chroms[c].length, 0);
  const gap = 1.6; const usable = 360 - gap * order.length;
  // Every chromosome is labelled, so the ring is pulled in to leave room for the ring of numbers.
  const cx = size / 2; const cy = size / 2; const r = size / 2 - 22;
  const start = {}; let a = -90;
  for (const c of order) { start[c] = a; a += (chroms[c].length / total) * usable + gap; }
  const pol = (deg, rr) => { const t = (deg * Math.PI) / 180; return [cx + rr * Math.cos(t), cy + rr * Math.sin(t)]; };

  let out = '';
  for (const c of order) {
    const sweep = (chroms[c].length / total) * usable;
    const [x1, y1] = pol(start[c], r); const [x2, y2] = pol(start[c] + sweep, r);
    const hot = colours[c];
    out += `<path d="M${x1.toFixed(1)} ${y1.toFixed(1)} A${r} ${r} 0 ${sweep > 180 ? 1 : 0} 1 ${x2.toFixed(1)} ${y2.toFixed(1)}" fill="none" stroke="${hot || '#2b3652'}" stroke-width="${hot ? 8 : 5}"/>`;
    // Every chromosome is named, so the reader can place the event without counting arcs. The ones
    // this junction does not touch stay grey and small, and the involved pair takes the panel colours.
    const [lx, ly] = pol(start[c] + sweep / 2, r + (hot ? 11 : 9));
    out += `<text x="${lx.toFixed(1)}" y="${(ly + 3).toFixed(1)}" text-anchor="middle" class="${hot ? 'circ-lab' : 'circ-lab-off'}"${hot ? ` fill="${hot}"` : ''}>${esc(c)}</text>`;
  }
  const angleAt = (c, pos) => start[c] + (Math.max(0, Math.min(pos, chroms[c].length)) / chroms[c].length) * (chroms[c].length / total) * usable;
  // The chord is a ribbon rather than a line: broad where it meets each chromosome, narrowing through
  // the middle. Both edges are quadratic curves sharing the centre as their control point, so the
  // taper falls out of the geometry — the width at each end is the only number chosen.
  const CHORD_HALF_WIDTH = 3.4;
  for (const j of junctions) {
    if (!chroms[j.aChr] || !chroms[j.bChr]) continue;
    const aDeg = angleAt(j.aChr, j.aPos); const bDeg = angleAt(j.bChr, j.bPos);
    const [x1, y1] = pol(aDeg, r - 6);
    const [x2, y2] = pol(bDeg, r - 6);
    // The width is laid along the arc at each end, so the ribbon meets the chromosome square on.
    const tangent = (deg) => { const t = (deg * Math.PI) / 180; return [-Math.sin(t), Math.cos(t)]; };
    const [ax, ay] = tangent(aDeg); const [bx, by] = tangent(bDeg);
    const p = (x, y, dx, dy, s) => `${(x + dx * s * CHORD_HALF_WIDTH).toFixed(1)} ${(y + dy * s * CHORD_HALF_WIDTH).toFixed(1)}`;
    out += `<path d="M${p(x1, y1, ax, ay, 1)} Q${cx} ${cy} ${p(x2, y2, bx, by, 1)}`
      + ` L${p(x2, y2, bx, by, -1)} Q${cx} ${cy} ${p(x1, y1, ax, ay, -1)} Z" fill="#ff5f57" fill-opacity="0.8"/>`
      + `<circle cx="${x1.toFixed(1)}" cy="${y1.toFixed(1)}" r="2.2" fill="#ff5f57"/><circle cx="${x2.toFixed(1)}" cy="${y2.toFixed(1)}" r="2.2" fill="#ff5f57"/>`;
  }
  return `<svg class="circos" viewBox="0 0 ${size} ${size}" role="img" aria-label="Genome overview with the translocation chord">${out}</svg>`;
}
