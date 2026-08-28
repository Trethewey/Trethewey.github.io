// gnomAD population-frequency lookup. Opt-in, network-only.
//
// Only a public genomic coordinate is sent — never patient data. Two public services are used:
//   1. Ensembl VEP HGVS endpoint converts the coding HGVS (e.g. NM_005373.3:c.1544G>T) to a GRCh38
//      genomic coordinate, because the app holds spliced coding sequence, not exon/genomic maps.
//      For an insertion or deletion a second Ensembl call fetches the reference bases to the left
//      of it, so the allele can be left-aligned to the position gnomAD files it under.
//   2. gnomAD's GraphQL API returns allele counts and frequencies for that coordinate.
//
// For a variant the app already has an hg38 coordinate for (Known/whitelist entries), step 1 is
// skipped.
//
// Step 1 also serves the Splice and ClinVar features through the main process, so it must handle
// insertions and deletions, not only substitutions. It does that by reading Ensembl's own VCF
// rendering of the variant rather than reconstructing one — see hgvsToGenomic.

import { leftAlign, needsLeftAlignment, variantId as idOf } from '../core/leftAlign.mjs';

const VEP_HGVS = 'https://rest.ensembl.org/vep/human/hgvs';
const ENSEMBL_SEQUENCE = 'https://rest.ensembl.org/sequence/region/human';
const GNOMAD_API = 'https://gnomad.broadinstitute.org/api';

// gnomAD v4 genetic-ancestry group labels.
const POP_LABELS = {
  afr: 'African / African-American', amr: 'Admixed American', asj: 'Ashkenazi Jewish',
  eas: 'East Asian', fin: 'Finnish', mid: 'Middle Eastern', nfe: 'Non-Finnish European',
  sas: 'South Asian', ami: 'Amish', remaining: 'Remaining',
};

async function getResponse(url, opts = {}) {
  const ctrl = new AbortController();
  const timer = setTimeout(() => ctrl.abort(), opts.timeoutMs || 30000);
  try {
    const res = await fetch(url, { ...opts, signal: ctrl.signal });
    if (!res.ok) throw new Error(`HTTP ${res.status} ${res.statusText}`);
    return res;
  } finally { clearTimeout(timer); }
}
const getJson = async (url, opts) => (await getResponse(url, opts)).json();
const getText = async (url, opts) => (await getResponse(url, opts)).text();

const COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A' };

// Ensembl's own VCF rendering of the variant, e.g. "7-140753336-A-T" or "17-7674968-AGAC-A".
// It is already on the plus strand and already carries the anchor base an insertion or deletion
// needs, so nothing here complements, reverses or invents a base.
function parseEnsemblVcfString(text) {
  const m = /^([0-9]{1,2}|X|Y|MT)-([0-9]+)-([ACGTN]+)-([ACGTN]+)$/.exec(String(text || '').trim().toUpperCase());
  if (!m) return null;
  const pos = Number(m[2]);
  if (!Number.isInteger(pos) || pos < 1) return null;
  return { chrom: m[1], pos, ref: m[3], alt: m[4] };
}

// Fallback for a single-base substitution when Ensembl sends no VCF rendering.
// Its HGVS endpoint reports the allele in the TRANSCRIPT strand, so on a minus-strand gene the two
// bases must be complemented to the plus-strand genome (gnomAD, SpliceAI and ClinVar all key on the
// plus strand). Without this, every minus-strand variant queries the wrong allele.
// Only plain single bases qualify: an insertion writes its missing side as "-", which is one
// character long but is not a base, and must not be passed off as one.
function substitutionFromAlleleString(r) {
  const [ref, alt] = String(r.allele_string || '').toUpperCase().split('/');
  if (!/^[ACGT]$/.test(ref || '') || !/^[ACGT]$/.test(alt || '')) return null;
  const flip = Number(r.strand) === -1;
  const pos = Number(r.start);
  if (!Number.isInteger(pos) || pos < 1) return null;
  return {
    chrom: String(r.seq_region_name),
    pos,
    ref: flip ? COMPLEMENT[ref] : ref,
    alt: flip ? COMPLEMENT[alt] : alt,
  };
}

const sameCoordinate = (a, b) => a.chrom === b.chrom && a.pos === b.pos && a.ref === b.ref && a.alt === b.alt;
const printCoordinate = (c) => `${c.chrom}:${c.pos} ${c.ref}>${c.alt}`;

/** Coding HGVS -> {chrom, pos, ref, alt} on GRCh38 PLUS strand, via Ensembl VEP.
 *
 *  Substitutions, insertions, deletions and combined deletion-insertions all resolve. The
 *  coordinate is read from Ensembl's own VCF rendering (`vcf_string`), which is why insertions and
 *  deletions work: that field carries the anchor base — the one unchanged base written into both
 *  alleles so neither side is empty — and it is already on the plus strand. Asking for it costs one
 *  extra field in the same reply to the same request; nothing further is sent off the machine. A
 *  5xx or a dropped connection is retried twice — Ensembl's REST service fails intermittently, and
 *  an unretried failure here loses the whole lookup rather than degrading it.
 *
 *  `placedBy: 'Ensembl'` records where the position came from. It matters for an insertion or
 *  deletion: Ensembl and other services do not always write the same allele at the same position
 *  inside a repeat (checked against VariantValidator — Ensembl puts CFTR p.Phe508del at
 *  7:117559591 TCTT>T, VariantValidator at 7:117559590 ATCT>A; both describe the same edit). The
 *  gnomAD lookup settles that with leftAlignForLookup below.
 */
export async function hgvsToGenomic(hgvs) {
  const url = `${VEP_HGVS}/${encodeURIComponent(hgvs)}?content-type=application/json&vcf_string=1`;
  const data = await retryTransient(() => getJson(url));
  const r = data && data[0];
  if (!r || !r.seq_region_name) throw new Error('Ensembl could not place this variant on the genome');
  if (r.assembly_name && String(r.assembly_name) !== 'GRCh38') {
    throw new Error(`Ensembl answered on ${r.assembly_name}; this app works in GRCh38 only`);
  }

  const fromVcf = parseEnsemblVcfString(r.vcf_string);
  const fromAlleles = substitutionFromAlleleString(r);
  // When both readings are available they must agree. If they do not, one of them is wrong and
  // there is no basis for choosing, so nothing is returned.
  if (fromVcf && fromAlleles && !sameCoordinate(fromVcf, fromAlleles)) {
    throw new Error(`Ensembl gave two different coordinates for this variant (${printCoordinate(fromVcf)} and ${printCoordinate(fromAlleles)}); refusing to choose between them`);
  }
  const coord = fromVcf || fromAlleles;
  if (!coord) throw new Error('Ensembl returned no usable genomic coordinate for this variant');
  return { ...coord, placedBy: 'Ensembl' };
}

/** Plus-strand GRCh38 reference bases for chrom:lo-hi (1-based, inclusive), from Ensembl.
 *  Ensembl's REST service answers 5xx and times out on and off; a server-side failure is retried
 *  twice before the window is given up, because giving up turns a findable allele into "absence not
 *  established". The browsers name a failed request differently from Node, so all four wordings are
 *  matched — the website runs this same code. */
const TRANSIENT = /HTTP 5\d\d|abort|fetch failed|failed to fetch|networkerror|load failed/i;
const RETRY_WAITS_MS = [500, 1000];
async function retryTransient(call) {
  for (let attempt = 0; ; attempt += 1) {
    try { return await call(); }
    catch (e) {
      if (!TRANSIENT.test(e.message) || attempt >= RETRY_WAITS_MS.length) throw e;
      await new Promise((r) => setTimeout(r, RETRY_WAITS_MS[attempt]));
    }
  }
}
async function referenceWindow(chrom, lo, hi) {
  const url = `${ENSEMBL_SEQUENCE}/${chrom}:${lo}..${hi}?content-type=text/plain;coord_system_version=GRCh38`;
  const seq = await retryTransient(async () => String(await getText(url)).trim().toUpperCase());
  if (!/^[ACGTN]+$/.test(seq) || seq.length !== hi - lo + 1) {
    throw new Error(`Ensembl returned ${seq.length} bases for ${chrom}:${lo}-${hi}, expected ${hi - lo + 1}`);
  }
  return { start: lo, seq };
}

// How far left of the variant to read. A polyglutamine tract is tens of bases and a long
// microsatellite hundreds; the window widens twice before the walk is abandoned.
const WINDOW_WIDTHS = [300, 3000, 30000];

/** Move an insertion or deletion to its left-most position (see core/leftAlign.mjs). A one-base
 *  substitution is returned as it came, with no request made. Throws when the walk cannot be
 *  completed — the caller decides what an unaligned allele means for its answer. */
export async function leftAlignForLookup(coord) {
  const chrom = String(coord.chrom ?? coord.chr).replace(/^chr/i, '');
  const c = { chrom, pos: Number(coord.pos), ref: String(coord.ref || '').toUpperCase(), alt: String(coord.alt || '').toUpperCase() };
  if (!needsLeftAlignment(c.ref, c.alt)) return leftAlign(c, null);
  let last = null;
  for (const width of WINDOW_WIDTHS) {
    const lo = Math.max(1, c.pos - width);
    const window = await referenceWindow(chrom, lo, c.pos + c.ref.length);
    last = leftAlign(c, window);
    if (last.ok) return last;
    if (last.code !== 'window' || lo === 1) break;
  }
  throw new Error(last && last.reason ? last.reason : 'left-alignment failed');
}

/** Query gnomAD for a genomic coordinate. Returns null if the variant is not in gnomAD. */
async function queryGnomad(chrom, pos, ref, alt, dataset = 'gnomad_r4') {
  const variantId = `${chrom.replace(/^chr/i, '')}-${pos}-${ref}-${alt}`;
  const query = `query($v:String!,$d:DatasetId!){
    variant(variantId:$v,dataset:$d){
      variant_id rsids
      flags
      genome{ ac an af populations{ id ac an } }
      exome{ ac an af populations{ id ac an } }
    }
  }`;
  const data = await getJson(GNOMAD_API, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ query, variables: { v: variantId, d: dataset } }),
    timeoutMs: 30000,
  });
  if (data.errors && (!data.data || !data.data.variant)) {
    // A "not found" comes back as an error message; treat genuinely-absent as null, re-throw others.
    const msg = data.errors.map((e) => e.message).join('; ');
    if (/not found|no variant/i.test(msg)) return { variantId, notInGnomad: true };
    throw new Error(msg);
  }
  const v = data.data && data.data.variant;
  if (!v) return { variantId, notInGnomad: true };

  // Merge exome + genome allele counts for a combined frequency, and per-population.
  const perPop = new Map();
  const addPops = (src) => { for (const p of (src && src.populations) || []) {
    const cur = perPop.get(p.id) || { id: p.id, ac: 0, an: 0 };
    cur.ac += p.ac || 0; cur.an += p.an || 0; perPop.set(p.id, cur);
  } };
  addPops(v.exome); addPops(v.genome);

  const populations = [...perPop.values()]
    .filter((p) => POP_LABELS[p.id])
    .map((p) => ({ id: p.id, label: POP_LABELS[p.id], ac: p.ac, an: p.an, af: p.an ? p.ac / p.an : 0 }))
    .sort((a, b) => b.af - a.af);

  const totalAc = (v.exome ? v.exome.ac : 0) + (v.genome ? v.genome.ac : 0);
  const totalAn = (v.exome ? v.exome.an : 0) + (v.genome ? v.genome.an : 0);
  const overallAf = totalAn ? totalAc / totalAn : 0;
  const popmax = populations.length ? populations[0] : null;

  return {
    variantId: v.variant_id || variantId,
    rsids: v.rsids || [],
    flags: v.flags || [],
    dataset,
    exome: v.exome ? { ac: v.exome.ac, an: v.exome.an, af: v.exome.af } : null,
    genome: v.genome ? { ac: v.genome.ac, an: v.genome.an, af: v.genome.af } : null,
    overall: { ac: totalAc, an: totalAn, af: overallAf },
    popmax,
    populations,
    notInGnomad: false,
  };
}

/** Full lookup from a coding HGVS (optionally with a known hg38 coordinate to skip VEP).
 *
 *  gnomAD files an insertion or deletion at one position only — the left-most, the VCF convention —
 *  and Ensembl writes the same allele at the 3′-most one, as HGVS does. Checked against the live
 *  services on 2026-08-28: gnomAD holds SMARCA2 c.705_707dup as 9-2039776-A-ACAG (about 1 in 30
 *  people carry it) and answers "Variant not found" for Ensembl's 9-2039814-A-AGCA; CFTR
 *  p.Phe508del and ASXL1 c.1934dup fail the same way. So the allele is left-aligned before it is
 *  queried, and the answer carries both spellings: `coord` is the coordinate as placed, `queried`
 *  the one gnomAD was asked about, `leftAligned` {from, to} when the two differ.
 */
export async function gnomadLookup(hgvs, knownCoord) {
  let coord = knownCoord;
  if (!coord) coord = await hgvsToGenomic(hgvs);
  const asPlaced = {
    chrom: String(coord.chrom ?? coord.chr).replace(/^chr/i, ''),
    pos: Number(coord.pos),
    ref: String(coord.ref || '').toUpperCase(),
    alt: String(coord.alt || '').toUpperCase(),
  };

  let queried = asPlaced;
  let leftAligned = null;
  let alignmentError = null;
  if (needsLeftAlignment(asPlaced.ref, asPlaced.alt)) {
    try {
      const a = await leftAlignForLookup(asPlaced);
      queried = { chrom: a.chrom, pos: a.pos, ref: a.ref, alt: a.alt };
      if (a.moved) leftAligned = { from: idOf(asPlaced), to: idOf(queried) };
    } catch (e) { alignmentError = e.message; }
  }

  // The position is worked out before gnomAD is asked, and the ClinVar card is keyed on it, so a
  // gnomAD outage must not take the position down with it: the failure is returned beside `queried`
  // rather than thrown.
  let result;
  try {
    result = await queryGnomad(queried.chrom, queried.pos, queried.ref, queried.alt);
  } catch (e) {
    return { unavailable: true, error: e.message, provisional: true, coord, queried, leftAligned, hgvs };
  }

  // A miss is a finding only when the allele was asked about in the form gnomAD keys on. When the
  // walk to that form failed, or the change swaps one run of bases for another (gnomAD's variant
  // table holds only single-base substitutions and anchored insertions and deletions, never a
  // multi-base substitution or deletion-insertion), "not found here" is all that has been
  // established.
  if (result.notInGnomad) {
    const multiBase = queried.ref.length > 1 && queried.alt.length > 1;
    if (alignmentError) {
      result.absenceUnconfirmed = true;
      result.absenceNote = `Not found as ${result.variantId}. Left-alignment failed: ${alignmentError}. Absence not established.`;
    } else if (multiBase) {
      result.absenceUnconfirmed = true;
      result.absenceNote = `Not found as ${result.variantId}. gnomAD does not file a change of one run of bases for another as a single record. Absence not established.`;
    }
  }
  // A provisional answer is one the position could not be settled for: it must not be cached, and it
  // must not be read as a finding about the variant.
  return { ...result, coord, queried, leftAligned, hgvs, provisional: !!alignmentError, alignmentError };
}
