// ClinVar, queried directly via NCBI E-utilities (not via GeneBe). Gives the modern three-way split:
// germline classification, somatic ONCOGENICITY, and somatic CLINICAL IMPACT (tier), each with a
// review-star count, plus the per-submission list.
//
// PRIVACY: querying by coordinate/rsID sends the variant location to NCBI — the same category as the
// GeneBe/Assessment feature, not the public-accession-only path. Gated and cached in the main process.
//
// Correctness note: a coordinate/rsID is multi-allelic (rs113488022 -> V600E AND V600G), so the
// candidate records are disambiguated by matching the variant's SPDI (ref/alt at position) before
// anything is shown — the same exactness rule the app applies to transcript versions.
//
// SPDI is NCBI's sequence-position-deletion-insertion notation. It is NOT the same shape as the
// VCF-style coordinate the rest of the app carries, so insertions and deletions have to be lined
// up before they can be compared — see spdiMatchesCoordinate below.

const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
const TOOL = 'tool=CodonCompass&email=codoncompass@local';

// Review-status text -> star count (ClinVar's documented scale). The API returns only the text.
const STARS = {
  'practice guideline': 4,
  'reviewed by expert panel': 3,
  'criteria provided, multiple submitters, no conflicts': 2,
  'criteria provided, multiple submitters': 2,
  'criteria provided, single submitter': 1,
  'criteria provided, conflicting classifications': 1,
  'no assertion criteria provided': 0,
  'no classification provided': 0,
};
export function reviewStars(status) {
  const key = String(status || '').trim().toLowerCase();
  return STARS[key] != null ? STARS[key] : 0;
}

// NCBI E-utilities rate-limit anonymous callers (3 requests a second) and return 429 when exceeded;
// they also emit transient 5xx. Retry those a couple of times with a short backoff so a burst
// doesn't surface as a spurious "not found". Every other HTTP status fails immediately: the request
// itself is what it objects to, so repeating it only makes the user wait.
//
// The status decision is made on the response, never by reading it back off an error message. The
// old version threw `HTTP <status>` for a non-retryable status and then tested the message with
// /^HTTP (?!429|5\d\d)/, which suppressed the rethrow for EVERY 5xx — so 501 and 505, deliberately
// left out of RETRY_STATUS, were retried three times and the error arrived 7 seconds late.
const RETRY_STATUS = new Set([429, 500, 502, 503, 504]);
const ATTEMPTS = 3;
async function fetchRetry(url, timeoutMs, parse) {
  let lastErr;
  for (let attempt = 0; attempt < ATTEMPTS; attempt += 1) {
    const ctrl = new AbortController();
    const timer = setTimeout(() => ctrl.abort(), timeoutMs);
    let fatal = null; // an error another attempt cannot fix
    try {
      const res = await fetch(url, { signal: ctrl.signal, headers: { 'User-Agent': 'CodonCompass' } });
      if (res.ok) return await parse(res);
      lastErr = new Error(`HTTP ${res.status}`);
      if (!RETRY_STATUS.has(res.status)) fatal = lastErr;
    } catch (e) {
      // A dropped connection, a timeout, or an unreadable body: worth another try.
      lastErr = e;
    } finally { clearTimeout(timer); }
    if (fatal) throw fatal;
    // Wait before trying again — but not after the last attempt, where the wait only delays the
    // error that is about to be thrown anyway.
    if (attempt < ATTEMPTS - 1) await new Promise((r) => setTimeout(r, 1200 * (attempt + 1)));
  }
  throw lastErr;
}
async function getJson(url, timeoutMs = 25000) { return fetchRetry(url, timeoutMs, (res) => res.json()); }
async function getText(url, timeoutMs = 30000) { return fetchRetry(url, timeoutMs, (res) => res.text()); }

function pack(c) {
  if (!c || !c.description) return null;
  return { description: c.description, reviewStatus: c.review_status || '', stars: reviewStars(c.review_status), lastEvaluated: c.last_evaluated || '' };
}

// ---- lining a VCF-style coordinate up with a ClinVar canonical SPDI ----
//
// The two notations describe the same edit in different ways:
//
//   VCF style (what the rest of the app carries): 1-based position, and an insertion or deletion
//     carries an ANCHOR BASE — one unchanged base written into both the reference and the
//     alternate allele so neither field is empty. NPM1's type A duplication is 5:171410539 C>CTCTG.
//
//   SPDI (sequence:position:deletion:insertion, what ClinVar publishes): 0-based, no anchor base,
//     and FULLY JUSTIFIED — inside a repeated stretch the deletion and insertion fields are both
//     expanded to cover the whole repeat. The same NPM1 duplication is
//     NC_000005.10:171410539:TCTG:TCTGTCTG.
//
// Comparing the fields literally therefore only ever works for a substitution, where there is no
// anchor base and nothing to expand.

const BASES = /^[ACGTN]*$/;

// The GRCh38 chromosome sequences, by their exact RefSeq accession INCLUDING the version.
//
// Every position this app handles is GRCh38, so an SPDI can only be lined up against one when it is
// written on one of these 24 sequences. The version is part of the test, not decoration: NC_000007.13
// is chromosome 7 of GRCh37, a different sequence with different coordinates, and treating it as
// chromosome 7 of GRCh38 would compare two positions that are not the same place. Anything else —
// another assembly, an unplaced scaffold, a gene or transcript accession, or the mitochondrion,
// which the app holds no transcript for — cannot be placed and is refused rather than guessed.
//
// The list is the set of distinct chromosome accessions in the GRCh38_chr column of the bundled
// NCBI MANE summary (output/refdb-build/MANE.summary.txt.gz). The control test re-derives it from
// that file and fails if the two ever differ.
const GRCH38_CHROMOSOMES = {
  'NC_000001.11': '1', 'NC_000002.12': '2', 'NC_000003.12': '3', 'NC_000004.12': '4',
  'NC_000005.10': '5', 'NC_000006.12': '6', 'NC_000007.14': '7', 'NC_000008.11': '8',
  'NC_000009.12': '9', 'NC_000010.11': '10', 'NC_000011.10': '11', 'NC_000012.12': '12',
  'NC_000013.11': '13', 'NC_000014.9': '14', 'NC_000015.10': '15', 'NC_000016.10': '16',
  'NC_000017.11': '17', 'NC_000018.10': '18', 'NC_000019.10': '19', 'NC_000020.11': '20',
  'NC_000021.9': '21', 'NC_000022.11': '22', 'NC_000023.11': 'X', 'NC_000024.10': 'Y',
};
export const GRCH38_CHROMOSOME_ACCESSIONS = GRCH38_CHROMOSOMES;

// RefSeq GRCh38 chromosome accession -> chromosome name. Null for everything else.
function chromosomeOfAccession(accession) {
  const key = String(accession || '').trim().toUpperCase();
  return Object.prototype.hasOwnProperty.call(GRCH38_CHROMOSOMES, key) ? GRCH38_CHROMOSOMES[key] : null;
}

/**
 * Trim a VCF-style allele down to the bases that actually change, dropping the anchor base.
 * Returns 0-based coordinates: `start` is the first changed base, `deleted` the reference bases
 * removed, `inserted` what replaces them. Either string may be empty (a pure insertion has no
 * deleted bases; a pure deletion has no inserted ones).
 * The shared tail is trimmed before the shared head, which places the change as far left as the
 * two alleles allow — the same end the fully justified SPDI position starts from.
 * @returns {{start:number,deleted:string,inserted:string}|null} null if the alleles are not plain bases.
 */
export function minimalAllele(pos, ref, alt) {
  let r = String(ref || '').toUpperCase();
  let a = String(alt || '').toUpperCase();
  if (!r || !a || !BASES.test(r) || !BASES.test(a)) return null;
  let start = Number(pos) - 1; // 1-based -> 0-based
  if (!Number.isInteger(start) || start < 0) return null;
  while (r.length && a.length && r[r.length - 1] === a[a.length - 1]) { r = r.slice(0, -1); a = a.slice(0, -1); }
  while (r.length && a.length && r[0] === a[0]) { r = r.slice(1); a = a.slice(1); start += 1; }
  return { start, deleted: r, inserted: a };
}

function parseSpdi(spdi) {
  const parts = String(spdi || '').split(':');
  if (parts.length !== 4) return null;
  const start = Number(parts[1]);
  if (!Number.isInteger(start) || start < 0) return null;
  return { sequence: parts[0], start, del: parts[2].toUpperCase(), ins: parts[3].toUpperCase() };
}

/**
 * Does a ClinVar canonical SPDI describe exactly the edit the caller asked about?
 *
 * Both edits are applied to the same stretch of reference sequence and the results compared. If
 * they read the same, the two notations describe the same allele; if they differ, they do not.
 * The reference bases come only from the two records themselves — the caller's reference allele,
 * and the SPDI's deletion field, which IS the reference over the span it covers — so nothing is
 * assumed about sequence the app has not been handed. Where the two spans leave a gap, the
 * reference in between is unknown and the answer is 'unknown', never a match.
 *
 * @param {string} spdi   ClinVar's canonical_spdi, e.g. 'NC_000007.14:117559590:TCTT:T'
 * @param {{chr:string,pos:number,ref:string,alt:string}} coord  VCF-style, 1-based, anchor base included
 * @returns {'match'|'different'|'unknown'}
 */
export function spdiMatchesCoordinate(spdi, coord) {
  const s = parseSpdi(spdi);
  if (!s) return 'unknown';

  const spdiChr = chromosomeOfAccession(s.sequence);
  if (!spdiChr) return 'unknown'; // not a genomic accession: cannot be placed against this position
  const wantChr = String(coord.chr || '').replace(/^chr/i, '').toUpperCase().replace(/^M$/, 'MT');
  if (spdiChr !== wantChr) return 'different';

  const q = minimalAllele(coord.pos, coord.ref, coord.alt);
  if (!q) return 'unknown';

  // The deletion field is allowed to be a plain length instead of the deleted sequence. Then the
  // reference bases are unknown and only an exact match of position, deleted length and inserted
  // sequence can be trusted; anything else is undecidable rather than a mismatch.
  if (/^\d+$/.test(s.del)) {
    return (s.start === q.start && Number(s.del) === q.deleted.length && s.ins === q.inserted) ? 'match' : 'unknown';
  }
  if (!BASES.test(s.del) || !BASES.test(s.ins)) return 'unknown';

  // Reference bases known from the caller: [q.start, q.start + q.deleted.length).
  // Reference bases known from the SPDI:   [s.start, s.start + s.del.length).
  // Build the window spanning both. A base known from neither means a gap we cannot fill.
  const lo = Math.min(q.start, s.start);
  const hi = Math.max(q.start + q.deleted.length, s.start + s.del.length);
  let reference = '';
  for (let i = lo; i < hi; i += 1) {
    const fromCoord = (i >= q.start && i < q.start + q.deleted.length) ? q.deleted[i - q.start] : null;
    const fromSpdi = (i >= s.start && i < s.start + s.del.length) ? s.del[i - s.start] : null;
    if (fromCoord === null && fromSpdi === null) return 'unknown';                    // gap
    if (fromCoord !== null && fromSpdi !== null && fromCoord !== fromSpdi) return 'different'; // references disagree
    reference += fromCoord !== null ? fromCoord : fromSpdi;
  }
  const apply = (start, deletedLength, inserted) =>
    reference.slice(0, start - lo) + inserted + reference.slice(start - lo + deletedLength);
  return apply(q.start, q.deleted.length, q.inserted) === apply(s.start, s.del.length, s.ins) ? 'match' : 'different';
}

// A record whose variation_set holds more than one allele is a haplotype or compound record — a
// COMBINATION of changes, so it is never the single allele that was asked about, even though its
// first allele can match. ClinVar's CFTR c.[1521_1523delCTT;3080T>C] shares its first allele with
// plain p.Phe508del and would otherwise be reported in its place.
function singleAlleleSet(r) {
  const set = (r && r.variation_set) || [];
  return set.length === 1 ? set[0] : null;
}

/**
 * Find and summarise the ClinVar record for a GRCh38 coordinate, disambiguated to the exact allele.
 * @param {{chr:string,pos:number,ref:string,alt:string}} coord
 * @param {string} [gene]  gene symbol, narrows the search
 * @returns {object|null}
 */
export async function clinvarByCoordinate(coord, gene) {
  const chr = String(coord.chr).replace(/^chr/i, '');
  // ClinVar files a record under the positions it actually changes. A VCF-style deletion is written
  // with an anchor base one place BEFORE the change, so searching that single position misses the
  // record entirely (CFTR p.Phe508del, written 7:117559590 ATCT>A, is filed at 117559591-117559593
  // and was invisible). Search the whole span the reference allele covers instead. For a
  // substitution or an insertion the reference allele is one base and this is the old single
  // position. If ClinVar ever returns more than retmax candidates the allele may go unfound — that
  // costs a lookup, never a wrong answer, because the allele match below is exact.
  const spanEnd = coord.pos + Math.max(String(coord.ref || '').length, 1) - 1;
  const posTerm = spanEnd > coord.pos ? `${coord.pos}:${spanEnd}[chrpos38]` : `${coord.pos}[chrpos38]`;
  const term = encodeURIComponent(`${chr}[chr] AND ${posTerm}${gene ? ` AND ${gene}[gene]` : ''}`);
  const search = await getJson(`${EUTILS}/esearch.fcgi?db=clinvar&term=${term}&retmode=json&retmax=200&${TOOL}`);
  const ids = (search.esearchresult && search.esearchresult.idlist) || [];
  if (!ids.length) return { notFound: true };

  const sum = await getJson(`${EUTILS}/esummary.fcgi?db=clinvar&id=${ids.join(',')}&retmode=json&${TOOL}`);
  const result = sum.result || {};

  // esummary can fail for one identifier inside an otherwise fine HTTP 200 reply: that identifier's
  // entry comes back as {"error":"cannot get document summary"} instead of a record, and an
  // identifier can be missing from the reply altogether. Neither is a record, and neither says
  // anything about the allele. Left in, an errored entry carries no variation_set, so it satisfied
  // the lone-candidate carve-out below and was reported as a real ClinVar record with nothing
  // classified — a transient NCBI failure shown to the user as "listed in ClinVar, unclassified".
  const isRecord = (id) => {
    const r = result[id];
    return !!r && typeof r === 'object' && !r.error;
  };
  const unreadable = ids.filter((id) => !isRecord(id));

  // Disambiguate: keep only the candidate whose canonical SPDI is the same edit as the queried
  // allele. spdiMatchesCoordinate handles the anchor base and SPDI's repeat expansion, so this
  // works for insertions and deletions as well as substitutions, and it stays allele-exact:
  // 'unknown' (two notations that cannot be lined up) is rejected exactly like 'different'.
  let chosen = null;
  for (const id of ids) {
    if (!isRecord(id)) continue;
    const entry = singleAlleleSet(result[id]);
    if (!entry || !entry.canonical_spdi) continue;
    if (spdiMatchesCoordinate(entry.canonical_spdi, coord) === 'match') { chosen = { id, r: result[id], entry }; break; }
  }
  // Accept a lone candidate ONLY when it has no SPDI to check against (can't disambiguate). If its
  // SPDI is present it already failed the allele match above, which means the ClinVar record at this
  // position is a DIFFERENT allele than the one queried — reporting it would show the wrong variant,
  // so treat this allele as absent from ClinVar. A haplotype record is rejected here too: it is a
  // combination of alleles, so it is not the single allele that was queried whether or not it
  // carries an SPDI.
  if (!chosen && ids.length === 1 && isRecord(ids[0])) {
    const r0 = result[ids[0]];
    const set = r0.variation_set || [];
    if (set.length <= 1 && !(set[0] && set[0].canonical_spdi)) chosen = { id: ids[0], r: r0, entry: set[0] || null };
  }
  if (!chosen) {
    // Nothing matched. If a candidate could not be read, we do not know what allele it described,
    // so "this allele is not in ClinVar" is a statement the lookup has not earned. Fail plainly
    // instead — the card then says ClinVar is unavailable, and nothing is written to the cache.
    if (unreadable.length) {
      const shown = unreadable.slice(0, 3).join(', ') + (unreadable.length > 3 ? ` and ${unreadable.length - 3} more` : '');
      throw new Error(`no summary returned for record${unreadable.length > 1 ? 's' : ''} ${shown}, so whether this allele is listed is unknown`);
    }
    return { notFound: true, candidates: ids.length };
  }

  const r = chosen.r;
  return {
    variationId: chosen.id,
    title: r.title,
    accession: r.accession,
    url: `https://www.ncbi.nlm.nih.gov/clinvar/variation/${chosen.id}/`,
    germline: pack(r.germline_classification),
    oncogenicity: pack(r.oncogenicity_classification),
    clinicalImpact: pack(r.clinical_impact_classification),
    submissionCount: (r.supporting_submissions && r.supporting_submissions.scv && r.supporting_submissions.scv.length) || 0,
    spdi: chosen.entry && chosen.entry.canonical_spdi,
    cdna: chosen.entry && chosen.entry.cdna_change,
  };
}

// ---- per-submission detail (lazy: only when the user expands the table) ----

function firstTag(block, tag) {
  const m = new RegExp(`<${tag}[^>]*>([\\s\\S]*?)</${tag}>`).exec(block);
  return m ? m[1].replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim() : '';
}
function attr(block, tag, name) {
  const m = new RegExp(`<${tag}\\b[^>]*\\b${name}="([^"]*)"`).exec(block);
  return m ? m[1] : '';
}

/** Per-submission rows for a VariationID, via efetch VCV XML (is_variationid is mandatory). */
export async function clinvarSubmissions(variationId) {
  const xml = await getText(`${EUTILS}/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=${encodeURIComponent(variationId)}&${TOOL}`);
  const rows = [];
  const re = /<ClinicalAssertion\b[\s\S]*?<\/ClinicalAssertion>/g;
  let m;
  while ((m = re.exec(xml))) {
    const block = m[0];
    const submitter = attr(block, 'ClinVarAccession', 'SubmitterName') || attr(block, 'ClinVarAccession', 'OrgID');
    const scv = attr(block, 'ClinVarAccession', 'Accession');
    // The classification tag inside <Classification> tells us germline from somatic. A submission
    // that carries none of the three tags is not a germline one: it is a submission with no
    // classification, and saying "germline" would invent the one fact the tag was there to give.
    // Such a row still reaches the table (it has a submitter, a review status and a condition), so
    // it is typed by what is actually there.
    let kind = 'not stated'; let value = '';
    const cls = /<Classification\b[\s\S]*?<\/Classification>/.exec(block);
    const c = cls ? cls[0] : block;
    if (/<SomaticClinicalImpact/.test(c)) { kind = 'somatic clinical impact'; value = firstTag(c, 'SomaticClinicalImpact'); }
    else if (/<OncogenicityClassification/.test(c)) { kind = 'somatic oncogenicity'; value = firstTag(c, 'OncogenicityClassification'); }
    else if (/<GermlineClassification/.test(c)) { kind = 'germline'; value = firstTag(c, 'GermlineClassification'); }
    const reviewStatus = firstTag(c, 'ReviewStatus');
    const condition = firstTag(block, 'ElementValue') || firstTag(block, 'Name');
    const date = attr(block, 'ClinicalAssertion', 'SubmissionDate') || attr(block, 'Description', 'DateLastEvaluated');
    if (submitter || value) rows.push({ submitter, scv, kind, value, reviewStatus, stars: reviewStars(reviewStatus), condition, date });
  }
  return rows;
}
