// Detailed SpliceAI via the Broad SpliceAI-lookup API. Returns the four donor/acceptor delta scores
// AND their genomic positions, the REF/ALT probabilities, the exon model and the per-base curve —
// everything needed to draw the splice impact on an exon diagram.
//
// PRIVACY: sends the variant coordinate to the Broad API — same category as the Assessment/GeneBe
// feature, not the public-accession-only path.
// LICENCE: SpliceAI scores are CC BY-NC 4.0 (non-commercial). Fine for research / internal use; a
// commercial release would need an Illumina licence.
// FRAGILE: the canonical host is often unreachable; a Cloud Run backend answers. Hosts can change,
// so the base list is ordered fallbacks and should be overridable from settings.

const DEFAULT_HOSTS = [
  'https://spliceailookup-api.broadinstitute.org',
  'https://spliceai-38-xwkwwwxdwq-uc.a.run.app',
];

async function tryGet(url, timeoutMs) {
  const ctrl = new AbortController();
  const timer = setTimeout(() => ctrl.abort(), timeoutMs);
  try {
    const res = await fetch(url, { signal: ctrl.signal, headers: { 'User-Agent': 'CodonCompass' } });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    return await res.json();
  } finally { clearTimeout(timer); }
}

const base = (t) => String(t || '').replace(/\.\d+$/, '');

/**
 * @param {{chr:string,pos:number,ref:string,alt:string}} coord  GRCh38
 * @param {string} transcript  the app's RefSeq accession, to pick the matching isoform
 * @param {object} [opts] { hosts?: string[], distance?: number }
 */
export async function spliceaiLookup(coord, transcript, opts = {}) {
  // This codebase carries the chromosome under two names — `chr` on the splice path, `chrom` on the
  // assessment and gnomAD path — and a caller that supplies the other one produced the string
  // "chrundefined-81511114-G-A", which the service rejects with a bare HTTP 400. Accept either rather
  // than depend on every caller remembering which is which.
  const rawChr = coord.chr ?? coord.chrom;
  if (rawChr == null || coord.pos == null) {
    throw new Error('No genomic coordinate for this variant (chromosome or position missing).');
  }
  const chr = String(rawChr).replace(/^chr/i, '');
  const variant = `chr${chr}-${coord.pos}-${String(coord.ref).toUpperCase()}-${String(coord.alt).toUpperCase()}`;
  const distance = opts.distance || 500;
  const hosts = (opts.hosts && opts.hosts.length ? opts.hosts : DEFAULT_HOSTS);

  let data = null; let lastErr = null; let usedHost = null;
  for (const host of hosts) {
    try {
      // First host gets a short timeout so a dead canonical host fails fast to the fallback.
      data = await tryGet(`${host}/spliceai/?hg=38&distance=${distance}&mask=0&variant=${encodeURIComponent(variant)}`, host === hosts[0] ? 3500 : 30000);
      usedHost = host; break;
    } catch (e) { lastErr = e; }
  }
  if (!data) throw new Error(`SpliceAI unreachable (${lastErr ? lastErr.message : 'no hosts'})`);
  if (data.error) throw new Error(String(data.error));

  const scores = data.scores || [];
  // Pick the score set for the app's transcript; fall back to MANE Select, then the first.
  const chosen = scores.find((s) => (s.t_refseq_ids || []).some((r) => base(r) === base(transcript)))
    || scores.find((s) => s.t_priority === 'MS')
    || scores[0];
  if (!chosen) return { variant, usedHost, notScored: true };

  const num = (v) => (v == null ? null : Number(v));
  const events = [
    { key: 'AG', label: 'Acceptor gain', ds: num(chosen.DS_AG), dp: num(chosen.DP_AG) },
    { key: 'AL', label: 'Acceptor loss', ds: num(chosen.DS_AL), dp: num(chosen.DP_AL) },
    { key: 'DG', label: 'Donor gain', ds: num(chosen.DS_DG), dp: num(chosen.DP_DG) },
    { key: 'DL', label: 'Donor loss', ds: num(chosen.DS_DL), dp: num(chosen.DP_DL) },
  ];
  const maxDs = Math.max(0, ...events.map((e) => e.ds || 0));

  return {
    variant, usedHost,
    transcript: (chosen.t_refseq_ids && chosen.t_refseq_ids[0]) || chosen.t_id,
    priority: chosen.t_priority,
    strand: chosen.t_strand,
    pos: coord.pos,
    events,
    maxDs,
    // Absolute genomic site of each event = variant pos + DP (signed, strand-independent).
    sites: events.filter((e) => e.dp != null).map((e) => ({ ...e, site: coord.pos + e.dp })),
    exonStarts: chosen.EXON_STARTS || [],
    exonEnds: chosen.EXON_ENDS || [],
    cdsStart: chosen.CDS_START,
    cdsEnd: chosen.CDS_END,
    curve: data.allNonZeroScores || [],
    interpretation: maxDs >= 0.8 ? 'high' : maxDs >= 0.5 ? 'moderate' : maxDs >= 0.2 ? 'low' : 'minimal',
  };
}
