import { STANDARD_CODE, translate, translateCodon, reverseComplement, toMrna } from '../core/geneticCode.mjs';
import { AMINO_ACIDS, aaByOne, three, fullName, POLARITY_COLOUR, AROMATICITY_NOTE } from '../core/aminoAcids.mjs';
import { parseHgvs } from '../core/hgvs.mjs';
import { clinvarByCoordinate, clinvarSubmissions } from '../net/clinvar.mjs';
import { resolveSubstitution, resolveIndel, consequenceLabel } from '../core/variant.mjs';
import { KNOWN_VARIANTS } from '../core/knownVariants.mjs';
import { maneToUniprot } from '../core/proteinMap.mjs';
import { resolveCodingCoordinate, spliceRegion, SPLICE_CONVENTION } from '../core/codingCoordinate.mjs';
import { AA_STRUCTURES } from '../core/aaStructures.mjs';
import { parseReport, detectEvents, describeDerivatives, partShort, parseGermlineAf } from '../core/fusion.mjs';
import { ideogram, derivativeMolecule, referenceLoci, miniCircos, chromosomeColours, chromosomeOrder } from './fusionDraw.mjs';
import { FUSION_HAEM_EXAMPLES, FUSION_HAEM_EXPECTATIONS } from '../core/haemFusions.mjs';

// Naming a textbook here was wrong: the four-way polarity grouping below is widely used, but the
// textbook it was credited to splits the side chains a different way (it gives the aromatic ones a
// class of their own). The grouping is described rather than attributed, so nothing is claimed about
// a source that was not checked.
const SCHEME = 'Polarity uses the common four-way biochemistry grouping (nonpolar / polar uncharged / acidic / basic); '
  + 'textbooks divide the side chains differently, some giving the aromatic ones a class of their own. '
  + 'Hydropathy is the Kyte-Doolittle index (1982). Average residue masses are the ExPASy values (residue = amino acid minus water). '
  + 'Charge is the net side-chain charge at physiological pH ~7.4.';

// The residue-mass convention, spelled out wherever a mass is shown outside the Amino acids tab —
// these are residue masses (the free amino acid minus water), so each is about 18 Da lighter than
// the free amino acid a reader may be comparing against.
const RESIDUE_MASS_HELP = 'Residue mass: the amino acid as it sits in a chain, i.e. the free amino '
  + 'acid minus water — about 18 Da lighter. ExPASy average values.';

const BASES4 = ['T', 'C', 'A', 'G'];
const DOMAIN_PALETTE = ['#4f83e2', '#37a89a', '#d99a4e', '#a86cff', '#e0567f', '#5bb98c', '#c98a3a', '#6c8fe0', '#3aa8a0', '#b06cff'];

// ---------- tiny DOM helper ----------
function el(tag, attrs = {}, ...kids) {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v == null) continue;
    if (k === 'class') n.className = v;
    else if (k === 'html') n.innerHTML = v;
    else if (k.startsWith('on') && typeof v === 'function') n.addEventListener(k.slice(2), v);
    else n.setAttribute(k, v);
  }
  for (const kid of kids.flat()) {
    if (kid == null || kid === false) continue;
    n.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return n;
}
const $ = (sel, root = document) => root.querySelector(sel);

// True in the hosted / single-file builds. GeneBe's API refuses browser (cross-origin) requests, so
// it is desktop-only; the renderer hides it here rather than auto-running it and showing an error.
const IS_WEB = typeof window !== 'undefined' && window.__CC_WEB === true;

function baseCell(b, { changed = false, big = false } = {}) {
  const u = String(b).toUpperCase();
  return el('div', { class: `base ${u}${changed ? ' changed' : ''}`, style: big ? '' : '' }, u);
}

// ---------- amino-acid "big" tile ----------
function aaTile(oneLetter, label) {
  const aa = aaByOne(oneLetter);
  const isStop = oneLetter === '*';
  const colour = aa ? (aa.colour || POLARITY_COLOUR[aa.polarity] || '#8a8f99') : '#8a8f99';
  return el('div', { class: 'aa-big' },
    el('div', { class: 'name muted' }, label),
    el('div', { class: 'letter', style: `color:${colour}` }, oneLetter),
    el('div', { class: 'three' }, isStop ? 'Ter (stop)' : `${aa ? aa.three : '?'}`),
    el('div', { class: 'name' }, isStop ? 'Termination' : (aa ? aa.name : 'Unknown')),
    aa && aa.polarity ? el('div', { class: 'pol', style: `background:${colour}` }, aa.polarity) : null,
    AA_STRUCTURES[oneLetter] ? aaStructureEl(oneLetter) : null,
  );
}

// ---------- state ----------
let transcriptIndex = [];

// ---------- one guard for every panel that fills itself from the network ----------
// Assessing a variant starts a new assessment and takes the next number. Any request that has to
// wait for an answer records that number before it starts and checks on arrival that it is still
// the current one; if it is not, the user has moved on and the answer is dropped rather than drawn,
// cached or attached under the new variant's name. This covers answers that come back from a cache
// too — those arrive almost at once, but still after the user could have typed the next variant.
let assessmentSeq = 0;
function beginAssessment() { assessmentSeq += 1; return assessmentSeq; }
function currentAssessment() { return assessmentSeq; }
function superseded(seq) { return seq !== assessmentSeq; }

// =====================================================================
//  ASSESS VIEW
// =====================================================================
function buildAssess() {
  const root = $('#view-assess');
  root.innerHTML = '';

  const input = el('input', { type: 'text', id: 'variant-input', placeholder: 'e.g. NM_006015.6:c.461A>C', spellcheck: 'false' });
  const assessBtn = el('button', { class: 'btn', onclick: () => runAssess() }, 'Assess');
  input.addEventListener('keydown', (e) => { if (e.key === 'Enter') runAssess(); });


  const inputCard = el('div', { class: 'card' },
    el('h2', {}, 'Assess a coding variant'),
    el('div', { class: 'sub' }, 'Paste an HGVS coding change — transcript, then c. position and base substitution. Everything is computed offline.'),
    el('div', { class: 'row' }, el('div', { style: 'flex:1;min-width:280px' }, input), assessBtn,
      // Clear empties the panels that answer for a variant as well as the box, so nothing is left
      // asserting a verdict for a variant that is no longer on screen.
      el('button', { class: 'btn small secondary', onclick: () => { input.value = ''; $('#assess-result').innerHTML = ''; forgetVariant(); } }, 'Clear')),
  );

  const result = el('div', { id: 'assess-result' });   // stays empty until a variant is entered

  root.append(inputCard, knownLookupCard(), result);
}

// Every panel that answers for "the variant" is emptied together. Called when nothing on screen has
// resolved: raising the assessment number only stops an answer still in flight from being drawn, it
// does not take down what an earlier variant already put on the page. The assessment pill and the
// Splice, gnomAD and Assessment tabs carry no variant name, so a leftover verdict there reads as a
// verdict on whatever is in the box now.
// `keepAnnotation` is for the one case where the input did go off to be annotated — a recognised but
// untranslatable change, say. The pill is then answering for the input on screen and must stay.
function forgetVariant({ keepAnnotation = false } = {}) {
  lastVariant = null;
  assessState = { status: 'none', hgvs: null, data: null };
  contextGene = null;
  renderAssessment();
  if (keepAnnotation) return;
  latestGeneBe = { status: 'none' };
  updateGenebePill();
  renderGeneBe();
}

async function runAssess() {
  const raw = $('#variant-input').value.trim();
  const box = $('#assess-result');
  box.innerHTML = '';
  // Empty input is not an error and not a prompt — just leave the result area blank, and stop the
  // panels answering for the variant that was there before it was cleared.
  if (!raw) { forgetVariant(); return; }
  // From here on this is the variant on screen. Anything still in flight for the previous one is
  // superseded; its answer must not land on this variant's panels.
  beginAssessment();
  const parsed = parseHgvs(raw);

  // GeneBe annotation (desktop only) runs independently of our local resolution — GeneBe does its own
  // HGVS->genomic conversion. This drives the personal "GeneBe" tab; the browser build has no GeneBe.
  const willAnnotate = parsed.ok && parsed.transcript && !IS_WEB;
  if (willAnnotate) {
    if (autoAnnotateOn()) annotateGeneBe(parsed.raw);
    else { latestGeneBe = { status: 'off' }; updateGenebePill(); renderGeneBe(); }
  }

  if (!parsed.ok) {
    forgetVariant();
    box.append(el('div', { class: 'card' }, el('div', { class: 'err' }, parsed.message)));
    return;
  }

  let lookup = { doc: null, from: null, otherVersions: [] };
  if (parsed.transcript) {
    try { lookup = await window.api.getTranscript(parsed.transcript); } catch (e) { /* ignore */ }
  }

  if (!parsed.supported) {
    // "Not translated" is true and used to be the whole answer, which made every intronic and
    // untranslated-region variant a dead end — including the splice-site variants this tool is most
    // often reached for. The protein consequence genuinely cannot be computed, but the variant can
    // still be PLACED on the chromosome from the bundled exon table, and once placed it gets a
    // splice-region assessment offline and the same population, clinical and predictor lookups that a
    // coding variant gets. So the card is built rather than refused.
    forgetVariant({ keepAnnotation: willAnnotate });
    box.append(await buildNonCodingCard(parsed, lookup));
    return;
  }

  if (!parsed.transcript) {
    forgetVariant({ keepAnnotation: willAnnotate });
    box.append(el('div', { class: 'card' }, el('div', { class: 'warn' }, 'Add the transcript accession, e.g. NM_006015.6:c.461A>C, so the reference codon can be read.')));
    return;
  }

  if (!lookup.doc) {
    forgetVariant({ keepAnnotation: willAnnotate });
    box.append(renderNotBundled(parsed, lookup));
    return;
  }

  // Remember the variant so the gnomAD and Assessment tabs can look it up without re-typing.
  lastVariant = { hgvs: parsed.raw, gene: lookup.doc.gene, coord: null };
  // A new variant invalidates any prior assessment; refresh it if the Assessment tab is already open.
  assessState = { status: 'none', hgvs: null, data: null };
  if (document.querySelector('#view-assessment.active')) openAssessment();

  renderResolved(box, parsed, lookup.doc, lookup.from);
}

/**
 * The card for a variant whose protein consequence cannot be computed — intronic, or in either
 * untranslated region.
 *
 * Until this existed the app printed "Not translated" and stopped, which made a splice-site variant —
 * one of the commonest reasons to reach for this tool — a dead end. The protein effect genuinely
 * cannot be worked out from the coding sequence, and this card does not pretend otherwise. What it
 * does is place the variant on the chromosome from the bundled exon table, which is enough to say
 * where it sits relative to the splice sites, and enough for every coordinate-keyed lookup the tool
 * already does for coding variants.
 */
async function buildNonCodingCard(parsed, lookup) {
  const card = el('div', { class: 'card' });
  const gene = lookup.doc && lookup.doc.gene ? lookup.doc.gene : '';
  const kindWord = parsed.region === 'intronic' ? 'Intronic'
    : (parsed.region === 'utr5' ? "5′ untranslated region" : (parsed.region === 'utr3' ? "3′ untranslated region" : 'Non-coding'));

  card.append(el('div', { class: 'row', style: 'justify-content:space-between; align-items:flex-start' },
    el('div', {}, el('h2', {}, parsed.raw),
      el('div', { class: 'sub' }, [gene, parsed.transcript].filter(Boolean).join(' · '))),
    el('div', { class: 'consequence c-unknown' }, el('span', { class: 'dot' }), `${kindWord} — no amino-acid change`)));

  // The exon table is the same one the Fusions tab uses; load it once and share it.
  if (!fusionTables) { try { fusionTables = await window.api.getFusionTables(); } catch { fusionTables = {}; } }
  const exonTable = fusionTables && fusionTables.exons && fusionTables.exons.genes;
  const g = exonTable && gene ? exonTable[gene] : null;

  if (!g) {
    card.append(el('div', { class: 'note warn', style: 'margin-top:10px' },
      `${parsed.message} The exon table holds no record for ${gene || 'this gene'}, so the position on the chromosome cannot be worked out either.`));
    return card;
  }
  if (g.tx && parsed.transcript && g.tx !== parsed.transcript) {
    card.append(el('div', { class: 'note muted', style: 'margin-top:8px; font-size:12px' },
      `Exon coordinates are held for ${g.tx}; you asked about ${parsed.transcript}. Positions below are on ${g.tx}, and the two transcripts may number their bases differently.`));
  }

  const res = resolveCodingCoordinate(g, parsed);
  if (!res.ok) {
    card.append(el('div', { class: 'note warn', style: 'margin-top:10px' }, res.reason));
    return card;
  }

  const sr = spliceRegion(g, res);
  const tone = sr.kind === 'essential-site' ? 'err' : (sr.kind === 'splice-region' ? 'warn' : 'muted');
  card.append(el('div', { class: `note ${tone}`, style: 'margin-top:10px' }, sr.label));

  const facts = el('dl', { class: 'fus-facts', style: 'margin-top:10px' });
  const fact = (k, v) => { facts.append(el('dt', {}, k), el('dd', {}, v)); };
  fact('Position', `chr${res.chr}:${res.pos.toLocaleString('en-GB')} (GRCh38)`);
  if (res.ref && res.alt) {
    fact('Change on the chromosome', `${res.ref}>${res.alt}${res.complemented ? ` — complemented, because ${gene} is on the minus strand` : ''}`);
  }
  if (res.exon) fact('Exon', `${res.exon.number} of ${res.exon.of}`);
  if (res.intron) fact('Intron', `between exons ${res.intron.after} and ${res.intron.before}`);
  if (res.region === 'intronic') fact('Distance from the exon', `${Math.abs(res.offset)} base${Math.abs(res.offset) === 1 ? '' : 's'}`);
  card.append(facts);

  card.append(el('div', { class: 'note muted', style: 'margin-top:8px; font-size:11.5px' },
    `${parsed.message} ${SPLICE_CONVENTION}`));

  // Everything keyed on the coordinate now works for this variant exactly as it does for a coding one.
  lastVariant = {
    hgvs: parsed.raw, gene, appProtein: null,
    // Both spellings: the splice path reads `chr`, the assessment and gnomAD path reads `chrom`.
    coord: { chr: res.chr, chrom: res.chr, pos: res.pos, ref: res.ref, alt: res.alt, source: 'bundled exon table' },
    nonCoding: { region: res.region, splice: sr.kind, offset: res.offset },
  };
  card.append(el('div', { class: 'row', style: 'gap:8px; margin-top:12px; flex-wrap:wrap' },
    el('button', { class: 'btn small', onclick: () => switchTab('splice') }, 'Splice scores'),
    el('button', { class: 'btn small secondary', onclick: () => switchTab('gnomad') }, 'Population frequency'),
    el('button', { class: 'btn small secondary', onclick: () => switchTab('assessment') }, 'Clinical significance')));
  return card;
}

function renderNotBundled(parsed, lookup) {
  const card = el('div', { class: 'card' });
  card.append(el('h2', {}, parsed.raw), el('div', { class: 'sub' }, `Transcript ${parsed.transcript} is not in the offline store.`));
  if (lookup.otherVersions && lookup.otherVersions.length) {
    card.append(el('div', { class: 'note', style: 'margin-top:6px' }, `Available offline: ${lookup.otherVersions.join(', ')}. Versions differ in numbering — use the exact one from your report.`));
  }
  const fetchBtn = el('button', { class: 'btn small', style: 'margin-top:10px' }, 'Fetch this transcript from NCBI (sends only the accession)');
  fetchBtn.addEventListener('click', async () => {
    fetchBtn.disabled = true; fetchBtn.textContent = 'Fetching…';
    try {
      const res = await window.api.fetchTranscript(parsed.transcript);
      const box = $('#assess-result');
      box.innerHTML = '';
      // This variant is now the resolved one. Without this, `renderResolved` stamps its protein
      // change onto whatever `lastVariant` still held — the variant assessed before this one — and
      // the Assessment tab then shows that variant's gene beside this variant's protein change.
      lastVariant = { hgvs: parsed.raw, gene: res.doc && res.doc.gene, coord: null };
      assessState = { status: 'none', hgvs: null, data: null };
      renderResolved(box, parsed, res.doc, 'fetched');
    } catch (e) {
      fetchBtn.disabled = false; fetchBtn.textContent = 'Fetch failed — try again';
      card.append(el('div', { class: 'err', style: 'margin-top:8px' }, `Fetch failed: ${e.message}`));
    }
  });
  card.append(fetchBtn);
  return card;
}

function propRow(label, refVal, altVal, changed, help) {
  return [el('div', { class: 'k', title: help || null }, label), el('div', { class: `v${changed ? ' warn' : ''}` }, `${refVal}  →  ${altVal}`)];
}

// Say where a transcript came from in words, because the four sources are not equally
// annotated: a live fetch carries the sequence only — no domains, no strand, no 3'UTR.
function sourceLabel(from) {
  if (from === 'bundled') return 'source: curated offline set';
  if (from === 'reference') return 'source: MANE reference set';
  if (from === 'cache') return 'source: previously fetched, cached on this machine';
  if (from === 'fetched') return 'source: fetched from NCBI just now (sequence only)';
  return `source: ${from}`;
}

function renderResolved(box, parsed, doc, from) {
  if (parsed.kind !== 'substitution') { renderIndelResolved(box, parsed, doc, from); return; }
  const r = resolveSubstitution(doc.cds, { position: parsed.position, ref: parsed.ref, alt: parsed.alt });
  if (!r.ok) { box.append(el('div', { class: 'card' }, el('div', { class: 'err' }, r.message))); return; }
  // Record this app-resolved protein change so the Assessment tab can cross-check GeneBe against it.
  if (lastVariant) { lastVariant.appProtein = r.proteinShort; lastVariant.appTranscript = parsed.transcript; lastVariant.appResidue = r.codonNumber; }

  const refAaRec = aaByOne(r.refAa);
  const altAaRec = aaByOne(r.altAa);

  const badge = el('div', { class: `consequence c-${r.consequence}` }, el('span', { class: 'dot' }), consequenceLabel(r.consequence));

  // Arm the Gene context tab for this gene, and offer a direct way in.
  contextGene = doc.gene || null;
  const contextBtn = el('button', { class: 'btn small', style: 'margin-top:8px' }, 'Gene context — pathway & reported role');
  contextBtn.addEventListener('click', () => { setGeneContext(doc.gene); switchTab('context'); });
  const gnomadBtn = el('button', { class: 'btn small secondary', style: 'margin-top:8px; margin-left:8px' }, 'gnomAD frequency');
  gnomadBtn.addEventListener('click', () => { lastGnomad = null; switchTab('gnomad'); });

  const header = el('div', { class: 'row', style: 'justify-content:space-between; align-items:flex-start; gap:16px' },
    el('div', {}, el('h2', {}, `${doc.gene || ''} ${r.proteinShort}`.trim()),
      el('div', { class: 'sub' }, `${parsed.transcript} · ${doc.proteinName || ''} · ${sourceLabel(from)}`),
      el('div', {}, contextBtn, gnomadBtn)),
    badge);

  // Block A — amino-acid change
  const blockA = el('div', { class: 'ib' }, el('div', { class: 'ib-h' }, 'Amino-acid change'),
    el('div', { class: 'aa-change' }, aaTile(r.refAa, 'Reference'), el('div', { class: 'arrow' }, '→'), aaTile(r.altAa, 'Altered')));

  // Block B — codon + coordinates
  const refCodonEl = el('div', { class: 'codon' }, r.refCodon.split('').map((b, i) => baseCell(b, { changed: i === r.posInCodon })));
  const altCodonEl = el('div', { class: 'codon' }, r.altCodon.split('').map((b, i) => baseCell(b, { changed: i === r.posInCodon })));
  const blockB = el('div', { class: 'ib' }, el('div', { class: 'ib-h' }, 'Codon & position'),
    el('div', { class: 'codon-change' },
      el('div', {}, el('div', { class: 'muted', style: 'font-size:12px;margin-bottom:4px' }, `Codon ${r.codonNumber}`), refCodonEl),
      el('div', { class: 'arrow' }, '→'),
      el('div', {}, el('div', { class: 'muted', style: 'font-size:12px;margin-bottom:4px' }, 'Altered'), altCodonEl)),
    el('div', { class: 'kv', style: 'margin-top:12px' },
      el('div', { class: 'k' }, 'Protein (3-letter)'), el('div', { class: 'v' }, r.proteinHgvs),
      el('div', { class: 'k' }, 'Protein (1-letter)'), el('div', { class: 'v' }, r.proteinShort),
      el('div', { class: 'k' }, 'Coding position'), el('div', { class: 'v' }, `c.${r.position} · codon ${r.codonNumber}, base ${r.posInCodon + 1}/3`),
      el('div', { class: 'k' }, 'Codon'), el('div', { class: 'v' }, `${r.refCodon} → ${r.altCodon}`)));

  // Block C — side-chain property comparison
  const props = el('div', { class: 'kv' });
  if (refAaRec && altAaRec && r.refAa !== '*' && r.altAa !== '*') {
    props.append(
      ...propRow('Polarity', refAaRec.polarity, altAaRec.polarity, refAaRec.polarity !== altAaRec.polarity),
      ...propRow('Charge (pH 7.4)', fmtCharge(refAaRec.charge), fmtCharge(altAaRec.charge), refAaRec.charge !== altAaRec.charge),
      // Which residues count as aromatic is a convention, like the polarity grouping and the mass —
      // histidine is the one that differs between tables — so it is named here too.
      ...propRow('Aromatic', refAaRec.aromatic ? 'yes' : 'no', altAaRec.aromatic ? 'yes' : 'no', refAaRec.aromatic !== altAaRec.aromatic, AROMATICITY_NOTE),
      ...propRow('Hydropathy', refAaRec.hydropathy, altAaRec.hydropathy, true),
      // "Mass (Da)" alone read as the free amino acid's mass; these are residue masses, 18 Da lighter.
      ...propRow('Residue mass (Da)', refAaRec.mass.toFixed(1), altAaRec.mass.toFixed(1), true, RESIDUE_MASS_HELP));
  } else {
    props.append(el('div', { class: 'note' }, 'Shown for amino-acid substitutions.'));
  }
  const blockC = el('div', { class: 'ib' }, el('div', { class: 'ib-h' }, 'Nature of the change'), props);

  const refCheck = r.refMatches
    ? el('div', { class: 'note', style: 'color:#8fe0b0; margin-top:8px' }, `Reference base check: c.${r.position} is ${r.refBaseInTranscript} in ${doc.accession} — matches the stated ${parsed.ref}.`)
    : el('div', { class: 'note err', style: 'margin-top:8px' }, `Reference base check: c.${r.position} is ${r.refBaseInTranscript} in ${doc.accession}, but the variant states ${parsed.ref}. Check the transcript version and strand before relying on this.`);

  const track = el('div', {},
    el('div', { class: 'ib-h', style: 'margin-bottom:6px' }, `Bases around the change — coding sequence c.${r.window.startCds}–${r.window.endCds} (5′→3′)`),
    buildTrack(r, doc));

  box.append(el('div', { class: 'card' },
    header,
    el('div', { class: 'inforow' }, blockA, blockB, blockC),
    refCheck,
    el('div', { class: 'divider' }),
    buildLollipop(doc, [{ residue: r.codonNumber ?? r.firstAffected, label: r.proteinShort, consequence: r.consequence }]),
    el('div', { class: 'divider' }),
    track));

  box.append(buildStructureCard(doc, r));
}

// ---------- indel / frameshift result ----------
function buildIndelTrack(r, doc) {
  const { ref, alt, eventIndex, eventLength } = r.window;
  const wrap = el('div', { class: 'track-wrap' });
  const track = el('div', { class: 'track' });
  const flank = (t) => el('div', { class: 'tflank' }, t);

  let altStart; let altLen;
  if (r.kind === 'deletion') { altStart = eventIndex; altLen = 0; }
  else if (r.kind === 'duplication') { altStart = eventIndex + eventLength; altLen = eventLength; }
  else if (r.kind === 'insertion') { altStart = eventIndex + 1; altLen = r.insertedNt.length; }
  else { altStart = eventIndex; altLen = r.insertedNt.length; }

  const makeRow = (label, seq, hlStart, hlLen) => {
    const row = el('div', { class: 'track-row' });
    row.append(el('div', { class: 'rowlabel' }, label), flank('5′'));
    for (let i = 0; i < seq.length; i += 1) {
      const inHl = hlLen > 0 && i >= hlStart && i < hlStart + hlLen;
      row.append(el('div', { class: `tbase ${seq[i]}${inHl ? ' hl' : ''}` }, seq[i]));
    }
    row.append(flank('3′'));
    return row;
  };

  track.append(makeRow('Reference', ref, eventIndex, r.kind === 'insertion' ? 0 : eventLength));
  track.append(makeRow('Altered', alt, altStart, altLen));
  wrap.append(track);

  const what = r.kind === 'deletion' ? `The outlined reference bases (${r.deletedNt}) are deleted.`
    : r.kind === 'duplication' ? `The outlined bases in the altered row are the duplicated copy (${r.insertedNt}).`
      : r.kind === 'insertion' ? `The outlined bases in the altered row are inserted (${r.insertedNt}).`
        : `The outlined reference bases (${r.deletedNt}) are replaced by ${r.insertedNt}.`;

  const strand = doc && doc.strand; const chr = doc && doc.chromosome;
  let orient = 'Shown on the coding (sense) strand, 5′ → 3′ in the direction of translation.';
  if (strand === 'minus') orient += ` ${doc.gene} is on the minus strand${chr ? ` of chromosome ${chr}` : ''}, so these coding bases are the reverse complement of the genomic (+) strand.`;
  else if (strand === 'plus') orient += ` ${doc.gene} is on the plus strand${chr ? ` of chromosome ${chr}` : ''}.`;

  return el('div', {}, wrap,
    el('div', { class: 'note', style: 'margin-top:8px' }, what),
    el('div', { class: 'note', style: 'margin-top:6px' }, orient));
}

function renderIndelResolved(box, parsed, doc, from) {
  const r = resolveIndel(doc.cds, parsed, { utr3: doc.utr3 });
  if (!r.ok) { box.append(el('div', { class: 'card' }, el('div', { class: 'err' }, r.message))); return; }
  if (lastVariant) { lastVariant.appProtein = r.proteinShort; lastVariant.appTranscript = parsed.transcript; lastVariant.appResidue = r.codonNumber ?? r.firstAffected; }

  const badge = el('div', { class: `consequence c-${r.consequence}` }, el('span', { class: 'dot' }), consequenceLabel(r.consequence));

  // Arm the Gene context tab for this gene, and offer a direct way in.
  contextGene = doc.gene || null;
  const contextBtn = el('button', { class: 'btn small', style: 'margin-top:8px' }, 'Gene context — pathway & reported role');
  contextBtn.addEventListener('click', () => { setGeneContext(doc.gene); switchTab('context'); });
  const gnomadBtn = el('button', { class: 'btn small secondary', style: 'margin-top:8px; margin-left:8px' }, 'gnomAD frequency');
  gnomadBtn.addEventListener('click', () => { lastGnomad = null; switchTab('gnomad'); });

  const header = el('div', { class: 'row', style: 'justify-content:space-between; align-items:flex-start; gap:16px' },
    el('div', {}, el('h2', {}, `${doc.gene || ''} ${r.proteinShort}`.trim()),
      el('div', { class: 'sub' }, `${parsed.transcript} · ${doc.proteinName || ''} · ${sourceLabel(from)}`),
      el('div', {}, contextBtn, gnomadBtn)),
    badge);

  const eventLabel = r.end !== r.start ? `c.${r.start}_${r.end}` : `c.${r.start}`;
  const kindWord = { deletion: 'deletion', duplication: 'duplication', insertion: 'insertion', delins: 'deletion-insertion' }[r.kind] || r.kind;
  const basesLabel = r.kind === 'deletion' ? 'Bases removed' : r.kind === 'insertion' ? 'Bases inserted' : r.kind === 'duplication' ? 'Bases duplicated' : 'Removed → inserted';
  const basesValue = r.kind === 'delins' ? `${r.deletedNt} → ${r.insertedNt}` : (r.deletedNt || r.insertedNt || '—');

  const blockA = el('div', { class: 'ib' }, el('div', { class: 'ib-h' }, 'Event'),
    el('div', { class: 'kv' },
      el('div', { class: 'k' }, 'Type'), el('div', { class: 'v' }, kindWord),
      el('div', { class: 'k' }, 'Coding position'), el('div', { class: 'v' }, eventLabel),
      el('div', { class: 'k' }, basesLabel), el('div', { class: 'v' }, basesValue),
      el('div', { class: 'k' }, 'Net change'), el('div', { class: 'v' }, `${r.netNt > 0 ? '+' : ''}${r.netNt} nt · ${r.frameshift ? 'frameshift' : 'in frame'}`)));

  const blockB = el('div', { class: 'ib' }, el('div', { class: 'ib-h' }, 'Protein consequence'),
    el('div', { class: 'kv' },
      el('div', { class: 'k' }, 'Protein (3-letter)'), el('div', { class: 'v' }, r.proteinHgvs),
      el('div', { class: 'k' }, 'Protein (1-letter)'), el('div', { class: 'v' }, r.proteinShort),
      el('div', { class: 'k' }, 'First affected'), el('div', { class: 'v' }, `${r.refAaThree || '?'}${r.firstAffected}`),
      el('div', { class: 'k' }, 'Protein length'), el('div', { class: 'v' }, `${r.proteinLength} aa → ${r.mutantProteinLength} aa`)));

  const fsCount = String(r.proteinShort || '').split('fs*')[1];
  const blockC = el('div', { class: 'ib' }, el('div', { class: 'ib-h' }, 'What this means'),
    el('div', { class: 'note' }, r.frameshift
      ? `The reading frame shifts from residue ${r.firstAffected}${fsCount && fsCount !== '?' ? `, and the new frame reaches a stop ${fsCount} residues later` : ''} — the protein is truncated from ${r.proteinLength} to ${r.mutantProteinLength} residues.`
      : 'The reading frame is preserved; only the listed residues are changed.'),
    r.stopNotFound ? el('div', { class: 'note warn', style: 'margin-top:6px' }, 'No new stop was found within the bundled 3′ untranslated region, so the stop position is shown as “?”.') : null,
    r.refMatches === false ? el('div', { class: 'note err', style: 'margin-top:6px' }, 'The stated reference bases do not match this transcript — check the transcript version.') : null,
    r.refMatches === true ? el('div', { class: 'note', style: 'margin-top:6px; color:#8fe0b0' }, 'The stated reference bases match the transcript.') : null);

  const track = el('div', {},
    el('div', { class: 'ib-h', style: 'margin-bottom:6px' }, `Bases around the change — coding sequence c.${r.window.startCds}–${r.window.endCds} (5′→3′)`),
    buildIndelTrack(r, doc));

  box.append(el('div', { class: 'card' },
    header,
    el('div', { class: 'inforow' }, blockA, blockB, blockC),
    el('div', { class: 'divider' }),
    buildLollipop(doc, [{ residue: r.codonNumber ?? r.firstAffected, label: r.proteinShort, consequence: r.consequence }]),
    el('div', { class: 'divider' }),
    track));

  box.append(buildStructureCard(doc, r));
}

// ---------- domain-level lollipop plot ----------
function consequenceColour(c) {
  return ({
    missense: '#e0a54e', nonsense: '#e05561', synonymous: '#46b877',
    'stop-loss': '#a86cff', 'start-loss': '#a86cff', 'stop-retained': '#46b877',
    frameshift: '#e05561',
    'inframe-deletion': '#a86cff', 'inframe-insertion': '#a86cff',
    'inframe-duplication': '#a86cff', 'inframe-delins': '#a86cff',
  })[c] || '#8a8f99';
}
function normDesc(d) { return String(d || '').replace(/\s*\d+$/, '').trim(); }
// Width of one character of the 11 px bold label type, in the diagram's own units. One constant, so
// the width used to place a label and the width used to decide whether it fits cannot disagree.
const LABEL_CHAR_PX = 6.4;
function labelWidthPx(text) { return String(text || '').length * LABEL_CHAR_PX; }
function truncateLabel(text, widthPx) {
  const max = Math.max(2, Math.floor((widthPx - 8) / LABEL_CHAR_PX));
  return text.length > max ? `${text.slice(0, max - 1)}…` : text;
}
function niceStep(L) {
  const target = L / 6;
  for (const s of [25, 50, 100, 200, 250, 500, 1000, 2000, 5000]) if (s >= target) return s;
  return 5000;
}

/**
 * Where one residue sits relative to the annotated domains, in words.
 *
 * `features` must be the domains as annotated, sorted by start — never the merged bars drawn above.
 * Neighbouring repeats of one domain are merged into a single bar for the picture, and that bar
 * spans the stretches between them, which are annotated as no domain at all. Reading membership off
 * the bar therefore places a residue in a domain the data does not put it in: XPO1 has HEAT repeats
 * at 602–639 and 775–813 and nothing in between, but the two are drawn as one bar, so residue 700
 * was reported as lying in "HEAT 3 ×5 (354–813)".
 */
function residueDomainSentence(residue, features) {
  const inside = features.filter((d) => residue >= d.start && residue <= d.end);
  if (inside.length) {
    return `Residue ${residue} lies in ${inside.map((d) => `${d.label} (${d.start}–${d.end})`).join('; ')}.`;
  }
  if (!features.length) return `Residue ${residue} lies outside the annotated domains.`;
  const before = features.filter((d) => d.end < residue).sort((a, b) => a.end - b.end).pop();
  const after = features.filter((d) => d.start > residue)[0];
  const name = (d) => `${d.label} (${d.start}–${d.end})`;
  if (before && after) {
    return `Residue ${residue} lies between ${name(before)} and ${name(after)}, in neither.`;
  }
  if (after) return `Residue ${residue} lies before the first annotated domain, ${name(after)}.`;
  return `Residue ${residue} lies after the last annotated domain, ${name(before)}.`;
}

/**
 * A marker position past the end of the protein, in words.
 *
 * A substitution in the stop codon resolves to codon L+1 of an L-residue protein. That is a position
 * in the coding sequence, not a residue, so it is named as what it is instead of being called a
 * residue and drawn on top of the last real one.
 */
function beyondProteinSentence(position, proteinLength) {
  return position === proteinLength + 1
    ? `Position ${position} is the stop codon, one past the last residue (${proteinLength}) — not a residue of the protein.`
    : `Position ${position} lies past the end of this protein (${proteinLength} residues).`;
}

/**
 * The coordinate frame the drawn positions are in, in words.
 *
 * Saying "MANE protein numbering" everywhere was wrong on the curated transcripts that are
 * deliberately not MANE: MYD88 NM_002468.4 is 309 residues and carries a domain ending at 306, which
 * cannot be a position in the 296-residue MANE protein at all. The numbering frame is always the
 * displayed transcript's own, so that is what is named; MANE is only mentioned when the record says
 * the transcript is MANE.
 */
function numberingFrame(doc) {
  if (!doc || !doc.accession) return 'positions in this transcript’s own protein numbering';
  return `positions in ${doc.accession}${doc.maneStatus ? ` (${doc.maneStatus})` : ''} protein numbering`;
}

/**
 * Domain lollipop. `markers` is an array of { residue, label, consequence } so several
 * variants in the same protein can share one diagram (labels are staggered to avoid
 * overlapping). Pass a single-element array for the normal one-variant view.
 */
function buildLollipop(doc, markers, opts = {}) {
  const marks = (Array.isArray(markers) ? markers : [markers]).filter((m) => m && Number.isInteger(m.residue));
  const L = doc.proteinLength || (doc.uniprot && doc.uniprot.length) || 0;
  const wrap = el('div', {});
  if (opts.heading !== false) {
    // Only claim MANE when the record actually says so. Six of the curated transcripts are
    // deliberately not MANE (MYD88 NM_002468.4, ASXL1 .5, SRSF2, SGK1, TNFAIP3, BIRC3), and a
    // live-fetched transcript is unknown, so an unconditional "(MANE)" here was simply wrong.
    const status = doc.maneStatus ? ` (${doc.maneStatus})` : '';
    wrap.append(el('div', { class: 'ib-h', style: 'margin-bottom:6px' },
      `Position in the protein — ${doc.gene || ''} ${L} aa${status}`));
  }
  if (!L) { wrap.append(el('div', { class: 'note' }, 'Protein length unknown; cannot draw the domain map.')); return wrap; }

  // Merge consecutive features that share a normalised label (e.g. EGF-like repeats).
  const raw = (doc.domains || []).filter((d) => Number.isInteger(d.start) && Number.isInteger(d.end) && d.end >= d.start).slice().sort((a, b) => a.start - b.start);
  const merged = [];
  for (const d of raw) {
    const nd = normDesc(d.label);
    const last = merged[merged.length - 1];
    if (last && normDesc(last.label) === nd && d.start <= last.end + Math.max(8, last.end - last.start)) {
      last.end = Math.max(last.end, d.end); last.count += 1;
    } else {
      merged.push({ label: d.label, ndesc: nd, start: d.start, end: d.end, count: 1, source: d.source, accession: d.accession });
    }
  }
  const colourMap = new Map(); let ci = 0;
  for (const d of merged) if (!colourMap.has(d.ndesc)) colourMap.set(d.ndesc, DOMAIN_PALETTE[ci++ % DOMAIN_PALETTE.length]);

  // Lane-pack overlapping domains.
  const lanes = [];
  for (const d of merged) {
    let placed = false;
    for (const lane of lanes) { if (d.start > lane[lane.length - 1].end) { lane.push(d); placed = true; break; } }
    if (!placed) lanes.push([d]);
  }
  const laneCount = Math.max(1, lanes.length);

  const NS = 'http://www.w3.org/2000/svg';
  const svg = (t, a, ...k) => { const n = document.createElementNS(NS, t); for (const [key, v] of Object.entries(a || {})) if (v != null) n.setAttribute(key, v); for (const kid of k.flat()) { if (kid == null) continue; n.append(kid.nodeType ? kid : document.createTextNode(String(kid))); } return n; };
  const VW = 1000, mL = 14, mR = 14, W = VW - mL - mR;
  const marksAt = marks.map((m) => ({ ...m })).sort((a, b) => a.residue - b.residue);
  // The drawing runs to the last residue, or past it when a marker does. A substitution in the stop
  // codon resolves to codon L+1; clamping that onto L put its circle exactly on the last residue's
  // pixel while the caption said position L+1 — the picture claimed the wrong residue.
  const XMAX = Math.max(L, ...marksAt.map((m) => m.residue));
  const x = (res) => mL + (Math.min(Math.max(res, 1), XMAX) - 1) / Math.max(1, XMAX - 1) * W;
  const circleR = 9, laneH = 20, laneGap = 6, markerRowH = 22, stemBase = 34, LABEL_GAP = 6;

  // Stagger markers whose labels would overlap: each gets a level (taller stem). The test is the
  // label's own width, not a fixed pixel gap — a frameshift label is nearly twice as wide as a
  // substitution one, so a single gap let two long labels share a level and overlap. Each label is
  // also held inside the drawing, so a label near either end is no longer cut off by the edge.
  const lastRightByLevel = [];
  for (const m of marksAt) {
    m._x = x(m.residue);
    m._label = truncateLabel(String(m.label || ''), W);
    const half = labelWidthPx(m._label) / 2;
    m._labelX = Math.min(Math.max(m._x, half + 2), VW - half - 2);
    let lvl = 0;
    while (lastRightByLevel[lvl] != null && m._labelX - half < lastRightByLevel[lvl] + LABEL_GAP) lvl += 1;
    lastRightByLevel[lvl] = m._labelX + half;
    m._level = lvl;
  }
  const maxLevel = marksAt.length ? Math.max(...marksAt.map((m) => m._level)) : 0;
  const backboneY = stemBase + maxLevel * markerRowH + 26;
  const domainsTop = backboneY - 9;
  const axisY = domainsTop + laneCount * (laneH + laneGap) + 4;
  const VH = axisY + 24;
  const markerY = (m) => backboneY - stemBase - m._level * markerRowH;

  const parts = [];
  for (const m of marksAt) {
    parts.push(svg('line', { x1: m._x, y1: markerY(m), x2: m._x, y2: backboneY, stroke: '#c9d3e6', 'stroke-width': 2 }));
  }
  parts.push(svg('rect', { x: mL, y: backboneY - 2, width: W, height: 4, rx: 2, fill: '#38445e' }));

  lanes.forEach((lane, li) => {
    const y = domainsTop + li * (laneH + laneGap);
    for (const d of lane) {
      const dx = x(d.start), dw = Math.max(3, x(d.end) - x(d.start));
      const rect = svg('rect', { x: dx, y, width: dw, height: laneH, rx: 4, fill: colourMap.get(d.ndesc), opacity: 0.92 });
      rect.append(svg('title', {}, `${d.label}${d.count > 1 ? ` ×${d.count}` : ''} (${d.start}–${d.end}) — ${d.source}`));
      parts.push(rect);
      if (dw > 46) parts.push(svg('text', { x: dx + dw / 2, y: y + laneH / 2 + 4, 'text-anchor': 'middle', fill: '#0c0f16', 'font-size': 11, 'font-weight': 700 }, truncateLabel(d.ndesc, dw)));
    }
  });

  // axis
  parts.push(svg('line', { x1: mL, y1: axisY, x2: mL + W, y2: axisY, stroke: '#38445e', 'stroke-width': 1 }));
  const ticks = [1]; const step = niceStep(L); for (let t = step; t < L - step / 2; t += step) ticks.push(t); ticks.push(L);
  for (const t of ticks) {
    parts.push(svg('line', { x1: x(t), y1: axisY, x2: x(t), y2: axisY + 4, stroke: '#6b7891', 'stroke-width': 1 }));
    parts.push(svg('text', { x: x(t), y: axisY + 16, 'text-anchor': 'middle', fill: '#6b7891', 'font-size': 10 }, t));
  }

  // Markers drawn last so they sit above the domains, and every circle before every label: drawing
  // each marker whole meant a promoted circle landed on top of the label already written next to it.
  for (const m of marksAt) {
    parts.push(svg('circle', { cx: m._x, cy: markerY(m), r: circleR, fill: consequenceColour(m.consequence), stroke: '#0b0f18', 'stroke-width': 1.5 }));
  }
  for (const m of marksAt) {
    parts.push(svg('text', { x: m._labelX, y: markerY(m) - 14, 'text-anchor': 'middle', fill: '#e7ecf6', 'font-size': 11, 'font-weight': 700 }, m._label));
  }

  wrap.append(el('div', { class: 'lolli-wrap' }, svg('svg', { viewBox: `0 0 ${VW} ${VH}`, width: '100%', preserveAspectRatio: 'xMidYMid meet', style: 'display:block' }, ...parts)));

  if (opts.caption !== false) {
    let text;
    if (!merged.length) {
      // Distinguish "this protein genuinely has no annotated domains" from "we have no domain
      // data for this transcript" — they mean very different things when writing a report.
      const src = doc.domainSource || '';
      let why;
      if (/could not be reconciled/.test(src)) why = 'Domains are known for this protein but could not be placed safely in this transcript’s numbering, so none are drawn.';
      else if (/no reviewed UniProt entry/.test(src)) why = 'No reviewed UniProt entry matched this transcript, so no domains are available.';
      else if (/no UniProt domain features/.test(src)) why = 'This protein has no domain features in UniProt.';
      else why = 'This transcript is not in the offline reference set, so no domain map is bundled for it.';
      text = `${why} Markers still show their residue of ${L}.`;
    } else {
      text = marksAt.map((m) => (m.residue > L
        ? beyondProteinSentence(m.residue, L)
        : residueDomainSentence(m.residue, raw))).join(' ');
    }
    wrap.append(el('div', { class: 'note', style: 'margin-top:6px' }, text));
  }

  // labelled legend — every domain named, with its range and source
  if (merged.length) {
    const legend = el('div', { class: 'domain-legend' });
    for (const d of merged) {
      legend.append(el('div', { class: 'dleg', title: `${d.source}${d.accession ? ' · ' + d.accession : ''}` },
        el('span', { class: 'dsw', style: `background:${colourMap.get(d.ndesc)}` }),
        el('span', { class: 'dlabel' }, `${d.label}${d.count > 1 ? ` ×${d.count}` : ''}`),
        el('span', { class: 'drange' }, `${d.start}–${d.end}`)));
    }
    wrap.append(legend);
    if (opts.source !== false) {
      wrap.append(el('div', { class: 'note muted', style: 'margin-top:6px;font-size:11.5px' },
        `Domains from ${doc.domainSource || 'InterPro/UniProt'} · ${numberingFrame(doc)}.`));
    }
  }
  return wrap;
}

function fmtCharge(c) { if (c === 1) return '+1'; if (c === -1) return '−1'; if (c === 0.1) return '≈0 (+)'; return '0'; }

function buildTrack(r, doc) {
  const { ref, alt, changeIndex, startCds } = r.window;
  const wrap = el('div', { class: 'track-wrap' });
  const track = el('div', { class: 'track' });
  const flank = (t) => el('div', { class: 'tflank' }, t);

  const makeRow = (label, seq) => {
    const row = el('div', { class: 'track-row' });
    row.append(el('div', { class: 'rowlabel' }, label));
    row.append(flank('5′'));
    for (let i = 0; i < seq.length; i++) {
      const cdsPos = startCds + i;
      const cell = el('div', { class: `tbase ${seq[i]}${i === changeIndex ? ' hl' : ''}` }, seq[i]);
      if (i > 0 && (cdsPos - 1) % 3 === 0) cell.style.marginLeft = '7px'; // codon boundary
      row.append(cell);
    }
    row.append(flank('3′'));
    return row;
  };

  // ruler
  const ruler = el('div', { class: 'track-row' });
  ruler.append(el('div', { class: 'rowlabel' }, ''));
  ruler.append(flank(''));
  for (let i = 0; i < ref.length; i++) {
    const cdsPos = startCds + i;
    const tick = el('div', { class: 'rtick' }, (i === changeIndex || cdsPos % 5 === 0) ? cdsPos : '');
    if (i > 0 && (cdsPos - 1) % 3 === 0) tick.style.marginLeft = '7px';
    ruler.append(tick);
  }
  ruler.append(flank(''));

  track.append(makeRow('Reference', ref), makeRow('Altered', alt), ruler);
  wrap.append(track);

  // base counts of the reference window
  const counts = { A: 0, C: 0, G: 0, T: 0 };
  for (const b of ref) if (counts[b] != null) counts[b]++;
  const countChips = el('div', { class: 'count-chips' },
    ['A', 'C', 'G', 'T'].map((b) => el('div', { class: 'count-chip' },
      el('div', { class: 'sw', style: `background:var(--base-${b})` }), `${b}: ${counts[b]}`)));

  // gene orientation
  const strand = doc && doc.strand;
  const chr = doc && doc.chromosome;
  let orient = 'Shown on the coding (sense) strand, 5′ → 3′ in the direction of translation — the same orientation HGVS c. numbering uses.';
  if (strand === 'minus') orient += ` ${doc.gene} is on the minus strand${chr ? ` of chromosome ${chr}` : ''}, so these coding bases are the reverse complement of the genomic (+) strand.`;
  else if (strand === 'plus') orient += ` ${doc.gene} is on the plus strand${chr ? ` of chromosome ${chr}` : ''}, so these coding bases match the genomic (+) strand.`;

  return el('div', {}, wrap, countChips, el('div', { class: 'note', style: 'margin-top:8px' }, orient));
}

// =====================================================================
//  PROTEIN STRUCTURE
// =====================================================================
let viewer = null;

function buildStructureCard(doc, r) {
  const card = el('div', { class: 'card' });
  const up = doc.uniprot;
  // Codon number for substitutions; first affected residue for indels/frameshifts.
  const resNum = r.codonNumber ?? r.firstAffected;
  // AlphaFold/RCSB are numbered by UniProt; map the MANE residue into that frame.
  const uni = maneToUniprot(doc.uniprotMap, resNum);
  const afResi = uni.pos;
  // maneToUniprot returns null where the two proteins genuinely cannot be lined up (KRAS 151-186,
  // EZH2 298-303, MYD88 2-14). Nothing is then highlighted on the model, so the card must neither
  // promise a highlight nor print the missing number.
  const unmapped = Boolean(up) && afResi == null;
  const head = el('div', { class: 'row', style: 'justify-content:space-between' },
    el('div', {}, el('h2', {}, 'Protein structure'),
      el('div', { class: 'sub' }, unmapped
        ? `${r.refAaThree || ''}${resNum} cannot be placed in this model’s numbering, so no residue is highlighted on it. A predicted model still loads automatically when available (public model by UniProt id, then offline).`
        : `Highlighting ${r.refAaThree || ''}${resNum}. A predicted model loads automatically when available (public model by UniProt id, then offline).`)),
    el('div', { class: 'chips' },
      up ? el('button', { class: 'btn small', onclick: () => loadStructure('alphafold', up.accession, afResi, msg) }, `AlphaFold model (${up.accession})`) : null,
      el('button', { class: 'btn small secondary', onclick: () => loadStructure('local', null, resNum, msg) }, 'Open local file')));

  const pdbRow = el('div', { class: 'row', style: 'margin-top:10px' },
    el('input', { type: 'text', id: 'pdbid', placeholder: 'PDB id, e.g. 1TUP', style: 'max-width:180px' }),
    el('button', { class: 'btn small secondary', onclick: () => loadStructure('rcsb', $('#pdbid').value, afResi, msg) }, 'Load PDB entry'),
    el('button', { class: 'btn small secondary', onclick: () => { if (viewer) { viewer.zoomTo(); viewer.render(); } } }, 'Reset view'));

  const viewerEl = el('div', { id: 'viewer' }, el('div', { class: 'hint' },
    up && up.accession ? 'Loading a predicted structure…' : 'No model available to load automatically. Enter a PDB id or open a local .pdb file.'));
  const msg = el('div', { class: 'note', style: 'margin-top:10px' });
  if (unmapped) msg.append(`This model is numbered by UniProt ${up.accession}, and MANE residue ${resNum} cannot be placed in that numbering${uni.reason ? `: ${uni.reason}` : ''}. No residue is highlighted on the model.`);
  else if (up && afResi !== resNum) msg.append(`This model is numbered by UniProt ${up.accession}: MANE residue ${resNum} corresponds to UniProt residue ${afResi}${uni.exact ? '' : ' (approximate)'}, which is what the highlight uses.`);
  else if (up) msg.append(`Numbering matches between MANE and UniProt ${up.accession}.`);

  card.append(head, pdbRow, viewerEl, msg);
  // Auto-load a predicted structure once the card is in the DOM.
  if (up && up.accession) setTimeout(() => loadStructure('alphafold', up.accession, afResi, msg, true), 30);
  return card;
}

async function loadStructure(kind, id, resi, msgEl, auto = false) {
  const viewerEl = $('#viewer');
  if (!viewerEl) return;
  if ((kind === 'alphafold' || kind === 'rcsb') && !String(id || '').trim()) {
    if (!auto) { viewerEl.innerHTML = ''; viewerEl.append(el('div', { class: 'hint err' }, 'Enter an identifier first.')); }
    return;
  }
  if (msgEl) for (const n of [...msgEl.querySelectorAll('.load-status')]) n.remove();
  viewerEl.innerHTML = ''; viewerEl.append(el('div', { class: 'hint' }, 'Loading…'));
  try {
    let res;
    if (kind === 'alphafold') res = await window.api.structureAlphaFold(id);
    else if (kind === 'rcsb') res = await window.api.structureRcsb(id);
    else res = await window.api.structureOpenLocal();
    if (!res) { viewerEl.innerHTML = ''; viewerEl.append(el('div', { class: 'hint' }, 'Cancelled.')); return; }
    const found = renderStructure(viewerEl, res, resi);
    const parts = [`Loaded: ${res.source}${res.predicted ? ' — a prediction, not an experimental structure.' : ''}`];
    if (Number.isInteger(resi) && found === 0) parts.push(`Residue ${resi} is not present in this model — a large protein may be split into fragments, or the numbering differs.`);
    if (msgEl) msgEl.append(el('div', { class: `note load-status${(Number.isInteger(resi) && found === 0) ? ' warn' : ''}`, style: 'margin-top:6px' }, parts.join(' ')));
  } catch (e) {
    viewerEl.innerHTML = '';
    if (auto) viewerEl.append(el('div', { class: 'hint' }, 'No structure loaded automatically (offline, or no model for this protein). Use the buttons above to load one.'));
    else viewerEl.append(el('div', { class: 'hint err' }, `Could not load structure: ${e.message}`));
  }
}

function renderStructure(container, res, resi) {
  container.innerHTML = '';
  const $3Dmol = window.$3Dmol;
  if (!$3Dmol) { container.append(el('div', { class: 'hint err' }, '3D library not available.')); return 0; }
  if (viewer) { try { viewer.clear(); } catch { /* */ } }
  viewer = $3Dmol.createViewer(container, { backgroundColor: '#0b0f18' });
  viewer.addModel(res.text, res.format || 'pdb');
  viewer.setStyle({}, { cartoon: { color: 'spectrum', opacity: 0.9 } });
  let found = 0;
  if (Number.isInteger(resi)) {
    try { const sel = viewer.selectedAtoms({ resi }); found = sel ? sel.length : 0; } catch { found = 0; }
    if (found > 0) {
      viewer.addStyle({ resi }, { stick: { radius: 0.35, colorscheme: 'whiteCarbon' } });
      viewer.addStyle({ resi }, { sphere: { scale: 0.5, color: '#ff4d6d' } });
      viewer.addResLabels({ resi }, { fontColor: 'white', backgroundColor: '#ff4d6d', fontSize: 12, showBackground: true });
      viewer.zoomTo({ resi });
    } else {
      viewer.zoomTo();
    }
  } else {
    viewer.zoomTo();
  }
  viewer.zoom(0.9);
  viewer.render();
  return found;
}

// =====================================================================
//  AMINO ACID REFERENCE VIEW
// =====================================================================
function polTag(polarity) {
  const c = POLARITY_COLOUR[polarity] || '#8a8f99';
  return el('span', { class: 'tag', style: `background:${c}` }, polarity);
}

function aaStructureEl(one) {
  const svg = AA_STRUCTURES[one];
  if (!svg) return el('span', { class: 'muted' }, '—');
  const d = el('div', { class: 'aa-struct' }); d.innerHTML = svg; return d;
}

// Reusable sortable table. columns: [{key,label,get?,render?,sortVal?,numeric?,className?,headClass?,sortable?}].
function sortableTable(columns, rows, initial = {}) {
  let sortKey = initial.key || null;
  let sortDir = initial.dir || 1;
  const thead = el('thead', {}); const headRow = el('tr', {}); const tbody = el('tbody', {}); const ths = [];

  function renderBody() {
    tbody.innerHTML = '';
    const data = rows.slice();
    if (sortKey) {
      const col = columns.find((c) => c.key === sortKey);
      data.sort((a, b) => {
        const va = col.sortVal ? col.sortVal(a) : col.get ? col.get(a) : '';
        const vb = col.sortVal ? col.sortVal(b) : col.get ? col.get(b) : '';
        const na = va == null || va === '' || va === '—'; const nb = vb == null || vb === '' || vb === '—';
        if (na && nb) return 0; if (na) return 1; if (nb) return -1;
        if (col.numeric) return (Number(va) - Number(vb)) * sortDir;
        return String(va).localeCompare(String(vb)) * sortDir;
      });
    }
    for (const row of data) {
      const tr = el('tr', {});
      for (const c of columns) {
        const td = el('td', { class: c.className || '' });
        const content = c.render ? c.render(row) : c.get ? c.get(row) : '';
        if (content != null && content.nodeType) td.append(content);
        else td.textContent = content == null ? '' : String(content);
        tr.append(td);
      }
      tbody.append(tr);
    }
  }
  function updateHeaders() { columns.forEach((c, i) => { ths[i]._arrow.textContent = sortKey === c.key ? (sortDir === 1 ? ' ▲' : ' ▼') : ''; }); }
  function setSort(key) {
    const col = columns.find((c) => c.key === key);
    if (!col || col.sortable === false) return;
    if (sortKey === key) sortDir = -sortDir; else { sortKey = key; sortDir = 1; }
    updateHeaders(); renderBody();
  }
  for (const c of columns) {
    const arrow = el('span', { class: 'sort-arrow' }, '');
    const th = el('th', { class: `${c.headClass || ''}${c.sortable === false ? '' : ' sortable'}` }, c.label, arrow);
    if (c.sortable !== false) th.addEventListener('click', () => setSort(c.key));
    th._arrow = arrow; ths.push(th); headRow.append(th);
  }
  thead.append(headRow);
  const table = el('table', { class: 'ref' }, thead, tbody);
  updateHeaders(); renderBody();
  return table;
}

function buildAmino() {
  const root = $('#view-amino');
  root.innerHTML = '';

  const legend = el('div', { class: 'legend' },
    Object.entries(POLARITY_COLOUR).map(([k, v]) => el('div', { class: 'item' },
      el('div', { class: 'sw', style: `background:${v}` }), k)));

  const columns = [
    { key: 'one', label: '1', headClass: 'mono', render: (a) => el('span', { class: 'aa-swatch', style: `background:${a.colour}` }, a.one), sortVal: (a) => a.one },
    { key: 'three', label: '3', className: 'mono', get: (a) => a.three },
    { key: 'name', label: 'Name', get: (a) => a.name },
    { key: 'polarity', label: 'Polarity', render: (a) => polTag(a.polarity), sortVal: (a) => a.polarity },
    { key: 'charge', label: 'Charge (pH 7.4)', className: 'mono', numeric: true, render: (a) => fmtCharge(a.charge), sortVal: (a) => a.charge },
    { key: 'aromatic', label: 'Aromatic', get: (a) => (a.aromatic ? 'yes' : 'no'), sortVal: (a) => (a.aromatic ? 1 : 0) },
    { key: 'hydropathy', label: 'Hydropathy (K-D)', className: 'mono', numeric: true, get: (a) => a.hydropathy },
    { key: 'mass', label: 'Mass (Da)', className: 'mono', numeric: true, get: (a) => a.mass.toFixed(2), sortVal: (a) => a.mass },
    { key: 'pka', label: 'Side-chain pKa', className: 'mono', numeric: true, get: (a) => (a.pka == null ? '—' : a.pka), sortVal: (a) => a.pka },
    { key: 'codons', label: 'Codons', className: 'mono', render: (a) => el('span', { style: 'font-size:12px' }, a.codons.join(' ')), sortVal: (a) => a.codons.length },
    { key: 'structure', label: 'Structure', sortable: false, render: (a) => aaStructureEl(a.one) },
  ];

  root.append(el('div', { class: 'card' },
    el('h2', {}, 'Amino-acid reference'),
    el('div', { class: 'sub' }, `${SCHEME} ${AROMATICITY_NOTE} Click a column heading to sort; the last column shows each residue’s chemical structure.`),
    legend,
    el('div', { style: 'overflow-x:auto' }, sortableTable(columns, AMINO_ACIDS.slice(), { key: 'one', dir: 1 }))));
}

// =====================================================================
//  CODON TABLE VIEW
// =====================================================================
function buildCodon() {
  const root = $('#view-codon');
  root.innerHTML = '';

  const legend = el('div', { class: 'legend' },
    [...Object.entries(POLARITY_COLOUR), ['stop', '#8a8f99']].map(([k, v]) => el('div', { class: 'item' },
      el('div', { class: 'sw', style: `background:${v}` }), k)));

  function buildGrid() {
    const grid = el('div', { class: 'gc' });
    for (const b1 of BASES4) {
      const block = el('div', { class: 'gc-block' });
      for (const b2 of BASES4) for (const b3 of BASES4) {
        const codon = b1 + b2 + b3; const one = STANDARD_CODE[codon]; const rec = aaByOne(one);
        const colour = one === '*' ? '#8a8f99' : (rec ? rec.colour : '#8a8f99');
        block.append(el('div', { class: 'gc-cell' },
          el('span', { class: 'codon' }, codon),
          el('span', { class: 'aa' },
            el('span', { class: 'sw', style: `background:${colour}` }, one),
            el('span', { class: 'nm' }, one === '*' ? 'Stop' : rec.three))));
      }
      grid.append(block);
    }
    return grid;
  }

  const rows = Object.keys(STANDARD_CODE).map((codon) => {
    const one = STANDARD_CODE[codon]; const rec = aaByOne(one);
    return { codon, one, three: one === '*' ? 'Ter' : rec.three, name: one === '*' ? 'Stop' : rec.name, polarity: one === '*' ? 'stop' : rec.polarity, colour: one === '*' ? '#8a8f99' : rec.colour };
  });
  const tableColumns = [
    { key: 'codon', label: 'Codon', className: 'mono', get: (r) => r.codon },
    { key: 'one', label: 'AA', render: (r) => el('span', { class: 'aa-swatch', style: `background:${r.colour}` }, r.one), sortVal: (r) => r.one },
    { key: 'three', label: '3-letter', className: 'mono', get: (r) => r.three },
    { key: 'name', label: 'Amino acid', get: (r) => r.name },
    { key: 'polarity', label: 'Polarity', render: (r) => el('span', { class: 'tag', style: `background:${r.colour}` }, r.polarity), sortVal: (r) => r.polarity },
  ];

  const holder = el('div', {});
  const gridBtn = el('button', { class: 'btn small' }, 'Classic grid');
  const tableBtn = el('button', { class: 'btn small secondary' }, 'Sortable table');
  const showGrid = () => { gridBtn.className = 'btn small'; tableBtn.className = 'btn small secondary'; holder.innerHTML = ''; holder.append(buildGrid()); };
  const showTable = () => { tableBtn.className = 'btn small'; gridBtn.className = 'btn small secondary'; holder.innerHTML = ''; holder.append(el('div', { style: 'overflow-x:auto' }, sortableTable(tableColumns, rows.slice(), { key: 'codon', dir: 1 }))); };
  gridBtn.addEventListener('click', showGrid);
  tableBtn.addEventListener('click', showTable);

  root.append(el('div', { class: 'card' },
    el('h2', {}, 'The standard genetic code'),
    el('div', { class: 'sub' }, 'NCBI translation table 1. DNA sense strand (T, not U). Coloured by the amino acid’s polarity. Switch to the sortable table and click a heading to sort by codon or amino acid.'),
    el('div', { class: 'row', style: 'gap:8px; margin-bottom:12px' }, gridBtn, tableBtn),
    legend,
    holder));
  showGrid();
}

// A two-strand codon figure for the Translate tool: the entered sequence 5′→3′ on top, its
// complementary template strand antiparallel 3′→5′ below, and the amino acid under each codon so the
// reading frame (which base begins each residue) is visible. Frame runs from the first base.
const NT_COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A' };
function buildTranslateFigure(seq) {
  const MAX = 300; // cap the drawn window so a full CDS paste stays light; the protein text shows all
  const shown = seq.slice(0, MAX);
  const fig = el('div', { class: 'codon-fig' });

  const gutter = el('div', { class: 'cf-gutter' },
    el('div', {}, '5′→3′'), el('div', {}, '3′→5′'), el('div', {}, 'aa'));
  const end = (top, bottom) => el('div', { class: 'cf-end' }, el('div', {}, top), el('div', {}, bottom));
  fig.append(gutter, end('5′', '3′'));

  for (let c = 0; c < shown.length; c += 3) {
    const codon = shown.slice(c, c + 3);
    const top = el('div', { class: 'cb-row' });
    const bot = el('div', { class: 'cb-row' });
    for (const b of codon) {
      top.append(el('div', { class: `tbase ${b}` }, b));
      const comp = NT_COMPLEMENT[b] || 'N';
      bot.append(el('div', { class: `tbase ${comp}` }, comp));
    }
    let aaCell;
    if (codon.length === 3) {
      const one = translateCodon(codon);
      const aa = one != null ? aaByOne(one) : null;
      const colour = one === '*' ? '#e0574a' : (aa ? (aa.colour || POLARITY_COLOUR[aa.polarity] || '#8a8f99') : '#8a8f99');
      aaCell = el('div', { class: 'cb-aa', title: one === '*' ? 'Ter (stop)' : (aa ? `${aa.three} — ${aa.name}` : 'unknown codon'), style: `color:${colour}` }, one == null ? '?' : one);
    } else {
      aaCell = el('div', { class: 'cb-aa muted', title: 'incomplete codon' }, '·');
    }
    fig.append(el('div', { class: 'codon-block' }, top, bot, aaCell));
  }
  fig.append(end('3′', '5′'));

  const note = el('div', { class: 'note', style: 'margin-top:6px' },
    'Top strand 5′→3′ is the sequence you entered; the complementary template strand is antiparallel below, 3′→5′. '
    + 'Each amino acid sits under its codon — the reading frame runs from the first base.'
    + (seq.length > MAX ? ` Figure shows the first ${MAX} of ${seq.length} bases; the protein above covers the whole sequence.` : ''));
  return el('div', {}, el('div', { class: 'track-wrap' }, fig), note);
}

// =====================================================================
//  SEQUENCE TOOLS VIEW
// =====================================================================
function buildTools() {
  const root = $('#view-tools');
  root.innerHTML = '';

  // Translate a sequence
  const seqInput = el('textarea', { placeholder: 'DNA coding sequence' });
  const seqOut = el('div', { class: 'kv', style: 'margin-top:12px' });
  const seqFig = el('div', { style: 'margin-top:12px' });
  const translateBtn = el('button', { class: 'btn', onclick: () => {
    const raw = seqInput.value.replace(/[^A-Za-z]/g, '').toUpperCase();
    seqOut.innerHTML = ''; seqFig.innerHTML = '';
    if (!raw) { seqOut.append(el('div', { class: 'note' }, 'Enter a sequence.')); return; }
    const t = translate(raw, { stopAtStop: true });
    const three3 = t.protein.split('').map((o) => three(o)).join('-');
    seqOut.append(
      el('div', { class: 'k' }, 'Length'), el('div', { class: 'v' }, `${raw.length} nt, ${t.protein.length} aa${t.stoppedAtStop ? ' (stop reached)' : ''}`),
      el('div', { class: 'k' }, 'Protein (1-letter)'), el('div', { class: 'v' }, t.protein || '—'),
      el('div', { class: 'k' }, 'Protein (3-letter)'), el('div', { class: 'v' }, three3 || '—'),
      el('div', { class: 'k' }, 'Reverse complement'), el('div', { class: 'v' }, reverseComplement(raw)),
      el('div', { class: 'k' }, 'mRNA (T→U)'), el('div', { class: 'v' }, toMrna(raw)));
    seqFig.append(buildTranslateFigure(raw));
  } }, 'Translate');

  const translateCard = el('div', { class: 'card' },
    el('h2', {}, 'Translate a sequence'),
    el('div', { class: 'sub' }, 'Standard genetic code, reading from the first base. Fully offline.'),
    seqInput, el('div', { style: 'margin-top:10px' }, translateBtn), seqOut, seqFig);

  // Manual codon assessor
  const refC = el('input', { type: 'text', placeholder: 'Reference codon', maxlength: '3', style: 'max-width:160px' });
  const altC = el('input', { type: 'text', placeholder: 'Altered codon', maxlength: '3', style: 'max-width:160px' });
  const codonOut = el('div', { style: 'margin-top:14px' });
  const codonBtn = el('button', { class: 'btn', onclick: () => {
    codonOut.innerHTML = '';
    const rc = refC.value.replace(/[^A-Za-z]/g, '').toUpperCase();
    const ac = altC.value.replace(/[^A-Za-z]/g, '').toUpperCase();
    if (rc.length !== 3 || ac.length !== 3) { codonOut.append(el('div', { class: 'err' }, 'Enter two three-base codons.')); return; }
    const ra = translateCodon(rc), aa = translateCodon(ac);
    if (ra == null || aa == null) { codonOut.append(el('div', { class: 'err' }, 'Codons must use only A, C, G, T.')); return; }
    let cons = 'missense';
    if (ra === aa) cons = 'synonymous';
    else if (aa === '*') cons = 'nonsense';
    else if (ra === '*') cons = 'stop-loss';
    codonOut.append(el('div', { class: 'aa-change' }, aaTile(ra, rc), el('div', { class: 'arrow' }, '→'), aaTile(aa, ac)),
      el('div', { class: `consequence c-${cons}`, style: 'margin-top:8px' }, el('span', { class: 'dot' }), consequenceLabel(cons)));
  } }, 'Compare codons');

  const codonCard = el('div', { class: 'card' },
    el('h2', {}, 'Compare two codons'),
    el('div', { class: 'sub' }, 'See the amino-acid effect of a codon change without a transcript — for example TAC → TCC. Fully offline.'),
    el('div', { class: 'row' }, refC, el('span', { class: 'arrow' }, '→'), altC, codonBtn),
    codonOut);

  root.append(el('div', { class: 'grid2' }, translateCard, codonCard));
}

// =====================================================================
//  ABOUT VIEW
// =====================================================================
async function buildAbout() {
  const root = $('#view-about');
  root.innerHTML = '';
  let info = {};
  try { info = await window.api.appInfo(); } catch { /* */ }
  if (!transcriptIndex.length) { try { transcriptIndex = await window.api.getIndex(); } catch { /* */ } }

  const src = (label, name, url) => el('li', {}, label + ' — ',
    el('a', { onclick: () => window.api.openExternal(url), style: 'cursor:pointer' }, name));
  const sources = el('ul', { class: 'sources' },
    src('Reference transcripts and proteins', 'NCBI RefSeq · MANE Select', 'https://www.ncbi.nlm.nih.gov/refseq/MANE/'),
    src('Gene → protein, function and canonical accession', 'UniProtKB', 'https://www.uniprot.org/'),
    src('Protein domains (Pfam, SMART, PROSITE, CDD, PANTHER), mapped to MANE', 'InterPro', 'https://www.ebi.ac.uk/interpro/'),
    src('Gene coordinates and strand', 'Ensembl', 'https://www.ensembl.org/'),
    src('Predicted 3D structures', 'AlphaFold DB', 'https://alphafold.ebi.ac.uk/'),
    src('Experimental 3D structures', 'RCSB PDB', 'https://www.rcsb.org/'),
    src('Clinical significance (germline and somatic)', 'ClinVar', 'https://www.ncbi.nlm.nih.gov/clinvar/'),
    src('Population allele frequency', 'gnomAD v4', 'https://gnomad.broadinstitute.org/'),
    src('Splice-impact prediction', 'SpliceAI (Broad lookup)', 'https://spliceailookup.broadinstitute.org/'),
    src('Experimental splice evidence', 'SpliceVarDB', 'https://splicevardb.org/'),
    src('One-shot clinical/population annotation engine', 'GeneBe', 'https://genebe.net/'),
    src('Curated haematology variants (COSMIC, GENIE, ClinVar, OncoKB, TCGA, CancerHotspots)', 'OncoKB', 'https://www.oncokb.org/'),
    src('Genetic code (standard translation table 1)', 'NCBI', 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi'),
    el('li', {}, 'Amino-acid properties — Kyte–Doolittle hydropathy (1982) and residue masses from ',
      el('a', { onclick: () => window.api.openExternal('https://web.expasy.org/protparam/'), style: 'cursor:pointer' }, 'ExPASy'), '.'));
  const sourcesNote = el('div', { class: 'sub', style: 'margin-top:8px' },
    IS_WEB
      ? 'ClinVar, SpliceAI, SpliceVarDB and GeneBe are network annotations available in the desktop app; this web build runs the offline core plus the gnomAD lookup.'
      : 'gnomAD sends only a public coordinate; GeneBe (opt-out) sends the variant. All other network sources take a public accession or coordinate, never patient data, and cache locally.');

  // Bundled transcripts — accurate scope for each build (web ships the haem panel; desktop the whole store).
  const geneCount = new Set(transcriptIndex.map((t) => t.gene)).size;
  let scope;
  if (IS_WEB) {
    scope = `Every MANE Select protein-coding gene is bundled — ${transcriptIndex.length.toLocaleString('en-GB')} transcripts across ${geneCount.toLocaleString('en-GB')} genes — each fetched on demand and resolved in your browser, nothing sent anywhere. A curated set of deliberately-versioned transcripts (e.g. non-MANE accessions such as MYD88 NM_002468.4) takes priority on lookup.`;
  } else if (info.referenceCount) {
    scope = `Every MANE Select protein-coding gene resolves offline — ${info.referenceCount.toLocaleString('en-GB')} transcripts bundled. A curated set of ${transcriptIndex.length} deliberately-versioned transcripts (e.g. non-MANE accessions such as MYD88 NM_002468.4) takes priority on lookup, and any other versioned accession can be fetched from NCBI on demand.`;
  } else {
    scope = `${transcriptIndex.length} transcripts are bundled and resolve offline.`;
  }

  root.append(
    el('div', { class: 'card' },
      el('h2', {}, 'Sources'), sources, sourcesNote),
    el('div', { class: 'card' },
      el('h2', {}, 'Bundled transcripts'),
      el('div', { class: 'note' }, scope)));

  // Build/diagnostics are meaningful only in the desktop app; on the web they are all blank.
  if (!IS_WEB) {
    root.append(el('div', { class: 'card' },
      el('h2', {}, 'Build'),
      el('div', { class: 'kv' },
        el('div', { class: 'k' }, 'Electron'), el('div', { class: 'v' }, info.versions?.electron || '—'),
        el('div', { class: 'k' }, 'Chromium'), el('div', { class: 'v' }, info.versions?.chrome || '—'),
        el('div', { class: 'k' }, 'Node'), el('div', { class: 'v' }, info.versions?.node || '—'),
        el('div', { class: 'k' }, 'Bundled data'), el('div', { class: 'v' }, info.bundledData || '—'),
        el('div', { class: 'k' }, 'Fetch cache'), el('div', { class: 'v' }, info.cacheDir || '—'))));
  }
}

// =====================================================================
//  LOLLIPOP BOARD — several variants drawn together for a workbook
// =====================================================================
function buildBoard() {
  const root = $('#view-board');
  root.innerHTML = '';
  const input = el('textarea', {
    id: 'board-input', spellcheck: 'false', style: 'min-height:132px',
    placeholder: 'NM_004333.6:c.1799T>A\nNM_015338.5:c.2945dup\nNM_000546.6:c.818G>A',
  });
  const status = el('span', { id: 'board-status', class: 'note', style: 'margin-left:10px' });

  root.append(
    el('div', { class: 'card' },
      el('h2', {}, 'Lollipop board'),
      el('div', { class: 'sub' }, 'Draw several variants as domain lollipop diagrams in one view, ready to screenshot for a variant assessment workbook. Variants on the same transcript share a diagram (labels stagger automatically); each gene gets its own.'),
      input,
      el('div', { class: 'row', style: 'margin-top:10px' },
        el('button', { class: 'btn', onclick: () => drawBoard() }, 'Draw diagrams'),
        el('button', { class: 'btn small secondary', onclick: () => {
          input.value = ['NM_004333.6:c.1799T>A', 'NM_000546.6:c.524G>A', 'NM_000546.6:c.743G>A', 'NM_000546.6:c.818G>A', 'NM_015338.5:c.2945dup', 'NM_002520.7:c.860_863dup'].join('\n');
        } }, 'Example set'),
        el('button', { class: 'btn small secondary', onclick: () => { input.value = ''; $('#board-output').innerHTML = ''; status.textContent = ''; } }, 'Clear'),
        status)),
    el('div', { id: 'board-output' }));
}

async function drawBoard() {
  const out = $('#board-output');
  const status = $('#board-status');
  out.innerHTML = '';
  const lines = $('#board-input').value.split('\n').map((s) => s.trim()).filter(Boolean);
  if (!lines.length) { out.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'Enter at least one variant, one per line.'))); return; }
  status.textContent = 'Drawing…';

  const groups = new Map();
  const problems = [];
  const missing = new Set();
  for (const line of lines) {
    const p = parseHgvs(line);
    if (!p.ok || !p.transcript || !p.supported) { problems.push(`${line} — ${p.message || 'could not be resolved'}`); continue; }
    let lookup = { doc: null };
    try { lookup = await window.api.getTranscript(p.transcript); } catch { /* offline */ }
    if (!lookup.doc) {
      missing.add(p.transcript);
      const others = (lookup.otherVersions || []).filter((v) => v !== p.transcript);
      problems.push(`${line} — transcript ${p.transcript} is not bundled${others.length ? `; you do have ${others.join(', ')} offline (versions can differ in numbering, so it is not substituted automatically)` : ''}`);
      continue;
    }
    const doc = lookup.doc;

    let residue; let label; let consequence;
    if (p.kind === 'substitution') {
      const r = resolveSubstitution(doc.cds, { position: p.position, ref: p.ref, alt: p.alt });
      if (!r.ok) { problems.push(`${line} — ${r.message}`); continue; }
      if (!r.refMatches) problems.push(`${line} — drawn, but the stated reference base does not match ${doc.accession}`);
      ({ codonNumber: residue, proteinShort: label, consequence } = r);
    } else {
      const r = resolveIndel(doc.cds, p, { utr3: doc.utr3 });
      if (!r.ok) { problems.push(`${line} — ${r.message}`); continue; }
      residue = r.firstAffected; label = r.proteinShort; consequence = r.consequence;
    }
    if (!groups.has(doc.accession)) groups.set(doc.accession, { doc, markers: [] });
    groups.get(doc.accession).markers.push({ residue, label, consequence });
  }

  // Compact panels: header + diagram + domain descriptors only. The per-residue
  // caption and the domain-source line are dropped here (one source note is added
  // once at the end) so a whole board screenshots cleanly into a workbook.
  for (const { doc, markers } of groups.values()) {
    out.append(el('div', { class: 'card board-card' },
      el('div', { class: 'row', style: 'justify-content:space-between; align-items:baseline; gap:12px' },
        el('h2', {}, `${doc.gene} — ${markers.map((m) => m.label).join(', ')}`),
        el('div', { class: 'sub' }, `${doc.accession} · ${doc.proteinLength} aa${doc.strand ? ` · ${doc.strand} strand` : ''}`)),
      buildLollipop(doc, markers, { heading: false, caption: false, source: false })));
  }
  if (groups.size) {
    // One footer for a board that can hold several transcripts, so it can name neither one
    // accession's numbering nor MANE: each card already carries its own accession above its diagram.
    const sources = [...new Set([...groups.values()].map(({ doc }) => doc.domainSource || 'InterPro/UniProt'))];
    out.append(el('div', { class: 'note muted', style: 'margin:2px 4px 14px; font-size:11.5px' },
      `Domains from ${sources.join('; ')} · positions are in each diagram's own transcript numbering, named above it.`));
  }

  if (missing.size) {
    const btn = el('button', { class: 'btn small' }, `Fetch ${missing.size} missing transcript${missing.size > 1 ? 's' : ''} from NCBI`);
    btn.addEventListener('click', async () => {
      btn.disabled = true; btn.textContent = 'Fetching…';
      for (const acc of missing) { try { await window.api.fetchTranscript(acc); } catch { /* reported on redraw */ } }
      drawBoard();
    });
    out.append(el('div', { class: 'card' },
      el('h2', {}, 'Missing transcripts'),
      el('div', { class: 'note' }, `Not bundled: ${[...missing].join(', ')}. Fetching sends only the accession to NCBI — never your variant.`),
      el('div', { class: 'note muted', style: 'margin-top:6px' }, 'A fetched transcript resolves the variant, but has no bundled domain map or 3′UTR, so its diagram shows no domains and a frameshift stop may read as “?”. Re-run the data scripts to bundle it fully.'),
      el('div', { class: 'row', style: 'margin-top:10px' }, btn)));
  }

  if (problems.length) {
    out.append(el('div', { class: 'card' }, el('h2', {}, 'Not drawn'),
      el('ul', {}, problems.map((t) => el('li', { class: 'note' }, t)))));
  }
  const drawn = [...groups.values()].reduce((n, g) => n + g.markers.length, 0);
  status.textContent = `${groups.size} diagram(s), ${drawn} variant(s).`;
}

// =====================================================================
//  KNOWN VARIANTS
// =====================================================================
// One searchable list across the curated reference set and the haematology variants imported from
// the whitelist, merged so one variant is one row. It feeds the lookup card on the Assess page: a
// search box that shows the best few matches, and a drop-down of the fifty most reported.
let knownRows = null;
// Why the haematology list is missing, when it is. "No match." is a statement about the list, so it
// may not stand in for a list that never arrived — the same rule the Fusions tab follows for its
// reference gene table.
let knownListError = null;

function knownRowsAll() {
  if (knownRows) return knownRows;
  knownRows = KNOWN_VARIANTS.map((v) => ({
    gene: v.gene,
    hgvsc: v.hgvsC,
    transcript: (v.hgvsC || '').split(':')[0],
    protein: v.expected,
    hg38: v.chrPos ? `${v.chrPos} ${v.genomicRefAlt || ''}`.trim() : '',
    samples: null,
    evidence: v.note || '',
    context: v.category || '',
  }));
  return knownRows;
}

const stripProteinPrefix = (x) => String(x || '').replace(/^p[.]/, '');

/** Fold one row's annotation into another. Nothing is dropped: sample counts take the larger, the
 *  evidence notes are unioned, and the first non-empty context is kept. */
function foldKnownRow(prim, sec) {
  if (sec.samples != null && (prim.samples == null || sec.samples > prim.samples)) prim.samples = sec.samples;
  prim.evidence = [...new Set([prim.evidence, sec.evidence].filter(Boolean).join(' · ').split(' · ').filter(Boolean))].join(' · ');
  if (!prim.context) prim.context = sec.context;
  return prim;
}

/**
 * The merged known-variant list (whitelist + curated). Cached.
 *
 * One variant must be one row, and every row must keep everything both sources said about it.
 * Three stages, in order, each with its own rule for what counts as "the same variant":
 *   1. the same HGVS coding string — the same change on the same transcript;
 *   2. the same genomic change (chromosome, position and alleles) in the same gene;
 *   3. the same coding change and protein change in the same gene, for a row that carries no
 *      genomic coordinate and so cannot reach stage 2.
 * Stage 3 claims only what the two rows' own strings say — that they write the same coding change —
 * and names the other transcript on the row rather than asserting a coordinate the row does not have.
 */
async function knownMergedRows() {
  if (knownRows && knownRows._whitelistMerged) return knownRows;
  const byKey = new Map();
  const keyOf = (hgvsc) => String(hgvsc || '').trim().toUpperCase();

  let wl = null;
  knownListError = null;
  try { wl = await window.api.getWhitelist(); } catch (e) { knownListError = e.message || String(e); }
  // The desktop provider now reports a failed read explicitly; the browser one, and any older build,
  // still turns it into an empty list. Treat both as a failure: the shipped file holds 2,141 variants
  // and is never legitimately empty.
  if (!knownListError && wl && wl.failed) knownListError = wl.reason || 'the list could not be read';
  if (!knownListError && !(wl && Array.isArray(wl.variants) && wl.variants.length)) {
    knownListError = 'the list came back empty';
  }
  for (const v of (wl && wl.variants) || []) {
    byKey.set(keyOf(v.hgvsc), {
      gene: v.gene,
      hgvsc: v.hgvsc,
      transcript: v.transcript,
      protein: v.proteinChange,
      hg38: v.hg38,
      samples: v.samples,
      evidence: [v.oncokb, v.clinvar].filter(Boolean).join(' · '),
      context: v.cancerTypes || '',
    });
  }
  // Stage 1. A curated row whose coding string already exists keeps the whitelist row's coordinate
  // and sample count, but its note and category are folded in rather than thrown away — 55 curated
  // rows lost both here, so JAK2 V617F's "MPN driver (PV, ET, PMF)" never reached the screen.
  for (const r of knownRowsAll()) {
    const key = keyOf(r.hgvsc);
    if (byKey.has(key)) foldKnownRow(byKey.get(key), r);
    else byKey.set(key, r);
  }
  // Stage 2. The same genomic change on two versions of one transcript is ONE variant: the newer
  // version's annotation wins (MYD88 L252P on NM_002468.5 over L265P on .4, same chr3:38,141,150
  // T>C) and the older protein name is kept as "previously called".
  // "Newer" only means anything within one accession series. Comparing the trailing version numbers
  // of two different series (NM_033360.4 against NM_004985.5) compares two unrelated counters, so
  // when the series differ neither row supersedes the other: the first seen stays primary and the
  // other name is recorded as an alternative rather than as a former name.
  const txVer = (t) => { const m = /[.](\d+)$/.exec(String(t || '')); return m ? Number(m[1]) : 0; };
  const txBase = (t) => String(t || '').replace(/[.]\d+$/, '');
  const gkey = (r) => {
    const g = String(r.hg38 || '').replace(/^chr/i, '').replace(/[\s,]/g, '').toUpperCase();
    return g ? `${r.gene}|${g}` : null;
  };
  const byG = new Map(); const merged = [];
  for (const r of byKey.values()) {
    const k = gkey(r);
    if (!k) { merged.push(r); continue; }
    const prev = byG.get(k);
    if (!prev) { byG.set(k, r); merged.push(r); continue; }
    const sameSeries = txBase(r.transcript) === txBase(prev.transcript);
    let prim = prev; let sec = r;
    if (sameSeries && txVer(r.transcript) > txVer(prev.transcript)) { prim = r; sec = prev; }
    foldKnownRow(prim, sec);
    if (stripProteinPrefix(sec.protein) && stripProteinPrefix(sec.protein) !== stripProteinPrefix(prim.protein)) {
      prim.alias = stripProteinPrefix(sec.protein);
      prim.aliasTx = sec.transcript;
      prim.aliasIsFormerName = sameSeries;
    }
    if (prim !== prev) { merged[merged.indexOf(prev)] = prim; byG.set(k, prim); }
  }
  // Stage 3. A row with no coordinate never reaches stage 2, so the eight curated KRAS hotspots on
  // NM_033360.4 sat beside the identical whitelist rows on NM_004985.5 — the same hotspot listed
  // twice, once with sample counts and once without. They are folded into the row that does carry a
  // coordinate when gene, coding change and protein change all match; the transcript that had no
  // coordinate is kept on the row so nothing is lost. The protein change is part of the test because
  // one coding position can carry more than one change: KRAS p.Q61H is c.183A>T in one whitelist row
  // and c.183A>C in another, two different alleles that must not be folded together.
  const codingChange = (h) => String(h || '').split(':').slice(1).join(':').trim().toUpperCase();
  const ckey = (r) => {
    const c = codingChange(r.hgvsc); const p = stripProteinPrefix(r.protein).toUpperCase();
    return c && p && r.gene ? `${r.gene}|${c}|${p}` : null;
  };
  const byC = new Map();
  for (const r of merged) {
    if (!gkey(r)) continue;
    const k = ckey(r);
    if (!k) continue;
    // A key claimed by two coordinate-carrying rows is ambiguous; fold into neither.
    byC.set(k, byC.has(k) ? null : r);
  }
  const kept = [];
  for (const r of merged) {
    const host = gkey(r) ? null : byC.get(ckey(r));
    if (!host) { kept.push(r); continue; }
    foldKnownRow(host, r);
    if (r.transcript && r.transcript !== host.transcript) host.alsoTranscript = r.transcript;
  }
  if (knownListError) {
    // Nothing is cached: the curated rows are returned so the search still answers for them, but a
    // retry has to be able to rebuild the list once the file arrives.
    knownRows = null;
    return kept;
  }
  knownRows = kept;
  knownRows._whitelistMerged = true;
  return knownRows;
}

// The other protein name a merged row carries, in words. Only a newer version of the SAME
// transcript can be said to have replaced an earlier name; two different accession series are just
// two names in use, and saying "previously" of those would invent a history.
function knownAliasText(r) {
  if (!r.alias) return '';
  return r.aliasIsFormerName ? `previously ${r.alias}` : `also called ${r.alias} on ${r.aliasTx}`;
}

// The merged evidence field is raw source text and cannot go on screen as it stands: the ClinVar
// strings carry underscores, the disease list carries US spellings, ICD-O fragments ("nos",
// "(c44._)") and unexplained short forms, and the curated notes are written in clinical
// abbreviations (MPN, PV, ET, PMF, WM, ABC-DLBCL).
// The two classification vocabularies are the part that can be stated properly: seven values across
// the whole shipped list, so each is spelled out here. Anything this map does not recognise is not
// shown — the search box still matches on it, so nothing is lost, but nothing raw is printed.
const CLASSIFICATION_WORDS = new Map([
  ['oncogenic', 'Oncogenic (OncoKB)'],
  ['likely oncogenic', 'Likely oncogenic (OncoKB)'],
  ['pathogenic', 'Pathogenic (ClinVar)'],
  ['likely_pathogenic', 'Likely pathogenic (ClinVar)'],
  ['pathogenic/likely_pathogenic', 'Pathogenic or likely pathogenic (ClinVar)'],
  ['pathogenic|other', 'Pathogenic, with other classifications also submitted (ClinVar)'],
  ['pathogenic/likely_pathogenic/pathogenic,_low_penetrance', 'Pathogenic or likely pathogenic, some submitted as low penetrance (ClinVar)'],
]);

/** The classification words a merged row carries, in readable form. Empty when it carries none. */
function knownEvidenceText(r) {
  const out = [];
  for (const token of String(r.evidence || '').split(' · ')) {
    const word = CLASSIFICATION_WORDS.get(token.trim().toLowerCase());
    if (word && !out.includes(word)) out.push(word);
  }
  return out.join(' · ');
}

/** Compact known-variant lookup on the Assess page: a search cell plus the canonical drop-down,
 *  ordered by how often each variant is reported (sample count), not alphabetically. */
function knownLookupCard() {
  const stripP = stripProteinPrefix;
  const kSearch = el('input', {
    type: 'search', placeholder: 'Search known variants — gene, protein change, transcript, hg38…',
    style: 'flex:1; min-width:240px; padding:8px 10px; border-radius:8px',
  });
  const placeholder = () => el('option', { value: '' }, 'Canonical variant…');
  const kSel = el('select', { id: 'canon-variant', class: 'fus-picker', style: 'max-width:360px' }, placeholder());
  const results = el('div', { class: 'known-results' });
  // A failed load is stated as soon as it happens, not only once something is typed. Before this, a
  // user who never typed saw a drop-down holding its placeholder alone and a card still offering to
  // search a list that had not arrived.
  const status = el('div', { class: 'known-status' });
  let rows = [];
  let loaded = false;
  // Named rather than inline so the "Try again" button below can run it after a failed load.
  const load = async () => {
    rows = await knownMergedRows();
    loaded = true;
    kSel.innerHTML = '';
    kSel.append(placeholder());
    const top = rows.filter((r) => r.samples != null)
      .sort((a, b) => b.samples - a.samples)
      .slice(0, 50);
    for (const r of top) {
      kSel.append(el('option', { value: r.hgvsc },
        `${r.gene} ${stripP(r.protein)}${r.alias ? ` (${knownAliasText(r)})` : ''} · ${r.samples.toLocaleString()} samples`));
    }
    // Say up front whether the list arrived, whatever is or is not typed in the search box.
    status.innerHTML = '';
    if (knownListError) {
      status.append(el('div', { class: 'note err', style: 'font-size:12px; margin-top:8px' },
        `The haematology variant list did not load, so only the ${rows.length} curated variants can be searched. Reason: ${knownListError}.`));
      status.append(el('button', { class: 'btn small', style: 'margin-top:6px', onclick: () => { loaded = false; status.innerHTML = ''; render(); load(); } }, 'Try again'));
    }
    // A search typed while the list was still loading was answered "No match." Redraw once it is in.
    render();
  };
  load();
  kSel.onchange = () => {
    if (!kSel.value) return;
    const v = kSel.value;
    kSel.value = '';
    loadVariant(v);
  };
  const render = () => {
    const q = kSearch.value.trim().toLowerCase();
    results.innerHTML = '';
    if (q.length < 2) return;
    // "No match" is a statement about the list. Until the list is in, the honest answer is that it
    // is still loading.
    if (!loaded) { results.append(el('div', { class: 'note muted', style: 'font-size:12px; margin-top:8px' }, 'Loading the known-variant list…')); return; }
    // The missing-list notice lives above, in `status`, so it is shown whether or not anything is
    // typed. Nothing to repeat here.
    const terms = q.split(/\s+/);
    const hits = rows.filter((r) => {
      const hay = `${r.gene} ${r.protein} ${r.alias || ''} ${r.hgvsc} ${r.alsoTranscript || ''} ${r.hg38} ${r.context || ''} ${r.evidence || ''}`.toLowerCase();
      return terms.every((t) => hay.includes(t));
    }).sort((a, b) => (b.samples || 0) - (a.samples || 0)).slice(0, 12);
    if (!hits.length) {
      results.append(el('div', { class: 'note muted', style: 'font-size:12px; margin-top:8px' },
        knownListError ? 'No match among the curated variants.' : 'No match.'));
      return;
    }
    for (const r of hits) {
      // Only the classification words go on screen. The curated note and the disease list are raw
      // source text — abbreviations, US spellings, underscores — so they stay in the search haystack
      // above rather than being printed as they are written.
      const note = knownEvidenceText(r);
      results.append(el('button', { class: 'known-hit', onclick: () => loadVariant(r.hgvsc) },
        el('strong', {}, r.gene), ' ',
        el('span', { class: 'mono' }, stripP(r.protein)),
        r.alias ? el('span', { class: 'muted', style: 'font-size:11px' }, ` (${knownAliasText(r)})`) : null,
        r.alsoTranscript ? el('span', { class: 'muted', style: 'font-size:11px' }, ` · same coding change on ${r.alsoTranscript}`) : null,
        el('span', { class: 'muted tail' }, `${r.samples != null ? `${r.samples.toLocaleString()} samples · ` : ''}${r.hgvsc}`),
        note ? el('span', { class: 'known-note' }, note) : null));
    }
  };
  kSearch.addEventListener('input', render);
  return el('div', { class: 'card' },
    el('h2', {}, 'Known variants'),
    el('div', { class: 'sub' }, 'Search the merged haematology list, or pick one of the fifty most frequently reported. Either loads straight into the assessor.'),
    el('div', { class: 'row', style: 'gap:10px; flex-wrap:wrap' }, kSearch, kSel),
    status,
    results);
}

function loadVariant(hgvs) {
  switchTab('assess');
  $('#variant-input').value = hgvs;
  runAssess();
}

// =====================================================================
//  GENEBE ANNOTATION VIEW  (sends the variant to genebe.net)
// =====================================================================
let settings = {};
let latestGeneBe = { status: 'none' };

function autoAnnotateOn() { return !(settings.genebe && settings.genebe.autoAnnotate === false); }

/**
 * The colour class for a classification word, across all three ClinVar vocabularies.
 *
 * The word "pathogenicity" contains "pathogenic", so ClinVar's very common aggregate description
 * "Conflicting classifications of pathogenicity" was matched by the pathogenic test and given the
 * solid red tag — a variant laboratories disagree about, coloured exactly like a settled Pathogenic
 * call. Conflicting is therefore tested first, and has a colour of its own.
 *
 * The somatic vocabularies used to match nothing and came out in the grey "unknown" style, so an
 * Oncogenic call looked the same as an unclassified one. Oncogenicity is the somatic counterpart of
 * pathogenicity and takes the same colours; the clinical-impact tiers are levels of evidence for
 * treatment rather than a pathogenicity call, so they get their own.
 */
function verdictClass(v) {
  const s = String(v || '').toLowerCase();
  if (/conflicting/.test(s)) return 'v-conflict';
  if (/(pathogenic|oncogenic)/.test(s) && !/benign/.test(s)) return 'v-path';
  if (/benign/.test(s)) return 'v-benign';                 // includes "Tier IV - Benign/Likely benign"
  if (/^tier\s/.test(s)) return 'v-tier';
  if (/uncertain|vus|risk|drug|other/.test(s)) return 'v-vus';
  return 'v-unknown';
}

function updateGenebePill() {
  const pill = $('#genebe-pill'); if (!pill) return;
  const s = latestGeneBe.status;
  const map = { none: 'Assessment: —', off: 'Assessment: off', loading: 'Assessment: …', error: 'Assessment: error' };
  if (s === 'ok') pill.textContent = `Assessment: ${latestGeneBe.data?.annotation?.acmg_classification || 'done'}`;
  else pill.textContent = map[s] || 'Assessment: —';
}

async function annotateGeneBe(hgvs) {
  const seq = currentAssessment();
  latestGeneBe = { status: 'loading', hgvs };
  updateGenebePill(); renderGeneBe();
  let res;
  try { res = await window.api.genebeAnnotate(hgvs); } catch (e) { res = { ok: false, error: e.message }; }
  // A slow answer for an earlier variant must not overwrite a newer one: the pill shows only the
  // classification, with no variant name, so a stale verdict there is unreadable as stale.
  if (superseded(seq)) return;
  latestGeneBe = res && res.ok ? { status: 'ok', hgvs, data: res } : { status: 'error', hgvs, error: (res && res.error) || 'unknown error' };
  // Reuse GeneBe's genomic coordinate for the gnomAD/ClinVar/Splice tabs so they don't re-resolve it.
  if (res && res.ok && res.genomic && lastVariant && lastVariant.hgvs === hgvs) {
    lastVariant.coord = { chr: res.genomic.chr, pos: res.genomic.pos, ref: res.genomic.ref, alt: res.genomic.alt };
  }
  updateGenebePill(); renderGeneBe();
}

function buildGeneBe() {
  const root = $('#view-genebe'); root.innerHTML = '';

  const auto = el('input', { type: 'checkbox', id: 'gb-auto' }); auto.checked = autoAnnotateOn();
  auto.addEventListener('change', async () => {
    settings = await window.api.setSettings({ genebe: { autoAnnotate: auto.checked } });
    updateGenebePill();
  });

  root.append(
    el('div', { class: 'card' },
      el('h2', {}, 'Variant assessment'),
      el('div', { class: 'sub' }, 'A quick comprehensive annotation — ACMG, ClinVar, population frequency, in-silico predictions and the RefSeq transcripts — cross-checked against this app’s own offline resolution.'),
      el('div', { class: 'note warn', style: 'margin-top:8px' }, 'This is the only tab that sends the variant off this machine — to GeneBe (genebe.net), the annotation engine. Everything else in Codon Compass stays offline.'),
      el('div', { class: 'row', style: 'align-items:center; gap:8px; margin-top:10px' }, auto, el('label', { for: 'gb-auto' }, 'Annotate each variant automatically when I assess it')),
      el('div', { class: 'disclaimer', style: 'margin-top:8px' }, 'The ACMG/pathogenicity calls are GeneBe’s automated computation — not a clinical determination. Verify against your laboratory’s validated interpretation.')),
    el('div', { id: 'genebe-result' }));
  renderGeneBe();
}

// GeneBe's 1-letter protein change (e.g. "V600E") for a given transcript, from its per-transcript
// consequences — used to cross-check against the app's own offline resolution.
function geneBeProteinChange(a, transcript) {
  const base = (t) => String(t || '').replace(/\.\d+$/, '');
  const cons = a.consequences || [];
  let c = cons.find((x) => x.feature === transcript)
    || cons.find((x) => base(x.feature) === base(transcript))
    || cons.find((x) => x.aa_ref && x.aa_alt);
  if (c && c.aa_ref && c.aa_start && c.aa_alt) return { text: `${c.aa_ref}${c.aa_start}${c.aa_alt}`, transcript: c.feature };
  return null;
}

// A compact tile for one in-silico predictor: label, score, colour-coded call, and a 0–1 bar where
// the score is a probability. Returns null when the predictor has no value (so nulls don't clutter).
function predTile(name, score, pred, opts = {}) {
  if ((score == null || score === '') && (pred == null || pred === '')) return null;
  const num = Number(score);
  const frac = opts.range ? Math.max(0, Math.min(1, (num - opts.range[0]) / (opts.range[1] - opts.range[0]))) : (isFinite(num) && num >= 0 && num <= 1 ? num : null);
  const cls = verdictClass(pred || '');
  const barColour = cls === 'v-path' ? '#e05561' : cls === 'v-benign' ? '#46b877' : '#e0a54e';
  return el('div', { style: 'flex:1 1 150px; min-width:140px; background:#152232; border:1px solid #26405c; border-radius:9px; padding:9px 11px' },
    el('div', { class: 'sub', style: 'font-size:11px; text-transform:uppercase; letter-spacing:.03em' }, name),
    el('div', { class: 'row', style: 'align-items:baseline; gap:8px; margin-top:2px' },
      el('div', { class: 'mono', style: 'font-size:17px; font-weight:700; color:#dbe9fb' }, score == null || score === '' ? '—' : String(score)),
      pred ? el('span', { class: `tag ${cls}`, style: 'font-size:10px' }, pred) : null),
    frac == null ? null : el('div', { style: 'margin-top:6px; height:5px; background:#0e1826; border-radius:3px; overflow:hidden' },
      el('div', { style: `width:${(frac * 100).toFixed(0)}%; height:100%; background:${barColour}` })));
}

function fmtAf(v) {
  if (v == null || v === '') return '—';
  const n = Number(v); if (!isFinite(n)) return String(v); if (n === 0) return '0';
  return `${n.toExponential(2)} (${(n * 100).toPrecision(3)}%)`;
}
function predRow(name, score, pred) {
  if ((score == null || score === '') && (pred == null || pred === '')) return null;
  return el('tr', {},
    el('td', {}, name),
    el('td', { class: 'mono' }, score == null || score === '' ? '—' : String(score)),
    el('td', {}, pred ? el('span', { class: `tag ${verdictClass(pred)}` }, pred) : '—'));
}

function renderGeneBe() {
  const box = $('#genebe-result'); if (!box) return;
  box.innerHTML = '';
  const s = latestGeneBe.status;
  if (s === 'none') { box.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'Assess a variant to see its assessment here.'))); return; }
  if (s === 'off') { box.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'Automatic annotation is off. Tick the box above to annotate the variants you assess.'))); return; }
  if (s === 'loading') { box.append(el('div', { class: 'card' }, el('div', { class: 'note' }, `Annotating ${latestGeneBe.hgvs}…`))); return; }
  if (s === 'error') { box.append(el('div', { class: 'card' }, el('h2', {}, latestGeneBe.hgvs || 'Assessment'), el('div', { class: 'err' }, `Could not annotate this variant: ${latestGeneBe.error}`))); return; }
  renderGeneBeAnnotation(box, latestGeneBe.data);
}

function starStr(n) { return '★'.repeat(Math.max(0, n)) + '☆'.repeat(Math.max(0, 4 - n)); }

// Where the ClinVar request goes.
//
// Desktop: through the main process, which is the only part of the app allowed to talk to NCBI and
// which caches the answer under userData. The renderer used to call NCBI itself, so nothing was
// cached and the two tabs that draw this card each fired their own search and summary for the same
// variant against NCBI's three-a-second limit for anonymous callers — making the "ClinVar
// unavailable" path more likely than it needed to be, and leaving the caching handler unreachable.
//
// Browser: there is no main process, so the page asks NCBI directly, as it always has. The web
// build's own stub answers "Desktop app only", which would take a working card away from it.
async function clinvarRequest(coord, gene) {
  if (!IS_WEB && window.api && typeof window.api.clinvar === 'function') {
    try { return await window.api.clinvar(coord, gene); } catch (e) { return { ok: false, error: e.message }; }
  }
  try { return { ok: true, ...(await clinvarByCoordinate(coord, gene)) }; } catch (e) { return { ok: false, error: e.message }; }
}
async function clinvarSubmissionRequest(variationId) {
  if (!IS_WEB && window.api && typeof window.api.clinvarSubmissions === 'function') {
    try { return await window.api.clinvarSubmissions(variationId); } catch (e) { return { ok: false, error: e.message }; }
  }
  try { return { ok: true, rows: await clinvarSubmissions(variationId) }; } catch (e) { return { ok: false, error: e.message }; }
}

// A ClinVar card fed from NCBI (not GeneBe): the germline / somatic-oncogenicity /
// somatic-clinical-impact split, each with review stars, and a lazy per-submission table.
function buildClinvarCard(coord, gene) {
  const body = el('div', {}, el('div', { class: 'note' }, 'Loading ClinVar…'));
  const card = el('div', { class: 'card ve-card ve-clinvar' }, el('h2', { class: 've-title' }, 'ClinVar'), body);
  (async () => {
    const res = await clinvarRequest(coord, gene);
    body.innerHTML = '';
    if (!res || !res.ok) { body.append(el('div', { class: 'note warn' }, `ClinVar unavailable: ${(res && res.error) || 'no response'}.`)); return; }
    if (res.notFound) { body.append(el('div', { class: 'note' }, 'This allele is not in ClinVar.')); return; }
    // Which record these classifications belong to. A record is accepted either because its own
    // allele identifier (the canonical SPDI) is the change that was asked about, or — when it is the
    // only record at the position and carries no allele identifier at all — on the position alone.
    // Without the record's name on the card, the second case was indistinguishable from the first,
    // and a record describing a different (often structural) change at the same base was shown under
    // this variant's heading with nothing on screen to reveal it.
    if (res.title) body.append(el('div', { class: 'note muted', style: 'font-size:11.5px; margin-bottom:6px' },
      `Record: ${res.title}${res.accession ? ` · ${res.accession}` : ''}`));
    body.append(res.spdi
      ? el('div', { class: 'note muted', style: 'font-size:11.5px; margin-bottom:6px' },
        `Matched on this record's own allele: ${res.spdi}${res.cdna ? ` · ${res.cdna}` : ''}.`)
      : el('div', { class: 'note warn', style: 'font-size:11.5px; margin-bottom:6px' },
        'This record carries no allele identifier in ClinVar, so it was matched on the position alone. '
        + 'It was the only record at this base; check above that it describes the change you entered.'));
    // ClinVar files three separate classifications for the same variant — a germline one and two
    // somatic ones — and they answer different questions. As stacked lines they sat in a column down
    // the left of a full-width card. As three tiles they read side by side, each carrying its own
    // review status in words rather than only as stars a reader has to decode.
    const cvTile = (label, c) => el('div', { class: 've-tile cv-tile' },
      el('div', { class: 've-tile-label' }, label),
      c
        ? el('div', { class: 'cv-tile-body' },
          el('span', { class: `tag ${verdictClass(c.description)}` }, c.description),
          el('div', { class: 'cv-stars', title: c.reviewStatus }, starStr(c.stars)),
          c.reviewStatus ? el('div', { class: 've-tile-sub' }, c.reviewStatus) : null)
        : el('div', { class: 've-tile-sub', style: 'margin-top:6px' }, 'not classified for this variant'));
    body.append(el('div', { class: 've-tiles' },
      cvTile('Germline', res.germline),
      cvTile('Somatic — oncogenicity', res.oncogenicity),
      cvTile('Somatic — clinical impact', res.clinicalImpact)));
    const subsHolder = el('div', { style: 'margin-top:8px' });
    const btn = el('button', { class: 'btn small secondary' }, `Show ${res.submissionCount} submissions`);
    // Opening the submission list used to be one-way: the button removed itself, and the only way back
    // was to leave the tab. It toggles now, and is fetched once.
    let loaded = false; let shown = false;
    btn.addEventListener('click', async () => {
      if (loaded) {
        shown = !shown;
        subsHolder.style.display = shown ? '' : 'none';
        btn.textContent = shown ? `Hide ${res.submissionCount} submissions` : `Show ${res.submissionCount} submissions`;
        return;
      }
      btn.disabled = true; btn.textContent = 'Loading…';
      const out = await clinvarSubmissionRequest(res.variationId);
      btn.disabled = false;
      if (out && out.ok) {
        subsHolder.append(renderClinvarSubs(out.rows || []));
        loaded = true; shown = true;
        btn.textContent = `Hide ${res.submissionCount} submissions`;
      } else {
        subsHolder.append(el('div', { class: 'note warn' }, `Submissions unavailable: ${(out && out.error) || 'no response'}.`));
        btn.remove();
      }
    });
    body.append(el('div', { class: 'row', style: 'align-items:center; gap:10px; margin-top:8px' },
      res.submissionCount ? btn : null,
      el('a', { onclick: () => window.api.openExternal(res.url) }, 'Open in ClinVar')));
    body.append(subsHolder);
  })();
  return card;
}
function renderClinvarSubs(rows) {
  const table = el('table', { class: 'ref' },
    el('thead', {}, el('tr', {}, ['Submitter', 'Type', 'Classification', 'Review', 'Condition', 'Date'].map((h) => el('th', {}, h)))),
    el('tbody', {}, ...rows.slice(0, 80).map((r) => el('tr', {},
      el('td', { style: 'font-size:12px' }, r.submitter),
      el('td', { style: 'font-size:12px' }, r.kind),
      el('td', {}, el('span', { class: `tag ${verdictClass(r.value)}`, style: 'font-size:10px' }, r.value || '—')),
      el('td', { style: 'font-size:11px; letter-spacing:1px' }, starStr(r.stars)),
      el('td', { style: 'font-size:12px' }, r.condition),
      el('td', { style: 'font-size:12px' }, r.date)))));
  return el('div', { style: 'overflow-x:auto; margin-top:8px' }, table);
}

// =====================================================================
//  ASSESSMENT (redesign) — open public sources, no GeneBe engine
//  Pulls ClinVar + dbNSFP (MyVariant.info), gnomAD, and Cancer Hotspots +
//  Mutation Assessor (Genome Nexus), anchored to the app's own resolution,
//  plus one-click deep-links to the analysis portals. Works web + desktop
//  (all hosts are CORS-open and whitelisted in the CSP).
// =====================================================================
let assessState = { status: 'none', hgvs: null, data: null };

const enc = (s) => encodeURIComponent(String(s || ''));
const firstOf = (x) => (Array.isArray(x) ? x[0] : x);

/** A figure tile: small capitalised label, the number large, one line of reading under it. Matches the
 *  gnomAD tab's own tiles so the same fact is presented the same way wherever it appears. */
function veTile(label, big, sub) {
  return el('div', { class: 've-tile' },
    el('div', { class: 've-tile-label' }, label),
    el('div', { class: 've-tile-big' }, big),
    sub ? el('div', { class: 've-tile-sub' }, sub) : null);
}

// dbNSFP returns each predictor as an ARRAY with one entry per transcript, labelled by Ensembl
// identifier in `dbnsfp.ensembl.transcriptid`. Taking the first entry showed whichever transcript the
// source happened to list first, unlabelled — for KRAS c.35G>A the four AlphaMissense values are
// 0.999, 0.9846, 0.9984 and 0.6429, so the number on screen was a coin toss between them. Pick the
// entry for the transcript the user actually typed, using the RefSeq-to-Ensembl pairing from
// data/mane-pairing.json. The Ensembl identifier is a lookup key only and is never displayed.
let manePairing = null;
async function loadManePairing() {
  if (manePairing) return manePairing;
  try { manePairing = await window.api.getManePairing(); } catch { manePairing = { byRefSeq: {} }; }
  return manePairing;
}
/** The array index dbNSFP holds for `transcript`, or null when it cannot be identified. */
function dbnsfpIndexFor(dn, transcript) {
  const ids = dn && dn.ensembl && dn.ensembl.transcriptid;
  if (!Array.isArray(ids) || !transcript) return null;
  const pair = manePairing && manePairing.byRefSeq
    && (manePairing.byRefSeq[transcript] || manePairing.byRefSeq[String(transcript).replace(/\.\d+$/, '')]);
  if (!pair || !pair.ensembl) return null;
  const want = String(pair.ensembl).replace(/\.\d+$/, '');
  const i = ids.findIndex((id) => String(id).replace(/\.\d+$/, '') === want);
  return i >= 0 ? i : null;
}
/** One predictor value, for the named transcript where that can be established. */
function forTranscript(value, index) {
  if (!Array.isArray(value)) return value;
  if (index != null && index < value.length) return value[index];
  return null;   // rather than an unlabelled arbitrary element
}
function num3(v) { if (v == null || v === '' || v === '.') return null; const n = Number(v); return isFinite(n) ? n.toFixed(3) : null; }
function scoreCall(v, hi, lo) { if (v == null || v === '' || v === '.') return null; const n = Number(v); if (!isFinite(n)) return null; return n >= hi ? 'damaging' : (n <= lo ? 'benign' : 'uncertain'); }
function amWord(p) { if (!p) return null; const s = String(p).toLowerCase(); if (s === 'p' || s.includes('patho')) return 'likely pathogenic'; if (s === 'b' || s.includes('benign')) return 'likely benign'; return 'ambiguous'; }
function siftWord(p) { if (!p) return null; return String(p).toUpperCase().startsWith('D') ? 'deleterious' : 'tolerated'; }
function ppWord(p) { if (!p) return null; const c = String(p).toUpperCase()[0]; return c === 'D' ? 'probably damaging' : c === 'P' ? 'possibly damaging' : 'benign'; }

async function fetchJson(url, timeoutMs = 15000) {
  const ctrl = new AbortController();
  const t = setTimeout(() => ctrl.abort(), timeoutMs);
  try { const r = await fetch(url, { signal: ctrl.signal }); if (!r.ok) return null; return await r.json(); }
  catch { return null; } finally { clearTimeout(t); }
}
function fetchMyVariant(c) {
  return fetchJson(`https://myvariant.info/v1/variant/chr${c.chr}:g.${c.pos}${c.ref}%3E${c.alt}?assembly=hg38&fields=clinvar,dbnsfp,cadd,cosmic,civic,dbsnp`);
}
function fetchGenomeNexus(c) {
  // GRCh38 is a separate Genome Nexus host; the default (www) is GRCh37 and would reject our hg38 coord.
  return fetchJson(`https://grch38.genomenexus.org/annotation/${c.chr}:g.${c.pos}${c.ref}%3E${c.alt}?fields=hotspots,mutation_assessor,annotation_summary`);
}
// Coding HGVS -> GRCh38 plus-strand coordinate via Ensembl VEP. The whole tab hangs off this, so it is
// resolved on its own. Ensembl VEP is genuinely flaky (frequent transient 500s and multi-second
// latency), so retry a few times with backoff before giving up.
async function resolveGenomicCoord(hgvs) {
  const url = `https://rest.ensembl.org/vep/human/hgvs/${encodeURIComponent(hgvs)}?content-type=application/json`;
  for (let attempt = 0; attempt < 4; attempt += 1) {
    if (assessState.hgvs !== hgvs) return null; // superseded — stop retrying
    const data = await fetchJson(url, 25000);
    const r = data && data[0];
    if (r && r.seq_region_name && r.allele_string) {
      let [ref, alt] = String(r.allele_string).split('/');
      if (!ref || !alt || ref.length !== 1 || alt.length !== 1) return null; // parsed, but not a single-base sub — real "can't", don't retry
      if (Number(r.strand) === -1) { ref = NT_COMPLEMENT[ref] || ref; alt = NT_COMPLEMENT[alt] || alt; } // minus-strand -> plus
      const chr = String(r.seq_region_name);
      return { chr, chrom: chr, pos: r.start, ref, alt };
    }
    if (attempt < 3) await new Promise((res) => setTimeout(res, 1000 * (attempt + 1))); // 1s, 2s, 3s
  }
  return null;
}

function buildAssessment() {
  const root = $('#view-assessment'); root.innerHTML = '';
  root.append(
    // One line, not a card. What it said before was three paragraphs restating the panel headings
    // below it, and the same sentence appeared again at the foot of the page.
    el('div', { class: 've-standing' },
      el('span', { class: 've-warn' }, 'Sends the variant coordinate to ClinVar, gnomAD, MyVariant.info and Genome Nexus.'),
      ' Database lookups, not a clinical determination.'),
    el('div', { id: 'assess-panel' }));
  renderAssessment();
}

function renderAssessment() {
  const box = $('#assess-panel'); if (!box) return;
  box.innerHTML = '';
  const s = assessState.status;
  if (s === 'none') { box.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'Assess a variant, then open this tab to see its assessment.'))); return; }
  if (s === 'loading') { box.append(el('div', { class: 'card' }, el('div', { class: 'note' }, `Assessing ${assessState.hgvs}…`))); return; }
  if (s === 'error') { box.append(el('div', { class: 'card' }, el('div', { class: 'err' }, assessState.error || 'Assessment failed.'))); return; }
  renderAssessmentData(box, assessState.data);
}

async function openAssessment() {
  if (!lastVariant || !lastVariant.hgvs) { assessState = { status: 'none' }; renderAssessment(); return; }
  if (assessState.status === 'ok' && assessState.hgvs === lastVariant.hgvs) { renderAssessment(); return; }
  await pullAssessment(lastVariant.hgvs);
}

async function pullAssessment(hgvs) {
  assessState = { status: 'loading', hgvs }; renderAssessment();
  // 1) Coordinate first, on its own — reuse a prior known coord, else resolve via Ensembl VEP (retried).
  const lv = lastVariant;
  let coord = (lv && lv.coord && lv.coord.pos)
    ? { chr: lv.coord.chr || lv.coord.chrom, chrom: lv.coord.chrom || lv.coord.chr, pos: lv.coord.pos, ref: lv.coord.ref, alt: lv.coord.alt }
    : await resolveGenomicCoord(hgvs);
  if (assessState.hgvs !== hgvs) return; // superseded while awaiting
  // Fallback: if VEP could not place it, let the gnomAD path try (a separate VEP call) and reuse its result.
  let gPre = null;
  if (!coord || !coord.chr) {
    try { gPre = await window.api.gnomad(hgvs, null); } catch (e) { gPre = { ok: false, error: e.message }; }
    if (assessState.hgvs !== hgvs) return;
    const m = gPre && gPre.variantId && /^(\w+)-(\d+)-([A-Za-z]+)-([A-Za-z]+)$/.exec(String(gPre.variantId));
    if (m) coord = { chr: m[1], chrom: m[1], pos: Number(m[2]), ref: m[3].toUpperCase(), alt: m[4].toUpperCase() };
  }
  if (!coord || !coord.chr) { assessState = { status: 'error', hgvs, error: `Could not map ${hgvs} to a genomic coordinate — Ensembl VEP is not responding (it is frequently flaky). Reassess in a moment.` }; renderAssessment(); return; }
  if (lastVariant && lastVariant.hgvs === hgvs) lastVariant.coord = coord;
  // 2) The annotation sources in parallel — independent, so one failing does not sink the rest.
  const [g, mv, gn] = await Promise.all([
    gPre || window.api.gnomad(hgvs, { chrom: coord.chrom, pos: coord.pos, ref: coord.ref, alt: coord.alt }).catch((e) => ({ ok: false, error: e.message })),
    fetchMyVariant(coord),
    fetchGenomeNexus(coord),
    loadManePairing(),   // needed to pick each predictor score for THIS transcript, not an arbitrary one
  ]);
  if (assessState.hgvs !== hgvs) return;
  const askedTranscript = (parseHgvs(hgvs) || {}).transcript || null;
  assessState = { status: 'ok', hgvs, data: { coord: { ...coord, transcript: askedTranscript }, gnomad: g, myvariant: mv, genomenexus: gn } };
  renderAssessment();
}

const VARIANT_LEVEL_PORTALS = new Set(['VarSome', 'Franklin', 'dbSNP', 'GeneBe', 'Ensembl', 'UCSC']);
function assessLinkBtn(label, url) {
  // Whether a portal opens THIS variant or searches the whole gene matters, and it used to be said in
  // a paragraph under the buttons. Each button now carries it itself, so the card needs no footnote.
  const scope = VARIANT_LEVEL_PORTALS.has(label) ? 'opens this exact variant' : `searches the gene in ${label}`;
  const b = el('button', { class: 'btn small secondary', title: `${scope} — ${url}` }, `${label} ↗`);
  b.addEventListener('click', () => window.api.openExternal(url));
  return b;
}

function renderAssessmentData(box, d) {
  // The chromosome arrives under either name — `chr` from the coordinate resolver, `chrom` from the
  // gnomAD path and from the offline resolution of a non-coding variant. Reading one of them printed
  // "GRCh38 undefined-25245350-C-T" whenever the other was supplied, and the same mismatch sent
  // "chrundefined-…" to the SpliceAI service earlier. Normalise once, here, and use `c` below.
  const c = { ...d.coord, chr: d.coord.chr ?? d.coord.chrom, chrom: d.coord.chrom ?? d.coord.chr };
  const vid = `${c.chr}-${c.pos}-${c.ref}-${c.alt}`;
  // The gene and protein change come from the app's own resolution, which is only this assessment's
  // if it is still the same variant. Otherwise this panel would head one variant's annotation with
  // another variant's gene; the coordinate below it is always the assessed one.
  const lv = (lastVariant && lastVariant.hgvs === assessState.hgvs) ? lastVariant : null;
  const gene = (lv && lv.gene) || '';
  const prot = (lv && lv.appProtein) || '';
  const trx = (lv && lv.appTranscript) || '';

  const rsid = (d.gnomad && d.gnomad.rsids && d.gnomad.rsids[0])
    || (d.myvariant && d.myvariant.dbsnp && d.myvariant.dbsnp.rsid)
    || (d.myvariant && d.myvariant.clinvar && d.myvariant.clinvar.rsid) || null;
  const rsidChip = () => (rsid ? el('span', { class: 've-dim mono' }, rsid) : null);

  // Everything that identifies the variant on ONE line. This was a full-height card carrying a name at
  // the left and a coordinate at the right with nothing in between.
  const anchor = el('div', { class: 've-ident' },
    el('strong', {}, `${gene} ${prot}`.trim() || assessState.hgvs || 'Variant'),
    trx ? el('span', { class: 've-dim' }, trx) : null,
    el('span', { class: 've-dim mono' }, `GRCh38 ${vid}`),
    rsidChip(),
    el('span', { class: 've-dim', style: 'margin-left:auto' }, trx ? 'resolved offline' : 'not resolved offline'));

  // Clinical — ClinVar directly from NCBI: germline / somatic-oncogenicity / somatic-clinical-impact,
  // allele-exact, with review stars, the per-submission table and the exact variation link.
  const clinicalCard = buildClinvarCard(c, gene);


  // Population — gnomAD
  const g = d.gnomad || {};
  let popBody;
  if (g.notInGnomad) popBody = el('div', { class: 'note' }, 'Not observed in gnomAD v4 — consistent with a somatic or very rare variant.');
  else if (g.ok && g.overall) {
    // The same three tiles the gnomAD tab shows, so the two places agree in both figures and wording.
    const o = g.overall;
    popBody = el('div', { class: 've-tiles' },
      veTile('Overall frequency (WES + WGS)', pctFmt(o.af), afReading(o.af)),
      veTile('Frequency as a ratio (AN ÷ AC)', oneInN(o.ac, o.an),
        `${o.ac == null ? '—' : o.ac.toLocaleString()} in ${o.an == null ? '—' : o.an.toLocaleString()} alleles`),
      g.popmax ? veTile('Highest population', pctFmt(g.popmax.af), g.popmax.label) : null);
  }
  else popBody = el('div', { class: 'note warn' }, `gnomAD unavailable${g.error ? `: ${g.error}` : ''}.`);
  const popCard = el('div', { class: 'card ve-card ve-gnomad' }, el('h2', { class: 've-title' }, 'Population frequency (gnomAD v4)'), popBody,
    gnomadLink(g) ? el('div', { class: 've-cap' }, gnomadLink(g)) : null);

  // Somatic — Genome Nexus (Cancer Hotspots + COSMIC colocated + Mutation Assessor)
  const gn = d.genomenexus || {};
  const hotAnn = gn.hotspots && gn.hotspots.annotation;
  const isHot = Array.isArray(hotAnn) && [].concat(...hotAnn).filter(Boolean).length > 0;
  const cosmicIds = (gn.colocatedVariants || []).map((v) => v && v.dbSnpId).filter((x) => /^COSV/i.test(String(x)));
  const ma = gn.mutation_assessor && (gn.mutation_assessor.functionalImpactPrediction || gn.mutation_assessor.functionalImpact);
  const somaticCard = el('div', { class: 'card ve-card ve-somatic' }, el('h2', { class: 've-title' }, 'Somatic context'),
    el('div', {}, el('strong', {}, 'Cancer hotspot: '), isHot ? el('span', { class: 'tag v-path' }, 'known hotspot residue') : el('span', { class: 'note muted' }, 'not annotated here — check cBioPortal / OncoKB')),
    cosmicIds.length ? el('div', { style: 'margin-top:6px' }, el('strong', {}, 'COSMIC: '), el('span', { class: 'mono', style: 'font-size:12px' }, cosmicIds.slice(0, 3).join(', '))) : null,
    ma ? el('div', { style: 'margin-top:6px' }, el('strong', {}, 'Mutation Assessor: '), el('span', {}, String(ma))) : null,
    // The gene-level links are only offered when a gene name is actually in hand. A link built from
    // an empty gene goes to a search for nothing, which looks like an answer.
    el('div', { class: 've-links', style: 'margin-top:10px' },
      gene ? assessLinkBtn('OncoKB', `https://www.oncokb.org/gene/${enc(gene)}/${enc(prot)}`) : null,
      gene ? assessLinkBtn('cBioPortal', `https://www.cbioportal.org/results/mutations?gene_list=${enc(gene)}&cancer_study_list=all&tab_index=tab_visualize&Action=Submit`) : null,
      (cosmicIds.length || gene) ? assessLinkBtn('COSMIC', cosmicIds.length ? `https://cancer.sanger.ac.uk/cosmic/search?q=${enc(cosmicIds[0])}` : `https://cancer.sanger.ac.uk/cosmic/search?q=${enc(gene)}`) : null,
      gene ? assessLinkBtn('CIViC', `https://civicdb.org/variants/home?gene=${enc(gene)}`) : null));

  // In-silico — dbNSFP via MyVariant.info
  const dn = (d.myvariant && d.myvariant.dbnsfp) || {};
  const tx = (d.coord && d.coord.transcript) || (lastVariant && lastVariant.appTranscript) || null;
  const ix = dbnsfpIndexFor(dn, tx);
  const v = (x) => forTranscript(x, ix);
  // CADD is genomic rather than per-transcript, so it has one value whichever transcript is asked about
  const caddPhred = (d.myvariant && d.myvariant.cadd && firstOf(d.myvariant.cadd.phred)) ?? firstOf(dn.cadd && dn.cadd.phred);
  const tiles = [
    predTile('REVEL', num3(v(dn.revel && dn.revel.score)), scoreCall(v(dn.revel && dn.revel.score), 0.7, 0.29), { range: [0, 1] }),
    predTile('AlphaMissense', num3(v(dn.alphamissense && dn.alphamissense.score)), amWord(v(dn.alphamissense && dn.alphamissense.pred)), { range: [0, 1] }),
    predTile('CADD (phred)', caddPhred != null ? Number(caddPhred).toFixed(1) : null, caddPhred != null ? (Number(caddPhred) >= 20 ? 'deleterious' : 'tolerated') : null),
    predTile('SIFT', num3(v(dn.sift && (dn.sift.converted_rankscore || dn.sift.score))), siftWord(v(dn.sift && dn.sift.pred))),
    predTile('PolyPhen-2', num3(v(dn.polyphen2 && dn.polyphen2.hdiv && dn.polyphen2.hdiv.score)), ppWord(v(dn.polyphen2 && dn.polyphen2.hdiv && dn.polyphen2.hdiv.pred))),
    predTile('PrimateAI', num3(v(dn.primateai && dn.primateai.score)), siftWord(v(dn.primateai && dn.primateai.pred))),
  ].filter(Boolean);
  const insilicoCard = el('div', { class: 'card ve-tight ve-predbar' }, el('h2', {}, 'In-silico predictions'),
    tiles.length ? el('div', { class: 'row', style: 'gap:10px; flex-wrap:wrap' }, ...tiles) : el('div', { class: 'note' }, 'No in-silico scores returned for this position.'),
    el('div', { class: 'note muted', style: 'margin-top:8px; font-size:11px' },
      ix != null && tx
        ? `dbNSFP via MyVariant.info · scores are for ${tx}. Other tools often show the highest score across all transcripts of the gene instead, which can differ.`
        : 'dbNSFP via MyVariant.info · the transcript these scores belong to could not be identified, so any that vary between transcripts are not shown.'));

  const portalBtns = [
    assessLinkBtn('VarSome', `https://varsome.com/variant/hg38/${vid}`),
    assessLinkBtn('Franklin', `https://franklin.genoox.com/clinical-db/variant/snp/chr${vid}-hg38`),
    rsid ? assessLinkBtn('dbSNP', `https://www.ncbi.nlm.nih.gov/snp/${enc(rsid)}`) : null,
    // GeneBe is the annotation engine this app already queries; the link opens the same variant on its
    // own site, where the full ACMG working is laid out.
    assessLinkBtn('GeneBe', `https://genebe.net/variant/hg38/${vid}`),
    // Ensembl's variant page (Variation/Explore) returns a server-side "Runtime Error … malformed JSON"
    // for ordinary rsIDs — checked 8 Aug 2026 on rs121913529 — so this points at the region view, which
    // works, and does not need an rsID either.
    assessLinkBtn('Ensembl', `https://www.ensembl.org/Homo_sapiens/Location/View?r=${enc(`${c.chr}:${c.pos}-${c.pos}`)}`),
    assessLinkBtn('UCSC', `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr${c.chr}%3A${c.pos}-${c.pos}`),
    // Gene-level searches, so they are only offered when the gene name is known.
    // LOVD and Mastermind were removed on 8 Aug 2026: LOVD blocks browser traffic from this network
    // outright ("your IP address or network has been blacklisted"), and Mastermind puts everything
    // behind a sign-up. A link that cannot be followed is worse than no link.
    gene ? assessLinkBtn('PubMed', `https://pubmed.ncbi.nlm.nih.gov/?term=${enc(`${gene} ${prot}`)}`) : null,
  ].filter(Boolean);
  const portals = el('div', { class: 'card ve-card ve-portals' }, el('h2', { class: 've-title' }, 'Analysis portals'),
    el('div', { class: 've-links' }, ...portalBtns));

  // Order follows how the panel is read: what the variant is, then the scores at a glance, then the
  // evidence behind them. gnomAD and the somatic context are short, so they share a column beside
  // ClinVar, which is the tall one — that is where the empty space was coming from.
  box.append(anchor,
    insilicoCard,
    // Balanced columns rather than two fixed ones: whichever of these is tallest varies by variant —
    // ClinVar is long when a variant has many submissions and two lines when it has none — so a fixed
    // grid leaves a ragged bottom under one column or the other. Letting the browser distribute them
    // keeps both columns ending at about the same place whatever the content turns out to be.
    // Three short boxes across, then ClinVar last and full width. ClinVar is the one that grows — its
    // submission list runs to dozens of rows — so it goes at the bottom, where opening it pushes
    // nothing else off the screen.
    el('div', { class: 've-three' }, popCard, somaticCard, portals),
    clinicalCard,
    el('div', { class: 've-foot' }, el('div', { class: 'note muted', style: 'font-size:11px' },
      'Live lookups: ClinVar & dbNSFP via MyVariant.info; gnomAD (Broad); Cancer Hotspots & Mutation Assessor via Genome Nexus. The variant coordinate is sent to these public services; results are database lookups, not a clinical determination.')));
}

function renderGeneBeAnnotation(box, data) {
  const a = data.annotation; const gp = data.genomic;
  const byGene = (a.acmg_by_gene && a.acmg_by_gene[0]) || {};
  const hgvsc = byGene.hgvs_c || ''; const hgvsp = byGene.hgvs_p || '';

  box.append(el('div', { class: 'card' },
    el('div', { class: 'row', style: 'justify-content:space-between; align-items:flex-start; gap:16px' },
      el('div', {},
        el('h2', {}, `${a.gene_symbol || byGene.gene_symbol || ''} ${hgvsp}`.trim()),
        el('div', { class: 'sub' }, `${a.transcript || byGene.transcript || ''} ${hgvsc} · ${a.effect || ''} · ${a.dbsnp || 'no dbSNP id'}`)),
      el('div', { class: `consequence ${verdictClass(a.acmg_classification)}` }, el('span', { class: 'dot' }), `ACMG: ${a.acmg_classification || '—'}`)),
    el('div', { class: 'note', style: 'margin-top:8px' }, `Genomic (hg38): chr${gp.chr}:${gp.pos} ${gp.ref}>${gp.alt}${data.cached ? ' · from local cache' : ''}`)));

  // Trust check: does GeneBe's protein change agree with this app's own offline resolution?
  // Only when both are the same variant. Comparing GeneBe's answer for one variant against the app's
  // resolution of another would print "GeneBe differs from the reference" about two different things.
  const sameVariant = Boolean(lastVariant && latestGeneBe.hgvs && lastVariant.hgvs === latestGeneBe.hgvs);
  const appP = sameVariant && lastVariant.appProtein;
  const isSub = /^[A-Z]\d+[A-Z]$/.test(appP || '');
  const gbP = geneBeProteinChange(a, lastVariant && lastVariant.appTranscript);
  if (appP && isSub && gbP) {
    const agree = appP.toUpperCase() === gbP.text.toUpperCase();
    box.append(el('div', { class: 'card', style: `border-color:${agree ? '#2f6d4f' : '#8a6a2a'}` },
      el('div', { class: 'row', style: 'align-items:center; gap:10px' },
        el('span', { style: `font-size:18px; color:${agree ? '#46b877' : '#e0a54e'}` }, agree ? '✓' : '⚠'),
        el('div', {}, agree
          ? el('div', {}, el('strong', {}, 'GeneBe agrees with the reference. '), el('span', { class: 'note' }, `Both resolve this to ${appP} on ${lastVariant.appTranscript}.`))
          : el('div', {}, el('strong', {}, 'GeneBe differs from the reference. '), el('span', { class: 'note' }, `This app resolves ${appP} on ${lastVariant.appTranscript}; GeneBe reports ${gbP.text}${gbP.transcript ? ' on ' + gbP.transcript : ''}. Check the transcript version.`))))));
  }

  const acmgCard = el('div', { class: 'card' }, el('h2', {}, 'ACMG classification (GeneBe)'),
    el('div', { class: 'kv' },
      el('div', { class: 'k' }, 'Classification'), el('div', { class: 'v' }, a.acmg_classification || '—'),
      el('div', { class: 'k' }, 'Score'), el('div', { class: 'v' }, String(a.acmg_score ?? '—')),
      el('div', { class: 'k' }, 'Combined'), el('div', { class: 'v' }, a.pathogenicity_classification_combined || '—')),
    el('div', { class: 'sub', style: 'margin:10px 0 6px' }, 'Criteria applied'),
    el('div', { class: 'chips' }, String(a.acmg_criteria || '').split(',').map((c) => c.trim()).filter(Boolean).map((c) => el('span', { class: 'chip' }, c)) || el('span', { class: 'note' }, 'none')));

  const clinvarCard = buildClinvarCard({ chr: gp.chr, pos: gp.pos, ref: gp.ref, alt: gp.alt }, a.gene_symbol || byGene.gene_symbol);

  const gnomadCard = el('div', { class: 'card' }, el('h2', {}, 'Population frequency (gnomAD)'),
    el('div', { class: 'kv' },
      el('div', { class: 'k' }, 'Exomes AF'), el('div', { class: 'v' }, fmtAf(a.gnomad_exomes_af)),
      el('div', { class: 'k' }, 'Exomes AC / hom'), el('div', { class: 'v' }, `${a.gnomad_exomes_ac ?? '—'} / ${a.gnomad_exomes_homalt ?? '—'}`),
      el('div', { class: 'k' }, 'Genomes AF'), el('div', { class: 'v' }, fmtAf(a.gnomad_genomes_af)),
      el('div', { class: 'k' }, 'Genomes AC / hom'), el('div', { class: 'v' }, `${a.gnomad_genomes_ac ?? '—'} / ${a.gnomad_genomes_homalt ?? '—'}`)),
    el('div', { class: 'note muted', style: 'margin-top:8px' }, 'Absent from gnomAD is shown as “—”. Rarity supports pathogenicity; a common allele argues against it.'));

  // In-silico predictions, grouped into a tile grid rather than one flat list.
  const missenseTiles = [
    predTile('AlphaMissense', a.alphamissense_score, a.alphamissense_prediction),
    predTile('REVEL', a.revel_score, a.revel_prediction),
    predTile('BayesDel', a.bayesdelnoaf_score, a.bayesdelnoaf_prediction, { range: [-1.3, 0.7] }),
    predTile('GeneBe consensus', typeof a.computational_score_selected === 'number' ? a.computational_score_selected.toFixed(3) : a.computational_score_selected, a.computational_prediction_selected),
  ].filter(Boolean);
  const conservTiles = [
    predTile('phyloP 100-way', a.phylop100way_score, a.phylop100way_prediction, { range: [-20, 10] }),
  ].filter(Boolean);
  const spliceTiles = [
    predTile('SpliceAI (max)', a.spliceai_max_score, a.spliceai_max_prediction),
    predTile('dbscSNV ADA', a.dbscsnv_ada_score, a.dbscsnv_ada_prediction),
  ].filter(Boolean);
  const mitoTiles = [
    predTile('MitoTip', a.mitotip_score, a.mitotip_prediction),
    predTile('APOGEE2', a.apogee2_score, a.apogee2_prediction),
  ].filter(Boolean);

  const group = (title, tiles) => tiles.length ? el('div', { style: 'margin-top:8px' },
    el('div', { class: 'sub', style: 'margin-bottom:6px' }, title),
    el('div', { class: 'row', style: 'gap:8px; flex-wrap:wrap' }, ...tiles)) : null;

  const predCard = el('div', { class: 'card' }, el('h2', {}, 'In-silico predictions'),
    group('Missense pathogenicity', missenseTiles),
    group('Conservation', conservTiles),
    group('Splice', spliceTiles),
    group('Mitochondrial', mitoTiles),
    el('div', { class: 'note muted', style: 'margin-top:10px; font-size:11.5px' }, 'For detailed SpliceAI donor/acceptor deltas and validated splice evidence, see the Splice tab.'));

  // RefSeq transcripts only (NM_/NR_/XM_/XR_): the table is titled RefSeq, yet GeneBe also returns the
  // Ensembl (ENST) rows and flags the ENST as MANE. MANE Select is the SAME transcript under a paired
  // NM_ and ENST accession, so we drop the ENST rows and put the MANE tick on the RefSeq one (its
  // mane_select field points to the paired ENST). The MANE Select transcript is listed first.
  const cons = (a.consequences || [])
    .filter((c) => c.feature && /^(NM_|NR_|XM_|XR_)/.test(c.feature))
    .sort((x, y) => ((y.mane_select ? 2 : 0) + (y.mane_plus ? 1 : 0)) - ((x.mane_select ? 2 : 0) + (x.mane_plus ? 1 : 0)));
  let transcriptCard = null;
  if (cons.length) {
    const rows = cons.slice(0, 20).map((c) => el('tr', {},
      el('td', { class: 'mono', style: 'font-size:12px' }, c.feature),
      el('td', { class: 'mono', style: 'font-size:12px' }, c.hgvs_c || ''),
      el('td', { class: 'mono', style: 'font-size:12px' }, c.hgvs_p || (c.aa_ref && c.aa_start ? `p.${c.aa_ref}${c.aa_start}${c.aa_alt || ''}` : '')),
      el('td', {}, (c.consequences || []).join(', ')),
      el('td', {}, c.exon_rank != null ? `${c.exon_rank}${c.exon_count ? '/' + c.exon_count : ''}` : ''),
      el('td', { title: c.mane_plus ? 'MANE Plus Clinical' : (c.mane_select ? 'MANE Select' : '') }, c.mane_select ? '✓' : (c.mane_plus ? '✓ Plus' : ''))));
    transcriptCard = el('div', { class: 'card' }, el('h2', {}, `RefSeq transcripts (${cons.length})`),
      el('div', { style: 'overflow-x:auto' }, el('table', { class: 'ref' },
        el('thead', {}, el('tr', {}, ['Transcript', 'HGVSc', 'HGVSp', 'Consequence', 'Exon', 'MANE'].map((h) => el('th', {}, h)))),
        el('tbody', {}, ...rows))));
  }

  box.append(el('div', { class: 'grid2' }, acmgCard, clinvarCard), el('div', { class: 'grid2' }, gnomadCard, predCard));
  if (transcriptCard) box.append(transcriptCard);
  box.append(el('div', { class: 'card' }, el('div', { class: 'note muted' },
    'Annotation by GeneBe — ', el('a', { onclick: () => window.api.openExternal('https://genebe.net') }, 'genebe.net'),
    '. The variant was sent to GeneBe’s servers. Population frequency, ClinVar and the ACMG call are GeneBe’s automated computation, not a clinical determination.')));
}

// =====================================================================
//  FUSIONS / TRANSLOCATIONS (BND) — visual breakend assessment
// =====================================================================
let fusionGeneLoci = null;  // {SYMBOL: {chr,strand,start,end,tx}} — lazy-loaded MANE gene spans
let fusionTables = null;    // { exons, cytobands, igLoci, igGenes } — lazy-loaded reference tables

// The drawing context and classifier options, built from whichever tables loaded. Everything
// degrades: with a table missing the figures lose that layer and the core falls back to span-only.
function fusionCtx() {
  const T = fusionTables || {};
  return { genes: fusionGeneLoci || {}, exons: T.exons, cytobands: T.cytobands, igLoci: T.igLoci, igGenes: T.igGenes };
}
function fusionCoreOpts() {
  const T = fusionTables || {};
  return { exons: (T.exons && T.exons.genes) || null, igLoci: (T.igLoci && T.igLoci.loci) || null };
}

/**
 * Pair each canonical fusion with its own breakend row.
 *
 * The picker used to take the row at the same position in the file as the expectation, and the
 * generated file's own control could not see the difference: `parseReport` skips a blank line
 * without recording it, so one stray blank line kept every count right while shifting every row
 * after it by one — the picker then loaded one fusion's junction under the next one's name and
 * disease. Each row is matched instead by the fusion name written in its own annotation column, so
 * the pairing depends on what the row says, not on where it sits.
 *
 * A name that matches no row, or more than one, yields no option at all: better a missing entry the
 * caller can report than an entry that loads the wrong junction.
 */
function pairFusionExamples(examplesText, expectations) {
  const rows = String(examplesText || '').split('\n').filter((r) => r.trim());
  const options = []; const unmatched = [];
  for (const e of expectations || []) {
    const hits = rows.filter((row) => String(row.split('\t')[3] || '').trim().startsWith(`${e.name} `));
    if (hits.length === 1) options.push({ name: e.name, disease: e.disease, row: hits[0] });
    else unmatched.push(e.name);
  }
  return { options, unmatched };
}

function buildFusions() {
  const root = $('#view-fusions'); root.innerHTML = '';
  const input = el('textarea', { placeholder: 'Paste BND rows from the SV report — one per line, e.g.\n3:187746165T]12:25056225]\t3q27.3;12p12.1\tBND\tN/A\t-\tPR-40/194;SR-23/159', style: 'min-height:120px; font-family:var(--mono); font-size:12px' });
  const out = el('div', { style: 'margin-top:14px' });
  const analyse = async () => {
    out.innerHTML = ''; out.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'Analysing…')));
    // Every gene named at a breakpoint comes from this table. If it does not load, the classifier
    // has nothing to look in and reports every breakpoint as "no MANE gene within 2 Mb" — a
    // statement about the genome standing in for a file that did not arrive. Nothing is analysed
    // until the table is there, so that sentence can never be read as a finding.
    let loadError = null;
    if (!fusionGeneLoci) {
      try {
        const g = (await window.api.getGeneLoci()).genes;
        if (g && Object.keys(g).length) fusionGeneLoci = g;
        else loadError = 'the gene table came back empty';
      } catch (e) { loadError = e.message || String(e); }
    }
    if (!fusionTables) { try { fusionTables = await window.api.getFusionTables(); } catch { fusionTables = {}; } }
    out.innerHTML = '';
    if (loadError) {
      const retry = el('button', { class: 'btn small', style: 'margin-top:10px', onclick: analyse }, 'Try again');
      out.append(el('div', { class: 'card', style: 'border-color:#6a2b30' },
        el('strong', { class: 'err' }, 'The reference gene table did not load'),
        el('div', { class: 'note muted', style: 'font-size:11.5px; margin-top:4px' },
          `Breakends cannot be analysed without it — no gene could be named at any breakpoint, and an empty table would make every breakpoint look like empty genome. Reason: ${loadError}.`),
        el('div', { class: 'row' }, retry)));
      return;
    }
    out.append(renderFusionEvents(input.value));
  };
  root.append(
    el('div', { class: 'card' },
      el('h2', {}, 'Fusions / translocations (BND)'),
      el('div', { class: 'sub' }, 'Paste breakend (BND) rows from the structural-variant report. Each junction is read in transcript direction — the gene that donates the promoter first, the gene that donates the body second — and reciprocal pairs are matched automatically. Fully offline; nothing leaves this machine.'),
      input,
      el('div', { class: 'row', style: 'gap:8px; margin-top:10px; align-items:center' },
        el('button', { class: 'btn', onclick: analyse }, 'Analyse breakends'),
        (() => {
          // Pick one canonical fusion at a time; the button loads the whole set.
          const { options, unmatched } = pairFusionExamples(FUSION_HAEM_EXAMPLES, FUSION_HAEM_EXPECTATIONS);
          const picker = el('select', { class: 'fus-picker', title: 'Canonical haematological fusions, constructed from the bundled reference tables. Not patient data.' },
            el('option', { value: '' }, 'Canonical fusion…'),
            ...options.map((o) => el('option', { value: o.row }, `${o.name} — ${o.disease}`)));
          picker.onchange = () => {
            if (picker.value === '') return;
            // Selecting a fusion shows that fusion alone.
            input.value = picker.value;
            picker.value = '';
            analyse();
          };
          // If the example set and the expectations ever drift apart, say which entries are missing
          // rather than offering an option that would load the wrong junction.
          return unmatched.length
            ? el('div', { class: 'row', style: 'gap:8px; align-items:center' }, picker,
              el('span', { class: 'note warn', style: 'font-size:11px' },
                `${unmatched.length} example${unmatched.length === 1 ? '' : 's'} could not be matched to a row and ${unmatched.length === 1 ? 'is' : 'are'} not offered: ${unmatched.join(', ')}.`))
            : picker;
        })(),
        el('button', { class: 'btn small secondary', onclick: () => { input.value = ''; out.innerHTML = ''; } }, 'Clear')),
      el('div', { class: 'note muted', style: 'margin-top:10px; font-size:11px' },
        'Examples are constructed from the bundled reference tables and verified in frame where a fusion protein is expected — not patient data. Rows without a cytoband column still read correctly, but the derivative goes unnamed.')),
    out);
}

// The product kind as the core named it, with the hyphens opened out. Formatting only — the words are
// the classifier's own (`product.kind`), never a gloss written here.
// Most kind slugs read as plain English once the hyphens go. The ones that do not get a phrase here.
const FUSION_KIND_LABEL = {
  // The product is promoter substitution either way; which side donates the promoter is already in
  // the reading column, and the productivity column carries the "possible" that the unnamed donor earns.
  'promoter-substitution-unnamed-donor': 'promoter substitution',
};
function fusionKindLabel(kind) {
  return FUSION_KIND_LABEL[kind] || String(kind || '').replace(/-/g, ' ');
}

// How a piece reads: gene plus the core's short part label; an IG/TCR locus by its name; else the
// bare chromosome. `fusionBareName` is the ledger variant — the name alone.
function fusionPieceLabel(seg) {
  if (seg.content.gene) return `${seg.content.gene} ${partShort(seg.content.part)}`;
  if (seg.content.igLocus) return `${seg.content.igLocus.symbol} locus`;
  return `chr${seg.chr}`;
}
function fusionBareName(seg) {
  return seg.content.gene || (seg.content.igLocus ? `${seg.content.igLocus.symbol} locus` : `chr${seg.chr}`);
}

// `short` is for the ledger, where the column is narrow; `word` is for the panel, which has the room.
const FUSION_VERDICT = {
  true: { word: 'productive junction', short: 'yes', cls: 'v-driver' },
  possible: { word: 'possibly productive', short: 'possible', cls: 'v-maybe' },
  false: { word: 'not productive', short: 'no', cls: 'v-unknown' },
};
function fusionVerdict(productive) { return FUSION_VERDICT[String(productive)] || FUSION_VERDICT.false; }

const FUSION_VERDICT_HELP = 'Whether a transcript can read across this junction, or a gene\'s promoter is changed by it.';

const FUSION_LEDGER_COLUMNS = [
  ['Derivative', 'The rearranged chromosome this junction builds. Select it to jump to the full reading.'],
  ['Reads 5′ → 3′', 'The two pieces in the order a transcript would read them, which is not always the order the report writes them.'],
  ['Product', 'What the junction makes, if anything.'],
  ['Productivity*', FUSION_VERDICT_HELP],
  ['Read pairs', 'Read pairs whose two ends land on opposite sides of the junction, exactly as the caller reported them.'],
  ['Split reads', 'Single reads crossing the breakpoint itself, part aligning to each side, exactly as the caller reported them.'],
];

/** Summary table: every junction in the paste, one row, productive first. Each row links to its panel. */
function drawFusionLedger(rows) {
  const head = el('tr', {}, ...FUSION_LEDGER_COLUMNS.map(([h, help], i) => el('th', { class: i > 3 ? 'num' : null, title: help }, h)));
  const body = rows.map(({ d, anchor }) => {
    const v = fusionVerdict(d.product.productive);
    const sup = (s) => (s ? `${s.alt} of ${s.total}` : '—');
    const jump = (e) => {
      e.preventDefault();
      const t = document.getElementById(anchor);
      if (!t) return;
      t.scrollIntoView({ block: 'center', behavior: 'smooth' });
      t.classList.add('flash');
      setTimeout(() => t.classList.remove('flash'), 1400);
    };
    return el('tr', { class: d.product.productive === true ? 'productive' : null },
      el('td', {}, el('a', { href: `#${anchor}`, onclick: jump }, d.derivative.name || `Junction ${d.index + 1}`)),
      el('td', {}, el('span', { class: 'g' }, fusionBareName(d.reading[0])),
        el('span', { class: 'arrow' }, d.product.productive === false ? '✕' : '→'),
        el('span', { class: 'g' }, fusionBareName(d.reading[1]))),
      el('td', { class: 'kind' }, fusionKindLabel(d.product.kind)),
      el('td', {}, el('span', { class: `tag ${v.cls}`, title: `${v.word} — ${fusionKindLabel(d.product.kind)}` }, v.short)),
      el('td', { class: 'num' }, sup(d.support && d.support.pr)),
      el('td', { class: 'num' }, sup(d.support && d.support.sr)));
  });
  // The verdict definitions live in a thesis-style footnote under the table, outside the card.
  return [
    el('div', { class: 'card' },
      el('div', { class: 'row', style: 'justify-content:space-between; align-items:baseline; gap:10px; flex-wrap:wrap' },
        el('h2', { style: 'margin:0; font-size:15px' }, 'Summary'),
        el('span', { class: 'note muted', style: 'font-size:11px' }, 'most consequential first · select a derivative to jump to it')),
      el('div', { class: 'track-wrap', style: 'margin-top:8px' },
        el('table', { class: 'fus-ledger' }, el('thead', {}, head), el('tbody', {}, ...body)))),
    el('p', { class: 'fus-tablenote' },
      '* Productivity — ', el('b', {}, 'yes'), ': a transcript reads across the junction, or an intact gene gains a new promoter or enhancer. ',
      el('b', {}, 'possible'), ': the arrangement allows a product; confirmation needs data beyond this report. ',
      el('b', {}, 'no'), ': no transcript can cross the junction.'),
  ];
}

function renderFusionEvents(text) {
  const breakends = parseReport(text);
  attachGermlineText(text, breakends);
  const skipped = breakends.unparsed || [];
  const frag = el('div', {});
  // Rows that were not understood are named. Silence here would read as "nothing to report".
  if (skipped.length) frag.append(el('div', { class: 'card', style: 'border-color:#7a4b28' },
    el('strong', { style: 'color:#f0c07a' }, `${skipped.length} row${skipped.length === 1 ? '' : 's'} not recognised as a breakend`),
    el('div', { class: 'note muted', style: 'font-size:11.5px; margin-top:4px' },
      'These were skipped. Copy-number and intra-chromosomal deletion/duplication rows are out of scope here; anything else is worth checking for a typo or a truncated cell.'),
    el('ul', { class: 'fus-caveats' }, ...skipped.slice(0, 12).map((s) => el('li', { class: 'mono' }, s))),
    skipped.length > 12 ? el('div', { class: 'note muted', style: 'font-size:11px' }, `…and ${skipped.length - 12} more.`) : null));

  if (!breakends.length) {
    frag.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'No BND rows recognised. Paste breakend rows such as 3:187746165T]12:25056225].')));
    return frag;
  }
  const { events } = detectEvents(breakends);
  const recip = events.filter((e) => e.reciprocal).length;
  const junctions = events.reduce((n, e) => n + e.breakends.length, 0);
  const collapsed = breakends.length - junctions;
  frag.append(el('div', { class: 'sub', style: 'margin-bottom:12px' },
    `${breakends.length} row${breakends.length === 1 ? '' : 's'} → ${junctions} junction${junctions === 1 ? '' : 's'} in ${events.length} event${events.length === 1 ? '' : 's'}${recip ? `, ${recip} reciprocal` : ''}.`
    + (collapsed > 0 ? ` ${collapsed} row${collapsed === 1 ? ' was' : 's were'} the same junction written from its other end and ${collapsed === 1 ? 'has' : 'have'} been merged.` : '')));
  // The ledger goes first: the whole paste at a glance, before any one event is opened up. Built from
  // the same describeDerivatives() call the panels use, so the two can never disagree.
  const ledgerRows = [];
  events.forEach((ev, ei) => {
    try {
      describeDerivatives(ev, fusionGeneLoci || {}, fusionCoreOpts()).forEach((d, di) => ledgerRows.push({ d, anchor: `fus-${ei}-${di}` }));
    } catch { /* the per-event catch below reports it */ }
  });
  if (ledgerRows.length) {
    const rank = (r) => (r.d.product.productive === true ? 0 : (r.d.product.productive === 'possible' ? 1 : 2));
    frag.append(...drawFusionLedger(ledgerRows.slice().sort((p, q) => rank(p) - rank(q))));
  }

  events.forEach((ev, ei) => {
    // One malformed event must not blank the whole panel with no explanation.
    try {
      frag.append(drawFusionEvent(ev, ei));
    } catch (err) {
      frag.append(el('div', { class: 'card', style: 'border-color:#6a2b30' },
        el('strong', { class: 'err' }, 'This event could not be drawn'),
        el('div', { class: 'note muted', style: 'font-size:11.5px; margin-top:4px' }, String((err && err.message) || err)),
        el('div', { class: 'mono', style: 'font-size:11px; margin-top:6px' }, ev.breakends.map((b) => b.raw).join('  ·  '))));
    }
  });
  return frag;
}

function fusionChip(label, value, cls, title) {
  return el('span', { class: `fus-chip${cls ? ` ${cls}` : ''}`, title: title || null }, el('b', {}, label), value);
}

// A number in the frequency column only counts as a population frequency when the report tagged it
// as one. `parseGermlineAf` marks a bare, untagged decimal `assumed`, because a somatic
// variant-allele fraction — the share of reads carrying the change — is written identically. An
// assumed number is shown as the unlabelled number it is and never earns the germline verdict.
function germlineIsStated(g) { return Boolean(g && g.populationDb && g.assumed !== true); }

// How an untagged frequency column is described wherever it appears, so the wording cannot drift.
const GERMLINE_ASSUMED_HELP = 'The report does not say what this column is. A germline population '
  + 'frequency and a somatic variant-allele fraction (the share of reads carrying the change) are '
  + 'written the same way here, so the number is left as an unlabelled number and counted as neither.';
function germlineAssumedText(g) {
  const others = g && g.otherCandidates ? g.otherCandidates : 0;
  // `asRead` is the field exactly as the row wrote it. Falling back to the parsed number would show
  // a row that wrote "1E-4" as 0.0001 — the same value, but not what the report said.
  return `${g.asRead || g.af} — column not labelled`
    + (others ? `; ${others} other unlabelled number${others === 1 ? '' : 's'} in the row, the first was taken` : '');
}

// The frequency field as the row wrote it, character for character. `parseGermlineAf` returns the
// value as a number, so the characters are gone by the time the panel is drawn. The row is split
// into fields and each field is offered back to that same parser, so the token found here is the one
// the parser read — there is no second copy of the column grammar in this file. If the row holds two
// different spellings of the same value, or none can be matched, nothing is claimed and the panel
// falls back to the parsed number.
function germlineAsRead(row, g) {
  if (!row || !g) return null;
  const hits = new Set();
  for (const field of String(row).split(/\t|\s{2,}/)) {
    const token = field.trim();
    const p = parseGermlineAf(token);
    if (p && p.af === g.af && Boolean(p.populationDb) === Boolean(g.populationDb)) hits.add(token);
  }
  return hits.size === 1 ? [...hits][0] : null;
}

// Pair each breakend with the row it came from, so the frequency field can be shown as written.
// `parseReport` takes one line at a time and yields at most one breakend per line, so running it
// line by line says exactly which row produced which breakend without re-reading the row grammar
// here. The germline record is the same object the derivative carries, so noting the text on it
// reaches the panel.
function attachGermlineText(text, breakends) {
  let i = 0;
  for (const line of String(text || '').split(/\r?\n/)) {
    if (!line.trim()) continue;
    if (!parseReport(line).length) continue;      // a row that yields no breakend
    const b = breakends[i]; i += 1;
    if (b && b.germline) b.germline.asRead = germlineAsRead(line, b.germline);
  }
}

// Chromosome order for display: numeric, with X and Y last, so t(3;12) is never shown as 12/3.
// Two names for the same chromosome must compare equal, or a caller that tie-breaks on position
// (`fusionChrOrder(a, b) || a.pos - b.pos`) never reaches its tie-break and the whole comparison is
// no longer a valid ordering — the result then depends on the sort's internals, not on the data.
function fusionChrOrder(a, b) {
  if (a === b) return 0;
  const rank = (c) => (c === 'X' ? 23 : (c === 'Y' ? 24 : (Number(c) || 99)));
  // Two different names that rank alike (an unplaced contig, say) fall back to alphabetical order.
  return rank(a) - rank(b) || (a < b ? -1 : 1);
}

// "der(14)" is standard cytogenetic shorthand but not self-explaining. Spell it out once per panel.
function fusionDerivativeTitle(name) {
  const m = /^der\((.+)\)$/.exec(name);
  return m ? `der(${m[1]}) — the rearranged chromosome ${m[1]}` : name;
}

/** One event: header with pills, ideograms + mini circos, then a Report-format panel per derivative. */
// The header tag for a derivative that is only POSSIBLY productive. The default is "possible fusion",
// which is wrong for the kinds where an unchanged protein is deregulated by a new upstream region —
// there is no fusion protein in any of these.
// It is also wrong for a truncated donor, where the summary says there is no fusion protein and the
// most that can form is a shortened copy of the donor. "possible fusion" over that summary is a
// contradiction on the same card. The only kind left on the default is 'chimeric-candidate', where a
// chimeric transcript really is what is expected.
const POSSIBLE_TAG = {
  'deregulation-candidate': 'possible deregulation',
  'promoter-substitution-unnamed-donor': 'possible deregulation',
  'unnamed-donor-intact': 'possible deregulation',
  'enhancer-adoption': 'possible deregulation',
  'truncated-donor': 'possible truncation',
};

function drawFusionEvent(ev, eventIndex = 0) {
  const ctx = fusionCtx();
  // An intra-chromosomal event (an inversion, say) has the same chromosome at both ends, so the
  // chromosome list holds it twice. Draw and colour each chromosome once.
  const chrs = [...new Set(ev.chrs)].sort(fusionChrOrder);
  const derivs = describeDerivatives(ev, ctx.genes, fusionCoreOpts());

  // Per chromosome: the breakpoint positions to mark on its ideogram.
  const perChr = {};
  for (const c of chrs) perChr[c] = { positions: [] };
  // Per breakpoint: what sits there. Keyed by position, not by chromosome — naming per chromosome
  // titles an inversion "MECOM × MECOM" and loses the partner the panel below names correctly.
  // A gene beats a nameless read of the same breakpoint, as before.
  const perBreak = new Map();
  for (const d of derivs) for (const seg of d.reading) {
    if (perChr[seg.chr]) perChr[seg.chr].positions.push(seg.pos);
    const key = `${seg.chr}:${seg.pos}`;
    const prev = perBreak.get(key);
    if (!prev || (!prev.content.gene && seg.content.gene)) perBreak.set(key, seg);
  }
  // Chromosome order, then position: t(3;12) is never titled 12 × 3, and an event whose breakpoints
  // share a chromosome is titled in position order rather than in whatever order the report happened
  // to write its rows. The two junctions of a reciprocal pair sit a few bases apart and name the same
  // partners, so repeats are dropped. Which piece is read first is a separate question, and is what
  // the reading line and the ledger below show.
  const ordered = [...perBreak.values()].sort((a, b) => fusionChrOrder(a.chr, b.chr) || a.pos - b.pos);
  const names = [];
  for (const seg of ordered) { const n = fusionBareName(seg); if (!names.includes(n)) names.push(n); }

  const title = names.join(' × ');
  const productive = derivs.filter((d) => d.product.productive === true);
  const possible = derivs.filter((d) => d.product.productive === 'possible');
  const highGermline = derivs.some((d) => germlineIsStated(d.germline) && d.germline.af >= 0.01);

  const header = el('div', { class: 'row', style: 'justify-content:space-between; align-items:baseline; gap:10px; flex-wrap:wrap' },
    el('div', {}, el('h2', { style: 'margin:0' }, title),
      el('div', { class: 'sub', style: 'font-size:12px; margin:0' },
        `chr${chrs.join(' ↔ chr')}${derivs.map((d) => d.derivative.name).filter(Boolean).length ? ` · ${derivs.map((d) => d.derivative.name).filter(Boolean).join(' + ')}` : ''}`)),
    el('div', { class: 'row', style: 'gap:6px; flex-wrap:wrap' },
      highGermline ? el('span', { class: 'tag v-maybe', title: 'High population germline frequency — a common germline structural variant rather than a somatic driver.' }, 'likely germline') : null,
      // `event-badge` is a hook, not a style: the same words appear in the ledger and in each
      // derivative panel, so a control that searches the whole card cannot tell what this badge says.
      productive.length ? el('span', { class: 'tag v-driver event-badge' }, productive[0].product.kind === 'promoter-substitution' ? 'promoter substitution' : 'productive')
        : (possible.length ? el('span', { class: 'tag v-driver event-badge', style: 'opacity:0.7' }, POSSIBLE_TAG[possible[0].product.kind] || 'possible fusion') : null),
      ev.reciprocal
        ? el('span', { class: 'tag v-unknown', title: `Breakpoint offsets between the two junctions: ${Object.entries(ev.offsets).map(([c, n]) => `chr${c} ${n.toLocaleString()} bp`).join(', ')}` }, `reciprocal · ${ev.confidence}`)
        : el('span', { class: 'tag v-unknown' }, 'single junction')));

  const card = el('div', { class: 'card' }, header);

  // Ideograms on the left, the whole-genome circos with the translocation chord on the right.
  const involved = {};
  Object.assign(involved, chromosomeColours(chrs));   // the figures' own rule, not a second copy
  const junctions = ev.breakends.map((b) => ({ aChr: b.localChr, aPos: b.localPos, bChr: b.mateChr, bPos: b.matePos }));
  const figRow = el('div', { class: 'fus-evfig' });
  const ideos = el('div', { class: 'fus-ideos' });
  ideos.innerHTML = chrs.map((c) => ideogram(ctx, c, perChr[c].positions, { width: 900 })).join('');
  const circ = el('div', { class: 'fus-circos' });
  circ.innerHTML = miniCircos(ctx, junctions, involved);
  figRow.append(ideos, circ);
  card.append(figRow);

  if (ev.reciprocal) card.append(el('div', { class: 'note muted', style: 'font-size:11.5px; margin-top:8px' },
    `The two junctions were matched as a reciprocal pair: they use opposite sides of the same cut on both chromosomes, offset by ${Object.entries(ev.offsets).map(([c, n]) => `${n.toLocaleString()} bp on chr${c}`).join(' and ')}.`
    + (ev.confidence === 'loose' ? ' The larger offset is consistent with a junction-associated deletion, so treat the pairing as probable rather than certain.' : '')));

  // ---- one Report-format panel per derivative, productive first ----
  derivs.forEach((d, di) => {
    const p = d.product;
    const panel = el('div', { class: 'fus-panel', id: `fus-${eventIndex}-${di}` });
    const v = fusionVerdict(p.productive);

    panel.append(el('div', { class: 'fus-head' },
      el('div', {},
        el('h3', {}, d.derivative.name ? fusionDerivativeTitle(d.derivative.name) : `Junction ${d.index + 1}`),
        el('div', { class: 'fus-reading' },
          el('span', { class: 'g' }, fusionPieceLabel(d.reading[0])),
          el('span', { class: 'arrow' }, p.productive === false ? '✕' : '→'),
          el('span', { class: 'g' }, fusionPieceLabel(d.reading[1])))),
      el('div', { class: 'fus-call' },
        el('span', { class: `tag ${v.cls}`, title: FUSION_VERDICT_HELP }, v.word),
        el('span', { class: 'mech' }, fusionKindLabel(p.kind)))));
    panel.append(el('div', { class: 'fus-sub' }, d.derivative.name ? d.derivative.basis : `Not named — ${d.derivative.basis}`));

    // the derivative molecule
    const fig = el('div', { class: 'track-wrap fus-figure' });
    fig.innerHTML = derivativeMolecule(ctx, d, 1080);
    panel.append(fig);

    // what it means — the call's own text, kept with the figure
    panel.append(el('div', { class: 'fus-headline' }, p.headline));
    panel.append(el('p', { class: 'fus-summary' }, p.summary));

    // warnings that need acting on stay visible
    const warns = el('div', { class: 'fus-chips' });
    if (!d.support) warns.append(fusionChip('Read support', 'not present in the pasted row', 'warn'));
    if (d.support && d.support.ambiguous) warns.append(fusionChip('Read support', 'this row carried more than one set of figures and the first was taken', 'warn', 'A normal column beside a tumour column, most likely. Check which sample these counts belong to.'));
    if (germlineIsStated(d.germline) && d.germline.af >= 0.01) warns.append(fusionChip('Population frequency', 'common enough to be an inherited structural variant rather than a tumour change', 'bad'));
    else if (d.germline && d.germline.assumed) warns.append(fusionChip('Frequency column', germlineAssumedText(d.germline), 'warn', GERMLINE_ASSUMED_HELP));
    if (d.derivative.warning) warns.append(fusionChip('Chromosome structure', d.derivative.warning === 'acentric'
      ? 'this junction on its own gives a fragment with no centromere, so it needs its reciprocal partner to survive'
      : 'both pieces carry a centromere, which makes an unstable dicentric chromosome', 'warn', d.derivative.basis));
    if (warns.childElementCount) panel.append(warns);

    // reference loci, boxed apart from the call
    const ref = el('div', { class: 'track-wrap fus-refsep' });
    ref.innerHTML = referenceLoci(ctx, d, 1080);
    panel.append(ref);

    // evidence behind the dark-green bar; any caveat rides inside it as a quiet note line
    panel.append(fusionEvidenceSection(d));
    card.append(panel);
  });
  return card;
}

/** Coordinates, read support and junction bookkeeping — folded behind the green bar. */
function fusionEvidenceSection(d) {
  const rows = el('dl', { class: 'fus-facts' });
  const fact = (term, value, title, cls) => {
    if (value === null || value === undefined || value === '') return;
    rows.append(el('dt', { title: title || null }, term), el('dd', { class: cls || null, title: title || null }, value));
  };
  d.reading.forEach((seg, i) => {
    const g = seg.content;
    const headName = g.gene || (g.igLocus ? `${g.igLocus.symbol} — ${g.igLocus.name}` : `Chromosome ${seg.chr}`);
    rows.append(el('dt', { class: `seg${i ? ' later' : ''}` },
      el('span', { class: 'g' }, headName),
      el('span', { class: 'part' }, g.gene ? `${g.strand} strand · ${partShort(g.part)}` : partShort(g.part))));
    fact('Breakpoint', `chr${seg.chr}:${seg.pos.toLocaleString()}`, null, 'm');
    if (d.bands && d.bands[seg.chr]) fact('Cytoband', d.bands[seg.chr].label, 'Taken from the pasted row, not computed here.');
    if (g.gene) {
      if (g.tx) fact('Transcript', g.tx, 'The MANE Select transcript this gene’s structure came from.', 'm');
      if (seg.exon) {
        fact('Structure', `${seg.exon.exons} exons`);
        fact('Breakpoint falls', seg.exon.label);
        fact('Coding kept', `${seg.exon.codingKept.toLocaleString()} of ${seg.exon.codingTotal.toLocaleString()} bases${seg.exon.codingIntact ? ' — all of it' : ''}`);
      } else if (g.keptFraction !== null && g.keptFraction !== undefined) {
        fact('Retained (span)', `${Math.round(g.keptFraction * 100)}% of ${g.span.toLocaleString()} bp`,
          'How much of the gene’s genomic span survives on the kept side. A span measurement — exon coordinates were not available.');
      }
    } else if (g.nearestGene && !g.igLocus) {
      fact('Nearest gene', `${g.nearestGene.gene}, ${(g.nearestGene.distance / 1000).toFixed(1)} kb away`,
        'Reported so a true gene desert can be told apart from a breakpoint that fell just outside the naming window.');
    }
  });

  const supportNote = 'Reproduced exactly as the caller wrote them, with no fraction derived here. '
    + 'Confirm whether your caller reports supporting/total or reference/supporting — the same pair of numbers means opposite things under the two conventions.';
  rows.append(el('dt', { class: 'seg later' }, el('span', { class: 'g' }, 'The junction itself')));
  if (d.support && d.support.pr) fact('Read pairs', `${d.support.pr.alt} of ${d.support.pr.total}`, `Read pairs whose two ends land on opposite sides of the junction. ${supportNote}`, 'm');
  if (d.support && d.support.sr) fact('Split reads', `${d.support.sr.alt} of ${d.support.sr.total}`, `Single reads that cross the breakpoint itself, part aligning to each side. ${supportNote}`, 'm');
  if (germlineIsStated(d.germline)) {
    fact('Population frequency', `${(d.germline.af * 100).toPrecision(3)}% (population database)`,
      'How often this junction is seen in people generally. A high figure points to an inherited variant rather than a tumour change.',
      d.germline.af >= 0.01 ? 'bad' : null);
  } else if (d.germline) {
    // Untagged: report the number, say plainly that its meaning is not stated, and read nothing into it.
    // No colour class — styles.css has no dd.warn rule, and the wording carries the caveat.
    fact('Frequency column', germlineAssumedText(d.germline), GERMLINE_ASSUMED_HELP, null);
  }
  if (d.product.frame) {
    const f = d.product.frame;
    fact('Reading frame', f.inFrame ? 'in frame' : `out of frame by ${f.shift} ${f.shift === 1 ? 'base' : 'bases'}`);
    fact('Fused coding', `${f.donorCodingBases.toLocaleString()} + ${f.acceptorCodingBases.toLocaleString()} = ${f.fusedCodingBases.toLocaleString()} bases`);
  }
  if (d.insertion) fact('Inserted bases', d.insertion, 'Novel or microhomologous sequence sitting at the join.', 'm');
  fact('Orientation', d.reverseComplement ? 'one piece is reverse-complemented' : 'both pieces in the same orientation',
    'Bookkeeping about how the two pieces are joined in the genome. The reading above is already in transcript order regardless.');
  if (d.mateRecords.length > 1) fact('Source rows', `${d.mateRecords.length}, merged`, `The same junction written from both ends:\n${d.mateRecords.join('\n')}`);
  fact('Breakend notation', d.notation, 'The row exactly as pasted.', 'm');

  const fold = el('details', { class: 'fus-fold evid' }, el('summary', {}, 'Coordinates and read support'), rows);
  // A genuine qualifier (the frame call assumes canonical splicing; an overlapping-gene alternative)
  // rides here as a quiet note, not as its own block.
  if (d.product.caveats && d.product.caveats.length) {
    fold.append(el('ul', { class: 'fus-caveats' }, ...d.product.caveats.map((c) => el('li', {}, c))));
  }
  return fold;
}

// =====================================================================
//  TABS + INIT
// =====================================================================
// ---- Pathway diagrams (drawn, not written) ----

const SVG_NS = 'http://www.w3.org/2000/svg';
function svg(tag, attrs = {}, ...kids) {
  const n = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) n.setAttribute(k, String(v));
  for (const k of kids) if (k !== null && k !== undefined) n.append(k);
  return n;
}

// One colour per kind of thing, so the shape of a pathway reads at a glance.
const NODE_STYLE = {
  ligand: { fill: '#1d3b4a', stroke: '#3d7d96', text: '#bfe6f2' },
  receptor: { fill: '#173f45', stroke: '#3aa4a0', text: '#b6f0e6' },
  kinase: { fill: '#2b2450', stroke: '#7d6bd0', text: '#d8cffb' },
  adaptor: { fill: '#1f3357', stroke: '#5b83c9', text: '#cbdcff' },
  tf: { fill: '#4a3517', stroke: '#c08a35', text: '#ffe0ad' },
  complex: { fill: '#2f2f38', stroke: '#7a7a8c', text: '#dcdce6' },
  process: { fill: '#123024', stroke: '#3f9469', text: '#b8f0d0' },
  other: { fill: '#26262e', stroke: '#6a6a78', text: '#d0d0da' },
};

const NODE_W = 118;
const NODE_H = 34;
const COL_W = 158;
const ROW_H = 76;
const PAD = 16;

/**
 * Draw one pathway as an SVG. The assessed gene's node is highlighted so the reader can see
 * immediately where their variant sits in the cascade.
 */
function buildPathwayDiagram(pathway, highlightGene) {
  const nodes = pathway.nodes || [];
  const byId = new Map(nodes.map((n) => [n.id, n]));
  const maxCol = Math.max(0, ...nodes.map((n) => n.col));
  const maxRow = Math.max(0, ...nodes.map((n) => n.row));
  const width = PAD * 2 + maxCol * COL_W + NODE_W;
  const height = PAD * 2 + maxRow * ROW_H + NODE_H;

  const xOf = (n) => PAD + n.col * COL_W;
  const yOf = (n) => PAD + n.row * ROW_H;
  const cx = (n) => xOf(n) + NODE_W / 2;
  const cy = (n) => yOf(n) + NODE_H / 2;

  const root = svg('svg', {
    viewBox: `0 0 ${width} ${height}`, width: '100%',
    style: `max-width:${width}px; height:auto; display:block`,
    role: 'img', 'aria-label': `${pathway.name} pathway diagram`,
  });

  const defs = svg('defs');
  defs.append(svg('marker', { id: 'arrow-act', viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 6, markerHeight: 6, orient: 'auto-start-reverse' },
    svg('path', { d: 'M 0 1 L 10 5 L 0 9 z', fill: '#8fb6d6' })));
  root.append(defs);

  // Edges first, so nodes sit on top of them.
  for (const e of pathway.edges || []) {
    const a = byId.get(e.from); const b = byId.get(e.to);
    if (!a || !b) continue;                       // never draw an edge to a node that is absent
    const x1 = cx(a); const y1 = cy(a); const x2 = cx(b); const y2 = cy(b);
    const dx = x2 - x1; const dy = y2 - y1;
    const len = Math.hypot(dx, dy) || 1;
    // Pull both ends back to the box edge so the line does not run under the label.
    const shrinkA = (Math.abs(dy) > Math.abs(dx) ? NODE_H / 2 : NODE_W / 2) + 2;
    const shrinkB = (Math.abs(dy) > Math.abs(dx) ? NODE_H / 2 : NODE_W / 2) + 8;
    const sx = x1 + (dx / len) * shrinkA; const sy = y1 + (dy / len) * shrinkA;
    const ex = x2 - (dx / len) * shrinkB; const ey = y2 - (dy / len) * shrinkB;
    const inhibit = e.kind === 'inhibits';
    root.append(svg('line', {
      x1: sx, y1: sy, x2: ex, y2: ey,
      stroke: inhibit ? '#d8737f' : '#8fb6d6',
      'stroke-width': 1.6,
      'stroke-dasharray': inhibit ? '4 3' : null,
      'marker-end': inhibit ? null : 'url(#arrow-act)',
    }));
    if (inhibit) {
      // A blunt crossbar is the convention for inhibition; an arrowhead would mislead.
      const nx = -(dy / len); const ny = dx / len;
      root.append(svg('line', {
        x1: ex - nx * 7, y1: ey - ny * 7, x2: ex + nx * 7, y2: ey + ny * 7,
        stroke: '#d8737f', 'stroke-width': 2.4,
      }));
    }
  }

  for (const n of nodes) {
    const st = NODE_STYLE[n.kind] || NODE_STYLE.other;
    const isHit = highlightGene && n.gene === highlightGene;
    const g = svg('g');
    g.append(svg('rect', {
      x: xOf(n), y: yOf(n), width: NODE_W, height: NODE_H, rx: n.kind === 'process' ? 16 : 6,
      fill: st.fill, stroke: isHit ? '#ffcf5c' : st.stroke, 'stroke-width': isHit ? 2.6 : 1.3,
    }));
    const t = svg('text', {
      x: cx(n), y: cy(n) + 4, 'text-anchor': 'middle',
      'font-size': 12.5, 'font-weight': isHit ? 700 : 500,
      fill: isHit ? '#ffcf5c' : st.text,
      'font-family': 'system-ui, -apple-system, Segoe UI, sans-serif',
    });
    t.textContent = n.label;
    g.append(t);
    root.append(g);
  }
  return root;
}

// ---- Gene context: what the protein does, and how it is reported to drive cancer ----

const MECHANISM_STYLE = {
  'gain-of-function': { label: 'Gain of function', bg: '#3a2130', fg: '#ffb3c1', border: '#a8546b' },
  'loss-of-function': { label: 'Loss of function', bg: '#152f3d', fg: '#9fd8ef', border: '#4a89a8' },
  'context-dependent': { label: 'Context-dependent', bg: '#3d3417', fg: '#f0d79a', border: '#a8873a' },
};

let contextGene = null;
let pathwayData = null;
let lastVariant = null;      // {hgvs, gene, coord} of the most recently resolved variant
let lastGnomad = null;       // cache of the last gnomAD result, keyed by hgvs

// The gnomAD dataset the app queries (gnomad_r4 = release v4.1.0).
const GNOMAD_VERSION = 'v4.1.0';

// ---- Splice analysis (SpliceAI deltas + a drawn exon model) ----

let lastSplice = null;   // cache of the last SpliceAI result, keyed by hgvs

const SPLICE_COLOUR = (ds) => (ds >= 0.8 ? '#e05561' : ds >= 0.5 ? '#e0a54e' : ds >= 0.2 ? '#d9b44e' : '#5b7089');

// SVG exon model: exons to scale, the variant marked, and each predicted splice site (variant
// position + SpliceAI delta position) marked and coloured by its delta score.
function buildSpliceDiagram(sp) {
  const starts = sp.exonStarts || []; const ends = sp.exonEnds || [];
  if (!starts.length) return el('div', { class: 'note muted' }, 'No exon model returned.');
  const sitePts = (sp.sites || []).filter((s) => s.ds >= 0.2).map((s) => s.site);
  const lo = Math.min(...starts, sp.pos, ...sitePts);
  const hi = Math.max(...ends, sp.pos, ...sitePts);
  const W = 860; const H = 96; const pad = 34; const midY = 58; const span = (hi - lo) || 1;
  const x = (g) => pad + ((g - lo) / span) * (W - 2 * pad);

  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px; height:auto`, role: 'img', 'aria-label': 'exon model with predicted splice sites' });
  // Intron baseline.
  root.append(svg('line', { x1: x(lo), y1: midY, x2: x(hi), y2: midY, stroke: '#33506f', 'stroke-width': 1.5 }));
  // Exons.
  for (let i = 0; i < starts.length; i++) {
    const x1 = Math.min(x(starts[i]), x(ends[i])); const x2 = Math.max(x(starts[i]), x(ends[i]));
    root.append(svg('rect', { x: x1, y: midY - 8, width: Math.max(1.5, x2 - x1), height: 16, rx: 2, fill: '#2b3f5c', stroke: '#4a6d96' }));
  }
  // 5' / 3' ends by strand (genomic left-to-right; minus strand reads right-to-left).
  const leftLbl = sp.strand === '-' ? "3'" : "5'"; const rightLbl = sp.strand === '-' ? "5'" : "3'";
  for (const [xx, lbl] of [[x(lo), leftLbl], [x(hi), rightLbl]]) {
    const t = svg('text', { x: xx, y: midY - 16, 'text-anchor': 'middle', 'font-size': 11, fill: '#8fa4bf', 'font-family': 'system-ui, sans-serif' }); t.textContent = lbl; root.append(t);
  }
  // Variant marker.
  root.append(svg('line', { x1: x(sp.pos), y1: midY - 20, x2: x(sp.pos), y2: midY + 20, stroke: '#ffcf5c', 'stroke-width': 2 }));
  const vt = svg('text', { x: x(sp.pos), y: midY + 34, 'text-anchor': 'middle', 'font-size': 11, 'font-weight': 700, fill: '#ffcf5c', 'font-family': 'system-ui, sans-serif' }); vt.textContent = 'variant'; root.append(vt);
  // Predicted splice sites.
  for (const s of (sp.sites || []).filter((e) => e.ds >= 0.2)) {
    const col = SPLICE_COLOUR(s.ds);
    root.append(svg('circle', { cx: x(s.site), cy: midY, r: 5, fill: col, stroke: '#0f1420', 'stroke-width': 1 }));
    const lt = svg('text', { x: x(s.site), y: midY - 22, 'text-anchor': 'middle', 'font-size': 10.5, fill: col, 'font-family': 'system-ui, sans-serif' }); lt.textContent = `${s.key} ${s.ds}`; root.append(lt);
  }
  return root;
}

async function renderSplice() {
  const view = $('#view-splice'); view.innerHTML = '';
  // The variant this render belongs to, captured once. Everything below reads `variant`, never the
  // module-level `lastVariant`, so an answer cannot be filed under a variant assessed since it was
  // asked for; `seq` then stops a superseded answer being drawn or cached at all.
  const variant = lastVariant;
  const seq = currentAssessment();
  if (!variant || !variant.hgvs) {
    view.append(el('div', { class: 'card' }, el('h2', {}, 'Splice'), el('div', { class: 'note' }, 'Assess a variant to see its splice analysis.')));
    return;
  }
  // The badge used to read "SpliceAI v1.3" as a fixed string. The scores come from a service outside
  // this app whose hosts and model can change, and the answer carries no version, so the name is
  // stated and the version only if one ever arrives with the scores.
  const versionBadge = el('span', { style: 'padding:3px 10px; border-radius:11px; font-size:12px; background:#2b2450; color:#d8cffb; border:1px solid #7d6bd0' }, 'SpliceAI');
  view.append(el('div', { class: 'card' },
    el('div', { class: 'row', style: 'align-items:center; gap:10px; flex-wrap:wrap' },
      el('h2', { style: 'margin:0' }, `Splice — ${variant.gene || ''} ${variant.appProtein || ''}`.trim()),
      versionBadge),
    el('div', { class: 'sub', style: 'margin-top:4px' }, IS_WEB
        ? 'Donor and acceptor delta scores. Sends the variant coordinate to the SpliceAI service — not offline.'
        : 'Donor/acceptor delta scores and experimentally validated splice evidence. Sends the variant coordinate to the SpliceAI service and SpliceVarDB — not offline.')));

  const holder = el('div', {});
  const svHolder = el('div', {});
  view.append(holder, svHolder);
  holder.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'Loading SpliceAI…')));

  // Resolve the plus-strand coordinate once (from the assessment, else via Ensembl) and share it.
  let coord = variant.coord || null;
  if (!coord) {
    try { const c = await window.api.resolveCoord(variant.hgvs); if (c && c.ok) coord = { chr: c.chr, pos: c.pos, ref: c.ref, alt: c.alt }; } catch { /* offline */ }
    // The coordinate belongs to the variant it was resolved for, never to whichever variant is on
    // screen when it arrives — the Assessment tab reuses this field for its ClinVar and gnomAD calls.
    if (coord) variant.coord = coord;
    if (superseded(seq)) return;
  }
  // SpliceVarDB loads in parallel with SpliceAI (independent network call).
  // SpliceVarDB sends no cross-origin header, so a browser can never read its reply. On the website the
  // card is left out entirely rather than shown as a standing "unavailable" notice, which tells the
  // reader nothing they can act on. The desktop application still has it.
  if (!IS_WEB) renderSpliceVarDb(svHolder, coord, variant, seq);

  let res = null;
  // Only a successful answer is worth keeping. A failure used to be cached under the variant too, so
  // once the host had been unreachable once, every later visit to this tab was served that failure
  // and the call was never made again — while the card invited the reader to expect a retry to work.
  if (lastSplice && lastSplice.hgvs === variant.hgvs && lastSplice.res && lastSplice.res.ok) res = lastSplice.res;
  else {
    try { res = await window.api.spliceai(variant.hgvs, variant.appTranscript, coord); }
    catch (e) { res = { ok: false, error: e.message }; }
    if (superseded(seq)) return;
    if (res && res.ok) lastSplice = { hgvs: variant.hgvs, res };
    else if (lastSplice && lastSplice.hgvs === variant.hgvs) lastSplice = null;
  }
  holder.innerHTML = '';
  // The version the scores actually came back with, if the service ever reports one.
  const reportedVersion = res && (res.version || res.modelVersion);
  if (reportedVersion) versionBadge.textContent = `SpliceAI ${reportedVersion}`;

  if (!res || !res.ok) {
    holder.append(el('div', { class: 'card' }, el('div', { class: 'note warn' }, `SpliceAI unavailable: ${(res && res.error) || 'no response'}.`),
      el('div', { class: 'note muted', style: 'margin-top:6px; font-size:11.5px' }, 'This needs the internet; the SpliceAI host can also be intermittently unreachable.'),
      el('div', { class: 'row', style: 'margin-top:8px' },
        el('button', { class: 'btn small', onclick: () => renderSplice() }, 'Try again'))));
  } else if (res.notScored) {
    holder.append(el('div', { class: 'card' }, el('div', { class: 'note' }, 'SpliceAI returned no score for this transcript.')));
  } else {
    const interpColour = { high: '#e05561', moderate: '#e0a54e', low: '#d9b44e', minimal: '#46b877' }[res.interpretation];
    const tiles = res.events.map((ev) => el('div', { style: 'flex:1 1 150px; min-width:150px; background:#152232; border:1px solid #26405c; border-radius:9px; padding:10px 12px' },
      el('div', { class: 'sub', style: 'font-size:11px; text-transform:uppercase' }, ev.label),
      el('div', { class: 'mono', style: 'font-size:19px; font-weight:700; color:#dbe9fb; margin-top:2px' }, ev.ds == null ? '—' : ev.ds.toFixed(2)),
      el('div', { style: 'margin-top:6px; height:5px; background:#0e1826; border-radius:3px; overflow:hidden' },
        el('div', { style: `width:${Math.round((ev.ds || 0) * 100)}%; height:100%; background:${SPLICE_COLOUR(ev.ds || 0)}` })),
      ev.dp != null ? el('div', { class: 'note muted', style: 'font-size:11px; margin-top:4px' }, `${ev.dp > 0 ? '+' : ''}${ev.dp} nt`) : null));

    holder.append(el('div', { class: 'card' },
      el('div', { class: 'row', style: 'justify-content:space-between; align-items:baseline; gap:12px' },
        el('h2', { style: 'margin:0' }, 'SpliceAI delta scores'),
        el('span', { style: `padding:3px 10px; border-radius:11px; font-size:12px; font-weight:700; background:${interpColour}22; color:${interpColour}; border:1px solid ${interpColour}` }, `${res.interpretation} impact (max ${res.maxDs.toFixed(2)})`)),
      el('div', { class: 'row', style: 'gap:10px; flex-wrap:wrap; margin-top:12px' }, ...tiles),
      el('div', { class: 'note muted', style: 'margin-top:10px; font-size:11.5px' }, `Transcript ${res.transcript} · Δ ≥ 0.2 low · ≥ 0.5 moderate · ≥ 0.8 high (Broad SpliceAI, CC BY-NC).`)));

    holder.append(el('div', { class: 'card' }, el('h2', { style: 'font-size:15px' }, 'Exon model'),
      el('div', { style: 'overflow-x:auto; margin-top:8px' }, buildSpliceDiagram(res)),
      el('div', { class: 'note muted', style: 'margin-top:6px; font-size:11.5px' }, 'Exons to scale; the amber line is the variant; coloured dots are SpliceAI-predicted gain/loss sites (Δ ≥ 0.2).')));
  }
}

// `variant` and `seq` are the same pair renderSplice captured, so this panel answers for the variant
// it was opened for and drops its answer if the user has moved on.
async function renderSpliceVarDb(host, coord, variant, seq) {
  const body = el('div', {}, el('div', { class: 'note' }, 'Loading…'));
  host.append(el('div', { class: 'card' }, el('h2', { style: 'font-size:15px' }, 'Experimental evidence (SpliceVarDB)'), body));
  if (!coord) {
    body.innerHTML = ''; body.append(el('div', { class: 'note muted' }, 'Could not resolve the genomic coordinate for this variant.'));
    return;
  }
  let res = null;
  try { res = await window.api.splicevardb(coord, variant.gene); } catch (e) { res = { ok: false, error: e.message }; }
  if (superseded(seq)) return;
  body.innerHTML = '';
  if (!res || !res.ok) { body.append(el('div', { class: 'note warn' }, `SpliceVarDB unavailable: ${(res && res.error) || 'no response'}.`)); return; }

  if (res.matches && res.matches.length) {
    const rows = res.matches.map((m) => el('tr', {},
      el('td', {}, el('span', { class: `tag ${/splice-altering/i.test(m.classification) ? 'v-path' : /normal/i.test(m.classification) ? 'v-benign' : 'v-vus'}` }, m.classification || '—')),
      el('td', {}, m.method || ''),
      el('td', {}, m.location || ''),
      el('td', { class: 'mono', style: 'font-size:12px' }, m.hgvs || ''),
      el('td', {}, el('a', { onclick: () => window.api.openExternal(m.url) }, 'view'))));
    body.append(
      el('div', { class: 'note', style: 'margin-bottom:8px' }, el('strong', {}, 'Experimentally tested for splicing. ')),
      el('div', { style: 'overflow-x:auto' }, el('table', { class: 'ref' },
        el('thead', {}, el('tr', {}, ['Result', 'Method', 'Location', 'HGVS', ''].map((h) => el('th', {}, h)))),
        el('tbody', {}, ...rows))));
  } else if (res.geneCovered) {
    body.append(el('div', { class: 'note' }, 'This variant has not been experimentally tested for splicing in SpliceVarDB.'),
      el('div', { class: 'note muted', style: 'margin-top:4px; font-size:11.5px' }, 'Absence here is not evidence the variant is splice-neutral — it simply has no experimental record.'));
  } else {
    body.append(el('div', { class: 'note muted' }, `${variant.gene || 'This gene'} is not covered by SpliceVarDB.`));
  }
  body.append(el('div', { class: 'note muted', style: 'margin-top:8px; font-size:11px' }, 'SpliceVarDB (Sullivan et al. 2024) — only the gene symbol was sent.'));
}

// ---- gnomAD population frequency ----

function pctFmt(af) {
  if (af == null) return '—';
  if (af === 0) return '0';
  if (af < 1e-4) return `${(af * 100).toExponential(1)} %`;
  return `${(af * 100).toPrecision(3)} %`;
}

// "1 in N", where N is exactly the allele number divided by the allele count (AN ÷ AC).
function oneInN(ac, an) {
  if (!ac || !an) return '—';
  return `1 : ${(an / ac).toLocaleString(undefined, { maximumFractionDigits: 1 })}`;
}

function afReading(overallAf) {
  // A plain-language read of what the frequency means for a somatic haem context.
  if (overallAf == null) return '';
  if (overallAf === 0) return 'Not observed';
  if (overallAf < 0.00005) return 'Ultra-rare in the population';
  if (overallAf < 0.001) return 'Rare in the population';
  if (overallAf < 0.01) return 'Uncommon in the population';
  return 'Common in the population — treat with caution as a driver';
}

async function renderGnomad() {
  const view = $('#view-gnomad');
  view.innerHTML = '';

  const input = el('input', {
    type: 'text', id: 'gnomad-input', spellcheck: 'false',
    placeholder: 'e.g. NM_004333.6:c.1799T>A',
    value: (lastVariant && lastVariant.hgvs) || '',
    style: 'flex:1; min-width:280px; padding:7px 10px; border-radius:6px',
  });
  const btn = el('button', { class: 'btn small' }, 'Look up in gnomAD');
  const result = el('div', {});

  const controls = el('div', { class: 'card' },
    el('div', { class: 'row', style: 'align-items:center; gap:10px; flex-wrap:wrap' },
      el('h2', { style: 'margin:0' }, 'gnomAD population frequency'),
      el('span', {
        style: 'display:inline-block; padding:3px 10px; border-radius:11px; font-size:13px; font-weight:700; background:#173f45; color:#b6f0e6; border:1px solid #3aa4a0',
      }, `gnomAD ${GNOMAD_VERSION}`)),
    el('div', { class: 'sub', style: 'margin:6px 0 10px' },
      'Sends only the public genomic coordinate (via Ensembl to place it, then gnomAD) — never patient data.'),
    el('div', { class: 'row', style: 'gap:10px; flex-wrap:wrap' }, input, btn,
      el('button', { class: 'btn small secondary', onclick: () => { input.value = ''; result.innerHTML = ''; } }, 'Clear')));
  view.append(controls, result);

  const lookup = async () => {
    const hgvs = input.value.trim();
    if (!hgvs) return;
    result.innerHTML = '';
    result.append(el('div', { class: 'card' }, el('div', { class: 'note' }, `Looking up ${hgvs} …`)));
    let res = null;
    try { res = await window.api.gnomad(hgvs, (lastVariant && lastVariant.hgvs === hgvs && lastVariant.coord) || null); }
    catch (e) { res = { ok: false, error: e.message }; }
    lastGnomad = { hgvs, res };
    renderGnomadResult(result, hgvs, res);
  };
  btn.addEventListener('click', lookup);
  input.addEventListener('keydown', (e) => { if (e.key === 'Enter') lookup(); });

  // If we arrived from a freshly assessed variant, look it up automatically.
  if (lastGnomad && lastGnomad.hgvs === input.value.trim()) renderGnomadResult(result, lastGnomad.hgvs, lastGnomad.res);
  else if (input.value.trim()) lookup();
}

// A clickable link to the variant's gnomAD page; opens in a new tab (system browser on desktop).
function gnomadLink(res) {
  if (!res || !res.variantId) return null;
  const url = `https://gnomad.broadinstitute.org/variant/${res.variantId}?dataset=${res.dataset || 'gnomad_r4'}`;
  const a = el('a', {
    style: 'cursor:pointer; color:#6b9fe0; text-decoration:underline; font-size:12.5px; white-space:nowrap',
    title: url,
  }, `Open ${res.variantId} in the gnomAD browser ↗`);
  a.addEventListener('click', () => window.api.openExternal(url));
  return a;
}

function renderGnomadResult(host, hgvs, res) {
  host.innerHTML = '';
  if (!res || !res.ok) {
    host.append(el('div', { class: 'card' },
      el('div', { class: 'note warn' }, `gnomAD lookup unavailable: ${(res && res.error) || 'no response'}.`),
      el('div', { class: 'note muted', style: 'margin-top:6px; font-size:11.5px' }, 'This feature needs the internet; a locked-down network may block it.')));
    return;
  }
  if (res.notInGnomad) {
    host.append(el('div', { class: 'card' },
      el('h2', { style: 'margin:0 0 6px' }, 'Not observed in gnomAD'),
      el('div', { class: 'note' }, `${res.variantId || hgvs} is absent from gnomAD v4 — consistent with a somatic or very rare variant.`),
      gnomadLink(res) ? el('div', { style: 'margin-top:8px' }, gnomadLink(res)) : null));
    return;
  }

  const tile = (label, big, sub) => el('div', {
    style: 'flex:1; min-width:150px; background:#152232; border:1px solid #26405c; border-radius:10px; padding:12px 14px',
  }, el('div', { class: 'sub', style: 'font-size:11.5px; text-transform:uppercase; letter-spacing:.04em' }, label),
    el('div', { style: 'font-size:22px; font-weight:700; margin-top:2px; color:#dbe9fb' }, big),
    sub ? el('div', { class: 'note muted', style: 'margin-top:2px; font-size:11.5px' }, sub) : null);

  const o = res.overall || {};
  // Explicit exomes + genomes split, so it is clear the headline is the combined figure.
  const seqLine = el('div', { class: 'note muted', style: 'margin-top:10px; font-size:12px' },
    `Combined exomes (WES) + genomes (WGS). ` +
    `Exomes: ${res.exome ? `${pctFmt(res.exome.af)} (${res.exome.ac}/${res.exome.an})` : 'not covered'} · ` +
    `Genomes: ${res.genome ? `${pctFmt(res.genome.af)} (${res.genome.ac}/${res.genome.an})` : 'not covered'}.`);

  const headline = el('div', { class: 'card' },
    el('div', { class: 'row', style: 'justify-content:space-between; align-items:baseline; gap:12px' },
      el('h2', { style: 'margin:0' }, res.variantId || hgvs),
      res.rsids && res.rsids.length ? el('div', { class: 'sub' }, res.rsids.join(', ')) : null),
    el('div', { class: 'row', style: 'gap:12px; flex-wrap:wrap; margin-top:12px' },
      tile('Overall frequency (WES + WGS)', pctFmt(o.af), afReading(o.af)),
      tile('Frequency as a ratio (AN ÷ AC)', oneInN(o.ac, o.an), `${o.ac == null ? '—' : o.ac.toLocaleString()} in ${o.an == null ? '—' : o.an.toLocaleString()} alleles`),
      res.popmax ? tile('Highest population', pctFmt(res.popmax.af), res.popmax.label) : null),
    seqLine);
  host.append(headline);

  // Per-population figure: a horizontal bar per genetic-ancestry group.
  if (res.populations && res.populations.length) {
    const maxAf = Math.max(...res.populations.map((p) => p.af), 1e-9);
    const rows = res.populations.map((p) => el('div', { class: 'row', style: 'align-items:center; gap:10px; margin:4px 0' },
      el('div', { style: 'width:180px; font-size:12.5px; color:#c7d6e8' }, p.label),
      el('div', { style: 'flex:1; background:#0e1826; border-radius:5px; height:16px; overflow:hidden' },
        el('div', { style: `width:${Math.max(2, (p.af / maxAf) * 100)}%; height:100%; background:linear-gradient(90deg,#3a7bd0,#6b9fe0)` })),
      el('div', { class: 'mono', style: 'width:130px; text-align:right; font-size:12px; color:#9fb6d0' }, `${pctFmt(p.af)} · ${p.ac}/${p.an}`)));
    host.append(el('div', { class: 'card' },
      el('div', { class: 'ib-h', style: 'margin-bottom:8px' }, 'By genetic-ancestry group'),
      ...rows));
  }

  const footKids = [];
  if (res.flags && res.flags.length) footKids.push(el('span', { class: 'note muted', style: 'font-size:11.5px' }, `Flags: ${res.flags.join(', ')}. `));
  const gl = gnomadLink(res);
  if (gl) footKids.push(gl);
  host.append(el('div', { class: 'row', style: 'align-items:center; gap:10px; margin:4px 4px 14px; flex-wrap:wrap' }, ...footKids,
    el('span', { class: 'note muted', style: 'font-size:11.5px' }, res.cached ? 'Cached on this machine.' : 'gnomAD v4.')));
}

// Break UniProt function prose into individual statements, one per bullet. Splits on sentence
// boundaries only, and rejoins fragments left by abbreviations (e.g. "TLR9." mid-clause) so a bullet
// is never a stray half-sentence.
function splitSentences(text) {
  return String(text || '')
    .split(/(?<=[.;])\s+(?=[A-Z0-9])/)
    .map((s) => s.trim())
    .filter((s) => s.length > 3);
}

// Only records which gene to draw. Rendering is left to switchTab, so that clicking through
// from the assess card does not start two overlapping async renders into the same view.
function setGeneContext(gene) {
  contextGene = gene || null;
}

async function renderGeneContext() {
  const view = $('#view-context');
  view.innerHTML = '';

  if (!pathwayData) {
    try { pathwayData = await window.api.getPathways(); } catch { pathwayData = { pathways: [], geneIndex: [] }; }
  }
  if (!contextGene) {
    view.append(el('div', { class: 'card' }, el('h2', {}, 'Pathway'),
      el('div', { class: 'note' }, 'Assess a variant to see its pathway.')));
    return;
  }

  // Two layers: the curated haem-onc pathway (drawn diagram + mechanism), where one exists, and
  // the whole-genome context (what the protein does + its canonical Reactome pathways), which
  // exists for almost every gene. Whichever is available is shown, so the tab is informative for
  // any gene, not only the ~68 on the curated maps.
  const entry = (pathwayData.geneIndex || []).find((g) => g.gene === contextGene);
  let ctx = null;
  try { ctx = await window.api.getGeneContext(contextGene); } catch { /* all-gene table absent */ }
  const info = ctx && ctx.entry;

  if (!entry && !info) {
    view.append(el('div', { class: 'card' }, el('h2', {}, contextGene),
      el('div', { class: 'note' }, 'No pathway or functional information is held for this gene.')));
    return;
  }

  // Header: gene name, and — only for a curated haem gene — the mechanism chip and clause.
  const headKids = [el('div', { class: 'row', style: 'align-items:center; gap:12px; flex-wrap:wrap' },
    el('h2', { style: 'margin:0' }, contextGene),
    entry ? (() => {
      const mech = MECHANISM_STYLE[entry.mechanism] || MECHANISM_STYLE['context-dependent'];
      return el('span', {
        style: `display:inline-block; padding:3px 10px; border-radius:11px; font-size:12px; font-weight:600;
                background:${mech.bg}; color:${mech.fg}; border:1px solid ${mech.border}`,
      }, mech.label);
    })() : null,
    info && info.uniprot ? el('span', { class: 'sub' }, `UniProt ${info.uniprot}`) : null)];
  if (entry) headKids.push(el('div', { class: 'note', style: 'margin-top:6px' }, entry.note));
  if (entry && entry.malignancies) headKids.push(el('div', { class: 'sub', style: 'margin-top:2px' }, entry.malignancies));
  else if (info && info.proteinName) headKids.push(el('div', { class: 'sub', style: 'margin-top:4px' }, info.proteinName));
  view.append(el('div', { class: 'card' }, ...headKids));

  // What the protein does — the UniProt function, broken into one scannable bullet per statement
  // rather than a wall of prose. The first sentence leads; the rest are bullets.
  if (info && info.function) {
    const sentences = splitSentences(info.function);
    const card = el('div', { class: 'card' }, el('div', { class: 'ib-h', style: 'margin-bottom:8px' }, 'What it does'));
    if (sentences.length) {
      card.append(el('div', { class: 'note', style: 'line-height:1.5; margin-bottom:8px' }, sentences[0]));
      if (sentences.length > 1) {
        const ul = el('ul', { style: 'margin:0; padding-left:18px; line-height:1.5' });
        for (const s of sentences.slice(1)) ul.append(el('li', { class: 'note', style: 'margin:3px 0' }, s));
        card.append(ul);
      }
    }
    view.append(card);
  }

  // Curated haem-onc diagrams, where this gene appears on one.
  if (entry) {
    const byId = new Map((pathwayData.pathways || []).map((p) => [p.id, p]));
    for (const pid of entry.pathwayIds || []) {
      const p = byId.get(pid);
      if (!p) continue;
      const card = el('div', { class: 'card' });
      card.append(el('div', { class: 'row', style: 'justify-content:space-between; align-items:baseline; gap:12px' },
        el('h2', { style: 'font-size:15px; margin:0' }, p.name),
        p.context ? el('div', { class: 'sub' }, p.context) : null));
      card.append(el('div', { style: 'overflow-x:auto; margin-top:10px' }, buildPathwayDiagram(p, contextGene)));
      view.append(card);
    }
    view.append(el('div', { class: 'note muted', style: 'margin:2px 4px 8px; font-size:11.5px' },
      'Arrow = activates · bar = inhibits · highlighted box = this gene. Curated haematology-oncology pathway.'));
  }

  // If this gene is on no drawn pathway figure, say so plainly — do NOT fall back to a list of
  // pathway names, which reads as clutter rather than a figure.
  if (!entry) {
    view.append(el('div', { class: 'note muted', style: 'margin:2px 4px 14px; font-size:11.5px' },
      'No pathway diagram for this gene yet. Function above is from UniProt (reviewed).'));
  } else {
    view.append(el('div', { class: 'note muted', style: 'margin:2px 4px 14px; font-size:11.5px' },
      'Function: UniProt (reviewed). Suggestive context, not a clinical determination.'));
  }
}

function switchTab(view) {
  const tabs = $('#tabs');
  for (const b of tabs.querySelectorAll('button')) b.classList.toggle('active', b.dataset.view === view);
  for (const v of document.querySelectorAll('.view')) v.classList.toggle('active', v.id === `view-${view}`);
  if (view === 'about') buildAbout();
  if (view === 'context') renderGeneContext();
  if (view === 'gnomad') renderGnomad();
  if (view === 'splice') renderSplice();
  if (view === 'genebe') renderGeneBe();       // original GeneBe assessment (desktop only)
  if (view === 'assessment') openAssessment(); // new redesigned Assessment (both platforms)
  if (view === 'assess' && viewer) setTimeout(() => { try { viewer.resize(); viewer.render(); } catch { /* */ } }, 50);
}

function setupTabs() {
  $('#tabs').addEventListener('click', (e) => {
    const btn = e.target.closest('button'); if (!btn) return;
    switchTab(btn.dataset.view);
  });
}

async function init() {
  // The desktop app hides the native title bar and reserves space for colour-matched controls; the
  // browser build keeps a normal full-width header. This class gates the difference.
  // The desktop app hides the native title bar; with the top bar now frozen, the window controls
  // always sit on the header, so no separate backdrop is needed.
  if (!IS_WEB) document.documentElement.classList.add('desktop-chrome');
  try { settings = (await window.api.getSettings()) || {}; } catch { settings = {}; }
  buildAssess();
  buildAssessment();                 // the new redesigned Assessment tab (both platforms)
  if (!IS_WEB) buildGeneBe();         // the original GeneBe assessment tab (desktop only)
  buildFusions();                    // BND fusion/translocation visualiser (offline, both platforms)
  buildBoard();
  buildAmino();
  buildCodon();
  buildTools();
  setupTabs();
  if (IS_WEB) {
    // GeneBe is desktop-only: its API blocks browser requests. Splice is NOT — the SpliceAI service
    // allows cross-origin calls, so the website scores splicing too; only SpliceVarDB inside that tab
    // is unavailable, and it says so itself.
    for (const v of ['genebe']) {
      const t = document.querySelector(`nav.tabs button[data-view="${v}"]`); if (t) t.remove();
    }
    const gbPill = $('#genebe-pill'); if (gbPill) gbPill.remove();
  } else {
    $('#genebe-pill').addEventListener('click', () => switchTab('genebe'));
    updateGenebePill();
  }
  // The transcript index (now the whole MANE store, ~1.9 MB) is only used by the Sources tab, so it is
  // loaded lazily there — startup fetches nothing extra.
  showBuildStamp();
  // Start blank — the assessor shows its neutral prompt until a variant is entered.
}

/**
 * Say which build this is, in the header.
 *
 * There are two live copies of Codon Compass — the application on this machine and the website — and
 * each is a copy of the source taken when it was built. A fix does not reach either until that copy is
 * rebuilt, and on 8 August 2026 both were a day behind while the source had a corrected classifier:
 * the packaged application still called the t(3;14) IGH-BCL6 event "No product" hours after it was
 * fixed, with nothing on screen to say so. This is that missing line.
 *
 * The stamp is written by the build (scripts/build_stamp.mjs). A copy running straight from source has
 * none, and says so rather than inventing one.
 */
async function showBuildStamp() {
  const pill = $('#build-pill');
  if (!pill) return;
  // The desktop page sits in src/renderer/, the website page at the root of dist-web.
  for (const path of ['../build-stamp.json', './build-stamp.json']) {
    try {
      const res = await fetch(path);
      if (!res.ok) continue;
      const s = await res.json();
      if (!s || !s.hash) continue;
      const when = new Date(s.builtAt);
      const built = Number.isNaN(when.getTime()) ? String(s.builtAt) : when.toLocaleString('en-GB');
      pill.textContent = `build ${s.hash}`;
      pill.title = `${s.target === 'website' ? 'Website build' : 'Desktop build'}, made ${built} from `
        + `${s.fileCount} source files. The identifier is a hash of that source: if it differs from `
        + 'another copy, the two are not running the same code.';
      return;
    } catch { /* try the next location */ }
  }
  pill.textContent = 'build: unstamped';
  pill.title = 'No build stamp found. This copy is running straight from the source tree rather than '
    + 'from a build, so there is nothing to compare against another copy.';
}

// Boot the whole app, unless a control test has imported this module to drive one function on a
// stubbed page. The test sets `window.__CC_TEST` before importing; the app itself never sets it.
if (!(typeof window !== 'undefined' && window.__CC_TEST === true)) init();

// ---- test surface ----
// Control tests (tests/rendererRaces.test.mjs) drive these real functions against a stubbed
// `window.api`, so a regression is caught in the shipped code rather than in a copy of it. The
// getters exist because the state they read is module-level and reassigned as the user works.
// Also hung on `window` so the running page can be driven from outside for a screenshot — the only
// honest way to check a layout is to look at it with real data in it.
if (typeof window !== 'undefined') { queueMicrotask(() => { window.__CC_T = __testing; }); }

export const __testing = {
  runAssess,
  beginAssessment,
  annotateGeneBe,
  renderSplice,
  pullAssessment,
  renderFusionEvents,
  renderGeneBe,
  renderAssessment,
  buildLollipop,
  fusionChrOrder,
  buildNonCodingCard,
  loadManePairing,
  dbnsfpIndexFor,
  forTranscript,
  forgetVariant,
  buildAssess,
  buildAmino,
  buildFusions,
  buildStructureCard,
  buildClinvarCard,
  drawBoard,
  verdictClass,
  numberingFrame,
  pairFusionExamples,
  knownMergedRows,
  knownEvidenceText,
  // The merged list is built once and cached in a module-level variable. A control that feeds it a
  // different whitelist has to be able to clear that cache first.
  resetKnownRows() { knownRows = null; knownListError = null; },
  get latestGeneBe() { return latestGeneBe; },
  set latestGeneBe(v) { latestGeneBe = v; },
  get lastVariant() { return lastVariant; },
  set lastVariant(v) { lastVariant = v; },
  get lastSplice() { return lastSplice; },
  set lastSplice(v) { lastSplice = v; },
  get assessState() { return assessState; },
  set assessState(v) { assessState = v; },
  get fusionGeneLoci() { return fusionGeneLoci; },
  set fusionGeneLoci(v) { fusionGeneLoci = v; },
  get fusionTables() { return fusionTables; },
  set fusionTables(v) { fusionTables = v; },
};
