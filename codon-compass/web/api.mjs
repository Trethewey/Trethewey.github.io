// Browser implementation of the window.api bridge.
//
// The Electron app provides window.api via preload.cjs -> IPC -> main process (file reads and
// network). This file provides the SAME 16-function surface for the hosted web build, using fetch()
// and browser APIs. The renderer and all of src/core are shared, byte-for-byte, between the two
// builds; only this adapter differs. See scripts/build_web.mjs.
//
// Offline-capable methods read static JSON served alongside the page. Network methods
// (live transcript fetch, 3D structures, GeneBe) are attempted directly and degrade with a clear
// message if the host's network or a proxy blocks them — the offline core keeps working.

import { gnomadLookup } from '../net/gnomad.mjs';

// Static data lives beside index.html, so all paths are RELATIVE (no leading slash): under
// christrethewey.dev/codon-compass/ they resolve to .../codon-compass/data/... correctly.
const DATA = 'data';

async function getJson(path) {
  const res = await fetch(path, { headers: { Accept: 'application/json' } });
  if (!res.ok) throw new Error(`HTTP ${res.status} for ${path}`);
  return res.json();
}

// A manifest lists every bundled transcript and which folder holds it, so getTranscript never has
// to probe with 404s and can compute the "other versions you do have" list.
let manifestPromise = null;
function manifest() {
  if (!manifestPromise) manifestPromise = getJson(`${DATA}/manifest.json`).catch(() => ({ byAccession: {} }));
  return manifestPromise;
}

const api = {
  async getIndex() {
    try { return await getJson(`${DATA}/transcripts/index.json`); } catch { return []; }
  },

  async getTranscript(accession) {
    const m = await manifest();
    const dir = m.byAccession && m.byAccession[accession];
    if (dir) {
      const doc = await getJson(`${DATA}/${dir}/${accession}.json`);
      return { doc, from: dir === 'transcripts' ? 'bundled' : 'reference', otherVersions: [] };
    }
    // Not held: report which other versions of the same base accession are available.
    const base = String(accession).replace(/\.\d+$/, '');
    const otherVersions = Object.keys((m && m.byAccession) || {})
      .filter((a) => a.replace(/\.\d+$/, '') === base);
    return { doc: null, from: null, otherVersions };
  },

  async getPathways() {
    try { return await getJson(`${DATA}/haemonc-pathways.json`); } catch { return { pathways: [], geneIndex: [] }; }
  },

  async getGeneContext(gene) {
    try {
      const entry = await getJson(`${DATA}/gene-context/${encodeURIComponent(gene)}.json`);
      return { entry, sources: entry && entry._sources ? entry._sources : null, unavailable: false };
    } catch {
      return { entry: null, sources: null, unavailable: false };
    }
  },

  async getWhitelist() {
    try { return await getJson(`${DATA}/haemonc-whitelist-variants.json`); } catch { return { variants: [], count: 0 }; }
  },

  async getGeneLoci() {
    try { return await getJson(`${DATA}/gene-loci.json`); } catch { return { genes: {} }; }
  },

  async getManePairing() {
    try { return await getJson(`${DATA}/mane-pairing.json`); }
    catch (e) { return { byRefSeq: {}, byEnsembl: {}, failed: true, reason: e.message }; }
  },

  async getFusionTables() {
    const grab = async (f) => { try { return await getJson(`${DATA}/${f}`); } catch { return null; } };
    return {
      exons: await grab('gene-exons.json'),
      cytobands: await grab('cytobands.json'),
      igLoci: await grab('ig-loci.json'),
      igGenes: await grab('ig-genes.json'),
    };
  },

  // ---- network features: attempt directly, degrade cleanly ----

  async fetchTranscript() {
    // Live NCBI fetch is not available in the hosted build (cross-origin / proxy). The offline
    // reference set already covers the bundled panel; anything outside it must be added to the set.
    throw new Error('Live transcript fetch is not available in the web version. This accession is not in the offline set.');
  },

  async structureAlphaFold(uniprot) {
    const acc = String(uniprot).trim().toUpperCase();
    try {
      const meta = await getJson(`https://alphafold.ebi.ac.uk/api/prediction/${acc}`);
      const url = meta && meta[0] && meta[0].pdbUrl;
      if (!url) throw new Error('no model');
      const res = await fetch(url);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const text = await res.text();
      return { format: 'pdb', text, source: `AlphaFold DB (${acc})`, predicted: true, url };
    } catch (e) {
      return { error: `Structure unavailable in the web version (${e.message}).` };
    }
  },

  async structureRcsb(pdbId) {
    const id = String(pdbId).trim().toUpperCase();
    try {
      const res = await fetch(`https://files.rcsb.org/download/${id}.pdb`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const text = await res.text();
      return { format: 'pdb', text, source: `RCSB PDB ${id}`, predicted: false, url: null };
    } catch (e) {
      return { error: `RCSB unavailable in the web version (${e.message}).` };
    }
  },

  structureOpenLocal() {
    // Read a local structure file through a hidden file input — no server involved.
    return new Promise((resolve) => {
      const input = document.createElement('input');
      input.type = 'file';
      input.accept = '.pdb,.cif,.ent';
      input.addEventListener('change', async () => {
        const file = input.files && input.files[0];
        if (!file) { resolve(null); return; }
        const text = await file.text();
        const format = file.name.toLowerCase().endsWith('.cif') ? 'cif' : 'pdb';
        resolve({ format, text, source: `Local file: ${file.name}`, predicted: null, url: null });
      }, { once: true });
      input.click();
    });
  },

  async appInfo() {
    return { versions: { web: '1' }, isPackaged: false, web: true, bundledData: DATA, cacheDir: null, structureCache: null };
  },

  async openExternal(url) { window.open(url, '_blank', 'noopener'); },

  async copyText(text) {
    try { await navigator.clipboard.writeText(String(text || '')); }
    catch { /* clipboard blocked; ignore, the UI shows its own confirmation */ }
  },

  async getSettings() {
    try { return JSON.parse(localStorage.getItem('codon-compass-settings') || '{}'); } catch { return {}; }
  },

  async setSettings(patch) {
    let cur = {};
    try { cur = JSON.parse(localStorage.getItem('codon-compass-settings') || '{}'); } catch { /* reset */ }
    const next = { ...cur, ...(patch || {}) };
    try { localStorage.setItem('codon-compass-settings', JSON.stringify(next)); } catch { /* storage full/blocked */ }
    return next;
  },

  async genebeAnnotate() {
    // GeneBe's API does not permit browser cross-origin calls, and it would send the variant off the
    // page. It is deliberately unavailable in the hosted build.
    return { ok: false, error: 'GeneBe annotation is available only in the desktop app.' };
  },

  async genebeWhoami() {
    return { ok: false, error: 'GeneBe is available only in the desktop app.' };
  },

  async gnomad(hgvs, knownCoord) {
    // gnomAD and Ensembl both allow cross-origin requests, so this works from the browser where the
    // network permits it, and degrades with a clear message where a proxy blocks it.
    try { const res = await gnomadLookup(hgvs, knownCoord || null); return { ok: true, cached: false, ...res }; }
    catch (e) { return { ok: false, error: e.message }; }
  },

  // The Assessment/Splice network annotations are desktop-only; the tabs are hidden in the web build,
  // but keep safe stubs so nothing can crash if a method is called.
  async genebeAnnotate2() { return { ok: false, error: 'Desktop app only.' }; },
  async clinvar() { return { ok: false, error: 'Desktop app only.' }; },
  async clinvarSubmissions() { return { ok: false, error: 'Desktop app only.' }; },
  async resolveCoord() { return { ok: false, error: 'Desktop app only.' }; },
  // SpliceAI IS reachable from a browser — the service sends a permissive cross-origin header — so the
  // website scores splicing exactly as the desktop application does. It needs a genomic coordinate;
  // for a coding variant that comes from the gnomAD lookup, for an intronic one from the bundled exon
  // table. Only the coordinate is sent, never anything else about the case.
  async spliceai(hgvs, transcript, knownCoord) {
    if (!knownCoord) return { ok: false, error: 'No genomic coordinate for this variant yet.' };
    try {
      const { spliceaiLookup } = await import('../net/spliceai.mjs');
      const res = await spliceaiLookup(knownCoord, transcript);
      return { ok: true, cached: false, ...res };
    } catch (e) { return { ok: false, error: e.message }; }
  },
  // SpliceVarDB sends no cross-origin header, so a browser cannot read its reply. Said plainly rather
  // than left as a silent absence.
  async splicevardb() {
    return { ok: false, error: 'SpliceVarDB does not allow browser requests, so it is available in the desktop application only.' };
  },
};

window.api = api;
window.__CC_WEB = true;   // renderer hides GeneBe (its API blocks browser requests) and adjusts labels
export default api;
