// Codon Compass — the pass/fail rules for the automated smoke run, and the crash-safe
// wrapper that runs it.
//
// Both live outside the Electron main process so `node --test` can exercise them. A gate
// that has quietly stopped gating is invisible otherwise: the smoke keeps printing PASS
// while the feature it was meant to protect is broken. The same goes for the crash safety
// at the foot of this file — main.mjs imports Electron, so nothing in it can be tested.
//
// Every gate must be provably present, not merely provably watching a field: each one is
// paired with a breakage in tests/smokeHarness.test.mjs that only it catches, so deleting a
// gate leaves that breakage uncaught and the controls go red naming it.
//
// What is gated, and what is not
// ------------------------------
// Gated: every check whose correct value is fixed and produced on this machine — the
// genetic-code and amino-acid tables, the bundled transcripts and their protein domains,
// the whole-genome reference store, the bundled pathway and gene-context files, the
// frameshift read-through path, the lollipop board, and the Fusions tab.
//
// Not gated: anything that needs a public web service (GeneBe, ClinVar, SpliceAI,
// SpliceVarDB, gnomAD, AlphaFold). This application is meant to work offline, so an
// offline run must not be reported as a failure. Those fields are still recorded in the
// report for a human to read.
//
// Screenshots sit between the two. A single lost image is not a failure: capturePage() throws
// a transient graphics error every few runs, and a lost picture does not change what the
// application computed, so those are listed as warnings with a count of how many of the run's
// attempts reached the disk. Saving NO image at all is different — the pictures are the
// artefacts a person actually reviews, and a run that produced none of them has not been
// reviewed, so that is gated. A picture that was taken and then could not be written is not a
// graphics flake at all (a full disk, a locked directory) and fails the run outright.

/** One rule. `fields` are the summary keys quoted back when it fails. */
const gate = (name, fields, expected, ok) => ({ name, fields, expected, ok });

const num = (v) => (typeof v === 'number' && Number.isFinite(v) ? v : NaN);
const text = (v) => (typeof v === 'string' ? v.trim() : '');
const has = (v, re) => re.test(text(v));

/**
 * The deterministic, offline checks. Every one of these must hold on a machine with no
 * network at all; if one cannot, it does not belong here.
 */
export const OFFLINE_GATES = [
  // --- the page itself ---
  gate('window title', ['title'], "the title is 'Codon Compass'",
    (s) => text(s.title) === 'Codon Compass'),
  gate('run reached the end', ['completed'], 'the smoke finished every step',
    (s) => s.completed === true),
  // Zero images means the capture path itself is broken — a graphics driver that never
  // recovers, an output directory that cannot be written — and the run has no artefacts for
  // anyone to look at. A run that lost some of its images still passes and says how many.
  gate('screenshots were saved', ['screenshotsWritten', 'screenshotsAttempted'],
    'at least one screenshot reached the disk',
    (s) => num(s.screenshotsWritten) >= 1),

  // --- assessment of the example variant, ARID1A NM_006015.6:c.461A>C ---
  gate('assessment rendered', ['hasConsequence'], 'a consequence badge is shown',
    (s) => s.hasConsequence === true),
  gate('assessment badge', ['badge'], 'the change is called a missense',
    (s) => has(s.badge, /missense/i)),
  gate('ARID1A protein change', ['proteinShort'], 'Y154S',
    (s) => has(s.proteinShort, /Y154S/)),
  gate('ARID1A domain legend', ['arid1aDomains'], 'at least one protein domain drawn',
    (s) => num(s.arid1aDomains) >= 1),

  // --- reference tables: fixed by biology, not by data files ---
  gate('amino-acid table', ['aminoRows'], '20 rows, one per amino acid',
    (s) => num(s.aminoRows) === 20),
  gate('amino-acid chemical structures', ['aminoStructures'], '20 drawn structures',
    (s) => num(s.aminoStructures) === 20),
  gate('codon grid', ['codonCells'], '64 cells, one per codon',
    (s) => num(s.codonCells) === 64),
  gate('codon sortable table', ['codonTableRows'], '64 rows, one per codon',
    (s) => num(s.codonTableRows) === 64),

  // --- the canonical known-variant picker ---
  gate('known-variant picker', ['knownCanonOptions'], '50 canonical variants',
    (s) => num(s.knownCanonOptions) === 50),
  gate('known-variant ordering', ['knownCanonTop'], 'KRAS G12D first, by sample count',
    (s) => has(s.knownCanonTop, /KRAS G12D/)),
  gate('MYD88 previous name surfaced', ['knownMyd88Alias'], 'the L252P entry shows "previously L265P"',
    (s) => s.knownMyd88Alias === true),
  gate('MYD88 search hits', ['knownMyd88Hits'], 'exactly one hit for MYD88 L252P',
    (s) => num(s.knownMyd88Hits) === 1),

  // --- BRAF V600E: bundled transcript plus bundled domains ---
  gate('BRAF protein change', ['brafProtein'], 'V600E',
    (s) => text(s.brafProtein) === 'V600E'),
  gate('BRAF domain caption', ['brafDomainCaption'], 'residue 600 reported inside a domain',
    (s) => has(s.brafDomainCaption, /Residue 600 lies in/)),

  // --- EZH2 Y646N: domains plus the MANE-to-UniProt residue mapping ---
  gate('EZH2 domain legend', ['ezh2Domains'], 'at least one protein domain drawn',
    (s) => num(s.ezh2Domains) >= 1),
  gate('EZH2 domain caption', ['ezh2DomainCaption'], 'residue 646 reported inside a domain',
    (s) => has(s.ezh2DomainCaption, /Residue 646 lies in/)),
  gate('EZH2 residue mapping', ['ezh2StructNote'], 'MANE residue 646 maps to UniProt residue 641',
    (s) => has(s.ezh2StructNote, /MANE residue 646 corresponds to UniProt residue 641/)),

  // --- MPL W515L: only the whole-genome reference store can resolve this one ---
  gate('MPL resolved from the reference store', ['mplProtein'], 'the heading names MPL and W515L',
    (s) => has(s.mplProtein, /MPL/) && has(s.mplProtein, /W515L/)),
  gate('MPL domain legend', ['mplDomains'], 'at least one protein domain drawn',
    (s) => num(s.mplDomains) >= 1),

  // --- drawn pathway figures and gene context, all from bundled files ---
  gate('pathway figure drawn', ['pathwayNodes', 'pathwayEdges'], 'at least one node and one edge',
    (s) => num(s.pathwayNodes) >= 1 && num(s.pathwayEdges) >= 1),
  gate('pathway mechanism named', ['pathwayMechanism'], 'a mechanism chip with text',
    (s) => text(s.pathwayMechanism).length > 0),
  gate('gene function text (MPL)', ['pathwayFunction'], 'the "What it does" text is present',
    (s) => text(s.pathwayFunction).length > 0),
  gate('gene function text (CFTR, non-haematological)', ['nonHaemFunction'], 'the "What it does" text is present',
    (s) => text(s.nonHaemFunction).length > 0),
  gate('ARID1A pathway figure', ['arid1aFigureNodes'], 'a drawn figure, at least one node',
    (s) => num(s.arid1aFigureNodes) >= 1),
  gate('ARID1A shows a figure, not a chip list', ['arid1aNoPathwayChipCard'], 'no pathway-name chip card',
    (s) => s.arid1aNoPathwayChipCard === true),

  // --- the frameshift path: ASXL1 NM_015338.5:c.2945dup ---
  // L983Afs*8 is what src/core/variant.mjs resolveIndel computes from the bundled coding
  // sequence read through into the 3-prime untranslated region; it is not a value typed in
  // by hand. tests/indelControl.test.mjs asserts the same change from the core side.
  gate('ASXL1 frameshift consequence', ['indelProtein'], 'L983Afs*8',
    (s) => text(s.indelProtein) === 'L983Afs*8'),
  gate('ASXL1 frameshift badge', ['indelBadge'], 'the change is called a frameshift',
    (s) => has(s.indelBadge, /frameshift/i)),

  // --- the lollipop board ---
  gate('board drew diagrams', ['boardDiagrams'], 'at least one diagram',
    (s) => num(s.boardDiagrams) >= 1),
  gate('board drew variant markers', ['boardMarkers'], 'at least one marker',
    (s) => num(s.boardMarkers) >= 1),
  gate('board status matches what was drawn', ['boardStatus', 'boardDiagrams'],
    'the status line counts the same number of diagrams as were drawn',
    (s) => {
      const m = /^\s*(\d+)\s+diagram/.exec(text(s.boardStatus));
      return Boolean(m) && Number(m[1]) === num(s.boardDiagrams);
    }),

  // --- Fusions: BCR::ABL1, then IGH::BCL2, then the t(3;14) IGH-BCL6 intron-1 break ---
  gate('fusion panel drew without error', ['fusionNoError', 'fusionUnrecognised'],
    'no draw error and the breakend row was recognised',
    (s) => s.fusionNoError === true && s.fusionUnrecognised === false),
  gate('fusion panel count', ['fusionPanels', 'fusionLedgerRows', 'fusionSectionsRendered'],
    'one panel, one ledger row, one section heading',
    (s) => num(s.fusionPanels) === 1 && num(s.fusionLedgerRows) === 1 && num(s.fusionSectionsRendered) === 1),
  gate('fusion detail folds', ['fusionFoldsRendered'], 'at least one expandable detail section',
    (s) => num(s.fusionFoldsRendered) >= 1),
  gate('BCR::ABL1 ledger row', ['fusionLedgerTop'], 'der(22) naming BCR',
    (s) => has(s.fusionLedgerTop, /der\(22\)/) && has(s.fusionLedgerTop, /BCR/)),
  gate('BCR::ABL1 driver verdict', ['fusionLedgerTopVerdict'], 'the driver tag',
    (s) => has(s.fusionLedgerTopVerdict, /v-driver/)),
  gate('BCR::ABL1 reading frame', ['fusionInFrame'], 'at least one headline reporting an in-frame junction',
    (s) => num(s.fusionInFrame) >= 1),
  gate('Philadelphia chromosome named', ['fusionPhiladelphia'], 'the panel names it',
    (s) => s.fusionPhiladelphia === true),
  gate('picking a second fusion replaces the view', ['fusionPickReplaces'], 'still one panel, not two',
    (s) => num(s.fusionPickReplaces) === 1),
  gate('IGH locus named', ['fusionIghNamed'], 'the immunoglobulin heavy chain locus is named',
    (s) => s.fusionIghNamed === true),
  // The t(3;14) case: BCL6's start codon sits in exon 3, so an intron-1 break leaves the
  // whole reading frame intact and the IGH locus supplies the promoter. Reported as
  // "No product" until 2026-08-08.
  gate('t(3;14) IGH-BCL6 is not called "No product"', ['fusionBcl6UtrNoProduct'],
    'no "No product" text',
    (s) => s.fusionBcl6UtrNoProduct === false),
  gate('t(3;14) headline names BCL6', ['fusionBcl6UtrHeadline'], 'BCL6 in the headline',
    (s) => has(s.fusionBcl6UtrHeadline, /BCL6/)),
  gate('t(3;14) called a promoter substitution', ['fusionBcl6UtrLedger'], 'promoter substitution in the ledger',
    (s) => has(s.fusionBcl6UtrLedger, /promoter substitution/)),
  gate('Clear empties input and output', ['fusionCleared'], 'no panels and an empty box',
    (s) => s.fusionCleared === true),
];

/**
 * Log prefixes that mean the run itself went wrong. A capture failure is deliberately not
 * here: screenshots are flaky and do not change what was computed. A screenshot that was
 * taken and could not then be SAVED is here, because that is the file system failing, and
 * whatever caused it will lose every later image too.
 */
export const FATAL_LOG_PREFIXES = [
  '[smoke-error]', '[smoke-fatal]', '[smoke-timeout]', '[render-gone]', '[preload-error]',
  '[screenshot-write-failed]',
];

/**
 * Log prefixes worth reporting but not worth failing on: a lost picture, and a wait that
 * reached its deadline. A wait that gave up is a warning rather than a failure because most
 * of them are on network answers, which an offline run never gets; if the thing waited for
 * was one this application computes for itself, its own gate fails a moment later and names
 * the field.
 */
export const WARNING_LOG_PREFIXES = ['[capture-failed]', '[capture-retry]', '[wait-timeout]'];

function seen(summary, fields) {
  return fields.map((f) => `${f}=${JSON.stringify(summary[f])}`).join(', ');
}

/** The verdict written before the run starts, so a crash cannot leave an older PASS on disk. */
export const NO_VERDICT_FAILURE =
  'the run did not finish: no verdict was reached (crash, hang or forced exit)';

/** The verdict written when the watchdog ends a hang. */
export const WATCHDOG_FAILURE = 'the run hung and was stopped by the watchdog';

/**
 * The record left in the project's usual output directory when a run was sent somewhere else
 * with --smoke-out=. Without it the old verdict stays there looking current, which is how a
 * twelve-day-old PASS survived two fresh runs.
 *
 * `screenshots` is empty because this notice is not a run and took no pictures. The pictures
 * the run before it left behind are named separately, so the next run in this directory can
 * still clear them: without that list they were orphaned, and a picture from two runs ago sat
 * in the directory looking like one of the current set.
 */
export function supersededNotice({ writtenTo, at, previousResult = null, previousScreenshots = [] }) {
  return {
    ok: false,
    result: 'SUPERSEDED',
    failures: [`this is not a verdict: the run of ${at} was written to ${writtenTo}, so no result for it exists here`],
    warnings: [],
    screenshots: [],
    summary: {},
    logs: [],
    supersededAt: at,
    resultWrittenTo: writtenTo,
    previousResult,
    previousScreenshots,
  };
}

/**
 * Apply every gate to a smoke summary.
 *
 * @param {object} summary the fields the smoke read from the live page
 * @param {string[]} logs the run's log lines
 * @returns {{ok:boolean, failures:string[], warnings:string[]}}
 */
export function evaluateSmoke(summary = {}, logs = []) {
  const s = summary || {};
  const lines = Array.isArray(logs) ? logs : [];
  const failures = [];
  for (const g of OFFLINE_GATES) {
    let passed = false;
    try { passed = Boolean(g.ok(s)); } catch { passed = false; }
    if (!passed) failures.push(`${g.name}: expected ${g.expected}; saw ${seen(s, g.fields)}`);
  }
  for (const line of lines) {
    if (FATAL_LOG_PREFIXES.some((p) => line.startsWith(p))) failures.push(`log: ${line}`);
  }
  const warnings = lines.filter((line) => WARNING_LOG_PREFIXES.some((p) => line.startsWith(p)));
  return { ok: failures.length === 0, failures, warnings };
}

// ---------------------------------------------------------------------------
// Waiting for the page to be ready, instead of guessing how long it takes
// ---------------------------------------------------------------------------
//
// The run used to advance on fixed sleeps — about a hundred seconds of them. A sleep is a
// guess about the machine: too short and the next step reads the page before it has finished
// drawing, so a field the gates check comes back empty and a healthy build is reported as a
// failure; too long and every run pays for the worst machine anyone might use it on.
//
// So each step now states what it is waiting FOR and how long it is prepared to wait. On a
// quick machine it goes on at the first poll; on a slow one it waits. The deadline is not a
// delay: reaching it does not fail the run by itself, it just stops the wait, and the gates
// then judge whatever the page really shows.

/**
 * Poll until something is true, or the deadline passes.
 *
 * The clock and the sleep are passed in so this can be driven by a fake clock in
 * tests/smokeHarness.test.mjs rather than by making the controls wait in real time.
 *
 * @param {() => (boolean|Promise<boolean>)} probe  asked repeatedly; a throw counts as "not yet"
 * @param {{timeoutMs:number, intervalMs?:number, sleep:(ms:number)=>Promise<void>, now:()=>number}} opts
 * @returns {Promise<{ok:boolean, timedOut:boolean, waitedMs:number, polls:number}>}
 */
export async function waitUntil(probe, { timeoutMs, intervalMs = 150, sleep, now }) {
  const startedAt = now();
  const waited = () => now() - startedAt;
  let polls = 0;
  for (;;) {
    let ready = false;
    // A probe that throws means the page is mid-render and the element is not there yet.
    // That is a "not yet", never an error to report.
    try { ready = Boolean(await probe()); } catch { ready = false; }
    polls += 1;
    if (ready) return { ok: true, timedOut: false, waitedMs: waited(), polls };
    const left = timeoutMs - waited();
    if (left <= 0) return { ok: false, timedOut: true, waitedMs: waited(), polls };
    await sleep(Math.min(intervalMs, left));
  }
}

// ---------------------------------------------------------------------------
// The crash-safe wrapper around a smoke run
// ---------------------------------------------------------------------------
//
// This lives here, not in main.mjs, for one reason: main.mjs imports Electron, so nothing in
// it can be driven by `node --test`. The pieces that decide what ends up on disk — the
// failure written before the run starts, the watchdog that ends a hang, the verdict written
// whatever happens, and the exit status — are therefore kept as one plain function with
// every outside effect passed in. tests/smokeHarness.test.mjs drives it with a fake file
// store and a fake timer.
//
// Nothing here touches the file system, the clock or the process directly.

/** A file name a run may delete: a plain name in its own output directory, nothing else. */
const isPlainFileName = (n) => typeof n === 'string' && n.length > 0 && !/[\\/]/.test(n) && n !== '.' && n !== '..';

/**
 * Run a smoke test with a verdict on disk at every moment.
 *
 * @param {object} deps
 *   - resultFile      where this run's verdict goes
 *   - defaultResultFile  the project's usual verdict file; if this run was sent elsewhere and
 *                        a verdict already sits there, it is replaced by a "superseded" record
 *   - mkdir, writeFile, readFile, remove, exists  file operations, all awaited
 *   - joinOut(name)   full path of an output file from its bare name
 *   - now()           an ISO timestamp
 *   - setTimer/clearTimer, exit(code), log(line)
 *   - timeoutMs, isPackaged
 * @param {(ctx:{summary:object, logs:string[], screenshots:string[], registerCleanup:Function}) => Promise<void>} body
 *   the steps that drive the application; it fills `summary` and pushes to `logs`.
 * @returns {Promise<object|null>} the report written, or null if even that failed
 */
export async function runGuardedSmoke(deps, body) {
  const {
    resultFile, defaultResultFile = null,
    mkdir, writeFile, readFile, remove, exists, joinOut,
    now, setTimer, clearTimer, exit, log,
    timeoutMs, isPackaged = false,
  } = deps;

  const startedAt = now();
  const logs = [];
  const screenshots = [];
  const cleanups = [];
  // `completed` is a gated field: it turns true only after the last step, so a run that stops
  // early cannot pass however good the fields it managed to collect look.
  const summary = { completed: false };

  const buildReport = (verdict) => ({
    ok: Boolean(verdict.ok),
    result: verdict.ok ? 'PASS' : 'FAIL',
    failures: verdict.failures || [],
    warnings: verdict.warnings || [],
    screenshots,
    summary,
    logs,
    // Which rules produced this verdict. A report that lists fewer rules than the current
    // set was written by older code, and saying so on the face of the file is the only way
    // a reader can tell a current result from one left behind.
    gates: { count: OFFLINE_GATES.length, names: OFFLINE_GATES.map((g) => g.name) },
    isPackaged,
    startedAt,
    at: now(),
  });

  const writeReport = async (verdict) => {
    const report = buildReport(verdict);
    await writeFile(resultFile, JSON.stringify(report, null, 2));
    return report;
  };

  await mkdir();

  // Remove the images the previous run recorded writing, so a picture this run fails to
  // capture cannot be mistaken for one it took. Only names that run listed, only bare file
  // names, only in this directory — never a pattern.
  try {
    const prior = JSON.parse(await readFile(resultFile));
    // Both lists: a normal report names the pictures that run took, and a superseded notice
    // names the ones the run before it took, which nothing else would ever clear.
    const stale = [...(prior.screenshots || []), ...(prior.previousScreenshots || [])];
    for (const name of stale.filter(isPlainFileName)) await remove(joinOut(name));
  } catch { /* no readable previous report: nothing of its to clear */ }

  await remove(resultFile);
  await writeReport({ ok: false, failures: [NO_VERDICT_FAILURE] });

  // A run sent elsewhere must not leave the usual place holding an older verdict.
  if (defaultResultFile && defaultResultFile !== resultFile) {
    try {
      if (await exists(defaultResultFile)) {
        let previousResult = null;
        let previousScreenshots = [];
        try {
          const prev = JSON.parse(await readFile(defaultResultFile));
          previousResult = prev.result || null;
          previousScreenshots = [...(prev.screenshots || []), ...(prev.previousScreenshots || [])];
        } catch { /* unreadable */ }
        const notice = supersededNotice({ writtenTo: resultFile, at: startedAt, previousResult, previousScreenshots });
        await writeFile(defaultResultFile, JSON.stringify(notice, null, 2));
      }
    } catch (e) {
      logs.push(`[smoke-error] could not mark the usual verdict file superseded: ${e.message}`);
    }
  }

  let timer = null;
  let report = null;
  try {
    timer = setTimer(() => {
      logs.push(`[smoke-timeout] the run passed ${Math.round(timeoutMs / 1000)} seconds and was stopped`);
      Promise.resolve(writeReport({ ok: false, failures: [WATCHDOG_FAILURE] }))
        .catch((e) => log('SMOKE could not write its timeout report: ' + e.message))
        .finally(() => { log('SMOKE_RESULT FAIL (timed out)'); exit(1); });
    }, timeoutMs);

    await body({ summary, logs, screenshots, registerCleanup: (fn) => cleanups.push(fn) });
  } catch (e) {
    // Anything the run threw — a page that never loaded, a missing element in the very first
    // step — lands here, and the verdict below is still written.
    logs.push(`[smoke-error] ${e && e.message ? e.message : String(e)}`);
  } finally {
    if (timer !== null) clearTimer(timer);
    try {
      report = await writeReport(evaluateSmoke(summary, logs));
    } catch (e) {
      log('SMOKE could not write its report: ' + e.message);
    }
    const ok = Boolean(report && report.ok);
    log('SMOKE_RESULT ' + (report ? report.result : 'FAIL') + ' ' + JSON.stringify(summary));
    if (report && report.failures.length) log('SMOKE_FAILURES\n  ' + report.failures.join('\n  '));
    for (const fn of cleanups) { try { await fn(); } catch { /* cleanup is best effort */ } }
    exit(ok ? 0 : 1);
  }
  return report;
}
