// Example breakend sets for the Fusions tab.
//
// Kept in the core, not the renderer, so `src/tests/fusion.test.mjs` can assert that each worked
// example still reads the way it is advertised to. An example that quietly stops demonstrating its
// point is worse than no example.

/** Six real rows, verbatim from the bundled structural-variant report. */
export const FUSION_REPORT_EXAMPLE = [
  '3:187746165T]12:25056225]\t3q27.3;12p12.1\tBND\tN/A\t-\tPR-40/194;SR-23/159',
  '3:187746173[12:25056242[T\t3q27.3;12p12.1\tBND\tN/A\t-\tPR-35/196;SR-26/161',
  '14:105864633ACCC[18:63126216[\t14q32.33;18q21.33\tBND\tN/A\t0.0001\tPR-28/150;SR-17/146',
  '14:105891702]18:63126209]GGG(+4)CTC\t14q32.33;18q21.33\tBND\tN/A\t-\tPR-37/178;SR-22/153',
  '2:15726277C]22:22343454]\t2p24.3;22q11.22\tBND\tN/A\t-\tPR-13/134',
  '12:124989860]16:28915184]A\t12q24.31;16p11.2\tBND\tN/A\t0.6408 (GE)\tPR-8/59',
].join('\n');

/**
 * A CONSTRUCTED teaching set — not patient data, and not lifted from any published case. Each row
 * exercises a different product class, so the tab can be understood without a real report to hand.
 *
 * No breakpoint here is invented. Every coordinate is derived from the bundled MANE Select spans in
 * data/gene-loci.json by the rule named beside it, and the two immunoglobulin-locus coordinates are
 * real breakpoints lifted from the report example above:
 *
 *   22:23,249,273   midpoint of BCR      23,180,509–23,318,037
 *    9:130,861,464  midpoint of ABL1    130,835,254–130,887,675
 *   12:11,772,525   midpoint of ETV6     11,649,674–11,895,377
 *   21:34,918,551   midpoint of RUNX1    34,787,801–35,049,302
 *    8:127,734,231  2 kb below MYC's span start (127,736,231)
 *   16:10,932,642   midpoint of the region where CIITA and DEXI overlap
 *   11:102,328,443  midpoint of BIRC3   102,317,484–102,339,403
 *   18:58,712,971   midpoint of MALT1    58,671,465–58,754,477
 *   17:7,677,955    midpoint of TP53      7,668,421–7,687,490
 *   14:105,864,633 and 2:15,726,277 — real, from the report example
 *
 * Cytobands are given to arm level only, which is all the derivative naming needs. Read-support and
 * germline figures are deliberately absent: inventing them would be worse than letting the panel say
 * they were not in the pasted row.
 */
export const FUSION_WORKED_EXAMPLES = [
  '22:23249273T[9:130861464[\t22q11;9q34\tBND\tBCR-ABL1 chimeric candidate, both plus strand',
  '12:11772525T]21:34918551]\t12p13;21q22\tBND\tETV6-RUNX1 chimeric candidate, opposite strands',
  '14:105864633A[8:127734231[\t14q32;8q24\tBND\tunnamed IG-locus donor into an intact MYC',
  '16:10932642T]14:105864633]\t16p13;14q32\tBND\tCIITA into an IG locus, breakpoint in two overlapping genes',
  '11:102328443T]18:58712971]\t\tBND\tBIRC3 and MALT1 joined the way round that cannot work, and no cytoband',
  '2:15726277T[17:7677955[\t2p24;17p13\tBND\tTP53 minus-strand donor, half the gene lost',
].join('\n');

/**
 * What each worked example is there to show. The control test asserts the tool still produces this,
 * so the set cannot drift into demonstrating something else.
 *   derivative: expected der() name, or null where the row deliberately carries no cytoband
 */
export const FUSION_WORKED_EXPECTATIONS = [
  { cell: '22:23249273T[9:130861464[', derivative: 'der(22)', kind: 'chimeric-candidate', reads: ['BCR', 'ABL1'],
    shows: 'the classic 5′-into-3′ chimera, both partners on the plus strand and joined in the same orientation' },
  { cell: '12:11772525T]21:34918551]', derivative: 'der(21)', kind: 'chimeric-candidate', reads: ['ETV6', 'RUNX1'],
    shows: 'the same shape with the partners on OPPOSITE strands and a reverse-complement join — it still reads forward, which is the case the old code got backwards' },
  { cell: '14:105864633A[8:127734231[', derivative: 'der(14)', kind: 'unnamed-donor-intact', reads: [null, 'MYC'],
    shows: 'an unnamed immunoglobulin-locus donor placed upstream of a gene that survives whole' },
  { cell: '16:10932642T]14:105864633]', derivative: 'der(14)', kind: 'deregulation-candidate', reads: ['CIITA', null],
    shows: 'a breakpoint inside two overlapping genes, where the alternative reading has to be surfaced' },
  { cell: '11:102328443T]18:58712971]', derivative: null, kind: 'no-read-through', reads: ['BIRC3', 'MALT1'],
    shows: 'two genes pointing at each other, so nothing transcribes across — and, with no cytoband, a derivative that stays unnamed rather than being guessed' },
  { cell: '2:15726277T[17:7677955[', derivative: 'der(17)', kind: 'truncated-donor', reads: ['TP53', null],
    shows: 'a minus-strand donor cut in half, kept apart from the near-whole case by the retained-span figure' },
];
