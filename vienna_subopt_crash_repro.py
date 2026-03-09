"""
Minimal reproducer for ViennaRNA 2.6.4 density_of_states[-1] out-of-bounds write.

SUMMARY
-------
Calling vrna_subopt_cb (or the Python equivalent fc.subopt()) on a cofold sequence
of a short mRNA fragment and an anti-Shine-Dalgarno rRNA sequence causes an
out-of-bounds write to density_of_states[-1] in ViennaRNA 2.6.4.

The word immediately before density_of_states[] in the BSS segment is the global
pointer `last_param_file`.  Writing +1 to it corrupts the pointer to a nonsense
address.  The *next* call to vrna_fold_compound() then crashes in strncpy() when
it attempts to copy from that invalid address.

AFFECTED VERSION : ViennaRNA 2.6.4 (confirmed crash)
FIXED IN         : ViennaRNA 2.7.x (correction term added; see Root Cause below)

ROOT CAUSE (src/ViennaRNA/subopt/subopt.c)
------------------------------------------
After backtracking a suboptimal structure, the code maps its energy to a
density_of_states bucket:

    # 2.6.4 (buggy) -- no lower-bound guard:
    e = (int)((structure_energy - min_en) * 10.)
    if e > MAXDOS:
        e = MAXDOS
    density_of_states[e] += 1          # crashes when e == -1

    # master / 2.7.x (partially fixed) -- uses a correction offset:
    correction = -0.1 if min_en < 0 else 0.1
    e = int((structure_energy - min_en) * 10. - correction)
    if e > MAXDOS:
        e = MAXDOS
    density_of_states[e] += 1          # still no lower-bound guard

In 2.6.4, the cofold backtracking accumulates a partial_energy (integer units of
0.01 kcal/mol) that can land below the stored MFE by enough that the conversion
to a 0.1 kcal/mol bucket index produces -1.  Without clamping, this writes one
word before the array.

SUGGESTED FIX
-------------
Add a lower-bound guard immediately after the upper-bound guard:

    if e < 0:
        e = 0
    elif e > MAXDOS:
        e = MAXDOS

CRASH SEQUENCE IN OSTIR
-----------------------
OSTIR's Rust binary (using ViennaRNA 2.6.4 via librna-sys) called:

    Step 1: vrna_subopt_cb("AUAAGGAGGU&ACCUCCUUA", delta=548, ...)
              -> density_of_states[-1]++  (last_param_file corrupted)
    Step 2: vrna_fold_compound("A", ...)
              -> strncpy(buf, last_param_file, ...)  -> SIGSEGV

The sequences are derived from:
    mRNA fragment : ATAAGGAGGT  (10 nt preceding the AUG start codon)
    ASD (rRNA 3') : ACCTCCTTA   (ViennaRNA's default anti-Shine-Dalgarno)
After T->U conversion and joining with '&':
    "AUAAGGAGGU&ACCUCCUUA"
"""

import sys

try:
    import RNA
except ImportError:
    sys.exit("ERROR: ViennaRNA Python bindings not found.  Install with: pip install viennarna")

print(f"ViennaRNA version: {RNA.__version__}")

# ---------------------------------------------------------------------------
# Model parameters matching OSTIR's defaults
# ---------------------------------------------------------------------------
md = RNA.md()
md.noLP = 1          # no lonely/isolated base pairs  (OSTIR default)
md.dangles = 2       # all dangles                    (OSTIR default for short sequences)
md.temperature = 37.0

# ---------------------------------------------------------------------------
# Step 1 – the triggering call
# ---------------------------------------------------------------------------
# These two strands are what OSTIR passes to vrna_subopt_cb when computing
# dG_mRNA_rRNA for the sequence "ATAAGGAGGTATG" with an ATG start codon at
# position 11 (1-indexed).
mrna = "ATAAGGAGGT"
asd  = "ACCTCCTTA"  # default anti-Shine-Dalgarno (3'->5' of 16S rRNA)

# ViennaRNA uses RNA alphabet with '&' as strand separator
seq = (mrna + "&" + asd).replace("T", "U").upper()
print(f"\nSequence passed to vrna_subopt: {seq!r}")

# Energy delta in OSTIR: ((3.0 + 2.481) * 100) rounded = 548 centi-kcal/mol
energy_delta = 548
print(f"Energy delta: {energy_delta} centi-kcal/mol ({energy_delta/100:.2f} kcal/mol)")

fc = RNA.fold_compound(seq, md)
(mfe_struct, mfe_en) = fc.mfe()
print(f"MFE: {mfe_struct}  dG = {mfe_en:.2f} kcal/mol")

# On ViennaRNA 2.6.4 this call writes density_of_states[-1]++ and
# corrupts the adjacent last_param_file global pointer.
# On 2.7.x it runs cleanly (the correction term was added).
print(f"\nCalling fc.subopt({energy_delta}) ...")
results = fc.subopt(energy_delta)
print(f"  -> returned {len(results)} structures (no crash on this version)")
for r in sorted(results, key=lambda x: x.energy):
    print(f"     {r.structure}  dG={r.energy:.2f}")

# ---------------------------------------------------------------------------
# Step 2 – the call that segfaults in 2.6.4 (reads the corrupted pointer)
# ---------------------------------------------------------------------------
# Any subsequent vrna_fold_compound call would crash.  We use "AU" (the
# shortest valid RNA) to keep the example minimal.
print(f"\nCalling vrna_fold_compound('AU', md) ...")
print("  On ViennaRNA 2.6.4 this SEGFAULTS inside strncpy()")
print("  because last_param_file was corrupted in the previous subopt call.")
fc2 = RNA.fold_compound("AU", md)
(ss2, mfe2) = fc2.mfe()
print(f"  -> structure={ss2!r}  dG={mfe2:.2f}  (clean on this version)")

print("\n=== Reproducer complete – no crash on ViennaRNA", RNA.__version__, "===")
print("Run this script with ViennaRNA 2.6.4 to observe the SIGSEGV.")
