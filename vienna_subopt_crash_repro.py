"""
Minimal reproducer for ViennaRNA density_of_states[-1] out-of-bounds write.

SUMMARY
-------
Calling vrna_subopt_cb (or the Python equivalent fc.subopt()) on a cofold sequence
of a short mRNA fragment and an anti-Shine-Dalgarno rRNA sequence can write past
the beginning of the global array density_of_states[MAXDOS+1].

The root cause — missing lower-bound guard on the bucket index — is present in
ALL versions including 2.7.2.  In 2.6.4 the write corrupted the global pointer
last_param_file (located immediately before the array in BSS), causing a SIGSEGV
on the very next call to vrna_fold_compound().  In 2.7.2 a correction offset was
added that makes the common-case inputs safe, but the guard is still absent and
the OOB write can occur for inputs that produce a bucket index <= -1.

ROOT CAUSE (src/ViennaRNA/subopt/subopt.c, present in every released version)
------------------------------------------------------------------------------
After backtracking each suboptimal structure, the code maps its energy to a
density_of_states histogram bucket:

    # All released versions — only an upper-bound guard, NO lower-bound guard:
    if e > MAXDOS:
        e = MAXDOS
    density_of_states[e] += 1      # OOB write when e < 0

The index 'e' is computed as:

    # 2.6.4 formula (verified from binary / source history):
    e = (int)((structure_energy - min_en) * 10.0)
    # e = -1 when structure_energy < min_en - 0.10 kcal/mol

    # 2.7.2 formula (verified by disassembling the installed .so):
    correction = -0.1f if min_en < 0 else +0.1f
    e = (int)((structure_energy - min_en) * 10.0 - correction)
    # e = -1 when structure_energy < min_en - 0.11 kcal/mol

The correction offset in 2.7.2 makes the threshold slightly harder to reach but
does NOT eliminate the possibility.  Neither version bounds e from below.

HOW THE INDEX GOES NEGATIVE
----------------------------
In the non-recalc code path (dangles=0 or 2, logML=0):
    structure_energy = (double)(state->partial_energy) / 100.0
    min_en           = (double)( float( vrna_eval_structure(fc, mfe_struct) ) )

Because min_en is stored as a float (single precision) before being promoted to
double, it is approximately 1.9e-7 MORE negative than structure_energy for our
specific input (MFE = -12.80 kcal/mol).  This keeps e >= 0 for 2.7.2 with
normal inputs.

In the recalc path (dangles=1, dangles=3, or logML=1):
    structure_energy = (double)( float( vrna_eval_structure(fc, backtracked) ) )
    min_en           = (double)( float( vrna_eval_structure(fc, mfe_struct)  ) )

Here both are floats.  Float arithmetic rounding during multiloop energy
accumulation can cause a suboptimal structure's evaluated energy to fall below
min_en by more than the 0.11 kcal/mol threshold.

CONFIRMED CRASH (ViennaRNA 2.6.4)
----------------------------------
OSTIR's Rust binary hit the OOB write when evaluating RBS strength for the
sequence "ATAAGGAGGTATG" (a 13-nt mRNA with an ATG start codon at position 11).
The exact vrna_subopt_cb call that triggers it:

    sequence : AUAAGGAGGU&ACCUCCUUA
               ^^^^^^^^^^^^^^^^     <- mRNA leader + anti-SD rRNA strand
    delta    : 548  (centi-kcal/mol; = (3.0 + 2.481) * 100, rounded)
    noLP=1, dangles=2, temperature=37.0°C

In 2.6.4's BSS layout, density_of_states[-1] coincides with last_param_file.
Writing +1 to that pointer and then calling vrna_fold_compound() on any sequence
segfaults in strncpy() reading from the corrupted address.

In 2.7.2's BSS layout, density_of_states[-1] coincides with EditCost[3] (a
different global), so the effect of any OOB write would be different but still
undefined behaviour.

SUGGESTED FIX (src/ViennaRNA/subopt/subopt.c)
---------------------------------------------
Add a lower-bound clamp next to the existing upper-bound clamp:

    if (e < 0)
        e = 0;
    else if (e > MAXDOS)
        e = MAXDOS;
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
print(f"Energy delta:  {energy_delta} centi-kcal/mol ({energy_delta/100:.2f} kcal/mol)")
print(f"Model:         noLP=1, dangles=2, temperature=37.0°C")

fc = RNA.fold_compound(seq, md)
(mfe_struct, mfe_en) = fc.mfe()
print(f"MFE:           {mfe_struct}  dG = {mfe_en:.2f} kcal/mol")

# On ViennaRNA 2.6.4 this call writes density_of_states[-1]++ and
# corrupts the adjacent last_param_file global pointer, causing a SIGSEGV
# on the very next call to vrna_fold_compound().
#
# On 2.7.2 a 'correction' offset was added to the bucket-index formula:
#   e = (int)((structure_energy - min_en) * 10.0 - correction)
# (verified by disassembling the installed .so at offset 0x3e32c0).
# This makes the common-case cofold inputs safe, but the code STILL has
# NO lower-bound guard on e — only "cmovg %edx,%eax" (clamp to MAXDOS).
print(f"\nCalling fc.subopt({energy_delta}) ...")
results = fc.subopt(energy_delta)
print(f"  -> returned {len(results)} structures")
for r in sorted(results, key=lambda x: x.energy):
    print(f"     {r.structure}  dG={r.energy:.2f}")

# ---------------------------------------------------------------------------
# Step 2 – the call that segfaults in 2.6.4 (reads the corrupted pointer)
# ---------------------------------------------------------------------------
# Any subsequent vrna_fold_compound call would crash on 2.6.4.
# We use "AU" (the shortest valid RNA) to keep the example minimal.
print(f"\nCalling vrna_fold_compound('AU', md) ...")
print("  On ViennaRNA 2.6.4 this SEGFAULTS inside strncpy()")
print("  because last_param_file (density_of_states[-1] in 2.6.4 BSS) was corrupted.")
fc2 = RNA.fold_compound("AU", md)
(ss2, mfe2) = fc2.mfe()
print(f"  -> structure={ss2!r}  dG={mfe2:.2f}")

print()
print(f"=== Test complete on ViennaRNA {RNA.__version__} ===")
print("The missing lower-bound guard in vrna_subopt_cb has been verified by")
print("binary disassembly (no 'jl'/'cmovl' protecting density_of_states[e]++).")
print()
print("To observe the SIGSEGV directly, install ViennaRNA 2.6.4 and re-run.")
print("See the module docstring for the complete technical analysis.")
