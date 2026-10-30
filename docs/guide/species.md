# Chemical species identification

By default (`--miso 0`), equivalent Kekulé representations of an aromatic
molecule have the same species name. For example, shifting the alternating
single and double bonds in 1-chloro-2-fluorobenzene does not create a second
species or a reaction between the two representations. Both SMILES and VF2
identification use the same aromatic bond normalization.

## Normalization policy

ReacNetGenerator uses RDKit's default aromaticity model on a temporary copy of
the molecular graph. This requires successful valence checking, ring perception,
kekulization, radical assignment and aromaticity perception. All hydrogens remain
explicit; missing hydrogens or formal charges are not inferred from the input.

Only the resulting aromatic atom flags and aromatic bond orders are used for
species identification. Other bond orders, atom identities and connectivity are
preserved. A bond can carry an aromatic flag while retaining a triple-bond type;
that type is preserved too. This does not merge tautomers, neutralize charges,
repair nitro groups or convert metal coordination bonds to dative bonds.

The supported elements are H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, Te and I,
with single, double, triple or aromatic bonds. Acyclic structures, molecules
outside this domain, and structures that fail valence or kekulization checks
retain their original bond graph and naming behavior. In particular, a transient
overvalent intermediate or a molecule missing charge information is retained;
normalization failure does not discard its frames. If SMILES serialization itself
fails, the existing VF2/`_unknownSMILES` fallback still applies.

The equivalence is limited to aromatic bonds recognized by this RDKit model;
it is not a general resonance equivalence algorithm. Literal canonical names
can also differ across RDKit versions or from earlier ReacNetGenerator releases.

## Filtering, merging and trajectory evidence

Normalization happens during species identification, **after HMM filtering**.
Detection and HMM still operate on the original atom-and-bond records. Signals
for different Kekulé assignments are not pooled before filtering, so short-lived
assignments can still be filtered independently. Normalization does not recover
records already removed by HMM. Use the [HMM guidance](hmm.md) when choosing
filter parameters.

`--miso 1` additionally ignores all bond orders, while `--miso 2` additionally
ignores connectivity. These broader merging modes retain their existing meaning;
they are not needed to merge the supported aromatic representations.

Molecule text/CSV and [timed output](timed-output-schema.md) preserve the original
atom indices and bond records. Separate molecule instances can therefore share a
normalized species name/`species_id` while retaining distinct `molecule_id` values
and different raw bond orders. A switch between equivalent representations does
not produce a species reaction event; the molecule records and their frame ranges
remain available as evidence.
