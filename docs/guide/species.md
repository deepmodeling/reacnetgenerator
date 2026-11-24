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

With the default `--miso 0`, HMM filters each raw atom-and-bond record
independently. Aromatic normalization happens during naming, so a short-lived
Kekulé assignment can be filtered out before it receives a species name.

With `--miso 1`, HMM pools observations of the **same atom indices** when their
bond graphs are isomorphic after ignoring bond orders. The isomorphism must
preserve input atom type IDs, even when two types have the same element name.
With `--miso 2`, the atom indices alone define the pool. Two O₂ molecules on
different atom pairs therefore have independent signals. A single O₂ molecule
that alternates between single and double bond assignments has one continuous
signal while those atoms stay bonded.

HMM filters the union of each pool's observed frames once. On an observed frame
retained by HMM, the raw structure for that frame owns the effective presence.
On a frame added by HMM, the analysis carries forward the preceding observed
variant in the pool, using the first observation for a leading gap. This
assignment supplies a graph for the inferred frame. The molecule CSV continues
to report actual observations. Atom routes treat variants in a pool as one
molecule, including across a change in inferred bond order.

Molecule text/CSV and [timed output](timed-output-schema.md) preserve the original
atom indices and bond records. Separate molecule instances can therefore share a
normalized species name/`species_id` while retaining distinct `molecule_id` values
and different raw bond orders. A switch between equivalent representations does
not produce a species reaction event; the molecule records and their frame ranges
remain available as evidence.
