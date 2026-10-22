# Reticulata Phase 2 — raw results and emitted certificates

Everything here is produced by `rust/reticulation`; nothing is hand-edited. The
interpretation is in [`docs/GENETICA_RETICULATA.md`](../../docs/GENETICA_RETICULATA.md)
§10. Regenerate with:

```bash
cd rust/reticulation
cargo build --release
./target/release/reticulation experiment   > .../experiment.json
./target/release/reticulation saturation   > .../saturation.json
./target/release/reticulation verify-tree --out <dir> > .../verify_tree.json
```

All three are deterministic given their seeds, so a reader who reruns them should
get these files back byte for byte.

| file | what it is |
| :--- | :--- |
| `experiment_n32.json` | The same grid at 32 leaves and three seeds (`--generations 5 --seeds 3 --max-quadruples 60000`), 35 960–46 376 quadruples per dataset searched exhaustively. The scale check where the triple statistic comes out ordered the wrong way round. |
| `experiment.json` | Seven simulated populations × five seeds. Pre-registered predictions with their criteria, then the measurements. Includes the block recording the two criteria that a first run showed were wrong, and why. |
| `saturation.json` | The length-scaling diagnostic: is the clonal four-point margin noise or saturation bias? Plus a divergence sweep giving the noise floor from mean p-distance 0.006 to 0.33. |
| `verify_tree.json` | Kernel verdicts for the two new certificate forms, in both directions. |
| `emitted/*.lean` | **The evidence itself.** Four generated Lean files, each embedding the sequences it is about. |

## Re-checking the certificates without rerunning the pipeline

This is the point of the architecture, so it should be easy. The two genuine
files must typecheck and the two corrupted ones must not:

```bash
cd <repo with formal/ANSE built>
lake env lean data/reticulata_phase2/emitted/reticulata_quartet.lean          # accepted, ~71 s
lake env lean data/reticulata_phase2/emitted/reticulata_quartet_bad.lean      # must FAIL
lake env lean data/reticulata_phase2/emitted/reticulata_complete_triple.lean  # accepted, ~42 s
lake env lean data/reticulata_phase2/emitted/reticulata_complete_triple_bad.lean  # must FAIL
```

The kernel computes every count from the embedded sequences by `rfl`; the files
hand it nothing but the sequences and the claimed counts. What remains trusted is
only that the embedded lists are the simulated population's sequences, and the
encoding is printed in each file's header (`A=0 C=1 G=2 T=3`, anything else `4`
and not comparable) so that can be checked too.

A file that fails where it should pass is a bug report about the Rust, not a
negative result.
