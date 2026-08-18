# FBKdock Example Data

A small, self-contained test set so you can verify that the pipeline works
end-to-end on your machine before running your own screening.

## Contents

```
examples/
├── protein/
│   └── 2Y9X.pdb              # Receptor: PPO3 tyrosinase + tropolone inhibitor
├── ligands/
│   ├── 3a.sdf ... 3e.sdf     # 5 ligand candidates (3D, with hydrogens)
└── config/
    ├── config_vinagpu.txt    # Grid template tuned for the 2Y9X binding site
    └── config_autodockvina.txt
```

- **2Y9X** is the crystal structure of PPO3, a tyrosinase from *Agaricus
  bisporus*, solved in complex with the inhibitor **tropolone (residue name
  `0TR`)**. It contains 8 chains (A–D: enzyme, E–H: lectin-like subunits),
  copper ions, and water molecules — a realistic PDB with the exact features
  the pipeline's interactive cleaning tools are designed for.
- **3a–3e.sdf** are five small organic molecules used as screening ligands.
- The **config templates** already contain the grid box computed around the
  tropolone binding site of chain A (`center_x/y/z = -9.923 / -26.885 /
  -43.059`, `size = 15 x 23 x 23`).

## How to use

1. Copy the receptor into the runtime directory:
   ```
   copy examples\protein\2Y9X.pdb Protein\        (Windows)
   cp examples/protein/2Y9X.pdb Protein/           (Linux)
   ```
   `Protein/` must contain **exactly one** PDB file.

2. Copy the ligands:
   ```
   copy examples\ligands\*.sdf Ligand\            (Windows)
   cp examples/ligands/*.sdf Ligand/               (Linux)
   ```

3. Launch the pipeline and follow the prompts (see the main README for the
   suggested answers for 2Y9X):
   ```
   conda activate FBKdock_env
   python main.py
   ```

4. The grid in `config_vinagpu.txt` / `config_autodockvina.txt` is already
   set for 2Y9X. **For your own protein** you must re-center the box on your
   binding site and update these files yourself.

## Re-centering the grid for your own protein

The docking box is defined by six numbers in the config template:

| Key | Meaning |
|---|---|
| `center_x/y/z` | Center of the search box (Angstrom) |
| `size_x/y/z` | Box size in each dimension (Angstrom) |

A practical workflow:

1. Open your receptor–ligand complex in PyMOL.
2. Select the binding-site ligand: `select lig, resn <LIG>`.
3. Get the centroid: `get_center lig` (PyMOL ≥ 2.4) or
   `print cmd.centerofmass("lig")`.
4. Use the printed coordinates as `center_x/y/z` and choose a box size that
   encloses the binding pocket (20–25 A per side is a typical starting
   point).
5. Copy the values into both config templates at the project root.

The box must fully enclose the co-crystallized ligand; otherwise the
redocking validation (Step 2) will fail on purpose — that is the
pipeline telling you the grid is wrong.
