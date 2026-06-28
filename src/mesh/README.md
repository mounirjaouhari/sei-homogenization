# `src/mesh/` — Mesh generation

Structured finite-element meshes for the electrode / SEI / electrolyte trilayer,
in the resolved (SEI explicitly meshed) and effective (SEI → interface) forms.

| File | Dim | Purpose |
|---|---|---|
| `build_trilayer_mesh.m` | 1D | FULL mesh: electrode (60 el.) + SEI (12 el.) + electrolyte (60 el.), 2-node linear elements |
| `build_eff_mesh.m` | 1D | EFFECTIVE mesh: two bulk meshes joined by a doubled interface node at the SEI mid-surface |
| `build_trilayer_mesh_2d.m` | 2D | FULL plane-strain Q4 mesh; the SEI band is explicitly resolved (`nys` elements through-thickness) |
| `build_eff_mesh_2d.m` | 2D | EFFECTIVE plane-strain mesh; SEI replaced by a zero-thickness interface line with **doubled** node pairs (`.iface_bot`, `.iface_top`) |

Convention (2D): `x` is tangential to the interface, `y` is the through-thickness
(normal) direction; the interface normal `n` points from electrode to electrolyte.

Each builder returns a struct with `.nodes`, `.elems`, `.eregion` (region id per
element: 1 = electrode, 2 = SEI, 3 = electrolyte) and boundary node lists
(`.bottom`, `.top`, `.left`, `.right`).
