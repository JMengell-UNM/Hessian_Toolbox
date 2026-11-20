## Write a simplified Hessian file from the parsed BDT components
from copy import deepcopy
import numpy as np




def write_hess(self, out_path, opp_dict=None):
    """Write a simplified .hess file (like Minimal_H2_vib.hess) from self.hess_dict.
    opp_dict: optional {key: new_value} overrides; sizes must match the original.
    """
    
    def _write_matrix(f, M, block=5):
        """ORCA-style matrix: blocks of 5 columns, each preceded by a header line of column indices."""
        n_rows, n_cols = M.shape
        for c0 in range(0, n_cols, block):
            c1 = min(c0 + block, n_cols)
            f.write("".join(f"{c:>19d}" for c in range(c0, c1)) + "\n")
            for r in range(n_rows):
                f.write(f"{r:>5d}" + "".join(f"{M[r, c]:>19.10E}" for c in range(c0, c1)) + "\n")
    
    d = deepcopy(self.hess_dict)
    if opp_dict:
        unknown = set(opp_dict) - set(d)
        if unknown:
            raise KeyError(f"Keys not in hess_dict: {unknown}")
        d.update(deepcopy(opp_dict))

    symbols = d["atomic_symbols"]
    n = len(symbols)
    n3 = 3 * n

    modes  = np.asarray(d["normal_modes"], float)
    freqs  = np.asarray(d["vibrational_frequencies"], float)
    coords = np.asarray(d["atomic_coords"], float).reshape(-1)
    masses = np.asarray(d["molecular_weights"], float)

    for name, got, want in [("normal_modes", modes.shape, (n3, n3)),
                            ("vibrational_frequencies", freqs.shape, (n3,)),
                            ("atomic_coords", coords.shape, (n3,)),
                            ("molecular_weights", masses.shape, (n,))]:
        if got != want:
            raise ValueError(f"{name}: shape {got}, expected {want}")

    with open(out_path, "w") as f:
        f.write("$orca_hessian_file\n\n")

        f.write(f"$hessian\n{n3}\n")
        _write_matrix(f, np.zeros((n3, n3)))      # zeros, like the minimal H2 file

        f.write(f"\n$vibrational_frequencies\n{n3}\n")
        for i, v in enumerate(freqs):
            f.write(f"{i:>5d}{v:>26.16f}\n")

        f.write(f"\n$normal_modes\n{n3} {n3}\n")
        _write_matrix(f, modes)

        f.write(f"\n$atoms\n{n}\n")
        for sym, m, (x, y, z) in zip(symbols, masses, coords.reshape(-1, 3)):
            f.write(f" {sym:<2s}{m:>12.5f}{x:>20.12f}{y:>20.12f}{z:>20.12f}\n")

        f.write(f"\n$frequency_scale_factor\n  {d.get('frequency_scale_factor', 1.0):.6f}\n")
        f.write("\n$end\n")                        # blank line before $end is needed by the parser


Hessian.write_hess = write_hess                    # attach to the existing class


# --- Parse H2 and write the simplified file ---
H2 = Hessian(r"C:\Users\apelr\OneDrive - University of New Mexico\Python Scripts\Josh_Scripts\H2_vib.hess")
out_path = r"C:\Users\apelr\OneDrive - University of New Mexico\Python Scripts\Josh_Scripts\H2_simplified_vib.hess"
H2.write_hess(out_path)

# --- Round-trip check: re-parse and compare every extracted component ---
H2_check = Hessian(out_path)
for key, a in H2.hess_dict.items():
    b = H2_check.hess_dict[key]
    if key == "atomic_symbols":
        same = (a == b)
    else:
        same = np.allclose(a, b, rtol=1e-6, atol=1e-9)
    print(f"{key:28s} {'OK' if same else 'MISMATCH'}")