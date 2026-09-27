"""Named, source-tagged nonbonded parameter sets for the ITO slab.

The slab is held rigid and its internal interactions are excluded, so only
SAM-slab cross terms matter.  Those come from OPLS geometric mixing
(``pair_modify mix geometric``) of the slab self terms below
with the LigParGen ligand self terms, exactly as ligand-ligand terms are
formed in the NiO workflow.  None of these sets is fitted to ITO adsorption
data; the pilot compares them as a sensitivity bracket.
"""
from __future__ import annotations

TWO_SIXTH = 2.0 ** (1.0 / 6.0)

# epsilon (kcal/mol), sigma (A)
PARAMETER_SETS: dict[str, dict] = {
    "uff-cation/clayff-anion": {
        "In": (0.599, 4.463 / TWO_SIXTH),
        "Sn": (0.567, 4.392 / TWO_SIXTH),
        "O": (0.1554, 3.1655),
        "Oh": (0.1554, 3.1655),
        "Hh": (0.0, 0.0),
        "sources": {
            "In,Sn": "UFF, Rappe et al., JACS 114, 10024 (1992); x1/D1 as tabulated in Open Babel data/UFF.prm (In3+3, Sn3); sigma = x1/2^(1/6)",
            "O,Oh,Hh": "CLAYFF ob/oh/ho, Cygan et al., J. Phys. Chem. B 108, 1255 (2004)",
        },
        "note": "Upper bracket: UFF metal well depths are large (0.6 kcal/mol).",
    },
    "clayff-cation/clayff-anion": {
        # CLAYFF octahedral Al (ao) as a bare-cation analogue for In/Sn.
        "In": (1.3298e-6, 4.2713),
        "Sn": (1.3298e-6, 4.2713),
        "O": (0.1554, 3.1655),
        "Oh": (0.1554, 3.1655),
        "Hh": (0.0, 0.0),
        "sources": {
            "In,Sn": "CLAYFF ao (octahedral Al) used as an electrostatics-dominated cation analogue; Cygan et al. 2004. NOT an In parameterization.",
            "O,Oh,Hh": "CLAYFF ob/oh/ho, Cygan et al. 2004",
        },
        "note": "Lower bracket: cation dispersion effectively off; adhesion is Coulomb + O dispersion.",
    },
}


def surface_pair_lines(type_ids: dict[str, int], parameter_set: str) -> list[str]:
    params = PARAMETER_SETS[parameter_set]
    lines = []
    for label, tid in sorted(type_ids.items(), key=lambda kv: kv[1]):
        eps, sigma = params[label]
        # sigma=0 is rejected by some LAMMPS mixing paths; any sigma with eps=0 is inert.
        sigma = sigma if sigma > 0 else 1.0
        lines.append(f"pair_coeff {tid} {tid} {eps:.7g} {sigma:.6f} # slab {label} [{parameter_set}]")
    return lines


HEADER = """# ITO/SAM extension force-field styles (rigid slab; ligand terms as in the NiO workflow)
bond_style      harmonic
angle_style     harmonic
dihedral_style  hybrid opls charmm
improper_style  cvff
pair_style      lj/cut/coul/long 10.0 8.0
pair_modify     mix geometric
kspace_style    pppm 1e-4
# ad differentiation: ~2x cheaper than ik for this rigid-slab geometry (measured)
kspace_modify   slab 3.0 diff ad
special_bonds   amber
"""
