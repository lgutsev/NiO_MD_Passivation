from __future__ import annotations
import argparse
from pathlib import Path
from .build import build
from .config import molecule_dir, molecule_manifest, missing_ligpargen
from .geometry import elements
from .lammps import charge, parse
from .validate import validate

TEMPLATE='''[molecule]\ndisplay_name = "{name}"\nslug = "{slug}"\nrole = "secondary"\nexpected_net_charge = 0.0\nphosphonic_acid_anchors = 1\nanchor_atom_ids = []\n\n[files]\nligpargen = "ligpargen.lmp"\ngeometry = ""\noverride = ""\n'''

def main(argv=None)->int:
    p=argparse.ArgumentParser(prog="nio-md-prep"); sub=p.add_subparsers(dest="command",required=True)
    i=sub.add_parser("init-molecule"); i.add_argument("name")
    i=sub.add_parser("inspect-molecule"); i.add_argument("name")
    i=sub.add_parser("build"); i.add_argument("config",type=Path); i.add_argument("--output",type=Path,required=True); i.add_argument("--packed-xyz",type=Path)
    i=sub.add_parser("validate"); i.add_argument("output",type=Path)
    i=sub.add_parser("prepare-sequential-stage2"); i.add_argument("config",type=Path); i.add_argument("--primary-final",type=Path,required=True); i.add_argument("--output",type=Path,required=True); i.add_argument("--packed-xyz",type=Path)
    a=p.parse_args(argv)
    try:
        if a.command=="init-molecule":
            slug=a.name.lower().replace("_","-").replace(" ","-"); folder=molecule_dir(slug); folder.mkdir(parents=True,exist_ok=True)
            manifest=folder/"molecule.toml"
            if not manifest.exists(): manifest.write_text(TEMPLATE.format(name=a.name,slug=slug),encoding="utf-8")
            print(f"Created {manifest}\nPlace the independently generated LigParGen LAMMPS output at:\n{(folder/'ligpargen.lmp').resolve()}")
        elif a.command=="inspect-molecule":
            folder,m=molecule_manifest(a.name); path=folder/m["files"].get("ligpargen","ligpargen.lmp")
            if not path.exists(): raise missing_ligpargen(path,f"nio-md-prep inspect-molecule {a.name}")
            d=parse(path); symbols=elements(d); print(f"{m['molecule']['display_name']}: {d.count('Atoms')} atoms, charge {charge(d):.8f}, elements {' '.join(symbols)}")
        elif a.command=="build": build(a.config,a.output,packed_xyz=a.packed_xyz); print(f"Build plan written to {a.output}")
        elif a.command=="validate": validate(a.output); print((a.output/"validation_report.txt").read_text(),end="")
        else:
            target=a.output/"stage2_secondary"
            build(a.config,target,primary_final=a.primary_final,packed_xyz=a.packed_xyz)
            print(f"Sequential stage 2 written to {target}")
    except (ValueError,FileNotFoundError,RuntimeError) as e:
        p.exit(2,f"error: {e}\n")
    return 0
