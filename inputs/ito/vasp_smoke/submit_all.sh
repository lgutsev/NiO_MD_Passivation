#!/bin/bash
# Submit every ITO VASP reference case. Review docs/ito/vasp-smoke.md first.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
(cd "$here/slab-bare" && sbatch run.sbatch)
(cd "$here/slab-oh" && sbatch run.sbatch)
(cd "$here/slab-oh-sn" && sbatch run.sbatch)
(cd "$here/mol-me-4pacz" && sbatch run.sbatch)
(cd "$here/ads-me-4pacz-phys" && sbatch run.sbatch)
(cd "$here/ads-me-4pacz-bidentate" && sbatch run.sbatch)
(cd "$here/mol-meo-2pacz" && sbatch run.sbatch)
(cd "$here/ads-meo-2pacz-phys" && sbatch run.sbatch)
(cd "$here/ads-meo-2pacz-bidentate" && sbatch run.sbatch)
