#!/usr/bin/env python3
# %autoindent

from pathlib import Path

import xarray as xr

raw_output_dir = Path("~/projects/sbg-uncertainty/isofit/examples/py-hypertrace/output/ht-trait-analysis/").expanduser()

ncfiles = list(raw_output_dir.glob("*.nc"))
ncfiles.sort()

ncf = ncfiles[0]

dat = xr.open_dataset(ncf)
