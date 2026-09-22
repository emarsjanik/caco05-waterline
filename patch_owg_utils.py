#!/usr/bin/env python3
"""
Patch the upstream OWG utils.py to accept any training CSV.

Upstream branches on the exact CSV filename in two places. Any new
name falls through every branch, so `df['path']` is never created and
`weights_path` is never assigned, and training crashes. This adds a
generic fallback to both, keyed on the CSV's own stem, without
altering behaviour for the original three datasets.

A backup is written to utils.py.orig before anything changes.
"""
import re, sys, shutil
from pathlib import Path

p = Path(sys.argv[1] if len(sys.argv) > 1 else "utils.py")
src = p.read_text()
if "GENERIC_CSV_FALLBACK" in src:
    print("already patched"); sys.exit(0)
shutil.copy2(p, str(p) + ".orig")

# 1. get_and_tidy_df: generic path construction
src = src.replace(
'''        df = df.rename(index=str, columns={" H": "H", " T": "T"})''',
'''        # GENERIC_CSV_FALLBACK: any other CSV names images by id + ".jpg"
        if 'path' not in df.columns:
                df['path'] = df['id'].map(lambda x: os.path.join(base_dir,
                                                                image_dir,'{}'.format(x)))+".jpg"

        df = df.rename(index=str, columns={" H": "H", " T": "T"})''')

# 2. get_weights_path: generic fallback on the CSV stem
src = src.replace(
'''        return weights_path''',
'''        # GENERIC_CSV_FALLBACK: derive a weights path from the CSV stem
        if 'weights_path' not in locals():
                stem = os.path.splitext(os.path.basename(input_csv_file))[0]
                kind = 'waveheight' if category == 'H' else 'waveperiod'
                weights_path = os.path.join(os.getcwd(), 'im'+str(imsize), 'res',
                        str(num_epochs)+'epoch', category, 'model'+str(counter),
                        'batch'+str(batch_size),
                        kind+'_weights_model'+str(counter)+'_'+str(batch_size)+'batch.best.'+stem+'.hdf5')
        return weights_path''')

p.write_text(src)
n = src.count("GENERIC_CSV_FALLBACK")
print(f"patched {p} ({n} fallback(s) added); original saved as {p}.orig")
