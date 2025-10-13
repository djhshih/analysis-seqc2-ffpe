#!/usr/bin/env python
import os
import glob

# %%
bam_paths = glob.glob("*.bam*")

for path in bam_paths:
    if "HCC1395_" in path:
        name = path.replace("HCC1395_", "")
        tokens = name.split(".")
        new_name = f"AMS_AB_T_{tokens[2]}.{tokens[1]}.{".".join(tokens[3:])}"
        os.rename(path, new_name)
        print(f"Renamed {path} to {new_name}")
        
    if "HCC1395BL_" in path:
        name = path.replace("HCC1395BL_", "")
        tokens = name.split(".")
        new_name = f"AMS_AB_N_{tokens[2]}.{tokens[1]}.{".".join(tokens[3:])}"
        os.rename(path, new_name)
        print(f"Renamed {path} to {new_name}")


