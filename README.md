# 0 - Dependencies
pip install dash plotly numpy

1 -  Input Data: How to Generate correl.dat ?

The .dat correlation files must be created using CPPTRAJ.

Example CPPTRAJ script:

```bash
parm ../../stripped_all.step3_pbcsetup.parm7
trajin "../../00_prodnoWAT_all.nc"

# Align trajectory to the first frame (recommended)
rms first

# Compute residue-by-residue cross-correlation on C and P atoms
matrix '@C,P' '@C,P' correl out correl.dat byres
```

This produces correl.dat, a square matrix with per-residue correlation values.

The Python app typically loads files named like:

correl_01.dat
correl_02.dat
correl_03.dat

# Editable User Settings

At the top of the script, several variables are designed to be modified by the user.
These control which files are loaded, the regions to plot, thresholds, and display options.

```python
filenames = ["correl_01_9mku.dat", "correl_01_9mku.dat", "correl_01_9mku.dat"]
labels = ["FILE1", "FILE2", "FILE3"]

max_arcs = 2500               # Max arcs drawn per region pair
min_corr_threshold = 0.1      # Only draw arcs above this |value|
index_correction = 3          # Shift residue numbers shown in labels
label_radius = 1.02           # Radius for residue-number labels
radius = 1.0                  # Base circle radius

showIndices = False           # Show residue index numbers on outer ring
showColorBar = True           # Show colorbar on last graph only
```

# Region Definitions

This block defines which residue intervals form each sector of the circular plot:
```python
regions = [
    (953, 1000),  # Region 1: BH-H1
    (1009,1021),  # Region 2: Lid
    (1296, 1303), # Region 3: DNA2
    (1304, 1311)  # Region 4: DNA
]

region_names = ['BH-H1', 'Lid', 'DNA2', 'DNA']
region_colors = ["#e181b0", "#00ff22", "#f01be2", "#f0db1b"]
```

You can add more regions as long as each has a name and color.
