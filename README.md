## Contact Angle Measurement
A Python tool for computing contact angles from molecular dynamics simulation frames of water droplets on surfaces.

### Input 
The script expects: PDB files obtained from MD analysis. For this example it is stored in frames_pure_water folder
- **PDB files**: Molecular frames located in `frames_pure_water/` directory
- **Naming convention**: `frame_<FRAME_NUMBER>.pdb` (e.g., `frame_100000.pdb`)

### Output
Output: The script produces: Average contact angle with standard error and `frames_and_degree.svg'.

### Notes
- Contact angles are computed from four directions of the water-surface interface (x-min, x-max, y-min, y-max)
- The algorithm uses the first 10 z-bins for fitting to focus on the primary interface region
