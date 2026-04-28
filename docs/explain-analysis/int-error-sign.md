Root Cause: NumPy 2.0 unsigned integer overflow

  The crash chain:

  1. Lib.binsize (cooler reading the hictk-produced mcool) returns self.resol as
   a numpy unsigned integer (likely np.uint32 or np.uint64), because hictk
  writes the bin-size HDF5 attribute as an unsigned type.
  2. NumPy 2.0 NEP 50 (required for Python 3.13) changed type promotion rules.
  Python scalars are now "weak" — they adapt to the numpy type instead of
  promoting it. So:
  # Line 314: test_region_start0 = test_region_start1 - (400 * self.resol)

  # Old NumPy 1.x behavior:
  1 - np.uint32(2000000)  →  np.int64(-1999999)   # promoted to signed

  # New NumPy 2.0 behavior:
  1 - np.uint32(2000000)  →  np.uint32(4292967297) # WRAPS! unsigned can't go 
  negative
  3. The wrapped value (~4.3 billion) is a huge positive number, so the safety
  check if test_region_start0 <= 1 at line 321 doesn't trigger (it's meant to
  catch negative values).
  4. The corrupted coordinate creates a region string like
  "chr1:4292967297-199000000" where start > end → ValueError: End coordinate 
  less than start.

  Fix

  Patch self.resol and chrsize to Python int in the installed stripenn, which
  has arbitrary precision and no unsigned issues. Run on HPC:

  conda activate mariner_env
  SITE_PKG=$(python -c "import stripenn; import os; 
  print(os.path.dirname(stripenn.__file__))")
  echo "Patching: ${SITE_PKG}/getStripe.py"

  # Fix 1: self.resol = int(resol) in __init__
  sed -i 's/        self.resol = resol/        self.resol = int(resol)/'
  "${SITE_PKG}/getStripe.py"

  # Fix 2: chrsize = int(...) in main_null_calc to prevent similar unsigned 
  issues
  sed -i 's/chrsize = self.chromnames2sizes\[chr\]/chrsize = 
  int(self.chromnames2sizes[chr])/' "${SITE_PKG}/getStripe.py"

  # Verify
  grep 'self.resol = ' "${SITE_PKG}/getStripe.py"
  grep 'chrsize = ' "${SITE_PKG}/getStripe.py"

  The int() cast converts numpy unsigned scalars to Python's arbitrary-precision
   integers. All downstream arithmetic then uses Python ints — no overflow
  possible. No other code changes needed; the fix propagates through every
  self.resol and chrsize usage.

rsync -avP expanse:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/logs/ /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/logs/

```
rsync -avhP expanse:/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/250402/{ctrl_merged,mut_merged}.hic \
  ~/sdsc/hic/250402/
```

rsync -avz --progress \
  expanse:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/outputs/250402/visualizations/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/outputs/250402/visualizations/ \
&& rsync -avz --progress \
  expanse:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/outputs/250831/visualizations/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/outputs/250831/visualizations/ \
&& rsync -avz --progress \
  expanse:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/outputs/combined/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/outputs/combined/


```
rsync -avzP \
  expanse:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/outputs/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/outputs/
```

```
rsync -avz --progress \
  expanse:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/outputs/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/outputs/
```

```c

 29 int* parse_rotor_indices(char* rotor_ind_str, int num_rotors) {
 30     int* indices = malloc(num_rotors * sizeof(int));
 31     int index = 0;
 32     for (int i = 0; rotor_ind_str[i] != '\0'; ++i) {
 33         if (rotor_ind_str[i] >= '0' && rotor_ind_str[i] <= '9') {
 34             indices[index] = rotor_ind_str[i] - '0';
 35             index++;
 36         }
 37     }
```

git reset --hard origin/main
