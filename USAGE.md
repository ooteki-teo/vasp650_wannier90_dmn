# Usage

## Generate and apply the patch

Use a clean, licensed VASP 6.5.0 source tree, Git, Ruby 2.6 or newer, and a `diff` implementation with `-u` and `-L` support. The generator checks the five affected source files against their VASP 6.5.0 fingerprints and writes the unified diff to standard output.

```bash
git clone https://github.com/ooteki-teo/vasp650_wannier90_dmn.git
patch_repo="$(cd vasp650_wannier90_dmn && pwd)"
vasp_src=/absolute/path/to/vasp.6.5.0
dmn_patch=/absolute/path/to/vasp.6.5.0-wannier90-dmn.patch

ruby "$patch_repo/make-patch.rb" "$vasp_src" > "$dmn_patch"
cd "$vasp_src"
git apply --check "$dmn_patch"
git apply "$dmn_patch"
```

A fingerprint mismatch indicates an already-patched file, a different release, or another source modification. Generate the diff against a clean 6.5.0 copy first. Keep the generated diff with your licensed source tree.

## Build

Enable the Wannier90 interface and link a compatible Wannier90 library in your existing `makefile.include`. Wannier90 3.1.0 was used for validation.

```makefile
CPP_OPTIONS += -DVASP2WANNIER90
LLIBS += /absolute/path/to/libwannier.a
```

```bash
make DEPS=1 -j 8 std
make DEPS=1 -j 8 ncl
```

Use `vasp_std` for scalar or collinear calculations and `vasp_ncl` for SOC spinors. Use CPU complex-plane-wave builds. The existing PEAD interface requires `NCORE = 1`.

## Export and run standalone Wannier90

Set the following tags together with the usual electronic-structure settings for your system.

```text
ISYM = 2
NCORE = 1
LWANNIER90 = .TRUE.
LWRITE_MMN_AMN = .TRUE.
LWRITE_DMN = .TRUE.
LWANNIER90_RUN = .FALSE.
```

Supply `NUM_WANN` and explicit projections in `WANNIER90_WIN` or the supported `wannier90.win` input. Include `site_symmetry = .true.` in the Wannier90 settings. The [Si INCAR](tests/Si/INCAR) provides a complete four-band example with projections at the bond centers. For an SOC example, see [TaAs](tests/TaAs-final/INCAR).

Choose a uniform k-point mesh closed under the VASP space-group operations. The `NUM_WANN` projection orbitals must likewise form a closed set, including their site and local axis choices. For disentanglement, select an outer window containing complete symmetry multiplets and the required target representations. Wannier90 3.1.0 supports the outer window in site-symmetry mode. `exclude_bands` is applied consistently to the exported DMN, AMN, MMN and EIG data.

After VASP finishes, keep the exported `.dmn`, `.amn`, `.mmn`, `.eig` and `.win` files in the working directory. For scalar and spinor runs, use

```bash
wannier90.x wannier90
```

For a collinear `ISPIN = 2` run, VASP exports separate matrices for the `wannier90.1` and `wannier90.2` seeds.

```bash
for spin in 1 2; do
    cp wannier90.win "wannier90.$spin.win"
    wannier90.x "wannier90.$spin"
done
```

`OUTCAR` lists the exported operations and their fractional translations. Its maximum band norm defect includes all exported bands. If the highest band cuts a degenerate multiplet, assess the norm defect within the actual disentanglement window. DMN stores the unitary space-group operations. The TaAs validation separately checks time reversal in the interpolated spectrum and Weyl-node orbits.

## Reproduce the regression cases

The test runner needs MPI, the patched VASP executables, standalone Wannier90, a VASP PAW-PBE library, and Python with NumPy, SciPy and spglib. Point the runner at an output directory that does not yet exist.

```bash
export VASP_STD=/absolute/path/to/vasp_std
export VASP_NCL=/absolute/path/to/vasp_ncl
export WANNIER90=/absolute/path/to/wannier90.x
export POTCAR_ROOT=/absolute/path/to/vasp-paw-pbe
export DMN_PYTHON=/absolute/path/to/python
bash "$patch_repo/tests/run-tests.sh" /absolute/path/to/new-output-directory
```

For the optional Si comparison, also set `DMN_QE_BIN` to a QE 7.5 executable directory and `ESPRESSO_PSEUDO` to an SSSP PBE directory containing `Si.us.pbe.z_4.ld1.psl.v1.0.0-high.upf`. The TaAs test uses coordinates from the [WannierTools TaAs example](https://github.com/quanshengwu/wannier_tools/tree/master/examples/TaAs) and searches near the W1 and W2 locations reported by [Garcia et al.](https://doi.org/10.1103/PhysRevResearch.2.013073)
