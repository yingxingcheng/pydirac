## [2024-6-3]

### Changed

- Remove `pydirac.analysis.constants` to `pydirac.analysis.utility`
- `get_svd_from_array` in `pydirac.analysis.polarizability.PolarizabilityCalculator` to accept an additional `rcond` parameter
- Updated `do_one_basis` function to include errors in the polarizability results
- Improved handling of residuals and calculation of covariance matrix in `get_svd_from_array`
- Updated `get_polarizability_from_output_list` and `do_one_basis` to better handle the structure and calculation of results
- Refactored test cases in `pydirac.tests.io.test_outputs`

### Fixed

- Corrected `np.alltrue` to `np.all` for checking if all roots are converged
- Fixed test assertions in `pydirac.tests.io.test_outputs` to use `pytest.approx` correctly

## [2023-4-5]

### Changed

- Remove `pydirac.analysis.constants` to `pydirac.analysis.utility`

### Added

- Test data of He atom

## [2023-4-4]

### Fixed

- Minor typos

### Changed

- move `setup.py` to `pyproject.toml`

### Added

- Documents
