# Development notes

## How to make releases

- Mark release in `CHANGELOG.md`
- Modify version in `pyproject.toml`
- Make a new commit and tag it with `vYear.Month.Day`
- Remove lingering `dist`
- Run `python3 -m build`
- Run `python3 -m twine upload dist/*`
