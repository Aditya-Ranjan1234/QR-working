# Architecture

The project has two independent execution paths:

## Browser application

`index.html`, `index.css`, and `app.js` form a static, dependency-free frontend. The JavaScript contains the GF(256), Reed–Solomon, and QR algorithms, so all interactive features work without a server. This is the version deployed by GitHub Pages.

## Optional Python API

`app.py` exposes the same educational operations as FastAPI endpoints under `/api`. It imports the reference implementations from `algorithms/` and serves the repository root for local development.

## Algorithm modules

- `algorithms/gf256.py`: finite-field arithmetic and lookup tables
- `algorithms/polynomial.py`: polynomial operations over GF(256)
- `algorithms/reedsolomon.py`: encoding, syndromes, error location, and correction
- `algorithms/qr_encoder.py`: QR matrix construction and masking
- `algorithms/qr_decoder.py`: matrix extraction, damage detection, and recovery
