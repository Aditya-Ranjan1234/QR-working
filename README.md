# QR-working

An educational, interactive QR-code and error-correction visualizer. The app lets you explore GF(256) arithmetic, polynomial operations, Reed–Solomon encoding/decoding, QR generation, damaged-codeword simulation, and recovery traces.

## Features

- GF(256) multiplication, polynomial representations, and primitive-polynomial table construction
- Reed–Solomon message encoding, error injection, Berlekamp–Massey tracing, and correction
- QR Code generation with version, error-correction, mask, data, and ECC details
- Scratch-and-recover simulation for experimenting with damaged QR matrices
- Fully client-side browser engine for free GitHub Pages hosting
- Optional FastAPI backend for local API experimentation

## Run locally

The GitHub Pages version needs no build step. Open `index.html` in a browser, or serve the repository with a local static server:

```powershell
python -m http.server 8000
```

For the optional FastAPI API, install the dependencies in a virtual environment and run:

```powershell
uvicorn app:app --reload
```

## Test

```powershell
python tests/test_qr.py
```

## Project layout

```text
QR-working/
├── index.html, index.css, app.js   # GitHub Pages frontend
├── algorithms/                     # GF(256), QR, and Reed–Solomon implementation
├── app.py                          # Optional FastAPI API
├── tests/                          # Python verification tests
└── docs/                           # Architecture and development notes
```

## Deploy to GitHub Pages

The repository root is already the publishable site. In GitHub, open **Settings → Pages**, choose **Deploy from a branch**, select `main` and `/ (root)`, and save. No paid hosting or server is required.

## License

This project is intended for learning and experimentation.
