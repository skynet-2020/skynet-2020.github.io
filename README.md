# Personal website

A minimal, one-page GitHub Pages site built around a custom spatial raster hero image.

## Repository structure

```text
.
├── index.html                    # Production webpage (added during site build)
├── styles.css                   # Production styles (added during site build)
├── assets/
│   ├── images/                  # Approved, web-optimized hero images
│   └── documents/               # CV and other downloadable documents
└── raster-generation/
    ├── README.md                # Raster workflow and conventions
    ├── scripts/                 # Reproducible R scripts
    └── outputs/
        ├── landscape/           # Landscape iterations and diagnostics
        └── portrait/            # Mobile iterations and diagnostics
```

Only final, web-optimized assets are copied into `assets/`. Experimental raster outputs remain local and are ignored by Git.
