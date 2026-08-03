# Raster generation

This workspace produces the site's hero raster as a reproducible hierarchy of spatial layers.

## Planned sequence

1. Broad, isotropic background field
2. Dominant central basin
3. Irregular basin geometry
4. Secondary regional minima and maxima
5. Medium-scale heterogeneity
6. Restrained fine-scale variation
7. Stabilized text-safe interior
8. Connectivity and value-distribution tuning
9. Truncated viridis color mapping
10. Final compositional adjustments

The landscape composition will be developed first. A coordinated portrait composition will follow after the landscape field and visual language are established.

## Output convention

Intermediate images should be numbered by stage, for example:

```text
01_broad-field.png
02_central-basin.png
03_irregular-basin.png
```

Generated iterations belong in `outputs/landscape/` or `outputs/portrait/`. These files are ignored by Git. Only selected, web-optimized final images should be copied to the repository-level `assets/images/` directory.
