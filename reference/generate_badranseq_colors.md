# Generate Categorical Color Palette

Internal function to generate publication-ready categorical colors
matching SCpubr's default palette. Uses colorspace "Dark 3" palette with
adjusted saturation and value for more vivid, slightly darker colors.

## Usage

``` r
generate_badranseq_colors(n, custom.colors = NULL)
```

## Arguments

- n:

  numeric. Number of colors needed.

- custom.colors:

  character. Optional custom color vector.

## Value

Character vector of hex colors.
