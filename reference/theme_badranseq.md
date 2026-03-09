# BadranSeq Publication-Ready Theme

**\[stable\]**

Theme function providing consistent SCpubr-like aesthetics for
publication-ready plots. Can be applied to any ggplot2 object.

## Usage

``` r
theme_badranseq(base_size = 14, plot.axes = TRUE)
```

## Arguments

- base_size:

  numeric. Base font size (default: 14).

- plot.axes:

  logical. Whether to show axes (default: TRUE).

## Value

ggplot2 theme object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
ggplot(mtcars, aes(mpg, wt)) +
  geom_point() +
  theme_badranseq()
} # }
```
