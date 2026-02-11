# Build a Single Violin Panel

Internal helper that constructs a single violin + boxplot + jitter +
centrality panel for one gene, using ggbetweenstats-style aesthetics.

## Usage

``` r
.build_violin_panel(
  df,
  y_col,
  group_col,
  title,
  colors,
  pt.size = 3,
  pt.alpha = 0.4,
  violin.alpha = 0.2,
  violin.width = 0.5,
  boxplot = TRUE,
  boxplot.width = 0.3,
  boxplot.alpha = 0.2,
  centrality = TRUE,
  centrality.fn = "median"
)
```

## Arguments

- df:

  data.frame. Long-format data with at least the y_col and group_col.

- y_col:

  character. Column name for the y-axis (expression values).

- group_col:

  character. Column name for the grouping variable.

- title:

  character. Plot title (typically the gene name).

- colors:

  character. Named vector of colors for groups.

- pt.size:

  numeric. Jitter point size.

- pt.alpha:

  numeric. Jitter point alpha.

- violin.alpha:

  numeric. Violin fill alpha.

- violin.width:

  numeric. Violin width.

- boxplot:

  logical. Whether to overlay a boxplot.

- boxplot.width:

  numeric. Boxplot width.

- boxplot.alpha:

  numeric. Boxplot alpha.

- centrality:

  logical. Whether to show centrality point and label.

- centrality.fn:

  character. Centrality function name ("median" or "mean").

## Value

A ggplot2 object.
