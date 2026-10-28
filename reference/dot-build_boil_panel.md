# Build a Single Boil Plot Panel

Internal helper that constructs a single boil plot panel for one gene: a
violin + boxplot + jitter + centrality layer for non-zero cells, and a
detection bar below zero showing the fraction of expressing cells.

## Usage

``` r
.build_boil_panel(
  df,
  group_col,
  title,
  x.title,
  colors,
  min.nonzero = 10,
  pt.size = 0.7,
  pt.alpha = 0.6,
  jitter.width = 0.1,
  violin.alpha = 0.8,
  violin.width = 0.9,
  boxplot = FALSE,
  boxplot.width = 0.3,
  boxplot.alpha = 0.2,
  centrality = TRUE,
  centrality.fn = "median"
)
```

## Arguments

- df:

  data.frame. Long-format data with `.value` and the group column.

- group_col:

  character. Column name for the grouping variable (factor).

- title:

  character. Plot title (typically the gene name).

- x.title:

  character. X-axis title (the grouping variable name).

- colors:

  character. Named vector of colors for groups.

- min.nonzero:

  numeric. Minimum expressing cells for a group to be kept in the panel.
  Groups below it are dropped entirely.

- pt.size:

  numeric. Jitter point size.

- pt.alpha:

  numeric. Jitter point alpha.

- jitter.width:

  numeric. Horizontal jitter width.

- violin.alpha:

  numeric. Violin fill alpha.

- violin.width:

  numeric. Maximum violin width; every violin is scaled to this width
  (`scale = "width"`).

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

A ggplot2 object, or NULL when no group reaches `min.nonzero`.
