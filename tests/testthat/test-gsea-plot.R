test_that("gsea_pbmc3k_cd8t dataset is available", {
  data(gsea_pbmc3k_cd8t, package = "BadranSeq")
  expect_true(exists("gsea_pbmc3k_cd8t"))
  expect_s4_class(gsea_pbmc3k_cd8t, "gseaResult")
  expect_true(nrow(gsea_pbmc3k_cd8t@result) >= 1)
})

test_that("do_GseaPlot returns a ggplot object", {
  data(gsea_pbmc3k_cd8t, package = "BadranSeq")

  p <- do_GseaPlot(
    gsea_pbmc3k_cd8t,
    analysis_name = "Test Analysis",
    simplify_go = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("do_GseaPlot handles empty results gracefully", {
  # Create a fake empty gseaResult-like object
  empty_result <- gsea_pbmc3k_cd8t
  empty_result@result <- gsea_pbmc3k_cd8t@result[0, ]

  expect_message(
    do_GseaPlot(empty_result, analysis_name = "Empty Test", simplify_go = FALSE),
    "No significant pathways"
  )

  p <- do_GseaPlot(empty_result, analysis_name = "Empty Test", simplify_go = FALSE)
  expect_null(p)
})

test_that("do_GseaPlot selection parameter overrides auto selection", {
  data(gsea_pbmc3k_cd8t, package = "BadranSeq")

  pathway_name <- gsea_pbmc3k_cd8t@result$Description[1]

  p <- do_GseaPlot(
    gsea_pbmc3k_cd8t,
    selection = pathway_name,
    simplify_go = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("do_GseaPlot warns on missing selected pathways (partial match)", {
  data(gsea_pbmc3k_cd8t, package = "BadranSeq")

  # Mix a real pathway with a fake one to trigger partial-match warning
  real_pathway <- gsea_pbmc3k_cd8t@result$Description[1]
  expect_warning(
    do_GseaPlot(
      gsea_pbmc3k_cd8t,
      selection = c(real_pathway, "Totally Fake Pathway"),
      simplify_go = FALSE
    ),
    "The following selected pathways were not found"
  )
})

test_that("do_GseaPlot returns NULL when no selected pathways match", {
  data(gsea_pbmc3k_cd8t, package = "BadranSeq")

  expect_message(
    do_GseaPlot(
      gsea_pbmc3k_cd8t,
      selection = c("Totally Fake Pathway"),
      simplify_go = FALSE
    ),
    "None of the selected pathways found"
  )

  p <- do_GseaPlot(
    gsea_pbmc3k_cd8t,
    selection = c("Totally Fake Pathway"),
    simplify_go = FALSE
  )
  expect_null(p)
})

test_that("do_GseaPlot supports viridis fill scale", {
  data(gsea_pbmc3k_cd8t, package = "BadranSeq")

  p <- do_GseaPlot(
    gsea_pbmc3k_cd8t,
    fill_scale = "viridis",
    simplify_go = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("do_GseaPlot respects number_of_pathways", {
  data(gsea_pbmc3k_cd8t, package = "BadranSeq")

  # This dataset only has 2 pathways, so n=1 should show 1 (if possible)
  p <- do_GseaPlot(
    gsea_pbmc3k_cd8t,
    number_of_pathways = 1,
    simplify_go = FALSE
  )

  expect_s3_class(p, "ggplot")
})
