# tests/testthat/test-visual_multi.R

skip_if_not_installed("ggcorrplot")

make_q3 <- function(J) {
  set.seed(1)
  m <- matrix(runif(J * J, -0.3, 0.3), J, J,
              dimnames = list(paste0("i", 1:J), paste0("i", 1:J)))
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m) <- 1
  m
}

# The two segment layers are the last two layers (after ggcorrplot's own layers)
segment_layers <- function(p) {
  is_seg <- vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1))
  p$layers[is_seg]
}

test_that("returns a ggplot with a horizontal and a vertical segment layer", {
  p <- visual_multi(make_q3(4), c(1, 1, 2, 2))
  expect_s3_class(p, "ggplot")
  expect_length(segment_layers(p), 2)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("segments are placed at the domain boundaries", {
  # domains 1,1 | 2,2,2 | 3  -> breaks at items 3 and 6
  Q3 <- make_q3(6)
  p <- visual_multi(Q3, c(1, 1, 2, 2, 2, 3))
  h <- segment_layers(p)[[1]]$data
  v <- segment_layers(p)[[2]]$data
  
  expect_equal(h$x,    c(2.5, 5.5))
  expect_equal(h$xend, c(5.5, 5.5))     # n - 0.5, here n = 6
  expect_equal(h$y,    c(3.5, 6.5))
  expect_equal(h$yend, h$y)             # horizontal
  
  expect_equal(v$x,    c(2.5, 5.5))
  expect_equal(v$xend, v$x)             # vertical
  expect_equal(v$y,    c(0.5, 0.5))
  expect_equal(v$yend, c(3.5, 6.5))
})

test_that("number of separator lines is the number of domains minus one", {
  for (itemD in list(c(1, 1, 2, 2), c(1, 2, 3, 4), c(1, 1, 1, 2, 2, 3))) {
    p <- visual_multi(make_q3(length(itemD)), itemD)
    expect_equal(nrow(segment_layers(p)[[1]]$data), length(unique(itemD)) - 1)
    expect_equal(nrow(segment_layers(p)[[2]]$data), length(unique(itemD)) - 1)
  }
})

test_that("domains are defined by consecutive runs, not by label", {
  # documents current behaviour: non-contiguous labels give extra lines
  p <- visual_multi(make_q3(4), c(1, 2, 1, 2))
  expect_equal(nrow(segment_layers(p)[[1]]$data), 3)
})

test_that("axis text, ticks and titles are removed", {
  th <- visual_multi(make_q3(4), c(1, 1, 2, 2))$theme
  for (el in c("axis.text.x", "axis.text.y", "axis.ticks.x",
               "axis.ticks.y", "axis.title.x", "axis.title.y")) {
    expect_s3_class(th[[el]], "element_blank")
  }
})

test_that("input matrix is not modified", {
  Q3 <- make_q3(4); before <- Q3
  visual_multi(Q3, c(1, 1, 2, 2))
  expect_identical(Q3, before)
})