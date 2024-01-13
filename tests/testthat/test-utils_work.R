test_that("pairwise_sim works", {
  mat1 = matrix(rnorm(100), nrow = 10)
  mat2 = matrix(rnorm(100), nrow = 10)

  sim_cols <- pairwise_sim(mat1, mat2, which = "cols")
  sim_rows <- pairwise_sim(mat1, mat2, which = "rows")

  expect_equal(
      sim_cols[1,2],
      lsa::cosine(mat1[,1], mat2[,2])[1,1]
  )

  expect_equal(
      sim_cols[2,1],
      lsa::cosine(mat1[,2], mat2[,1])[1,1]
  )

  expect_equal(
      sim_rows[1,2],
      lsa::cosine(mat1[1,], mat2[2,])[1,1]
  )

  expect_equal(
      sim_rows[2,1],
      lsa::cosine(mat1[2,], mat2[1,])[1,1]
  )
})
