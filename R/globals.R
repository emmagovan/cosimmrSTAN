if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    names = c(
      "n_samples",
      "Proportion",
      "Source",
      "density",
      "Beta"
    ),
    package = "cosimmrSTAN",
    add = FALSE
  )
}
