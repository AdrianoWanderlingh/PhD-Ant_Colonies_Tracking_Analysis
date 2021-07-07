rm()

set.seed(20200618)

# Create example data frame
newDF <- function() {
  observations <- rpois(1, lambda = 10)
  data.frame(group = sample(letters[1:3], observations, replace = TRUE),
             measure= rnorm(observations))
}
newDF()

# Example nested list
ex1 <- list(
  a1 = list(
    b1 = newDF()
  ),
  a2 = list(
    b2 = newDF(),
    b3 = newDF(),
    b4 = newDF()
  ),
  a3 = list(
    b5 = newDF()
  )
)
str(ex1)

# Example function that needs to be applied to the entire data frame, not its columns

sum1 <- function(x) x$measure + 1

ex1_a2 <- ex1$a2
result2 <- lapply(ex1_a2, FUN= sum1)
head(result2[[1]])

