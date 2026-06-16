# A function that adds two numbers, adds 2 and then multiplies by c
funky_maths <- function (a, b, c){
  # sum a and b
  sum_ab <- a + b
  # add 2 to the sum
  sum_ab_plus_2 <- sum_ab + 2
  # multiply the result by c
  result <- sum_ab_plus_2 * c
  return(result)
}
