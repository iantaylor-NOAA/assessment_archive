# make simple barplot showing assessments per year

library(ggplot2)
info <- read.csv("inst/extdata/Stock Assessment History - main 24 Feb 2025.csv")

unique(info$Type) |> t() |> t()
#       [,1]
#  [1,] "Full"
#  [2,] "Update"
#  [3,] "DR"
#  [4,] "DL"
#  [5,] "DM"
#  [6,] "Full (Rejected)"
#  [7,] "DM (not reviewed)"
#  [8,] "DL (Rejected)"
#  [9,] "Full (Withdrawn)"
# [10,] "COU"

info$Type2 <- factor(
  info$Type,
  levels =
    rev(c(
      "Full",
      "Full (Rejected)",
      "Full (Withdrawn)",
      "Update",
      "DM",
      "DM (not reviewed)",
      "DL",
      "DL (Rejected)",
      "COU",
      "DR"
    ))
)

info |>
  dplyr::filter(grepl("Full", Type) |
    grepl("Update", Type) |
    grepl("DM", Type)) |>
  ggplot(aes(Year, fill = Type2)) +
  geom_histogram(breaks = seq(2002.5, 2024.5, 1))
