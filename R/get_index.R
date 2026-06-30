# code modified from initial response by Gemini to Google Search
# "read list of files in directory within a branch of a github repository using an R script"
# then modified by Copilot to fix a bug and deal with nested directories
# then modified by Ian from there

list_est_by_area_files <- function() {
  url <- "https://api.github.com/repos/pfmc-assessments/auto-indexwc/git/trees/autogen-results?recursive=1"

  request_obj <- request(url) |> req_perform()
  content_json <- resp_body_json(request_obj)

  all_tree_items <- content_json$tree
  file_items <- keep(all_tree_items, ~ .x$type == "blob")
  file_paths <- map_chr(file_items, ~ .x$path)

  grep("est_by_area.csv", file_paths, value = TRUE)
}

get_est_by_area_csv <- function(common_name) {
  common_name <- gsub(" ", "_", common_name) |> tolower()
  est_files <- list_est_by_area_files()

  matched_files <- grep(common_name, est_files, value = TRUE)
  if (length(matched_files) == 0) {
    cli::cli_abort("No file found for common name: {common_name}")
  }
  if (length(matched_files) > 1) {
    cli::cli_abort("Multiple files found for common name: {common_name}")
  }

  raw_url <- paste0("https://raw.githubusercontent.com/pfmc-assessments/auto-indexwc/autogen-results/", matched_files)
  read.csv(raw_url, stringsAsFactors = FALSE)
}
