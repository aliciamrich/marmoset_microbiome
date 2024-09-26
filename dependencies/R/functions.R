sort.taxa <- function(df) {
  df <- df %>%
    arrange(superkingdom,
            phylum,
            class,
            order,
            family,
            genus,
            species)
  return(df)
}

fix.strings <-  function(df) {
  df <- df %>%
    mutate(across(where(is.character), ~str_remove_all(.x, fixed("'")))) %>%
    mutate(across(where(is.character), ~str_remove_all(.x, fixed("[")))) %>%
    mutate(across(where(is.character), ~str_remove_all(.x, fixed("]")))) %>%
    mutate(across(where(is.character), ~str_trim(.x, "both"))) %>%
    mutate(across(where(is.character), ~str_squish(.x)))
  return(df)
}

export.list <- function(df, filename) {
  write.table(df,
              paste0(params$local, "dataframes/", filename, ".txt"),
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}


backup.df  <- function(df, filename) {
  write.table(df,
              paste0(params$local, "dataframes/", filename, "_", Sys.Date(), ".tsv"),
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE)
}


since.start <- function(date.col, units) {
  as.numeric(as.period(interval(ymd(params$day1), date.col), unit = units), units)
}


check.duplicates <- function(df2, df, group) {
  df2 <- df %>%
    group_by(group) %>%
    filter(n() > 1) %>%
    ungroup()
}

group.taxonomy <- function(df) {
  group_by(.data = df,
           superkingdom,
           phylum,
           class,
           order,
           family,
           genus,
           species)
}

read.tables <- function(file) {
  data <- read.table(file, 
                     header = TRUE, 
                     sep = "\t", 
                     stringsAsFactors = FALSE) %>%
    tibble()
  return(data)
}


read.recent.version.csv <- function(directory, pattern) {
  files             <- list.files(path       = paste0(params$local, "/", directory, "/"), 
                                  pattern    = paste0(pattern, "\\d{4}-\\d{1,2}-\\d{1,2}\\.csv"), 
                                  full.names = TRUE)
  dates             <- gsub(".*_(\\d{4}-\\d{1,2}-\\d{1,2})\\.csv", "\\1", files)
  dates             <- as.Date(dates, format = "%Y-%m-%d")
  most_recent_index <- which.max(dates)
  most_recent_file  <- files[most_recent_index]
  data              <- read.csv(most_recent_file, header = TRUE)
  
  return(data)
}


read.recent.version.tsv <- function(directory, pattern) {
  files             <- list.files(path       = paste0(params$local, "/", directory, "/"), 
                                  pattern    = paste0(pattern, "\\d{4}-\\d{1,2}-\\d{1,2}\\.tsv"), 
                                  full.names = TRUE)
  dates             <- gsub(".*_(\\d{4}-\\d{1,2}-\\d{1,2})\\.tsv", "\\1", files)
  dates             <- as.Date(dates, format = "%Y-%m-%d")
  most_recent_index <- which.max(dates)
  most_recent_file  <- files[most_recent_index]
  data              <- read.table(most_recent_file, sep = "\t", header = T)
  
  return(data)
}

#Terminal Code as Text
open.job <- function(name, mem, hrs, cpus) {
  cat("paste to terminal:\n\nsrun --partition=guest --nodes=1 --ntasks-per-node=1", paste0(" --job-name=", name, " --mem=", mem, "GB --time=", hrs, ":00:00 --cpus-per-task=", cpus), "--pty $SHELL\n\n")
}

conda.env <- function(env) {
  cat("paste to terminal:\n\ncd", params$work_dir, "\nmodule load anaconda\nconda activate", env, "\n\n")
}

load.pkg <- function(pkg, wd) {
  cat("paste to terminal:\n\ncd", paste0(params$work_dir, wd), "\nmodule load", pkg, "\n\n")
}

sbatch <- function(script) {
  cat("paste to terminal:\n\ncd", params$work_dir, "\n\nsbatch", paste0(script, ".sh"), "\n\nsqueue -u aliciarich\n\n")
}

name.subdir <- function(name, subdir) {
  cat("paste to terminal:\n\n", paste0(name, "=", params$work_dir, subdir), "\n\n", paste0("mkdir -p", " $", name), "\n\n")
}

out.dir <- function(name) {
  cat("paste to terminal:\n\n", paste0(name, "=", params$work_dir, "/", name, "/", params$seqrun, "/", params$date), "\n\n", paste0("mkdir -p", " $", name), "\n\n")
}

read.tables <- function(file) {
  data <- read.table(file, 
                     header = TRUE, 
                     sep = "\t", 
                     stringsAsFactors = FALSE) %>%
    tibble()
  return(data)
}

read_alignment_file  <- function(file) {
  read.csv(file, stringsAsFactors = FALSE, fill = TRUE)
}

