# IN CASE THE LOCATION OF THE SAVE CHANGES, USE THESE TWO SHELL COMMANDS TO CHANGE THE PATHS IN THE R PROGRAMS
# filelist=($(find ~/Desktop/mouseproject/software | xargs grep '/share/adl/pnlong/software/useful_functions.Rsave' | cut -d : -f 1 | tr '\n' ' '))
# for file in "${filelist[@]}"; do sed -i.bak 's+/share/adl/pnlong/software/useful_functions.Rsave+/share/adl/pnlong/mouseproject/software/useful_functions.Rsave+g' $file; done

# a function to cleanly extract the arguments of an R script, returning a named vector where the names are argument names and values are argument values
get_args <- function() {
  args <- t(matrix(unlist(strsplit(x = commandArgs(trailingOnly = TRUE),
                                   split = "=")),
                   nrow = 2))
  
  vec <- args[,2]
  names(vec) <- args[,1]
  return(vec)
}

# a function to separate a vector into n parts that are (somewhat) equally spaced apart
vector_separate_at_indices <- function(x, n) {
  sequence <- round(
    seq(from = (length(x) %/% n),
        to = length(x),
        length.out = n))
  
  return(x[sequence])
}

# a new pipe operator: saves the current values in a pipe to a given name while allowing the pipe to continue (2ND_NEW VARIABLE <- x %>% filter() %->% NEW_VARIABLE %>% pull())
`%->%` <- function(value, x)
{
  library(lazyeval)
  x <- lazyeval::lazy(x)
  assign(deparse(x$expr), value, x$env)
  value
}

# a function to convert a 2D vector (list of vectors/lists) into a tibble
two_dimensional_vector_to_tibble <- function(x, names = NA) {
  output <- x %>%
    unlist() %>%
    matrix(nrow = length(x),
           ncol = length(x[[1]]),
           byrow = TRUE) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    as_tibble()
  
  if (!mean(is.na(names))) colnames(output) <- names
  
  return(output)
}

wc_l <- function(file) {
  if (!file.exists(file)) {
    return(0)
  } else {
    x <- file(file)
    on.exit(close(x))
    count = length(readLines(con = x))
    
    return(count)
  }
}

correct_LOD_values <- function(x) {
  # x is a LOD_table tibble with 3 columns: TEST_SNP, LOD_METHOD, LOD_VALUE
  
  # check if it is snpbased if the LOD_METHOD is only ever equal to 1
  if (mean(unique(x$LOD_METHOD) == c(1)) == 1) { # so if it is snpbased LOD table
    return(x)
  } else if (mean(sort(unique(x$LOD_METHOD)) == c(1, 2)) == 1) { # so if it is hapbased LOD table
    return(x %>%
             pivot_wider(names_from = LOD_METHOD, values_from = LOD_VALUE) %>%
             mutate("1" = map2_dbl(.x = `1`,
                                   .y = `2`,
                                   .f = ~ if_else(condition = is.na(.y), # if LOD_2 is NA
                                                  true = as.numeric(.x), # use LOD_1
                                                  false = as.numeric(.y)))) %>% # otherwise, use LOD_2
             pivot_longer(cols = `1`, names_to = "LOD_METHOD", values_to = "LOD_VALUE") %>%
             mutate(LOD_METHOD = as.numeric(LOD_METHOD)) %>%
             select(TEST_SNP, LOD_METHOD, LOD_VALUE) %>%
             filter(LOD_METHOD == 1)
    )
  }
}

# SAVE TO BE ACCESSED BY ANY PROGRAM
save(list = c("get_args", "vector_separate_at_indices", "%->%", "two_dimensional_vector_to_tibble", "wc_l", "correct_LOD_values"), file = "/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
