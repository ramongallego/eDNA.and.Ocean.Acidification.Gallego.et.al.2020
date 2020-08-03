# A function that takes a long dataframe and turns it into a number of different objects:

# community table tibble.to.table
# 
# distance matrix tibble.to.matrix
# 
# Permanova object tibble.to.permanova
# 
# Environmental table tibble.to.env

tibble.to.table <- function (long.table,taxon,Abundance, sample.name, transformation, ...) {
  
  sample.name = rlang::enquo(sample.name)
  taxon = rlang::enquo(taxon) # quote taxon so it works with select
 
  Abundance = rlang::enquo(Abundance)
  
  
  # So now the eDNA index works let's see if we can use decostand
  
  
  
  cols = long.table %>% ungroup %>% dplyr::select(!!taxon) %>% distinct() %>% pull() # capture the unique values of the taxa
  
  
  long.table %>%
    ungroup %>% 
    select(!!taxon, !!Abundance, !!sample.name, ...) %>% # select at least three columns: nReads, taxa and the sampleID
    spread (key = !!taxon, value = !!Abundance, fill = 0) -> matrix_1 # Make a wide table, like vegan likes
  
  env <- colnames(matrix_1)[!colnames(matrix_1) %in% cols] # Which column has all the sample info
  
  dplyr::select (matrix_1, env) -> env # The values of the sample info in the order they appear in the wide table
  
  samples <- pull (matrix_1, !!sample.name)
  dplyr::select (matrix_1, cols) -> matrix_1 # select only the spp
  data.matrix(matrix_1) -> matrix_1 # make it a data matrix
  
  dimnames(matrix_1)[[1]] <- samples # name it with the sample info
  return (matrix_1)
  
  }

tibble_to_matrix <- function (long.table,taxon,Abundance, sample.name, distance = "bray", transformation, ...) {
  
  sample.name = rlang::enquo(sample.name)
  taxon = rlang::enquo(taxon) # quote taxon so it works with select
  choice = rlang::enquo(distance)
  Abundance = rlang::enquo(Abundance)
  
  
  # So now the eDNA index works let's see if we can use decostand
  
  
  
  cols = long.table %>% ungroup %>% dplyr::select(!!taxon) %>% distinct() %>% pull() # capture the unique values of the taxa
  
  
  long.table %>%
    ungroup %>% 
    select(!!taxon, !!Abundance, !!sample.name, ...) %>% # select at least three columns: nReads, taxa and the sampleID
    spread (key = !!taxon, value = !!Abundance, fill = 0) -> matrix_1 # Make a wide table, like vegan likes
  
  env <- colnames(matrix_1)[!colnames(matrix_1) %in% cols] # Which column has all the sample info
  
  dplyr::select (matrix_1, env) -> env # The values of the sample info in the order they appear in the wide table
  
  samples <- pull (matrix_1, !!sample.name)
  dplyr::select (matrix_1, cols) -> matrix_1 # select only the spp
  data.matrix(matrix_1) -> matrix_1 # make it a data matrix
  
  dimnames(matrix_1)[[1]] <- samples # name it with the sample info
  
  # Now that the matrix is ready to go, we can apply decostand if needed
  
  good.transformations <- c("total", "max", "frequency", "normalize", "range", 
                            "rank", "rrank", "standardize", "pa", "chi.square", "hellinger", 
                            "log")
  if (transformation == "sq") {
    sqrt(matrix_1) -> matrix_1
  }else{
  if(transformation %in%  good.transformations ){
    decostand(matrix_1, method = transformation) -> matrix_1
  }
  }
  #simper(matrix_1, env$Site)
  vegdist(matrix_1, method = distance) -> matrix_1 
  return(matrix_1)
  # calculate the dissimilarity matrix
  
#   adonis(matrix_1 ~  Site*Date, data= env)
  # as.tibble(subset(melt(as.matrix(matrix_1))))
}

tibble_to_permanova <- function (long.table,taxon,Abundance, sample.name, distance = "bray", transformation, formula = "",...) {
  
  sample.name = rlang::enquo(sample.name)
  taxon = rlang::enquo(taxon) # quote taxon so it works with select
  choice = rlang::enquo(distance)
  Abundance = rlang::enquo(Abundance)
  
  
  # So now the eDNA index works let's see if we can use decostand
  
  
  
  cols = long.table %>% ungroup %>% dplyr::select(!!taxon) %>% distinct() %>% pull() # capture the unique values of the taxa
  
  
  long.table %>%
    ungroup %>% 
    select(!!taxon, !!Abundance, !!sample.name, ...) %>% # select at least three columns: nReads, taxa and the sampleID
    spread (key = !!taxon, value = !!Abundance, fill = 0) -> matrix_1 # Make a wide table, like vegan likes
  
  env <- colnames(matrix_1)[!colnames(matrix_1) %in% cols] # Which column has all the sample info
  
  dplyr::select (matrix_1, env) -> env # The values of the sample info in the order they appear in the wide table
  
  samples <- pull (matrix_1, !!sample.name)
  dplyr::select (matrix_1, cols) -> matrix_1 # select only the spp
  data.matrix(matrix_1) -> matrix_1 # make it a data matrix
  
  dimnames(matrix_1)[[1]] <- samples # name it with the sample info
  
  # Now that the matrix is ready to go, we can apply decostand if needed
  
  good.transformations <- c("total", "max", "frequency", "normalize", "range", 
                            "rank", "rrank", "standardize", "pa", "chi.square", "hellinger", 
                            "log")
  if (transformation == "sq") {
    sqrt(matrix_1) -> matrix_1
  }else{
    if(transformation %in%  good.transformations ){
      decostand(matrix_1, method = transformation) -> matrix_1
    }
  }
  #simper(matrix_1, env$Site)
  vegdist(matrix_1, method = distance) -> matrix_1 
  return(matrix_1)
  # calculate the dissimilarity matrix
  
     adonis(matrix_1 ~  Site*Date, data= env)
  # as.tibble(subset(melt(as.matrix(matrix_1))))
}


tibble_to_env <- function (long.table,taxon,Abundance,sample.name, ...) {  # A function that returns the accompaining environmental matrix  
  
  # Enquos to allow me to play with the data
  
  sample.name = rlang::enquo(sample.name)
  taxon = rlang::enquo(taxon) # quote taxon so it works with select
  #choice = rlang::enquo(distance)
  Abundance = rlang::enquo(Abundance)
  cols = long.table %>% dplyr::select(!!taxon) %>% distinct() %>% pull() 
  #return(cols)
  long.table %>%
    ungroup %>% 
    select(!!taxon, !!Abundance, !!sample.name, ...) %>% 
    spread (key = !!taxon, value = !!Abundance, fill = 0) -> matrix_1
  # return (matrix_1)
  env <- colnames(matrix_1)[!colnames(matrix_1) %in% cols]
  # return(env)
  dplyr::select (matrix_1, env) -> env
  return(env)
}

