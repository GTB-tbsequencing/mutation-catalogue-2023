## This function executes the original SOLO algorithm, following Timothy Walker's implementation (PMID 26116186).
## The inputTab has columns sample_id, mutation, phenotype and [SOnly or het]; maxIter is a positive integer or Inf.
## If removeSOnly is FALSE, the first stage of classifying and removing the variants occurring only in S is skipped.
originalSOLO = function(inputTab, maxIter = MAX_ITER, removeSOnly = TRUE) {
  useHet = ("het" %in% colnames(inputTab))
  ## Construct bipartite directed graph with each edge having a sample as tail and a mutation it contains as head
  G = graph_from_data_frame(inputTab, directed = TRUE)
  ## The mutations are in column 2 
  mutations = na.omit(unique(inputTab[[2]]))
  ## The output matrix contains, for each mutation, its classification, iteration classified, SOLO R/total counts
  output = matrix(NA, length(mutations), 4, dimnames = list(mutations, c("class", "iter", "Tcnt", "Rcnt")))
  iter = 1L
  currentSetS = c()
  ## The loop proceeds till either maxIter is reached, all mutations are classified or no new ones are classified S
  while (iter <= maxIter && length(mutations) > 0 && (iter == 1L || length(currentSetS) > 0)) {
    ## Remove mutations classified S from the list of active mutations and the graph, along with incident edges
    G %<>% delete.vertices(currentSetS)
    mutations %<>% setdiff(currentSetS)
    ## Compute the total out-degree, which is 0 for mutations and is the number of mutations contained for samples 
    Tdegree = igraph::degree(G, mode = "out")
    ## The singleton samples are those containing a single mutation
    currentSingletons = names(Tdegree)[Tdegree == 1]
    ## The initial set of mutations classified as S is empty
    currentSetS = c()
    ## Start by classifying all the mutations occurring only in S isolates as S at the first iteration, if necessary
    if (iter == 1L && removeSOnly) {
      if (useHet) {
        ## Use what we know; extract the subgraph containing only the R samples and non-het mutations
        Rsubgraph = subgraph.edges(G, eids = E(G)[phenotype == "R" & !het], delete.vertices = FALSE)
        ## Compute the degrees in this subgraph; degree 0 means a mutation never occurs in R samples as a non-het
        Rdegree = igraph::degree(Rsubgraph)
        ## Also extract the subgraph containing only the S samples and non-het mutations
        Ssubgraph = subgraph.edges(G, eids = E(G)[phenotype == "S" & !het], delete.vertices = FALSE)
        ## Compute the degrees in this subgraph; degree > 0 means a mutation occurs in some S samples as a non-het
        Sdegree = igraph::degree(Ssubgraph)
        ## Classify all the mutations satisfying both criteria as S
        currentSetS = mutations[Rdegree[mutations] == 0 & Sdegree[mutations] > 0]
      } else { 
        ## Extract the mutations marked SOnly; some are excluded due to hets
        currentSetS = inputTab %>% 
          dplyr::filter(SOnly) %>%
          pull(2) %>%
          unique()
      }
    }
    ## Compute the edges connecting the singleton samples to their mutations, and get their endpoints as a matrix
    if (useHet) {
      currentEdgeSeq = E(G)[.from(currentSingletons) & !het]
    } else {
      currentEdgeSeq = E(G)[.from(currentSingletons)]
    }
    currentEdges   = ends(G, currentEdgeSeq)
    ## The full counts are obtained by looking at all the mutations occurring in singletons (i.e. solo mutations)
    currentTCounts  = table(currentEdges[,2])
    currentSet      = names(currentTCounts)
    ## Identify the R edges among them (they originate from R samples) and classify their mutation endpoints as R
    currentREdges  = (currentEdgeSeq$phenotype == "R")
    currentSetR    = currentEdges[currentREdges, 2]
    ## The mutations found in R edges are counted as SOLO_R; note that after tabulation, the names become unique
    currentRCounts  = table(currentSetR)
    currentSetR     = names(currentRCounts)
    ## The remaining mutations tare found in singletons are only found in S samples and so get classified as S
    currentSetS %<>% c(setdiff(currentSet, currentSetR))
    ## Record the classification of the newly classified mutations as S or R, accordingly; leave the rest alone
    output[currentSetS, "class"] %<>% ifelse(is.na(.), SCODE, .)
    output[currentSetR, "class"] %<>% ifelse(is.na(.), RCODE, .)
    ## Record the iteration of the newly classified mutations, but leave the rest alone; increment iteration counter
    output[c(currentSetR, currentSetS), "iter"] %<>% ifelse(is.na(.), iter, .)
    iter %<>% add(1)
    ## Record the SOLO counts of the newly classified mutations as the current SOLO counts; leave the rest alone
    output[currentSet,  "Tcnt"] %<>% ifelse(is.na(.), currentTCounts, .)
    output[currentSetR, "Rcnt"] %<>% ifelse(is.na(.), currentRCounts, .)
  }
  ## Convert the output into a tibble, replace NA values with 0 and add the original mutation names as 1st column
  initMutations = rownames(output)
  output %<>%
    as_tibble() %>%
    mutate_all(~{replace_na(., 0) %>% as.integer()}) %>%
    mutate(mutation = initMutations, .before = 1)
  output
}

## This function executes preprocessing, the SOLO algorithm (for maxIter steps), and postprocessing of inputTab.
## The column names of inputTab are required to be drug, sample_id, variant, and phenotype; an optional fifth
## column can contain information about those mutations that only occur in S isolates, which allows het removal;
## alternatively, to get the most general case, it can contain information about which of the mutations are hets.
## If the stage is also specified, then it is added as an additional variable at the end of the output table.
## If removeSOnly is FALSE, the first stage of classifying and removing the variants occurring only in S is skipped.
runSOLOPipeline = function(inputTab, maxIter, stage = NULL, removeSOnly = TRUE) {
  stopifnot(ncol(inputTab) %in% 4:5 && all(colnames(inputTab)[1:4] == c("drug", "sample_id", "variant", "phenotype")))
  ## Initialise the output table
  outputs = tibble()
  ## For each drug, process the first maxIter steps of the original solo algorithm
  for (drugName in unique(inputTab$drug)) {
    curTab = inputTab %>%
      dplyr::filter(drug == drugName) %>%
      select(-drug)
    outputs %<>% 
      bind_rows(originalSOLO(inputTab = curTab, maxIter = maxIter, removeSOnly = removeSOnly) %>%
                  mutate(drug = drugName, .before = 1))
  }
  ## Update the outputs for compatibility with the original format and return them
  outputs %<>%
    rename(variant = "mutation") %>%
    mutate_at("class", ~{recode(as.factor(.), !!!CODE_KEY) %>% as.character()}) %>%
    mutate(Scnt = Tcnt - Rcnt, Tcnt = NULL)
  if (!is.null(stage)) {
    outputs %<>%
      mutate(stage = stage)
  }
  outputs
}

testSOLO = function() {
  inputTabAlt = tibble(drug = rep("Dummy", 7),
                    sample_id = rep(1:4, c(3, 1, 2, 1)), 
                    variant = LETTERS[c(1, 2, 3, 1, 3, 4, 2)], 
                    phenotype = rep(c("R", "S"), c(4, 3)),
                    SOnly = c(rep(FALSE, 5), TRUE, FALSE))
  outputTabAlt = runSOLOPipeline(inputTabAlt, 2L)
  stopifnot(all(outputTabAlt$class == c("R", "S", "S", "S")))
  stopifnot(all(outputTabAlt$iter  == c(1, 1, 2, 1)))
  stopifnot(all(outputTabAlt$Scnt  == c(0, 1, 1, 0)))
  stopifnot(all(outputTabAlt$Rcnt  == c(1, 0, 0, 0)))
  inputTabHet = inputTabAlt
  inputTabHet$SOnly = FALSE
  outputTabHet = runSOLOPipeline(inputTabHet, 2L)
  stopifnot(all(outputTabHet$class == c("R", "S", "U", "U")))
  stopifnot(all(outputTabHet$iter  == c(1, 1, 0, 0)))
  stopifnot(all(outputTabHet$Scnt  == c(0, 1, 0, 0)))
  stopifnot(all(outputTabHet$Rcnt  == c(1, 0, 0, 0)))
  trickyCase = tibble(drug = rep("Dummy", 6),
                      sample_id = rep(1:3, each = 2),
                      variant = rep(c("x", "y"), 3),
                      phenotype = c(rep("S", 4), rep("R", 2)),
                      het = c(rep(FALSE, 5), TRUE))
  outputTricky = runSOLOPipeline(trickyCase, 2L)
  stopifnot(all(outputTricky$class == c("R", "S")))
  stopifnot(all(outputTricky$iter  == c(2, 1)))
  stopifnot(all(outputTricky$Scnt  == c(2, 0)))
  stopifnot(all(outputTricky$Rcnt  == c(1, 0)))
  return(TRUE)
}
