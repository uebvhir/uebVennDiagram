###################################################
###################################################
#' A createVennEuler Function
#'
#' Generates Venn and Euler diagram from lists of elements (top tables) with the corresponding lists of shared and non-shared elements
#' @param topTabs = list of data frames (which are expected to be in "topTable" shape) to be used as input data
#' @param compNames = labels (character vector), corresponding with each input list, to be used in the plots and output table
#' @param label = name (character) to be included as global label in all output files created
#' @param colFeat = name (character) of the column where to look for the elements to be compared (such as "Name", "Symbol" or "GeneID")
#' @param colPVal = name (character) of the column with p-values to be used as filtering criteria, such as "Adj.p.val" or "P.value"
#' @param pval = threshold (numeric) value to be applied for filtering by p-value (0.05 means to filter out cases with p-val>0.05)
#' @param pltR = TRUE/FALSE, to determine if the plots should be done directly in R
#' @param pltPdf = TRUE/FALSE, to determine if the plots should be done in an output PDF file
#' @param eul = TRUE/FALSE, to determine whether or not to plot the Euler diagram
#' @param venn = TRUE/FALSE, to determine whether or not to plot the Venn diagram
#' @param csv= TRUE/FALSE, to determine whether or not to write the table of shared elements as csv output file
#' @export createVennEuler
#' @import VennDiagram
#' @import venneuler
#' @import grid
#' @author Miriam Mota <miriam.mota@vhir.org> and Ferran Brianso <ferran.brianso@vhir.org>
#' @examples
#' ## load("topTableList.RData") # this loads the topTableList data object
#' ## load("compNamesList.RData") # this loads the topTableList data object
#' ## compNamesList <- c("A1vsB1","A1vsC1","B1vsC1","A2vsB2","A2vsC2")
#' sharedElems <- createVennEuler(topTabs = topTableList,
#'                                compNames = compNamesList,
#'                                label = "demo",
#'                                colFeat = "Name",
#'                                colPVal = "P.Value",
#'                                pval = 0.01,
#'                                pltR = TRUE,
#'                                pltPdf = TRUE,
#'                                eul = TRUE,
#'                                venn = TRUE,
#'                                csv = TRUE)
#' @return sharedElements (data.frame): table with elements shared by each combination of lists provided
#' @keywords Venn Euler diagram sets union intersection lists VennDiagram EulerDiagram venneuler
#' @references Hanbo Chen (2014). VennDiagram: Generate high-resolution Venn and Euler plots. R package version 1.6.9. https://CRAN.R-project.org/package=VennDiagram
#' @references Lee Wilkinson (2011). venneuler: Venn and Euler Diagrams. R package version 1.1-0. https://CRAN.R-project.org/package=venneuler



#####################################################################
###### MAIN FUNCTION createVennEuler()
#####################################################################
createVennEuler <- function(topTabs, compNames, 
                            label = "selected", 
                            colFeat = "X", colPVal = "P.Value", pval = 0.05, 
                            pltR = TRUE, pltPdf = TRUE, venn = TRUE, eul = TRUE, csv = TRUE){

  ## Initializing lists
  list_genes_sel <- list()
  
  ## Reading input data
  for (i in 1:length(topTabs)) {
    colpval <- which(names(topTabs[[i]]) == colPVal)
    colFeature <- which(names(topTabs[[i]]) == colFeat)
    list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colpval] < pval])
  }

  ## Creating Venn Diagram
  if (venn) {
    venn.plot <- venn.diagram(list_genes_sel,
                              category.names = compNames,
                              fill = rainbow(length(compNames)),
                              #fill = c("tomato", "orchid4", "turquoise3"),
                              alpha = 0.50,
                              resolution = 600,
                              cat.cex = 0.9,
                              main = paste0("Venn diagram (" , colPVal, " < ", pval,")"),
                              filename = NULL)
    if (pltPdf) {
      pdf(paste0("VennDiagram.", label, ".", colPVal, pval, ".pdf"))
      grid.draw(venn.plot)
      dev.off()
    }
    if (pltR) {grid.draw(venn.plot)}
  }
  
  ## Creating Euler Diagram
  if (eul) {
    set <- NULL
    for (i in 1:length(compNames)) {
      set <- c(set, rep(compNames[i],length(list_genes_sel[[i]])))
    }
    v <- venneuler(data.frame(elements = c(unlist(list_genes_sel)),
                              sets = set))
    if (pltPdf) {
      pdf(paste0("EulerDiagram.", label, ".", colPVal, pval, ".pdf"))
      plot(v, main = paste0("Euler diagram (" , colPVal, " < ", pval,")"))
      dev.off()
    }
    if (pltR) {plot(v, main = paste0("Euler diagram (" , colPVal, " < ", pval,")"))}
  }

  ## Obtaining lists of combined elements and computing shared elements
  names(list_genes_sel) <- paste0(compNames, "_")
  combs <-  unlist(lapply(1:length(list_genes_sel),
                          function(j) combn(names(list_genes_sel), j, simplify = FALSE)),
                   recursive = FALSE)
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
  #str(combs)

  elements <- lapply(combs, function(i) Setdiff(list_genes_sel[i], list_genes_sel[setdiff(names(list_genes_sel), i)]))
  n.elements <- sapply(elements, length)
  list_res <- list(elements = elements, n.elements = n.elements)

  seq.max <- seq_len(max(n.elements))
  mat <- sapply(elements, "[", i = seq.max)
  mat[is.na(mat)] <- ""
  sharedElements <- rbind(t(data.frame(n.elements)),data.frame(mat))

  ## Writing table of shared elements as csv
  if (csv) {
    write.csv(sharedElements,
              file = paste0("sharedElements.", label, ".", colPVal, pval, ".csv"), 
              row.names = FALSE)
  }
  
  ## Returning table as data.frame
  return(sharedElements)
}
#####################################################################


###################################################################
###### Internal functions used to extract the lists of shared elements 
###### both x and y are, in all cases, lists of character strings
###################################################################
Intersect <- function(x) {
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function(x) {
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function(x, y) {
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}
#####################################################################


#####################################################################
###### Example call of createVennEuler() function
#####################################################################
#
##topTableList <- list()
##topTableList[[1]] <- topTable.A1vsB1 
##topTableList[[2]] <- topTable.A1vsC1 
##topTableList[[3]] <- topTable.B1vsC1 
##topTableList[[4]] <- topTable.A2vsB2
##topTableList[[5]] <- topTable.A2vsC2
##save(topTableList, file = "topTableList.RData")
#
#require(VennDiagram)
#require(venneuler)
#
#load("topTableList.RData")
#compNamesList <- c("A1vsB1","A1vsC1","B1vsC1","A2vsB2","A2vsC2")
#
#sharedElems <- createVennEuler(topTabs = topTableList,
#                               compNames = compNamesList,
#                               label = "test",
#                               colFeat = "Name",
#                               colPVal = "P.Value",
#                               pval = 0.01,
#                               pltR = TRUE,
#                               pltPdf = TRUE,
#                               eul = TRUE,
#                               venn = TRUE,
#                               csv = TRUE)
#

#####################################################################

