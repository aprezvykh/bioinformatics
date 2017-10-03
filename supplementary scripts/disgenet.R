library("xlsx")
down <- read.xlsx('~/mouse significant trends.xlsx', header = TRUE, sheetIndex = 1)
up <- read.xlsx('~/mouse significant trends.xlsx', header = TRUE, sheetIndex = 2)
DisGeNET <- function(input, entity = c("gene", "disease"), identifier) {
  loadNamespace("httr")
  
  stopifnot(is.character(input))
  entity <- match.arg(entity)
  stopifnot(is.character(identifier), length(identifier) == 1L)
  
  STR <- switch(entity, gene={
    if (!identifier %in% c("entrez", "hgnc"))
      stop("entity='gene' 'identifier' must be 'entrez' or 'hgnc'")
    if (identifier == "entrez")
      "c2.geneId"
    else                            # identifier = 'hgnc'
      "c2.name"
  }, disease={
    if (!identifier %in% c("cui", "mesh", "omim"))
      stop("entity='disease' 'identifier' must be 'cui', 'mesh' or 'omim'")
    paste0("c1.", identifier)
  })
  
  url <- "http://www.disgenet.org/oql"
  terms <- paste(sprintf("'%s'", input), collapse=", ")
  
  oql <- paste0(
    "DEFINE
            c0='/data/gene_disease_score_onexus',
            c1='/data/diseases',
            c2='/data/genes',
            c3='/data/sources'
        ON
            'http://bitbucket.org/janis_pi/disgenet_onexus.git'
        SELECT
            c1 (cui, name, diseaseClassName, STY, MESH, omimInt),
            c2 (geneId, name, uniprotId, description, pathName, pantherName),
            c0 (score, pmids)
        FROM
            c0
        WHERE
            (c3 = 'ALL' AND ", STR, " IN (",  terms, ")
        ORDER BY ",
    STR, ", c0.score DESC")
  
  response <- httr::POST(url, body=oql)
  httr::stop_for_status(response)
  tbl <- read.csv(text=httr::content(response), header=TRUE, sep="\t")
  
  bad <- !input %in% tbl$c2.name
  if (any(bad))
    warning("entitites not in DisGeNET:\n  ",
            paste(sQuote(input[bad]), collapse=", "),
            call.=FALSE)
  
  tbl
}

input <- c("CDK1", "CDK1A", "CDK2")
result <- DisGeNET(input, 'gene', 'hgnc')
head(result)
  