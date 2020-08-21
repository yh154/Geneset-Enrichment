
log <- file(snakemake@log[[1]])
sink(log, append = TRUE)
sink(log, append = TRUE, type = "message")

# prepare required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

list.of.pkgs <- c("clusterProfiler","pathview", "enrichplot","DOSE","ggplot2", "msigdbr")
new.pkgs <- list.of.pkgs[!(list.of.pkgs %in% installed.packages()[,"Package"])]
if(length(new.pkgs))
  BiocManager::install(new.pkgs)

sapply(list.of.pkgs, library, character.only = TRUE)

# prepare input gene list (named vector of LogFCs)
fun_prep_input <- function(gene_list){
  gene_list = na.omit(gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
}

# GO term enrichment
fun_gseGO <- function(gene_list, ont, keyType, minGS,maxGS, pvalueCutoff, orgdb, pAdjMethod){
      gse <- clusterProfiler::gseGO(
             geneList=gene_list,
             ont =ont,
             keyType = keyType,
             nPerm = 10000,
             minGSSize = minGS,
             maxGSSize = maxGS,
             pvalueCutoff = pvalueCutoff,
             verbose = TRUE,
             OrgDb = orgdb,
             pAdjustMethod = pAdjMethod
       )
}


# KEGG enrichment
fun_convert_geneid <- function(gene_list, fromType, toType, orgdb){
    #  Convert gene IDs for gseKEGG function
    new.ids<-clusterProfiler::bitr(names(gene_list),fromType = fromType,toType = toType,OrgDb=orgdb)

    new.gene.list <- merge(new.ids, as.data.frame(gene_list), by.x =1 , by.y='row.names')
    new.gene.list <- new.gene.list[order(abs(new.gene.list[['gene_list']]),decreasing=T),]
    new.gene.list <- new.gene.list[!duplicated(new.gene.list[c(fromType)]),]
    new.gene.list <- new.gene.list[!duplicated(new.gene.list[c(toType)]),]

    kegg_gene_list <- new.gene.list[['gene_list']]
    names(kegg_gene_list) <- new.gene.list[[toType]]
    kegg_gene_list<-na.omit(kegg_gene_list)
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
}

fun_gseKEGG <- function(kegg_gene_list, keyType, kegg_organism, minGS, maxGS, pvalueCutoff, pAdjustMethod){
   gse <- clusterProfiler::gseKEGG(
           geneList = kegg_gene_list,
           keyType = keyType,
           organism  = kegg_organism,
           nPerm  = 10000,
           minGSSize = minGS,
           maxGSSize = maxGS,
           pvalueCutoff = pvalueCutoff,
           pAdjustMethod = pAdjustMethod
    )
}

fun_MsigDB <- function(msigdb_gene_list, keyType, msigdb, minGS, maxGS, pvalueCutoff, pAdjustMethod ){
  gse <- clusterProfiler::GSEA(
         msigdb_gene_list,
         TERM2GENE = msigdb,
         nPerm = 1000,
         minGSSize = minGS,
         maxGSSize = maxGS,
         pvalueCutoff = pvalueCutoff,
         pAdjustMethod = pAdjustMethod
         )
}

####### main #########
set.seed(1000)
output = snakemake@output[[1]]
# prepare organism
organism = snakemake@params[["organism_OrgDb"]]
kegg_organism = snakemake@params[["organism_KEGG"]]
msigdb_organism = snakemake@params[["organism_MsigDB"]]
if(!(organism %in% installed.packages()[,"Package"]))
    BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

input<-read.table(snakemake@input[["gene_list"]],sep="\t",header=FALSE)
cols=colnames(input)
vals = input[[cols[2]]]
names(vals)=input[[cols[1]]]
if(!is.numeric(vals) | !is.character(names(vals)))
  stop("Error found in gene list. Please make sure gene list is correctly formatted." )
gene_list <- fun_prep_input(vals)

# GO
gse_go <- fun_gseGO(
    gene_list = gene_list,
    ont=snakemake@params[["GO_ontology"]],
    keyType=snakemake@params[["GO_keyType"]],
    minGS=snakemake@params[["GO_min_geneset_size"]],
    maxGS=snakemake@params[["GO_max_geneset_size"]],
    pvalueCutoff = snakemake@params[["GO_pval_cutoff"]],
    orgdb = organism,
    pAdjMethod=snakemake@params[["GO_padj_method"]]
)


# KEGG
keyType = snakemake@params[["KEGG_keyType"]]
if(!(tolower(keyType) %in% c('kegg', 'ncbi-geneid','entrezid', 'ncib-proteinid', 'uniprot'))){
    kegg_gene_list <- fun_convert_geneid(
        gene_list = gene_list,
        fromType = toupper(keyType),
        toType = "ENTREZID",
        orgdb=organism
    )
    keyType = "ncbi-geneid"
}else{
  kegg_gene_list = gene_list
}

gse_kegg = fun_gseKEGG(
  kegg_gene_list = kegg_gene_list,
  keyType = tolower(keyType),
  kegg_organism = kegg_organism,
  minGS = snakemake@params[["KEGG_min_geneset_size"]],
  maxGS = snakemake@params[["KEGG_max_geneset_size"]],
  pvalueCutoff = snakemake@params[["KEGG_pval_cutoff"]],
  pAdjustMethod = snakemake@params[["KEGG_padj_method"]]
)

# MsigDB
keyType = toupper(snakemake@params[["MsigDB_keyType"]])
if(!(keyType %in% c('ENTREZID', 'SYMBOL'))){
  msigdb_gene_list <- fun_convert_geneid(
    gene_list = gene_list,
    fromType = keyType,
    toType = "SYMBOL",
    orgdb=organism
  )
  keyType = "SYMBOL"
}else{
  msigdb_gene_list = gene_list
}

cols = c("gs_name","gene_symbol")
if(keyType == "ENTREZID"){
  cols = c("gs_name","entrez_gene")
}
db=msigdbr::msigdbr(species = msigdb_organism,
           category = snakemake@params[["MsigDB_category"]],
           subcategory = snakemake@params[["MsigDB_subcategory"]])[,cols]

gse_msigdb = fun_MsigDB(
  msigdb_gene_list = msigdb_gene_list,
  keyType = keyType,
  msigdb = db,
  minGS = snakemake@params[["MsigDB_min_geneset_size"]],
  maxGS = snakemake@params[["MsigDB_max_geneset_size"]],
  pvalueCutoff = snakemake@params[["MsigDB_pval_cutoff"]],
  pAdjustMethod = snakemake@params[["MsigDB_padj_method"]]
)


gses = list("GO"=gse_go, "KEGG" = gse_kegg, "MsigDB" = gse_msigdb)
genes = list("GO"=gene_list, "KEGG" = kegg_gene_list, "MsigDB" = msigdb_gene_list)

idx = sapply(gses,nrow)>0
gses <- gses[idx]
genes <- genes[idx]


saveRDS(gses, gsub(".pdf",".Rds",output))

p1 <- lapply(gses, function(gse){
           clusterProfiler::dotplot(gse,
           showCategory=snakemake@params[['dotplot_showCategory']],
           color=snakemake@params[['dotplot_color']],
           font.size=snakemake@params[['dotplot_fontSize']],
           title=snakemake@params[['dotplot_title']],
           split=".sign") + facet_grid(.~.sign)})

p2 <- lapply(gses, function(gse){
           clusterProfiler::emapplot(gse,
           showCategory = snakemake@params[['emaplot_showCategory']],
           color = snakemake@params[['emaplot_color']],
           layout = snakemake@params[['emaplot_layout']],
           pie_scale = snakemake@params[['emaplot_pie_scale']],
           line_scale = snakemake@params[['emaplot_line_scale']])})

p3=list()
for(i in 1:length(gses)) {
  p3[[names(gses)[i]]] = clusterProfiler::cnetplot(gses[[i]],
             foldChange=genes[[names(gses)[i]]],
             showCategory = snakemake@params[['cnetplot_showCategory']],
             layout = snakemake@params[['cnetplot_layout']],
             colorEdge = snakemake@params[['cnetplot_colorEdge']],
             circular = snakemake@params[['cnetplot_circular']],
             node_label = snakemake@params[['cnetplot_node_label']])}

p4 <- lapply(gses, function(gse){
           clusterProfiler::ridgeplot(gse,
           showCategory = snakemake@params[['redgeplot_showCategory']],
           fill = snakemake@params[['redgeplot_fill']],
           core_enrichment = snakemake@params[['redgeplot_core_enrichment']])+theme(axis.text.y=element_text(size=8))
           })

p5 <- lapply(gses, function(gse){
           enrichplot::gseaplot(gse,
           title = snakemake@params[['gseaplot_title']],
           geneSetID = snakemake@params[['gseaplot_geneSetID']])})

if("KEGG" %in% names(gses)){
    cwd = getwd()
    outdir = dirname(output)
    setwd(outdir)
    pathway.id = as.character(unlist(strsplit(snakemake@params[['pathview_pathway_id']],",")))
    pathway.id = gsub("\\s+","",pathway.id)
    lapply(pathway.id, function(pid){
              pathview(gene.data=genes[['KEGG']],
              pathway.id=pid,
              species = kegg_organism)})
    setwd(cwd)
}
pdf(output,height=8, width=11)
for(i in 1:length(gses)){
  print(p1[[i]]);print(p2[[i]]);print(p3[[i]]);print(p4[[i]]);print(p5[[i]])
}
dev.off()
message("Done.")
