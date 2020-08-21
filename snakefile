configfile: "config.yaml"

rule all:
    input:
        config["output"]

rule enrichment:
    input:
        gene_list=config["input"]["gene_list"],
    output:
        config["output"]
    log:
        "/Users/yuanhao/Analysis/Poirier/workflow/enrichment/enrichment.log"
    #conda:
    #    "/Users/yuanhao/Analysis/Poirier/workflow/preseq/preseqr.yaml"
    params:
        organism_OrgDb=config["organism"]["OrgDb"],
        organism_KEGG=config["organism"]["KEGG"],
        GO_ontology = config["GO"]["ontology"],
        GO_keyType = config["GO"]["key_type"],
        GO_min_geneset_size = config["GO"]["min_geneset_size"],
        GO_max_geneset_size = config["GO"]["max_geneset_size"],
        GO_pval_cutoff = config["GO"]["pvalue_cutoff"],
        GO_padj_method = config["GO"]["padj_method"],
        KEGG_keyType = config["KEGG"]["keyType"],
        KEGG_min_geneset_size = config["KEGG"]["min_geneset_size"],
        KEGG_max_geneset_size = config["KEGG"]["max_geneset_size"],
        KEGG_pval_cutoff = config["KEGG"]["pvalue_cutoff"],
        KEGG_padj_method = config["KEGG"]["padj_method"],
        dotplot_showCategory = config["dotplot"]["showCategory"],
        dotplot_color = config["dotplot"]["color"],
        dotplot_fontSize = config["dotplot"]["fontSize"],
        dotplot_title = config["dotplot"]["title"],
        emaplot_showCategory = config["emaplot"]["showCategory"],
        emaplot_color = config["emaplot"]["color"],
        emaplot_layout = config["emaplot"]["layout"],
        emaplot_pie_scale = config["emaplot"]["pie_scale"],
        emaplot_line_scale = config["emaplot"]["line_scale"],
        cnetplot_showCategory = config["cnetplot"]["showCategory"],
        cnetplot_layout = config["cnetplot"]["layout"],
        cnetplot_colorEdge = config["cnetplot"]["colorEdge"],
        cnetplot_circular = config["cnetplot"]["circular"],
        cnetplot_node_label = config["cnetplot"]["node_label"],
        redgeplot_showCategory = config["redgeplot"]["showCategory"],
        redgeplot_fill = config["redgeplot"]["fill"],
        redgeplot_core_enrichment = config["redgeplot"]["core_enrichment"],
        gseaplot_geneSetID = config["gseaplot"]["geneSetID"],
        gseaplot_title = config["gseaplot"]["title"],
        pathview_pathway_id = config["pathvew"]["pathway_id"]

    script:
        "script/enrichment.R"
