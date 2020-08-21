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
        "enrichment.log"
    params:
        organism_OrgDb=config["organism"]["OrgDb"],
        organism_KEGG=config["organism"]["KEGG"],
        organism_MsigDB=config["organism"]["MsigDB"],
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
        MsigDB_keyType = config["MsigDB"]["keyType"],
        MsigDB_category = config["MsigDB"]["category"],
        MsigDB_subcategory = config["MsigDB"]["subcategory"],
        MsigDB_min_geneset_size = config["MsigDB"]["min_geneset_size"],
        MsigDB_max_geneset_size = config["MsigDB"]["max_geneset_size"],
        MsigDB_pval_cutoff = config["MsigDB"]["pvalue_cutoff"],
        MsigDB_padj_method = config["MsigDB"]["padj_method"],
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
