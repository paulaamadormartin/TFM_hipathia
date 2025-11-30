if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)


pathways_signor_mmu=readRDS("/Users/paula/Documents/Bioinformatica/tfm/25_26/signor2Hipathia/tmp/hipathia_mmu/2025_30_nov._dom._10_54Vhi_2.13.1_mgi.RDS")
pathways_signor_mmu$species
subG_signor=pathways_signor_mmu$pathigraphs$mmuSIGXOR0AC$effector.subgraphs$`P-mmuSIGXOR0AC-14`
V(subG_signor)[degree(subG_signor, mode = "out")==0]$genesList
igraph::as_data_frame(subG_signor, what = "vertices") %>% View
get_entrez_effector=function(subG_signor){
  V(subG_signor)[degree(subG_signor, mode = "out")==0]$genesList %>% unlist()  
}
get_entrez_effector(pathways_signor_mmu$pathigraphs$mmuSIGXOR0AC$effector.subgraphs$`P-mmuSIGXOR0AC-14`)

get_effector_gene_pathway_map = function(pathways_signor_mmu) {
  result = list()
  
  for (pw_id in names(pathways_signor_mmu$pathigraphs)) {
    pw = pathways_signor_mmu$pathigraphs[[pw_id]]
    pw_name = pw$path.name
    
    for (subgraph_name in names(pw$effector.subgraphs)) {
      subG = pw$effector.subgraphs[[subgraph_name]]
      genes = get_entrez_effector(subG)
      
      if (length(genes) > 0) {
        df = data.frame(
          gene = genes,
          pathway_id = pw_id,
          pathway_name = pw_name,
          effector_circuit = subgraph_name,
          stringsAsFactors = FALSE
        )
        result[[paste(pw_id, subgraph_name, sep = "_")]] <- df
      }
    }
  }
  
  final_df = do.call(rbind, result)
  return(final_df)
}
effector_gene_pathway_map=get_effector_gene_pathway_map(pathways_signor_mmu)

#los entrez id están asociados a humanos, hacemos una conversión para pasarlos a ratón
entrez_ids = unique(effector_gene_pathway_map$gene)
entrez_ids = as.character(entrez_ids)  # por si acaso

mart_hs <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mart_mm <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")


human_basic <- getBM(
  attributes = c(
    "entrezgene_id",   # Entrez humano
    "hgnc_symbol",     # símbolo humano
    "ensembl_gene_id"  # ID Ensembl humano
  ),
  filters  = "entrezgene_id",
  values   = entrez_ids,
  mart     = mart_hs
)

head(human_basic)

human_homologs <- getBM(
  attributes = c(
    "ensembl_gene_id",              # Ensembl humano (para unir)
    "mmusculus_homolog_ensembl_gene" # Ensembl ortólogo de ratón
  ),
  filters  = "ensembl_gene_id",
  values   = unique(human_basic$ensembl_gene_id),
  mart     = mart_hs
)

head(human_homologs)

human_full <- human_basic %>%
  inner_join(human_homologs, by = "ensembl_gene_id") %>%
  filter(mmusculus_homolog_ensembl_gene != "")

head(human_full)

mouse_annot <- getBM(
  attributes = c(
    "ensembl_gene_id",  # Ensembl ratón
    "entrezgene_id",    # Entrez ratón
    "mgi_symbol"        # símbolo ratón
  ),
  filters  = "ensembl_gene_id",
  values   = unique(human_full$mmusculus_homolog_ensembl_gene),
  mart     = mart_mm
)

colnames(mouse_annot) <- c("ensembl_gene_mm", "entrez_mm", "symbol_mm")

head(mouse_annot)

mapping_hs_mm <- human_full %>%
  inner_join(
    mouse_annot,
    by = c("mmusculus_homolog_ensembl_gene" = "ensembl_gene_mm")
  )

colnames(mapping_hs_mm)[1:4] <- c(
  "entrez_hs", "symbol_hs", "ensembl_hs", "ensembl_mm"
)

head(mapping_hs_mm)

mapping_hs_mm$entrez_hs <- as.character(mapping_hs_mm$entrez_hs)
effector_gene_pathway_map$gene <- as.character(effector_gene_pathway_map$gene)

effector_gene_pathway_map_mmu_filtered = effector_gene_pathway_map_mmu %>%
  filter(!is.na(entrez_mm))

effector_gene_pathway_csv = data.frame(
  lapply(effector_gene_pathway_map_mmu_filtered, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) paste(x, collapse = ";"))
    } else {
      col
    }
  }),
  stringsAsFactors = FALSE
)
head(effector_gene_pathway_csv)
write.csv(effector_gene_pathway_csv, file = "/Users/paula/Documents/Bioinformatica/tfm/25_26/SIGNOR_circuits_annotated_mmu.csv", row.names = FALSE)
