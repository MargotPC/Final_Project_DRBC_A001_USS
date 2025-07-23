# ─── CONFIGURACIÓN ──────────────────────────────────────────────────────────────
# Instalar paquetes necesarios (solo una vez)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "readr", "dplyr"))

# Cargar librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(readr)
library(dplyr)
library(ggplot2)

# ─── CARGAR DATOS ────────────────────────────────────────────────────────────────
# Leer archivo procesado con la columna 'regulation'
data <- read_csv("C:/Users/margo/Documents/USS/bioinformatica/trabajo_final/RNA-Seq-expression-Norilsk2019.csv")

# Reemplazar padj == 0 por 1e-300
data$`Adjusted p-value`[data$`Adjusted p-value` == 0] <- 1e-300

# Clasificar según regulación
data <- data %>%
  filter(
    !is.na(`log_2 fold change`),
    is.finite(`log_2 fold change`),
    !is.na(`Adjusted p-value`),
    is.finite(`Adjusted p-value`),
    `Adjusted p-value` > 0
  ) %>%
  mutate(
    regulation = case_when(
      `log_2 fold change` > 2 & `Adjusted p-value` < 0.01 ~ "Upregulated",
      `log_2 fold change` < -2 & `Adjusted p-value` < 0.01 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# ─── EXTRAER GENES SOBREEXPRESADOS ───────────────────────────────────────────────
genes_up <- data %>%
  filter(regulation == "Upregulated") %>%
  pull(Gene) %>%
  unique()

# ─── CONVERTIR A ENTREZ ID ───────────────────────────────────────────────────────
gene_map_up <- bitr(genes_up,
                    fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db)

# ─── ANÁLISIS GO ─────────────────────────────────────────────────────────────────
go_up <- enrichGO(gene          = gene_map_up$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "ALL",  # También "BP", "MF", "CC"
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

# ─── GUARDAR RESULTADOS ──────────────────────────────────────────────────────────
# Exportar tabla como CSV
go_df <- as.data.frame(go_up)
write.csv(go_df, "GO_upregulated_results.csv", row.names = FALSE)

# ─── GRAFICAR ────────────────────────────────────────────────────────────────────
# Barplot de los 10 términos más representativos
#barplot(go_up, showCategory = 10)

# También puedes usar:
dotplot(go_up,split ="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + theme(plot.margin = margin(5, 5, 5, 25))
# Crear Bubble Plot con personalización extra
ggsave("go_new1.png", width = 10, height = 20, dpi = 300)

# Cruzar data con SYMBOL y ENTREZID
fc_df <- data %>%
  inner_join(gene_map_up, by = c("Gene" = "ENSEMBL")) %>%
  select(ENTREZID, SYMBOL, `log_2 fold change`) %>%
  filter(!is.na(ENTREZID) & !is.na(SYMBOL)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

# Crear vector con nombres = SYMBOL (para etiquetas) y valores = log2FC
fc_vector_named <- setNames(fc_df$`log_2 fold change`, fc_df$SYMBOL)

# Asignar los símbolos directamente al objeto enriquecido
go_up_readable <- setReadable(go_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")

# Graficar heatplot
png("heatmap_GOnew.png", width = 1600, height = 1000, res = 150)
heatplot(go_up_readable, showCategory = 10, foldChange = fc_vector_named)
dev.off()


