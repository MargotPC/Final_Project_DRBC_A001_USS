# Cargar librerías necesarias
library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)  # Para etiquetas que no se sobrepongan
library(biomaRt)
# Leer el archivo CSV
data <- read_csv("C:/Users/margo/Documents/USS/bioinformatica/trabajo_final/RNA-Seq-expression-Norilsk2019.csv")

# Reemplazar padj == 0 por 1e-300
data$`Adjusted p-value`[data$`Adjusted p-value` == 0] <- 1e-300

# Filtrar y preparar datos
data <- data %>%
  filter(
    !is.na(`log_2 fold change`),
    is.finite(`log_2 fold change`),
    !is.na(`Adjusted p-value`),
    is.finite(`Adjusted p-value`),
    `Adjusted p-value` > 0
  ) %>%
  mutate(
    neg_log10_padj = -log10(`Adjusted p-value`),
    regulation = case_when(
      `log_2 fold change` > 2 & `Adjusted p-value` < 0.01 ~ "Sobreexpresados",
      `log_2 fold change` < -2 & `Adjusted p-value` < 0.01 ~ "Subexpresados",
      TRUE ~ "No significativos"
    )
  )
# Conectar a Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Mapear los nombres
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = data$Gene,
  mart = ensembl
)

# Combinar los nombres cortos
data <- left_join(data, gene_map, by = c("Gene" = "ensembl_gene_id"))

# Seleccionar genes a etiquetar (top 5 de cada grupo significativo)
genes_etiquetados <- data %>%
  filter(regulation != "No significativos") %>%
  arrange(`Adjusted p-value`) %>%
  group_by(regulation) %>%
  slice_head(n = 7)

# Crear volcano plot con etiquetas
ggplot(data, aes(x = `log_2 fold change`, y = neg_log10_padj)) +
  geom_point(aes(color = regulation), alpha = 0.7) +
  geom_text_repel(
    data = genes_etiquetados,
    aes(label = hgnc_symbol),
    size = 3,
    max.overlaps = 15,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "gray50"
  ) +
  scale_color_manual(values = c("Sobreexpresados" = "red", "Subexpresados" = "blue", "No significativos" = "gray")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  labs(
    x = expression(Fold~Change~(Log[2])),
    y = expression(-Log[10]~(p-value)),
    color = ""
  ) +
  theme_minimal()
ggsave("volcano_plot.png", width = 8, height = 5, dpi = 300)

