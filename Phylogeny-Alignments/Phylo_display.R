# Load libraries
library(ggtree)
library(ggplot2)
library(ape)
library(phytools)
library(phangorn)
library(tibble)
library(tidytree)
library(treeio)
library(dplyr)
library(extrafont)
loadfonts()

# Read in your tree file (replace "treefile.nwk" with the path to your tree file)
tree <- read.tree("FastTreeMatrix_v1.tree")
tree$tip.label
str(tree)

# Define the outgroup
#outgroup <- ("DQ174731_Chromera_velia")
#outgroup <- "OP390085.1_Stenophora_sp."
#outgroup <- ("KX351970_Myosporidium_spraguei")
#outgroup <- ("MF278562_Mitosporidium_daphniae")
#outgroup <- ("MZ580904_Capitella_capitata")
#outgroup <- ("JX242610_Basidiobolus_microsporus_OG")
outgroup <- ("Ancora_sagittata")

# Reroot the tree
## Setting resolve.root to true adds a node along the branch connecting the root taxon and the rest of the tree. Edgelabel set to true would allow root function to account for correct replacement of node labels.
tree2 <- ape::root(tree, outgroup, edgelabel = TRUE, resolve.root = TRUE)
#tree2 <- midpoint(tree)

# Add root lenght
#tr <- root(tree, outgroup, resolve.root = TRUE, edgelabel = TRUE)
#theroot <- Ntip(tr) + 1
#i <- which(tr$edge[, 1] == theroot & tr$edge.length == 0)

#adjust root length with tr$edge.length[i] <- ..., for instance:
#basal <- which(tr$edge[, 1] == theroot)
#tr$edge.length[i] <- mean(tr$edge.length[basal])
##The method minimally fixes a zero-length root branch

#or modify root branches with
## make the root perfectly midpoint
#tree2$edge.length[which(!(tree2$edge[,1] %in% tree2$edge[,2]))] <- sum(tree2$edge.length[which(!(tree2$edge[,1] %in% tree2$edge[,2]))])/2

# This shortens your tree to fit tip labels. Adjust the factor for a better fit.
xlim_adj <- max(ggtree(tree2)$data$x) * 2.5

# Extend the length of your branches by multiplying the edge lengths by a factor (e.g., 1.5)
#tree$edge.length <- tree$edge.length * 1

# Convert node labels to percentages and filter out values below 50%
tree2$node.label
tree2$node.label <- as.numeric(tree2$node.label) * 100
tree2$node.label <- round(tree2$node.label, 0)
tree2$node.label

# Create a ggtree object
p <- ggtree(tree2, ladderize = TRUE, layout="rectangular")
str(p)

# Modify the tip labels (Fungi)
p$data <- p$data %>%
  mutate(label = case_when(
    label == "id_3_total_supporting_reads_17_Run2_Barcode1_Pieris" ~ "Run 2 Barcode1 Pieris brassicae (17) UnonMet",
    label == "id_0_total_supporting_reads_15_Run2_Barcode1_Pieris" ~ "Run 2 Barcode1 Pieris brassicae (15) UnonMet",
    TRUE ~ label
  ))

# Modify the tip labels (Microsporidia)
p$data <- p$data %>%
  mutate(label = case_when(
    label == "id_5_total_supporting_reads_205_Run3_Barcode7_Gryllus_assimilis_V1f_530R" ~ "Run3 Barcode 7 Gryllus assimilis DK V1F/530R (205)",
    label == "id_33_total_supporting_reads_172_Run3_Barcode6_Gryllus_bimaculatus_V1f_530R" ~ "Run3 Barcode 6 Gryllus bimaculatus V1F/530R (172)",
    label == "MG663123.1_Pleistophoridae_sp_YST-2017a" ~ "MG663123 Albopleistophora grylli",
    label == "id_84_total_supporting_reads_80_Run3_Barcode5_Gryllus_bimaculatus_V1f_530R" ~ "Run3 Barcode 5 Gryllus bimaculatus DK V1F/530R (80)",
    label == "id_51_total_supporting_reads_81_Run3_Barcode6_Gryllus_bimaculatus_V1f_530R" ~ "Run3 Barcode 6 Gryllus bimaculatus V1F/530R (81)",
    label == "id_9_total_supp_69_Run3_Barcode5_Gryllus_bimaculatus_V1f_530R" ~ "Run3 Barcode 5 Gryllus bimaculatus DK V1F/530R (69)",
    label == "id_84_total_supp_reads_80_Run3_Barcode5_Gryllus_bimaculatus_V1f_530R" ~ "Run 3 Barcode 5 Gryllus bimaculatus DK V1F/530R (80)",
    TRUE ~ label
  ))

# Modify the tip labels (gregarines 18S)
p$data <- p$data %>%
  mutate(label = case_when(
    label == "id_1_total_supporting_reads_11_Run1_Barcode1_Aphomia_sociella_GregF_GregR" ~ "Run1 Barcode1 Aphomia sociella GregF-R (11)",
    label == "id_6_total_supporting_reads_18_Run1_Barcode8_Gryllus_bimaculatus_UnonMet" ~ "Run1 Barcode8 Gryllus bimaculatus UnonMet (18)",
    label == "id_3_total_supporting_reads_130_Run3_Barcode1_Acheta_domesticus_gregarine_UnonMet" ~ "Run3 Barcode1 Acheta domesticus gregarine UnonMet (130)",
    label == "id_0_total_supporting_reads_279_Run3_Barcode1_Acheta_domesticus_gregarine_WL1_EukP3" ~ "Run3 Barcode1 Acheta domesticus gregarine WL1 EukP3 (279)",
    label == "id_12_total_supporting_reads_29_Run1_Barcode10_Gryllus_assimilis_UnonMet" ~ "Run1 Barcode10 Gryllus assimilis UnonMet (29)",
    label == "id_28_total_supporting_reads_195_Run3_Barcode2_Blaptica_dubia_GregF_GregR" ~ "Run3 Barcode2 Blaptica dubia GregF-R (195)",
    label == "id_4_total_supporting_reads_100_Run3_Barcode7_Alphitobius_diaperinus_gregarine_GregF_GregR" ~ "Run3 Barcode7 Alphitobius diaperinus gregarine GregF-R (100)",
    label == "id_17_total_supporting_reads_290_Run3_Barcode3_Tenebrio_molitor_GregF_GregR" ~ "Run3 Barcode3 Tenebrio molitor GregF-R (290)",
    label == "id_1_total_supporting_reads_24_Run1_Barcode5_Gryllus_assimilis_GregF_GregR" ~ "Run1 Barcode5-2-7 Gryllus assimilis GregF-R (24)",
    label == "id_64_total_supporting_reads_81_Run3_Barcode6_Gryllus_bimaculatus_GregF_GregR" ~ "Run3 Barcode6 Gryllus bimaculatus GregF-R (81)",
    label == "id_0_total_supporting_reads_114_Run3_Barcode5_Gryllus_bimaculatus_GregF_GregR" ~ "Run3 Barcode5 Gryllus bimaculatus DK GregF-R (114)",
    label == "id_916_total_supporting_reads_99_Run3_Barcode7_Gryllus_assimilis_GregF_GregR" ~ "Run3 Barcode7 Gryllus assimilis DK GregF-R (99)",
    label == "id_10_total_supporting_reads_235_Run3_Barcode1_Acheta_domesticus_gregarine_GregF_GregR" ~ "Run3 Barcode1 Acheta domesticus gregarine GregF-R (235)",
    label == "id_14_total_supporting_reads_284_Run3_Barcode4_Acheta_domesticus_GregF_GregR" ~ "Run3 Barcode4 Acheta domesticus GregF-R (284)",
    label == "id_2087_total_supporting_reads_45_Run3_Barcode8_Gryllus_assimilis_GregF_GregR" ~ "Run3 Barcode8 Gryllus assimilis GregF-R (45)",
    label == "id_21_total_supporting_reads_2_Run1_Barcode5_Gryllus_assimilis_KinetF" ~ "Run1 Barcode5 Gryllus assimilis KintetF-R (2)",
    label == "id_3_total_supporting_reads_2_Run1_Barcode1_Aphomia_sociella_Kinet" ~ "Run1 Barcode1 Aphomia sociella KintetF-R (2)",
    label == "id_1_total_supporting_reads_3_Run1_Barcode2_Gryllus_bimaculatus_Kinet" ~ "Run1 Barcode2 Gryllus bimaculatus KintetF-R (3)",
    label == "id_2_total_supporting_reads_2_Run1_Barcode7_Gryllus_sigillatus_Kinet" ~ "Run1 Barcode7 Gryllus sigillatus KintetF-R (2)",
    TRUE ~ label
  ))

# Modify the tip labels (gregarines 28S)
p$data <- p$data %>%
  mutate(label = case_when(
    label == "id_0_total_supporting_reads_4_Run2_Barcode3_Gryllus_bimaculatus_GregRC_iLSUr" ~ "Run2 Barcode3 Gryllus bimaculatus GregRC iLSUr (4)",
    label == "id_1_total_supporting_reads_3_Run1_Barcode5_Gryllus_assimilis_KinetRC_iLSUr" ~ "Run1 Barcode5 Gryllus assimilis KinetRC_iLSUr (3)",
    label == "id_0_total_supporting_reads_3_Run1_Barcode7_Gryllus_sigillatus_KinetRC_iLSUr" ~ "Run1_Barcode7_Gryllus_sigillatus_KinetRC_iLSUr (3)",
    label == "id_0_total_supporting_reads_6_Run2_Barcode4_Gryllus_assimilis_GregRC_iLSUr" ~ "Run2_Barcode4 & Run1_Barcode5_Gryllus_assimilis_GregRC_iLSUr (6)",
    label == "id_0_total_supporting_reads_4_Run1_Barcode2_Gryllus_bimaculatus_GregRC_iLSUr" ~ "Run1_Barcode2_Gryllus_bimaculatus_GregRC_iLSUr (4)",
    label == "id_1_total_supporting_reads_8_Run1_Barcode5_Gryllus_assimilis_GregRC_iLSUr" ~ "Run1_Barcode5_Gryllus_assimilis_GregRC_iLSUr (8)",
    label == "JF412715_2041_4815_Gregarina_sp__Phaedon_Daegwallyeong_2010" ~ "JF412715_Gregarina_sp__Phaedon",
    label == "id_0_total_supporting_reads_2_Run1_Barcode3_Galleria_mellonella_KinetRC_iLSUr" ~ "Run1 Barcode3 Galleria mellonella KinetRC/iLSUr (2)",
    label == "id_0_total_supporting_reads_7_Run1_Barcode7_Gryllus_sigillatus_GregRC_iLSUr" ~ "Run1 Barcode7-5 Gryllus sigillatus & G. assimilis GregRC/iLSUr (7)",
    label == "id_0_total_supporting_reads_2_Run1_Barcode5_Gryllus_assimilis_GregRC_iLSUr" ~ "Run1 Barcode5 Gryllus assimilis GregRC/iLSUr (2)",
    label == "id_0_total_supporting_reads_2_Run1_Barcode5_Gryllus_assimilis_KinetRC_iLSUr" ~ "Run1 Barcode5 Gryllus assimilis KinetRC/iLSUr (2)",
    label == "id_1_total_supporting_reads_2_Run1_Barcode1_Aphomia_sociella_GregRC_iLSUr" ~ "Run1 Barcode1 Aphomia sociella GregRC/iLSUr (2)",
    label == "id_5_total_supporting_reads_2_Run1_Barcode7_Gryllus_sigillatus_GregRC_iLSUr" ~ "Run1 Barcode7 Gryllus sigillatus GregRC/iLSUr (2)",
    label == "id_15_total_supporting_reads_2_Run1_Barcode5_Gryllus_assimilis_GregRC_iLSUr" ~ "Run1 Barcode5 Gryllus assimilis GregRC/iLSUr (2)",
    label == "GJRY01020788_TSA_PRJNA784797_Gryllodes_sigillatus_TRINITY_DN5437_c0_g1_i1" ~ "TSA_Gryllodes_sigillatus",
    label == "GHUU01044847_TSA_PRJNA485997_Acheta_domesticus_cl_883774_16_RC" ~ "TSA_Acheta_domesticus",
    label == "" ~ "",
    TRUE ~ label
  ))

p$data
str(p)

#p <- ggtree(tree2, ladderize = TRUE, layout="rectangular", size = 0.3)
# Plot the tree with new labels
p <- p + 
  geom_tiplab(aes(label = label), hjust = 0, size = 4, linesize = 0.1, offset = 0.001, fontface = "italic", family = "Times New Roman") + 
  geom_treescale(y = -0.95, fontsize = 3.9) +
  geom_text2(aes(label = round(as.numeric(label), 2), 
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 0 & as.numeric(label) <= 100), 
             vjust = -0.5, hjust = 1.2, size = 3, check_overlap = TRUE) + 
  #geom_text(aes(label=ifelse(!isTip, node, '')), hjust=-0.5, vjust=0.5, size=3) + # see the nodes number
  theme(legend.text = element_text(size = 8)) + 
  xlim(0, xlim_adj) +
  scale_fill_identity(guide = "none")

# Display the tree
p

# Rotate for COI
p <- ggtree::rotate(p, 25)
#p <- ggtree::rotate(p, 13)
p

ggsave("Greg_18S.png", p, dpi = 900, width = 6.5, height = 4.5, units = "in")
