library(tidyverse)
library(readr)
library(Kmisc)
library(readxl)

#------------------------- 1.read data -------------------------------------------

geneCount <- read.table("Orthogroups.GeneCount.tsv", header=T)		# orthofinder output file.
species <- read_excel("SpeciesInfo.xlsx")		
# The file includes species names(in order on the species tree), abbreviations(corresponds to the orthofinder output file), and information belonging to branch, order, and family.


#------------------------- 2.view data -------------------------------------------
# view the distribution of the data by scatter plot.
ggplot(geneCount,aes(x=Orthogroup,y=Total))+
  geom_point(position_jitter(width=0.1,height=0),alpha=0.5,size =0.5) +
  theme_bw()+theme(panel.grid.major.x = element_blank(),axis.title.x = element_blank()) +
  ggtitle("GeneCount") +
  ylab("gene count in each Orthogroup") +
  ggsave("Orthogroups_GeneCount_scatter_plot.pdf", device = "pdf", dpi = 600, width = 2, height = 2.5)

# view the distribution of the data by line plot.
ggplot(geneCount,aes(x=Orthogroup,y=Total))+
  geom_line(group = 1) +
  geom_vline(xintercept = "OG0000099", colour="#BB0000", linetype="dashed") +
  geom_vline(xintercept = "OG0000199", colour="#8d4bbb", linetype="dashed") +
  geom_vline(xintercept = "OG0000299", colour="#1ee3cf", linetype="dashed") +
  theme_bw() + theme(panel.grid.major.x = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  annotate(geom = "text", x = "OG0000450", y = 3000, colour="#BB0000", label = "top 100 (gene count = 134)") +
  annotate(geom = "text", x = "OG0000550", y = 2000, colour="#8d4bbb", label = "top 200 (gene count = 42)") +
  annotate(geom = "text", x = "OG0000650", y = 1000, colour="#1ee3cf", label = "top 300 (gene count = 23)") +
  # annotate x,y is where the comment tags.
  ggtitle("GeneCount") +
  ylab("gene count in each Orthogroup") +
  ggsave("Orthogroups_GeneCount_line_plot.pdf", device = "pdf", dpi = 600, width = 6, height = 2.5)


#------------------------- 3.transform the data to analyse. --------------------------------------

geneCount_clean <- geneCount %>%
  # select(-Total) %>%
  column_to_rownames("Orthogroup") %>%		# column_to_rownames(): converts a column to a row name.
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Abbreviation") %>%		# rownames_to_column(): converts a row name to a column.
  left_join(species) %>%		# join the species information data.
  gather(starts_with("OG"), key = Orthogroup, value = geneCount) %>%
  write_tsv(file = "Orthogroups_GeneCount_clean.xls")


#------------------------- 4. highly conserved gene family --------------------------------------

# "highly conserved" families are those present in more than 80% of species in each evolutionary branch.
# get the branch level data.
geneCount_clean.branch <- geneCount_clean %>%
  select(branch, Orthogroup, geneCount) %>%
  group_by(branch, Orthogroup) %>%
  summarise(geneCount = sum(as.numeric(geneCount))) %>%
  spread(key = Orthogroup, value = geneCount) %>%
  column_to_rownames("branch") %>%
  t() %>%
  as.data.frame() %>%
  select(-Total, Total)

# Count the proportion of species that contain the corresponding Orthogroup in each branch.
# Orthogroups will be retained only if they are present in more than 0.8 species in a branch.
Orthogroups_proportion.branch <- geneCount_clean %>%
  select(Abbreviation, branch, Orthogroup, geneCount) %>%
  mutate(geneExistence = case_when(geneCount == 0 ~ 0, geneCount >= 1 ~ 1)) %>%
  group_by(Orthogroup, branch) %>%
  summarise(speciesCount = n(), geneExistenceCount = sum(geneExistence)) %>%
  mutate(proportion = geneExistenceCount/speciesCount) %>%
  filter(proportion >= 0.8) %>%
  filter(branch != "Total")

### highly conserved gene family
highlyConserved <-
  geneCount_clean.branch %>%
  filter_all(all_vars(. > 0)) %>%
  rownames_to_column("Orthogroup") %>%
  select(Orthogroup) %>%
  mutate(Type = "highlyConserved") %>%
  semi_join(Orthogroups_proportion.branch) %>%
  write_tsv(file = "highlyConservedGeneFamily.xls")


#------------------------- 5. species-specific gene family --------------------------------------

# "species-specific" families are those present only in a particular species not in other species;
# get the species level data
geneCount_clean.species <- geneCount_clean %>%
  select(Abbreviation, Orthogroup, geneCount) %>%
  spread(key = Orthogroup, value = geneCount) %>%
  column_to_rownames("Abbreviation") %>%
  t() %>%
  as.data.frame() %>%
  select(-Total, Total)

### speciesSepcificFamily
speciesSepcificFamily <-
  geneCount_clean.species %>%
  # rownames_to_column("Orthogroup") %>%
  # rowwise() %>%
  # mutate(Total = sum(across(where(is.numeric)))) %>%
  # column_to_rownames("Orthogroup") %>%
  mutate_all(~(./Total)) %>%
  select(-Total) %>%
  filter_all(any_vars(. == 1)) %>%
  rownames_to_column("Orthogroup") %>%
  select(Orthogroup) %>%
  mutate(Type = "speciesSepcific") %>%
  write_tsv(file = "speciesSepcificFamily_geneId.xls")

## Convert "proportion" to "gene count"
# speciesSepcificOrthogroups_count <-
#   speciesSepcificOrthogroups_proportion %>%
#   select(Orthogroup) %>%
#   left_join(geneCount) %>%
#   write_tsv(file = "speciesSepcificOrthogroups_count.xls")


#------------------------- 5. Summary output --------------------------------------

Orthogroups_Type <- rbind(branchConservedOrthogroups,speciesSepcificOrthogroups)

geneCount %>% select(Orthogroup) %>%
  left_join(Orthogroups_Type) %>%
  group_by(Type) %>%
  summarise(Orthogroup = n()) %>%
  write_tsv(file = "Orthogroups_Type.xls", na = "other")


#------------------------- 6. Pick out the representatives gene family and plot --------------------------------------

data <- read_excel("subfamily_number.xlsx",sheet=1)

used_data <- data %>% as.data.frame() %>% select(Abbreviation,Chitin_bind_1,Nodulin_late)

used_data$Abbreviation <- factor(used_data$Abbreviation,levels=used_data[,1])

plot_chitin <- ggplot(used_data,aes(Abbreviation,Chitin_bind_1,group=1))+
  geom_line(color="#f8766d")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
 
plot_nodulin <- ggplot(used_data,aes(Abbreviation,Nodulin_late,group=1))+
  geom_line(color="#00bfc4")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  
pdf("Chitin_bind_1-line.pdf",width = 40,height = 24)
print(plot_chitin)
dev.off()

pdf("Nodulin_late-line.pdf",width = 40,height = 24)
print(plot_nodulin)
dev.off()

tr_data <- gather(used_data,key="family",value="GeneCount",Chitin_bind_1,Nodulin_late)

tr_plot <- ggplot(tr_data,aes(x=Abbreviation,y=GeneCount,group=family,color=family))+
  geom_line()+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
pdf("line-Chitin-Nodulin.pdf",width = 40,height = 24)
print(tr_plot)
dev.off()

