library(dplyr)

SRBC_site_macros <- read.csv("data/raw/SRBC_macros.csv")
PA_taxonomy <- read.csv("data/raw/PADEP_macro_taxonomy.csv")
SRBC_macros <- data.frame(SRBC_site_macros)

genus_lookup <- PA_taxonomy %>%
  select(ITIS_PHYLUM, ITIS_CLASS, ITIS_ORDER, ITIS_FAMILY, ITIS_GENUS) %>%
  distinct()
names(genus_lookup)[names(genus_lookup) == "ITIS_PHYLUM"] <- "Phylum"
names(genus_lookup)[names(genus_lookup) == "ITIS_CLASS"] <- "Class"
names(genus_lookup)[names(genus_lookup) == "ITIS_ORDER"] <- "Order"
names(genus_lookup)[names(genus_lookup) == "ITIS_FAMILY"] <- "Family"
names(genus_lookup)[names(genus_lookup) == "ITIS_GENUS"] <- "Genus"

#merged with taxonomy table, but some of the taxa listed weren't actually genera,
#they were a higher taxonomic level
SRBC_macro_tax <- SRBC_macros %>%
 left_join(genus_lookup, by = "Genus")

write.csv(SRBC_macro_tax, "data/processed/SRBC_macro_taxonomy_wPhylum.csv")
