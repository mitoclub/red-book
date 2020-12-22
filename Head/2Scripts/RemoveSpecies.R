library(dplyr)
library(ape)

data = read.csv('../../Body/2Derived/IucnGenbank.csv')
data = data[, -1]

length(unique(data$SpeciesName)) # 956

unique_data = distinct(data)

table(unique_data$category)

# CR  DD  EN  EW  EX  LC  NT  VU 
# 41  50 128   2   5 545  72 113

remove_categories = unique_data %>% 
  mutate(as.character(category)) %>%
  filter(SpeciesName != 'Homo sapiens') %>%
  filter(category %in% c('CR', 'EN', 'LC', 'NT', 'VU'))

table(remove_categories$category)

species_list = paste(sub(' ', '_', remove_categories$SpeciesName), collapse = ',')

species_list

# check why the tree doesn't merge with species list

tree = read.tree('../../Body/1Raw/FromKG/FcC_supermatrix.part.treefile')

species_list

"Tarsius_syrichta" %in% tree$tip.label

View(as.data.frame(tree$tip.label))

not_found = remove_categories %>%
  filter(!SpeciesName %in% c('Tarsius syrichta', 'Peromyscus maniculatus', 
                             'Macropus robustus', 'Sorex gracillimus', 'Viverricula indica',
                             'Tarsius bancanus', 'Castor fiber', 'Myotis bombinus',
                             'Herpestes javanicus'))

paste(sub(' ', '_', not_found$SpeciesName), collapse = ' ')

table(not_found$category)

# CR  DD  EN  EW  EX  LC  NT  VU 
# 41   0 128   0   0 538  70 112

##### checking the pruned tree

pruned_tree = read.tree('../../Body/2Derived/FcC_supermatrix.part.treefile-v')

length(pruned_tree$tip.label)

setdiff(sub(' ', '_', not_found$SpeciesName), pruned_tree$tip.label)

### take species from the tree

not_found$SpeciesName = sub(' ', '_', not_found$SpeciesName)
final = not_found[not_found$SpeciesName %in% pruned_tree$tip.label,]

table(final$category)
# CR  DD  EN  EW  EX  LC  NT  VU 
# 41   0 127   0   0 535  70 111

write.csv(final, '../../Body/2Derived/IucnGenbank_reduced.csv')
