library(dplyr)

iucn = read.csv('../../Body/1Raw/IUCN.csv', sep = ';')
genbank = read.csv('../../Body/2Derived/taxonomy_from_KG_genbank.csv')

table(iucn$class) # MAMMALIA 6295

mammIucn = filter(iucn, class == 'MAMMALIA')
mammIucn$SpeciesName = mammIucn$specieName

# Mustela eversmannii / eversmanii
mammIucn$SpeciesName = sub('eversmanii', 'eversmannii', mammIucn$SpeciesName)
# mammIucn = select(mammIucn, SpeciesName, category)

table(mammIucn$category)
# CR    DD    EN    EW    EX    LC LR/cd LR/lc LR/nt    NT    VU 
# 281   902   607     2    93  3430     1     0     1   375   603 

mammGenbank = filter(genbank, grepl("Eukaryota>Metazoa>Chordata>Craniata>Vertebrata>Euteleostomi>Mammalia>", FullTaxonomy)) %>%
  select(SpeciesName, FullTaxonomy) %>%
  mutate(SpeciesName = trimws(SpeciesName))

withNA = full_join(mammIucn, mammGenbank)

mamm = inner_join(mammIucn, mammGenbank)

# a = merge(mammIucn, mammGenbank, by='SpeciesName')

table(mamm$category)

# CR    DD    EN    EW    EX    LC LR/cd LR/lc LR/nt    NT    VU 
# 41    50   128     2     5   549     0     0     0    73   115

# species from genbank which weren't merged

setdiff(mammGenbank$SpeciesName, mamm$SpeciesName)

# [1] "Tursiops australis"                   "Petaurista hainana"                  
# [3] "Saimiri oerstedii citrinellus"        "Mus musculus domesticus"             
# [5] "Cervus nippon kopschi"                "Cercopithecus lhoesti"               
# [7] "Equus przewalskii"                    "Equus ovodovi"                       
# [9] "Equus burchellii"                     "Chiropotes israelita"                
# [11] "Callicebus lugens"                    "Uropsilus sp. 1 FT-2014"             
# [13] "Tamias sibiricus"                     "Manis temminckii"                    
# [15] "Megaladapis edwardsi"                 "Bison priscus"                       
# [17] "Manis tricuspis"                      "Canis anthus"                        
# [19] "Myotis petax"                         "Uropsilus sp. 4 FT-2015"             
# [21] "Profelis aurata"                      "Puma yagouaroundi"                   
# [23] "Ovis orientalis breed Asian mouflon"  "Ovis vignei breed Urial"             
# [25] "Arctotherium sp."                     "Capricornis sp. YZ-2016"             
# [27] "Aotus azarai"                         "Callithrix pygmaea"                  
# [29] "Cercopithecus albogularis"            "Lagothrix lagotricha"                
# [31] "Cebus xanthosternos"                  "Macaca leucogenys"                   
# [33] "Bison schoetensacki"                  "Petaurista yunanensis"               
# [35] "Homo heidelbergensis"                 "Lepilemur grewcocki"                 
# [37] "Lepilemur wrighti"                    "Lepilemur ahmansoni"                 
# [39] "Lepilemur sahamalazensis"             "Lepilemur jamesi"                    
# [41] "Elephas antiquus"                     "Microcebus tanosi"                   
# [43] "Bos frontalis"                        "Mammut americanum"                   
# [45] "Aonyx cinerea"                        "Manis gigantea"                      
# [47] "Rhinolophus yunnanensis"              "Sundamys annandalei"                 
# [49] "Hsunycteris thomasi"                  "Tamias quadrivittatus"               
# [51] "Tamias rufus"                         "Tamias canipes"                      
# [53] "Tamias dorsalis"                      "Tamias cinereicollis"                
# [55] "Tamias umbrinus"                      "Toromys albiventris"                 
# [57] "Mylodon darwinii"                     "Macaca mulatta vestita"              
# [59] "Capreolus pygargus tianschanicus"     "Proechimys gularis"                  
# [61] "Proechimys semispinosus calidior"     "Ovis nivicola lydekkeri"             
# [63] "Proechimys poliopus"                  "Budorcas taxicolor tibetana"         
# [65] "Neodon fuscus"                        "Homo sapiens neanderthalensis"       
# [67] "Nycticebus coucang insularis"         "Megalonyx jeffersonii"               
# [69] "Megatherium americanum"               "Pantherina griselda"                 
# [71] "Blarinella cf."                       "Acratocnus ye"                       
# [73] "Nothrotheriops shastensis"            "Parocnus serus"                      
# [75] "Myotragus balearicus"                 "Arctictis binturong albifrons"       
# [77] "Budorcas taxicolor taxicolor"         "Miniopterus fuliginosus"             
# [79] "Bootherium bombifrons"                "Rucervus duvaucelii branderi"        
# [81] "Aotus azarai azarai"                  "Monachus schauinslandi"              
# [83] "Phoca fasciata"                       "Phoca groenlandica"                  
# [85] "Cervus nippon taiouanus"              "Muntiacus reevesi micrurus"          
# [87] "Ursus thibetanus mupinensis"          "Mammuthus primigenius"               
# [89] "Anomalurus sp. GP-2005"               "Ursus thibetanus formosanus"         
# [91] "Bos grunniens"                        "Ailurus fulgens styani"              
# [93] "Canis lupus lupus"                    "Camelus bactrianus"                  
# [95] "Camelus dromedarius"                  "Mus musculus musculus"               
# [97] "Canis lupus chanco"                   "Canis lupus familiaris"              
# [99] "Uncia uncia"                          "Ursus spelaeus"                      
# [101] "Arctodus simus"                       "Ursus thibetanus ussuricus"          
# [103] "Ursus thibetanus thibetanus"          "Gorilla gorilla gorilla"             
# [105] "Rusa unicolor swinhoei"               "Manis tetradactyla"                  
# [107] "Sus scrofa domesticus"                "Giraffa camelopardalis angolensis"   
# [109] "Lama glama"                           "Canis lupus laniger"                 
# [111] "Mus musculus castaneus"               "Coelodonta antiquitatis"             
# [113] "Eulemur fulvus fulvus"                "Eulemur fulvus mayottensis"          
# [115] "Eulemur macaco macaco"                "Varecia variegata variegata"         
# [117] "Cervus nippon hortulorum"             "Cervus elaphus xanthopygus"          
# [119] "Cervus elaphus yarkandensis"          "Homo sp. Altai"                      
# [121] "Sus scrofa taiwanensis"               "Rucervus eldi"                       
# [123] "Cervus elaphus songaricus"            "Panthera tigris amoyensis"           
# [125] "Capra hircus"                         "Elephantulus sp. VB001"              
# [127] "Microtus fortis fortis"               "Microtus fortis calamorum"           
# [129] "Mammuthus columbi"                    "Odobenus rosmarus rosmarus"          
# [131] "Physeter catodon"                     "Rhinolophus ferrumequinum korai"     
# [133] "Prionailurus bengalensis euptilurus"  "Apodemus chejuensis"                 
# [135] "Pseudois schaeferi"                   "Przewalskium albirostris"            
# [137] "Plecotus rafinesquii"                 "Platanista minor"                    
# [139] "Hydropotes inermis argyropus"         "Panthera leo persica"                
# [141] "Rhinopithecus bieti 1 RL-2012"        "Rhinopithecus bieti 2 RL-2012"       
# [143] "Pygathrix cinerea 1 RL-2012"          "Pygathrix cinerea 2 RL-2012"         
# [145] "Saimiri boliviensis boliviensis"      "Eospalax baileyi"                    
# [147] "Rhinolophus pumilus"                  "Rhinolophus monoceros"               
# [149] "Bos indicus"                          "Bubalus bubalis"                     
# [151] "Zaglossus bruijni"                    "Equus caballus"                      
# [153] "Felis catus"                          "Equus asinus"                        
# [155] "Ovis aries"                           "Cavia porcellus"                     
# [157] "Bos taurus"                           "Mus musculus molossinus"             
# [159] "Balaenoptera brydei"                  "Cervus nippon yesoensis"             
# [161] "Cervus nippon centralis"              "Cervus nippon yakushimae"            
# [163] "Echymipera rufescens australis"       "Cricetulus griseus"                  
# [165] "Lama pacos"                           "Cervus nippon sichuanicus"           
# [167] "Callicebus donacophilus"              "Rhinolophus ferrumequinum quelpartis"
# [169] "Hexaprotodon liberiensis"             "Gazella erlangeri"                   
# [171] "Mazama gouazoupira"                   "Neotragus moschatus"                 
# [173] "Nannospalax galili"                   "Nannospalax judaei"                  
# [175] "Spalax carmeli"                       "Nannospalax golani"                  
# [177] "Loxodonta cyclotis"                   "Taurotragus derbianus"               
# [179] "Hemitragus jayakari"                  "Eospalax cansus"                     
# [181] "Callicebus cupreus"                 

write.csv(mamm, '../../Body/2Derived/IucnGenbank.csv')
write.csv(withNA, '../../Body/2Derived/IucnGenbankWithNA.csv')
