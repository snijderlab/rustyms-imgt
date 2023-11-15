use serde::{Deserialize, Serialize};
use std::fmt::Display;
#[allow(dead_code)]

macro_rules! species {
    ($($identifier:ident, $common:expr, $imgt:expr)*) => {
        #[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Ord, PartialOrd)]
        pub enum Species {
            $($identifier,)*
        }

        impl Display for Species {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "{}", self.common_name())
            }
        }


        impl Species {
            pub fn common_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $common,)*
                }
            }
            pub fn imgt_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $imgt,)*
                }
            }
            pub fn scientific_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => &$imgt[..$imgt.find(" (").unwrap_or($imgt.len())],)*
                }
            }
            pub fn ident(&self) -> &'static str {
                match self {
                    $(Self::$identifier => stringify!($identifier),)*
                }
            }

            /// Get the species name from IMGT name tag, or None if it is not a proper species
            pub fn from_str(s: &str) -> Result<Option<Self>, ()> {
                match s {
                    $($imgt => Ok(Some(Self::$identifier)),)*
                    "synthetic construct (synthetic construct)" |
                    "unidentified" | "unclassified sequences" | "unidentified cloning vector" |
                    "Cloning vector AbVec-hIgG1" |
                    "Cloning vector AbVec-hIgKappa" |
                    "Cloning vector pASK88-huHRS3-VH-EP3/1" |
                    "Cloning vector pchiIGHG1" |
                    "Cloning vector pchiIGKC" |
                    "Cloning vector pCL" |
                    "Cloning vector pCLZip" |
                    "Cloning vector pMAB136" |
                    "Cloning vector pUR4546" |
                    "Cloning vector pUR4585" |
                    "Expression vector p28BIOH-LIC4" |
                    "Expression vector pFUSE-HEAVY" |
                    "Expression vector pFUSE-hFc2-adapt-scFv" |
                    "Expression vector pFUSE-LIGHT" |
                    "Expression vector pFUSE-mFc2-adapt-scFv" |
                    "Expression vector pFUSE-rFc2-adapt-scFv" |
                    "Expression vector pHIN-PEP" |
                    "Expression vector pHIN-TRI" |
                    "Expression vector pSFV4" |
                    "Expression vector pTH-HIN" |
                    "Phagemid vector pGALD7" |
                    "Phagemid vector pGALD7DL" |
                    "Phagemid vector pGALD7DLFN" |
                    "Phagemid vector pGALD9" |
                    "Phagemid vector pGALD9DL" |
                    "Phagemid vector pGALD9DLFN" |
                    "Phagemid vector pMID21" | "Enterobacteria phage M13 vector DY3F63"
    => Ok(None),
                    _ => Err(()),
                }
            }
        }
    };
}

species!(
    AcanthopagrusSchlegelii, "Black porgy", "Acanthopagrus schlegelii (black porgy)"
    AcipenserBaerii, "Siberian sturgeon", "Acipenser baerii (Siberian sturgeon)"
    AcipenserGueldenstaedtii, "Russian sturgeon", "Acipenser gueldenstaedtii (Russian sturgeon)"
    AcipenserRuthenus, "Sterlet", "Acipenser ruthenus (sterlet)"
    AcipenserSchrenckii, "Amur sturgeon", "Acipenser schrenckii (Amur sturgeon)"
    AcipenserSinensis, "Chinese sturgeon", "Acipenser sinensis (Chinese sturgeon)"
    AiluropodaMelanoleuca, "Giant panda", "Ailuropoda melanoleuca (giant panda)"
    AlligatorSinensis, "Chinese alligator", "Alligator sinensis (Chinese alligator)"
    AmblyrajaGeorgiana, "Antarctic starry skate", "Amblyraja georgiana (Antarctic starry skate)"
    AmblyrajaHyperborea, "Arctic skate", "Amblyraja hyperborea (Arctic skate)"
    AmbystomaMexicanum, "Axolotl", "Ambystoma mexicanum (axolotl)"
    AmeivaAmeiva, "Jungle runners", "Ameiva ameiva"
    AmiaCalva, "Bowfin", "Amia calva (bowfin)"
    AmphiprionClarkii, "Yellowtail clownfish", "Amphiprion clarkii (yellowtail clownfish)"
    AnarhichasMinor, "Spotted wolffish", "Anarhichas minor (spotted wolffish)"
    AnasPlatyrhynchos, "Mallard", "Anas platyrhynchos (mallard)"
    AnguillaAnguilla, "European eel", "Anguilla anguilla (European eel)"
    AnguillaJaponica, "Japanese eel", "Anguilla japonica (Japanese eel)"
    AnolisCarolinensis, "Green anole", "Anolis carolinensis (green anole)"
    AnoplopomaFimbria, "Sablefish", "Anoplopoma fimbria (sablefish)"
    AnserAnser, "Domestic goose", "Anser anser (Domestic goose)"
    AnserAnserDomesticus, "Domestic goose", "Anser anser domesticus"
    AnserCaerulescens, "Snow goose", "Anser caerulescens (Snow goose)"
    AnserSpGIGHV2011, "Geese GIGHV2011", "Anser sp. GIGHV2011"
    AnserSpGIGLV2009, "Geese GIGLV2009", "Anser sp. GIGLV2009"
    AotusAzarai, "Azara's night monkey", "Aotus azarai (Azara's night monkey)"
    AotusNancymaae, "Ma's night monkey", "Aotus nancymaae (Ma's night monkey)"
    AotusTrivirgatus, "Douroucouli", "Aotus trivirgatus (douroucouli)"
    ArgyropelecusHemigymnus, "Half-naked hatchetfish", "Argyropelecus hemigymnus (half-naked hatchetfish)"
    AtelesBelzebuth, "White-bellied spider monkey", "Ateles belzebuth (white-bellied spider monkey)"
    AtelesGeoffroyi, "Black-handed spider monkey", "Ateles geoffroyi (black-handed spider monkey)"
    AtherinaBoyeri, "Big-scale sand smelt", "Atherina boyeri (big-scale sand smelt)"
    BalaenopteraAcutorostrata, "Minke whale", "Balaenoptera acutorostrata (minke whale)"
    BalaenopteraOmurai, "Omura's baleen whale", "Balaenoptera omurai (Omura's baleen whale)"
    BathyrajaAlbomaculata, "White-dotted skate", "Bathyraja albomaculata (white-dotted skate)"
    BathyrajaBrachyurops, "Broadnose skate", "Bathyraja brachyurops (broadnose skate)"
    BathyrajaEatonii, "Eaton's skate", "Bathyraja eatonii (Eaton's skate)"
    BosGaurus, "Gaur", "Bos gaurus (gaur)"
    BosIndicus, "Domestic zebu", "Bos indicus (zebu cattle)"
    BosJavanicus, "Banteng", "Bos javanicus (banteng)"
    BosTaurus, "Domestic bovine", "Bos taurus (bovine)"
    BosTaurusXBosIndicus, "Bos Tauris and Bos indicus cross", "Bos taurus x Bos indicus"
    BovichtusDiacanthus, "Tristan clipfish", "Bovichtus diacanthus"
    BubalusBubalis, "Water buffalo", "Bubalus bubalis (water buffalo)"
    BuergeriaBuergeri, "Buerger's frog", "Buergeria buergeri (Buerger's frog)"
    CaimanCrocodilus, "Spectacled caiman", "Caiman crocodilus (spectacled caiman)"
    CairinaMoschata, "Muscovy duck", "Cairina moschata (Muscovy duck)"
    CallithrixJacchus, "white-tufted-ear marmoset", "Callithrix jacchus (white-tufted-ear marmoset)"
    CallorhinchusMilii, "Elephant shark", "Callorhinchus milii (elephant shark)"
    Camelidae, "Camels", "Camelidae"
    CamelusBactrianus, "Bactrian camel", "Camelus bactrianus (Bactrian camel)"
    CamelusDromedarius, "Arabian camel", "Camelus dromedarius (Arabian camel)"
    CanisLupus, "Gray wolf", "Canis lupus (gray wolf)"
    CanisLupusFamiliaris, "Domestic dog", "Canis lupus familiaris (dog)"
    CanisSp, "Dogs", "Canis sp."
    CapraHircus, "Domestic goat", "Capra hircus (goat)"
    CarassiusAuratus, "Goldfish", "Carassius auratus (goldfish)"
    CarassiusLangsdorfii, "Japanese silver crucian carp", "Carassius langsdorfii (Japanese silver crucian carp)"
    CarcharhinusLeucas, "Bull shark", "Carcharhinus leucas (bull shark)"
    CarcharhinusPlumbeus, "Sandbar shark", "Carcharhinus plumbeus (sandbar shark)"
    CarolliaPerspicillata, "Seba's short-tailed bat", "Carollia perspicillata (Seba's short-tailed bat)"
    CaviaPorcellus, "Domestic guinea pig", "Cavia porcellus (domestic guinea pig)"
    CephalopachusBancanus, "Horsfield's tarsier", "Cephalopachus bancanus (Horsfield's tarsier)"
    CeratotheriumSimum, "White rhinoceros", "Ceratotherium simum (white rhinoceros)"
    CercocebusAtys, "Sooty mangabey", "Cercocebus atys (sooty mangabey)"
    CercocebusTorquatus, "Collared mangabey", "Cercocebus torquatus (collared mangabey)"
    CervusElaphusHispanicus, "Spanish red deer", "Cervus elaphus hispanicus"
    ChaenocephalusAceratus, "Blackfin icefish", "Chaenocephalus aceratus (blackfin icefish)"
    ChampsocephalusEsox, "Pike icefish", "Champsocephalus esox (pike icefish)"
    ChannaArgus, "Northern snakehead", "Channa argus (northern snakehead)"
    ChannaStriata, "Snakehead murrel", "Channa striata (snakehead murrel)"
    ChaunaTorquata, "Southern screamer", "Chauna torquata (southern screamer)"
    ChelonAuratus, "Golden grey mullet", "Chelon auratus (golden grey mullet)"
    ChionodracoHamatus, "Antarctic icefish", "Chionodraco hamatus (Antarctic icefish)"
    ChionodracoRastrospinosus, "Ocellated icefish", "Chionodraco rastrospinosus (ocellated icefish)"
    ChlorocebusAethiops, "Grivet", "Chlorocebus aethiops (grivet)"
    ClupeaPallasii, "Pacific herring", "Clupea pallasii (Pacific herring)"
    ColobusGuereza, "Mantled guereza", "Colobus guereza (mantled guereza)"
    ColobusPolykomos, "King colobus", "Colobus polykomos (king colobus)"
    CricetinaeSp, "Hamster", "Cricetinae gen. sp. (Hamster)"
    CricetulusMigratorius, "Armenian hamster", "Cricetulus migratorius (Armenian hamster)"
    CrocodylusSiamensis, "Siamese crocodile", "Crocodylus siamensis (Siamese crocodile)"
    CtenopharyngodonIdella, "Grass carp", "Ctenopharyngodon idella (grass carp)"
    CygnodracoMawsoni, "Mawson's dragonfish", "Cygnodraco mawsoni (Mawson's dragonfish)"
    CynoglossusSemilaevis, "Tongue sole", "Cynoglossus semilaevis (tongue sole)"
    CynopterusSphinx, "Indian short-nosed fruit bat", "Cynopterus sphinx (Indian short-nosed fruit bat)"
    CyprinusCarpio, "Common carp", "Cyprinus carpio (common carp)"
    DanioRerio, "Zebrafish", "Danio rerio (zebrafish)"
    DaubentoniaMadagascariensis, "Aye-aye", "Daubentonia madagascariensis (aye-aye)"
    DelphinapterusLeucas, "Beluga whale", "Delphinapterus leucas (beluga whale)"
    DelphinusCapensis, "Long-beaked common dolphin", "Delphinus capensis (long-beaked common dolphin)"
    DicentrarchusLabrax, "European seabass", "Dicentrarchus labrax (European seabass)"
    DissostichusMawsoni, "Antarctic toothfish", "Dissostichus mawsoni (Antarctic toothfish)"
    DrosophilaMelanogaster, "Fruit fly", "Drosophila melanogaster (fruit fly)"
    ElapheTaeniura, "Beauty snake", "Elaphe taeniura (beauty snake)"
    EleginopsMaclovinus, "Patagonian blennie", "Eleginops maclovinus (Patagonian blennie)"
    ElopsAaurus, "Ladyfish", "Elops saurus (ladyfish)"
    EpinephelusAkaara, "Hong Kong grouper", "Epinephelus akaara (Hong Kong grouper)"
    EpinephelusCoioides, "Orange-spotted grouper", "Epinephelus coioides (orange-spotted grouper)"
    EptatretusBurgeri, "Inshore hagfish", "Eptatretus burgeri (inshore hagfish)"
    EptesicusFuscus, "Big brown bat", "Eptesicus fuscus (big brown bat)"
    EquusAsinus, "Ass", "Equus asinus (ass)"
    EquusBurchelliiAntiquorum, "Burchell's zebra", "Equus burchellii antiquorum"
    EquusCaballus, "Domestic horse", "Equus caballus (horse)"
    ErythrocebusPatas, "Red guenon", "Erythrocebus patas (red guenon)"
    EscherichiaColi, "E. coli", "Escherichia coli (E. coli)"
    EsoxLucius, "Northern pike", "Esox lucius (northern pike)"
    EublepharisMacularius, "Leopard gecko", "Eublepharis macularius (Leopard gecko)"
    EulemurFulvus, "Brown lemur", "Eulemur fulvus (brown lemur)"
    FelineLeukemiaVirus, "Feline leukemia virus", "Feline leukemia virus"
    FelisCatus, "Domestic cat", "Felis catus (domestic cat)"
    FelisSp, "Cats", "Felis sp."
    GadusMorhua, "Atlantic cod", "Gadus morhua (Atlantic cod)"
    GalagoSenegalensis, "Senegal galago", "Galago senegalensis (Senegal galago)"
    GallusGallus, "Domestic chicken", "Gallus gallus (chicken)"
    GasterosteusAculeatus, "Three-spined stickleback", "Gasterosteus aculeatus (three-spined stickleback)"
    GinglymostomaCirratum, "Nurse shark", "Ginglymostoma cirratum (nurse shark)"
    GobionotothenGibberifrons, "Humped rockcod", "Gobionotothen gibberifrons (humped rockcod)"
    GorillaGorilla, "Western gorilla", "Gorilla gorilla (western gorilla)"
    GorillaGorillaGorilla, "Western lowland gorilla", "Gorilla gorilla gorilla (western lowland gorilla)"
    GrampusGriseus, "Risso's dolphin", "Grampus griseus (Risso's dolphin)"
    GymnodracoAcuticeps, "Ploughfish", "Gymnodraco acuticeps"
    GymnogypsCalifornianus, "California condor", "Gymnogyps californianus (California condor)"
    HaemorhousMexicanus, "House finch", "Haemorhous mexicanus (house finch)"
    HemibagrusMacropterus, "Largefin longbarbel catfish", "Hemibagrus macropterus"
    HepacivirusC, "Hepacivirus C", "Hepacivirus C"
    HeterocephalusGlaber, "Naked mole-rat", "Heterocephalus glaber (naked mole-rat)"
    HeterodontusFrancisci, "Horn shark", "Heterodontus francisci (horn shark)"
    HippoglossusHippoglossus, "Atlantic halibut", "Hippoglossus hippoglossus (Atlantic halibut)"
    HistiodracoVelifer, "Histiodraco velifer", "Histiodraco velifer"
    HomoSapiens, "Human", "Homo sapiens (human)"
    HoolockHoolock, "Hoolock gibbon", "Hoolock hoolock (hoolock gibbon)"
    HusoHuso, "Beluga", "Huso huso (beluga)"
    HydrolagusColliei, "Spotted ratfish", "Hydrolagus colliei (spotted ratfish)"
    HylobatesLar, "Common gibbon", "Hylobates lar (common gibbon)"
    IctalurusPunctatus, "Channel catfish", "Ictalurus punctatus (channel catfish)"
    IsoodonMacrourus, "Northern brown bandicoot", "Isoodon macrourus (northern brown bandicoot)"
    KogiaSima, "Dwarf sperm whale", "Kogia sima (dwarf sperm whale)"
    LabeobarbusIntermedius, "Labeobarbus intermedius", "Labeobarbus intermedius"
    LabeoRohita, "Rohu", "Labeo rohita (rohu)"
    LamaGlama, "Llama", "Lama glama (llama)"
    LarimichthysCrocea, "Large yellow croaker", "Larimichthys crocea (large yellow croaker)"
    LatimeriaChalumnae, "Coelacanth", "Latimeria chalumnae (coelacanth)"
    LatimeriaMenadoensis, "Menado coelacanth", "Latimeria menadoensis (Menado coelacanth)"
    LatrisLineata, "Striped trumpeter", "Latris lineata (striped trumpeter)"
    LemurCatta, "Ring-tailed lemur", "Lemur catta (Ring-tailed lemur)"
    LeontopithecusRosalia, "Golden lion tamarin", "Leontopithecus rosalia (golden lion tamarin)"
    LepilemurRuficaudatus, "Red-tailed sportive lemur", "Lepilemur ruficaudatus (red-tailed sportive lemur)"
    LepisosteusOsseus, "Longnose gar", "Lepisosteus osseus (longnose gar)"
    LepusAmericanus, "Snowshoe hare", "Lepus americanus (snowshoe hare)"
    LepusCalifornicus, "Black-tailed jackrabbit", "Lepus californicus (black-tailed jackrabbit)"
    LepusCallotis, "White-sided jackrabbit", "Lepus callotis (white-sided jackrabbit)"
    LepusCapensis, "Brown hare", "Lepus capensis (brown hare)"
    LepusCastroviejoi, "Broom Hare", "Lepus castroviejoi (Broom Hare)"
    LepusEuropaeus, "European hare", "Lepus europaeus (European hare)"
    LepusGranatensis, "Granada hare", "Lepus granatensis (Granada hare)"
    LepusSaxatilis, "Scrub hare", "Lepus saxatilis (scrub hare)"
    LepusTimidus, "Mountain hare", "Lepus timidus (Mountain hare)"
    LeucorajaErinacea, "Little skate", "Leucoraja erinacea (little skate)"
    LipotesVexillifer, "Yangtze River dolphin", "Lipotes vexillifer (Yangtze River dolphin)"
    LutjanusSanguineus, "Humphead snapper", "Lutjanus sanguineus (humphead snapper)"
    MacacaArctoides, "Stump-tailed macaque", "Macaca arctoides (stump-tailed macaque)"
    MacacaAssamensis, "Assam macaque", "Macaca assamensis (Assam macaque)"
    MacacaCyclopis, "Taiwan macaque", "Macaca cyclopis (Taiwan macaque)"
    MacacaFascicularis, "Crab-eating macaque", "Macaca fascicularis (crab-eating macaque)"
    MacacaMulatta, "Rhesus monkey", "Macaca mulatta (Rhesus monkey)"
    MacacaNemestrina, "Pig-tailed macaque", "Macaca nemestrina (pig-tailed macaque)"
    MacacaSilenus, "Liontail macaque", "Macaca silenus (liontail macaque)"
    MacacaThibetana, "Pere David's macaque", "Macaca thibetana (Pere David's macaque)"
    MarecaStrepera, "Gadwall", "Mareca strepera (gadwall)"
    MarmotaHimalayana, "Himalayan marmot", "Marmota himalayana (Himalayan marmot)"
    MarmotaMonax, "Woodchuck", "Marmota monax (woodchuck)"
    MauremysMutica, "Yellowpond turtle", "Mauremys mutica (yellowpond turtle)"
    MelanogrammusAeglefinus, "Haddock", "Melanogrammus aeglefinus (haddock)"
    MeleagrisGallopavo, "Turkey", "Meleagris gallopavo (turkey)"
    MerionesUnguiculatus, "Mongolian gerbil", "Meriones unguiculatus (Mongolian gerbil)"
    MesocricetusAuratus, "Golden hamster", "Mesocricetus auratus (golden hamster)"
    MicrocebusMurinus, "Gray mouse lemur", "Microcebus murinus (gray mouse lemur)"
    MonodelphisDomestica, "Gray short-tailed opossum", "Monodelphis domestica (gray short-tailed opossum)"
    MurinaeSp, "Old world rats and mice", "Murinae gen. sp."
    Mus, "mouse", "Mus (mouse)"
    MusCookii, "Cook's mouse", "Mus cookii (Cook's mouse)"
    MusMinutoides, "Southern African pygmy mouse", "Mus minutoides (Southern African pygmy mouse)"
    MusMusculus, "House mouse", "Mus musculus (house mouse)"
    MusMusculusCastaneus, "Southeastern Asian house mouse", "Mus musculus castaneus (southeastern Asian house mouse)"
    MusMusculusDomesticus, "Western European house mouse", "Mus musculus domesticus (western European house mouse)"
    MusMusculusMolossinus, "Japanese wild mouse", "Mus musculus molossinus (Japanese wild mouse)"
    MusMusculusMusculus, "Eastern European house mouse", "Mus musculus musculus (eastern European house mouse)"
    MusPahari, "Shrew mouse", "Mus pahari (shrew mouse)"
    MusSaxicola, "Spiny mouse", "Mus saxicola (spiny mouse)"
    MusSp, "Mice", "Mus sp. (mice)"
    MusSpretus, "Western wild mouse", "Mus spretus (western wild mouse)"
    MustelaPutoriusFuro, "Domestic ferret", "Mustela putorius furo (domestic ferret)"
    MustelaSp, "Ferret", "Mustela sp."
    MyotisLucifugus, "Little brown bat", "Myotis lucifugus (little brown bat)"
    NeophocaenaPhocaenoides, "Indo-Pacific finless porpoise", "Neophocaena phocaenoides (Indo-Pacific finless porpoise)"
    NeovisonVison, "American mink", "Neovison vison (American mink)"
    NomascusConcolor, "Black crested gibbon", "Nomascus concolor (Black crested gibbon)"
    NotamacropusEugenii, "Tammar wallaby", "Notamacropus eugenii (tammar wallaby)"
    NototheniaCoriiceps, "Black rockcod", "Notothenia coriiceps (black rockcod)"
    NycticebusCoucang, "Slow loris", "Nycticebus coucang (slow loris)"
    OncorhynchusGorbuscha, "Pink salmon", "Oncorhynchus gorbuscha (pink salmon)"
    OncorhynchusMykiss, "Rainbow trout", "Oncorhynchus mykiss (rainbow trout)"
    OncorhynchusTshawytscha, "Chinook salmon", "Oncorhynchus tshawytscha (Chinook salmon)"
    OrectolobusMaculatus, "Spotted wobbegong", "Orectolobus maculatus (spotted wobbegong)"
    OreochromisNiloticus, "Nile tilapia", "Oreochromis niloticus (Nile tilapia)"
    OrnithorhynchusAnatinus, "Platypus", "Ornithorhynchus anatinus (platypus)"
    OryctolagusCuniculus, "Rabbit", "Oryctolagus cuniculus (rabbit)"
    OryctolagusCuniculusAlgirus, "European rabbit", "Oryctolagus cuniculus algirus"
    OryctolagusCuniculusCuniculus, "Rabbit", "Oryctolagus cuniculus cuniculus"
    OryziasLatipes, "Japanese medaka", "Oryzias latipes (Japanese medaka)"
    OryziasMelastigma, "Indian medaka", "Oryzias melastigma (Indian medaka)"
    OtolemurCrassicaudatus, "Thick-tailed bush baby", "Otolemur crassicaudatus (thick-tailed bush baby)"
    OvisAries, "Domestic sheep", "Ovis aries (sheep)"
    OvisSp, "Sheep", "Ovis sp."
    PacifastacusLeniusculus, "Signal crayfish", "Pacifastacus leniusculus (signal crayfish)"
    PagetopsisMacropterus, "Pagetopsis macropterus", "Pagetopsis macropterus"
    PagrusMajor, "Red seabream", "Pagrus major (red seabream)"
    PangasianodonHypophthalmus, "Striped catfish", "Pangasianodon hypophthalmus (striped catfish)"
    PanPaniscus, "Pygmy chimpanzee", "Pan paniscus (pygmy chimpanzee)"
    PantheraPardus, "Leopard", "Panthera pardus (leopard)"
    PanTroglodytes, "Chimpanzee", "Pan troglodytes (chimpanzee)"
    PanTroglodytesVerus, "Western chimpanzee", "Pan troglodytes verus"
    PapioAnubis, "Olive baboon", "Papio anubis (olive baboon)"
    PapioAnubisAnubis, "Olive baboon anubis", "Papio anubis anubis"
    PapioHamadryas, "Hamadryas baboon", "Papio hamadryas (hamadryas baboon)"
    PapioPapio, "Guinea baboon", "Papio papio (Guinea baboon)"
    ParalichthysOlivaceus, "Japanese flounder", "Paralichthys olivaceus (Japanese flounder)"
    PelodiscusSinensis, "Chinese soft-shelled turtle", "Pelodiscus sinensis (Chinese soft-shelled turtle)"
    PelteobagrusFulvidraco, "Yellow catfish", "Pelteobagrus fulvidraco (yellow catfish)"
    PerdixPerdix, "Grey partridge", "Perdix perdix (grey partridge)"
    PeromyscusManiculatus, "North American deer mouse", "Peromyscus maniculatus (North American deer mouse)"
    PetromyzonMarinus, "Sea lamprey", "Petromyzon marinus (sea lamprey)"
    PhascogaleCalura, "Red-tailed phascogale", "Phascogale calura (red-tailed phascogale)"
    PhasianusColchicus, "Ring-necked pheasant", "Phasianus colchicus (Ring-necked pheasant)"
    PhyseterCatodon, "Sperm whale", "Physeter catodon (sperm whale)"
    PitheciaPithecia, "White-faced saki", "Pithecia pithecia (white-faced saki)"
    PlataleaAjaja, "Roseate spoonbil", "Platalea ajaja"
    Platyrrhini, "New World monkeys", "Platyrrhini (New World monkeys)"
    PlecoglossusAltivelisAltivelis, "Ayu sweetfish", "Plecoglossus altivelis altivelis"
    PleurodelesWaltl, "Iberian ribbed newt", "Pleurodeles waltl (Iberian ribbed newt)"
    PogonophryneScotti, "Pogonophryne scotti", "Pogonophryne scotti"
    PolyprionOxygeneios, "Hāpuku", "Polyprion oxygeneios"
    PongoAbelii, "Sumatran orangutan", "Pongo abelii (Sumatran orangutan)"
    PongoPygmaeus, "Bornean orangutan", "Pongo pygmaeus (Bornean orangutan)"
    PresbytisComata, "Grizzled leaf monkey", "Presbytis comata (grizzled leaf monkey)"
    PresbytisFemoralis, "Banded leaf monkey", "Presbytis femoralis (banded leaf monkey)"
    PresbytisMelalophos, "Mitred leaf monkey", "Presbytis melalophos (mitred leaf monkey)"
    PropithecusVerreauxi, "White sifaka", "Propithecus verreauxi (white sifaka)"
    ProtopterusAethiopicus, "Marbled lungfish", "Protopterus aethiopicus (marbled lungfish)"
    PseudobatosProductus, "Shovelnose guitarfish", "Pseudobatos productus (shovelnose guitarfish)"
    PteropusAlecto, "Black flying fox", "Pteropus alecto (black flying fox)"
    PythonBivittatus, "Burmese python", "Python bivittatus (Burmese python)"
    RachycentronCanadum, "Cobia", "Rachycentron canadum (cobia)"
    RajaEglanteria, "Clearnose skate", "Raja eglanteria (clearnose skate)"
    RattusFuscipes, "Bush rat", "Rattus fuscipes (bush rat)"
    RattusLeucopus, "Mottle-tailed rat", "Rattus leucopus (mottle-tailed rat)"
    RattusNorvegicus, "Norway rat", "Rattus norvegicus (Norway rat)"
    RattusRattus, "Black rat", "Rattus rattus (black rat)"
    RattusSordidus, "Australian dusky field rat", "Rattus sordidus (Australian dusky field rat)"
    RattusSp, "Rats", "Rattus sp. (rats)"
    RattusTunneyi, "Tunney's rat", "Rattus tunneyi (Tunney's rat)"
    RattusVillosissimus, "Long-haired rat", "Rattus villosissimus (long-haired rat)"
    RhinocerosUnicornis, "Greater Indian rhinoceros", "Rhinoceros unicornis (greater Indian rhinoceros)"
    RousettusLeschenaultii, "Leschenault's rousette", "Rousettus leschenaultii (Leschenault's rousette)"
    SaguinusLabiatus, "Red-chested mustached tamarin", "Saguinus labiatus (red-chested mustached tamarin)"
    SaguinusMidas, "Midas tamarin", "Saguinus midas (Midas tamarin)"
    SaguinusOedipus, "Cotton-top tamarin", "Saguinus oedipus (cotton-top tamarin)"
    SaimiriBoliviensisBoliviensis, "Bolivian squirrel monkey", "Saimiri boliviensis boliviensis (Bolivian squirrel monkey)"
    SaimiriSciureus, "Common squirrel monkey", "Saimiri sciureus (common squirrel monkey)"
    SalmoMarmoratus, "Salmo marmoratus", "Salmo marmoratus"
    SalmoSalar, "Atlantic salmon", "Salmo salar (Atlantic salmon)"
    SalmoTrutta, "River trout", "Salmo trutta (river trout)"
    SalvelinusAlpinus, "Arctic char", "Salvelinus alpinus (Arctic char)"
    SanderVitreus, "Walleye", "Sander vitreus (walleye)"
    SapajusApella, "Tufted capuchin", "Sapajus apella (Tufted capuchin)"
    SciaenopsOcellatus, "Red drum", "Sciaenops ocellatus (red drum)"
    ScophthalmusMaximus, "Turbot", "Scophthalmus maximus (turbot)"
    ScyliorhinusCanicula, "Smaller spotted catshark", "Scyliorhinus canicula (smaller spotted catshark)"
    SeriolaQuinqueradiata, "Japanese amberjack", "Seriola quinqueradiata (Japanese amberjack)"
    SilurusAsotus, "Amur catfish", "Silurus asotus (Amur catfish)"
    SilurusMeridionalis, "Silurus meridionalis", "Silurus meridionalis"
    SinipercaChuatsi, "Mandarin fish", "Siniperca chuatsi (mandarin fish)"
    SousaChinensis, "Indo-pacific humpbacked dolphin", "Sousa chinensis (Indo-pacific humpbacked dolphin)"
    SparusAurata, "Gilthead seabream", "Sparus aurata (gilthead seabream)"
    SphoeroidesNephelus, "Southern puffer", "Sphoeroides nephelus (southern puffer)"
    SqualusAcanthias, "Spiny dogfish", "Squalus acanthias (spiny dogfish)"
    StegastesLeucostictus, "Beaugregory", "Stegastes leucostictus (beaugregory)"
    StegastesPartitus, "Bicolor damselfish", "Stegastes partitus (bicolor damselfish)"
    StenellaAttenuata, "Bridled dolphin", "Stenella attenuata (bridled dolphin)"
    StenellaCoeruleoalba, "Striped dolphin", "Stenella coeruleoalba (striped dolphin)"
    StreptomycesViridochromogenes, "Streptomyces viridochromogenes", "Streptomyces viridochromogenes"
    StruthioCamelus, "African ostrich", "Struthio camelus (African ostrich)"
    SuncusMurinus, "House shrew", "Suncus murinus (house shrew)"
    SusScrofa, "Domestic pig", "Sus scrofa (pig)"
    SylvilagusCunicularis, "Mexican cottontail", "Sylvilagus cunicularis (Mexican cottontail)"
    SylvilagusFloridanus, "Eastern cottontail", "Sylvilagus floridanus (eastern cottontail)"
    SymphalangusSyndactylus, "Siamang", "Symphalangus syndactylus (siamang)"
    TachyglossusAculeatus, "Australian echidna", "Tachyglossus aculeatus (Australian echidna)"
    TachysurusFulvidraco, "Yellow catfish", "Tachysurus fulvidraco (yellow catfish)"
    TachysurusVachellii, "Tachysurus vachellii", "Tachysurus vachellii"
    TaeniopygiaGuttata, "Zebra finch", "Taeniopygia guttata (zebra finch)"
    TakifuguRubripes, "Torafugu", "Takifugu rubripes (torafugu)"
    TarsiusDentatus, "Dian's tarsier", "Tarsius dentatus (Dian's tarsier)"
    TarsiusLariang, "Lariang tarsier", "Tarsius lariang (Lariang tarsier)"
    TarsiusSyrichta, "Philippine tarsier", "Tarsius syrichta (Philippine tarsier)"
    TetraodonNigroviridis, "Spotted green pufferfish", "Tetraodon nigroviridis (spotted green pufferfish)"
    TrachemysScripta, "Red-eared slider turtle", "Trachemys scripta (red-eared slider turtle)"
    TrachemysScriptaElegans, "Red-eared slider turtle elegans", "Trachemys scripta elegans"
    TrachypithecusCristatus, "Silvery lutung", "Trachypithecus cristatus (Silvery lutung)"
    TrachypithecusObscurus, "Dusky leaf-monkey", "Trachypithecus obscurus (Dusky leaf-monkey)"
    TrematomusBernacchii, "Emerald rockcod", "Trematomus bernacchii (emerald rockcod)"
    TrematomusHansoni, "Striped rockcod", "Trematomus hansoni (striped rockcod)"
    TrematomusLoennbergii, "Deepwater notothen", "Trematomus loennbergii (deepwater notothen)"
    TrematomusNewnesi, "Dusky notothen", "Trematomus newnesi (dusky notothen)"
    TrematomusPennellii, "Sharp-spined notothen", "Trematomus pennellii (sharp-spined notothen)"
    TriakisScyllium, "Banded houndshark", "Triakis scyllium (banded houndshark)"
    TrichechusManatusLatirostris, "Florida manatee", "Trichechus manatus latirostris (Florida manatee)"
    TrichosurusVulpecula, "Common brushtail", "Trichosurus vulpecula (common brushtail)"
    TursiopsAduncus, "Indo-pacific bottlenose dolphin", "Tursiops aduncus (Indo-pacific bottlenose dolphin)"
    TursiopsTruncatus, "Common bottlenose dolphin", "Tursiops truncatus (common bottlenose dolphin)"
    VicugnaPacos, "Alpaca", "Vicugna pacos (alpaca)"
    Xenopus, "Xenopus", "Xenopus"
    XenopusLaevis, "African clawed frog", "Xenopus laevis (African clawed frog)"
    XenopusLaevisOrGilli, "African or Cape clawed frog", "Xenopus laevis/gilli"
    XenopusSp, "Clawed frog", "Xenopus sp. (clawed frog)"
    XenopusTropicalis, "Tropical clawed frog", "Xenopus tropicalis (tropical clawed frog)"
);
