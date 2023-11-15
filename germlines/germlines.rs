use std::{sync::OnceLock, collections::HashMap};
use crate::shared::{Germlines, Species};
/// Get the germlines for any of the available species. See the tables below for which species have which data available.
///
/// ## Bos taurus / Domestic bovine
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|12/17|3/10|7/9|
/// |IGKV|6/6|1/1|0/0|
/// |IGLV|25/41|5/9|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Camelus dromedarius / Arabian camel
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|0/0|0/0|0/0|
/// |IGKV|8/8|4/4|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Canis lupus familiaris / Domestic dog
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|35/35|5/5|4/9|
/// |IGKV|13/24|4/8|0/0|
/// |IGLV|71/71|5/5|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Capra hircus / Domestic goat
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|0/0|0/0|0/0|
/// |IGKV|5/5|1/1|0/0|
/// |IGLV|24/24|1/1|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Danio rerio / Zebrafish
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|1/1|0/0|0/0|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Equus caballus / Domestic horse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|19/32|6/12|10/19|
/// |IGKV|17/33|4/12|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Felis catus / Domestic cat
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|0/0|0/0|0/0|
/// |IGKV|12/12|5/5|0/0|
/// |IGLV|32/32|10/10|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Gallus gallus / Domestic chicken
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|2/3|1/3|0/0|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|2/3|1/3|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Gorilla gorilla gorilla / Western lowland gorilla
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|53/112|5/20|9/15|
/// |IGKV|26/57|5/15|0/0|
/// |IGLV|27/54|6/13|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Homo sapiens / Human
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|55/275|6/29|10/63|
/// |IGKV|40/74|5/10|0/0|
/// |IGLV|26/65|4/5|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Ictalurus punctatus / Channel catfish
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|0/0|8/8|0/0|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Lemur catta / Ring-tailed lemur
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|56/56|4/4|3/3|
/// |IGKV|11/28|4/8|0/0|
/// |IGLV|119/119|10/10|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Macaca fascicularis / Crab-eating macaque
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|1/1|6/6|1/1|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Macaca mulatta / Rhesus monkey
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|122/361|7/28|8/17|
/// |IGKV|58/83|5/10|0/0|
/// |IGLV|55/107|7/12|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Mus musculus / House mouse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|114/252|4/14|7/12|
/// |IGKV|98/211|4/8|0/0|
/// |IGLV|3/8|3/14|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Mus musculus domesticus / Western European house mouse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|1/1|0/0|0/0|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|1/1|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Mus spretus / Western wild mouse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|0/0|0/0|0/0|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|2/2|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Mustela putorius furo / Domestic ferret
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|0/0|0/0|0/0|
/// |IGKV|37/37|4/4|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Oncorhynchus mykiss / Rainbow trout
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|42/46|18/27|1/2|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Ornithorhynchus anatinus / Platypus
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|34/34|10/10|6/6|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Oryctolagus cuniculus / Rabbit
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|9/11|6/15|7/8|
/// |IGKV|5/7|5/10|0/0|
/// |IGLV|20/20|1/1|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Ovis aries / Domestic sheep
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|0/0|0/0|0/0|
/// |IGKV|5/5|1/1|0/0|
/// |IGLV|18/22|1/1|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Pongo abelii / Sumatran orangutan
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|47/47|3/3|6/6|
/// |IGKV|30/30|4/4|0/0|
/// |IGLV|25/25|4/4|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Pongo pygmaeus / Bornean orangutan
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|53/53|3/3|9/9|
/// |IGKV|34/34|4/4|0/0|
/// |IGLV|28/28|5/5|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Rattus norvegicus / Norway rat
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|97/109|4/5|7/8|
/// |IGKV|20/20|6/6|0/0|
/// |IGLV|7/7|2/3|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Salmo salar / Atlantic salmon
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|66/80|19/26|2/4|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Sus scrofa / Domestic pig
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|6/6|5/10|10/20|
/// |IGKV|11/19|5/10|0/0|
/// |IGLV|10/17|3/4|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
/// ## Vicugna pacos / Alpaca
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |IGHV|4/4|6/6|5/5|
/// |IGKV|0/0|0/0|0/0|
/// |IGLV|0/0|0/0|0/0|
/// |IGIV|0/0|0/0|0/0|
/// 
/// _Number of genes / number of alleles_
/// 
pub fn germlines(species: Species) -> Option<&'static Germlines> {match species {
Species::BosTaurus => Some(lock_BosTaurus()),
Species::CamelusDromedarius => Some(lock_CamelusDromedarius()),
Species::CanisLupusFamiliaris => Some(lock_CanisLupusFamiliaris()),
Species::CapraHircus => Some(lock_CapraHircus()),
Species::DanioRerio => Some(lock_DanioRerio()),
Species::EquusCaballus => Some(lock_EquusCaballus()),
Species::FelisCatus => Some(lock_FelisCatus()),
Species::GallusGallus => Some(lock_GallusGallus()),
Species::GorillaGorillaGorilla => Some(lock_GorillaGorillaGorilla()),
Species::HomoSapiens => Some(lock_HomoSapiens()),
Species::IctalurusPunctatus => Some(lock_IctalurusPunctatus()),
Species::LemurCatta => Some(lock_LemurCatta()),
Species::MacacaFascicularis => Some(lock_MacacaFascicularis()),
Species::MacacaMulatta => Some(lock_MacacaMulatta()),
Species::MusMusculus => Some(lock_MusMusculus()),
Species::MusMusculusDomesticus => Some(lock_MusMusculusDomesticus()),
Species::MusSpretus => Some(lock_MusSpretus()),
Species::MustelaPutoriusFuro => Some(lock_MustelaPutoriusFuro()),
Species::OncorhynchusMykiss => Some(lock_OncorhynchusMykiss()),
Species::OrnithorhynchusAnatinus => Some(lock_OrnithorhynchusAnatinus()),
Species::OryctolagusCuniculus => Some(lock_OryctolagusCuniculus()),
Species::OvisAries => Some(lock_OvisAries()),
Species::PongoAbelii => Some(lock_PongoAbelii()),
Species::PongoPygmaeus => Some(lock_PongoPygmaeus()),
Species::RattusNorvegicus => Some(lock_RattusNorvegicus()),
Species::SalmoSalar => Some(lock_SalmoSalar()),
Species::SusScrofa => Some(lock_SusScrofa()),
Species::VicugnaPacos => Some(lock_VicugnaPacos()),
_=>None}}
pub fn all_germlines() -> impl std::iter::Iterator<Item = &'static Germlines> {
std::iter::once(lock_BosTaurus())
.chain(std::iter::once(lock_CamelusDromedarius()))
.chain(std::iter::once(lock_CanisLupusFamiliaris()))
.chain(std::iter::once(lock_CapraHircus()))
.chain(std::iter::once(lock_DanioRerio()))
.chain(std::iter::once(lock_EquusCaballus()))
.chain(std::iter::once(lock_FelisCatus()))
.chain(std::iter::once(lock_GallusGallus()))
.chain(std::iter::once(lock_GorillaGorillaGorilla()))
.chain(std::iter::once(lock_HomoSapiens()))
.chain(std::iter::once(lock_IctalurusPunctatus()))
.chain(std::iter::once(lock_LemurCatta()))
.chain(std::iter::once(lock_MacacaFascicularis()))
.chain(std::iter::once(lock_MacacaMulatta()))
.chain(std::iter::once(lock_MusMusculus()))
.chain(std::iter::once(lock_MusMusculusDomesticus()))
.chain(std::iter::once(lock_MusSpretus()))
.chain(std::iter::once(lock_MustelaPutoriusFuro()))
.chain(std::iter::once(lock_OncorhynchusMykiss()))
.chain(std::iter::once(lock_OrnithorhynchusAnatinus()))
.chain(std::iter::once(lock_OryctolagusCuniculus()))
.chain(std::iter::once(lock_OvisAries()))
.chain(std::iter::once(lock_PongoAbelii()))
.chain(std::iter::once(lock_PongoPygmaeus()))
.chain(std::iter::once(lock_RattusNorvegicus()))
.chain(std::iter::once(lock_SalmoSalar()))
.chain(std::iter::once(lock_SusScrofa()))
.chain(std::iter::once(lock_VicugnaPacos()))
}
static LOCK_BosTaurus: OnceLock<Germlines> = OnceLock::new();
fn lock_BosTaurus()->&'static Germlines{LOCK_BosTaurus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic bovine.bin")).unwrap()})}
static LOCK_CamelusDromedarius: OnceLock<Germlines> = OnceLock::new();
fn lock_CamelusDromedarius()->&'static Germlines{LOCK_CamelusDromedarius.get_or_init(|| {bincode::deserialize(include_bytes!("Arabian camel.bin")).unwrap()})}
static LOCK_CanisLupusFamiliaris: OnceLock<Germlines> = OnceLock::new();
fn lock_CanisLupusFamiliaris()->&'static Germlines{LOCK_CanisLupusFamiliaris.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic dog.bin")).unwrap()})}
static LOCK_CapraHircus: OnceLock<Germlines> = OnceLock::new();
fn lock_CapraHircus()->&'static Germlines{LOCK_CapraHircus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic goat.bin")).unwrap()})}
static LOCK_DanioRerio: OnceLock<Germlines> = OnceLock::new();
fn lock_DanioRerio()->&'static Germlines{LOCK_DanioRerio.get_or_init(|| {bincode::deserialize(include_bytes!("Zebrafish.bin")).unwrap()})}
static LOCK_EquusCaballus: OnceLock<Germlines> = OnceLock::new();
fn lock_EquusCaballus()->&'static Germlines{LOCK_EquusCaballus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic horse.bin")).unwrap()})}
static LOCK_FelisCatus: OnceLock<Germlines> = OnceLock::new();
fn lock_FelisCatus()->&'static Germlines{LOCK_FelisCatus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic cat.bin")).unwrap()})}
static LOCK_GallusGallus: OnceLock<Germlines> = OnceLock::new();
fn lock_GallusGallus()->&'static Germlines{LOCK_GallusGallus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic chicken.bin")).unwrap()})}
static LOCK_GorillaGorillaGorilla: OnceLock<Germlines> = OnceLock::new();
fn lock_GorillaGorillaGorilla()->&'static Germlines{LOCK_GorillaGorillaGorilla.get_or_init(|| {bincode::deserialize(include_bytes!("Western lowland gorilla.bin")).unwrap()})}
static LOCK_HomoSapiens: OnceLock<Germlines> = OnceLock::new();
fn lock_HomoSapiens()->&'static Germlines{LOCK_HomoSapiens.get_or_init(|| {bincode::deserialize(include_bytes!("Human.bin")).unwrap()})}
static LOCK_IctalurusPunctatus: OnceLock<Germlines> = OnceLock::new();
fn lock_IctalurusPunctatus()->&'static Germlines{LOCK_IctalurusPunctatus.get_or_init(|| {bincode::deserialize(include_bytes!("Channel catfish.bin")).unwrap()})}
static LOCK_LemurCatta: OnceLock<Germlines> = OnceLock::new();
fn lock_LemurCatta()->&'static Germlines{LOCK_LemurCatta.get_or_init(|| {bincode::deserialize(include_bytes!("Ring-tailed lemur.bin")).unwrap()})}
static LOCK_MacacaFascicularis: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaFascicularis()->&'static Germlines{LOCK_MacacaFascicularis.get_or_init(|| {bincode::deserialize(include_bytes!("Crab-eating macaque.bin")).unwrap()})}
static LOCK_MacacaMulatta: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaMulatta()->&'static Germlines{LOCK_MacacaMulatta.get_or_init(|| {bincode::deserialize(include_bytes!("Rhesus monkey.bin")).unwrap()})}
static LOCK_MusMusculus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculus()->&'static Germlines{LOCK_MusMusculus.get_or_init(|| {bincode::deserialize(include_bytes!("House mouse.bin")).unwrap()})}
static LOCK_MusMusculusDomesticus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculusDomesticus()->&'static Germlines{LOCK_MusMusculusDomesticus.get_or_init(|| {bincode::deserialize(include_bytes!("Western European house mouse.bin")).unwrap()})}
static LOCK_MusSpretus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusSpretus()->&'static Germlines{LOCK_MusSpretus.get_or_init(|| {bincode::deserialize(include_bytes!("Western wild mouse.bin")).unwrap()})}
static LOCK_MustelaPutoriusFuro: OnceLock<Germlines> = OnceLock::new();
fn lock_MustelaPutoriusFuro()->&'static Germlines{LOCK_MustelaPutoriusFuro.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic ferret.bin")).unwrap()})}
static LOCK_OncorhynchusMykiss: OnceLock<Germlines> = OnceLock::new();
fn lock_OncorhynchusMykiss()->&'static Germlines{LOCK_OncorhynchusMykiss.get_or_init(|| {bincode::deserialize(include_bytes!("Rainbow trout.bin")).unwrap()})}
static LOCK_OrnithorhynchusAnatinus: OnceLock<Germlines> = OnceLock::new();
fn lock_OrnithorhynchusAnatinus()->&'static Germlines{LOCK_OrnithorhynchusAnatinus.get_or_init(|| {bincode::deserialize(include_bytes!("Platypus.bin")).unwrap()})}
static LOCK_OryctolagusCuniculus: OnceLock<Germlines> = OnceLock::new();
fn lock_OryctolagusCuniculus()->&'static Germlines{LOCK_OryctolagusCuniculus.get_or_init(|| {bincode::deserialize(include_bytes!("Rabbit.bin")).unwrap()})}
static LOCK_OvisAries: OnceLock<Germlines> = OnceLock::new();
fn lock_OvisAries()->&'static Germlines{LOCK_OvisAries.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic sheep.bin")).unwrap()})}
static LOCK_PongoAbelii: OnceLock<Germlines> = OnceLock::new();
fn lock_PongoAbelii()->&'static Germlines{LOCK_PongoAbelii.get_or_init(|| {bincode::deserialize(include_bytes!("Sumatran orangutan.bin")).unwrap()})}
static LOCK_PongoPygmaeus: OnceLock<Germlines> = OnceLock::new();
fn lock_PongoPygmaeus()->&'static Germlines{LOCK_PongoPygmaeus.get_or_init(|| {bincode::deserialize(include_bytes!("Bornean orangutan.bin")).unwrap()})}
static LOCK_RattusNorvegicus: OnceLock<Germlines> = OnceLock::new();
fn lock_RattusNorvegicus()->&'static Germlines{LOCK_RattusNorvegicus.get_or_init(|| {bincode::deserialize(include_bytes!("Norway rat.bin")).unwrap()})}
static LOCK_SalmoSalar: OnceLock<Germlines> = OnceLock::new();
fn lock_SalmoSalar()->&'static Germlines{LOCK_SalmoSalar.get_or_init(|| {bincode::deserialize(include_bytes!("Atlantic salmon.bin")).unwrap()})}
static LOCK_SusScrofa: OnceLock<Germlines> = OnceLock::new();
fn lock_SusScrofa()->&'static Germlines{LOCK_SusScrofa.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic pig.bin")).unwrap()})}
static LOCK_VicugnaPacos: OnceLock<Germlines> = OnceLock::new();
fn lock_VicugnaPacos()->&'static Germlines{LOCK_VicugnaPacos.get_or_init(|| {bincode::deserialize(include_bytes!("Alpaca.bin")).unwrap()})}
