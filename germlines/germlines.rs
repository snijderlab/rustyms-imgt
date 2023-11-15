use std::{sync::OnceLock, collections::HashMap};
use crate::shared::{Germlines, Species};
/// Get the germlines for any of the available species. See the tables below for which species have which data available.
///
/// Germlines for Bos taurus Domestic bovine
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|12|3|7|
/// |Light|31|6|-|
/// 
/// Germlines for Camelus dromedarius Arabian camel
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|0|0|0|
/// |Light|8|4|-|
/// 
/// Germlines for Canis lupus familiaris Domestic dog
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|35|5|4|
/// |Light|84|9|-|
/// 
/// Germlines for Capra hircus Domestic goat
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|0|0|0|
/// |Light|29|2|-|
/// 
/// Germlines for Danio rerio Zebrafish
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|1|0|0|
/// |Light|0|0|-|
/// 
/// Germlines for Equus caballus Domestic horse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|19|6|10|
/// |Light|17|4|-|
/// 
/// Germlines for Felis catus Domestic cat
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|0|0|0|
/// |Light|44|15|-|
/// 
/// Germlines for Gallus gallus Domestic chicken
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|2|1|0|
/// |Light|2|1|-|
/// 
/// Germlines for Gorilla gorilla gorilla Western lowland gorilla
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|53|5|9|
/// |Light|53|11|-|
/// 
/// Germlines for Homo sapiens Human
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|55|6|10|
/// |Light|66|9|-|
/// 
/// Germlines for Ictalurus punctatus Channel catfish
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|0|8|0|
/// |Light|0|0|-|
/// 
/// Germlines for Lemur catta Ring-tailed lemur
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|56|4|3|
/// |Light|130|14|-|
/// 
/// Germlines for Macaca fascicularis Crab-eating macaque
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|1|6|1|
/// |Light|0|0|-|
/// 
/// Germlines for Macaca mulatta Rhesus monkey
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|122|7|8|
/// |Light|113|12|-|
/// 
/// Germlines for Mus musculus House mouse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|114|4|7|
/// |Light|101|7|-|
/// 
/// Germlines for Mus musculus domesticus Western European house mouse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|1|0|0|
/// |Light|0|1|-|
/// 
/// Germlines for Mus spretus Western wild mouse
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|0|0|0|
/// |Light|2|0|-|
/// 
/// Germlines for Mustela putorius furo Domestic ferret
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|0|0|0|
/// |Light|37|4|-|
/// 
/// Germlines for Oncorhynchus mykiss Rainbow trout
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|42|18|1|
/// |Light|0|0|-|
/// 
/// Germlines for Ornithorhynchus anatinus Platypus
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|34|10|6|
/// |Light|0|0|-|
/// 
/// Germlines for Oryctolagus cuniculus Rabbit
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|9|6|7|
/// |Light|25|6|-|
/// 
/// Germlines for Ovis aries Domestic sheep
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|0|0|0|
/// |Light|23|2|-|
/// 
/// Germlines for Pongo abelii Sumatran orangutan
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|47|3|6|
/// |Light|55|8|-|
/// 
/// Germlines for Pongo pygmaeus Bornean orangutan
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|53|3|9|
/// |Light|62|9|-|
/// 
/// Germlines for Rattus norvegicus Norway rat
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|97|4|7|
/// |Light|27|8|-|
/// 
/// Germlines for Salmo salar Atlantic salmon
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|66|19|2|
/// |Light|0|0|-|
/// 
/// Germlines for Sus scrofa Domestic pig
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|6|5|10|
/// |Light|21|8|-|
/// 
/// Germlines for Vicugna pacos Alpaca
            ///
/// | Kind | V | J | C |
/// |------|---|---|---|
/// |Heavy|4|6|5|
/// |Light|0|0|-|
/// 
pub fn germlines(species: Species) -> Option<&'static Germlines> {match species {
Species::BosTaurus => Some(LOCK0.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic bovine.bin")).unwrap()})),
Species::CamelusDromedarius => Some(LOCK1.get_or_init(|| {bincode::deserialize(include_bytes!("Arabian camel.bin")).unwrap()})),
Species::CanisLupusFamiliaris => Some(LOCK2.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic dog.bin")).unwrap()})),
Species::CapraHircus => Some(LOCK3.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic goat.bin")).unwrap()})),
Species::DanioRerio => Some(LOCK4.get_or_init(|| {bincode::deserialize(include_bytes!("Zebrafish.bin")).unwrap()})),
Species::EquusCaballus => Some(LOCK5.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic horse.bin")).unwrap()})),
Species::FelisCatus => Some(LOCK6.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic cat.bin")).unwrap()})),
Species::GallusGallus => Some(LOCK7.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic chicken.bin")).unwrap()})),
Species::GorillaGorillaGorilla => Some(LOCK8.get_or_init(|| {bincode::deserialize(include_bytes!("Western lowland gorilla.bin")).unwrap()})),
Species::HomoSapiens => Some(LOCK9.get_or_init(|| {bincode::deserialize(include_bytes!("Human.bin")).unwrap()})),
Species::IctalurusPunctatus => Some(LOCK10.get_or_init(|| {bincode::deserialize(include_bytes!("Channel catfish.bin")).unwrap()})),
Species::LemurCatta => Some(LOCK11.get_or_init(|| {bincode::deserialize(include_bytes!("Ring-tailed lemur.bin")).unwrap()})),
Species::MacacaFascicularis => Some(LOCK12.get_or_init(|| {bincode::deserialize(include_bytes!("Crab-eating macaque.bin")).unwrap()})),
Species::MacacaMulatta => Some(LOCK13.get_or_init(|| {bincode::deserialize(include_bytes!("Rhesus monkey.bin")).unwrap()})),
Species::MusMusculus => Some(LOCK14.get_or_init(|| {bincode::deserialize(include_bytes!("House mouse.bin")).unwrap()})),
Species::MusMusculusDomesticus => Some(LOCK15.get_or_init(|| {bincode::deserialize(include_bytes!("Western European house mouse.bin")).unwrap()})),
Species::MusSpretus => Some(LOCK16.get_or_init(|| {bincode::deserialize(include_bytes!("Western wild mouse.bin")).unwrap()})),
Species::MustelaPutoriusFuro => Some(LOCK17.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic ferret.bin")).unwrap()})),
Species::OncorhynchusMykiss => Some(LOCK18.get_or_init(|| {bincode::deserialize(include_bytes!("Rainbow trout.bin")).unwrap()})),
Species::OrnithorhynchusAnatinus => Some(LOCK19.get_or_init(|| {bincode::deserialize(include_bytes!("Platypus.bin")).unwrap()})),
Species::OryctolagusCuniculus => Some(LOCK20.get_or_init(|| {bincode::deserialize(include_bytes!("Rabbit.bin")).unwrap()})),
Species::OvisAries => Some(LOCK21.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic sheep.bin")).unwrap()})),
Species::PongoAbelii => Some(LOCK22.get_or_init(|| {bincode::deserialize(include_bytes!("Sumatran orangutan.bin")).unwrap()})),
Species::PongoPygmaeus => Some(LOCK23.get_or_init(|| {bincode::deserialize(include_bytes!("Bornean orangutan.bin")).unwrap()})),
Species::RattusNorvegicus => Some(LOCK24.get_or_init(|| {bincode::deserialize(include_bytes!("Norway rat.bin")).unwrap()})),
Species::SalmoSalar => Some(LOCK25.get_or_init(|| {bincode::deserialize(include_bytes!("Atlantic salmon.bin")).unwrap()})),
Species::SusScrofa => Some(LOCK26.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic pig.bin")).unwrap()})),
Species::VicugnaPacos => Some(LOCK27.get_or_init(|| {bincode::deserialize(include_bytes!("Alpaca.bin")).unwrap()})),
_=>None}}
static LOCK0: OnceLock<Germlines> = OnceLock::new();
static LOCK1: OnceLock<Germlines> = OnceLock::new();
static LOCK2: OnceLock<Germlines> = OnceLock::new();
static LOCK3: OnceLock<Germlines> = OnceLock::new();
static LOCK4: OnceLock<Germlines> = OnceLock::new();
static LOCK5: OnceLock<Germlines> = OnceLock::new();
static LOCK6: OnceLock<Germlines> = OnceLock::new();
static LOCK7: OnceLock<Germlines> = OnceLock::new();
static LOCK8: OnceLock<Germlines> = OnceLock::new();
static LOCK9: OnceLock<Germlines> = OnceLock::new();
static LOCK10: OnceLock<Germlines> = OnceLock::new();
static LOCK11: OnceLock<Germlines> = OnceLock::new();
static LOCK12: OnceLock<Germlines> = OnceLock::new();
static LOCK13: OnceLock<Germlines> = OnceLock::new();
static LOCK14: OnceLock<Germlines> = OnceLock::new();
static LOCK15: OnceLock<Germlines> = OnceLock::new();
static LOCK16: OnceLock<Germlines> = OnceLock::new();
static LOCK17: OnceLock<Germlines> = OnceLock::new();
static LOCK18: OnceLock<Germlines> = OnceLock::new();
static LOCK19: OnceLock<Germlines> = OnceLock::new();
static LOCK20: OnceLock<Germlines> = OnceLock::new();
static LOCK21: OnceLock<Germlines> = OnceLock::new();
static LOCK22: OnceLock<Germlines> = OnceLock::new();
static LOCK23: OnceLock<Germlines> = OnceLock::new();
static LOCK24: OnceLock<Germlines> = OnceLock::new();
static LOCK25: OnceLock<Germlines> = OnceLock::new();
static LOCK26: OnceLock<Germlines> = OnceLock::new();
static LOCK27: OnceLock<Germlines> = OnceLock::new();
