use std::{
    collections::HashMap,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::ParseIntError,
    ops::RangeInclusive,
    str::FromStr,
};

#[path = "../../shared/mod.rs"]
mod shared;

use crate::shared::*;
use flate2::read::GzDecoder;
use itertools::Itertools;
use rustyms::AminoAcid;

fn main() {
    let file = File::open("../data/imgt.dat.Z").unwrap();
    let mut output = BufWriter::new(File::create("../germlines/germlines.rs").unwrap());
    let mut docs = BufWriter::new(File::create("../germlines/germlines.md").unwrap());
    let mut error = BufWriter::new(File::create("errors.dat").unwrap());
    let data = parse_dat(BufReader::new(GzDecoder::new(file)));
    let mut grouped = HashMap::new();
    for element in data {
        let species = element.as_ref().unwrap().species;
        for gene in element.unwrap().genes {
            match gene.clone().finish(species.to_owned()) {
                Ok(gene) => grouped
                    .entry(species)
                    .or_insert(Germlines::new(species))
                    .insert(gene),
                Err(err) => {
                    writeln!(error, "ERROR FOR GENE:\n{species}\t{gene}{err}\n",).unwrap();
                }
            }
        }
        //writeln!(output, "{}", element.unwrap()).unwrap();
    }

    writeln!(
        output,
        "#![allow(non_snake_case,non_upper_case_globals)]\nuse std::sync::OnceLock;\nuse crate::shared::{{Germlines, Species}};"
    )
    .unwrap();
    writeln!(output, "/// Get the germlines for any of the available species. See the tables below for which species have which data available.").unwrap();
    writeln!(output, "///").unwrap();
    let mut found_species = Vec::new();
    let mut found_germlines: Vec<(Species, Germlines)> = grouped.into_iter().collect();
    found_germlines.sort_unstable_by_key(|g| g.0);
    for (species, germlines) in found_germlines {
        writeln!(
            docs,
            "## {} / {}

| Kind | V | J | C |
|------|---|---|---|
|IGHV{}
|IGKV{}
|IGLV{}
|IGIV{}

_Number of genes / number of alleles_
",
            species.scientific_name(),
            species.common_name(),
            germlines.h.doc_row(),
            germlines.k.doc_row(),
            germlines.l.doc_row(),
            germlines.i.doc_row(),
        )
        .unwrap();
        found_species.push(species);

        let mut file = std::fs::File::create(format!("../germlines/{species}.bin")).unwrap();
        file.write_all(&bincode::serialize::<Germlines>(&germlines).unwrap())
            .unwrap();
    }
    writeln!(
        output,
        "pub fn germlines(species: Species) -> Option<&'static Germlines> {{match species {{"
    )
    .unwrap();

    for species in &found_species {
        writeln!(output, "Species::{0} => Some(lock_{0}()),", species.ident()).unwrap();
    }
    writeln!(output, "_=>None}}}}").unwrap();
    writeln!(
        output,
        "/// Get all germlines in one iterator, see [`germlines()`] for more information about the available germlines\npub fn all_germlines() -> impl std::iter::Iterator<Item = &'static Germlines> {{"
    )
    .unwrap();
    let mut first = true;
    for species in &found_species {
        if first {
            first = false;
            writeln!(output, "std::iter::once(lock_{}())", species.ident()).unwrap();
        } else {
            writeln!(
                output,
                ".chain(std::iter::once(lock_{}()))",
                species.ident()
            )
            .unwrap();
        }
    }
    writeln!(output, "}}").unwrap();

    for species in &found_species {
        writeln!(
            output,
            "static LOCK_{0}: OnceLock<Germlines> = OnceLock::new();\nfn lock_{0}()->&'static Germlines{{LOCK_{0}.get_or_init(|| {{bincode::deserialize(include_bytes!(\"{species}.bin\")).unwrap()}})}}",
            species.ident(),
        )
        .unwrap();
    }
}

fn parse_dat<T: std::io::Read>(
    reader: BufReader<T>,
) -> impl Iterator<Item = Result<DataItem, String>> {
    reader
        .lines()
        .batching(|f| {
            let mut data = PreDataItem::default();
            let mut next = f.next();
            while let Some(Ok(line)) = next {
                if line != "//" {
                    if line.starts_with("ID") {
                        data.id = line;
                    } else if line.starts_with("KW") {
                        data.kw.extend(
                            line[5..]
                                .split(';')
                                .map(|s| s.trim().to_string())
                                .filter(|s| !s.is_empty()),
                        )
                    } else if line.starts_with("FH   Key") {
                        data.ft_key_width = line.find("Location").expect("Incorrect FH line") - 5;
                    } else if line.starts_with("FT") {
                        data.ft.push(line);
                    } else if line.starts_with("OS") && data.os.is_none() {
                        data.os = Species::from_imgt(line[5..].trim()).unwrap_or_else(|()| {
                            println!("Not a species name: `{line}`");
                            None
                        });
                    } else if line.starts_with("  ") {
                        data.sq.extend(
                            line.chars()
                                .filter(|c| *c == 'c' || *c == 'a' || *c == 't' || *c == 'g'),
                        )
                    }
                } else {
                    return Some(data);
                }
                next = f.next();
            }
            None
        })
        .filter(|pre| {
            pre.kw.contains(&"immunoglobulin (IG)".to_string())
                && pre.kw.contains(&"functional".to_string())
                && pre.os.is_some()
        })
        .map(DataItem::new)
}

#[derive(Default, Debug)]
struct PreDataItem {
    id: String,
    kw: Vec<String>,
    ft_key_width: usize,
    ft: Vec<String>,
    os: Option<Species>,
    sq: String,
}

#[derive(Debug)]
struct DataItem {
    id: String,
    genes: Vec<IMGTGene>,
    regions: Vec<Region>,
    species: Species,
    sequence: String,
}

#[derive(Clone, Debug)]
struct IMGTGene {
    key: String,
    location: Location,
    allele: String,
    regions: HashMap<String, Region>,
}

#[derive(Clone, Debug)]
struct Region {
    key: String,
    location: Location,
    reported_seq: String,
    found_seq: Option<(String, Vec<AminoAcid>)>,
    allele: String,
    functional: bool,
    partial: bool,
}

#[derive(Clone, Debug)]
enum Location {
    Normal(RangeInclusive<usize>),
    Complement(RangeInclusive<usize>),
    SingleNormal(usize),
    SingleComplement(usize),
}

impl Location {
    fn contains(&self, other: &Location) -> bool {
        match (self, other) {
            (Self::Complement(s), Self::Complement(o)) | (Self::Normal(s), Self::Normal(o)) => {
                s.start() <= o.start() && s.end() >= o.end()
            }
            (Self::Complement(s), Self::SingleComplement(o)) => s.contains(o),
            (Self::Normal(s), Self::SingleNormal(o)) => s.contains(o),
            _ => false,
        }
    }

    fn get_aa_loc(&self, inner: &Self) -> RangeInclusive<usize> {
        match (self, inner) {
            (Self::Complement(s), Self::Complement(o)) | (Self::Normal(s), Self::Normal(o)) => {
                (o.start() - s.start()) / 3..=(o.end() - s.start()) / 3
            }
            (Self::Normal(s), Self::SingleNormal(o))
            | (Self::Complement(s), Self::SingleComplement(o)) => {
                (o - s.start()) / 3..=(o - s.start())
            }
            _ => panic!("Invalid get_aa_loc outer/inner pairing"),
        }
    }
}

impl Display for Location {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Complement(range) => write!(f, "c{}..{}", range.start(), range.end()),
            Self::Normal(range) => write!(f, "{}..{}", range.start(), range.end()),
            Self::SingleComplement(loc) => write!(f, "c{}", loc),
            Self::SingleNormal(loc) => write!(f, "{}", loc),
        }
    }
}

impl FromStr for Location {
    type Err = ParseIntError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(tail) = s.strip_prefix("complement(") {
            tail.trim_end_matches(')')
                .split_once("..")
                .map(|(start, end)| {
                    Ok(Self::Complement(
                        start.trim_start_matches('<').parse::<usize>()? - 1
                            ..=end.trim_start_matches('>').parse::<usize>()? - 1,
                    ))
                })
                .unwrap_or_else(|| Ok(Self::SingleComplement(tail.trim_end_matches(')').parse()?)))
        } else {
            s.split_once("..")
                .map(|(start, end)| {
                    Ok(Self::Normal(
                        start.trim_start_matches('<').parse::<usize>()? - 1
                            ..=end.trim_start_matches('>').parse::<usize>()? - 1,
                    ))
                })
                .unwrap_or_else(|| Ok(Self::SingleNormal(s.parse()?)))
        }
    }
}

impl DataItem {
    pub fn new(data: PreDataItem) -> Result<Self, String> {
        let mut result = Self {
            id: data.id[5..].split(';').next().unwrap().to_string(),
            species: data.os.ok_or("No species found")?,
            sequence: data.sq,
            genes: Vec::new(),
            regions: Vec::new(),
        };
        let mut current = None;
        let mut sequence = false;
        for line in data.ft {
            let line = &line[5..];
            if !line.starts_with(' ') || current.is_none() {
                if let Some(region) = current {
                    result.add_region(region);
                }
                let (key, location) = (&line[..data.ft_key_width], &line[data.ft_key_width..]);
                let location = location
                    .trim()
                    .parse()
                    .unwrap_or_else(|_| panic!("`{}` not a valid location", location));
                let seq = result.get_sequence(&location);
                current = Some(Region {
                    key: key.trim().to_string(),
                    location,
                    reported_seq: String::new(),
                    found_seq: seq,
                    allele: String::new(),
                    functional: false,
                    partial: false,
                });
                continue;
            }
            if let Some(current) = &mut current {
                let trimmed = line.trim();
                if sequence {
                    current.reported_seq = trimmed.trim_end_matches('\"').to_string();
                    if trimmed.ends_with('\"') {
                        sequence = false;
                    }
                } else if let Some(tail) = trimmed.strip_prefix("/translation=\"") {
                    current.reported_seq = tail.trim_end_matches('\"').to_string();
                    if !trimmed.ends_with('\"') {
                        sequence = true;
                    }
                } else if let Some(tail) = trimmed.strip_prefix("/IMGT_allele=\"") {
                    current.allele = tail.trim_end_matches('\"').to_string();
                } else if trimmed.starts_with("/functional") {
                    current.functional = true;
                } else if trimmed.starts_with("/partial") {
                    current.partial = true;
                }
            }
        }
        if let Some(region) = current {
            result.add_region(region);
        }

        Ok(result)
    }

    fn add_region(&mut self, region: Region) {
        if ["V-GENE", "C-GENE", "J-GENE"].contains(&region.key.as_str()) // , "D-GENE"
            && region.functional
            && !region.partial
            && region.allele.starts_with("IG")
        {
            self.genes.push(IMGTGene {
                key: region.key,
                location: region.location,
                allele: region.allele,
                regions: HashMap::new(),
            });
        } else if [
            "FR1-IMGT",
            "FR2-IMGT",
            "FR3-IMGT",
            "CDR1-IMGT",
            "CDR2-IMGT",
            "CDR3-IMGT",
            "1st-CYS",
            "2nd-CYS",
            "CONSERVED-TRP",
            "J-PHE",
            "CH1",
            "CH2",
            "CH3",
            "CH4",
            "CH5",
            "CH6",
            "CH7",
            "CH8",
            "CH9",
            "CHS",
            "J-REGION",
            //"D-REGION",
        ]
        .contains(&region.key.as_str())
        {
            if let Some(gene) = self
                .genes
                .iter_mut()
                .find(|g| g.location.contains(&region.location))
            {
                gene.regions.insert(region.key.clone(), region);
            } else {
                self.regions.push(region)
            }
        }
    }

    fn get_sequence(&self, slice: &Location) -> Option<(String, Vec<AminoAcid>)> {
        Some(translate(match slice {
            Location::Normal(range) => self.sequence.get(range.clone())?.to_string(),
            Location::SingleNormal(index) => {
                char::from(*self.sequence.as_bytes().get(*index)?).to_string()
            }
            Location::Complement(range) => {
                complement(self.sequence.get(range.clone())?.to_string())
            }
            Location::SingleComplement(index) => {
                complement(char::from(*self.sequence.as_bytes().get(*index)?).to_string())
            }
        }))
    }
}

fn complement(s: String) -> String {
    let map = HashMap::from([(b'a', b't'), (b't', b'a'), (b'c', b'g'), (b'g', b'c')]);
    String::from_utf8(
        s.as_bytes()
            .iter()
            .map(|c| {
                *map.get(c)
                    .unwrap_or_else(|| panic!("Invalid sequence: {} in `{s}`", char::from(*c)))
            })
            .rev()
            .collect(),
    )
    .unwrap()
}

fn translate(s: String) -> (String, Vec<AminoAcid>) {
    if s.len() < 3 {
        (s, Vec::new())
    } else {
        (
            s.clone(),
            (0..=s.len() - 3)
                .step_by(3)
                .filter_map(|chunk| {
                    AminoAcid::from_dna(&s[chunk..chunk + 3]).unwrap_or_else(|()| {
                        panic!("Not a valid codon: `{}`", &s[chunk..chunk + 3])
                    })
                })
                .collect(),
        )
    }
}

impl Display for DataItem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{}\t{}\n{}", self.id, self.species, self.sequence)?;
        for gene in &self.genes {
            writeln!(f, "G {gene}")?;
        }
        for region in &self.regions {
            writeln!(f, "R {region}")?;
        }
        Ok(())
    }
}

impl Display for IMGTGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{}\t{}\t{}", self.key, self.location, self.allele)?;
        for region in self.regions.values().sorted_by_key(|r| &r.key) {
            writeln!(f, "  R {region}")?;
        }
        Ok(())
    }
}

impl Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}",
            self.key,
            self.location,
            // self.sequence,
            // dna,
            // self.found_seq.0,
            self.found_seq
                .as_ref()
                .map(|seq| seq.1.iter().map(|a| a.char()).collect::<String>())
                .unwrap_or("<NO SEQ!>".to_string()),
        )
    }
}

impl IMGTGene {
    fn finish(self, species: Species) -> Result<Germline, String> {
        let get = |key| -> Result<Vec<AminoAcid>, String> {
            Ok(self
                .regions
                .get(key)
                .ok_or(format!("Could not find {key}"))?
                .found_seq
                .as_ref()
                .ok_or(format!("{key} does not have a sequence"))?
                // .unwrap_or_else(|| panic!("No sequence for `{key}` `{}`", self.allele))
                .1
                .clone())
        };
        let regions = if self.key == "V-GENE" {
            vec![
                (shared::Region::FR1, get("FR1-IMGT")?),
                (shared::Region::CDR1, get("CDR1-IMGT")?),
                (shared::Region::FR2, get("FR2-IMGT")?),
                (shared::Region::CDR2, get("CDR2-IMGT")?),
                (shared::Region::FR3, get("FR3-IMGT")?),
                (shared::Region::CDR3, get("CDR3-IMGT")?),
            ]
        } else if self.key == "C-GENE" {
            let mut seq = vec![
                (shared::Region::CH1, get("CH1")?),
                (shared::Region::CH2, get("CH2")?),
                (shared::Region::CH3, get("CH3")?),
            ];
            let mut possibly_add = |region, key: &str| -> Result<(), String> {
                if self.regions.contains_key(key) {
                    seq.push((
                        region,
                        self.regions
                            .get(key)
                            .ok_or(format!("Could not find {key}"))?
                            .found_seq
                            .as_ref()
                            .ok_or(format!("{key} does not have a sequence"))?
                            // .unwrap_or_else(|| panic!("No sequence for `{key}` `{}`", self.allele))
                            .1
                            .clone(),
                    ))
                }
                Ok(())
            };
            possibly_add(shared::Region::CH4, "CH4")?;
            possibly_add(shared::Region::CH5, "CH5")?;
            possibly_add(shared::Region::CH6, "CH6")?;
            possibly_add(shared::Region::CH7, "CH7")?;
            possibly_add(shared::Region::CH8, "CH8")?;
            possibly_add(shared::Region::CH9, "CH9")?;
            seq.push((shared::Region::CHS, get("CHS")?)); // TODO: what if only the combined CHX-CHS is present in the database
            seq
        } else if self.key == "J-GENE" {
            vec![(shared::Region::FR4, get("J-REGION")?)] // TODO: not fully correct right, has some CDR3 as well, and has quite some conserved residues
        } else if self.key == "D-GENE" {
            vec![(shared::Region::CDR3, get("D-REGION")?)]
        } else {
            Vec::new()
        };
        let sequence = regions.iter().flat_map(|r| r.1.clone()).collect();
        let regions = regions.iter().map(|r| (r.0, r.1.len())).collect();
        let conserved_map = HashMap::from([
            ("1st-CYS", Annotation::Cysteine1),
            ("2nd-CYS", Annotation::Cysteine2),
            ("CONSERVED-TRP", Annotation::Tryptophan),
            ("J-PHE", Annotation::Phenylalanine),
        ]);
        let conserved = self
            .regions
            .iter()
            .filter(|(key, _)| {
                ["1st-CYS", "2nd-CYS", "CONSERVED-TRP", "J-PHE"].contains(&key.as_str())
            })
            .map(|(key, _)| (conserved_map[key.as_str()], 0))
            .collect();
        let (name, allele) = Gene::from_imgt_name_with_allele(self.allele.as_str())?;
        Ok(Germline {
            name,
            alleles: vec![(
                allele,
                AnnotatedSequence {
                    sequence,
                    regions,
                    conserved,
                },
            )],
        })
    }
}
