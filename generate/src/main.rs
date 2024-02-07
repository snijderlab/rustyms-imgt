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
use itertools::Itertools;
use rustyms::{spectrum::AnnotatedPeak, AminoAcid};

fn main() {
    let file = File::open("../data/imgt.dat")
        .expect("Please provide the 'imgt.dat' file in the 'data' directory.");
    let mut output = BufWriter::new(File::create("../germlines/germlines.rs").unwrap());
    let mut docs = BufWriter::new(File::create("../germlines/germlines.md").unwrap());
    let mut error = BufWriter::new(File::create("errors.dat").unwrap());
    let data = parse_dat(BufReader::new(file));
    let mut grouped = HashMap::new();
    let mut errors = Vec::new();
    for element in data {
        let species = element.as_ref().unwrap().species;
        for gene in element.unwrap().genes {
            match gene.clone().finish() {
                Ok(gene) => grouped
                    .entry(species)
                    .or_insert(Germlines::new(species))
                    .insert(gene),
                Err(err) => {
                    errors.push((species, gene, err));
                }
            }
        }
        //writeln!(output, "{}", element.unwrap()).unwrap();
    }

    for (species, errors) in errors
        .iter()
        .map(|(species, gene, err)| (species, (gene, err)))
        .into_group_map()
        .into_iter()
        .map(|(species, genes)| {
            (
                species,
                genes
                    .into_iter()
                    .map(|(gene, err)| (gene.key.clone(), (gene, err)))
                    .into_group_map(),
            )
        })
    {
        writeln!(error, "SPECIES: {species}").unwrap();
        for (gene, errors) in errors {
            writeln!(error, "GENE: {gene}").unwrap();
            for (gene, err) in errors {
                writeln!(error, "ERROR FOR GENE:\n{species}\t{gene}\t{err}\n").unwrap();
            }
        }
    }

    writeln!(
        output,
        "#![allow(non_snake_case,non_upper_case_globals)]\nuse std::sync::OnceLock;\nuse crate::shared::{{Germlines, Species}};"
    )
    .unwrap();
    writeln!(output, "/// Get the germlines for any of the available species. See the main documentation for which species have which data available.").unwrap();
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
    // germlines
    writeln!(
        output,
        "pub fn germlines(species: Species) -> Option<&'static Germlines> {{match species {{"
    )
    .unwrap();

    for species in &found_species {
        writeln!(output, "Species::{0} => Some(lock_{0}()),", species.ident()).unwrap();
    }
    writeln!(output, "_=>None}}}}").unwrap();
    // all_germlines
    writeln!(
        output,
"/// Get all germlines in one iterator, see the main documentation for more information about the available germlines
pub fn all_germlines() -> impl std::iter::Iterator<Item = &'static Germlines> {{"
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
    // par_germlines
    writeln!(
        output,
"/// Get all germlines in one parallel iterator, see the main documentation for more information about the available germlines
#[cfg(feature = \"rayon\")]
use rayon::prelude::*;
#[cfg(feature = \"rayon\")]
pub fn par_germlines() -> impl rayon::prelude::ParallelIterator<Item = &'static Germlines> {{"
    )
    .unwrap();
    let mut first = true;
    for species in &found_species {
        if first {
            first = false;
            writeln!(output, "rayon::iter::once(lock_{}())", species.ident()).unwrap();
        } else {
            writeln!(
                output,
                ".chain(rayon::iter::once(lock_{}()))",
                species.ident()
            )
            .unwrap();
        }
    }
    writeln!(output, "}}").unwrap();

    for species in &found_species {
        writeln!(
            output,
"static LOCK_{0}: OnceLock<Germlines> = OnceLock::new();
fn lock_{0}()->&'static Germlines{{LOCK_{0}.get_or_init(|| {{bincode::deserialize(include_bytes!(\"{species}.bin\")).unwrap()}})}}",
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

#[derive(Clone, Debug, Eq, PartialEq)]
struct IMGTGene {
    key: String,
    location: Location,
    allele: String,
    regions: HashMap<String, Region>,
}

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
struct Region {
    key: String,
    location: Location,
    reported_seq: String,
    found_seq: Option<(String, Vec<AminoAcid>)>,
    allele: String,
    functional: bool,
    partial: bool,
    shift: usize,
    splice_aa: Option<AminoAcid>,
}

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
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

    fn get_aa_loc(&self, inner: &Self) -> Option<RangeInclusive<usize>> {
        if !self.contains(inner) {
            None
        } else {
            match (self, inner) {
                (Self::Complement(s), Self::Complement(o)) | (Self::Normal(s), Self::Normal(o)) => {
                    Some((o.start() - s.start()) / 3..=(o.end() - s.start()) / 3)
                }
                (Self::Normal(s), Self::SingleNormal(o))
                | (Self::Complement(s), Self::SingleComplement(o)) => {
                    Some((o - s.start()) / 3..=(o - s.start()) / 3)
                }
                _ => None,
            }
        }
    }

    /// Break the location around the given amino acid index in the location. If the position is outside the range or this location is a single it returns None.
    fn splice(&self, position: usize) -> Option<(Self, Self)> {
        match self {
            Self::Normal(s) => {
                let mid_point = *s.start() + position * 3;
                if mid_point >= *s.end() {
                    None
                } else {
                    Some((
                        Self::Complement((*s.start())..=mid_point),
                        Self::Complement(mid_point..=*s.end()),
                    ))
                }
            }
            Self::Complement(s) => {
                let mid_point = *s.end() - position * 3;
                if mid_point <= *s.start() {
                    None
                } else {
                    Some((
                        Self::Complement((*s.start())..=mid_point),
                        Self::Complement(mid_point..=*s.end()),
                    ))
                }
            }
            _ => None,
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
        let mut current: Option<Region> = None;
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
                current = Some(Region {
                    key: key.trim().to_string(),
                    location,
                    reported_seq: String::new(),
                    found_seq: None,
                    allele: String::new(),
                    functional: false,
                    partial: false,
                    shift: 0,
                    splice_aa: None,
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
                } else if let Some(tail) = trimmed.strip_prefix("/codon_start=") {
                    current.shift = tail
                        .parse::<usize>()
                        .map_err(|_| format!("Not a valid codon_start: '{tail}'"))?
                        - 1;
                } else if let Some(tail) = trimmed.strip_prefix("/splice-expectedcodon=") {
                    if let Some(i) = tail.find(']') {
                        current.splice_aa = AminoAcid::try_from(tail.as_bytes()[i - 1]).ok();
                    }
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

    fn add_region(&mut self, mut region: Region) {
        // Get the actual sequence
        region.found_seq = self.get_sequence(&region.location, region.shift);

        // Determine if what this region is and if is warrants keeping
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
            "FR4-IMGT",
            "CDR1-IMGT",
            "CDR2-IMGT",
            "CDR3-IMGT",
            "1st-CYS",
            "2nd-CYS",
            "CONSERVED-TRP",
            "J-REGION",
            // "J-TRP",
            // "J-PHE",
            "J-MOTIF",
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
            "H", //"D-REGION",
            "M",
            "M1",
            "M2",
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

    fn get_sequence(&self, slice: &Location, shift: usize) -> Option<(String, Vec<AminoAcid>)> {
        let (inner_shift, shift) = if shift == 2 { (1, 0) } else { (0, shift) };

        Some(translate(
            &match slice {
                Location::Normal(range) => {
                    if *range.start() < inner_shift {
                        return None;
                    }
                    self.sequence
                        .get(range.start() - inner_shift..=*range.end())?
                        .to_string()
                }
                Location::SingleNormal(index) => {
                    char::from(*self.sequence.as_bytes().get(*index)?).to_string()
                }
                Location::Complement(range) => complement(
                    self.sequence
                        .get(*range.start()..=*range.end() + inner_shift)?
                        .to_string(),
                ),
                Location::SingleComplement(index) => {
                    complement(char::from(*self.sequence.as_bytes().get(*index)?).to_string())
                }
            }[shift..],
        ))
        .map(|(s, v)| (s.to_owned(), v))
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

fn translate(s: &str) -> (&str, Vec<AminoAcid>) {
    if s.len() < 3 {
        (s, Vec::new())
    } else {
        (
            s,
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
        for region in self.regions.values().sorted_by_key(|reg| &reg.key) {
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
    fn finish(self) -> Result<Germline, String> {
        let get = |key| -> Result<(Vec<AminoAcid>, Location), String> {
            self.regions
                .get(key)
                .ok_or(format!("Could not find {key}"))
                .and_then(|region| {
                    region
                        .found_seq
                        .as_ref()
                        .map(|seq| {
                            let mut final_seq = region
                                .splice_aa
                                .map(|aa| vec![aa])
                                .filter(|_| region.shift != 2)
                                .unwrap_or_default();
                            final_seq.extend(seq.1.clone());
                            (final_seq, region.location.clone())
                        })
                        .ok_or(format!("{key} does not have a sequence"))
                })
        };
        let mut additional_annotations = Vec::new();
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
            let mut seq = Vec::new();
            let mut possibly_add = |region, key: &str| -> Result<(), String> {
                if self.regions.contains_key(key) {
                    seq.push((
                        region,
                        self.regions
                            .get(key)
                            .ok_or(format!("Could not find {key}"))
                            .and_then(|region| {
                                region
                                    .found_seq
                                    .as_ref()
                                    .map(|seq| {
                                        let mut final_seq = region
                                            .splice_aa
                                            .map(|aa| vec![aa])
                                            .filter(|_| region.shift != 2)
                                            .unwrap_or_default();
                                        final_seq.extend(seq.1.clone());
                                        (final_seq, region.location.clone())
                                    })
                                    .ok_or(format!("{key} does not have a sequence"))
                            })?,
                    ))
                }
                Ok(())
            };
            possibly_add(shared::Region::CH1, "CH1")?;
            possibly_add(shared::Region::H, "H")?;
            possibly_add(shared::Region::CH2, "CH2")?;
            possibly_add(shared::Region::CH3, "CH3")?;
            possibly_add(shared::Region::CH4, "CH4")?;
            possibly_add(shared::Region::CH5, "CH5")?;
            possibly_add(shared::Region::CH6, "CH6")?;
            possibly_add(shared::Region::CH7, "CH7")?;
            possibly_add(shared::Region::CH8, "CH8")?;
            possibly_add(shared::Region::CH9, "CH9")?;
            possibly_add(shared::Region::CHS, "CHS")?; // TODO: what if only the combined CHX-CHS is present in the database
                                                       // possibly_add(shared::Region::M, "M")?; // TODO: Figure out if support for membrane bound is needed, if so provide a way to switch between the two versions
                                                       // possibly_add(shared::Region::M1, "M1")?;
                                                       // possibly_add(shared::Region::M2, "M2")?;
            if seq.is_empty() {
                return Err("Empty C sequence".to_string());
            }
            seq
        } else if self.key == "J-GENE" {
            // if self.regions.contains_key("FR4-IMGT") {
            //     let fr4 = get("FR4-IMGT")?;
            //     let j = get("J-REGION")?;
            //     let cdr3_len = j.0.len() - fr4.0.len();
            //     let j = fix_j(j, cdr3_len);
            //     additional_annotations.extend(j.1);
            //     j.0
            // } else if self.regions.contains_key("J-MOTIF") {
            //     let motif = get("J-MOTIF")?;
            //     let j = get("J-REGION")?;
            //     let loc =
            //         j.1.get_aa_loc(&motif.1)
            //             .ok_or("J-MOTIF does not fall into J-REGION")?;
            //     let j = fix_j(j, *loc.start());
            //     additional_annotations.extend(j.1);
            //     j.0
            // } else {
            let j = get("J-REGION")?;
            let motif = j.0.iter().tuple_windows().position(|(a, b, _, d)| {
                (*a == AminoAcid::W || *a == AminoAcid::F)
                    && *b == AminoAcid::G
                    && *d == AminoAcid::G
            });
            if let Some(motif_start) = motif {
                let j = fix_j(j, motif_start);
                additional_annotations.extend(j.1);
                j.0
            } else {
                vec![(shared::Region::FR4, j)] // TODO: not fully correct right, has some CDR3 as well, and has quite some conserved residues
            }
            // }
        } else if self.key == "D-GENE" {
            vec![(shared::Region::CDR3, get("D-REGION")?)]
        } else {
            Vec::new()
        };
        let sequence: Vec<AminoAcid> = regions.iter().flat_map(|reg| reg.1 .0.clone()).collect();
        let region_lengths = regions.iter().map(|reg| (reg.0, reg.1 .0.len())).collect();
        let conserved_map = HashMap::from([
            ("1st-CYS", Annotation::Cysteine1),
            ("2nd-CYS", Annotation::Cysteine2),
            ("CONSERVED-TRP", Annotation::Tryptophan),
            ("J-PHE", Annotation::Phenylalanine),
            ("J-TRP", Annotation::Tryptophan),
        ]);
        let mut conserved = self
            .regions
            .iter()
            .filter(|(key, _)| {
                ["1st-CYS", "2nd-CYS", "CONSERVED-TRP", "J-PHE", "J-TRP"].contains(&key.as_str())
            })
            .map(|(key, region)| {
                find_aa_location(&region.location, &regions)
                    .map(|index| (conserved_map[key.as_str()], index))
                    .ok_or(format!("Cannot find location of '{key}' '{region}'"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        conserved.extend(
            find_possible_n_glycan_locations(&sequence)
                .iter()
                .map(|i| (Annotation::NGlycan, *i)),
        );
        conserved.extend(additional_annotations);
        let (name, allele) = Gene::from_imgt_name_with_allele(self.allele.as_str())?;
        Ok(Germline {
            name,
            alleles: vec![(
                allele,
                AnnotatedSequence {
                    sequence: sequence.into(),
                    regions: region_lengths,
                    conserved,
                },
            )],
        })
    }
}

fn find_aa_location(
    location: &Location,
    sections: &[(shared::Region, (Vec<AminoAcid>, Location))],
) -> Option<usize> {
    let mut start = 0;
    for section in sections {
        if let Some(index) = section.1 .1.get_aa_loc(location) {
            return Some(start + index.start());
        }
        start += section.1 .0.len();
    }
    None
}

fn find_possible_n_glycan_locations(sequence: &[AminoAcid]) -> Vec<usize> {
    let mut result = Vec::new();
    for (index, aa) in sequence.windows(3).enumerate() {
        if aa[0] == AminoAcid::N && (aa[2] == AminoAcid::S || aa[2] == AminoAcid::T) {
            result.push(index);
        }
    }
    result
}

fn fix_j(
    j: (Vec<AminoAcid>, Location),
    cdr3_length: usize,
) -> (
    Vec<(shared::Region, (Vec<AminoAcid>, Location))>,
    Vec<(Annotation, usize)>,
) {
    let (cdr3_loc, fr4_loc) =
        j.1.splice(cdr3_length)
            .expect("CDR3 should fit in full FR4 of J gene");
    let cdr3 = (j.0[..cdr3_length].to_vec(), cdr3_loc);
    let fr4 = (j.0[cdr3_length..].to_vec(), fr4_loc);

    let mut annotations = Vec::new();
    if fr4.0[0] == AminoAcid::W {
        annotations.push((Annotation::Tryptophan, cdr3_length));
    } else if fr4.0[0] == AminoAcid::F {
        annotations.push((Annotation::Phenylalanine, cdr3_length));
    }
    if fr4.0[1] == AminoAcid::G {
        annotations.push((Annotation::Glycine, cdr3_length + 1));
    }
    if fr4.0[3] == AminoAcid::G {
        annotations.push((Annotation::Glycine, cdr3_length + 3));
    }

    (
        vec![(shared::Region::CDR3, cdr3), (shared::Region::FR4, fr4)],
        annotations,
    )
}
