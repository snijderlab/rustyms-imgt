use itertools::Itertools;
use rustyms::{AminoAcid, LinearPeptide};
use serde::{Deserialize, Serialize};
use std::{fmt::Display, str::FromStr};

use super::species::Species;

/// A selection of germlines from a single species. Use the [`Self::get`] method to retrieve the sequences you are interested in.
#[derive(Serialize, Deserialize, Debug)]
pub struct Germlines {
    species: Species,
    pub(crate) h: Chain,
    pub(crate) k: Chain,
    pub(crate) l: Chain,
    pub(crate) i: Chain,
}

impl Germlines {
    /// Get the species for which this are the germlines
    pub fn species(&self) -> Species {
        self.species
    }

    pub(crate) fn new(species: Species) -> Self {
        Self {
            species,
            h: Chain::default(),
            k: Chain::default(),
            l: Chain::default(),
            i: Chain::default(),
        }
    }

    pub(crate) fn insert(&mut self, germline: Germline) {
        match &germline.name.kind {
            Kind::Heavy => self.h.insert(germline),
            Kind::LightKappa => self.k.insert(germline),
            Kind::LightLambda => self.l.insert(germline),
            Kind::I => self.i.insert(germline),
        };
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub(crate) struct Chain {
    pub variable: Vec<Germline>,
    pub joining: Vec<Germline>,
    pub constant: Vec<Germline>,
}

impl Chain {
    pub(crate) fn insert(&mut self, germline: Germline) {
        let db = match &germline.name.segment {
            Segment::V => &mut self.variable,
            Segment::J => &mut self.joining,
            Segment::C(_) => &mut self.constant,
        };

        match db.binary_search_by_key(&germline.name, |g| g.name.clone()) {
            Ok(index) => db[index].alleles.extend(germline.alleles),
            Err(index) => db.insert(index, germline),
        }
    }

    pub(crate) fn doc_row(&self) -> String {
        format!(
            "|{}/{}|{}/{}|{}/{}|",
            self.variable.len(),
            self.variable.iter().map(|g| g.alleles.len()).sum::<usize>(),
            self.joining.len(),
            self.joining.iter().map(|g| g.alleles.len()).sum::<usize>(),
            self.constant.len(),
            self.constant.iter().map(|g| g.alleles.len()).sum::<usize>(),
        )
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct Germline {
    pub name: Gene,
    pub alleles: Vec<(usize, AnnotatedSequence)>,
}

#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct AnnotatedSequence {
    pub sequence: LinearPeptide,
    pub regions: Vec<(Region, usize)>,
    pub conserved: Vec<(Annotation, usize)>,
}

/// A germline gene name, broken up in its constituent parts.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct Gene {
    pub kind: Kind,
    pub segment: Segment,
    pub number: Option<usize>,
    pub family: Vec<(usize, String)>,
}

impl Display for Gene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fn to_roman(n: usize) -> &'static str {
            [
                "0", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
            ][n]
        }

        write!(
            f,
            "IG{}{}{}{}",
            self.kind,
            self.segment,
            if let Some(n) = &self.number {
                format!("({})", to_roman(*n))
            } else {
                String::new()
            },
            if self.number.is_some() && !self.family.is_empty() {
                "-"
            } else {
                ""
            }
        )?;

        let mut first = true;
        let mut last_str = false;
        for element in &self.family {
            if !first && !last_str {
                write!(f, "-")?;
            }
            write!(f, "{}{}", element.0, element.1)?;
            last_str = !element.1.is_empty();
            first = false;
        }
        Ok(())
    }
}

impl Gene {
    pub fn from_imgt_name(s: &str) -> Result<(Self, usize), String> {
        fn parse_name(s: &str) -> (Option<(usize, String)>, &str) {
            let num = s
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect::<String>();
            if num.is_empty() {
                return (None, s);
            }
            let tail = s
                .chars()
                .skip(num.len())
                .take_while(|c| c.is_ascii_alphabetic())
                .collect::<String>();
            let rest = &s[num.len() + tail.len()..];
            (Some((num.parse().unwrap(), tail)), rest)
        }

        fn from_roman(s: &str) -> Option<usize> {
            match s {
                "I" => Some(1),
                "II" => Some(2),
                "III" => Some(3),
                "IV" => Some(4),
                "V" => Some(5),
                "VI" => Some(6),
                "VII" => Some(7),
                "VIII" => Some(8),
                "IX" => Some(9),
                "X" => Some(10),
                _ => None,
            }
        }

        if s.starts_with("IG") {
            let kind = s[2..3]
                .parse()
                .map_err(|()| format!("Invalid kind: `{}`", &s[2..3]))?;
            let segment = s[3..4]
                .parse()
                .map_err(|()| format!("Invalid segment: `{}`", &s[3..4]))?;
            let mut start = 4;
            let number = if &s[4..5] == "(" {
                let end = s[5..].find(')').ok_or(format!(
                    "Invalid segment number `{}` out of `{}`",
                    &s[4..],
                    s
                ))?;
                start += end + 2;
                Some(from_roman(&s[5..5 + end]).ok_or(format!(
                    "Invalid roman numeral (or too big) `{}`",
                    &s[5..5 + end]
                ))?)
            } else {
                None
            };
            let tail = &s[start..];
            let mut tail = tail.trim_start_matches('-');
            let mut family = Vec::new();
            while let (Some(branch), t) = parse_name(tail) {
                family.push(branch);
                tail = t.trim_start_matches('-');
            }

            let allele = if let Some(tail) = tail.strip_prefix('*') {
                tail.parse()
                    .map_err(|_| format!("Invalid allele spec: `{}`", &tail))
            } else {
                Err(format!("Invalid allele spec: `{tail}`"))
            }?;

            Ok((
                Self {
                    kind,
                    segment,
                    number,
                    family,
                },
                allele,
            ))
        } else {
            Err("Gene name does not start with IG")?
        }
    }
}

/// Any kind of germline
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Hash, Debug)]
pub enum Kind {
    Heavy = 0,
    LightKappa,
    LightLambda,
    /// Fish I kind
    I,
}

impl TryFrom<usize> for Kind {
    type Error = ();
    fn try_from(i: usize) -> Result<Self, Self::Error> {
        match i {
            0 => Ok(Self::Heavy),
            1 => Ok(Self::LightKappa),
            2 => Ok(Self::LightLambda),
            3 => Ok(Self::I),
            _ => Err(()),
        }
    }
}

impl FromStr for Kind {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "K" => Ok(Self::LightKappa),
            "L" => Ok(Self::LightLambda),
            "H" => Ok(Self::Heavy),
            "I" => Ok(Self::I),
            _ => Err(()),
        }
    }
}

impl Display for Kind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Heavy => "H",
                Self::LightKappa => "K",
                Self::LightLambda => "L",
                Self::I => "I",
            }
        )
    }
}

/// Any segment in a germline, eg variable, joining
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord, Clone, Hash, Copy)]
pub enum Segment {
    /// Variable
    V,
    /// Joining
    J,
    /// Constant, potentially with the type of constant given as well
    C(Option<Constant>),
}

/// Any type of constant segment
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord, Clone, Copy, Hash)]
pub enum Constant {
    A,
    D,
    E,
    G,
    M,
}

impl FromStr for Segment {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "V" => Ok(Self::V),
            "J" => Ok(Self::J),
            "C" => Ok(Self::C(None)),
            "A" => Ok(Self::C(Some(Constant::A))),
            "D" => Ok(Self::C(Some(Constant::D))),
            "E" => Ok(Self::C(Some(Constant::E))),
            "G" => Ok(Self::C(Some(Constant::G))),
            "M" => Ok(Self::C(Some(Constant::M))),
            _ => Err(()),
        }
    }
}

impl Display for Segment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::V => "V",
                Self::J => "J",
                Self::C(None) => "C",
                Self::C(Some(Constant::A)) => "A",
                Self::C(Some(Constant::D)) => "D",
                Self::C(Some(Constant::E)) => "E",
                Self::C(Some(Constant::G)) => "G",
                Self::C(Some(Constant::M)) => "M",
            }
        )
    }
}

/// Any region in a germline, eg FR1, CDR1
#[derive(Debug, Serialize, Deserialize, Copy, Clone)]
pub enum Region {
    CDR1,
    CDR2,
    CDR3,
    FR1,
    FR2,
    FR3,
    FR4,
    CH1,
    CH2,
    CH3,
    CH4,
    CH5,
    CH6,
    CH7,
    CH8,
    CH9,
    CHS,
}

impl Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::CDR1 => "CDR1",
                Self::CDR2 => "CDR2",
                Self::CDR3 => "CDR3",
                Self::FR1 => "FR1",
                Self::FR2 => "FR2",
                Self::FR3 => "FR3",
                Self::FR4 => "FR4",
                Self::CH1 => "CH1",
                Self::CH2 => "CH2",
                Self::CH3 => "CH3",
                Self::CH4 => "CH4",
                Self::CH5 => "CH5",
                Self::CH6 => "CH6",
                Self::CH7 => "CH7",
                Self::CH8 => "CH8",
                Self::CH9 => "CH9",
                Self::CHS => "CHS",
            }
        )
    }
}

/// Any annotation in a germline, eg conserved residues
#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
pub enum Annotation {
    Cysteine1,
    Cysteine2,
    Tryptophan,
    Phenylalanine,
}

impl Display for Annotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Cysteine1 => "Cys1",
                Self::Cysteine2 => "Cys2",
                Self::Tryptophan => "Trp",
                Self::Phenylalanine => "Phe",
            }
        )
    }
}
