use itertools::Itertools;
use rustyms::{AminoAcid, LinearPeptide};
use serde::{Deserialize, Serialize};
use std::{fmt::Display, str::FromStr};

use super::species::Species;

#[derive(Serialize, Deserialize)]
pub struct Germlines {
    pub species: Species,
    pub heavy_variable: Vec<Germline>,
    pub heavy_joining: Vec<Germline>,
    pub heavy_constant: Vec<Germline>,
    pub light_variable: Vec<Germline>,
    pub light_joining: Vec<Germline>,
}

impl Germlines {
    pub fn new(species: Species) -> Self {
        Self {
            species,
            heavy_variable: Vec::new(),
            heavy_joining: Vec::new(),
            heavy_constant: Vec::new(),
            light_variable: Vec::new(),
            light_joining: Vec::new(),
        }
    }

    pub fn insert(&mut self, germline: Germline) {
        let name = &germline.name;
        let db = match (name.kind, &name.segment) {
            (Kind::Heavy, Segment::V) => &mut self.heavy_variable,
            (Kind::Heavy, Segment::J) => &mut self.heavy_joining,
            (Kind::Heavy, Segment::C(_)) => &mut self.heavy_constant,
            (Kind::LightKappa | Kind::LightLambda, Segment::V) => &mut self.light_variable,
            (Kind::LightKappa | Kind::LightLambda, Segment::J) => &mut self.light_joining,
            _ => panic!("Unknown combination of kind + segment"),
        };

        match db.binary_search_by_key(&germline.name, |g| g.name.clone()) {
            Ok(index) => db[index].alleles.extend(germline.alleles),
            Err(index) => db.insert(index, germline),
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct Germline {
    pub name: Gene,
    pub alleles: Vec<(usize, Allele)>,
}

#[derive(Serialize, Deserialize)]
pub struct Allele {
    pub sequence: Vec<AminoAcid>,
    pub regions: Vec<(Region, usize)>,
    pub conserved: Vec<(Annotation, usize)>,
}

#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone)]
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
    pub fn from_str(s: &str) -> Result<(Self, usize), String> {
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

#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum Kind {
    Heavy,
    LightKappa,
    LightLambda,
}

impl FromStr for Kind {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "K" => Ok(Self::LightKappa),
            "L" => Ok(Self::LightLambda),
            "H" => Ok(Self::Heavy),
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
            }
        )
    }
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord, Clone)]
pub enum Segment {
    V,
    D,
    J,
    C(Option<Constant>),
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord, Clone, Copy)]
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
            //"D" => Ok(Self::D),
            "J" => Ok(Self::J),
            "C" => Ok(Self::C(None)),
            "A" => Ok(Self::C(Some(Constant::A))),
            "D" => Ok(Self::C(Some(Constant::D))), // TODO: How is this used??
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
                Self::D => "D",
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

#[derive(Clone, Copy, Serialize, Deserialize)]
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
