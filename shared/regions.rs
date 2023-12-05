use rustyms::LinearPeptide;
use serde::{Deserialize, Serialize};
use std::{fmt::Display, str::FromStr};

use super::species::Species;

/// A selection of germlines from a single species. Use the [`Self::get`] method to retrieve the sequences you are interested in.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct Germlines {
    pub(crate) species: Species,
    pub(crate) h: Chain,
    pub(crate) k: Chain,
    pub(crate) l: Chain,
    pub(crate) i: Chain,
}

impl Germlines {
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

impl<'a> IntoIterator for &'a Germlines {
    type IntoIter = std::array::IntoIter<(Kind, &'a Chain), 4>;
    type Item = (Kind, &'a Chain);

    fn into_iter(self) -> Self::IntoIter {
        [
            (Kind::Heavy, &self.h),
            (Kind::LightKappa, &self.k),
            (Kind::LightLambda, &self.l),
            (Kind::I, &self.i),
        ]
        .into_iter()
    }
}

#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(feature = "rayon")]
impl<'a> IntoParallelIterator for &'a Germlines {
    type Iter = rayon::array::IntoIter<(Kind, &'a Chain), 4>;
    type Item = (Kind, &'a Chain);

    fn into_par_iter(self) -> Self::Iter {
        [
            (Kind::Heavy, &self.h),
            (Kind::LightKappa, &self.k),
            (Kind::LightLambda, &self.l),
            (Kind::I, &self.i),
        ]
        .into_par_iter()
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub(crate) struct Chain {
    pub variable: Vec<Germline>,
    pub joining: Vec<Germline>,
    pub constant: Vec<Germline>,
}

impl Chain {
    pub(crate) fn insert(&mut self, mut germline: Germline) {
        let db = match &germline.name.segment {
            Segment::V => &mut self.variable,
            Segment::J => &mut self.joining,
            Segment::C(_) => &mut self.constant,
        };

        match db.binary_search_by_key(&germline.name, |g| g.name.clone()) {
            // If there are multiple copies of the same region keep the one with the most annotations + regions
            Ok(index) => {
                match db[index]
                    .alleles
                    .binary_search_by_key(&germline.alleles[0].0, |a| a.0)
                {
                    Ok(allele_index) => {
                        if germline.alleles[0].1.conserved.len()
                            + germline.alleles[0].1.regions.len()
                            > db[index].alleles[allele_index].1.conserved.len()
                                + db[index].alleles[allele_index].1.regions.len()
                        {
                            db[index].alleles[allele_index] = germline.alleles.pop().unwrap()
                        }
                    }
                    Err(allele_index) => db[index]
                        .alleles
                        .insert(allele_index, germline.alleles.pop().unwrap()),
                }
            }
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

impl<'a> IntoIterator for &'a Chain {
    type IntoIter = std::array::IntoIter<(Segment, &'a [Germline]), 3>;
    type Item = (Segment, &'a [Germline]);

    fn into_iter(self) -> Self::IntoIter {
        [
            (Segment::V, self.variable.as_slice()),
            (Segment::J, self.joining.as_slice()),
            (Segment::C(None), self.constant.as_slice()),
        ]
        .into_iter()
    }
}

#[cfg(feature = "rayon")]
impl<'a> IntoParallelIterator for &'a Chain {
    type Iter = rayon::array::IntoIter<(Segment, &'a [Germline]), 3>;
    type Item = (Segment, &'a [Germline]);

    fn into_par_iter(self) -> Self::Iter {
        [
            (Segment::V, self.variable.as_slice()),
            (Segment::J, self.joining.as_slice()),
            (Segment::C(None), self.constant.as_slice()),
        ]
        .into_par_iter()
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct Germline {
    pub name: Gene,
    pub alleles: Vec<(usize, AnnotatedSequence)>,
}

impl<'a> IntoIterator for &'a Germline {
    type IntoIter = std::slice::Iter<'a, (usize, AnnotatedSequence)>;
    type Item = &'a (usize, AnnotatedSequence);

    fn into_iter(self) -> Self::IntoIter {
        self.alleles.iter()
    }
}

#[cfg(feature = "rayon")]
impl<'a> IntoParallelIterator for &'a Germline {
    type Iter = rayon::slice::Iter<'a, (usize, AnnotatedSequence)>;
    type Item = &'a (usize, AnnotatedSequence);

    fn into_par_iter(self) -> Self::Iter {
        self.alleles.par_iter()
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct AnnotatedSequence {
    pub sequence: LinearPeptide,
    /// The different regions in the sequence, defined by their name and length
    pub regions: Vec<(Region, usize)>,
    /// 0 based locations of single amino acid annotations, overlapping with the regions defined above
    pub conserved: Vec<(Annotation, usize)>,
}

/// A germline gene name, broken up in its constituent parts.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct Gene {
    pub kind: Kind,
    pub segment: Segment,
    pub number: Option<usize>,
    pub family: Vec<(Option<usize>, String)>,
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
            write!(
                f,
                "{}{}",
                element.0.map(|i| i.to_string()).unwrap_or_default(),
                element.1
            )?;
            last_str = !element.1.is_empty();
            first = false;
        }
        Ok(())
    }
}

impl Gene {
    /// Get an IMGT name with allele, eg IGHV3-23*03
    pub fn from_imgt_name_with_allele(s: &str) -> Result<(Self, usize), String> {
        let s = s.split(" or ").next().unwrap(); // Just ignore double names
        let (gene, tail) = Self::from_imgt_name_internal(s)?;
        if tail.is_empty() {
            return Ok((gene, 1));
        }
        let allele = if let Some(tail) = tail.strip_prefix('*') {
            tail.parse()
                .map_err(|_| format!("Invalid allele spec: `{}`", &tail))
        } else {
            Err(format!("Invalid allele spec: `{tail}`"))
        }?;
        Ok((gene, allele))
    }

    /// Get an IMGT name, eg IGHV3-23
    pub fn from_imgt_name(s: &str) -> Result<Self, String> {
        Self::from_imgt_name_internal(s).map(|(gene, _)| gene)
    }

    fn from_imgt_name_internal(s: &str) -> Result<(Self, &str), String> {
        fn parse_name(s: &str) -> (Option<(Option<usize>, String)>, &str) {
            let num = s
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect::<String>();
            let tail = s
                .chars()
                .skip(num.len())
                .take_while(|c| c.is_ascii_alphabetic())
                .collect::<String>();
            let rest = &s[num.len() + tail.len()..];
            if num.is_empty() && tail.is_empty() {
                return (None, s);
            }
            let num = if num.is_empty() {
                None
            } else {
                Some(num.parse().unwrap())
            };
            (Some((num, tail)), rest)
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

            Ok((
                Self {
                    kind,
                    segment,
                    number,
                    family,
                },
                tail,
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
    O,
    // DD,
    // MD,
    T,
}

impl FromStr for Segment {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "V" => Ok(Self::V),
            "J" => Ok(Self::J),
            "C" => Ok(Self::C(None)),
            "A" => Ok(Self::C(Some(Constant::A))),
            // "DD" => Ok(Self::C(Some(Constant::DD))),
            // "MD" => Ok(Self::C(Some(Constant::MD))),
            "D" => Ok(Self::C(Some(Constant::D))),
            "E" => Ok(Self::C(Some(Constant::E))),
            "G" => Ok(Self::C(Some(Constant::G))),
            "M" => Ok(Self::C(Some(Constant::M))),
            "O" => Ok(Self::C(Some(Constant::O))),
            "T" => Ok(Self::C(Some(Constant::T))),
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
                Self::C(Some(Constant::O)) => "O",
                // Self::C(Some(Constant::DD)) => "DD",
                // Self::C(Some(Constant::MD)) => "MD",
                Self::C(Some(Constant::T)) => "T",
            }
        )
    }
}

/// Any region in a germline, eg FR1, CDR1
#[derive(Debug, Serialize, Deserialize, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Region {
    CDR1,
    CDR2,
    CDR3,
    FR1,
    FR2,
    FR3,
    FR4,
    CH1,
    H,
    CH2,
    CH3,
    CH4,
    CH5,
    CH6,
    CH7,
    CH8,
    CH9,
    CHS,
    M,
    M1,
    M2,
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
                Self::H => "H",
                Self::CH2 => "CH2",
                Self::CH3 => "CH3",
                Self::CH4 => "CH4",
                Self::CH5 => "CH5",
                Self::CH6 => "CH6",
                Self::CH7 => "CH7",
                Self::CH8 => "CH8",
                Self::CH9 => "CH9",
                Self::CHS => "CHS",
                Self::M => "M",
                Self::M1 => "M1",
                Self::M2 => "M2",
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
    NGlycan,
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
                Self::NGlycan => "NGly",
            }
        )
    }
}

#[test]
fn imgt_names() {
    assert_eq!(
        Gene::from_imgt_name_with_allele("IGHV3-23*03")
            .map(|(g, a)| (g.to_string(), a))
            .unwrap(),
        ("IGHV3-23".to_string(), 3)
    );
    assert_eq!(
        Gene::from_imgt_name_with_allele("IGKV6-d*01")
            .map(|(g, a)| (g.to_string(), a))
            .unwrap(),
        ("IGKV6-d".to_string(), 1)
    );
}
