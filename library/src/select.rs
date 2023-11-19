#[cfg(feature = "rayon")]
use rayon::prelude::*;
use std::collections::HashSet;

use rustyms::LinearPeptide;

pub use crate::shared::*;

/// Get a specific germline
pub fn get_germline(
    species: Species,
    gene: Gene,
    allele: Option<usize>,
) -> Option<Allele<'static>> {
    crate::germlines(species).and_then(|g| g.find(species, gene, allele))
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// The selection rules for iterating over all alleles in a selection of germlines.
pub struct Selection {
    /// The species you want, None allows all, otherwise only the species specified will be returned
    pub species: Option<HashSet<Species>>,
    /// The kind of alleles you want, None allows all, otherwise only the kinds specified will be returned
    pub kinds: Option<HashSet<Kind>>,
    /// The kind of segments you want, None allows all, otherwise only the segments specified will be returned
    pub segments: Option<HashSet<Segment>>,
    /// The way of handling alleles you want
    pub allele: AlleleSelection,
}

impl Selection {
    /// Builder pattern method to add a species selection, will replace any previously set selection
    pub fn species(self, species: impl Into<HashSet<Species>>) -> Self {
        Self {
            species: Some(species.into()),
            ..self
        }
    }
    /// Builder pattern method to add a kind selection, will replace any previously set selection
    pub fn kind(self, kinds: impl Into<HashSet<Kind>>) -> Self {
        Self {
            kinds: Some(kinds.into()),
            ..self
        }
    }
    /// Builder pattern method to add a segment selection, will replace any previously set selection
    pub fn segment(self, segments: impl Into<HashSet<Segment>>) -> Self {
        Self {
            segments: Some(segments.into()),
            ..self
        }
    }
    /// Builder pattern method to add an allele selection, will replace any previously set selection
    pub fn allele(self, allele: AlleleSelection) -> Self {
        Self { allele, ..self }
    }
    /// Get the selected alleles
    pub fn germlines(&self) -> impl Iterator<Item = Allele<'_>> {
        crate::all_germlines()
            .filter(|g| {
                self.species
                    .as_ref()
                    .map(|s| s.contains(&g.species))
                    .unwrap_or(true)
            })
            .flat_map(|g| g.into_iter().map(|c| (g.species, c.0, c.1)))
            .filter(|(_, kind, _)| {
                self.kinds
                    .as_ref()
                    .map(|k| k.contains(kind))
                    .unwrap_or(true)
            })
            .flat_map(|(species, _, c)| c.into_iter().map(move |g| (species, g.0, g.1)))
            .filter(|(_, segment, _)| {
                self.segments
                    .as_ref()
                    .map(|s| s.contains(segment))
                    .unwrap_or(true)
            })
            .flat_map(|(species, _, germlines)| germlines.iter().map(move |a| (species, a)))
            .flat_map(move |(species, germline)| {
                germline
                    .into_iter()
                    .take(self.allele.take_num())
                    .map(move |(a, seq)| (species, &germline.name, *a, seq))
            })
            .map(Into::into)
    }
    #[cfg(feature = "rayon")]
    /// Get the selected alleles in parallel fashion, only available if you enable the feature "rayon" (on by default)
    pub fn par_germlines(&self) -> impl ParallelIterator<Item = Allele<'_>> {
        crate::par_germlines()
            .filter(|g| {
                self.species
                    .as_ref()
                    .map(|s| s.contains(&g.species))
                    .unwrap_or(true)
            })
            .flat_map(|g| g.into_par_iter().map(|c| (g.species, c.0, c.1)))
            .filter(|(_, kind, _)| {
                self.kinds
                    .as_ref()
                    .map(|k| k.contains(kind))
                    .unwrap_or(true)
            })
            .flat_map(|(species, _, c)| c.into_par_iter().map(move |g| (species, g.0, g.1)))
            .filter(|(_, segment, _)| {
                self.segments
                    .as_ref()
                    .map(|s| s.contains(segment))
                    .unwrap_or(true)
            })
            .flat_map(|(species, _, germlines)| {
                germlines.into_par_iter().map(move |a| (species, a))
            })
            .flat_map(move |(species, germline)| {
                germline
                    .into_par_iter()
                    .take(self.allele.take_num())
                    .map(move |(a, seq)| (species, &germline.name, *a, seq))
            })
            .map(Into::into)
    }
}

impl Default for Selection {
    /// Get a default selection, which gives all kinds and segments but only returns the first allele
    fn default() -> Self {
        Self {
            species: None,
            kinds: None,
            segments: None,
            allele: AlleleSelection::First,
        }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
/// The way of handling alleles you want
pub enum AlleleSelection {
    /// Return all alleles
    All,
    /// Only return the first allele. It can have a number higher than 1 if the previous alleles are not functional.
    First,
}

impl AlleleSelection {
    fn take_num(&self) -> usize {
        match self {
            Self::First => 1,
            Self::All => usize::MAX,
        }
    }
}

/// A returned allele
#[non_exhaustive]
#[derive(Debug)]
pub struct Allele<'a> {
    /// The species where this gene originates from
    pub species: Species,
    /// The gene where this is the sequence for, eg `IGHV3-23`
    pub gene: std::borrow::Cow<'a, Gene>,
    /// The allele number, in IMGT this follows the name, eg `*01` is the allele in `IGHV3-23*01`
    pub allele: usize,
    /// The actual sequence, the sequences are flat amino acids, no modification of any way shape or form
    pub sequence: &'a LinearPeptide,
    /// The regions in the sequence, every region has an annotation and a length, all lengths together are the same length as the full sequence
    pub regions: &'a [(Region, usize)],
    /// Any additional annotations, every annotation has beside the kind it is also it location, as index in the sequence
    pub annotations: &'a [(Annotation, usize)],
}

impl<'a> Allele<'a> {
    /// Get the IMGT name for this allele
    pub fn name(&self) -> String {
        format!("{}*{:02}", self.gene, self.allele)
    }

    /// Get the region for a specific index into the sequence, None if outside range,
    /// the additional bool indicates if this is the starting position for the region
    pub fn region(&self, index: usize) -> Option<(Region, bool)> {
        let mut left = index;
        let mut regions_index = 0;
        let mut next = self.regions[regions_index];
        while left > next.1 {
            left -= next.1;
            regions_index += 1;
            if regions_index == self.regions.len() {
                return None;
            }
            next = self.regions[regions_index];
        }
        Some((next.0, left == 1))
    }

    /// Get all annotations for this position
    pub fn annotations(&self, index: usize) -> impl Iterator<Item = Annotation> + 'a {
        self.annotations
            .iter()
            .filter(move |a| a.1 == index)
            .map(|a| a.0)
    }
}

impl<'a> From<(Species, &'a Gene, usize, &'a AnnotatedSequence)> for Allele<'a> {
    fn from(value: (Species, &'a Gene, usize, &'a AnnotatedSequence)) -> Self {
        Self {
            species: value.0,
            gene: std::borrow::Cow::Borrowed(value.1),
            allele: value.2,
            sequence: &value.3.sequence,
            regions: &value.3.regions,
            annotations: &value.3.conserved,
        }
    }
}

impl Germlines {
    pub fn find(&self, species: Species, gene: Gene, allele: Option<usize>) -> Option<Allele<'_>> {
        let chain = match gene.kind {
            Kind::Heavy => &self.h,
            Kind::LightKappa => &self.k,
            Kind::LightLambda => &self.l,
            Kind::I => &self.i,
        };
        let genes = match gene.segment {
            Segment::V => &chain.variable,
            Segment::J => &chain.joining,
            Segment::C(_) => &chain.constant,
        };
        genes
            .binary_search_by(|g| g.name.cmp(&gene))
            .ok()
            .and_then(|g| {
                let g = &genes[g];
                allele
                    .map(|a| g.alleles.iter().find(|(ga, _)| a == *ga))
                    .unwrap_or(g.alleles.first())
            })
            .map(move |(a, seq)| Allele {
                species,
                gene: std::borrow::Cow::Owned(gene),
                allele: *a,
                sequence: &seq.sequence,
                regions: &seq.regions,
                annotations: &seq.conserved,
            })
    }
}

#[cfg(test)]
mod tests {
    use crate::Selection;
    use crate::{Kind, Segment, Species};

    #[test]
    fn try_first_human() {
        let selection = Selection::default()
            .species([Species::HomoSapiens])
            .kind([Kind::Heavy])
            .segment([Segment::V]);
        let first = selection.germlines().next().unwrap();
        assert_eq!(first.name(), "IGHV1-2*01");
    }
}
