use std::collections::HashSet;

use rustyms::LinearPeptide;

pub use crate::shared::*;

/// The selection rules for iterating over all alleles in a selection of germlines.
pub struct Selection {
    /// The kind of alleles you want, None allows all, otherwise only the kinds specified will be returned
    pub kinds: Option<HashSet<Kind>>,
    /// The kind of segments you want, None allows all, otherwise only the segments specified will be returned
    pub segments: Option<HashSet<Segment>>,
    /// The way of handling alleles you want
    pub allele: AlleleSelection,
}

impl Selection {
    /// Builder pattern method to add a kind selection
    pub fn kind(self, kinds: impl Into<HashSet<Kind>>) -> Self {
        Self {
            kinds: Some(kinds.into()),
            ..self
        }
    }
    /// Builder pattern method to add a segment selection
    pub fn segment(self, segments: impl Into<HashSet<Segment>>) -> Self {
        Self {
            segments: Some(segments.into()),
            ..self
        }
    }
    /// Builder pattern method to add an allele selection
    pub fn allele(self, allele: AlleleSelection) -> Self {
        Self { allele, ..self }
    }
}

impl Default for Selection {
    /// Get a default selection, which gives all kinds and segments but only returns the first allele
    fn default() -> Self {
        Self {
            kinds: None,
            segments: None,
            allele: AlleleSelection::First,
        }
    }
}

/// The way of handling alleles you want
pub enum AlleleSelection {
    /// Return all alleles
    All,
    /// Only return the first allele. It can have a number higher than 1 if the previous alleles are not functional.
    First,
}

/// A returned allele
#[non_exhaustive]
#[derive(Debug)]
pub struct Allele<'a> {
    /// The gene where this is the sequence for, eg `IGHV3-23`
    pub gene: &'a Gene,
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
    pub fn name(&self) -> String {
        format!("{}*{:02}", self.gene, self.allele)
    }
}

impl Germlines {
    /// Retrieve the actual sequences fitting to the given selection rules.
    pub fn get<'a>(
        &'static self,
        selection: &'a Selection,
    ) -> Option<impl Iterator<Item = Allele<'a>> + 'a> {
        GermlinesIterator::new(selection, self)
    }

    // TODO: Make a find method which finds a specific germline (does binary search on the name and then either gives the first allele or a specific one)

    fn index(&self, index: Kind) -> &Chain {
        match index {
            Kind::Heavy => &self.h,
            Kind::LightKappa => &self.k,
            Kind::LightLambda => &self.l,
            Kind::I => &self.i,
        }
    }
}

impl Chain {
    fn index<'a>(&'a self, selection: &'a Selection, index: Segment) -> AllelesIterator<'a> {
        AllelesIterator {
            selection,
            index: 0,
            elements: match index {
                Segment::V => &self.variable,
                Segment::J => &self.joining,
                Segment::C(_) => &self.constant,
            },
            sub_index: 0,
            sub_elements: None,
        }
    }
}

impl Kind {
    fn next(&self) -> Option<Self> {
        match self {
            Self::Heavy => Some(Self::LightKappa),
            Self::LightKappa => Some(Self::LightLambda),
            Self::LightLambda => Some(Self::I),
            Self::I => None,
        }
    }
}

struct GermlinesIterator<'a> {
    selection: &'a Selection,
    germlines: &'static Germlines,
    iter: ChainsIterator<'a>,
    index: Kind,
}

impl<'a> GermlinesIterator<'a> {
    fn new(selection: &'a Selection, germlines: &'static Germlines) -> Option<Self> {
        let mut index = Kind::Heavy;
        let (iter, index) = (|| {
            if let Some(kinds) = &selection.kinds {
                if kinds.contains(&index) {
                    return Some((
                        ChainsIterator::new(selection, germlines.index(index))?,
                        index,
                    ));
                }
                while let Some(i) = index.next() {
                    index = i;
                    if kinds.contains(&index) {
                        return Some((
                            ChainsIterator::new(selection, germlines.index(index))?,
                            index,
                        ));
                    }
                }
                None
            } else {
                Some((
                    ChainsIterator::new(selection, germlines.index(index))?,
                    index.next().unwrap_or(Kind::I),
                ))
            }
        })()?;
        Some(Self {
            selection,
            germlines,
            iter,
            index,
        })
    }
}

impl<'a> Iterator for GermlinesIterator<'a> {
    type Item = Allele<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut next = self.iter.next();
        while next.is_none() {
            if self.selection.kinds.is_none() {
                self.iter = ChainsIterator::new(self.selection, self.germlines.index(self.index))?;
                self.index = Kind::try_from(self.index as usize + 1).ok()?;
            } else if let Some(kinds) = &self.selection.kinds {
                for index in self.index as usize..=Kind::I as usize {
                    let index = Kind::try_from(index).ok()?;
                    if kinds.contains(&index) {
                        self.iter =
                            ChainsIterator::new(self.selection, self.germlines.index(self.index))?;
                        self.index =
                            Kind::try_from((index as usize + 1).max(Kind::I as usize)).ok()?; // +1
                        break;
                    }
                }
                return None;
            }
            next = self.iter.next();
        }
        next
    }
}

impl Segment {
    fn next(&self) -> Option<Self> {
        match self {
            Self::V => Some(Self::J),
            Self::J => Some(Self::C(None)),
            Self::C(_) => None,
        }
    }
}

struct ChainsIterator<'a> {
    selection: &'a Selection,
    chain: &'static Chain,
    iter: AllelesIterator<'a>,
    index: Segment,
}

impl<'a> ChainsIterator<'a> {
    fn new(selection: &'a Selection, chain: &'static Chain) -> Option<Self> {
        let mut index = Segment::V;
        let (iter, index) = (|| {
            if let Some(segments) = &selection.segments {
                if segments.contains(&index) {
                    return Some((chain.index(selection, index), index));
                }
                while let Some(i) = index.next() {
                    index = i;
                    if segments.contains(&index) {
                        return Some((chain.index(selection, index), index));
                    }
                }
                None
            } else {
                Some((
                    chain.index(selection, index),
                    index.next().unwrap_or(Segment::C(None)),
                ))
            }
        })()?;
        Some(Self {
            selection,
            chain,
            iter,
            index,
        })
    }
}

impl<'a> Iterator for ChainsIterator<'a> {
    type Item = Allele<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut next = self.iter.next();
        while next.is_none() {
            if self.selection.segments.is_none() {
                self.iter = self.chain.index(self.selection, self.index);
                self.index = self.index.next()?;
            } else if let Some(segments) = &self.selection.segments {
                while let Some(index) = self.index.next() {
                    self.index = index;
                    if segments.contains(&index) {
                        // TODO: Allow Constant(None) to select all constant chains
                        self.iter = self.chain.index(self.selection, self.index);
                        break;
                    }
                }
                return None;
            }
            next = self.iter.next();
        }
        next
    }
}

struct AllelesIterator<'a> {
    selection: &'a Selection,
    index: usize,
    elements: &'a [Germline],
    sub_index: usize,
    sub_elements: Option<(&'a Germline, &'a [(usize, AnnotatedSequence)])>,
}

impl<'a> Iterator for AllelesIterator<'a> {
    type Item = Allele<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.elements.len() {
            return None;
        }
        if self.sub_elements.is_none() || self.sub_index >= self.sub_elements.unwrap().1.len() {
            self.sub_elements = match self.selection.allele {
                AlleleSelection::First => Some((
                    &self.elements[self.index],
                    &self.elements[self.index].alleles[0..1],
                )),
                AlleleSelection::All => Some((
                    &self.elements[self.index],
                    &self.elements[self.index].alleles[..],
                )),
            };
            self.sub_index = 0;
            self.index += 1;
        }

        if let Some((germline, sub)) = self.sub_elements {
            let res = &sub[self.sub_index];
            self.sub_index += 1;
            Some(Allele {
                gene: &germline.name,
                allele: res.0,
                sequence: &res.1.sequence,
                regions: &res.1.regions,
                annotations: &res.1.conserved,
            })
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::germlines;
    use crate::Selection;
    use crate::{Kind, Segment, Species};

    #[test]
    fn try_first_human() {
        let germlines = germlines(Species::HomoSapiens).unwrap();
        let selection = Selection::default()
            .kind([Kind::Heavy])
            .segment([Segment::V]);
        let first = germlines.get(&selection).unwrap().next().unwrap();
        assert_eq!(first.name(), "IGHV1-2*01");
    }
}
