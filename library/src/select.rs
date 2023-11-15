use std::collections::HashSet;

pub use crate::shared::*;

pub struct Selection {
    pub kind: Option<HashSet<Kind>>,
    pub segment: Option<HashSet<Segment>>,
    pub allele: AlleleSelection,
}

pub enum AlleleSelection {
    All,
    First,
}

#[non_exhaustive]
#[derive(Debug)]
pub struct FoundAllele<'a> {
    pub gene: &'a Gene,
    pub allele: usize,
    pub sequence: &'a AnnotatedSequence,
}

impl Germlines {
    pub fn get<'a>(
        &'static self,
        selection: &'a Selection,
    ) -> Option<impl Iterator<Item = FoundAllele<'a>> + 'a> {
        GermlinesIterator::new(selection, self)
    }

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
            if let Some(kinds) = &selection.kind {
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
    type Item = FoundAllele<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut next = self.iter.next();
        while next.is_none() {
            if self.selection.kind.is_none() {
                self.iter = ChainsIterator::new(self.selection, self.germlines.index(self.index))?;
                self.index = Kind::try_from(self.index as usize + 1).ok()?;
            } else if let Some(kinds) = &self.selection.kind {
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
            if let Some(segments) = &selection.segment {
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
    type Item = FoundAllele<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut next = self.iter.next();
        while next.is_none() {
            if self.selection.segment.is_none() {
                self.iter = self.chain.index(self.selection, self.index);
                self.index = self.index.next()?;
            } else if let Some(segments) = &self.selection.segment {
                while let Some(index) = self.index.next() {
                    self.index = index;
                    if segments.contains(&index) {
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
    type Item = FoundAllele<'a>;
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
            Some(FoundAllele {
                gene: &germline.name,
                allele: res.0,
                sequence: &res.1,
            })
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use crate::germlines;
    use crate::Selection;
    use crate::{Kind, Segment, Species};

    #[test]
    fn try_first_human() {
        let germlines = germlines(Species::HomoSapiens).unwrap();
        let selection = Selection {
            kind: Some(HashSet::from([Kind::Heavy])),
            segment: Some(HashSet::from([Segment::V])),
            allele: crate::AlleleSelection::First,
        };
        let first = germlines.get(&selection).unwrap().next().unwrap();
        assert_eq!(first.gene.to_string(), "IGHV1-2");
        assert_eq!(first.allele, 2);
    }
}
