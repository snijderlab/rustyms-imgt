#![doc = include_str!("../../germlines/germlines.md")]

mod fancy;
#[path = "../../germlines/germlines.rs"]
mod germlines;
mod itertools_extension;
mod select;
#[path = "../../shared/mod.rs"]
mod shared;

use itertools_extension::*;
use ordered_float::OrderedFloat;
use std::collections::HashSet;

pub use fancy::*;
use germlines::{all_germlines, germlines, par_germlines};
use itertools::Itertools;
use rustyms::{
    align::{AlignType, Alignment},
    AminoAcid, LinearPeptide, Tolerance,
};
pub use select::*;
pub use shared::*;

/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
pub fn consecutive_align<const STEPS: u16>(
    sequence: LinearPeptide,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species>>,
    chains: Option<HashSet<ChainType>>,
    allele: AlleleSelection,
    tolerance: Tolerance,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
) -> Vec<Vec<(Allele<'static>, Alignment)>> {
    assert!(genes.len() >= 2);
    let mut output: Vec<Vec<(Allele<'static>, Alignment)>> = Vec::with_capacity(genes.len());

    let mut prev = 0;
    for n in 0..genes.len() {
        let (left_sequence, use_species) = if n == 0 {
            (sequence.clone(), species.clone())
        } else {
            prev += output[n - 1][0].1.start_b + output[n - 1][0].1.len_b();
            let mut left_sequence: LinearPeptide =
                sequence.clone().sequence.into_iter().skip(prev).collect();
            left_sequence.c_term = sequence.c_term.clone();
            (left_sequence, Some([output[n - 1][0].0.species].into()))
        };

        if left_sequence.is_empty() {
            break;
        }

        output.push(
            Selection {
                species: use_species,
                chains: chains.clone(),
                allele: allele.clone(),
                genes: Some([genes[n].0].into()),
            }
            .germlines()
            .map(|seq| {
                let a = seq.sequence.clone();
                (
                    seq,
                    rustyms::align::align::<STEPS>(
                        a,
                        left_sequence.clone(),
                        matrix,
                        tolerance,
                        genes[n].1,
                    ),
                )
            })
            .k_largest_by_key(return_number, |i| OrderedFloat(i.1.normalised_score))
            .collect_vec(),
        );
    }
    output
}

/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
#[cfg(feature = "rayon")]
pub fn par_consecutive_align<const STEPS: u16>(
    sequence: LinearPeptide,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species>>,
    chains: Option<HashSet<ChainType>>,
    allele: AlleleSelection,
    tolerance: Tolerance,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
) -> Vec<Vec<(Allele<'static>, Alignment)>> {
    use rayon::iter::ParallelIterator;

    assert!(genes.len() >= 2);
    let mut output: Vec<Vec<(Allele<'static>, Alignment)>> = Vec::with_capacity(genes.len());

    let mut prev = 0;
    for n in 0..genes.len() {
        let (left_sequence, use_species) = if n == 0 {
            (sequence.clone(), species.clone())
        } else {
            prev += output[n - 1][0].1.start_b + output[n - 1][0].1.len_b();
            let mut left_sequence: LinearPeptide =
                sequence.clone().sequence.into_iter().skip(prev).collect();
            left_sequence.c_term = sequence.c_term.clone();
            (left_sequence, Some([output[n - 1][0].0.species].into()))
        };

        if left_sequence.is_empty() {
            break;
        }

        output.push(
            Selection {
                species: use_species,
                chains: chains.clone(),
                allele: allele.clone(),
                genes: Some([genes[n].0].into()),
            }
            .par_germlines()
            .map(|seq| {
                let a = seq.sequence.clone();
                (
                    seq,
                    rustyms::align::align::<STEPS>(
                        a,
                        left_sequence.clone(),
                        matrix,
                        tolerance,
                        genes[n].1,
                    ),
                )
            })
            .collect::<Vec<_>>()
            .into_iter()
            .k_largest_by_key(return_number, |i| OrderedFloat(i.1.normalised_score))
            .collect_vec(),
        );
    }
    output
}
