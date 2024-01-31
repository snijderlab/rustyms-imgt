#![doc = include_str!("../../germlines/germlines.md")]

mod fancy;
#[path = "../../germlines/germlines.rs"]
mod germlines;
mod select;
#[path = "../../shared/mod.rs"]
mod shared;

use std::collections::HashSet;

pub use fancy::*;
use germlines::{all_germlines, germlines, par_germlines};
use itertools::Itertools;
use rustyms::{
    align::{Alignment, Type},
    AminoAcid, LinearPeptide, Tolerance,
};
pub use select::*;
pub use shared::*;

/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
pub fn consecutive_align<const STEPS: usize>(
    sequence: LinearPeptide,
    genes: &[(GeneType, bool)],
    species: Option<HashSet<Species>>,
    chains: Option<HashSet<ChainType>>,
    allele: AlleleSelection,
    tolerance: Tolerance,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
) -> Vec<(Allele<'static>, Alignment)> {
    assert!(genes.len() >= 2);
    const GLOBAL_A_LEFT: Type = Type::new(true, true, true, false);
    let mut output: Vec<(Allele, Alignment)> = Vec::with_capacity(genes.len());

    for n in 1..genes.len() {
        let left_sequence = if n == 0 {
            sequence.clone()
        } else {
            let mut left_sequence: LinearPeptide = sequence
                .clone()
                .sequence
                .into_iter()
                .skip(output[n - 1].1.start_b + output[n - 1].1.len_b())
                .collect();
            left_sequence.c_term = sequence.c_term.clone();
            left_sequence
        };

        if left_sequence.is_empty() {
            break;
        }

        output[n] = Selection {
            species: species.clone(),
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
                    if genes[n].1 {
                        GLOBAL_A_LEFT
                    } else {
                        Type::GLOBAL_A
                    },
                ),
            )
        })
        .max_by(|a, b| b.1.normalised_score.total_cmp(&a.1.normalised_score))
        .unwrap_or_else(|| {
            panic!(
                "Germline sequences should be available for gene {} in consecutive alignment",
                genes[n].0
            )
        });
    }
    output
}

/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
#[cfg(feature = "rayon")]
pub fn par_consecutive_align<const STEPS: usize>(
    sequence: LinearPeptide,
    genes: &[(GeneType, bool)],
    species: Option<HashSet<Species>>,
    chains: Option<HashSet<ChainType>>,
    allele: AlleleSelection,
    tolerance: Tolerance,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
) -> Vec<(Allele<'static>, Alignment)> {
    use rayon::iter::ParallelIterator;

    assert!(genes.len() >= 2);
    const GLOBAL_A_LEFT: Type = Type::new(true, true, true, false);
    let mut output: Vec<(Allele, Alignment)> = Vec::with_capacity(genes.len());

    for n in 1..genes.len() {
        let left_sequence = if n == 0 {
            sequence.clone()
        } else {
            let mut left_sequence: LinearPeptide = sequence
                .clone()
                .sequence
                .into_iter()
                .skip(output[n - 1].1.start_b + output[n - 1].1.len_b())
                .collect();
            left_sequence.c_term = sequence.c_term.clone();
            left_sequence
        };

        if left_sequence.is_empty() {
            break;
        }

        output[n] = Selection {
            species: species.clone(),
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
                    if genes[n].1 {
                        GLOBAL_A_LEFT
                    } else {
                        Type::GLOBAL_A
                    },
                ),
            )
        })
        .max_by(|a, b| b.1.normalised_score.total_cmp(&a.1.normalised_score))
        .unwrap_or_else(|| {
            panic!(
                "Germline sequences should be available for gene {} in consecutive alignment",
                genes[n].0
            )
        });
    }
    output
}
