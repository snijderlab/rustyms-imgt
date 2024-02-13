//! This crate handles parsing the [IMGT LIGM-DB database](https://www.imgt.org/) into structures compatible with rustyms.
//! It additionally stores all regions and annotations. There are two main ways of selecting germline(s), specified by name
//! [`get_germline`] or by building a query over the data [`Selection`].
//!
//! <details><summary>Data present per species</summary>
//!
#![doc = include_str!("../../germlines/germlines.md")]
//!
//! </details>
//!
//! ```
//! use rustyms_imgt::*;
//! let selection = Selection::default()
//!                           .species([Species::HomoSapiens])
//!                           .chain([ChainType::Heavy])
//!                           .gene([GeneType::V]);
//! let first = selection.germlines().next().unwrap();
//! assert_eq!(first.name(), "IGHV1-2*01");
//! ```

#![warn(clippy::all, clippy::pedantic, clippy::nursery, missing_docs)]
#![allow(
    clippy::must_use_candidate,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::wildcard_imports,
    clippy::module_name_repetitions,
    clippy::suboptimal_flops,
    clippy::too_many_lines
)]

mod fancy;
#[path = "../../germlines/germlines.rs"]
mod germlines;
mod itertools_extension;
mod select;
mod shared;

use itertools_extension::*;
use std::collections::HashSet;

pub use fancy::*;
use germlines::{all_germlines, germlines, par_germlines};
use itertools::Itertools;
use rustyms::{
    align::{AlignType, Alignment, OwnedAlignment},
    AminoAcid, LinearPeptide, Tolerance,
};
pub use select::*;
pub use shared::*;

/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
pub fn consecutive_align<const STEPS: u16>(
    sequence: &LinearPeptide,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species>>,
    chains: Option<HashSet<ChainType>>,
    allele: AlleleSelection,
    tolerance: Tolerance,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
) -> Vec<Vec<(Allele<'static>, OwnedAlignment)>> {
    assert!(genes.len() >= 2);
    let mut output: Vec<Vec<(Allele<'static>, OwnedAlignment)>> = Vec::with_capacity(genes.len());

    let mut prev = 0;
    for n in 0..genes.len() {
        let (left_sequence, use_species) = if n == 0 {
            (sequence.clone(), species.clone())
        } else {
            prev += output[n - 1][0].1.start_b() + output[n - 1][0].1.len_b();
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
                let alignment = rustyms::align::align::<STEPS>(
                    seq.sequence,
                    &left_sequence,
                    matrix,
                    tolerance,
                    genes[n].1,
                )
                .to_owned();
                (seq, alignment)
            })
            .k_largest_by(return_number, |a, b| a.1.cmp(&b.1))
            .collect_vec(),
        );
    }
    output
}

/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
#[cfg(feature = "rayon")]
pub fn par_consecutive_align<const STEPS: u16>(
    sequence: &LinearPeptide,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species>>,
    chains: Option<HashSet<ChainType>>,
    allele: AlleleSelection,
    tolerance: Tolerance,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
) -> Vec<Vec<(Allele<'static>, OwnedAlignment)>> {
    use rayon::iter::ParallelIterator;

    assert!(genes.len() >= 2);
    let mut output: Vec<Vec<(Allele<'static>, OwnedAlignment)>> = Vec::with_capacity(genes.len());

    let mut prev = 0;
    for n in 0..genes.len() {
        let (left_sequence, use_species) = if n == 0 {
            (sequence.clone(), species.clone())
        } else {
            prev += output[n - 1][0].1.start_b() + output[n - 1][0].1.len_b();
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
                let alignment = rustyms::align::align::<STEPS>(
                    seq.sequence,
                    &left_sequence,
                    matrix,
                    tolerance,
                    genes[n].1,
                );
                (seq, alignment.to_owned())
            })
            .collect::<Vec<_>>()
            .into_iter()
            .k_largest_by(return_number, |a, b| a.1.cmp(&b.1))
            .collect_vec(),
        );
    }
    output
}
