#![doc = include_str!("../../germlines/germlines.md")]

mod fancy;
#[path = "../../germlines/germlines.rs"]
mod germlines;
mod select;
#[path = "../../shared/mod.rs"]
mod shared;

pub use fancy::*;
use germlines::{all_germlines, germlines, par_germlines};
pub use select::*;
pub use shared::*;
