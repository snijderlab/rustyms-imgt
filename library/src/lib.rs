#![doc = include_str!("../../germlines/germlines.md")]

#[path = "../../germlines/germlines.rs"]
mod germlines;
mod select;
#[path = "../../shared/mod.rs"]
mod shared;

use germlines::{all_germlines, germlines};
pub use select::*;
pub use shared::*;
