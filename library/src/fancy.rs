pub use crate::shared::*;
use std::fmt::Display;
use std::fmt::Write;

pub trait FancyDisplay: Display {
    fn to_fancy_string(&self) -> String;
}

impl FancyDisplay for Gene {
    fn to_fancy_string(&self) -> String {
        let mut f = String::new();
        fn to_roman(n: usize) -> &'static str {
            ["0", "Ⅰ", "Ⅱ", "Ⅲ", "Ⅳ", "Ⅴ", "Ⅵ", "Ⅶ", "Ⅷ", "Ⅸ", "Ⅹ"][n]
        }

        write!(
            f,
            "IG{}{}{}{}",
            self.chain.to_fancy_string(),
            self.gene.to_fancy_string(),
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
        )
        .unwrap();

        let mut first = true;
        let mut last_str = false;
        for element in &self.family {
            if !first && !last_str {
                write!(f, "-").unwrap();
            }
            write!(
                f,
                "{}{}",
                element.0.map(|i| i.to_string()).unwrap_or_default(),
                element.1
            )
            .unwrap();
            last_str = !element.1.is_empty();
            first = false;
        }
        f
    }
}

impl FancyDisplay for ChainType {
    fn to_fancy_string(&self) -> String {
        match self {
            Self::Heavy => "H",
            Self::LightKappa => "Κ",
            Self::LightLambda => "Λ",
            Self::Iota => "Ι",
        }
        .to_string()
    }
}

impl FancyDisplay for GeneType {
    fn to_fancy_string(&self) -> String {
        match self {
            Self::V => "V",
            Self::J => "J",
            Self::C(None) => "C",
            Self::C(Some(Constant::A)) => "Α",
            Self::C(Some(Constant::D)) => "Δ",
            Self::C(Some(Constant::E)) => "Ε",
            Self::C(Some(Constant::G)) => "Ɣ",
            Self::C(Some(Constant::M)) => "Μ",
            Self::C(Some(Constant::O)) => "Ο",
            Self::C(Some(Constant::T)) => "Τ",
        }
        .to_string()
    }
}
