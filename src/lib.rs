//! 库接口（供集成测试复用）。

pub mod chem;
pub mod engine;
pub mod model;
pub mod parser;
pub mod storage;
pub mod web;

pub const FIXTURE: &str = include_str!("../fixtures/sample.rxn");
pub const FIXTURE_NAME: &str = "sample.rxn";
