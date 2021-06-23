<h1 align="center">ark-bpr20 (CRH-based)</h1>

<p align="center">
    <img src="https://github.com/arkworks-rs/groth16/workflows/CI/badge.svg?branch=master">
    <a href="https://github.com/arkworks-rs/groth16/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/arkworks-rs/groth16/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="https://deps.rs/repo/github/arkworks-rs/groth16"><img src="https://deps.rs/repo/github/arkworks-rs/groth16/status.svg"></a>
</p>

The arkworks ecosystem consist of Rust libraries for designing and working with __zero knowledge succinct non-interactive arguments (zkSNARKs)__. This repository contains an efficient implementation of the simulation extractable zk-SNARK published in CANS20 paper [[BPR20a]](https://link.springer.com/chapter/10.1007/978-3-030-65411-5_22) which is also avalable in Section 4 of [[BPR20b]](https://eprint.iacr.org/2020/1306), and is done by Oussama Amine (University of Oslo) and Karim Baghery (KU Leuven).

This library is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone -b CRH-based https://github.com/arkworks-rs/bpr20.git CRH-based
cd CRH-based
cargo build --release
```

This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test
```
and for benchmarking the scheme with `RAYON_NUM_THREADS=4` threads, run the following command,  
```bash
RAYON_NUM_THREADS=4 cargo bench --no-default-features --features "std parallel" -- --nocapture
```

## Empirical performance

Below is the empricial performance of several weak and strong simulation extractable zk-SNARKs in `Arkworks`. Note that Groth's zk-SNARK is proven to achieve weak simulation extractability [[BKSV20]](https://eprint.iacr.org/2020/811).  
We benchmark the zk-SNARKs on an R1CS instance for different curves and report proving and verifying times for each constrain with 100 iterations for the prover and 100.000 iterations for the verification. 
The benchmarks were obtained using a 2.50 GHz Intel Core i5-7200U CPU, in multi-threaded mode with N=4, with 16 GB RAM, using the Bls12_381, MNT4_298, MNT6_298, MNT4_753 and MNT6_753 curves.

Abbreviations used: <i>SE</i> = simulation extractable, <i>ns</i> = nanoseconds, <i>RO</i> = Random Oracle, <i>CRH</i> = Collision Resistant Hash.

| Curve | zk-SNARK | Secuiry | Per-constraint proving time, ns | Per-constraint verifying time, ns |
| :--- | :---: | :---: | :---: | :---: |
| BLS12-381 | [Groth16](arkworks-rs/groth16)               | Weak SE   | 28513 | ??? |
| BLS12-381 | [GM17](arkworks-rs/gm17)                     | Strong SE | ???   | ??? |
| BLS12-381 | [BPR20-RO](arkworks-rs/bpr20/RO-based)       | Strong SE | ???   | ??? |
| BLS12-381 | [BPR20-CRH](arkworks-rs/bpr20/CRH-based)     | Strong SE | ???   | ??? |
| **MNT4-298** | [Groth16](arkworks-rs/groth16)            | Weak SE   | 29462 | ??? |
| **MNT4-298** | [GM17](arkworks-rs/gm17)                  | Strong SE | ???   | ??? |
| **MNT4-298** | [BPR20-RO](arkworks-rs/bpr20/RO-based)    | Strong SE | ???   | ??? |
| **MNT4-298** | [BPR20-CRH](arkworks-rs/bpr20/CRH-based)  | Strong SE | ???   | ??? |
| MTN6-298 | [Groth16](arkworks-rs/groth16)                | Weak SE   | 26478 | ??? |
| MTN6-298 | [GM17](arkworks-rs/gm17)                      | Strong SE | ???   | ??? |
| MTN6-298 | [BPR20-RO](arkworks-rs/bpr20/RO-based)        | Strong SE | ???   | ??? |
| MTN6-298 | [BPR20-CRH](arkworks-rs/bpr20/CRH-based)      | Strong SE | ???   | ??? |
| **MNT4-753** | [Groth16](arkworks-rs/groth16)            | Weak SE   | 255144 | ??? |
| **MNT4-753** | [GM17](arkworks-rs/gm17)                  | Strong SE | ???   | ??? |
| **MNT4-753** | [BPR20-RO](arkworks-rs/bpr20/RO-based)    | Strong SE | ???   | ??? |
| **MNT4-753** | [BPR20-CRH](arkworks-rs/bpr20/CRH-based)  | Strong SE | ???   | ??? |
| MTN6-753 | [Groth16](arkworks-rs/groth16)                | Weak SE   | 198195 | ??? |
| MTN6-753 | [GM17](arkworks-rs/gm17)                      | Strong SE | ???   | ??? |
| MTN6-753 | [BPR20-RO](arkworks-rs/bpr20/RO-based)        | Strong SE | ???   | ??? |
| MTN6-753 | [BPR20-CRH](arkworks-rs/bpr20/CRH-based)      | Strong SE | ???   | ??? |



## License

This library is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in this library by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

## Acknowledgements

This work was supported by:
the Defense Advanced Research Projects Agency (DARPA) under Contract No. HR001120C0085; 
a Google Faculty Award;
the National Science Foundation;
the UC Berkeley Center for Long-Term Cybersecurity;
and donations from the Ethereum Foundation, the Interchain Foundation, and Qtum.

An earlier version of this library was developed as part of the paper *"[ZEXE: Enabling Decentralized Private Computation][zexe]"*.

[zexe]: https://ia.cr/2018/962

