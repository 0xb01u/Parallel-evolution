# About
This is a sequential re-implementation of the evolution simulation using the Rust programming language. It tries to replicate the original, unoptimized source code's behavior and "intention" as close as possible. Thus, it doesn't present any special optimizations, by parallelization or otherwise. The purpose of this implementation was just to practice Rust programming.

## Building
For debugging (no compiler optimizations, debugging symbols, debugging-purposed special functions):
```
$ cargo build --features "debug"
```
The binary will be placed inside `./target/dev/`

For production (compiler optimizations):
```
$ cargo build --release
```
The binary will be placed inside `./target/release/`

You have to have installed rust and cargo on your system (they are usually bundled together).

