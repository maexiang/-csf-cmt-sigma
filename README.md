# csf-cmt-sigma

Affine rescaling analysis of lattice QCD symmetry energy σ(T, Φ).

## Overview

This repository contains the analysis code for testing whether the
QCD symmetry energy σ(T) at finite gravitational potential Φ is an
affine rescaling of the zero-potential lattice result:

    σ(T, Φ) = σ_lat( T / (1 + β Φ / c²) )

If the scaling holds, the shape of σ(T) is environment-independent
and only the critical temperature shifts.

## Repository structure

    data/       lattice QCD input (T/T_c, σ/σ0)
    src/        analysis code
    figures/    output figures

## Reproduce

    python src/rescaling.py

## License

MIT License. See LICENSE.
