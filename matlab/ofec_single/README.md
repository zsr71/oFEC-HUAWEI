# MATLAB oFEC Single Pipeline

`run_ofec_single.m` is a MATLAB translation of `apps/ofec_single.cpp`. It wires
oFEC bit generation, encoding, QPSK modulation, AWGN, log-likelihood recovery,
and Chase-256 decoding to reproduce the same pre-/post-FEC BER reporting as the
C++ utility.

## Usage

```matlab
result = run_ofec_single(); % uses the default debug_L6 scenario
```

Optional name/value pairs mirror the configurable constants in the C++ entry
point, for example:

```matlab
result = run_ofec_single('Label','debug_L4', ...
                         'EbN0dB', 3.0, ...
                         'ChaseL', 4, ...
                         'AlphaFill', 0.15, ...
                         'BetaFill', 0.75);
```

The returned struct contains the scenario label, the Eb/N0 setting, BER
statistics before and after decoding, and per-tile early-stop percentages.

> **Note:** The current implementation focuses on the float (LLR\_BITS=16) path,
> mirroring the configuration used by the single-run workflow.
