
# LDPC Encoder/Decoder for IEEE 802.16e (Rate 1/2 QC-LDPC)

This project implements an **LDPC (Low-Density Parity-Check) Encoder and Decoder** in **C++** based on the **IEEE 802.16e WiMAX standard** for the **Rate-1/2 Quasi-Cyclic LDPC (QC-LDPC) code**.  
It simulates the full digital communication chain:

```
Data → Encode → Modulation (BPSK) → Channel (AWGN) → Find LLR → Decode → Data
````



## Features

- **Standard**: IEEE 802.16e (WiMAX) Rate-1/2 QC-LDPC  
  - Base matrix: `12 × 24`  
  - Lifting size: `Z = 96`  
  - Codeword length: `N = 2304`  
  - Message length: `K = 1152`  
- **Encoder**: Systematic encoding based on the parity-check matrix `H`.
- **Decoder**:  
  - Implements **Sum-Product Algorithm (SPA)** in log-domain using the **φ-function**.  
  - Supports configurable iteration count (default = **20**).
- **Simulation**:  
  - BPSK modulation  
  - AWGN channel with adjustable SNR  
  - BER (Bit Error Rate) performance evaluation across SNR values.


## Project Structure

- **LDPC Matrix Construction**  
  - Builds the full `H` matrix from the IEEE 802.16e base matrix.
- **Encoding**  
  - Computes parity bits using Gaussian elimination and systematic encoding.
- **Channel Model**  
  - Adds **AWGN noise** using Box-Muller method for Gaussian samples.
- **Decoding**  
  - Implements iterative **SPA decoding** in the log domain with the φ-function.
- **Simulation Loop**  
  - Runs Monte Carlo simulations across SNR range `0.0 → 2.0 dB` in `0.2 dB` steps.  
  - Stops when `5000` codewords are processed **or** `100` bit errors are observed.  
  - Prints results in CSV format:
    ```
    SNR, sets, total_bits, bit_errors, BER
    ```


## Build & Run

### Requirements
- C++11 or newer
- Standard libraries only (no external dependencies)

### Compile
```bash
g++ -O2 -std=c++11 ldpc.cpp -o ldpc
````

### Run

```bash
./ldpc
```

The output will be printed to `stdout` in CSV format, e.g.:

```
SNR,sets,total_bits,bit_errors,BER
0,5000,5760000,xxxx,xxxx
0.2,5000,5760000,xxxx,xxxx
...
2.0,5000,5760000,xxxx,xxxx
```



## Example Result (BER vs. SNR)

<img width="1387" height="567" alt="performance" src="https://github.com/user-attachments/assets/37384471-05ed-4a54-895b-8227ea993f4e" />

---

## References

* IEEE 802.16e WiMAX Standard – LDPC coding specifications.
* Gallager, R. G. *Low-Density Parity-Check Codes*, 1962.
* Richardson, T. & Urbanke, R. *Modern Coding Theory*, 2008.

