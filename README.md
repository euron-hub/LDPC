# Using C++ to write LDPC Encoder/Decoder for the IEEE 802.16e Rate 1/2 QC-LDPC H-Matrix

This project simulates a process of "Data → Encode → Modulation (BPSK) → Channel (AWGN) → Find LLR → Decode → Data".

## Features
- **Standard**: IEEE 802.16e (WiMAX) **Rate-1/2 QC-LDPC**, base matrix 12×24, lifting size **Z=96**  
- **Encoder**: Systematic encoding from H 
- **Decoder**: Sum-Product (SPA) with **φ-function**  
- **Iteration**: Iteration is default to be 20 times

#### Performance
<img width="1387" height="567" alt="performance" src="https://github.com/user-attachments/assets/37384471-05ed-4a54-895b-8227ea993f4e" />


