# Lightning

Paper: Lightning, Field-Agnostic Super-Efficient
Polynomial Commitment Scheme (`https://eprint.iacr.org/2026/258.pdf`).


Prototype C++ code that implements the **Lightning PCS**, it uses **Brakedown** as the inner code.



## Dependencies

- [mcl](https://github.com/herumi/mcl) 
- GMP
- OpenSSL (for `SHA256`)

## Build

From the project root:

```bash
g++ -O2 -std=c++17 -o lightning brakedown.cpp lightning.cpp -lmcl -lgmp -lcrypto
```

## Run

```bash
./lightning 24 0.004
```


