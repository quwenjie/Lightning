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
g++ -O2 -std=c++17 -o lightning brakedown.cpp lightning_brakedown.cpp -lmcl -lgmp -lcrypto
```

## Run

```bash
./lightning
```

## Parameters

Lightning-related parameters currently live as constants near the top of `lightning_brakedown.cpp`:

- `poly_var`: number of variables (default `24`)
- `delta`: the target Lightning code distance
