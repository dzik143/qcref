# QCRef
  - Reference implementation of **molecular integrals** over **Gaussian-Type Orbitals** (GTO).

# Overview:
  - **Naive and slow** (not optimized), but **clean and easy-to-understand** implementation.
  - Node.js and Python code.
  - May be a starting point for learning and/or research related to molecular integrals.

# Which integrals are implemented?
  - `<a|b>`   = overlap integral (S),
  - `<a|T|b>` = kinetic energy integral (KEI),
  - `<a|V|b>` = nuclear attraction integral (NAI),
  - `<ab|1/r|cd>` = electron repulsion integral (ERI).

# How does it work?
  - Integrals are calculated using **Obara-Saika algorithm**.
  - Original paper:
    Ref. Obara, S.; Saika, A. J Chem Phys 1986, 84, 3963.
    https://doi.org/10.1063/1.450106

# Example: One electron integrals (python):
  ```python
    from ObaraSaika_1E  import ObaraSaika_Overlap, ObaraSaika_Nuclear, ObaraSaika_Kinetic

    # XYZ coords of A,B functions centers.
    ra = [0.1, 0.5, 1.2]
    rb = [0.2, 0.7, 1.5]

    # XYZ coords of nuclei center (nuclear attraction integral only).
    rc = [0.3, 0.9, 1.8]

    # Exponential parameters of A,B functions (so-called zeta or alpha).
    za = 0.6
    zb = 0.8

    # Angular momentum for A and B functions.
    #     a func  b func
    #     x y z   x y z
    q = [ 0,0,0,  1,0,2 ]
    #     s       fxzz

    # Calculate one-electron integrals.
    s_ab   = ObaraSaika_Overlap(za, zb, ra, rb, q)
    nai_ab = ObaraSaika_Nuclear(za, zb, ra, rb, rc, q)
    kei_ab = ObaraSaika_Kinetic(za, zb, ra, rb, q)
  ```

# Example: Two electrons integral (python):
  ```python
    from ObaraSaika_ERI import ObaraSaika_ERI

    # XYZ coords of A,B,C,D functions centers.
    ra = [0.1, 0.5, 1.2]
    rb = [0.2, 0.7, 1.5]
    rc = [0.3, 0.9, 1.8]
    rd = [0.4, 1.1, 2.1]

    # Exponential parameters of A,B,C,D functions (so-called zeta or alpha).
    za = 0.6
    zb = 0.8
    zc = 1.0
    zd = 1.3

    # Angular momentum for A,B,C,D functions.
    #     a func  b func  c func  d func
    #     x y z   x y z   x y z   x y z
    q = [ 0,0,0,  1,0,1,  0,1,2,  1,0,0 ]
    #     s       dxz     fxzz    px

    # Calculate one-electron integrals.
    g_ab = ObaraSaika_ERI(za, zb, zc, zd, ra, rb, rc, rd, q)
  ```
