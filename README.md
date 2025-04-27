# DipoleFluctuations

This is a repository for an example project in the Oxford Dyson Program 2024-2025. I attempt to reestablish results achieved in Lau & Binney (2019), where multipole components of an Isocchrone N-body snapshots were analyzed to determine their contributions to the relaxation of the system. The N-Body simulation uses Falcon code created by Dehnen, and modified by Professor James Binney for parallel processing. The ReadSnapShot code attached delivers the files "PowerRadialVariation#.csv" and "DipoleVar#.csv" at time steps 0, 25, 50, ... 250. 

"PowerRadialVariation#.csv" delivers in the following format:
radius | |c1|^2 | |c2|^2 
0      | 0      | 0      
...

Where c1 and c2 are the dipole and quadrupole coefficients as defined in Lau & Binney (2019)

"DipoleVar#.csv"
radius | <cos(q)>
...

Where <cos(q)> Is the cosine angle between inner and outer dipoles, as specified in Lau & Binney (2019)


1. Dehnen, W. (2000) ‘A very fast and momentum-Conserving tree code’, The Astrophysical Journal, 536(1). doi:10.1086/312724.
2. Dehnen, W. (2002) ‘A hierarchical O(n) force calculation algorithm’, Journal of Computational Physics, 179(1), pp. 27–42. doi:10.1006/jcph.2002.7026. 
3. Lau, J.Y. and Binney, J. (2019) ‘Relaxation of spherical stellar systems’, Monthly Notices of the Royal Astronomical Society, 490(1), pp. 478–490. doi:10.1093/mnras/stz2567
