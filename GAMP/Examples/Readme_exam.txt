The directory of "2017244/results_reference" is the processing results for reference.
NOTE:    G    --> GPS
         GR   --> GPS + GLONASS
         GRCE --> GPS + GLONASS + BDS + Galileo
         kin  --> Kinematic mode
         SF   --> Single-frequency mode
         DF   --> Dual-frequency mode
         GIM  --> Ionosphere-constrained PPP
      no_GIM  --> Standard PPP
         wum  --> The adopted final precise products

e.g.,
./sh_ppp_1site 2017 244 1 G kin SF Y
The command line indicates GPS-only kinematic single-frequency ionosphere-constrained PPP.

./sh_ppp_1site 2016 183 1 GRC sta DF N
The command line indicates GPS + GLONASS + BDS static single-frequency standard PPP.