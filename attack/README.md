# Attacks

- Implementation of our attacks. Sage and pari are needed.
- To compile the pari implementation in pari_tools run `make`.
- To run tests for pari run `make test` and then run `test`.
- Experiments are in `experiments.py` with results in `results` and graph visualisation in `results_visual`.
- ZVP-GLV attack on Shamir's trick with GLV-SAC is in `zvp_glv_sac`
- ZVP-GLV attack on the interleaving algorithm (alternative version) in `zvp_glv_inter_easy_prec`. It uses precomputed DCP points in `results/interleaving_secp256k1_remapped_*`. 
- The "hard" ZVP-GLV attack on the interleaving algorithm is done only on the first window and so the results can be reconstructed from the DCP points from above.
- The classical ZVP attack on the signed LTR is in `zvp_ltr_signed`
- The `dcp_analysis` jupyter notebook contains a few short scripts for the analysis of dcp.