# Attacks

Implementation of our attack. Sage and pari are needed.
To compile the pari implementation in pari_tools run `make`.
To run tests for pari run `make test` and then run `test`.
Experiments are in `experiments.py` with results in `results`.


## ZVP-GLV

Implementation of ZVP-GLV on SM/MSM in `zvp_glv_sm.py` and `zvp_glv_msm.py`

- Example of input/output for SM:

    Input: `{"curve": [0, 3], "field": 290551549417814750367658356818917722163, "lam": 113295953501134880362935974281619428327, "beta": 10580401330595003717415621069543210058, "order": 290551549417814750333932559803390773607, "k": 44532035905438573032946033634342206470, "k0": 10927427345391088509, "k1": 10219904264712657963, "zvp_target": 4}`

    Output:
    `{"params": , "k0_upper": 4, "zvp_recovered_0": 3, "nguesses_0": [1, 1, 2, 2], "k1_upper": 4, "zvp_recovered_1": 3, "nguesses_1": [1, 1, 2, 2], "time_zvp": 1.2122609615325928}`

- see `results/zvp_glv_*` for more examples

## 2-step

Implementation of the 2-step attack in `twostep.py`

- Example of input/output for 2-step on SM:

    Input: `{"params": {"curve": [0, 3], "field": 4181204319369823, "lam": 4064216550687834, "beta": 2676826888224834, "order": 4181204318045707, "k": 2585698485112112, "k0": 36224342, "k1": 40038981, "zvp_target": 6}`
    
    Output: `{"k0_upper": 34, "zvp_recovered_0": 6, "nguesses_0": [1, 1, 2, 3, 5, 1], "k1_upper": 38, "zvp_recovered_1": 6, "nguesses_1": [1, 1, 2, 3, 5, 1], "time_zvp": 5.867018222808838, "time_bsgs": 3.5167787075042725}`

- see `results/twostep_*` for more examples



## Notes
- All DCP instances are hardcoded for the polynomial 
`x1+x2+2` only. If needed, `pari_tools.dcp.c`,`zvp_glv_sm`, `msm.py` can be changed to other polynomial.
- For convenience we assume that the decomposition of the private scalar results in positive scalars.
