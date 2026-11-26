# LGE test data (not shipped)

LGE material data and deck-run configs used only by `tests/test_gui_case_model.cpp`
(GUI acceptance tests 1 and 2). They stay out of the release folder and out of
`fluids/`.

| File | Source |
|---|---|
| `fw68d_pair/solution.csv` | copy of `data/R410A_FW68D/fw68d.csv` |
| `fw68d_pair/pair.txt` | header written for the GUI fluid format |
| `os_load_130hz.csv` | `D:\Seongho\Sim\pancake\lge-compressor-transient-20261119-02a4b6c\inputs\os_load_130hz.csv` |
| `configs/o4_tester_fbns_thermal.txt` | `lge-thermal-v2-20261029\oil\O4\O4\0003_equilibrium_thermal3_806c7579` |
| `configs/x4_tester_coupled_thermal.txt` | `lge-thermal-v3-c274-20261114\mix\gas_on\X4\runs\X4\0003_equilibrium_thermal3_474e01a5` |
| `configs/optical_rig_coupled.txt` | `silicone-coupled-20261114-b271002-local\runs\A-G_002_pinned_720x192` |
| `configs/compressor_gcvc_steady.txt` | `lge-compressor-transient-20261119-02a4b6c\b0\mean\009_newton6` |
| `configs/compressor_vc_steady.txt` | `lge-compressor-vc-steady-20261123-00e5d1a\runs\V-mean\007_newton4` |
| `configs/compressor_transient.txt` | `lge-compressor-transient-20261119-02a4b6c\b2\rev03` |
| `configs/p5g_tester_fbns_gas.txt` | `lge-fbns-tester-20261117\runs\P5g\runs\0007_newton4_0755bbea` |

All config paths are relative to `D:\Seongho\Sim\pancake`. The configs are byte
copies (CRLF, afml-r6 paths kept), because the round-trip test checks that
saving keeps those bytes.
