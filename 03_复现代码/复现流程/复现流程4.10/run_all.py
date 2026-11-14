
# -*- coding: utf-8 -*-
"""run_all.py — 复现 §4.10（金属度通道）：量、正文比对、Table 5 逐格、Figure 9 重画。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
rcs = {}
steps = ("step1_build.py", "step2_check_410.py", "step3_table5.py", "step4_make_figure9.py")
for s in steps:
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time() - t0, rc))
    rcs[s] = rc
    if s == "step1_build.py" and rc != 0:
        sys.exit(rc)
sys.exit(any(rcs.get(s, 1) for s in steps[1:]) and 1 or 0)
