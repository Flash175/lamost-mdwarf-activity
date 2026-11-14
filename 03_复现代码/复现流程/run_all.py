# -*- coding: utf-8 -*-
"""run_all.py — 依次执行复现流程五步并汇报。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = ["step1_sample_chain.py","step2_bins_and_band.py","step3_decompositions.py",
         "step4_break_scan.py","step5_check_all.py"]
want = sys.argv[1:] or ["1","2","3","4","5"]; last = 0
for s in STEPS:
    if s[4] not in want and not s.startswith("step5"): continue
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step5"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)
