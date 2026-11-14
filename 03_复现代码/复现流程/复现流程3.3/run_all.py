# -*- coding: utf-8 -*-
"""run_all.py — 复现 §3.3（6D 运动学 + 附录 Table 5 + LAMOST-RV 分片与交叉核对）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_build.py", "step2_bins.py", "step3_errors.py", "step4_lamost.py",
          "step5_crosscheck.py", "step6_table5.py", "step7_check_33.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step7"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)
