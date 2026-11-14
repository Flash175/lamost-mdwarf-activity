# -*- coding: utf-8 -*-
"""run_all.py — 复现 §3.4（冷端 dM5-dM8：活动率、Wilson、探测比例、S/N 匹配、区间重叠）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_build.py", "step2_coolend.py", "step3_matched.py",
          "step4_overlap.py", "step5_check_34.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step5"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)
