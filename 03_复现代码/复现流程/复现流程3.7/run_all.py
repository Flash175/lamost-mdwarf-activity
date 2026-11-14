# -*- coding: utf-8 -*-
"""run_all.py — 复现 §3.7（清理后的 Teff 轴：K 5,047 谱 + M 642,596 星）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_teffaxis.py", "step2_check_37.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step2"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)
