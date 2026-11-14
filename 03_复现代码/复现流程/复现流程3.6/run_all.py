# -*- coding: utf-8 -*-
"""run_all.py — 复现 §3.6（热端 K 型：逐型比例、χ²、Teff 轴、S/N 匹配、M 侧对照）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_kpilot.py", "step2_mside.py", "step4_check_36.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step4"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)