# -*- coding: utf-8 -*-
"""run_all.py — 复现 §2.6（独立 Hβ 与 L_Hα/L_bol 判据）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_hbeta_tracks.py", "step2_efficiency.py", "step3_lbol.py", "step4_check_26.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step4"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)
