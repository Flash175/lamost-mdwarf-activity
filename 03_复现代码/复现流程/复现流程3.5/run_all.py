# -*- coding: utf-8 -*-
"""run_all.py — 复现 §3.5（周期样本、Ro 构造、面板/|z| 分层拟合、Figure 2）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_sample.py", "step2_fits.py", "step2b_tau_variant.py", "step3_figure.py", "step4_check_35.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step4"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)
