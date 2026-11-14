# -*- coding: utf-8 -*-
"""run_all.py — 复现 §3.2（探测比例、S/N 匹配、阈值不完备、单历元、变幅示例、曝光采样）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_build.py", "step2_detection.py", "step3_snmatch.py",
          "step4_threshold.py", "step5_epoch.py", "step6_variability.py",
          "step7_exposure.py", "step8_check_32.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time()-t0, rc))
    if s.startswith("step8"): last = rc
    elif rc != 0: sys.exit(rc)
sys.exit(last)
