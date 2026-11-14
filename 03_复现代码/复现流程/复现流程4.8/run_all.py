
# -*- coding: utf-8 -*-
"""run_all.py — 复现 §4.8（热端试点）。"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); last = 0
for s in ("step1_build.py", "step2_check_48.py"):
    print("\n=== %s ===" % s); t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(HERE, s)])
    print("--- %s finished in %.0fs (rc=%d)" % (s, time.time() - t0, rc))
    if s.startswith("step2"):
        last = rc
    elif rc != 0:
        sys.exit(rc)
sys.exit(last)
