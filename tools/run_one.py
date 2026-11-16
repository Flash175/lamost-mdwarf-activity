# -*- coding: utf-8 -*-
"""run_one.py — 由 Snakemake 以 `script:` 方式调用，执行单个复现脚本。

约定：
- 子进程用 config 里指定的解释器（默认 py311），不使用 PATH 里的 python。
- stdout/stderr 全部落到 logs/raw/<tag>.log。
- 无论成败都写 <output>（JSON 元信息：rc / status / 耗时 / 日志路径）。
- rc != 0 时以非零码退出，让 Snakemake 立即停在该规则上并报告是哪个子流程失败。

也可以脱离 Snakemake 手工调用（便于排查）：
    python tools/run_one.py --script 03_复现代码/复现流程/step1_sample_chain.py \\
        --cwd 03_复现代码/复现流程 --tag step1 --done logs/step1.done
"""
import argparse
import datetime
import json
import os
import pathlib
import subprocess
import sys
import time

DEFAULT_PY = r"C:\Users\lingz\anaconda3\envs\py311\python.exe"


def run(py: str, script: str, cwd: str, tag: str, out_path: str, extra=None) -> int:
    out = pathlib.Path(out_path)
    log = pathlib.Path("logs") / "raw" / f"{tag}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not os.path.isabs(script):
        script = os.path.abspath(script)
    if not os.path.exists(script):
        print(f"[{tag}] 脚本不存在: {script}", file=sys.stderr)
        return 2

    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"

    t0 = time.time()
    argv = [py, script] + [str(x) for x in (extra or [])]
    with open(log, "w", encoding="utf-8") as fh:
        fh.write(f"# start {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"# cwd    {os.path.abspath(cwd)}\n")
        fh.write(f"# python {py}\n")
        fh.write(f"# argv   {' '.join(argv)}\n\n")
        fh.flush()
        rc = subprocess.call(argv, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)

    elapsed = round(time.time() - t0, 1)
    status = "PASS" if rc == 0 else "FAIL"
    out.write_text(
        json.dumps(
            {
                "tag": tag,
                "script": script,
                "workdir": os.path.abspath(cwd),
                "python": py,
                "rc": rc,
                "status": status,
                "elapsed_s": elapsed,
                "log": str(log),
                "finished": datetime.datetime.now().isoformat(timespec="seconds"),
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"[{tag}] {status} rc={rc} {elapsed}s -> {log}")
    return rc


def main():
    # 1) Snakemake 注入模式
    try:
        sm = snakemake  # noqa: F821  (由 Snakemake 注入)
        rc = run(
            str(sm.params.py),
            str(sm.params.script),
            str(sm.params.cwd),
            str(sm.params.tag),
            str(sm.output[0]),
            extra=list(getattr(sm.params, "extra", []) or []),
        )
    except NameError:
        # 2) 手工调用模式
        ap = argparse.ArgumentParser(description="运行单个复现脚本并记录日志")
        ap.add_argument("--script", required=True)
        ap.add_argument("--cwd", required=True)
        ap.add_argument("--tag", required=True)
        ap.add_argument("--done", default=None,
                        help="完成标记 JSON 路径（被调脚本自己的 --out 会原样透传）")
        ap.add_argument("--py", default=DEFAULT_PY)
        a, unknown = ap.parse_known_args()
        rc = run(a.py, a.script, a.cwd, a.tag, a.done or f"logs/{a.tag}.done", extra=unknown)

    sys.exit(rc)


if __name__ == "__main__":
    main()
