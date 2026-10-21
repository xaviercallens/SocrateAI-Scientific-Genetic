#!/usr/bin/env python3
"""Build the textbook, or any single section of it, in isolation.

The book is deliberately split one-file-per-section so that a section can be
drafted, compiled and reviewed without the other nine in play. Building a single
section reuses `main.tex`'s preamble verbatim -- a section that compiles alone
therefore also compiles in the book, which is the property that makes
section-at-a-time review trustworthy.

    python build.py                 # whole book
    python build.py --section 3     # just section 3, into build/sec03.pdf
    python build.py --list          # section inventory and build status

Exit status is non-zero if LaTeX reported an error, an undefined reference, an
undefined citation or an overfull box, so this doubles as the quality gate.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SECTIONS = ROOT / "sections"
BUILD = ROOT / "build"
MAIN = ROOT / "main.tex"


def section_files() -> list[Path]:
    """Section sources in reading order (sec01..sec10), appendices excluded."""
    return sorted(SECTIONS.glob("sec[0-9][0-9]_*.tex"))


def preamble() -> str:
    """Everything in main.tex before \\begin{document}."""
    text = MAIN.read_text(encoding="utf-8")
    cut = text.index(r"\begin{document}")
    return text[:cut]


def solo_wrapper(target: Path) -> str:
    """A standalone document containing exactly one section, same preamble."""
    body = "\n".join([
        r"\begin{document}",
        r"\setcounter{section}{%d}" % (section_number(target) - 1),
        r"\input{%s}" % target.with_suffix("").as_posix(),
        r"\end{document}",
        "",
    ])
    return preamble() + body


def section_number(path: Path) -> int:
    m = re.match(r"sec(\d\d)_", path.name)
    if m is None:
        raise ValueError(f"not a numbered section file: {path.name}")
    return int(m.group(1))


def run_latex(tex: Path, workdir: Path) -> tuple[int, str]:
    """Run latexmk to completion and return (exit code, combined log)."""
    proc = subprocess.run(
        ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
         f"-outdir={workdir}", str(tex)],
        # latexmk echoes bytes from fonts and logs that are not valid UTF-8, so
        # decode permissively rather than letting the build crash on them.
        cwd=str(ROOT), capture_output=True, text=True, errors="replace", timeout=1800,
    )
    # Grade the .log only. latexmk reruns LaTeX until references settle, and its
    # stdout carries every pass; counting those would report first-pass "undefined
    # reference" warnings that the final pass has already resolved. The .log file
    # is overwritten each pass, so it is exactly the last one.
    log = workdir / (tex.stem + ".log")
    text = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    return proc.returncode, text if text else (proc.stdout + proc.stderr)


def grade(log: str) -> dict[str, int]:
    """The four defects we refuse to ship, counted from the LaTeX log."""
    return {
        "errors": len(re.findall(r"^! ", log, re.M)),
        "undefined_refs": len(re.findall(r"Reference .* undefined", log)),
        "undefined_cites": len(re.findall(r"Citation .* undefined", log)),
        "overfull": len(re.findall(r"Overfull \\hbox", log)),
    }


def report(label: str, code: int, log: str, *, strict: bool = True) -> bool:
    """Print the grade and say whether it passes.

    A solo section build cannot resolve `\\ref`s into other sections, so undefined
    references are expected there and are reported without failing. The full-book
    build is the authoritative gate and fails on all four defects.
    """
    g = grade(log)
    fatal = dict(g) if strict else {k: v for k, v in g.items()
                                    if k not in {"undefined_refs", "undefined_cites"}}
    ok = code == 0 and not any(fatal.values())
    mark = "PASS" if ok else "FAIL"
    print(f"[{mark}] {label}: exit={code} " + " ".join(f"{k}={v}" for k, v in g.items()))
    if not ok:
        for line in re.findall(r"^! .*$", log, re.M)[:5]:
            print(f"        {line}")
        for line in re.findall(r"^.*(?:undefined|Overfull \\hbox).*$", log, re.M)[:5]:
            print(f"        {line.strip()[:150]}")
    return ok


def build_section(n: int) -> bool:
    matches = [p for p in section_files() if section_number(p) == n]
    if not matches:
        print(f"no section file for {n}; have "
              f"{[section_number(p) for p in section_files()]}")
        return False
    target = matches[0]
    BUILD.mkdir(exist_ok=True)
    wrapper = ROOT / f"_solo_sec{n:02d}.tex"
    wrapper.write_text(solo_wrapper(target), encoding="utf-8")
    try:
        code, log = run_latex(wrapper, BUILD)
        ok = report(f"section {n:02d} ({target.name})", code, log, strict=False)
    finally:
        wrapper.unlink(missing_ok=True)
    produced = BUILD / f"_solo_sec{n:02d}.pdf"
    if produced.exists():
        shutil.move(str(produced), str(BUILD / f"sec{n:02d}.pdf"))
    return ok


def build_book() -> bool:
    BUILD.mkdir(exist_ok=True)
    code, log = run_latex(MAIN, BUILD)
    return report("full book", code, log)


def inventory() -> None:
    files = section_files()
    print(f"{'sec':>4}  {'file':44}  {'lines':>6}  {'built':>5}")
    for p in files:
        n = section_number(p)
        built = "yes" if (BUILD / f"sec{n:02d}.pdf").exists() else "-"
        lines = len(p.read_text(encoding='utf-8').splitlines())
        print(f"{n:>4}  {p.name:44}  {lines:>6}  {built:>5}")
    missing = sorted({1, 2, 3, 4, 5, 6, 7, 8, 9, 10} - {section_number(p) for p in files})
    if missing:
        print(f"\nnot yet written: {missing}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--section", type=int, help="build only this section (1-10)")
    ap.add_argument("--all-sections", action="store_true",
                    help="build every section separately, then report")
    ap.add_argument("--list", action="store_true", help="show the inventory")
    args = ap.parse_args()

    if args.list:
        inventory()
        return 0
    if args.section is not None:
        return 0 if build_section(args.section) else 1
    if args.all_sections:
        results = {section_number(p): build_section(section_number(p))
                   for p in section_files()}
        print("\n" + " ".join(f"{n}:{'ok' if v else 'FAIL'}" for n, v in sorted(results.items())))
        return 0 if all(results.values()) else 1
    return 0 if build_book() else 1


if __name__ == "__main__":
    sys.exit(main())
