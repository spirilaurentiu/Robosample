"""Mechanical LaTeX balance check for $$...$$ blocks in a markdown file.

Not a pipeline - a single tool the paper-ingestor subagent runs after writing equations.md, so it never has to count braces by eye. Exit 1 if any issue.

    python3 references/tools/check_latex.py references/papers/<key>/equations.md
"""

import re
import sys

ENVS = ("align", "aligned", "cases", "bmatrix", "pmatrix", "matrix", "array", "split")


def issues(latex: str) -> list[str]:
    out = []
    if latex.count("{") != latex.count("}"):
        out.append(f"brace mismatch {{={latex.count('{')} }}={latex.count('}')}")
    if latex.count(r"\left") != latex.count(r"\right"):
        out.append(rf"\left={latex.count(r'\left')} \right={latex.count(r'\right')}")
    for env in ENVS:
        b = len(re.findall(rf"\\begin\{{{env}\*?\}}", latex))
        e = len(re.findall(rf"\\end\{{{env}\*?\}}", latex))
        if b != e:
            out.append(rf"\begin/\end {env}: {b}/{e}")
    return out


def main(path: str) -> int:
    text = open(path, encoding="utf-8", errors="replace").read()
    blocks = re.findall(r"\$\$(.*?)\$\$", text, flags=re.S)
    bad = 0
    for i, blk in enumerate(blocks, 1):
        # locate the nearest preceding <!-- eq:... --> anchor for a useful label
        anchor = "?"
        m = list(re.finditer(r"<!--\s*(eq:[^\s>]+)", text[: text.find(blk)]))
        if m:
            anchor = m[-1].group(1)
        errs = issues(blk)
        if errs:
            bad += 1
            print(f"[FAIL] block {i} ({anchor}): {'; '.join(errs)}")
    print(f"\n{len(blocks)} equations, {bad} with issues.")
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
