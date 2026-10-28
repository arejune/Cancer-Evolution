"""Text with numbers taken from the ledger, for the tables. Tokens:
  {L:id}  the rendered text of ledger entry id
  {M:x}   a label that contains digits but is not a number (for example A0)
  $...$   inline LaTeX math
Plain text may not contain digits.
"""
import re, json
from threshold_transfer.paths import RESULTS
LEDGER = json.load(open(RESULTS / "ledger.json"))["entries"]
TOK = re.compile(r"\{L:([^}]+)\}|\{M:([^}]+)\}|\$([^$]+)\$")
UNI = {"τ": r"$\tau$", "Δ": r"$\Delta$", "ρ": r"$\rho$", "φ": r"$\phi$", "κ": r"$\kappa$", "θ": r"$\theta$",
       "≤": r"$\leq$", "≥": r"$\geq$", "×": r"$\times$", "−": r"$-$", "→": r"$\rightarrow$", "–": "--", "·": r"$\cdot$",
       "%": r"\%", "&": r"\&", "#": r"\#", "_": r"\_"}


def parts(s):
    """Yield (kind, value): kind in text / num / label / math."""
    pos = 0
    for m in TOK.finditer(s):
        if m.start() > pos: yield "text", s[pos:m.start()]
        if m.group(1): yield "num", m.group(1)
        elif m.group(2): yield "label", m.group(2)
        else: yield "math", m.group(3)
        pos = m.end()
    if pos < len(s): yield "text", s[pos:]


def num_text(i):
    """Rendered text of ledger entry i, with a true minus sign before digits."""
    t = LEDGER[i]["text"]
    return re.sub(r"(?<![\w.])-(?=\d)", "−", t)          # true minus sign before a digit


def tex_text(s):
    """Replace Unicode symbols and LaTeX special characters by their LaTeX forms."""
    return "".join(UNI.get(c, c) for c in s)


def to_tex(s):
    """Render a rich string as LaTeX."""
    out = []
    for k, v in parts(s):
        if k == "text": out.append(tex_text(v))
        elif k == "num": out.append(tex_text(num_text(v)))
        elif k == "label": out.append(tex_text(v))
        else: out.append(f"${v}$")
    return "".join(out)


def ledger_ids(s):
    """Ledger ids used in a rich string."""
    return [v for k, v in parts(s) if k == "num"]


def check_no_digits(s, where=""):
    """Raise ValueError if the plain text of a rich string contains a digit."""
    for k, v in parts(s):
        if k == "text" and re.search(r"\d", v):
            raise ValueError(f"typed digit in plain text {where!r}: {v!r}")
