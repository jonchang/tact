"""Small helpers for writing the documentation figures as SVG, with one palette for light and dark mode."""

from __future__ import annotations

import re
from pathlib import Path
from xml.sax.saxutils import escape

FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif'

# Each class: (light declarations, dark declarations). Only the classes a figure uses go into its stylesheet.
STYLES: dict[str, tuple[str, str]] = {
    # backgrounds and text
    "bg": ("fill: #ffffff;", "fill: #1e2129;"),
    "fg": ("fill: #1f2328;", "fill: #e3e5e8;"),
    "muted": ("fill: #6b7280;", "fill: #9aa0a6;"),
    "ok-text": ("fill: #2e7d32;", "fill: #66bb6a;"),
    "bad-text": ("fill: #c62828;", "fill: #ef5350;"),
    "new-text": ("fill: #1565c0;", "fill: #64b5f6;"),
    "sole-text": ("fill: #a15c00;", "fill: #ffb74d;"),
    # branches
    "branch": ("fill: none; stroke: #1f2328; stroke-width: 1.6;", "stroke: #e3e5e8;"),
    "edge": ("fill: none; stroke: #6b7280; stroke-width: 2;", "stroke: #9aa0a6;"),
    "unsampled": ("fill: none; stroke: #9ca3af; stroke-width: 1.5; stroke-dasharray: 3 2;", "stroke: #6b7280;"),
    "placed": ("fill: none; stroke: #1565c0; stroke-width: 1.5; stroke-dasharray: 3 2;", "stroke: #64b5f6;"),
    "unlocked": ("fill: none; stroke: #2e7d32; stroke-width: 3.5;", "stroke: #66bb6a;"),
    "new-edge": ("fill: none; stroke: #1565c0; stroke-width: 2.5;", "stroke: #64b5f6;"),
    # clades
    "clade": ("fill: #e5e7eb; stroke: #9ca3af; stroke-width: 1.5;", "fill: #343a46; stroke: #6b7280;"),
    "new-clade": (
        "fill: #1565c0; fill-opacity: 0.10; stroke: #1565c0; stroke-width: 1.5; stroke-dasharray: 5 3;",
        "fill: #64b5f6; fill-opacity: 0.14; stroke: #64b5f6;",
    ),
    "sole": ("fill: #f5d9a8; stroke: #c27c0e; stroke-width: 1.5;", "fill: #4a3b22; stroke: #ffb74d;"),
    "intruder": ("fill: #f6c9c9; stroke: #c62828; stroke-width: 1.5;", "fill: #4a2628; stroke: #ef5350;"),
    # nodes
    "dot": ("fill: #6b7280;", "fill: #9aa0a6;"),
    "new-dot": ("fill: #1565c0;", "fill: #64b5f6;"),
    "tip": ("fill: #1f2328;", "fill: #e3e5e8;"),
    "crown": ("fill: #1f2328; stroke: #ffffff; stroke-width: 1.5;", "fill: #e3e5e8; stroke: #1e2129;"),
    "true-crown": ("fill: #ffffff; stroke: #6b7280; stroke-width: 2;", "fill: #1e2129; stroke: #9aa0a6;"),
    # annotations
    "band": ("fill: #2e7d32; fill-opacity: 0.10;", "fill: #66bb6a; fill-opacity: 0.12;"),
    "sole-box": (
        "fill: #c27c0e; fill-opacity: 0.06; stroke: #c27c0e; stroke-width: 1.5; stroke-dasharray: 6 4;",
        "fill: #ffb74d; stroke: #ffb74d;",
    ),
    "axis": ("fill: none; stroke: #6b7280; stroke-width: 1.5;", "stroke: #9aa0a6;"),
    "bracket": ("fill: none; stroke: #6b7280; stroke-width: 1.5;", "stroke: #9aa0a6;"),
    "ok-line": ("fill: none; stroke: #2e7d32; stroke-width: 1.5;", "stroke: #66bb6a;"),
    "bad-line": ("fill: none; stroke: #c62828; stroke-width: 1.5;", "stroke: #ef5350;"),
    "arrowhead": ("fill: #c62828;", "fill: #ef5350;"),
    "shackle": ("fill: none; stroke: #6b7280; stroke-width: 1.5;", "stroke: #9aa0a6;"),
    # modifiers
    "dashed": ("stroke-dasharray: 4 3;", ""),
    "small": ("font-size: 12px;", ""),
    "italic": ("font-style: italic;", ""),
    "bold": ("font-weight: bold;", ""),
}

PADLOCK = """<g id="padlock">
      <path class="shackle" d="M-3.5,-1 v-3 a3.5,3.5 0 0 1 7,0 v3"/>
      <rect class="muted" x="-5.5" y="-1" width="11" height="8" rx="1.5"/>
    </g>"""

ARROW = """<marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" \
orient="auto-start-reverse">
      <path class="arrowhead" d="M0,0 L10,5 L0,10 z"/>
    </marker>"""


def num(value: float) -> str:
    """Format a coordinate compactly and reproducibly."""
    text = f"{value:.1f}"
    return text[:-2] if text.endswith(".0") else text


def i(text: str) -> str:
    """Italicize part of a label, such as a genus name."""
    return f'<tspan class="italic">{escape(text)}</tspan>'


def text_width(markup: str, size: float = 12) -> float:
    """Roughly estimate the rendered width of a label, for laying out legends."""
    return len(re.sub(r"<[^>]+>", "", markup)) * size * 0.5


class Figure:
    """An SVG figure built up from elements, with a stylesheet holding only the classes it uses."""

    def __init__(self, width: float, height: float, title: str, desc: str, comment: str = "") -> None:
        """Start an empty figure."""
        self.width, self.height = width, height
        self.title, self.desc, self.comment = title, desc, comment
        self.defs: list[str] = []
        self.body: list[str] = []

    def add(self, element: str) -> None:
        """Add a raw SVG element."""
        self.body.append(element)

    def path(self, cls: str, d: str, extra: str = "") -> None:
        """Add a path."""
        self.add(f'<path class="{cls}" d="{d}"{extra}/>')

    def text(self, x: float, y: float, markup: str, cls: str = "fg", anchor: str | None = None) -> None:
        """Add a label; `markup` may contain tspans, such as from `i`."""
        anchor_attr = f' text-anchor="{anchor}"' if anchor else ""
        self.add(f'<text x="{num(x)}" y="{num(y)}"{anchor_attr} class="{cls}">{markup}</text>')

    def circle(self, cls: str, x: float, y: float, r: float) -> None:
        """Add a circle."""
        self.add(f'<circle class="{cls}" cx="{num(x)}" cy="{num(y)}" r="{num(r)}"/>')

    def rect(self, cls: str, x: float, y: float, width: float, height: float, rx: float = 0) -> None:
        """Add a rectangle."""
        rx_attr = f' rx="{num(rx)}"' if rx else ""
        self.add(
            f'<rect class="{cls}" x="{num(x)}" y="{num(y)}" width="{num(width)}" height="{num(height)}"{rx_attr}/>'
        )

    def triangle(self, cls: str, x: float, y: float, x_end: float, half: float) -> None:
        """Add a collapsed clade: a triangle from its crown at (x, y) to a base of height 2 * half at x_end."""
        d = f"M{num(x)},{num(y)} L{num(x_end)},{num(y - half)} L{num(x_end)},{num(y + half)} Z"
        self.path(cls, d)

    def padlock(self, x: float, y: float) -> None:
        """Add a padlock icon, marking a locked clade."""
        if PADLOCK not in self.defs:
            self.defs.append(PADLOCK)
        self.add(f'<use href="#padlock" x="{num(x)}" y="{num(y)}"/>')

    def arrow(self, cls: str, d: str) -> None:
        """Add a path ending in an arrowhead."""
        if ARROW not in self.defs:
            self.defs.append(ARROW)
        self.path(cls, d, ' marker-end="url(#arrow)"')

    def legend(self, x: float, y: float, items: list[tuple[str, str, str]]) -> None:
        """Add a row of legend items, each (kind, class, label) where kind is "line", "triangle" or "tip"."""
        for kind, cls, label in items:
            if kind == "triangle":
                self.path(cls, f"M{num(x)},{num(y - 5)} L{num(x + 22)},{num(y - 10)} L{num(x + 22)},{num(y)} Z")
            else:
                self.path(cls, f"M{num(x)},{num(y - 5)} H{num(x + 22)}")
                if kind == "tip":
                    self.circle("tip", x + 22, y - 5, 3.5)
            self.text(x + 30, y - 1, label, "small muted")
            x += 30 + text_width(label) + 28

    def stylesheet(self) -> str:
        """CSS for the classes this figure uses, with dark mode overrides."""
        markup = "\n".join(self.defs + self.body)
        used = {cls for attr in re.findall(r'class="([^"]+)"', markup) for cls in attr.split()} | {"bg"}
        names = [name for name in STYLES if name in used]
        light = "\n".join(f"    .{name} {{ {STYLES[name][0]} }}" for name in names)
        dark = "\n".join(f"      .{name} {{ {STYLES[name][1]} }}" for name in names if STYLES[name][1])
        return (
            f"  <style>\n{light}\n    @media (prefers-color-scheme: dark) {{\n{dark}\n    }}\n"
            f"    text {{ font-family: {FONT}; font-size: 13px; }}\n  </style>"
        )

    def write(self, path: Path) -> None:
        """Write the figure to `path`."""
        w, h = num(self.width), num(self.height)
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            'role="img" aria-labelledby="title desc">',
            f'  <title id="title">{escape(self.title)}</title>',
            f'  <desc id="desc">{escape(self.desc)}</desc>',
            "  <!-- Generated by docs/figures/make_figures.py; edit that instead of this file. -->",
        ]
        if self.comment:
            parts.append(f"  <!-- {self.comment} -->")
        parts.append(self.stylesheet())
        if self.defs:
            parts.append("  <defs>\n    " + "\n    ".join(self.defs) + "\n  </defs>")
        parts.append(f'  <rect class="bg" width="{w}" height="{h}" rx="8"/>')
        parts.extend("  " + element for element in self.body)
        parts.append("</svg>\n")
        path.write_text("\n".join(parts), encoding="utf-8")
