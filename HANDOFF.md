# Handoff: UCL AI Agents Seminar

## What this is

Materials for a 1-hour seminar at UCL Physics: "AI Coding Agents: A
Reductive Account for Physicists." Demystifies AI coding agents by
reducing them to composable primitives physicists already have intuition
for.

## Current state

**Slides are fully rebuilt** with a custom LaTeX Beamer system replacing
the original reportlab-based `build_slides.py`.

### Build system

- **Engine**: LuaLaTeX (TeX Live 2023)
- **Theme**: Custom `TJO` Beamer theme (`slides/beamer*TJO.sty`)
- **Fonts**: Whitney (text), MathTime Pro 2 (math), Iosevka Custom (code)
- **Colors**: Whitney Teal palette (dk2=#335B74, accent1=#1CADE4)
- **Build**: `make` at project root; `make s00` through `make s04` for
  individual sections

### Sections (5 decks, 40 slides total)

| File | Slides | Content |
|------|--------|---------|
| `section00.tex` | 10 | WHY: timeline, audience poll, personal journey, project showcase |
| `section01.tex` | 8 | LLM as stateless nondeterministic function |
| `section02.tex` | 7 | Chat illusion, HTTP reality, context window |
| `section03.tex` | 7 | Agent loop, tools, live demo |
| `section04.tex` | 8 | Audience vote on what to build next |

### Key files

```
slides/
  beamerthemeTJO.sty          # Master theme (loads 4 sub-themes)
  beamercolorthemeTJO.sty     # Whitney Teal color palette
  beamerfontthemeTJO.sty      # Whitney + MathTime Pro 2 + Iosevka
  beamerinnerthemeTJO.sty     # Title page, \sectiondivider, \statementslide, codeblock
  beamerouterthemeTJO.sty     # Frame title + cyan accent line, footline
  preamble.tex                # Shared packages, TikZ styles, math shortcuts
  section0[0-4].tex           # Slide content
Makefile                      # Build system
build_slides.py               # Old reportlab system (kept for reference)
```

### Demo code (unchanged from initial commit)

- `01-*/demo.py` — Temperature sweep showing nondeterminism
- `02-*/single_call.py` — Raw API call with JSON payload
- `02-*/chat_loop.py` — Terminal chat showing history growth
- `03-*/agent.py` — 124-line agent with read/write tools
- `03-*/tools.py` — Sandboxed file tools

### Known issues

1. **Section 00, slide 4** (survey): 4th bullet item clips at bottom.
   Fix: reduce `\fontsize` from 22 to 20.
2. **Section 00, slide 8** (projects): 9th table row (alethfeld) pushed
   off bottom. Fix: reduce table font or drop to 8 rows.
3. **Section 03, slide 3** (agent loop diagram): TikZ diagram slightly
   clipped at bottom. Fix: reduce `node distance` or scale tikzpicture.

### Style guide

A comprehensive style guide was extracted from 15+ years of TJO's
presentations and lives at `~/Projects/presentations/` (private repo
`tobiasosborne/presentations`). Reference it for any future slide work.

## What to do next

- Fix the three minor clipping issues noted above
- Consider adding screenshots/renders from actual projects (Lyr.jl
  volume renders, vectorfeld UI) to section 00 for visual impact
- Rehearse with timing — 40 slides for 60 minutes is comfortable
- The old `build_slides.py` can be removed once the new system is stable

## Technical notes

- `codeblock` uses `tcblisting` (verbatim) — frames containing it need
  `[fragile]`, and `\begin{codeblock}` / `\end{codeblock}` must start
  at column 0 (no leading whitespace)
- Font loading order in `beamerfontthemeTJO.sty` is critical:
  `luatex85` → `mtpro2[lite]` → `fontspec[no-math]`
- mtpro2 `[lite]` is used; `[complete]` may work with all 136 pfb files
