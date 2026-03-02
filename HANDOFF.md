# Handoff: AI Agents Seminar

## What this is

Materials for a 1-hour seminar: "AI Coding Agents: A
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

### Sections (5 decks, 45 slides total)

| File | Slides | Content |
|------|--------|---------|
| `section00.tex` | 12 | WHY: timeline, audience poll, personal journey, project showcase (2 pages) |
| `section01.tex` | 9 | LLM as stateless nondeterministic function, temperature formula + limits |
| `section02.tex` | 8 | Chat illusion, HTTP reality, no-memory statement, context window |
| `section03.tex` | 8 | Agent loop, tools, "LLM never executes" statement, live demo |
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

All three previously reported clipping issues have been fixed:

1. ~~Section 00, slide 4 (survey)~~: font reduced 22→20pt.
2. ~~Section 00, slide 8 (projects)~~: split into two slides (5 + 4 rows),
   font bumped 13→14pt.
3. ~~Section 03, slide 3 (agent loop diagram)~~: node distance reduced,
   boxes compacted.

Additionally, dense slides were split per the style guide ("one idea per
slide, max 3–4 short lines"):

- **Section 00, slide 6** "The problem" → quote statement slide +
  evidence slide
- **Section 01, slide 5** "Temperature" → formula slide + limits slide
- **Section 02, slide 4** "How chat works" → code block + "No memory"
  statement slide
- **Section 03, slide 4** "Tool definition" → explanation + "LLM never
  executes" statement slide

No known issues remain.

### Style guide

A comprehensive style guide was extracted from 15+ years of TJO's
presentations and lives at `~/Projects/presentations/` (private repo
`tobiasosborne/presentations`). Reference it for any future slide work.

## What to do next

- Consider adding screenshots/renders from actual projects (Lyr.jl
  volume renders, vectorfeld UI) to section 00 for visual impact
- Rehearse with timing — 45 slides for 60 minutes is comfortable
- The old `build_slides.py` can be removed once the new system is stable

## Technical notes

- `codeblock` uses `tcblisting` (verbatim) — frames containing it need
  `[fragile]`, and `\begin{codeblock}` / `\end{codeblock}` must start
  at column 0 (no leading whitespace)
- Font loading order in `beamerfontthemeTJO.sty` is critical:
  `luatex85` → `mtpro2[lite]` → `fontspec[no-math]`
- mtpro2 `[lite]` is used; `[complete]` may work with all 136 pfb files
- Font paths in `beamerfontthemeTJO.sty` are absolute to
  `/home/tobiasosborne/`; update if building on a different machine
