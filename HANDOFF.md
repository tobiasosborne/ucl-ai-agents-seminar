# Handoff: AI Agents Seminar

## What this is

Materials for a 1-hour seminar: "Large Language Models: A Physicist's
Perspective." Demystifies LLMs and AI coding agents by rebuilding the
entire stack from first principles for a physics audience.

## Current state

**Slides rebuilt twice** — first from v1 review notes, then refined
from v2 annotations (annotated PDF + two voice memo transcripts).

### Build system

- **Engine**: LuaLaTeX (TeX Live 2023)
- **Theme**: Custom `TJO` Beamer theme (`slides/beamer*TJO.sty`)
- **Fonts**: Whitney (text), MathTime Pro 2 (math), Iosevka Custom (code)
- **Colors**: Whitney Teal palette (dk2=#335B74, accent1=#1CADE4)
- **Build**: `make` at project root builds `slides/seminar.pdf`

### Sections (5 sections + freestyle, 54 slides total)

| Section | Slides | Content |
|---------|--------|---------|
| 00: Why this talk? | 14 | GPT-5.2 headline, timeline, survey, personal confession, skepticism, evidence, mental model ("unhinged MSc student"), the gap, three requirements, formalization, accelerating research, transition |
| 01: What is an LLM? | 12 | From outside looking in, stateless, nondeterministic, probability distribution, temperature/Boltzmann, limits, tokens, inference, Transformers references, closing statement |
| 02: The illusion of chat | 10 | Single API call, cURL reality, HTTP reality, how chat works, no memory, context window (bar chart), context rot, system prompt, closing statement |
| 03: From function to agent | 11 | Two ways in (web vs terminal), coding agents intro, agent loop (code + diagram), tools (+ Thorsten Ball ref), LLM never executes, two tools, production equivalence, live demo |
| 04: Beyond the loop | 8 | Real workflow, bash polling loop, multi-agent orchestration, natural endpoint, limitations (3 points), automation spectrum (Labour→Cognition) |
| 05: Live demo — FeO | — | Audience challenge: compute FeO potential energy curves (see `05-live-demo-feo/README.md`) |
| Freestyle | 1 | "Let's try something / No promises" |

### Key files

```
slides/
  beamerthemeTJO.sty          # Master theme (loads 4 sub-themes)
  beamercolorthemeTJO.sty     # Whitney Teal color palette
  beamerfontthemeTJO.sty      # Whitney + MathTime Pro 2 + Iosevka
  beamerinnerthemeTJO.sty     # Title page (with logos), section dividers, statement slides, codeblock
  beamerouterthemeTJO.sty     # Frame title + cyan accent line, footline
  preamble.tex                # Shared packages, TikZ styles, math shortcuts, logo paths
  seminar.tex                 # Unified deck (all sections)
  assets/
    gpt52pr.PNG               # GPT-5.2 press release screenshot
    Innovailia_logo_2.png     # InnovAILia logo
    Logo-Paket/               # LUH logos (RGB, CMYK, S-W, Pantone)
Makefile                      # Build system
notes/
  old/                        # v1 review materials (archived)
    REMEDIATION_PLAN.md
    REVIEW_SYNTHESIS.md
    voice_memo_transcript.txt
    annotatedslides.pdf
    seminarnotes.pdf
  seminar.pdf                 # Current annotated slides (v2)
  Seminar v2.m4a              # Voice memo (v2)
  Seminar_v2_transcript.md    # Whisper transcript (tiny model)
  Seminar v2 transcripts .txt # Second transcript (alternative ASR)
```

### Demo code

- `01-*/demo.py` — Temperature sweep showing nondeterminism (fun prompt about fictional particles)
- `02-*/single_call.py` — Raw API call with JSON payload
- `02-*/chat_loop.py` — Terminal chat showing history growth
- `02-*/curl_call.sh` — Standalone cURL demo matching the slide (`KEY` env var required)
- `03-*/agent.py` — 124-line agent with read/write tools
- `03-*/tools.py` — Sandboxed file tools
- `03-*/workspace/data.csv` — Sample data file for agent demo
- `05-*/smoke_test_feo.py` — Single-point CASSCF+NEVPT2 validation
- `05-*/compute_feo_pec.py` — Production PEC scan (SA-CASSCF(12,12) + SC-NEVPT2, X2C, DK basis)
- `05-*/plot_feo_pec.py` — PEC plotting and spectroscopic constants extraction

### Changes in v2 (from v1)

Applied from annotated slides PDF + two voice memo transcripts:

- **Survey slide**: Removed "used" prefix, added question marks
- **Personal confession**: "tried" → "experimented", added "for research", comma → colon
- **Untrustworthy**: Removed italic and quote marks
- **Evidence**: Left-aligned bullet text, teal conclusion stays centered
- **Mental model**: "deranged" → "unhinged", removed trailing comma
- **Project showcase slide**: Deleted (felt like bragging)
- **Formalization**: TODO placeholder for auto-formalization screenshot
- **From outside, looking in**: Removed "From" and "entire", left-aligned text
- **Probability distribution**: "over tokens" merged into first line, added "one after another"
- **Two ways in**: Moved "same API underneath" label below arrow to avoid overlap
- **New "Coding agents" slide**: Intro slide before agent loop ("Empowering LLMs with autonomous tool use")
- **Tools slide**: Added Thorsten Ball "How to Build an Agent" reference (ampcode.com)
- **Be skeptical**: Removed title bar, kept 3 content points as plain slide
- **New "Automation" slide**: Labour→Cognition spectrum (Calculator, CAS, Automated proving, AI tools)
- **Removed "Live coding" section divider**: Kept "Let's try something" slide
- **Context rot**: Left-aligned the three body statements
- **Timeline**: Rebalanced box spacing
- **Transformers slide**: New references slide (Vaswani, Karpathy, Alammar, PicoGPT.jl)

### Known issues / TODOs

- Insert screenshot of auto-formalization project after the formalization slide
- Add picture of transformer architecture to the Transformers references slide
- Check Nondeterministic slide for text overflow (P(next token|context) may clip)
- Minor overfull hbox warnings on some code blocks (inherent to verbatim in beamer)

### Style guide

A comprehensive style guide lives at `~/Projects/presentations/STYLE_GUIDE.md`.

### LaTeX talk script

A full written script of the talk, typeset as lecture notes using the
same `amsart` + custom style template as `~/Projects/quantum-noise-and-decoherence/`.

- **Engine**: pdfLaTeX (TeX Live 2023), Computer Modern fonts (portable)
- **Master doc**: `latex/LLMSeminar.tex`
- **Style**: `latex/talk-style.sty` (adapted from `qnd-style.sty` — same color palette, boxes, theorem environments, plus `audience` box and `listings` for code)
- **Macros**: `latex/talk-macros.sty` (LLM-specific: `\Str`, `\Tok`, `\Prob`, TikZ component styles)
- **Bibliography**: `latex/references.bib`
- **Build**: `cd latex && pdflatex LLMSeminar && bibtex LLMSeminar && pdflatex LLMSeminar && pdflatex LLMSeminar`

| Section file | Title | Content |
|---|---|---|
| `sections/sec00.tex` | Why this talk? | GPT-5.2, timeline, survey, skeptic's confession, the paradox, mental model, the gap, requirements |
| `sections/sec01.tex` | What is an LLM? | f: String→String, stateless, nondeterministic, temperature/Boltzmann, tokens, inference, summary |
| `sections/sec02.tex` | The illusion of chat | API calls, how chat works, context window, context rot, system prompt |
| `sections/sec03.tex` | From function to agent | Agent loop, tools, key insight, feedback, live demo narrative |
| `sections/sec04.tex` | Beyond the loop | Bash orchestration, multi-agent, limitations, automation spectrum |
| `sections/sec05.tex` | Applications and discussion | Examples, FeO computation, Q&A highlights (lit reviews, security, reasoning, correctness, formalisation) |

### Transcript

YouTube auto-captions downloaded and cleaned: `transcripts/talk_clean.txt`
(90 paragraphs, 87 minutes). Raw VTT at `transcripts/talk_audio.en.vtt`.

## What to do next

- Provide auto-formalization screenshot and transformer architecture image
- Rehearse with timing — 54 slides for 60 minutes is comfortable
- Prepare live demo environment (Claude Code terminal + web interface)

## Technical notes

- `codeblock` uses `tcblisting` (verbatim) — frames containing it need
  `[fragile]`, and `\begin{codeblock}` / `\end{codeblock}` must start
  at column 0 (no leading whitespace)
- Font loading order in `beamerfontthemeTJO.sty` is critical:
  `luatex85` → `mtpro2[lite]` → `fontspec[no-math]`
- Whitney font path: `/usr/share/fonts/opentype/whitney/` (system install)
- Iosevka: loaded via fontconfig name lookup (no hardcoded path)
- All demo scripts use model ID `claude-sonnet-4-6` (no date suffix)
- Logo paths defined in `preamble.tex` via `\luhlogo` and `\innovailialogo`
- TikZ decorations library (`decorations.pathreplacing`) loaded for context window brace
