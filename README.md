# AI Coding Agents: A Reductive Account for Physicists

Materials for a 1-hour seminar at UCL Physics. The goal is to demystify
AI coding agents by showing that every part of the stack reduces to
simple, composable primitives that physicists already have intuition for.

## Setup

```bash
git clone https://github.com/<you>/ucl-ai-agents-seminar.git
cd ucl-ai-agents-seminar
pip install -r requirements.txt
export ANTHROPIC_API_KEY=sk-ant-...
```

## Regenerate slides

```bash
make
```

This compiles all four LaTeX Beamer decks (via LuaLaTeX) and copies the
PDFs to each section directory. Requires LuaLaTeX with Whitney, MathTime
Pro 2, and Iosevka fonts installed.

To build a single section: `make s01`, `make s02`, `make s03`, `make s04`.

The slide source lives in `slides/` with a custom Beamer theme (`TJO`)
matching the Whitney Teal design system.

## Structure

| Directory | What it covers |
|-----------|---------------|
| `01-stateless-nondeterministic-function/` | An LLM is a stateless nondeterministic function f : String -> String |
| `02-illusion-of-chat/` | Chat is a client-side loop. Memory is an illusion. |
| `03-primitive-agent/` | The agent loop: LLM + tools + history |
| `04-audience-vote/` | Live audience vote on what to build next |

Each section directory has its own README with instructions for running
its demos and what the audience should observe.

## Requirements

- Python 3.10+
- An Anthropic API key (for the demo scripts)
- LuaLaTeX (TeX Live 2023+) with Beamer, tcolorbox, TikZ
- Whitney font family (OTF)
- MathTime Pro 2 (Type 1)
- Iosevka Custom (TTF, monospace)
