# Accelerating Physics Discovery with Machine Intelligence

**1st Grand Challenge Conference** --- Berlin University Alliance
23--25 March 2026, Berlin

30-minute talk adapted from the AI Agents Seminar material.

## Build

```bash
cd slides
make
```

## Structure

27 slides in five parts:

1. **The Landscape** --- GPT-5.2, timeline, thesis (capability vs trust)
2. **What is an LLM?** --- f : String -> String, nondeterminism, temperature, context
3. **From Function to Agent** --- the agent loop, tool use
4. **Accelerating Discovery** --- multi-agent research, adversarial verification, convergence of AI + quantum + formal methods
5. **The Harder Problem** --- why LLMs fail at research, real AI harms, automation spectrum, human in the loop

## Differences from the seminar

- Cut from 54 to 27 slides (~50% reduction)
- Removed: cURL/HTTP demos, bash polling loop, live coding section, audience vote, detailed code demos
- Added: adversarial verification (from blog), convergence diagram (AI + quantum + formal methods), social impact slide (real AI harms), summary slide
- Reframed multi-agent diagram for research (conjectures, verification, experiments)
- Merged context window + context rot into one slide
- Title and framing aligned to conference abstract
