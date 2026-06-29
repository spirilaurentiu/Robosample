For each .md file in references/papers/*.md, read the full paper and fill in references/index.yaml (or whatever format you chose to reference all papers faster - this will be consumed by you). Use ASCII only. Here are some tags with explainations:

- `summary`: 2-4 sentences, the gist and why it matters to our algorithm
- `tags`: short topic keywords
- `key_terms`: the concepts/definitions the paper introduces
- `key_equations`: the named equations our program would need to obey, in short form. Some papers are OCR-garbled. LaTeX suffers the most from it, mainly from unclosed equations.
- `invariants`: the obligations those equations impose on the code (for Rule 8 tests)
- `path`: path to full .md file

Example:

---
`title`: "Hybrid Monte Carlo"
`authors`: ["Duane, S.", "Kennedy, A. D.", "Pendleton, B. J.", "Roweth, D."]
`year`: 1987
`venue`: "Physics Letters B"
`doi`: "10.1016/0370-2693(87)91197-X"
`tags`: [HMC, molecular-dynamics, sampling, detailed-balance, acceptance]
`path`: `references/papers/0021-9991_2877_2990098-5.md`
`key_terms`:

- Hybrid Monte Carlo (MD trajectory as a global MH proposal)
- Guidance vs acceptance Hamiltonian
- `key_equations`:
  - Acceptance: A = min(1, exp(-beta * dH))
- `invariants`:
  - Proposal must be reversible and volume-preserving for exact MH.
  - Guidance Hamiltonian may differ from the acceptance Hamiltonian without bias.
- summary:
  Introduces hybrid Monte Carlo: MD trajectories as global proposals accepted by a Metropolis test
  on the full Hamiltonian. Source of the distinct guidance/acceptance Hamiltonian idea Robosample
  relies on. Relevant to: proposal reversibility, volume preservation, exact-dH acceptance.

---
