# Scientific Document Review for LLM-Assisted Molecular Simulation System

You are reviewing a technical document that defines the assumptions, components, and reasoning structure of a molecular simulation and sampling software system with integrated LLM-assisted coding and scientific reasoning. The same document will also serve as reference to future onboarding developers (engineers, biochemists, physicists etc).

## Objective

Evaluate and improve the document with respect to:

* Scientific correctness and internal consistency
* Clarity of definitions and conceptual precision
* Suitability for LLM-assisted reasoning and decomposition
* Robustness under limited context and competing information sources
* Resistance to prompt dilution and information overload
* Known System Failure Modes (must actively assess)

The system is known to suffer from:

* Limited context handling: important dependencies may be lost across sections
* Unstructured content drift: definitions and assumptions may not remain stable across modules
* Information competition: conflicting definitions across components
* Model overwhelm: excessive simultaneous constraints or overloaded reasoning steps
* Prompt dilution: key constraints become less influential as document length increases

Explicitly identify where each of these occurs or could occur.

## Review Tasks

1. Structural Audit: Identify whether the document has a clear hierarchical structure:

    * Axioms / assumptions
    * Definitions
    * Computational components
    * Interfaces between components

    Flag any sections where structure is implicit or ambiguous.

2. Conceptual Consistency Check: Detect contradictions in:

    * Definitions
    * Coordinate systems
    * Probabilistic assumptions
    * Sampling rules

    Highlight competing or duplicated definitions.

3. LLM-Compatibility Assessment: Evaluate whether the document can be reliably used by an LLM to:

    * Decompose tasks into subproblems
    * Maintain correct variable/state tracking
    * Avoid cross-contamination between components
    * Preserve invariants across reasoning steps

4. Identify where ambiguity would lead to reasoning failure:

    * Information Density and Cognitive Load

        * Identify overly dense sections
        * Detect where too many constraints are introduced simultaneously
        * Flag where simplification or modularization is needed
        * Improvement Requirements

5. For every issue found:

    * Provide a precise correction
    * Rewrite the problematic section in a cleaner form
    * If necessary, propose a modular decomposition of the concept

6. Additionally:

    * Suggest a minimal “core axioms” version of the document
    * Suggest a LLM execution-friendly version (highly structured, stepwise, low ambiguity)
    * Output Format
    * Executive Summary (max 10 lines)
    * Critical Issues (grouped by failure mode)
    * Line-by-line or section-by-section fixes (if applicable)
    * Rewritten key sections
    * Proposed improved structure of the document
    * Minimal core axioms version
    * LLM-optimized version for execution
    * Additional Instruction

Do not be permissive. Assume subtle ambiguity is a failure. Prefer restructuring over patching when uncertainty propagates across sections.
