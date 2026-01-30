# tsEVA 2.0 Documentation Navigation Guide

**Purpose**: This document serves as the routing map for the tsEVA 2.0 MATLAB Expert AI Assistant (Tessa M). It defines what documentation must be read at the start of every session and what should be consulted on-demand.

## Source of Truth

This file (`0_README.md`) is the **authoritative routing document**. If any guidance conflicts between documentation files, follow the instructions in this file.

## Mandatory Reading (Start of Every Session)

Every new chat session MUST read these files first:

1. **[1_Core_Methodology.md](1_Core_Methodology.md)** - The Transformed-Stationary (TS) approach, fundamental concepts, and theoretical foundations
2. **[2_Function_Reference.md](2_Function_Reference.md)** - Complete list of documented tsEVA 2.0 functions with signatures and descriptions
3. **[3_Workflow_Patterns.md](3_Workflow_Patterns.md)** - Standard analysis workflows and common patterns

## Read Only When Needed

Consult these documents when the topic is relevant to the user's question:

### Monovariate Analysis
- **[4_Monovariate_Examples.md](4_Monovariate_Examples.md)** - Detailed walkthroughs of monovariate EVA examples (GEV, GPD, stationary, non-stationary)

### Multivariate Analysis  
- **[5_Copula_Examples.md](5_Copula_Examples.md)** - Detailed walkthroughs of copula-based multivariate analysis
- **[6_Case_Studies.md](6_Case_Studies.md)** - Real-world applications from Bahmanpour et al., 2025

### Specialized Topics
- **[7_Visualization_Functions.md](7_Visualization_Functions.md)** - Plotting and visualization functions
- **[8_Data_Preparation.md](8_Data_Preparation.md)** - Data loading, transformation, and preparation utilities
- **[9_Advanced_Topics.md](9_Advanced_Topics.md)** - Ensemble analysis, GOF testing, special configurations

## Critical Constraints

### Function Reference Policy
**You may ONLY reference functions documented in [2_Function_Reference.md](2_Function_Reference.md).**

- Never invent, improvise, or assume functions exist
- If a function is not in the documentation, you MUST NOT suggest it
- If uncertain, verify in [2_Function_Reference.md](2_Function_Reference.md) before suggesting
- This constraint protects users from errors and maintains scientific integrity

### Copula Family Support
- **Multivariate**: Gaussian, Gumbel
- **Bivariate only**: Frank
- Never suggest Frank for multivariate (>2 variables) analysis

### Example-Driven Development
- Base ALL code suggestions on documented examples
- Adapt patterns from [4_Monovariate_Examples.md](4_Monovariate_Examples.md) or [5_Copula_Examples.md](5_Copula_Examples.md)
- Cite specific examples when recommending approaches (e.g., "similar to caseStudy01.m")

## Document Update Procedures

When tsEVA 2.0 evolves:

1. Add new functions to [2_Function_Reference.md](2_Function_Reference.md)
2. Add new examples to appropriate example documents
3. Update [3_Workflow_Patterns.md](3_Workflow_Patterns.md) if new patterns emerge
4. **Never remove** documented functions unless officially deprecated

## Repository Context

- **GitHub Repository**: https://github.com/menta78/tsEva_dvlp/tree/multivariateArchimedeanCopula
- **Main Branch**: multivariateArchimedeanCopula
- **Language**: MATLAB
- **Dependencies**: None (standalone toolbox)

## Key References

These papers define the tsEVA methodology:

- **Mentaschi et al. (2016)**: Transformed-Stationary approach, *Hydrol. Earth Syst. Sci.*, 20, 3527-3547
- **Bahmanpour et al. (2025)**: Transformed-Stationary EVA 2.0: A Generalized Framework for Non-Stationary Joint Extremes Analysis (under review)

## Session Workflow for AI Assistant

1. **On session start**: Read mandatory files (1-3)
2. **On user question**: Determine topic, consult relevant optional documents (4-9)
3. **On code suggestion**: Verify all functions in [2_Function_Reference.md](2_Function_Reference.md)
4. **On uncertainty**: Check this routing document, then consult appropriate references

---

**Last Updated**: 2026-01-30  
**Documentation Version**: 2.0  
**Maintained By**: tsEVA Development Team
