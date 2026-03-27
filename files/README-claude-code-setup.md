# Claude Code Setup for temperheic

This directory contains the files to add to your temperheic repository root for a Claude Code workflow. Copy these into your existing package directory.

## File Tree

```
temperheic/                          # Your existing package root
│
├── CLAUDE.md                        # ★ Main project context — loaded every session
│                                    #   Project overview, build commands, conventions,
│                                    #   @imports to docs/ for on-demand reference
│
├── .claude/
│   ├── rules/
│   │   ├── r-package-conventions.md # Always loaded: R code style, testing, roxygen rules
│   │   └── thermal-hydrology-domain.md # Always loaded: physical parameters, sign conventions,
│   │                                #   Meacham Creek site context
│   └── commands/
│       ├── check-package.md         # /check-package — run test + check pipeline
│       └── review-function.md       # /review-function fit_ols — verify against source paper
│
├── docs/                            # Project knowledge (referenced via @imports in CLAUDE.md)
│   ├── design-philosophy.md         # Package design principles
│   ├── tidyverse-fp-guidelines.md   # R coding conventions (copy from project knowledge)
│   ├── methods-review.md            # Full theoretical framework (copy from project knowledge)
│   ├── goals-checklist.md           # Current development phases (copy from project knowledge)
│   ├── architecture.md              # R toolkit vision, Layers 1-5 (copy from project knowledge)
│   ├── project-instructions.md      # ReacTran patterns (copy from project knowledge)
│   └── refs/
│       ├── luce-2013-reference.md   # ★ Reference card: η method, key equations, UQ
│       ├── luce-2017-reference.md   # ★ Reference card: BC independence, multi-freq diagnostic
│       ├── van-kampen-2022-reference.md # ★ Reference card: LPMLEn method overview
│       └── bertagnoli-2024-reference.md # ★ Reference card: iFLOW, benchmarking test cases
│
├── refs/                            # Full PDFs — read by Claude Code on demand
│   ├── luce_et_al_2013.pdf
│   ├── luce_et_al_2017.pdf
│   ├── van_kampen_et_al_2022.pdf
│   └── bertagnoli_et_al_2024.pdf
│
├── R/                               # (existing) Package source code
├── tests/                           # (existing) testthat tests
├── inst/extdata/                    # (existing) Meacham Creek datasets
├── man/                             # (existing) roxygen-generated docs
├── DESCRIPTION                      # (existing)
└── NAMESPACE                        # (existing)
```

## How it works

### Always in context (every session):
- `CLAUDE.md` — project overview, conventions, build commands
- `.claude/rules/*.md` — R style rules, physical parameter ranges, sign conventions

### Loaded on demand (when Claude needs them):
- `docs/*.md` — referenced via `@` imports in CLAUDE.md; Claude reads when relevant
- `docs/refs/*.md` — lightweight reference cards with key equations and implementation mapping

### Loaded when needed for exact details:
- `refs/*.pdf` — full research papers; Claude reads these when it needs exact derivations, figures, or table values that aren't in the reference cards

### Slash commands:
- `/check-package` — run the full test + check pipeline
- `/review-function {name}` — verify a function against source paper equations

## Setup instructions

1. Copy `CLAUDE.md` to your temperheic package root
2. Copy the `.claude/` directory to your package root
3. Copy the `docs/` directory to your package root (or adapt paths if you already have a docs/ dir)
4. Create a `refs/` directory and place your PDF files there (rename to match the reference cards)
5. Copy your existing project knowledge markdown files into `docs/`:
   - `temperheic-design-philosophy.md` → `docs/design-philosophy.md` (provided)
   - `Tidyverse___Functional_Programming_Guidelines_for_temperheic.md` → `docs/tidyverse-fp-guidelines.md`
   - `streambed-flux-methods-review.md` → `docs/methods-review.md`
   - `temperheic-goals-checklist-2026-02-final.md` → `docs/goals-checklist.md`
   - `r-applied-thermal-hydrology-toolkit-architecture.md` → `docs/architecture.md`
   - `PROJECT_INSTRUCTIONS_temp_signal_modeling.md` → `docs/project-instructions.md`
6. Add to `.gitignore`: `.claude/settings.local.json` (if it appears)
7. Run `/init` in Claude Code to let it learn your codebase, then review/merge with the provided CLAUDE.md

## What to add to .gitignore

```
# Claude Code local settings (personal, not shared)
.claude/settings.local.json
```

## Token budget comparison

| Content | Tokens (approx) | When loaded |
|---------|-----------------|-------------|
| CLAUDE.md | ~800 | Every session |
| .claude/rules/ (2 files) | ~600 | Every session |
| Reference card (each) | ~500-800 | On demand |
| Full project knowledge doc (each) | ~2,000-4,000 | On demand |
| Full PDF paper (each) | ~30,000-60,000 | On demand |

The reference cards are the key innovation: ~500 tokens vs ~40,000 tokens for a full PDF, covering ~90% of what you need in day-to-day development. Claude reads the full PDF only when you need exact details.
