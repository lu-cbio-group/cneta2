# Using AI tools

AI tools are welcome in CNETA: for code, tests, documentation and help
with writing. This page is the policy for the people who use them. The
instructions for the tools themselves are in `AGENTS.md` at the
repository root, and most tools read that file on their own (see
[Setting up your tool](#setting-up-your-tool)).

There are six rules. Each comes with the reason for it and, where it
helps, what other projects do.

## The rules

### 1. You answer for every change

Read, run and understand everything you submit, whoever or whatever
wrote it. In review, you should be able to explain any line.

- **Why.** CNETA's results go into papers. If the author cannot explain
  a change, a reviewer cannot tell a typo from a misunderstanding of the
  method.
- **Elsewhere.** Nearly every published policy starts here: the Linux
  kernel, LLVM, Fedora, NumPy, CPython.

### 2. Say so in the pull request

When an AI tool wrote a substantial part of a change (code, tests or
text), add one line to the pull request description:

```text
AI assistance: <tool> drafted <the code / the tests / this text>; I have read, run and checked it.
```

Agents that follow `AGENTS.md` add this line themselves; check that it
is true before you post. Commit messages need nothing.

- **Why.** Reviewers read AI-written code differently: they look harder
  for plausible-looking mistakes. It also records how the software was
  made.
- **Elsewhere.** NumPy, SciPy and scikit-learn ask for the same in the
  pull request description. Fedora and LLVM ask for a commit trailer
  (`Assisted-by:`). The Linux kernel used such a trailer during 2026 and
  then moved the note to the cover letter of a patch series.

### 3. Keep patient data and unpublished work out of AI tools

Do not paste, upload or let an agent open patient-derived data,
unpublished data or unpublished manuscript text. For examples and tests,
use the data in `tests/data/` or data simulated with `cnets`.

- **Why.** Most AI tools send what they read to the provider's servers.
  Patient-derived data falls under data-protection law and the group's
  data agreements; unpublished results and text are not ours to share
  yet.

### 4. Small pull requests; core changes start with an issue

One topic per pull request. A substantial change to core functionality
(the models, the likelihood, the optimiser, the statistics, input and
output formats) is tied to an issue, and the pull request links it.
Agents are told to stop and ask when a change grows beyond about 300
changed lines.

- **Why.** AI tools make large changes cheap to write, but not cheap to
  review. A small pull request can be checked against the method; a
  large one built on a misunderstanding has to be redone.

### 5. Agents open draft pull requests, people do the rest

An agent (for example the GitHub Copilot cloud agent) may open a
**draft** pull request. A person reviews it, marks it ready, writes the
comments and merges.

- **Elsewhere.** nf-core's `AGENTS.md` has the same rule. pandas, NumPy,
  scikit-learn and Astropy go further and close pull requests that an
  agent opened by itself.

### 6. AI help with writing is fine, if you stand behind it

You may let an AI tool draft or polish a pull request description, a
review reply or a documentation page, as long as you have read it and
agree with every sentence.

- **Why.** Most of us write in a second language. The content still has
  to be yours.
- **Elsewhere.** Some projects decide the other way: NumPy asks
  contributors not to use AI tools for pull request descriptions and
  comments, and pandas asks that AI-written text be quoted and never
  speak for the contributor.

## What `AGENTS.md` tells the tools

So that you know what to expect from an agent in this repository:

- the commands to build, test and check the documentation, and to run
  them before calling a change done;
- to fix a bug with a test that fails before the fix;
- to stop and ask before changing results, test reference values, CI,
  the build or dependencies, or before removing a known-failure marker;
- the words that must be used precisely (site, bin, segment, model,
  `--cn_type`, ...);
- branch names, the commit-message shape, and pull request titles,
  including "(changes results)" when the same input and seed give
  different output;
- to add the AI-assistance line of rule 2 to the pull request
  description;
- to write in plain, short English.

## Setting up your tool

| Tool | Reads | Note |
|---|---|---|
| Claude Code | `CLAUDE.md`, which imports `AGENTS.md` | Personal notes go in `CLAUDE.local.md`, which git ignores |
| GitHub Copilot: cloud agent, code review, VS Code, CLI | `AGENTS.md` | |
| GitHub Copilot Chat on github.com and in JetBrains IDEs | `.github/copilot-instructions.md` | a short version that points to `AGENTS.md` |
| Codex, Cursor | `AGENTS.md` | |
| Gemini CLI | `GEMINI.md` | add `AGENTS.md` under `context.fileName` in its settings |
| A chat in the browser | nothing | paste only what rule 3 allows |

## Good practice

- **Plan first.** Ask the tool what it intends to change, and why, before
  it changes anything.
- **One task per session.** Start a new session for a new topic.
- **Test first.** For a bug, ask for a test that fails, then the fix.
- **Read the diff.** `git diff` and `git diff --staged` before every
  commit; don't commit what you haven't read.
- **Run the tests yourself**: [Testing and CI](testing.md).
- **Check the science.** AI tools explain methods confidently and are
  sometimes wrong. Check against the user guide and the papers, and ask
  when unsure.
- **Teach the file.** When an agent makes the same mistake twice, add a
  line to `AGENTS.md` in a pull request.

## What other projects do

Checked on 6 November 2026.

| Project | AI use | Disclosure |
|---|---|---|
| [Linux kernel](https://docs.kernel.org/process/generated-content.html) | allowed; the person who submits is responsible | note in the patch series' cover letter |
| [LLVM](https://llvm.org/docs/AIToolPolicy.html) | allowed with a person in the loop; no agents acting alone | label substantial tool output, for example `Assisted-by:` |
| [Fedora](https://docs.fedoraproject.org/en-US/council/policy/ai-contribution-policy/) | allowed; the contributor is accountable | required for significant parts; `Assisted-by:` |
| [NumPy](https://numpy.org/devdocs/dev/ai_policy.html), [SciPy](https://docs.scipy.org/doc/scipy/dev/conduct/ai_policy.html) | allowed; you must understand it; no pull requests by agents | required in the pull request description |
| [scikit-learn](https://scikit-learn.org/dev/developers/contributing.html) | no fully automated pull requests | in the pull request description |
| [pandas](https://pandas.pydata.org/docs/dev/development/contributing.html) | allowed responsibly; AI must not speak for you | tool and model named |
| [nf-core/modules](https://github.com/nf-core/modules/blob/HEAD/AGENTS.md) | agents may open draft pull requests only | "Generated by" note |
| [Bioconductor](https://contributions.bioconductor.org/ai-policy-third-party.html) | allowed; you are responsible | in the submission issue and pull request |
| [JOSS](https://joss.readthedocs.io/en/latest/policies.html) | allowed | required section in the paper |
| [QEMU](https://www.qemu.org/docs/master/devel/code-provenance.html), [Gentoo](https://wiki.gentoo.org/wiki/Project:Council/AI_policy) | not accepted | — |

Two studies shaped `AGENTS.md`:

- GitHub looked at more than 2,500 `AGENTS.md` files
  ([blog post, November 2025](https://github.blog/ai-and-ml/github-copilot/how-to-write-a-great-agents-md-lessons-from-over-2500-repositories/)).
  The useful ones put the commands first and sort their rules into
  "always", "ask first" and "never". They start small and add a rule
  when an agent makes a mistake.
- Gloaguen and colleagues (ETH Zurich and LogicStar.ai) measured whether
  such files help ([arXiv 2602.11988](https://arxiv.org/abs/2602.11988), 2026).
  Files written by people helped a little; files generated by an AI tool
  made results slightly worse and cost more; overviews of the repository
  did not help agents find the right files. So `AGENTS.md` stays short
  and hand-written, and leaves out what the README and these pages
  already say.

The format itself is maintained by the
[Agentic AI Foundation](https://agents.md/) under the Linux Foundation.

## Changing these rules

Propose a change in a pull request that edits this page and `AGENTS.md`
together, so the policy and the instructions stay in step.
